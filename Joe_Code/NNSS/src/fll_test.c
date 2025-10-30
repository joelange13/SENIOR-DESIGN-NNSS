// fll_stream_apply.c
// Streaming band-edge FLL that outputs a frequency-corrected file.
// Build:  g++ -O2 -std=c++17 fll_stream_apply.c -lm -o fll_stream
// Usage:  ./fll_stream "<path-to.cfile>" [Fs_Hz]
//
// Notes:
// - Input/Output: interleaved float32 (I,Q,I,Q,...).
// - Parameters: SPS=8, alpha=0.5, Ntaps=401 (tweak if desired).
// - Error: e = (Pu - Pl) / (Pu + Pl + eps)  [matches your probe polarity].
// - Loop: PI with decimated updates and warm-up to avoid transient bias.

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

typedef struct { float re, im; } cplx;

static inline cplx cexpj(float th){ return (cplx){cosf(th), sinf(th)}; }
static inline cplx cmul(cplx a, cplx b){ return (cplx){a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re}; }
static inline float cabs2(cplx a){ return a.re*a.re + a.im*a.im; }
static inline float sinc(float x){ return (fabsf(x)<1e-8f)?1.0f:sinf((float)M_PI*x)/((float)M_PI*x); }

/* ---------- Hamming-windowed low-pass prototype (DC-normalized) ---------- */
static void lp_sinc_hamming(float fc_cyc_per_samp, int Ntaps, float *h){
    int M = (Ntaps-1)/2;
    for(int n=0;n<Ntaps;n++){
        float m = (float)(n - M);
        float w = 0.54f - 0.46f*(float)cos(2.0*M_PI*(double)n/(double)(Ntaps-1));
        h[n] = 2.f*fc_cyc_per_samp * sinc(2.f*fc_cyc_per_samp*m) * w;
    }
    float s=0.f; for(int i=0;i<Ntaps;i++) s += h[i];
    if (s != 0.f) for(int i=0;i<Ntaps;i++) h[i] /= s;
}

/* ---------- Complex band-edge filters at +/- f_edge ---------- */
static void design_bandedges(int sps, float alpha, int Ntaps, cplx *hu, cplx *hl){
    float df = alpha/(2.f*sps);              // slice roll-off width
    float *hlp = (float*)calloc((size_t)Ntaps, sizeof(float));
    lp_sinc_hamming(df, Ntaps, hlp);

    float f_edge = 1.f/(2.f*sps);            // half symbol rate (cycles/sample)
    int M = (Ntaps-1)/2;
    for(int n=0;n<Ntaps;n++){
        float th = 2.f*(float)M_PI*f_edge*(float)(n - M);
        float c = cosf(th), s = sinf(th);
        hu[n].re = hlp[n]*c;  hu[n].im = hlp[n]*s;   // +f_edge
        hl[n].re = hlp[n]*c;  hl[n].im =-hlp[n]*s;   // -f_edge
    }
    free(hlp);
}

/* ---------- Streaming complex FIR ---------- */
typedef struct {
    int Nt;
    cplx *taps;
    cplx *d;     // delay line [Nt]
} cfir_t;

static void cfir_init(cfir_t *f, const cplx *taps, int Nt){
    f->Nt = Nt;
    f->taps = (cplx*)malloc((size_t)Nt*sizeof(cplx));
    f->d    = (cplx*)calloc((size_t)Nt,sizeof(cplx));
    memcpy(f->taps, taps, (size_t)Nt*sizeof(cplx));
}
static inline cplx cfir_push(cfir_t *f, cplx x){
    memmove(&f->d[1], &f->d[0], (size_t)(f->Nt-1)*sizeof(cplx));
    f->d[0] = x;
    float re=0.f, im=0.f;
    for(int k=0;k<f->Nt;k++){
        re += f->d[k].re*f->taps[k].re - f->d[k].im*f->taps[k].im;
        im += f->d[k].re*f->taps[k].im + f->d[k].im*f->taps[k].re;
    }
    return (cplx){re, im};
}
static void cfir_free(cfir_t *f){ free(f->taps); free(f->d); }

/* ---------- Load & save interleaved cf32 ---------- */
static int load_cfile(const char *path, cplx **out){
    FILE *f = fopen(path,"rb");
    if(!f){ perror("open"); return -1; }
    fseek(f,0,SEEK_END); long bytes = ftell(f); fseek(f,0,SEEK_SET);
    if(bytes < 0){ fclose(f); return -1; }
    size_t Nfloat = (size_t)bytes/sizeof(float);
    if((Nfloat & 1u) != 0u){ fclose(f); fprintf(stderr,"File length not even floats.\n"); return -1; }
    size_t Nsamp = Nfloat / 2;
    float *buf = (float*)malloc(Nfloat*sizeof(float));
    size_t rd = fread(buf, sizeof(float), Nfloat, f); fclose(f);
    if(rd != Nfloat){ free(buf); fprintf(stderr,"Read error.\n"); return -1; }
    cplx *x = (cplx*)malloc(Nsamp*sizeof(cplx));
    for(size_t n=0;n<Nsamp;n++){ x[n].re = buf[2*n+0]; x[n].im = buf[2*n+1]; }
    free(buf); *out = x; return (int)Nsamp;
}
static int save_cfile(const char *path, const cplx *x, int N){
    FILE *f = fopen(path,"wb"); if(!f){ perror("write"); return -1; }
    for(int n=0;n<N;n++){
        float r=x[n].re, i=x[n].im;
        if(fwrite(&r,sizeof(float),1,f)!=1 || fwrite(&i,sizeof(float),1,f)!=1){ fclose(f); return -1; }
    }
    fclose(f); return 0;
}

/* ---------- FLL state ---------- */
typedef struct {
    float Kp, Ki;       // loop gains
    float freq_cps;     // cycles/sample
    float phase;        // radians
    float Pu, Pl;       // smoothed edge powers
    float beta;         // power smoother
} fll_t;

static void fll_init(fll_t *l, float Kp, float Ki, float beta){
    l->Kp=Kp; l->Ki=Ki; l->beta=beta;
    l->freq_cps=0.f; l->phase=0.f;
    l->Pu=0.f; l->Pl=0.f;
}

int main(int argc, char **argv){
    if(argc < 2){
        fprintf(stderr, "Usage: %s <path-to.cfile> [Fs_Hz]\n", argv[0]);
        return 1;
    }
    const char *inpath = argv[1];
    double Fs_Hz = (argc >= 3) ? atof(argv[2]) : 0.0;

    // Band-edge params (match your setup)
    const int   SPS   = 8;
    const float ALPHA = 0.5f;
    const int   NTAPS = 401;

    // Loop parameters (stable defaults)
    const float Kp   = 0.004f;     // proportional
    const float Ki   = 0.00004f;   // integral
    const float beta = 0.0015f;    // power smoothing
    const int   DECIM= 8;          // PI updates every DECIM samples
    const int   warm = NTAPS + 12*SPS;

    // Expected CFO guard (cycles/sample). Adjust if needed.
    const float CPS_CLAMP = 0.05f;

    // Load input
    cplx *x; int N = load_cfile(inpath, &x);
    if(N <= 0){ fprintf(stderr,"Failed to load: %s\n", inpath); return 1; }
    printf("Loaded %d complex samples.\n", N);

    // Design band-edge filters
    cplx *hu = (cplx*)malloc((size_t)NTAPS*sizeof(cplx));
    cplx *hl = (cplx*)malloc((size_t)NTAPS*sizeof(cplx));
    design_bandedges(SPS, ALPHA, NTAPS, hu, hl);

    // Init streaming FIRs
    cfir_t fu, fl; cfir_init(&fu, hu, NTAPS); cfir_init(&fl, hl, NTAPS);

    // Init FLL
    fll_t loop; fll_init(&loop, Kp, Ki, beta);
    // Optional seed for faster pull-in on this known file:
    loop.freq_cps = 0.0100f;

    // Output (corrected)
    cplx *y = (cplx*)malloc((size_t)N*sizeof(cplx));

    // Decimated update accumulators
    int decim_cnt = 0;
    double e_accum = 0.0;

    // For reporting
    double mean_en = 0.0; int cnt_en = 0;

    for(int n=0; n<N; n++){
        // Mix by current NCO to get corrected sample
        cplx ejm = cexpj(-loop.phase);
        cplx ycorr = cmul(x[n], ejm);

        // Update band-edge filters on corrected sample
        cplx yu = cfir_push(&fu, ycorr);
        cplx yl = cfir_push(&fl, ycorr);

        // Smooth powers always
        float pu = cabs2(yu), pl = cabs2(yl);
        loop.Pu = (1.f - loop.beta)*loop.Pu + loop.beta*pu;
        loop.Pl = (1.f - loop.beta)*loop.Pl + loop.beta*pl;

        // Normalized error (Pu - Pl) to match your probe
        float num = loop.Pu - loop.Pl;
        float den = loop.Pu + loop.Pl + 1e-12f;
        float e   = num / den;

        if (n < warm){
            // No PI updates during warmup; just advance by freq
            loop.phase += 2.f*(float)M_PI * loop.freq_cps;
        } else {
            // Accumulate error for decimated PI update
            e_accum += (double)e;
            decim_cnt++;

            if (decim_cnt == DECIM){
                float e_avg = (float)(e_accum / (double)DECIM);
                e_accum = 0.0; decim_cnt = 0;

                // Soft clamp e to avoid spikes
                if (e_avg >  0.5f) e_avg =  0.5f;
                if (e_avg < -0.5f) e_avg = -0.5f;

                // PI update
                loop.freq_cps += loop.Ki * e_avg;            // integral
                if (loop.freq_cps >  CPS_CLAMP) loop.freq_cps =  CPS_CLAMP;
                if (loop.freq_cps < -CPS_CLAMP) loop.freq_cps = -CPS_CLAMP;

                loop.phase    += 2.f*(float)M_PI * loop.freq_cps
                                 + loop.Kp * e_avg;           // proportional
            } else {
                // Between updates, advance only by freq
                loop.phase += 2.f*(float)M_PI * loop.freq_cps;
            }

            // Post-warmup mean error (for report)
            mean_en += (double)e;
            cnt_en++;
        }

        // Keep phase bounded (optional)
        if (loop.phase >  1e6f) loop.phase -= 1e6f;
        if (loop.phase < -1e6f) loop.phase += 1e6f;

        y[n] = ycorr;

        if((n % 8000) == 0){
            double freq_cps = loop.freq_cps;
            if(Fs_Hz > 0.0){
                printf("n=%7d  freq_est = %+0.6f cps  (%+0.3f Hz)\n",
                       n, freq_cps, freq_cps*Fs_Hz);
            } else {
                printf("n=%7d  freq_est = %+0.6f cps\n", n, freq_cps);
            }
        }
    }

    // Save corrected file
    char outpath[1024];
    snprintf(outpath, sizeof(outpath), "%s_derot.cfile", inpath);
    if(save_cfile(outpath, y, N)==0){
        printf("Wrote corrected file: %s\n", outpath);
    } else {
        fprintf(stderr,"Failed to write corrected file.\n");
    }

    // Final report
    double final_cps = loop.freq_cps;
    if(Fs_Hz > 0.0){
        printf("\nFinal estimates:\n  freq_est_final = %+0.9f cps  (%+0.6f Hz)\n",
               final_cps, final_cps*Fs_Hz);
    } else {
        printf("\nFinal estimates:\n  freq_est_final = %+0.9f cps\n", final_cps);
    }
    if(cnt_en) printf("  mean normalized edge error (post-warmup) = %0.6g\n", mean_en / (double)cnt_en);

    // Cleanup
    cfir_free(&fu); cfir_free(&fl);
    free(hu); free(hl);
    free(x); free(y);
    return 0;
}
