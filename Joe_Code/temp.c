// pfb_sync.c
// Build:  gcc -O2 -std=c11 -lm -o pfb_sync pfb_sync.c
// Usage:  ./pfb_sync in.cf32 out.cf32 [--sps N] [--span K] [--alpha A] [--L P]
//                         [--loopbw BW] [--damp ZETA] [--orelim R] [--noagc] [--notrim] [--quiet]

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

/* ---------------------- your RRC generator (as given) ---------------------- */
void rrc_generate_taps(float alpha, int sps, int span, float *taps_out){
    int num_taps = span * sps + 1; 
    int mid = num_taps/2;
    double a = (double)alpha;
    double T = 1.0;
    double Ts = T/(double)sps;
    double eps = 1e-12;
    double sum = 0.0;
    for(int n = 0; n < num_taps; n++){
        int k = n - mid;
        double t = k * Ts;
        double x = 4.0 * a * t/T;
        double h;
        if(fabs(t) < eps){
            h = 1.0 + a * (4.0/M_PI - 1.0);
        }else if(a > 0.0 && fabs(fabs(t) - T/(4.0*a)) < (2.0*Ts)){
            double theta = M_PI/(4.0*a);
            h = (a/sqrt(2.0)) * ((1.0 + 2.0/M_PI)*sin(theta) + (1.0 - 2.0/M_PI)*cos(theta));
        }else{
            double num1 = cos((1.0 + a)*M_PI*t/T);
            double num2 = (T/(4.0*a*t)) * sin((1.0 - a)*M_PI*t/T);
            double den1 = M_PI*t/T;
            double den2 = 1.0 - x*x;
            h = (4.0*a) * (num1 + num2)/(den1 * den2); 
        }

        taps_out[n] = (float)h;
        sum += h*h; 
    }
    double scale;
    if (sum > 0.0) scale = 1.0/sqrt(sum);
    else           scale = 1.0;
    for(int n = 0; n < num_taps; n++){
        taps_out[n] = (float)(taps_out[n]*scale);
    } 
}

/* ------------------------- polyphase bank utilities ------------------------ */
typedef struct {
    int L;      // number of phases
    int tpp;    // taps per phase
    float *b;   // flattened [L][tpp], phase-major
} PolyBank;

// Build bank by padding proto to multiple of L, then de-interleaving: h_p[k] = h[p + k*L]
static PolyBank pfb_from_proto(const float *proto, int proto_len, int L){
    PolyBank pb = {0,0,NULL};
    if (L <= 0 || proto_len <= 0) return pb;

    int padded = proto_len;
    if (padded % L != 0) {
        padded = proto_len + (L - (proto_len % L));
    }
    int tpp = padded / L;

    float *buf = (float*)calloc((size_t)L * tpp, sizeof(float));
    if (!buf) return pb;

    for (int p = 0; p < L; p++){
        for (int k = 0; k < tpp; k++){
            int src = p + k*L;
            float v = 0.0f;
            if (src < proto_len) v = proto[src];
            buf[p*tpp + k] = v;
        }
    }
    pb.L = L; pb.tpp = tpp; pb.b = buf;
    return pb;
}

static void pfb_free(PolyBank *pb){
    if (!pb) return;
    if (pb->b){ free(pb->b); pb->b = NULL; }
    pb->L = 0; pb->tpp = 0;
}

static inline void pfb_dot_phase(const PolyBank *pb,
                                 const float *xIQ, size_t center_idx,
                                 int phase, float *outI, float *outQ)
{
    const float *taps = &pb->b[phase * pb->tpp];
    int half = (pb->tpp - 1)/2;

    double accI = 0.0, accQ = 0.0;
    int k;
    for (k = -half; k <= half; k++){
        long ii = (long)center_idx + k;
        if (ii < 0) continue;
        // caller must ensure bounds for upper end; we protect below
        const float *s = &xIQ[(size_t)ii * 2];
        double hk = (double)taps[k + half];
        accI += hk * (double)s[0];
        accQ += hk * (double)s[1];
    }
    *outI = (float)accI;
    *outQ = (float)accQ;
}

/* -------------------------- PFB clock sync state --------------------------- */
typedef struct {
    // loop state
    double omega, omega_min, omega_max;  // samples/symbol estimate
    double mu;                           // fractional timing [0,1)
    double gain_mu, gain_omega;          // loop gains (alpha, beta)

    // bank
    int sps_nom;         // nominal sps
    int span_syms;       // RRC span (symbols)
    int L;               // phases
    int proto_len;       // original taps length
    float *proto;        // owned taps
    PolyBank bank;

    // memory for Gardner
    float yk_prev_I, yk_prev_Q;

    // options
    int verbose;
} PFBClockSync;

static PFBClockSync* pfbcs_create(float alpha, int sps_nom, int span_syms, int L,
                                  double loop_bw, double damping, double omega_rel_lim,
                                  int verbose)
{
    if (sps_nom <= 0 || span_syms <= 0 || L <= 0) return NULL;

    PFBClockSync *cs = (PFBClockSync*)calloc(1, sizeof(PFBClockSync));
    if (!cs) return NULL;

    cs->sps_nom = sps_nom;
    cs->span_syms = span_syms;
    cs->L = L;
    cs->verbose = verbose;

    // build prototype RRC taps (allocate here)
    int N = span_syms * sps_nom + 1;
    cs->proto = (float*)malloc(sizeof(float) * (size_t)N);
    if (!cs->proto){ free(cs); return NULL; }
    rrc_generate_taps(alpha, sps_nom, span_syms, cs->proto);
    cs->proto_len = N;

    // polyphase bank
    cs->bank = pfb_from_proto(cs->proto, cs->proto_len, L);

    // 2nd order loop filter (Gardner slope constants)
    // These constants work well in practice; you can tweak.
    double Kp = 2.7;
    double Ki = 3.0;
    double denom = 1.0 + 2.0*damping*loop_bw + loop_bw*loop_bw;
    cs->gain_mu    = (4.0*damping*loop_bw) / (denom * Kp);
    cs->gain_omega = (4.0*loop_bw*loop_bw) / (denom * Ki);

    cs->omega = (double)sps_nom;
    cs->omega_min = cs->omega * (1.0 - omega_rel_lim);
    cs->omega_max = cs->omega * (1.0 + omega_rel_lim);

    cs->mu = 0.0;
    cs->yk_prev_I = 0.0f;
    cs->yk_prev_Q = 0.0f;

    return cs;
}

static void pfbcs_free(PFBClockSync *cs){
    if (!cs) return;
    if (cs->proto) free(cs->proto);
    pfb_free(&cs->bank);
    free(cs);
}

// returns number of output symbols
static size_t pfbcs_process(PFBClockSync *cs,
                            const float *xIQ, size_t nin_complex,
                            float *outIQ)
{
    const int tpp = cs->bank.tpp;
    const int half = (tpp - 1)/2;

    // leave room for center and mid windows
    size_t in_idx = (size_t)(half + (int)floor(0.5 * cs->omega) + 2);
    size_t end_idx = nin_complex - (size_t)half - 2;
    size_t out_count = 0;

    while (in_idx < end_idx) {
        // choose phase for current mu
        double pf = cs->mu * (double)cs->L;
        int phase = (int)floor(pf + 0.5);
        if (phase < 0) phase = 0;
        if (phase >= cs->L) phase = cs->L - 1;

        // y_k: symbol-center sample (fractional-delay matched filter)
        float ykI=0.f, ykQ=0.f;
        pfb_dot_phase(&cs->bank, xIQ, in_idx, phase, &ykI, &ykQ);

        // y_mid: true half-symbol earlier sample
        int half_sps = (int)floor(0.5 * cs->omega);
        size_t mid_idx = in_idx - (size_t)half_sps;

        double mu_mid = cs->mu - 0.5;
        while (mu_mid < 0.0) mu_mid += 1.0;
        int phase_mid = (int)floor(mu_mid * (double)cs->L + 0.5);
        if (phase_mid < 0) phase_mid = 0;
        if (phase_mid >= cs->L) phase_mid = cs->L - 1;

        float ymidI=0.f, ymidQ=0.f;
        pfb_dot_phase(&cs->bank, xIQ, mid_idx, phase_mid, &ymidI, &ymidQ);

        // Gardner error: Re{ ymid * conj( yk - yk_prev ) }
        float dyI = ykI - cs->yk_prev_I;
        float dyQ = ykQ - cs->yk_prev_Q;
        double e = (double)ymidI * (double)dyI + (double)ymidQ * (double)dyQ;

        // loop update
        cs->omega += cs->gain_omega * e;
        if (cs->omega < cs->omega_min) cs->omega = cs->omega_min;
        if (cs->omega > cs->omega_max) cs->omega = cs->omega_max;

        cs->mu += cs->omega + cs->gain_mu * e;

        // output sample
        outIQ[2*out_count + 0] = ykI;
        outIQ[2*out_count + 1] = ykQ;
        out_count++;

        // advance input by floor(mu)
        int adv = (int)floor(cs->mu);
        if (adv < 1) adv = 1;
        cs->mu -= (double)adv;
        in_idx += (size_t)adv;

        cs->yk_prev_I = ykI;
        cs->yk_prev_Q = ykQ;

        if (cs->verbose && (out_count % 1000 == 0)) {
            printf("[clk] out=%zu  mu=%.4f  omega=%.4f  e=%+.6g  phase=%d\n",
                   out_count, cs->mu, cs->omega, e, phase);
        }
    }
    return out_count;
}

/* ------------------------------- CLI helpers ------------------------------- */
static bool flag_present(int argc, char **argv, const char *flag){
    int i;
    for (i = 1; i < argc; i++){
        if (strcmp(argv[i], flag) == 0) return true;
    }
    return false;
}
static const char* get_opt(int argc, char **argv, const char *key, const char *defv){
    int i;
    for (i = 1; i+1 < argc; i++){
        if (strcmp(argv[i], key) == 0) return argv[i+1];
    }
    return defv;
}
static void usage(const char *p){
    fprintf(stderr,
        "Usage: %s in.cf32 out.cf32 [--sps N] [--span K] [--alpha A] [--L P]\n"
        "          [--loopbw BW] [--damp ZETA] [--orelim R]\n"
        "          [--noagc] [--notrim] [--quiet]\n", p);
}

/* ----------------------------------- main ---------------------------------- */
int main(int argc, char **argv)
{
    if (argc < 3) { usage(argv[0]); return 1; }
    const char *in_path  = argv[1];
    const char *out_path = argv[2];

    // defaults tuned for your chain
    int    SPS    = atoi(get_opt(argc, argv, "--sps",    "4"));
    int    SPAN   = atoi(get_opt(argc, argv, "--span",   "15"));
    float  ALPHA  = (float)atof(get_opt(argc, argv, "--alpha", "0.22"));
    int    L      = atoi(get_opt(argc, argv, "--L",      "16"));
    double LOOPBW = atof(get_opt(argc, argv, "--loopbw", "0.0628"));
    double DAMP   = atof(get_opt(argc, argv, "--damp",   "0.707"));
    double ORELIM = atof(get_opt(argc, argv, "--orelim", "0.375"));

    int do_agc  = flag_present(argc, argv, "--noagc")  ? 0 : 1;
    int do_trim = flag_present(argc, argv, "--notrim") ? 0 : 1;
    int quiet   = flag_present(argc, argv, "--quiet")  ? 1 : 0;

    // read cf32 input
    FILE *fi = fopen(in_path, "rb");
    if (!fi){ perror("open in"); return 1; }
    if (fseek(fi, 0, SEEK_END) != 0){ perror("fseek"); fclose(fi); return 1; }
    long nbytes = ftell(fi);
    if (nbytes < 0){ perror("ftell"); fclose(fi); return 1; }
    if (fseek(fi, 0, SEEK_SET) != 0){ perror("fseek"); fclose(fi); return 1; }

    size_t nfloats = (size_t)nbytes / sizeof(float);
    if (nfloats % 2 != 0){ fprintf(stderr,"Input float count not even (I/Q)\n"); fclose(fi); return 1; }
    size_t nin_c = nfloats / 2;

    float *xIQ = (float*)malloc(sizeof(float) * nfloats);
    if (!xIQ){ fprintf(stderr,"malloc xIQ failed\n"); fclose(fi); return 1; }
    size_t rd = fread(xIQ, sizeof(float), nfloats, fi);
    fclose(fi);
    if (rd != nfloats){ fprintf(stderr,"short read\n"); free(xIQ); return 1; }

    // build clock sync
    PFBClockSync *cs = pfbcs_create(ALPHA, SPS, SPAN, L, LOOPBW, DAMP, ORELIM, quiet ? 0 : 1);
    if (!cs){ fprintf(stderr,"pfbcs_create failed\n"); free(xIQ); return 1; }
    if (!quiet){
        printf("PFB: taps=%d -> phases=%d, taps/phase=%d\n", cs->proto_len, cs->bank.L, cs->bank.tpp);
    }

    // process
    float *yIQ = (float*)malloc(sizeof(float) * 2 * nin_c);
    if (!yIQ){ fprintf(stderr,"malloc yIQ failed\n"); pfbcs_free(cs); free(xIQ); return 1; }
    size_t nout = pfbcs_process(cs, xIQ, nin_c, yIQ);

    // Optional AGC to ~unit RMS
    if (do_agc && nout > 100){
        size_t used = nout;
        if (used > 4000) used = 4000;
        double s2 = 0.0;
        size_t i;
        for (i=0;i<used;i++){
            double I = yIQ[2*i+0], Q = yIQ[2*i+1];
            s2 += I*I + Q*Q;
        }
        double rms = 0.0;
        if (used > 0) rms = sqrt(s2 / (double)(2*used));
        if (rms > 1e-12){
            double g = 1.0 / rms;
            for (i=0;i<nout;i++){
                yIQ[2*i+0] = (float)(yIQ[2*i+0] * g);
                yIQ[2*i+1] = (float)(yIQ[2*i+1] * g);
            }
            if (!quiet) printf("AGC scale: x%.3f (RMS->1)\n", (double)g);
        }
    }

    // optional trimming of head/tail transients
    size_t start = 0, stop = nout;
    if (do_trim){
        size_t head = (size_t)(2 * SPAN);
        size_t tail = (size_t)(2 * SPAN);
        if (nout > head) start = head;
        if (nout > tail) stop  = nout - tail;
        if (stop < start){ start = 0; stop = nout; }
    }
    size_t valid = 0;
    if (stop > start) valid = stop - start;

    // write cf32 output
    FILE *fo = fopen(out_path, "wb");
    if (!fo){ perror("open out"); free(yIQ); pfbcs_free(cs); free(xIQ); return 1; }
    size_t wr = fwrite(&yIQ[2*start], sizeof(float), 2*valid, fo);
    fclose(fo);
    if (wr != 2*valid){ fprintf(stderr,"short write: %zu/%zu floats\n", wr, 2*valid); }

    if (!quiet){
        printf("Clock sync done. In=%zu complex, Out=%zu symbols, Wrote=%zu (trim=%s, agc=%s)\n",
               nin_c, nout, valid, do_trim ? "yes":"no", do_agc ? "yes":"no");
    }

    free(yIQ);
    pfbcs_free(cs);
    free(xIQ);
    return 0;
}
