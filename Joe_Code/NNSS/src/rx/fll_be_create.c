// src/rx/fll_be_create.c
#include "api.h"
#include "fll_be_internal.h"

/* Hamming-windowed LP prototype (DC normalized) */
static void lp_sinc_hamming(float fc_cyc_per_samp, int Ntaps, float *h){
    int M = (Ntaps-1)/2;
    for(int n=0; n<Ntaps; n++){
        float m = (float)(n - M);
        float w = 0.54f - 0.46f*(float)cos(2.0*M_PI*(double)n/(double)(Ntaps-1));
        h[n] = 2.f*fc_cyc_per_samp * sinc_f(2.f*fc_cyc_per_samp*m) * w;
    }
    float s=0.f; for(int i=0;i<Ntaps;i++) s += h[i];
    if (s != 0.f) for(int i=0;i<Ntaps;i++) h[i] /= s;
}

/* Design complex band-edge filters at +/- f_edge */
static void design_bandedges(int sps, float alpha, float be_frac, int Ntaps, fll_cplx *hu, fll_cplx *hl){
    float df = alpha/(2.f*sps);
    float *hlp = (float*)calloc((size_t)Ntaps, sizeof(float));
    lp_sinc_hamming(df, Ntaps, hlp);

    float f_edge = (1.f/(2.f*sps)) * (be_frac > 0.f ? be_frac : 1.f);
    int M = (Ntaps-1)/2;
    for(int n=0; n<Ntaps; n++){
        float th = 2.f*(float)M_PI*f_edge*(float)(n - M);
        float c = cosf(th), s = sinf(th);
        hu[n].re = hlp[n]*c;  hu[n].im = hlp[n]*s;   /* +f_edge */
        hl[n].re = hlp[n]*c;  hl[n].im =-hlp[n]*s;   /* -f_edge */
    }
    free(hlp);
}

FLLBandEdge* fll_be_create(float alpha, int sps, int span,
                           float be_frac, double loop_bw, double damping,
                           double pwr_beta, int verbose)
{
    (void)damping;
    FLLBandEdge *f = (FLLBandEdge*)calloc(1, sizeof(FLLBandEdge));
    if(!f) return NULL;

    f->alpha   = alpha;
    f->sps     = sps;
    f->span    = span;
    f->be_frac = (be_frac <= 0.f) ? 1.0f : be_frac;
    f->Ntaps   = span * sps + 1;
    if (f->Ntaps < 401)  f->Ntaps = 401;     /* ensure robust averaging */
    if (f->Ntaps > 1201) f->Ntaps = 1201;

    /* Heuristic mapping from loop_bw to Kp/Ki (conservative defaults) */
    float bw = (float)((loop_bw > 0.0) ? loop_bw : 0.004);
    f->Kp = 0.004f;
    f->Ki = 0.00004f;
    if (f->Kp < 0.001f) f->Kp = 0.001f;
    if (f->Ki < 1e-6f)  f->Ki = 1e-6f;

    f->beta = (float)((pwr_beta > 0.0) ? pwr_beta : 0.0015);

    f->decim = 8;
    f->warm  = f->Ntaps + 12*f->sps;
    f->cps_clamp = 0.05f;

    f->freq_cps = 0.f;
    f->phase    = 0.f;
    f->Pu = 0.f; f->Pl = 0.f;
    f->decim_cnt = 0; f->e_accum = 0.0;
    f->verbose = verbose;

    /* Allocate and design taps */
    f->hu = (fll_cplx*)malloc((size_t)f->Ntaps*sizeof(fll_cplx));
    f->hl = (fll_cplx*)malloc((size_t)f->Ntaps*sizeof(fll_cplx));
    if(!f->hu || !f->hl){ free(f->hu); free(f->hl); free(f); return NULL; }
    design_bandedges(f->sps, f->alpha, f->be_frac, f->Ntaps, f->hu, f->hl);

    /* Init FIRs */
    cfir_init(&f->fu, f->hu, f->Ntaps);
    cfir_init(&f->fl, f->hl, f->Ntaps);

    if (f->verbose){
        fprintf(stderr,
            "FLL: SPS=%d span=%d Ntaps=%d alpha=%.3f be_frac=%.3f  Kp=%.5f Ki=%.6f beta=%.5f decim=%d warm=%d clamp=%.3f\n",
            f->sps, f->span, f->Ntaps, f->alpha, f->be_frac, f->Kp, f->Ki, f->beta, f->decim, f->warm, f->cps_clamp);
    }
    return f;
}
