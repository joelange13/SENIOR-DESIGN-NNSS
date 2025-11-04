// src/rx/fll_be_internal.h
#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct { float re, im; } fll_cplx;

static inline fll_cplx cexpj_f(float th){ fll_cplx z = { (float)cos(th), (float)sin(th) }; return z; }
static inline fll_cplx cmul_f(fll_cplx a, fll_cplx b){
    fll_cplx c; c.re = a.re*b.re - a.im*b.im; c.im = a.re*b.im + a.im*b.re; return c;
}
static inline float   cabs2_f(fll_cplx a){ return a.re*a.re + a.im*a.im; }
static inline float   sinc_f(float x){ return (fabsf(x) < 1e-8f)? 1.0f : sinf((float)M_PI*x)/((float)M_PI*x); }

/* Complex FIR */
typedef struct {
    int      Nt;
    fll_cplx *taps;
    fll_cplx *d;
} cfir_t;

static inline void cfir_init(cfir_t *f, const fll_cplx *taps, int Nt){
    f->Nt = Nt;
    f->taps = (fll_cplx*)malloc((size_t)Nt*sizeof(fll_cplx));
    f->d    = (fll_cplx*)calloc((size_t)Nt, sizeof(fll_cplx));
    memcpy(f->taps, taps, (size_t)Nt*sizeof(fll_cplx));
}
static inline fll_cplx cfir_push(cfir_t *f, fll_cplx x){
    memmove(&f->d[1], &f->d[0], (size_t)(f->Nt-1)*sizeof(fll_cplx));
    f->d[0] = x;
    float re=0.f, im=0.f;
    for(int k=0;k<f->Nt;k++){
        re += f->d[k].re*f->taps[k].re - f->d[k].im*f->taps[k].im;
        im += f->d[k].re*f->taps[k].im + f->d[k].im*f->taps[k].re;
    }
    fll_cplx y = { re, im };
    return y;
}
static inline void cfir_free(cfir_t *f){
    if (f->taps) free(f->taps);
    if (f->d)    free(f->d);
    f->taps=NULL; f->d=NULL; f->Nt=0;
}

/* Forward-declared public handle lives here so layout is shared by all .c files */
struct FLLBandEdge {
    /* Config */
    int    sps;
    int    span;
    float  alpha;
    float  be_frac;
    int    Ntaps;
    int    decim;
    int    warm;
    float  cps_clamp;

    /* Loop gains / smoothing */
    float  Kp, Ki;
    float  beta;

    /* NCO + powers */
    float  freq_cps;
    float  phase;
    float  Pu, Pl;

    /* Band-edge FIRs */
    cfir_t   fu, fl;
    fll_cplx *hu, *hl;

    /* Accumulators */
    int      decim_cnt;
    double   e_accum;

    /* Verbose once */
    int verbose;
};
