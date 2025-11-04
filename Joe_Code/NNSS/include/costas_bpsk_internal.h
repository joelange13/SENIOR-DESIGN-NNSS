#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct { float re, im; } cb_cplx;
static inline cb_cplx cmul(cb_cplx a, cb_cplx b){
    cb_cplx c; c.re = a.re*b.re - a.im*b.im; c.im = a.re*b.im + a.im*b.re; return c;
}
static inline cb_cplx cexpj_f(float th){ cb_cplx z = { (float)cos(th), (float)sin(th) }; return z; }

struct CostasBPSK {
    /* loop */
    float Kp, Ki;
    float freq_cps;   /* cycles/sample (symbol-rate domain) */
    float phase;      /* radians */

    /* simple AGC (RMS tracker) */
    float p_hat;      /* power estimate */
    float agc_beta;   /* e.g., 1e-3 */
    float gain;       /* 1/sqrt(p_hat) style */

    int   verbose;
};
