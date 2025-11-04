#define _USE_MATH_DEFINES
#include "api.h"
#include "costas_bpsk_internal.h"

size_t costas_bpsk_process(CostasBPSK *c,
                           const float *xIQ, size_t nin_complex,
                           float *yIQ)
{
    for (size_t n = 0; n < nin_complex; n++){
        /* input at ~1 sps (post-PFB) */
        cb_cplx x = { xIQ[2*n+0], xIQ[2*n+1] };

        /* derotate by current NCO */
        cb_cplx ejm = cexpj_f(-c->phase);
        cb_cplx y   = cmul(x, ejm);

        /* tiny internal AGC for stability (decision-directed prefers ~unit power) */
        float p = y.re*y.re + y.im*y.im;
        c->p_hat = (1.0f - c->agc_beta)*c->p_hat + c->agc_beta*p;
        c->gain  = (c->p_hat > 1e-12f) ? (1.0f / sqrtf(c->p_hat)) : 1.0f;
        y.re *= c->gain; y.im *= c->gain;

        /* BPSK Costas DD phase detector:
           For BPSK, a simple & robust choice is e = I * Q (after derotation + AGC).
           (When constellation is near ±1 on I and Q≈0, sign(I)*Q ≈ I*Q anyway.) */
        float e = y.re * y.im;

        /* 2nd-order PI update in symbol-rate domain */
        c->freq_cps += c->Ki * e;                    /* integral */
        /* optional clamp if you want: e.g., +/-0.05 cps */
        if (c->freq_cps >  0.05f) c->freq_cps =  0.05f;
        if (c->freq_cps < -0.05f) c->freq_cps = -0.05f;
        c->phase     += 2.0f*(float)M_PI * c->freq_cps
                        + c->Kp * e;

        /* write corrected sample */
        yIQ[2*n+0] = y.re;
        yIQ[2*n+1] = y.im;
    }
    /* keep phase bounded (optional) */
    if (c->phase >  1e6f) c->phase -= 1e6f;
    if (c->phase < -1e6f) c->phase += 1e6f;

    return nin_complex;
}
