// src/rx/fll_be_process.c

#include "api.h"
#include "fll_be_internal.h"

size_t fll_be_process(FLLBandEdge *S,
                      const float *xIQ, size_t nin_complex,
                      float *yIQ)
{
    int decim_cnt = S->decim_cnt;
    double e_acc  = S->e_accum;

    size_t nout = 0;
    for(size_t n=0; n<nin_complex; n++){
        /* Mix by current NCO (derotate) */
        fll_cplx x = { xIQ[2*n+0], xIQ[2*n+1] };
        fll_cplx ejm = cexpj_f(-S->phase);
        fll_cplx y   = cmul_f(x, ejm);

        /* Edge filters on corrected sample */
        fll_cplx yu = cfir_push(&S->fu, y);
        fll_cplx yl = cfir_push(&S->fl, y);

        /* Smoothed powers */
        float pu = cabs2_f(yu), pl = cabs2_f(yl);
        S->Pu = (1.f - S->beta)*S->Pu + S->beta*pu;
        S->Pl = (1.f - S->beta)*S->Pl + S->beta*pl;

        /* Normalized error (Pu - Pl) to match your probe polarity */
        float num = S->Pu - S->Pl;
        float den = S->Pu + S->Pl + 1e-12f;
        float e   = num / den;

        if ((int)n < S->warm){
            /* No PI updates during warmup */
            S->phase += 2.f*(float)M_PI * S->freq_cps;
        } else {
            /* Decimated PI updates */
            e_acc += (double)e;
            decim_cnt++;
            if (decim_cnt == S->decim){
                float e_avg = (float)(e_acc / (double)S->decim);
                e_acc = 0.0; decim_cnt = 0;

                if (e_avg >  0.5f) e_avg =  0.5f;
                if (e_avg < -0.5f) e_avg = -0.5f;

                S->freq_cps += S->Ki * e_avg;                    /* I */
                if (S->freq_cps >  S->cps_clamp) S->freq_cps =  S->cps_clamp;
                if (S->freq_cps < -S->cps_clamp) S->freq_cps = -S->cps_clamp;

                S->phase    += 2.f*(float)M_PI * S->freq_cps      /* advance */
                               + S->Kp * e_avg;                   /* P */
            } else {
                S->phase += 2.f*(float)M_PI * S->freq_cps;
            }
        }

        /* Write corrected sample */
        yIQ[2*n+0] = y.re;
        yIQ[2*n+1] = y.im;
        nout++;
    }

    S->decim_cnt = decim_cnt;
    S->e_accum   = e_acc;
    return nout;
}
