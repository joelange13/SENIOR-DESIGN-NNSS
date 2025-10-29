#include "api.h"
#include "pfb_internal.h"



/* Local dot-product helper for a given phase */
static inline void pfb_dot_phase(const PolyBank *pb, const float *xIQ, size_t center_idx,
                                 int phase, float *outI, float *outQ)
{
    const float *taps = &pb->b[phase * pb->tpp];
    int half = (pb->tpp - 1)/2;
    double accI = 0.0, accQ = 0.0;
    for (int k = -half; k <= half; k++){
        long ii = (long)center_idx + k;
        if (ii < 0) continue;
        const float *s = &xIQ[(size_t)ii * 2];
        double hk = (double)taps[k + half];
        accI += hk * (double)s[0];
        accQ += hk * (double)s[1];
    }
    *outI = (float)accI;
    *outQ = (float)accQ;
}

/*
 * pfbcs_process(cs, xIQ, nin_complex, outIQ)
 *  Input : PFBClockSync*, input interleaved IQ (nin_complex)
 *  Output: interleaved symbol-rate matched-filter outputs; returns #symbols
 */
size_t pfbcs_process(PFBClockSync *cs, const float *xIQ, size_t nin_complex, float *outIQ){
    const int tpp  = cs->bank.tpp;
    const int half = (tpp - 1)/2;

    size_t in_idx  = (size_t)(half + cs->sps_nom + 2);
    size_t end_idx = nin_complex - (size_t)half - (size_t)cs->sps_nom - 2;
    size_t out_count = 0;

    double env_hat = 1.0, env_beta = 1e-3;

    while (in_idx < end_idx) {
        double pf  = cs->mu * (double)cs->L;
        int phase  = (int)(pf);
        if (phase < 0) phase = 0;
        if (phase >= cs->L) phase = cs->L - 1;

        float ykI=0.f, ykQ=0.f;
        pfb_dot_phase(&cs->bank, xIQ, in_idx, phase, &ykI, &ykQ);

        double half_omega = cs->omega * 0.5;
        int half_sps      = (int)floor(half_omega);
        size_t mid_idx    = in_idx - (size_t)half_sps;

        double mu_mid = cs->mu + 0.5;
        if (mu_mid >= 1.0) mu_mid -= 1.0;

        double pf_mid = mu_mid * (double)cs->L;
        int phase_mid = (int)(pf_mid);
        if (phase_mid < 0) phase_mid = 0;
        if (phase_mid >= cs->L) phase_mid = cs->L - 1;

        float ymidI=0.f, ymidQ=0.f;
        pfb_dot_phase(&cs->bank, xIQ, mid_idx, phase_mid, &ymidI, &ymidQ);

        float dyI = ykI - cs->yk_prev_I;
        float dyQ = ykQ - cs->yk_prev_Q;
        double e  = (double)ymidI * (double)dyI + (double)ymidQ * (double)dyQ;

        double p  = (double)ykI*ykI + (double)ykQ*ykQ;
        env_hat = (1.0 - env_beta)*env_hat + env_beta * (p > 1e-12 ? p : env_hat);
        double e_norm = e / (env_hat + 1e-12);

        cs->mu    += cs->gain_mu    * e_norm;
        cs->omega += cs->gain_omega * e_norm;

        if (cs->omega < cs->omega_min) cs->omega = cs->omega_min;
        if (cs->omega > cs->omega_max) cs->omega = cs->omega_max;

        outIQ[2*out_count + 0] = ykI;
        outIQ[2*out_count + 1] = ykQ;
        out_count++;

        cs->mu += cs->omega;
        int adv = (int)floor(cs->mu);
        if (adv < 1) adv = 1;
        cs->mu -= (double)adv;
        in_idx += (size_t)adv;

        cs->yk_prev_I = ykI;
        cs->yk_prev_Q = ykQ;

        if (cs->verbose && (out_count % 1000 == 0)) {
            printf("[clk] out=%zu  mu=%.4f  omega=%.4f  e=%+.6g\n", out_count, cs->mu, cs->omega, e_norm);
        }
    }
    return out_count;
}
