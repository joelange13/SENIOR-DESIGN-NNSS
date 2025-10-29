#include "api.h"
/*
 * agc_apply_inplace(a, xIQ, nin_complex)
 *  Input : AGC*, interleaved IQ, sample count
 *  Output: scales xIQ in-place to target power
 */
void agc_apply_inplace(StreamAGC *a, float *xIQ, size_t nin_complex){
    if (!a || !xIQ || nin_complex == 0) return;
    const double eps = 1e-12;
    for (size_t i = 0; i < nin_complex; i++){
        double I = xIQ[2*i+0];
        double Q = xIQ[2*i+1];
        double p = I*I + Q*Q;
        a->p_hat = (1.0 - a->beta)*a->p_hat + a->beta*p;
        double g_target = sqrt(a->target_power / (a->p_hat + eps));
        a->gain = (1.0 - a->gamma)*a->gain + a->gamma*g_target;
        if (a->gain < a->min_gain) a->gain = a->min_gain;
        if (a->gain > a->max_gain) a->gain = a->max_gain;
        xIQ[2*i+0] = (float)(a->gain * I);
        xIQ[2*i+1] = (float)(a->gain * Q);
    }
}
