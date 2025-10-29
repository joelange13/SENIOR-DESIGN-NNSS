#include "api.h"
/*
 * agc_bootstrap(a, xIQ, nin_complex, sample_cap)
 *  Input : AGC*, interleaved IQ, sample count, max samples to average
 *  Output: adjusts initial gain/power estimate
 */
void agc_bootstrap(StreamAGC *a, const float *xIQ, size_t nin_complex, size_t sample_cap){
    if (!a || !xIQ || nin_complex == 0) return;
    size_t N = nin_complex < sample_cap ? nin_complex : sample_cap;
    double sum = 0.0; size_t cnt = 0;
    for (size_t i = 0; i < N; i++){
        double I = xIQ[2*i+0], Q = xIQ[2*i+1];
        double p = I*I + Q*Q;
        if (p > 0) { sum += p; cnt++; }
    }
    if (cnt > 0){
        double p_avg = sum / (double)cnt;
        if (p_avg > 1e-12){
            double g = sqrt(a->target_power / p_avg);
            if (g < a->min_gain) g = a->min_gain;
            if (g > a->max_gain) g = a->max_gain;
            a->gain = g;
            a->p_hat = p_avg;
        }
    }
}
