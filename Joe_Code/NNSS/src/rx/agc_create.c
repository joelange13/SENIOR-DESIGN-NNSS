#include "api.h"
/*
 * agc_create(beta, gamma, target_rms, init_gain)
 *  Input : IIR betas, desired RMS, initial gain
 *  Output: heap-allocated StreamAGC* (or NULL)
 */
StreamAGC* agc_create(double beta, double gamma, double target_rms, double init_gain){
    StreamAGC *a = (StreamAGC*)calloc(1, sizeof(StreamAGC));
    if (!a) return NULL;
    a->beta = (beta > 0.0 && beta < 1.0) ? beta : 1e-3;
    a->gamma = (gamma > 0.0 && gamma < 1.0) ? gamma : 1e-3;
    a->target_power = (target_rms > 0.0 ? target_rms : 1.0);
    a->target_power *= a->target_power;
    a->p_hat = a->target_power;
    a->gain  = (init_gain > 0.0 ? init_gain : 1.0);
    a->min_gain = 0.05;
    a->max_gain = 20.0;
    return a;
}
