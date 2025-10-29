#include "api.h"
/*
 * rrc_generate_taps(alpha, sps, span, taps_out)
 *  Input : alpha (0..1), sps (>=2), span (symbols), taps_out length = span*sps+1
 *  Output: taps_out filled with unit-energy RRC impulse response
 */
void rrc_generate_taps(float alpha, int sps, int span, float *taps_out){
    const int num_taps = span * sps + 1;
    const int mid      = num_taps / 2;
    const double a     = (double)alpha;
    const double T     = 1.0;
    const double Ts    = T / (double)sps;
    const double eps   = 1e-12;

    for (int n = 0; n < num_taps; n++) {
        const int k   = n - mid;
        const double t = k * Ts;
        double h;
        if (fabs(t) < eps) {
            h = 1.0 + a * (4.0/M_PI - 1.0);
        } else if (a > 0.0 && fabs(fabs(t) - T/(4.0*a)) < 1e-9) {
            const double theta = M_PI/(4.0*a);
            h = (a / sqrt(2.0)) * ((1.0 + 2.0/M_PI)*sin(theta) + (1.0 - 2.0/M_PI)*cos(theta));
        } else {
            const double num = sin(M_PI * t * (1.0 - a) / T)
                             + (4.0 * a * t / T) * cos(M_PI * t * (1.0 + a) / T);
            const double den = (M_PI * t / T) * (1.0 - pow(4.0 * a * t / T, 2.0));
            h = num / den;
        }
        taps_out[n] = (float)h;
    }
    /* Normalize to unit energy */
    double e = 0.0;
    for (int n = 0; n < num_taps; n++) e += (double)taps_out[n] * (double)taps_out[n];
    const double scale = (e > 0.0) ? 1.0 / sqrt(e) : 1.0;
    for (int n = 0; n < num_taps; n++) taps_out[n] = (float)(taps_out[n] * scale);
}
