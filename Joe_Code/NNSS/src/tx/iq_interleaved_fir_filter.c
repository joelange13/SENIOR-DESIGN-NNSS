#include "api.h"
/*
 * iq_interleaved_fir_filter(x, in_floats, h, num_taps, y)
 *  Input : x = interleaved IQ (len in_floats = 2*Nin), real taps h[num_taps]
 *  Output: y = interleaved IQ, length 2*(Nin + num_taps - 1)
 */
void iq_interleaved_fir_filter(float *x, size_t in_floats, float *h, int num_taps, float *y){
    size_t nin  = in_floats / 2;
    size_t nout = nin + (size_t)num_taps - 1;
    memset(y, 0, sizeof(float) * 2 * nout);

    for(size_t n = 0; n < nout; n++){
        double accI = 0.0, accQ = 0.0;
        for(int k = 0; k < num_taps; k++){
            long idx = (long)n - k;
            if(idx >= 0 && (size_t)idx < nin){
                size_t ii = (size_t)idx * 2;
                accI += (double)h[k]*(double)x[ii + 0];
                accQ += (double)h[k]*(double)x[ii + 1];
            }
        }
        y[2*n + 0] = (float)accI;
        y[2*n + 1] = (float)accQ;
    }
}
