#include "api.h"
/*
 * tx_pulse_shape(out_iq, nsamps_unshaped, sps, alpha, span, &shaped_out_iq, &len)
 *  Input : unshaped interleaved IQ (len nsamps_unshaped), RRC params
 *  Output: allocates *shaped_out_iq (caller must free), sets *len to float count
 */
void tx_pulse_shape(float *out_iq, size_t nsamps_unshaped, int sps, float alpha, int span,
                    float **shaped_out_iq, size_t *nsamps_shaped_out){
    int num_taps = span * sps + 1;
    float *taps = (float*)malloc(sizeof(float) * num_taps);
    rrc_generate_taps(alpha, sps, span, taps);

    size_t nin_complex  = nsamps_unshaped / 2;
    size_t nout_complex = nin_complex + (size_t)num_taps - 1;
    float *y = (float*)malloc(sizeof(float) * 2 * nout_complex);

    iq_interleaved_fir_filter(out_iq, nsamps_unshaped, taps, num_taps, y);

    *shaped_out_iq       = y;
    *nsamps_shaped_out   = 2 * nout_complex;

    free(taps);
}
