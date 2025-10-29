#include "api.h"
/*
 * bits_to_float_bpsk_iq(packet_bits, num_bits, out_iq, sps)
 *  Input : bit array, number of bits, out_iq interleaved (I,Q), samples/symbol
 *  Output: number of complex samples written (num_bits * sps)
 *          I = +1 for 0, I = -1 for 1, Q = 0
 */
size_t bits_to_float_bpsk_iq(bool *packet_bits, size_t num_bits, float* out_iq, int sps){
    size_t num_complex = num_bits * (size_t)sps;
    memset(out_iq, 0, 2 * num_complex * sizeof(float));
    for(size_t i = 0; i < num_bits; i++){
        size_t n = i * (size_t)sps;
        size_t w = 2 * n;
        out_iq[w + 0] = packet_bits[i] ? -1.0f : +1.0f;
        out_iq[w + 1] = 0.0f;
    }
    return (int)num_complex;
}
