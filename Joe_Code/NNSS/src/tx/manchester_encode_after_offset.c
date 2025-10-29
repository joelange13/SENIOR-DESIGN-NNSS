#include "api.h"
/*
 * manchester_encode_after_offset(in_bits, n_bits, offset, &out_bits)
 *  Input : in_bits[], total n_bits, do not encode first 'offset' bits
 *  Output: allocates *out_bits with length offset + 2*(n_bits-offset); returns length
 */
size_t manchester_encode_after_offset(const bool *in_bits, size_t n_bits,
                                      size_t offset, bool **out_bits){
    size_t out_n = offset + 2*(n_bits - offset);
    bool *out = (bool*)malloc(out_n * sizeof(bool));
    if (!out) return 0;

    memcpy(out, in_bits, offset * sizeof(bool));
    size_t w = offset;
    for (size_t k = offset; k < n_bits; ++k) {
        bool b = in_bits[k];
        if (!b) { out[w++] = 0; out[w++] = 1; }
        else    { out[w++] = 1; out[w++] = 0; }
    }
    *out_bits = out;
    return out_n;
}
