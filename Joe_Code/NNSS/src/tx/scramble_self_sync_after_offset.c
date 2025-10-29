#include "api.h"
/*
 * scramble_self_sync_after_offset(bits, n_bits, offset)
 *  Input : bits[] in-place scrambled after 'offset' (preamble preserved)
 *  Output: bits[] modified using s[k] = d[k] XOR s[k-3] XOR s[k-5]
 */
void scramble_self_sync_after_offset(bool *bits, size_t n_bits, size_t offset){
    bool s1=0, s2=0, s3=0, s4=0, s5=0;
    for (size_t k = 0; k < n_bits; ++k) {
        if (k < offset) continue;
        bool d   = bits[k];
        bool s_k = d ^ s3 ^ s5;
        bits[k]  = s_k;
        s5 = s4; s4 = s3; s3 = s2; s2 = s1; s1 = s_k;
    }
}
