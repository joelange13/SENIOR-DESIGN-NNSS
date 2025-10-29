#include "api.h"
/*
 * descramble_self_sync_after_offset(bits, n, offset)
 *  Input : bit array scrambled with s[k]=d[k]^d[k-3]^d[k-5], skip first 'offset'
 *  Output: in-place descramble to original bits after offset
 */
void descramble_self_sync_after_offset(bool *s, size_t n, size_t offset){
    bool d1=0, d2=0, d3=0, d4=0, d5=0;
    for (size_t k = 0; k < n; ++k) {
        if (k < offset) continue;
        bool s_k = s[k];
        bool d_k = s_k ^ d3 ^ d5;
        s[k] = d_k;
        d5=d4; d4=d3; d3=d2; d2=d1; d1=d_k;
    }
}
