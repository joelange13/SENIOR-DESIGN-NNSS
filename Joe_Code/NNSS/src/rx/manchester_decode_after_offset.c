#include "api.h"
/*
 * manchester_decode_after_offset(in, n, offset, &out)
 *  Input : in-bits with Manchester code after offset
 *  Output: allocates out-bits decoded; returns length or 0 on error
 */
size_t manchester_decode_after_offset(const bool *in, size_t n, size_t offset, bool **outp){
    if ((n < offset) || ((n-offset) & 1)) return 0;
    size_t out_n = offset + (n - offset)/2;
    bool *out = (bool*)malloc(out_n * sizeof(bool));
    if (!out) return 0;
    memcpy(out, in, offset * sizeof(bool));
    size_t w = offset;
    for (size_t i = offset; i < n; i += 2) {
        bool a = in[i], b = in[i+1];
        if (a==0 && b==1) out[w++] = 0;
        else if (a==1 && b==0) out[w++] = 1;
        else { free(out); return 0; }
    }
    *outp = out;
    return out_n;
}
