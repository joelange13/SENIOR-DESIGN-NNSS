#include "api.h"
/*
 * write_cf32(file, shaped_out_iq, shaped_out_iq_length)
 *  Input : file (binary), interleaved floats (I,Q) length = shaped_out_iq_length
 *  Output: writes to file; prints a confirmation line
 */
void write_cf32(FILE *file, float *shaped_out_iq, size_t shaped_out_iq_length){
    (void)fwrite(shaped_out_iq, sizeof(float), shaped_out_iq_length, file);
    printf("Wrote cf32 IQ to file.\n");
}
