#include "api.h"
/*
 * agc_free(a)
 *  Input : AGC pointer (nullable)
 *  Output: frees memory
 */
void agc_free(StreamAGC *a){ if (a) free(a); }
