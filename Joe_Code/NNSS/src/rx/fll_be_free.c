// src/rx/fll_be_free.c

#include "api.h"
#include "fll_be_internal.h"

void fll_be_free(FLLBandEdge *S){
    if(!S) return;
    cfir_free(&S->fu);
    cfir_free(&S->fl);
    if (S->hu) free(S->hu);
    if (S->hl) free(S->hl);
    free(S);
}

double fll_be_get_freq_cps(const FLLBandEdge *S){ return (double)S->freq_cps; }
double fll_be_get_freq_rad(const FLLBandEdge *S){ return 2.0*M_PI*(double)S->freq_cps; }
