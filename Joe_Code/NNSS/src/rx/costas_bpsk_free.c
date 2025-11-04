#define _USE_MATH_DEFINES
#include "api.h"
#include "costas_bpsk_internal.h"

void   costas_bpsk_free(CostasBPSK *c){ if(c) free(c); }
double costas_bpsk_get_phase_rad(const CostasBPSK *c){ return c ? (double)c->phase : 0.0; }
double costas_bpsk_get_freq_cps (const CostasBPSK *c){ return c ? (double)c->freq_cps : 0.0; }
