// include/pfb_internal.h
#ifndef PFB_INTERNAL_H
#define PFB_INTERNAL_H

#include "api.h"

// Polyphase bank (phase-major layout)
typedef struct {
    int   L;      // number of phases
    int   tpp;    // taps per phase (odd)
    float *b;     // flattened [L][tpp], phase-major: b[p*tpp + k]
} PolyBank;

// Opaque in api.h, concrete here for all rx*.c
struct PFBClockSync {
    // loop state
    double omega, omega_min, omega_max;
    double mu;
    double gain_mu, gain_omega;

    // bank/taps
    int    sps_nom;
    int    span_syms;
    int    L;
    int    proto_len;
    float *proto;
    PolyBank bank;

    // Gardner memory
    float  yk_prev_I, yk_prev_Q;

    // inline AGC bookkeeping (even if disabled)
    double agc_gain, agc_ref, agc_rate;
    int    agc_enabled;

    // options
    int    verbose;
};

#endif // PFB_INTERNAL_H
