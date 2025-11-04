#define _USE_MATH_DEFINES
#include "api.h"
#include "costas_bpsk_internal.h"

CostasBPSK* costas_bpsk_create(double loop_bw, double damping, double agc_beta, int verbose){
    (void)damping; /* reserved for PLL design mapping; we use a conservative heuristic */
    CostasBPSK *c = (CostasBPSK*)calloc(1, sizeof(CostasBPSK));
    if(!c) return NULL;

    /* Heuristic: small bandwidth at 1 sps; tweak if you want a flag later */
    float bw = (float)((loop_bw > 0.0) ? loop_bw : 0.01f);
    c->Kp = 2.0f * bw;          /* proportional */
    c->Ki = 0.25f * bw * bw;    /* integral; gentle so it doesn’t wander */
    if (c->Kp < 1e-4f) c->Kp = 1e-4f;
    if (c->Ki < 1e-6f) c->Ki = 1e-6f;

    c->freq_cps = 0.0f;
    c->phase    = 0.0f;

    c->agc_beta = (float)((agc_beta > 0.0) ? agc_beta : 1e-3);
    c->p_hat    = 1.0f;   /* start sane */
    c->gain     = 1.0f;

    c->verbose  = verbose;
    if (c->verbose){
        fprintf(stderr, "Costas(BPSK): Kp=%.6f Ki=%.6f  bw=%.4g  agc_beta=%.4g\n",
                c->Kp, c->Ki, bw, c->agc_beta);
    }
    return c;
}
