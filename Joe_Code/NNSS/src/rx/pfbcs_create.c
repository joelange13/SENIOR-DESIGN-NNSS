#include "api.h"
#include "pfb_internal.h"



/* Local helpers (private) */
static PolyBank pfb_from_proto(const float *proto, int proto_len, int L){
    PolyBank pb = (PolyBank){0,0,NULL};
    if (L <= 0 || proto_len <= 0) return pb;
    int padded = proto_len;
    if (padded % L != 0) padded = proto_len + (L - (proto_len % L));
    int tpp = padded / L;
    if ((tpp % 2) == 0) { padded += L; tpp += 1; }
    float *buf = (float*)calloc((size_t)L * tpp, sizeof(float));
    if (!buf) return pb;
    for (int p = 0; p < L; p++){
        for (int k = 0; k < tpp; k++){
            int src = p + k*L;
            float v = 0.0f;
            if (src < proto_len) v = proto[src];
            buf[p*tpp + k] = v;
        }
    }
    pb.L = L; pb.tpp = tpp; pb.b = buf;
    return pb;
}




void pfbcs_free(PFBClockSync *cs); /* forward so we can reuse helpers */

/*
 * pfbcs_create(alpha, sps_nom, span_syms, L, loop_bw, damping, omega_rel_lim, enable_agc, agc_rate, verbose)
 *  Output: allocated PFBClockSync* or NULL
 */
PFBClockSync* pfbcs_create(float alpha, int sps_nom, int span_syms, int L,
                           double loop_bw, double damping, double omega_rel_lim,
                           int enable_agc, double agc_rate, int verbose)
{
    if (sps_nom <= 0 || span_syms <= 0 || L <= 0) return NULL;

    PFBClockSync *cs = (PFBClockSync*)calloc(1, sizeof(PFBClockSync));
    if (!cs) return NULL;

    cs->sps_nom = sps_nom; cs->span_syms = span_syms; cs->L = L; cs->verbose = verbose;

    int N = span_syms * sps_nom + 1;
    cs->proto = (float*)malloc(sizeof(float) * (size_t)N);
    if (!cs->proto){ free(cs); return NULL; }
    rrc_generate_taps(alpha, sps_nom, span_syms, cs->proto);
    cs->proto_len = N;

    cs->bank = pfb_from_proto(cs->proto, cs->proto_len, L);

    double denom = 1.0 + 2.0*damping*loop_bw + loop_bw*loop_bw;
    double Kp = 1.0, Ki = 2.0;
    cs->gain_mu    = (4.0*damping*loop_bw) / (denom * Kp);
    cs->gain_omega = (4.0*loop_bw*loop_bw) / (denom * Ki);

    cs->omega = (double)sps_nom;
    cs->omega_min = cs->omega * (1.0 - omega_rel_lim);
    cs->omega_max = cs->omega * (1.0 + omega_rel_lim);
    cs->mu = 0.5;

    cs->agc_enabled = enable_agc;
    cs->agc_gain = 1.0;
    cs->agc_ref  = 1.0;
    cs->agc_rate = agc_rate;

    if (verbose){
        printf("PFB: taps=%d -> phases=%d, taps/phase=%d\n", cs->proto_len, cs->bank.L, cs->bank.tpp);
        printf("Loop: gain_mu=%.6f, gain_omega=%.6f\n", cs->gain_mu, cs->gain_omega);
    }
    return cs;
}
