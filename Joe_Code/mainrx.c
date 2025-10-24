#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

/* ----------------------  RRC generator (as given) ---------------------- */
static void rrc_generate_taps(float alpha, int sps, int span, float *taps_out)
{
    const int num_taps = span * sps + 1;     // odd length
    const int mid      = num_taps / 2;
    const double a     = (double)alpha;
    const double T     = 1.0;                 // symbol period
    const double Ts    = T / (double)sps;
    const double eps   = 1e-12;

    for (int n = 0; n < num_taps; n++) {
        const int k   = n - mid;
        const double t = k * Ts;              // time in symbol units
        double h;

        // t = 0 → closed form
        if (fabs(t) < eps) {
            h = 1.0 + a * (4.0/M_PI - 1.0);
        }
        // t = ±T/(4α) → closed form
        else if (a > 0.0 && fabs(fabs(t) - T/(4.0*a)) < 1e-9) {
            const double theta = M_PI/(4.0*a);
            h = (a / sqrt(2.0)) * ((1.0 + 2.0/M_PI)*sin(theta) + (1.0 - 2.0/M_PI)*cos(theta));
        }
        // general case
        else {
            const double num = sin(M_PI * t * (1.0 - a) / T)
                             + (4.0 * a * t / T) * cos(M_PI * t * (1.0 + a) / T);
            const double den = (M_PI * t / T) * (1.0 - pow(4.0 * a * t / T, 2.0));
            h = num / den;
        }

        taps_out[n] = (float)h;
    }

    // Normalize for unity energy (so TX+RX ≈ RC with unity peak)
    double e = 0.0;
    for (int n = 0; n < num_taps; n++) e += (double)taps_out[n] * (double)taps_out[n];
    const double scale = (e > 0.0) ? 1.0 / sqrt(e) : 1.0;
    for (int n = 0; n < num_taps; n++) taps_out[n] = (float)(taps_out[n] * scale);
}

/* ===================== Stream-level AGC (pre-PFB) ===================== */

typedef struct {
    double p_hat;         // smoothed power estimate (I^2+Q^2)
    double gain;          // current linear gain
    double beta;          // power IIR smoothing (1e-3 .. 1e-2)
    double gamma;         // gain smoothing (1e-3 .. 5e-3)
    double target_power;  // desired output power after AGC (e.g., 1.0)
    double min_gain, max_gain;
} StreamAGC;

static StreamAGC* agc_create(double beta, double gamma,
                             double target_rms, double init_gain)
{
    StreamAGC *a = (StreamAGC*)calloc(1, sizeof(StreamAGC));
    if (!a) return NULL;
    a->beta = (beta > 0.0 && beta < 1.0) ? beta : 1e-3;
    a->gamma = (gamma > 0.0 && gamma < 1.0) ? gamma : 1e-3;
    a->target_power = (target_rms > 0.0 ? target_rms : 1.0);
    a->target_power *= a->target_power;     // store as power
    a->p_hat = a->target_power;             // start reasonable
    a->gain  = (init_gain > 0.0 ? init_gain : 1.0);
    a->min_gain = 0.05;
    a->max_gain = 20.0;
    return a;
}

static void agc_free(StreamAGC *a){
    if (a) free(a);
}

// Quick one-time bootstrap so we start near the right gain.
static void agc_bootstrap(StreamAGC *a, const float *xIQ,
                          size_t nin_complex, size_t sample_cap)
{
    if (!a || !xIQ || nin_complex == 0) return;
    size_t N = nin_complex < sample_cap ? nin_complex : sample_cap;
    double sum = 0.0;
    size_t cnt = 0;
    for (size_t i = 0; i < N; i++){
        double I = xIQ[2*i+0], Q = xIQ[2*i+1];
        double p = I*I + Q*Q;
        if (p > 0) { sum += p; cnt++; }
    }
    if (cnt > 0){
        double p_avg = sum / (double)cnt;
        if (p_avg > 1e-12){
            double g = sqrt(a->target_power / p_avg);
            if (g < a->min_gain) g = a->min_gain;
            if (g > a->max_gain) g = a->max_gain;
            a->gain = g;
            a->p_hat = p_avg;
        }
    }
}

// Apply AGC in place to the raw input stream (I/Q interleaved).
static void agc_apply_inplace(StreamAGC *a, float *xIQ, size_t nin_complex)
{
    if (!a || !xIQ || nin_complex == 0) return;
    const double eps = 1e-12;
    for (size_t i = 0; i < nin_complex; i++){
        double I = xIQ[2*i+0];
        double Q = xIQ[2*i+1];

        // update power IIR
        double p = I*I + Q*Q;
        a->p_hat = (1.0 - a->beta)*a->p_hat + a->beta*p;

        // target gain from smoothed power
        double g_target = sqrt(a->target_power / (a->p_hat + eps));

        // smooth the gain
        a->gain = (1.0 - a->gamma)*a->gain + a->gamma*g_target;
        if (a->gain < a->min_gain) a->gain = a->min_gain;
        if (a->gain > a->max_gain) a->gain = a->max_gain;

        // scale sample
        xIQ[2*i+0] = (float)(a->gain * I);
        xIQ[2*i+1] = (float)(a->gain * Q);
    }
}


/* ------------------------- polyphase bank utilities ------------------------ */
typedef struct {
    int L;      // number of phases
    int tpp;    // taps per phase
    float *b;   // flattened [L][tpp], phase-major
} PolyBank;

// Build bank by padding proto to multiple of L, then de-interleaving: h_p[k] = h[p + k*L]
static PolyBank pfb_from_proto(const float *proto, int proto_len, int L){
    PolyBank pb = {0,0,NULL};
    if (L <= 0 || proto_len <= 0) return pb;

    int padded = proto_len;
    if (padded % L != 0) {
        padded = proto_len + (L - (proto_len % L));
    }
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

static void pfb_free(PolyBank *pb){
    if (!pb) return;
    if (pb->b){ free(pb->b); pb->b = NULL; }
    pb->L = 0; pb->tpp = 0;
}

static inline void pfb_dot_phase(const PolyBank *pb,
                                 const float *xIQ, size_t center_idx,
                                 int phase, float *outI, float *outQ)
{
    const float *taps = &pb->b[phase * pb->tpp];
    int half = (pb->tpp - 1)/2;

    double accI = 0.0, accQ = 0.0;
    int k;
    for (k = -half; k <= half; k++){
        long ii = (long)center_idx + k;
        if (ii < 0) continue;
        const float *s = &xIQ[(size_t)ii * 2];
        double hk = (double)taps[k + half];
        accI += hk * (double)s[0];
        accQ += hk * (double)s[1];
    }
    *outI = (float)accI;
    *outQ = (float)accQ;
}

/* -------------------------- PFB clock sync state --------------------------- */
typedef struct {
    // loop state
    double omega, omega_min, omega_max;  // samples/symbol estimate
    double mu;                           // fractional timing [0,1)
    double gain_mu, gain_omega;          // loop gains (alpha, beta)

    // bank
    int sps_nom;         // nominal sps
    int span_syms;       // RRC span (symbols)
    int L;               // phases
    int proto_len;       // original taps length
    float *proto;        // owned taps
    PolyBank bank;

    // memory for Gardner
    float yk_prev_I, yk_prev_Q;

    // inline AGC
    double agc_gain;
    double agc_ref;      // reference amplitude
    double agc_rate;     // adaptation rate
    int agc_enabled;

    // options
    int verbose;
} PFBClockSync;

static PFBClockSync* pfbcs_create(float alpha, int sps_nom, int span_syms, int L,
                                  double loop_bw, double damping, double omega_rel_lim,
                                  int enable_agc, double agc_rate, int verbose)
{
    if (sps_nom <= 0 || span_syms <= 0 || L <= 0) return NULL;

    PFBClockSync *cs = (PFBClockSync*)calloc(1, sizeof(PFBClockSync));
    if (!cs) return NULL;

    cs->sps_nom = sps_nom;
    cs->span_syms = span_syms;
    cs->L = L;
    cs->verbose = verbose;

    // build prototype RRC taps (allocate here)
    int N = span_syms * sps_nom + 1;
    cs->proto = (float*)malloc(sizeof(float) * (size_t)N);
    if (!cs->proto){ free(cs); return NULL; }
    rrc_generate_taps(alpha, sps_nom, span_syms, cs->proto);
    cs->proto_len = N;

    // polyphase bank
    cs->bank = pfb_from_proto(cs->proto, cs->proto_len, L);

    // Loop filter (back to working values)
    double Kp = 1.0;
    double Ki = 2.0;
    double denom = 1.0 + 2.0*damping*loop_bw + loop_bw*loop_bw;
    cs->gain_mu    = (4.0*damping*loop_bw) / (denom * Kp);
    cs->gain_omega = (4.0*loop_bw*loop_bw) / (denom * Ki);

    cs->omega = (double)sps_nom;
    cs->omega_min = cs->omega * (1.0 - omega_rel_lim);
    cs->omega_max = cs->omega * (1.0 + omega_rel_lim);

    cs->mu = 0.5;
    cs->yk_prev_I = 0.0f;
    cs->yk_prev_Q = 0.0f;

    // AGC setup
    cs->agc_enabled = enable_agc;
    cs->agc_gain = 1.0;
    cs->agc_ref = 1.0;      // target amplitude
    cs->agc_rate = agc_rate;

    return cs;
}

static void pfbcs_free(PFBClockSync *cs){
    if (!cs) return;
    if (cs->proto) free(cs->proto);
    pfb_free(&cs->bank);
    free(cs);
}

// returns number of output symbols
// ===================== REPLACEMENT: pfbcs_process (no AGC) =====================
static size_t pfbcs_process(PFBClockSync *cs,
    const float *xIQ, size_t nin_complex,
    float *outIQ)
{
const int tpp  = cs->bank.tpp;
const int half = (tpp - 1)/2;

// guard so dot products always have taps available
size_t in_idx  = (size_t)(half + cs->sps_nom + 2);
size_t end_idx = nin_complex - (size_t)half - (size_t)cs->sps_nom - 2;
size_t out_count = 0;

// optional: normalize TED by a smoothed envelope to de-sensitize to level
double env_hat = 1.0;          // simple IIR of |y_k|^2
const double env_beta = 1e-3;  // small smoothing

while (in_idx < end_idx) {
// ---- center sample y_k ----
double pf  = cs->mu * (double)cs->L;
int phase  = (int)(pf);
if (phase < 0) phase = 0;
if (phase >= cs->L) phase = cs->L - 1;

float ykI=0.f, ykQ=0.f;
pfb_dot_phase(&cs->bank, xIQ, in_idx, phase, &ykI, &ykQ);

// ---- mid sample y_{k-1/2} (coarse: integer index + phase) ----
double half_omega = cs->omega * 0.5;
int half_sps      = (int)floor(half_omega);
size_t mid_idx    = in_idx - (size_t)half_sps;

double mu_mid = cs->mu + 0.5;
if (mu_mid >= 1.0) mu_mid -= 1.0;

double pf_mid = mu_mid * (double)cs->L;
int phase_mid = (int)(pf_mid);
if (phase_mid < 0) phase_mid = 0;
if (phase_mid >= cs->L) phase_mid = cs->L - 1;

float ymidI=0.f, ymidQ=0.f;
pfb_dot_phase(&cs->bank, xIQ, mid_idx, phase_mid, &ymidI, &ymidQ);

// ---- Gardner error on these matched-filter outputs ----
float dyI = ykI - cs->yk_prev_I;
float dyQ = ykQ - cs->yk_prev_Q;
double e  = (double)ymidI * (double)dyI + (double)ymidQ * (double)dyQ;

// optional: normalize TED by smoothed envelope (helps if input level drifts)
double p  = (double)ykI*ykI + (double)ykQ*ykQ;
env_hat = (1.0 - env_beta)*env_hat + env_beta * (p > 1e-12 ? p : env_hat);
double e_norm = e / (env_hat + 1e-12);

// ---- loop updates ----
cs->mu    += cs->gain_mu    * e_norm;
cs->omega += cs->gain_omega * e_norm;

if (cs->omega < cs->omega_min) cs->omega = cs->omega_min;
if (cs->omega > cs->omega_max) cs->omega = cs->omega_max;

// ---- output the symbol-rate sample (no AGC scaling) ----
outIQ[2*out_count + 0] = ykI;
outIQ[2*out_count + 1] = ykQ;
out_count++;

// ---- advance input pointer by whole samples from μ ----
cs->mu += cs->omega;
int adv = (int)floor(cs->mu);
if (adv < 1) adv = 1;
cs->mu -= (double)adv;
in_idx += (size_t)adv;

// store previous (for dy)
cs->yk_prev_I = ykI;
cs->yk_prev_Q = ykQ;

if (cs->verbose && (out_count % 1000 == 0)) {
printf("[clk] out=%zu  mu=%.4f  omega=%.4f  e=%+.6g\n",
out_count, cs->mu, cs->omega, e_norm);
}
}
return out_count;
}


/* ------------------------------- CLI helpers ------------------------------- */
static bool flag_present(int argc, char **argv, const char *flag){
    int i;
    for (i = 1; i < argc; i++){
        if (strcmp(argv[i], flag) == 0) return true;
    }
    return false;
}
static const char* get_opt(int argc, char **argv, const char *key, const char *defv){
    int i;
    for (i = 1; i+1 < argc; i++){
        if (strcmp(argv[i], key) == 0) return argv[i+1];
    }
    return defv;
}
static void usage(const char *p){
    fprintf(stderr,
        "Usage: %s in.cf32 out.cf32 [--sps N] [--span K] [--alpha A] [--L P]\n"
        "          [--loopbw BW] [--damp ZETA] [--orelim R] [--agcrate R]\n"
        "          [--noagc] [--notrim] [--quiet]\n", p);
}

/* ----------------------------------- main ---------------------------------- */
int main(int argc, char **argv)
{
    if (argc < 3) { usage(argv[0]); return 1; }
    const char *in_path  = argv[1];
    const char *out_path = argv[2];

    // defaults
    int    SPS     = atoi(get_opt(argc, argv, "--sps",     "8"));
    int    SPAN    = atoi(get_opt(argc, argv, "--span",    "11"));
    float  ALPHA   = (float)atof(get_opt(argc, argv, "--alpha",  "0.5"));
    int    L       = atoi(get_opt(argc, argv, "--L",       "16"));
    double LOOPBW  = atof(get_opt(argc, argv, "--loopbw",  "0.0628"));
    double DAMP    = atof(get_opt(argc, argv, "--damp",    "0.707"));
    double ORELIM  = atof(get_opt(argc, argv, "--orelim",  "0.375"));
    double AGCRATE = atof(get_opt(argc, argv, "--agcrate", "0.01"));  // Faster AGC for startup

    int do_agc  = flag_present(argc, argv, "--noagc")  ? 0 : 1;
    int do_trim = flag_present(argc, argv, "--notrim") ? 0 : 1;
    int quiet   = flag_present(argc, argv, "--quiet")  ? 1 : 0;

    // read cf32 input
    FILE *fi = fopen(in_path, "rb");
    if (!fi){ perror("open in"); return 1; }
    if (fseek(fi, 0, SEEK_END) != 0){ perror("fseek"); fclose(fi); return 1; }
    long nbytes = ftell(fi);
    if (nbytes < 0){ perror("ftell"); fclose(fi); return 1; }
    if (fseek(fi, 0, SEEK_SET) != 0){ perror("fseek"); fclose(fi); return 1; }

    size_t nfloats = (size_t)nbytes / sizeof(float);
    if (nfloats % 2 != 0){ fprintf(stderr,"Input float count not even (I/Q)\n"); fclose(fi); return 1; }
    size_t nin_c = nfloats / 2;

    float *xIQ = (float*)malloc(sizeof(float) * nfloats);
    if (!xIQ){ fprintf(stderr,"malloc xIQ failed\n"); fclose(fi); return 1; }
    size_t rd = fread(xIQ, sizeof(float), nfloats, fi);
    fclose(fi);
    if (rd != nfloats){ fprintf(stderr,"short read\n"); free(xIQ); return 1; }

    // --- AGC on the input stream (do this BEFORE PFB/clock sync) ---
    StreamAGC *agc = agc_create(/*beta=*/1e-3, /*gamma=*/1e-3,
        /*target_rms=*/1.0, /*init_gain=*/1.0);
    if (!agc){ fprintf(stderr,"AGC alloc failed\n"); return 1; }

    // Optional one-shot bootstrap with first ~1000 symbols:
    agc_bootstrap(agc, xIQ, nin_c, 1000);

    // In-place normalize the whole input stream:
    agc_apply_inplace(agc, xIQ, nin_c);

    // (You can free now or keep it if you plan further calls)
    agc_free(agc);

    // Build clock sync with AGC disabled
    PFBClockSync *cs = pfbcs_create(ALPHA, SPS, SPAN, L, LOOPBW, DAMP, ORELIM,
            /*enable_agc=*/0, /*agc_rate=*/0.0,
            quiet ? 0 : 1);

    if (!cs){ fprintf(stderr,"pfbcs_create failed\n"); free(xIQ); return 1; }
    
    // Pre-estimate AGC gain from input signal to avoid startup transient
    if (do_agc && nin_c > 100) {
        size_t sample_count = nin_c < 1000 ? nin_c : 1000;
        double sum_mag = 0.0;
        size_t valid = 0;
        for (size_t i = 0; i < sample_count; i++) {
            double I = xIQ[2*i+0], Q = xIQ[2*i+1];
            double mag = sqrt(I*I + Q*Q);
            if (mag > 1e-9) {
                sum_mag += mag;
                valid++;
            }
        }
        if (valid > 10) {
            double avg_mag = sum_mag / (double)valid;
            // Set initial gain to get roughly unit amplitude after matched filter
            // The matched filter has gain ~sqrt(sps) for symbols
            double target_output = 1.0;
            double expected_mf_gain = sqrt((double)SPS) * 0.5; // approximate
            cs->agc_gain = target_output / (avg_mag * expected_mf_gain);
            if (cs->agc_gain < 0.01) cs->agc_gain = 0.01;
            if (cs->agc_gain > 100.0) cs->agc_gain = 100.0;
            if (!quiet) printf("Initial AGC gain from input: %.3f (input avg mag: %.6f)\n", 
                              cs->agc_gain, avg_mag);
        }
    }
    
    if (!quiet){
        printf("PFB: taps=%d -> phases=%d, taps/phase=%d\n", cs->proto_len, cs->bank.L, cs->bank.tpp);
        printf("Loop: gain_mu=%.6f, gain_omega=%.6f\n", cs->gain_mu, cs->gain_omega);
        if (do_agc) printf("AGC: enabled, rate=%.6f\n", AGCRATE);
        
        // Warn if taps/phase is too low
        if (cs->bank.tpp < 5) {
            printf("WARNING: Only %d taps/phase. Consider increasing --span or reducing --L\n", cs->bank.tpp);
            printf("         Recommended: span*sps/L >= 5, so span >= %d for L=%d\n", 
                   (5*L + SPS-1)/SPS, L);
        }
    }

    // process
    float *yIQ = (float*)malloc(sizeof(float) * 2 * nin_c);
    if (!yIQ){ fprintf(stderr,"malloc yIQ failed\n"); pfbcs_free(cs); free(xIQ); return 1; }
    size_t nout = pfbcs_process(cs, xIQ, nin_c, yIQ);

    // optional trimming of head/tail transients
    size_t start = 0, stop = nout;
    if (do_trim){
        size_t head = (size_t)(2 * SPAN);
        size_t tail = (size_t)(2 * SPAN);
        if (nout > head) start = head;
        if (nout > tail) stop  = nout - tail;
        if (stop < start){ start = 0; stop = nout; }
    }
    size_t valid = 0;
    if (stop > start) valid = stop - start;

    // write cf32 output
    FILE *fo = fopen(out_path, "wb");
    if (!fo){ perror("open out"); free(yIQ); pfbcs_free(cs); free(xIQ); return 1; }
    size_t wr = fwrite(&yIQ[2*start], sizeof(float), 2*valid, fo);
    fclose(fo);
    if (wr != 2*valid){ fprintf(stderr,"short write: %zu/%zu floats\n", wr, 2*valid); }

    if (!quiet){
        printf("Clock sync done. In=%zu complex, Out=%zu symbols, Wrote=%zu (trim=%s, agc=%s)\n",
               nin_c, nout, valid, do_trim ? "yes":"no", do_agc ? "yes":"no");
        if (do_agc) printf("Final AGC gain: %.3f\n", cs->agc_gain);
    }

    free(yIQ);
    pfbcs_free(cs);
    free(xIQ);
    return 0;
}