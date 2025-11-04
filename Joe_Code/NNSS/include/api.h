#ifndef API_H
#define API_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#define _USE_MATH_DEFINES
#include <math.h>

/* -------- TX build-time flags -------- */
#define USE_SCRAMBLER  0
#define USE_MANCHESTER 0
#define PREAMBLE_BYTES 8u
#define PREAMBLE_BITS  (PREAMBLE_BYTES * 8u)

/* -------- RX build-time flags -------- */


/* -------- Packet (TX) -------- */
typedef struct{
    uint8_t  preamble[8];
    uint16_t sequence;
    uint8_t  length;
    uint8_t  payload[253];
    uint16_t crc;
} Packet;

/* -------- CRC16 globals (TX) -------- */
extern uint16_t CRC16_table[256];
extern bool     crc_table_built;

/* ===================== TX prototypes (existing) ===================== */
FILE*   OpenFile(const char* path, bool is_input_file);
void    build_CRC16_table(void);
uint16_t CRC16(Packet pkt);
int     byte_to_bits(uint8_t byte, bool* packet_bits, size_t bits_cap, size_t* bit_pos);
size_t  packet_to_bits(const Packet *pkt, bool *packet_bits, size_t bits_cap);
size_t     bits_to_float_bpsk_iq(bool *packet_bits, size_t num_bits, float* out_iq, int sps);
void    rrc_generate_taps(float alpha, int sps, int span, float *taps_out);
void    iq_interleaved_fir_filter(float *x, size_t in_floats, float *h, int num_taps, float *y);
void    tx_pulse_shape(float *out_iq, size_t nsamps_unshaped, int sps, float alpha, int span,
                       float **shaped_out_iq, size_t *nsamps_shaped_out);
void    scramble_self_sync_after_offset(bool *bits, size_t n_bits, size_t offset);
size_t  manchester_encode_after_offset(const bool *in_bits, size_t n_bits,
                                       size_t offset, bool **out_bits);
void    debugPrints(Packet pkt, bool *packet_bits, size_t num_bits,
                    float *out_iq, int num_iq_samps, float *shaped_out_iq, size_t shaped_out_iq_length);
void    write_cf32(FILE *file, float *shaped_out_iq, size_t shaped_out_iq_length);

/* ===================== RX types ===================== */
typedef struct {
    double p_hat;
    double gain;
    double beta;
    double gamma;
    double target_power;
    double min_gain, max_gain;
} StreamAGC;

typedef struct PFBClockSync PFBClockSync; /* opaque */

/* ===================== RX prototypes (one-per-file) ===================== */
/* AGC */
StreamAGC* agc_create(double beta, double gamma, double target_rms, double init_gain);
void       agc_free(StreamAGC *a);
void       agc_bootstrap(StreamAGC *a, const float *xIQ, size_t nin_complex, size_t sample_cap);
void       agc_apply_inplace(StreamAGC *a, float *xIQ, size_t nin_complex);

/* PFB clock sync */
PFBClockSync* pfbcs_create(float alpha, int sps_nom, int span_syms, int L,
                           double loop_bw, double damping, double omega_rel_lim,
                           int enable_agc, double agc_rate, int verbose);
void          pfbcs_free(PFBClockSync *cs);
/* Returns number of output symbols written to outIQ (I/Q interleaved) */
size_t        pfbcs_process(PFBClockSync *cs, const float *xIQ, size_t nin_complex, float *outIQ);

/* Line coding (RX) */
void   descramble_self_sync_after_offset(bool *s, size_t n, size_t offset);
size_t manchester_decode_after_offset(const bool *in, size_t n, size_t offset, bool **outp);

/* CLI helpers (shared) */
bool        flag_present(int argc, char **argv, const char *flag);
const char* get_opt(int argc, char **argv, const char *key, const char *defv);

/* ===================== FLL Band-Edge ===================== */
/* Opaque handle */
typedef struct FLLBandEdge FLLBandEdge;

/* Create a band-edge FLL.
 * alpha   : RRC rolloff (0..1)
 * sps     : samples per symbol
 * span    : RRC span in symbols (FIR length ~= span*sps+1, clamped to >=201)
 * be_frac : 1.0 = place filters at the theoretical band-edges;
 *           <1.0 pulls them slightly inward (e.g., 0.95)
 * loop_bw : normalized loop BW per-sample (heuristic mapping to Kp/Ki)
 * damping : currently unused in mapping (reserved)
 * pwr_beta: single-pole IIR for power smoothing (0..1), e.g., 0.0015
 * verbose : nonzero -> prints settings on create
 */
FLLBandEdge* fll_be_create(float alpha, int sps, int span,
                           float be_frac, double loop_bw, double damping,
                           double pwr_beta, int verbose);

/* Process IQ (I/Q interleaved floats); writes frequency-corrected IQ.
 * Returns number of complex samples written (== nin_complex).
 */
size_t       fll_be_process(FLLBandEdge *fll,
                            const float *xIQ, size_t nin_complex,
                            float *yIQ);

/* Destroy */
void         fll_be_free(FLLBandEdge *fll);

/* Frequency getters:
 * cycles/sample and rad/sample (set Hz = cps * Fs separately if desired).
 */
double       fll_be_get_freq_cps(const FLLBandEdge *fll);
double       fll_be_get_freq_rad(const FLLBandEdge *fll);

/* ===================== Costas Loop (BPSK, DD) ===================== */
typedef struct CostasBPSK CostasBPSK;

/* loop_bw: normalized per-sample (post-PFB, ~1 sample/symbol)
   damping: typical 0.707
   agc_beta: small power smoother (e.g., 1e-3) for internal RMS tracking
   verbose: print settings once
*/
CostasBPSK* costas_bpsk_create(double loop_bw, double damping, double agc_beta, int verbose);

/* Process nin complex samples (I/Q interleaved) -> derotated output in yIQ.
   Returns nin (same length). */
size_t      costas_bpsk_process(CostasBPSK *c,
                                const float *xIQ, size_t nin_complex,
                                float *yIQ);

/* Destroy */
void        costas_bpsk_free(CostasBPSK *c);

/* Inspectors */
double      costas_bpsk_get_phase_rad(const CostasBPSK *c);   /* current NCO phase */
double      costas_bpsk_get_freq_cps(const CostasBPSK *c);    /* cycles/sample */



#endif /* API_H */
