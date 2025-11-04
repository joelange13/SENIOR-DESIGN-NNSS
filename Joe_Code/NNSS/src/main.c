#include "api.h"

/* Defaults (unchanged) */
#define DEFAULT_TX_IN   "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\hello.txt"
#define DEFAULT_TX_OUT  "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\goodbye.cfile"
#define DEFAULT_RX_IN   "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\goodbye_offset.cfile"
#define DEFAULT_RX_OUT  "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\rx_out.cfile"

/* Extra diagnostics */
#define DEFAULT_FLL_DIAG   "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\fll_corrected.cfile"
#define DEFAULT_COSTAS_OUT "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\costas_out.cfile"

static void usage_tx(const char *p){
    fprintf(stderr, "Usage (TX): %s --mode tx [--in path] [--out path]\n", p);
}
static void usage_rx(const char *p){
    fprintf(stderr, "Usage (RX): %s --mode rx [--in path] [--out path] "
                    "[--sps N] [--span K] [--alpha A] [--L P] "
                    "[--loopbw BW] [--damp ZETA] [--orelim R] "
                    "[--fll_bw B] [--be_frac R] "
                    "[--costas_bw B] [--costas_agc A] [--nocostas] "
                    "[--noagc] [--notrim] [--quiet]\n", p);
}

static const char* get_opt_str(int argc, char **argv, const char *key, const char *defv){
    for (int i = 1; i+1 < argc; i++) if (strcmp(argv[i], key) == 0) return argv[i+1];
    return defv;
}

int main(int argc, char **argv){
    if (!flag_present(argc, argv, "--mode")){
        fprintf(stderr, "Missing --mode {tx|rx}\n");
        usage_tx(argv[0]); usage_rx(argv[0]);
        return 1;
    }
    const char *mode = get_opt(argc, argv, "--mode", "");

    if (strcmp(mode, "tx") == 0){
        /* ----- TX PATH (defaults preserved) ----- */
        const char* in_path  = get_opt_str(argc, argv, "--in",  DEFAULT_TX_IN);
        const char* out_path = get_opt_str(argc, argv, "--out", DEFAULT_TX_OUT);

        FILE* inputFile  = OpenFile(in_path,  true);
        FILE* outputFile = OpenFile(out_path, false);
        if (!inputFile || !outputFile) {
            if (inputFile)  fclose(inputFile);
            if (outputFile) fclose(outputFile);
            return 1;
        }

        build_CRC16_table();
        Packet pkt;
        const uint8_t preamble[8] = {0xAA,0xAA,0xAA,0xAA,0xAA,0xAA,0xAA,0xD5};
        uint16_t sequence = 0;

        const int   sps   = 8;
        const float alpha = 0.5f;
        const int   span  = 11;

        for (;;) {
            size_t nread = fread(pkt.payload, 1, sizeof(pkt.payload), inputFile);
            if (nread == 0) break;

            memcpy(pkt.preamble, preamble, sizeof(pkt.preamble));
            pkt.sequence = sequence++;
            pkt.length   = (uint8_t)nread;
            pkt.crc      = CRC16(pkt);

            const size_t need_bits = 8u*sizeof(pkt.preamble) + 16u + 8u +
                                     (size_t)pkt.length * 8u + 16u;
            bool *packet_bits = (bool*)malloc(need_bits * sizeof(bool));
            if (!packet_bits){ fprintf(stderr, "malloc packet_bits failed\n"); break; }

            size_t num_bits = packet_to_bits(&pkt, packet_bits, need_bits);

            #if USE_SCRAMBLER
                scramble_self_sync_after_offset(packet_bits, num_bits, PREAMBLE_BITS);
            #elif USE_MANCHESTER
                bool *coded = NULL;
                size_t coded_bits = manchester_encode_after_offset(packet_bits, num_bits, PREAMBLE_BITS, &coded);
                if (!coded_bits) { fprintf(stderr,"Manchester malloc failed\n"); free(packet_bits); break; }
                free(packet_bits); packet_bits = coded; num_bits = coded_bits;
            #endif

            debugPrints(pkt, packet_bits, (int)num_bits, NULL, 0, NULL, 0);
            if (num_bits == 0){ free(packet_bits); continue; }

            size_t out_iq_floats = 2u * (size_t)num_bits * (size_t)sps;
            float *out_iq = (float*)malloc(out_iq_floats * sizeof(float));
            if (!out_iq){ fprintf(stderr,"malloc out_iq failed\n"); free(packet_bits); break; }

            size_t num_iq_samps = bits_to_float_bpsk_iq(packet_bits, num_bits, out_iq, sps);
            float *shaped_out_iq=NULL; size_t shaped_len=0;
            if (num_iq_samps > 0){
                tx_pulse_shape(out_iq, 2*(size_t)num_iq_samps, sps, alpha, span, &shaped_out_iq, &shaped_len);
                write_cf32(outputFile, shaped_out_iq, shaped_len);
            }
            free(shaped_out_iq); free(out_iq); free(packet_bits);
        }
        fclose(inputFile); fclose(outputFile);
        return 0;
    }

    if (strcmp(mode, "rx") == 0){
        /* ----- RX PATH ----- */
        const char *in_path   = get_opt_str(argc, argv, "--in",  DEFAULT_RX_IN);
        const char *out_path  = get_opt_str(argc, argv, "--out", DEFAULT_RX_OUT);

        int    SPS     = atoi(get_opt(argc, argv, "--sps",     "8"));
        int    SPAN    = atoi(get_opt(argc, argv, "--span",    "11"));
        float  ALPHA   = (float)atof(get_opt(argc, argv, "--alpha",  "0.5"));
        int    L       = atoi(get_opt(argc, argv, "--L",       "16"));
        double LOOPBW  = atof(get_opt(argc, argv, "--loopbw",  "0.0628"));
        double DAMP    = atof(get_opt(argc, argv, "--damp",    "0.707"));
        double ORELIM  = atof(get_opt(argc, argv, "--orelim",  "0.375"));

        /* FLL tuning (optional CLI) */
        double FLL_BW    = atof(get_opt(argc, argv, "--fll_bw",   "0.004"));  /* per-sample, pre-PFB */
        float  BE_FRAC   = (float)atof(get_opt(argc, argv, "--be_frac", "1.0"));

        /* Costas tuning (symbol-rate domain) */
        double COSTAS_BW  = atof(get_opt(argc, argv, "--costas_bw",  "0.01"));
        double COSTAS_AGC = atof(get_opt(argc, argv, "--costas_agc", "0.001"));

        int do_agc    = flag_present(argc, argv, "--noagc")    ? 0 : 1;
        int do_trim   = flag_present(argc, argv, "--notrim")   ? 0 : 1;
        int do_costas = flag_present(argc, argv, "--nocostas") ? 0 : 1;
        int quiet     = flag_present(argc, argv, "--quiet")    ? 1 : 0;

        FILE *fi = fopen(in_path, "rb");
        if (!fi){ perror("open in"); return 1; }
        fseek(fi, 0, SEEK_END); long nbytes = ftell(fi);
        if (nbytes < 0){ perror("ftell"); fclose(fi); return 1; }
        fseek(fi, 0, SEEK_SET);
        size_t nfloats = (size_t)nbytes / sizeof(float);
        if (nfloats % 2 != 0){ fprintf(stderr,"Input float count not even (I/Q)\n"); fclose(fi); return 1; }
        size_t nin_c = nfloats / 2;

        float *xIQ = (float*)malloc(sizeof(float) * nfloats);
        if (!xIQ){ fprintf(stderr,"malloc xIQ failed\n"); fclose(fi); return 1; }
        size_t rd = fread(xIQ, sizeof(float), nfloats, fi);
        fclose(fi);
        if (rd != nfloats){ fprintf(stderr,"short read\n"); free(xIQ); return 1; }

        /* --- AGC (optional) --- */
        if (do_agc){
            StreamAGC *agc = agc_create(/*beta*/1e-3, /*gamma*/2.0, /*target_rms*/0.5, /*init_gain*/1.0);
            if (!agc){ fprintf(stderr,"agc_create failed\n"); free(xIQ); return 1; }
            agc_bootstrap(agc, xIQ, nin_c, (size_t)(8*SPS));
            agc_apply_inplace(agc, xIQ, nin_c);
            agc_free(agc);
            if (!quiet) fprintf(stderr,"AGC applied.\n");
        }

        /* --- FLL Band-Edge (diagnostic output + feed PFB) --- */
        FLLBandEdge *fll = fll_be_create(ALPHA, SPS, SPAN, BE_FRAC, FLL_BW, /*damping*/0.707, /*pwr_beta*/0.0015, quiet?0:1);
        if (!fll){ fprintf(stderr,"fll_be_create failed\n"); free(xIQ); return 1; }

        float *xIQ_corr = (float*)malloc(sizeof(float) * 2 * nin_c);
        if (!xIQ_corr){ fprintf(stderr,"malloc xIQ_corr failed\n"); fll_be_free(fll); free(xIQ); return 1; }

        size_t nfix = fll_be_process(fll, xIQ, nin_c, xIQ_corr);
        if (nfix != nin_c) { fprintf(stderr,"fll_be_process wrote %zu of %zu\n", nfix, nin_c); }

        if (!quiet){
            const char *fll_out_path = DEFAULT_FLL_DIAG;
            FILE *fll_out = fopen(fll_out_path, "wb");
            if (fll_out) {
                size_t fll_wr = fwrite(xIQ_corr, sizeof(float), 2*nfix, fll_out);
                fclose(fll_out);
                printf("Wrote %zu complex samples to: %s\n", nfix, fll_out_path);
            } else {
                fprintf(stderr, "Warning: could not write FLL output file\n");
            }
            printf("FLL: est freq = %.6g cps (%.6g rad/spl)\n",
                   fll_be_get_freq_cps(fll), fll_be_get_freq_rad(fll));
        }

        /* Replace input with corrected for downstream timing recovery */
        free(xIQ); xIQ = xIQ_corr; nin_c = nfix;

        /* --- PFB timing recovery --- */
        PFBClockSync *cs = pfbcs_create(ALPHA, SPS, SPAN, L, LOOPBW, DAMP, ORELIM,
                                        /*enable_agc=*/0, /*agc_rate=*/0.0,
                                        quiet ? 0 : 1);
        if (!cs){ fprintf(stderr,"pfbcs_create failed\n"); fll_be_free(fll); free(xIQ); return 1; }

        float *yIQ_pfb = (float*)malloc(sizeof(float) * 2 * nin_c);
        if (!yIQ_pfb){ fprintf(stderr,"malloc yIQ_pfb failed\n"); pfbcs_free(cs); fll_be_free(fll); free(xIQ); return 1; }
        size_t nout = pfbcs_process(cs, xIQ, nin_c, yIQ_pfb);

        /* Optional transient trim (applies to both PFB and Costas files) */
        size_t start = 0, stop = nout;
        if (do_trim){
            size_t head = (size_t)(2 * SPAN);
            size_t tail = (size_t)(2 * SPAN);
            if (nout > head) start = head;
            if (nout > tail) stop  = nout - tail;
            if (stop < start){ start = 0; stop = nout; }
        }
        size_t valid = (stop > start) ? (stop - start) : 0;

        /* --- Write PFB-only output (unchanged path) --- */
        {
            FILE *fo = fopen(out_path, "wb");
            if (!fo){ perror("open out"); free(yIQ_pfb); pfbcs_free(cs); fll_be_free(fll); free(xIQ); return 1; }
            size_t wr = fwrite(&yIQ_pfb[2*start], sizeof(float), 2*valid, fo);
            fclose(fo);
            if (wr != 2*valid){ fprintf(stderr,"short write: %zu/%zu floats (PFB)\n", wr, 2*valid); }
        }

        /* --- Costas (BPSK) on PFB output -> new diagnostic file --- */
        if (do_costas){
            CostasBPSK *cl = costas_bpsk_create(COSTAS_BW, 0.707, COSTAS_AGC, quiet?0:1);
            if (!cl){
                fprintf(stderr,"costas_bpsk_create failed\n");
                free(yIQ_pfb); pfbcs_free(cs); fll_be_free(fll); free(xIQ);
                return 1;
            }
            float *yIQ_costas = (float*)malloc(sizeof(float) * 2 * nout);
            if (!yIQ_costas){
                fprintf(stderr,"malloc yIQ_costas failed\n");
                costas_bpsk_free(cl); free(yIQ_pfb); pfbcs_free(cs); fll_be_free(fll); free(xIQ);
                return 1;
            }
            size_t nz = costas_bpsk_process(cl, yIQ_pfb, nout, yIQ_costas);

            /* Write Costas-corrected stream */
            const char *costas_path = DEFAULT_COSTAS_OUT;
            FILE *fc = fopen(costas_path, "wb");
            if (!fc){
                perror("open costas");
            } else {
                /* Trim with same start/stop so time windows match */
                size_t wr = fwrite(&yIQ_costas[2*start], sizeof(float), 2*valid, fc);
                fclose(fc);
                if (wr != 2*valid){ fprintf(stderr,"short write: %zu/%zu floats (Costas)\n", wr, 2*valid); }
                if (!quiet){
                    printf("Costas: phase=%.6f rad  freq=%.6g cps  (wrote %zu sym) -> %s\n",
                        (float)costas_bpsk_get_phase_rad(cl),
                        (float)costas_bpsk_get_freq_cps(cl),
                        valid, costas_path);
                }
            }

            free(yIQ_costas);
            costas_bpsk_free(cl);
        }

        if (!quiet){
            printf("Clock sync done. In=%zu complex, Out=%zu symbols, Wrote=%zu (trim=%s, agc=%s, costas=%s)\n",
                   nin_c, nout, valid, do_trim ? "yes":"no", do_agc ? "yes":"no", do_costas ? "yes":"no");
        }

        /* Clean up */
        free(yIQ_pfb);
        pfbcs_free(cs);
        fll_be_free(fll);
        free(xIQ);
        return 0;
    }

    fprintf(stderr, "Unknown --mode %s\n", mode);
    usage_tx(argv[0]); usage_rx(argv[0]);
    return 1;
}
