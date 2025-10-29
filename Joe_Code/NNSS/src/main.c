#include "api.h"

/* Defaults (your originals for TX) */
#define DEFAULT_TX_IN  "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\hello.txt"
#define DEFAULT_TX_OUT "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\goodbye.cfile"

/* RX sensible defaults (use TX output as RX input, write a new file) */
#define DEFAULT_RX_IN  "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\goodbye.cfile"
#define DEFAULT_RX_OUT "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\rx_out.cfile"

static void usage_tx(const char *p);

static void usage_tx(const char *p){
    fprintf(stderr, "Usage (TX): %s --mode tx [--in path] [--out path]\n", p);
}

static void usage_rx(const char *p){
    fprintf(stderr, "Usage (RX): %s --mode rx [--in path] [--out path] "
                    "[--sps N] [--span K] [--alpha A] [--L P] "
                    "[--loopbw BW] [--damp ZETA] [--orelim R] "
                    "[--noagc] [--notrim] [--quiet]\n", p);
}


/* Small helpers so we can use flags for in/out too */
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
        /* ----- RX PATH (defaults aligned to your TX output) ----- */
        const char *in_path  = get_opt_str(argc, argv, "--in",  DEFAULT_RX_IN);
        const char *out_path = get_opt_str(argc, argv, "--out", DEFAULT_RX_OUT);

        int    SPS     = atoi(get_opt(argc, argv, "--sps",     "8"));
        int    SPAN    = atoi(get_opt(argc, argv, "--span",    "11"));
        float  ALPHA   = (float)atof(get_opt(argc, argv, "--alpha",  "0.5"));
        int    L       = atoi(get_opt(argc, argv, "--L",       "16"));
        double LOOPBW  = atof(get_opt(argc, argv, "--loopbw",  "0.0628"));
        double DAMP    = atof(get_opt(argc, argv, "--damp",    "0.707"));
        double ORELIM  = atof(get_opt(argc, argv, "--orelim",  "0.375"));

        int do_agc  = flag_present(argc, argv, "--noagc")  ? 0 : 1;
        int do_trim = flag_present(argc, argv, "--notrim") ? 0 : 1;
        int quiet   = flag_present(argc, argv, "--quiet")  ? 1 : 0;

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

        if (do_agc){
            StreamAGC *agc = agc_create(1e-3, 1e-3, 1.0, 1.0);
            if (!agc){ fprintf(stderr,"AGC alloc failed\n"); free(xIQ); return 1; }
            agc_bootstrap(agc, xIQ, nin_c, 1000);
            agc_apply_inplace(agc, xIQ, nin_c);
            agc_free(agc);
        }

        PFBClockSync *cs = pfbcs_create(ALPHA, SPS, SPAN, L, LOOPBW, DAMP, ORELIM,
                                        /*enable_agc=*/0, /*agc_rate=*/0.0,
                                        quiet ? 0 : 1);
        if (!cs){ fprintf(stderr,"pfbcs_create failed\n"); free(xIQ); return 1; }

        float *yIQ = (float*)malloc(sizeof(float) * 2 * nin_c);
        if (!yIQ){ fprintf(stderr,"malloc yIQ failed\n"); pfbcs_free(cs); free(xIQ); return 1; }
        size_t nout = pfbcs_process(cs, xIQ, nin_c, yIQ);

        size_t start = 0, stop = nout;
        if (do_trim){
            size_t head = (size_t)(2 * SPAN);
            size_t tail = (size_t)(2 * SPAN);
            if (nout > head) start = head;
            if (nout > tail) stop  = nout - tail;
            if (stop < start){ start = 0; stop = nout; }
        }
        size_t valid = (stop > start) ? (stop - start) : 0;

        FILE *fo = fopen(out_path, "wb");
        if (!fo){ perror("open out"); free(yIQ); pfbcs_free(cs); free(xIQ); return 1; }
        size_t wr = fwrite(&yIQ[2*start], sizeof(float), 2*valid, fo);
        fclose(fo);
        if (wr != 2*valid){ fprintf(stderr,"short write: %zu/%zu floats\n", wr, 2*valid); }

        if (!quiet){
            printf("Clock sync done. In=%zu complex, Out=%zu symbols, Wrote=%zu (trim=%s, agc=%s)\n",
                   nin_c, nout, valid, do_trim ? "yes":"no", do_agc ? "yes":"no");
        }

        free(yIQ); pfbcs_free(cs); free(xIQ);
        return 0;
    }

    fprintf(stderr, "Unknown --mode '%s' (use 'tx' or 'rx')\n", mode);
    usage_tx(argv[0]); usage_rx(argv[0]);
    return 1;
}
