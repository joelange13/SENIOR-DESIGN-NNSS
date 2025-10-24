// slice_cf32.c
// Slice cf32 symbols to bits (BPSK I/Q) and ALSO scan bit-alignments
// for special word 0xD5, then attempt to parse seq/len/payload/crc.
// Writes the same sliced bitfile as before.
//
// Usage (examples):
//   ./slice_cf32 symbols.cf32 sliced.bits --mod bpsk --axis i --rotate 180
//
// Notes:
// - Threshold = 0.0 (change if you want).
// - Rotation 180 flips sign, 0 = as-is.
// - Axis: i or q.
// - CRC16: CCITT-FALSE (poly 0x1021, init 0xFFFF), payload-only.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>


static void die(const char* msg){ fprintf(stderr, "%s\n", msg); exit(1); }

typedef struct { float i,q; } c32;

static uint16_t CRC16_CCITT_FALSE(const uint8_t* buf, size_t len) {
    const uint16_t poly = 0x1021;
    uint16_t crc = 0xFFFF;
    for (size_t i = 0; i < len; i++) {
        crc ^= (uint16_t)buf[i] << 8;
        for (int b = 0; b < 8; b++) {
            crc = (crc & 0x8000) ? (uint16_t)((crc << 1) ^ poly) : (uint16_t)(crc << 1);
        }
    }
    return crc;
}

static void write_bits_msb(const uint8_t* bits, size_t nbits, const char* outpath){
    FILE* f = fopen(outpath, "wb");
    if(!f) die("Failed to open output bit file");
    uint8_t acc=0; int k=0;
    for(size_t i=0;i<nbits;i++){
        acc = (uint8_t)((acc<<1) | (bits[i]?1:0));
        if(++k==8){ fwrite(&acc,1,1,f); acc=0; k=0; }
    }
    if(k){ acc <<= (8-k); fwrite(&acc,1,1,f); }
    fclose(f);
}

int main(int argc, char** argv){
    if(argc < 3) {
        fprintf(stderr, "Usage: %s <in.cf32> <out.bits> [--mod bpsk] [--axis i|q] [--rotate 0|180]\n", argv[0]);
        return 1;
    }

    const char* inpath  = argv[1];
    const char* outbits = argv[2];

    // defaults
    const char* mod = "bpsk";
    char axis = 'i';
    double rotate_deg = 0.0;

    // parse simple flags
    for(int i=3;i<argc;i++){
        if(strcmp(argv[i],"--mod")==0 && i+1<argc) { mod = argv[++i]; }
        else if(strcmp(argv[i],"--axis")==0 && i+1<argc) { axis = (char)tolower(argv[++i][0]); }
        else if(strcmp(argv[i],"--rotate")==0 && i+1<argc) { rotate_deg = atof(argv[++i]); }
    }

    if (strcasecmp(mod,"bpsk")!=0) {
        die("Only --mod bpsk is supported in this build.");
    }
    if(axis!='i' && axis!='q') die("axis must be i or q");

    // read cf32
    FILE* f = fopen(inpath, "rb");
    if(!f) die("Failed to open input cf32");
    fseek(f,0,SEEK_END);
    long bytes = ftell(f);
    fseek(f,0,SEEK_SET);
    if(bytes < 0 || bytes % (int)sizeof(float) != 0) die("Bad file size");
    size_t nfloats = (size_t)bytes / sizeof(float);
    if(nfloats % 2 != 0) die("cf32 must have even number of floats");
    size_t nsamps = nfloats/2;

    float* buf = (float*)malloc(bytes);
    if(!buf) die("OOM");
    if(fread(buf,sizeof(float),nfloats,f) != nfloats) die("Read error");
    fclose(f);

    c32* x = (c32*)buf;

    // slice to bits
    const double rot = fmod(rotate_deg, 360.0);
    int flip = (fabs(rot-180.0) < 1e-6) ? -1 : +1; // only 0 or 180 supported

    uint8_t* bits = (uint8_t*)malloc(nsamps);
    if(!bits) die("OOM bits");

    size_t nbits = 0;
    for(size_t n=0;n<nsamps;n++){
        float v = (axis=='i') ? x[n].i : x[n].q;
        v = (float)(flip * v);
        bits[nbits++] = (v < 0.0f) ? 1 : 0; // BPSK: map sign to bit
    }

    printf("Sliced %zu symbols (BPSK). Produced %zu bits.\n", nsamps, nbits);
    printf("Options: rotate=%.1f°, axis=%c, MSB-first\n", rotate_deg, axis);

    // always write the sliced bits (packed MSB-first) to the requested file
    write_bits_msb(bits, nbits, outbits);

    // -------- bit alignment scan for special word 0xD5 --------
    // We’ll try all 8 bit offsets; at each offset, repack into bytes, scan for 0xD5.
    // If found, try to parse: seq(lo,hi), len, payload[len], crc(lo,hi).
    // Note: capture may not include a whole packet; we’ll report if incomplete.

    int total_hits = 0;
    for(int bitoff = 0; bitoff < 8; bitoff++){
        size_t usable_bits = (nbits > (size_t)bitoff) ? (nbits - (size_t)bitoff) : 0;
        size_t nbytes = usable_bits / 8;
        if(nbytes == 0) continue;

        uint8_t* bytes8 = (uint8_t*)malloc(nbytes);
        if(!bytes8) die("OOM bytes8");
        // pack MSB-first starting at bitoff
        for(size_t B=0; B<nbytes; B++){
            uint8_t acc = 0;
            for(int b=0;b<8;b++){
                size_t bi = (size_t)bitoff + B*8 + b;
                acc = (uint8_t)((acc<<1) | (bits[bi] ? 1 : 0));
            }
            bytes8[B] = acc;
        }

        // scan for 0xD5
        size_t hits = 0;
        for(size_t i=0;i<nbytes;i++){
            if(bytes8[i] == 0xD5){
                // Optional: check for preceding AA run (commented)
                // int aa_ok = 0;
                // if(i>=7){
                //     aa_ok = 1;
                //     for(int k=1;k<=7;k++) if(bytes8[i-k] != 0xAA) { aa_ok = 0; break; }
                // }

                // Try to parse a packet starting *after* D5:
                // [D5][seq_lo][seq_hi][len][payload...][crc_lo][crc_hi]
                size_t pos = i + 1;
                if(pos + 3 > nbytes){ // need at least seq_lo, seq_hi, len
                    // incomplete header
                } else {
                    uint16_t seq = (uint16_t)bytes8[pos] | ((uint16_t)bytes8[pos+1] << 8);
                    uint8_t  len = bytes8[pos+2];
                    size_t need = pos + 3 + (size_t)len + 2; // header + payload + crc
                    if(need <= nbytes){
                        const uint8_t* payload = &bytes8[pos+3];
                        uint16_t rx_crc = (uint16_t)payload[len] | ((uint16_t)payload[len+1] << 8);
                        uint16_t calc = CRC16_CCITT_FALSE(payload, len);

                        printf("\n[hit] bit_offset=%d, byte_index=%zu  (global bit=%zu)\n",
                               bitoff, i, (size_t)bitoff + i*8);
                        printf("      D5 found. seq=%u (0x%04X), len=%u\n", seq, seq, len);
                        printf("      CRC rx=0x%04X calc=0x%04X  %s\n",
                               rx_crc, calc, (rx_crc==calc)?"OK":"FAIL");

                        // Show first up to 32 payload bytes
                        size_t show = (len < 32) ? len : 32;
                        printf("      payload[0:%zu]:", show);
                        for(size_t t=0;t<show;t++){
                            printf(" %02X", payload[t]);
                        }
                        printf("%s\n", (show<len)?" ...":"");

                        hits++;
                        total_hits++;
                    } else {
                        // header OK but not enough bytes for full payload+crc
                        printf("\n[hit] bit_offset=%d, byte_index=%zu  (global bit=%zu)\n",
                               bitoff, i, (size_t)bitoff + i*8);
                        printf("      D5 found but capture ends before full packet. ");
                        printf("Have %zu bytes after D5; need %zu.\n",
                               (nbytes - (i+1)), (size_t)3 + (size_t)len + 2);
                        hits++;
                        total_hits++;
                    }
                }
            }
        }

        if(hits==0){
            // keep quiet per offset; we’ll summarize after loop
        }
        free(bytes8);
    }

    if(total_hits==0){
        printf("No D5-based packets parsed. ");
        printf("Either alignment still off wrt packet start or capture too short.\n");
    } else {
        printf("\nDone. Found %d D5 markers across bit alignments.\n", total_hits);
    }

    free(bits);
    free(buf);
    return 0;
}
