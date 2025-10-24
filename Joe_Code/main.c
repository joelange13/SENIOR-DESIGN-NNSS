#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#define _USE_MATH_DEFINES
#include <math.h>

typedef struct{
    uint8_t preamble[8];
    uint16_t sequence;
    uint8_t length;
    uint8_t payload[253];
    uint16_t crc;
} Packet;

FILE* OpenFile(const char* path, bool is_input_file){
    if(is_input_file){
        FILE* ioFile = fopen(path, "rb");
        if(!ioFile){
            perror("Error Opening Input File\n");
            return NULL;
        }else{
            printf("Opened Input File.\n");
            return ioFile;
        }
    }else{
        FILE* ioFile = fopen(path, "wb");
        if(!ioFile){
            perror("Error Opening Output File\n");
            return NULL;
        }else{
            printf("Opened Output File.\n");
            return ioFile;
        }
    }
}

static uint16_t CRC16_table[256];
bool crc_table_built = false;

static void build_CRC16_table(){
    uint16_t crc_polynomial = 0x1021;
    uint16_t crc;
    for(uint16_t i = 0; i < 256; ++i){
        crc = (uint16_t)(i << 8);
        for(int b = 0; b < 8; ++b){
            if((crc & 0x8000) !=0){
                crc = (uint16_t)((crc << 1)^crc_polynomial);
            }else{
                crc = (uint16_t)(crc << 1);
            }
        }
        CRC16_table[i] = crc;
    }
    crc_table_built = true;
}

uint16_t CRC16(Packet pkt){
    uint16_t crc = 0xFFFF;
    for(int i = 0; i < pkt.length; i++){
        uint8_t idx = (uint8_t)(((crc >> 8) & 0xFF)^pkt.payload[i]);
        crc = (uint16_t)((crc << 8) ^ CRC16_table[idx]);
    }
    return crc;
}

int byte_to_bits(uint8_t byte, bool* packet_bits, size_t bits_cap, size_t* bit_pos){
    int b = 0; 
    while(b < 8){
        if(*bit_pos >= bits_cap){
            return -1;
        }
        uint8_t mask = (uint8_t)(0x80u >> b);
        if((byte & mask) != 0){
            packet_bits[*bit_pos] = true;
        }else{
            packet_bits[*bit_pos] = false;
        }
        (*bit_pos)++;
        b++;
    }
    return 0;
}

size_t packet_to_bits(const Packet *pkt, bool *packet_bits, size_t bits_cap ){
    size_t need_bits = 64u + 16u + 8u + (size_t)pkt->length * 8u + 16u;
    if(bits_cap < need_bits){
        return 0;
    }
    size_t bit_pos = 0;
    int rc;
    for(int i = 0; i < 8; i++){
        rc = byte_to_bits(pkt->preamble[i], packet_bits, bits_cap, &bit_pos);
        if(rc != 0){
            return 0;
        }
    }
    rc = byte_to_bits((uint8_t)(pkt->sequence & 0xFF), packet_bits, bits_cap, &bit_pos);
    if(rc != 0){
        return 0;
    }
    rc = byte_to_bits((uint8_t)((pkt->sequence >> 8) & 0xFF), packet_bits, bits_cap, &bit_pos);
    if(rc != 0){
        return 0;
    }
    rc = byte_to_bits(pkt->length, packet_bits, bits_cap, &bit_pos);
    if(rc != 0){
        return 0;
    }
    size_t i = 0;
    while(i < pkt->length){
        rc = byte_to_bits(pkt->payload[i], packet_bits, bits_cap, &bit_pos);
        if(rc != 0){
            return 0;
        }
        i++;
    }
    rc = byte_to_bits((uint8_t)(pkt->crc & 0xFF), packet_bits, bits_cap, &bit_pos);
    if(rc != 0){
        return 0;
    }
    rc = byte_to_bits((uint8_t)((pkt->crc >> 8) & 0xFF), packet_bits, bits_cap, &bit_pos);
    if(rc != 0){
        return 0;
    }
    return bit_pos;
}

int bits_to_float_bpsk_iq(bool *packet_bits, size_t num_bits, float* out_iq, int sps){
    size_t num_complex = num_bits * (size_t)sps;
    memset(out_iq, 0, 2*num_complex*sizeof(float));
    for(size_t i =0; i < num_bits; i++){
        size_t n = i * (size_t)sps;
        size_t w = 2 * n;
        if(packet_bits[i]){
            out_iq[w] = -1.0f;
            out_iq[w + 1] = 0.0f;
        }else{
            out_iq[w] = +1.0f;
            out_iq[w + 1] = 0.0f;
        } 
    }
    return (int)num_complex;
} 

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


static void iq_interleaved_fir_filter(float *x, size_t in_floats, float *h, int num_taps, float *y){
    size_t nin = in_floats/2;
    size_t nout = nin + (size_t)num_taps - 1;
    memset(y, 0, sizeof(float) * 2 * nout);

    for(size_t n = 0; n < nout; n++){
        double accI = 0.0;
        double accQ = 0.0;
        for(int k = 0; k < num_taps; k++){
            long idx = (long)n - k;
            if(idx >= 0 && (size_t)idx < nin){
                size_t ii = (size_t)idx * 2;
                accI += (double)h[k]*(double)x[ii + 0];
                accQ += (double)h[k]*(double)x[ii + 1];
            }
        }
        y[2*n + 0] = (float)accI;
        y[2*n + 1] = (float)accQ;
    }
}

void tx_pulse_shape(float *out_iq, size_t nsamps_unshaped, int sps, float alpha, int span, float **shaped_out_iq, size_t *nsamps_shaped_out){
    int num_taps = span * sps + 1;
    float *taps = (float*)malloc(sizeof(float)*num_taps);
    rrc_generate_taps(alpha, sps, span, taps);
    size_t nin_complex = nsamps_unshaped/2;
    size_t nout_complex = nin_complex + (size_t)num_taps - 1;
    float *y = (float*)malloc(sizeof(float) * 2 * nout_complex);
    iq_interleaved_fir_filter(out_iq, nsamps_unshaped, taps, num_taps, y);
    *shaped_out_iq = y;
    *nsamps_shaped_out = 2 * nout_complex;
    free(taps);
}

void debugPrints(Packet pkt, bool *packet_bits, size_t num_bits, float *out_iq, int num_iq_samps, float *shaped_out_iq, size_t shaped_out_iq_length){
    printf("\n"); 
    printf("Packet Contents \n");
    printf("Preamble: \n");
    for(int i = 0; i < sizeof(pkt.preamble); i++){
        printf("%02X", pkt.preamble[i]);
        printf(" ");
    }
    printf("\n");
    printf("Sequence: \n");
    printf("%d\n", pkt.sequence);
    printf("Length: \n");
    printf("%d\n", pkt.length);
    printf("Payload: \n");
    for(int i = 0; i < sizeof(pkt.payload); i++){
        printf("%02X", pkt.payload[i]);
        printf(" ");
    }
    printf("\n");
    printf("CRC: \n");
    printf("%04X", pkt.crc);
    printf("\n");
    printf("\n"); 
    printf("Size: %d\n",sizeof(pkt));
    printf("\n");
    printf("Packet Bits: \n");
    for(int i = 0; i < num_bits; i++){
        if(packet_bits[i]){
            printf("1");
        }else{
            printf("0");
        }
    }
    printf("\n");
    printf("Number of Bits: %d\n", num_bits);
    printf("\n");
    //printf("IQ Data (unfiltered): \n");
    //for(int i = 0; i < num_iq_samps; i++){
    //    int idx = 2 * i; 
    //    printf("%+.1f, %+.1f | ", out_iq[idx], out_iq[idx + 1]);
    //}
    //printf("\n");
    printf("Number of IQ Samples: %d\n", num_iq_samps);
    printf("\n");
    //printf("Shaped IQ Data: \n");
    //for(size_t i = 0; i < shaped_out_iq_length/2; i++){
    //    int idx = 2 * i; 
    //    printf("%+.1f, %+.1f | ", shaped_out_iq[idx], shaped_out_iq[idx + 1]);
    //}
    //printf("\n");
    printf("Shaped IQ size: %zu\n", shaped_out_iq_length);
}

void write_cf32(FILE *file, float *shaped_out_iq, size_t shaped_out_iq_length){
    size_t wrote = fwrite(shaped_out_iq, sizeof(float), shaped_out_iq_length, file);
    printf("Wrote cf32 IQ to file.\n");
}

int main(void){
    printf("Hello World\n");

    const char* in_path  = "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\hello.txt";
    const char* out_path = "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\goodbye.cfile";

    FILE* inputFile  = OpenFile(in_path,  true);
    FILE* outputFile = OpenFile(out_path, false);
    if (!inputFile || !outputFile) {
        if (inputFile)  fclose(inputFile);
        if (outputFile) fclose(outputFile);
        return 1;
    }

    build_CRC16_table();
    printf("%s\n", crc_table_built ? "CRC Table Built" : "CRC Table Error");

    Packet pkt;
    const uint8_t preamble[8] = {0xAA,0xAA,0xAA,0xAA,0xAA,0xAA,0xAA,0xD5};
    uint16_t sequence = 0;

    // TX params
    const int   sps   = 8;
    const float alpha = 0.5f;
    const int   span  = 11;

    for (;;) {
        // Read up to 253 bytes from the input file
        size_t nread = fread(pkt.payload, 1, sizeof(pkt.payload), inputFile);
        if (nread == 0) break;

        // Fill packet fields for this block
        memcpy(pkt.preamble, preamble, sizeof(pkt.preamble));
        pkt.sequence = sequence++;
        pkt.length   = (uint8_t)nread;

        // Compute CRC over actual payload length
        pkt.crc = CRC16(pkt);

        // ---- Build bitstream (compute exact bit budget) ----
        const size_t need_bits =
            8u*sizeof(pkt.preamble) +   // preamble (8 bytes -> 64 bits)
            16u +                        // sequence (2 bytes)
            8u  +                        // length   (1 byte)
            (size_t)pkt.length * 8u +    // payload
            16u;                         // CRC (2 bytes)

        bool *packet_bits = (bool*)malloc(need_bits * sizeof(bool));
        if (!packet_bits) {
            fprintf(stderr, "malloc failed for packet_bits (%zu bits)\n", need_bits);
            break;
        }

        size_t num_bits = packet_to_bits(&pkt, packet_bits, need_bits);

        // Debug dump (optional)
        debugPrints(pkt, packet_bits, (int)num_bits, NULL, 0, NULL, 0);

        if (num_bits == 0) {
            // Nothing to send for this packet; clean up and continue
            free(packet_bits);
            continue;
        }

        // ---- Map bits to BPSK (I-axis) ----
        size_t out_iq_floats = 2u * (size_t)num_bits * (size_t)sps; // I/Q interleaved
        float *out_iq = (float*)malloc(out_iq_floats * sizeof(float));
        if (!out_iq) {
            fprintf(stderr, "malloc failed for out_iq (%zu floats)\n", out_iq_floats);
            free(packet_bits);
            break;
        }

        int num_iq_samps = bits_to_float_bpsk_iq(packet_bits, num_bits, out_iq, sps);
        printf("Number of Bits: %zu\n", num_bits);
        printf("Number of IQ Samples: %d\n\n", num_iq_samps);

        // ---- RRC pulse shaping ----
        float *shaped_out_iq = NULL;
        size_t shaped_out_iq_length = 0;
        if (num_iq_samps > 0) {
            tx_pulse_shape(out_iq, 2*(size_t)num_iq_samps, sps, alpha, span,
                           &shaped_out_iq, &shaped_out_iq_length);
            printf("Shaped IQ size: %zu\n", shaped_out_iq_length);
            write_cf32(outputFile, shaped_out_iq, shaped_out_iq_length);
        }

        // Cleanup
        free(shaped_out_iq);
        free(out_iq);
        free(packet_bits);
    }

    fclose(inputFile);
    fclose(outputFile);
    return 0;
}
