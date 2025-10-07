#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#define _USE_MATH_DEFINES
#include <math.h>

typedef struct{
    uint8_t preamble[4];
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
    for(int i = 0; i < sizeof(pkt.length); i++){
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
    size_t need_bits = 32u + 16u + 8u + (size_t)pkt->length * 8u + 16u;
    if(bits_cap < need_bits){
        return 0;
    }
    size_t bit_pos = 0;
    int rc;
    for(int i = 0; i < 4; i++){
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

static void rrc_generate_taps(float alpha, int sps, int span, float *taps_out){
    int num_taps = span * sps;
    int mid = num_taps/2;
    double a = (double)alpha;
    double T = 1.0;
    double Ts = T/(double)sps;
    double eps = 1e-12;
    double sum = 0.0;
    for(int n = 0; n < num_taps; n++){
        int k = n - mid;
        double t = k * Ts;
        double x = 4.0 * a * t/T;
        double h;
        if(fabs(t) < eps){
            h = 1.0 + a * (4.0/M_PI - 1.0);
        }else if(a > 0.0 && fabs(fabs(t) - T/(4.0*a)) < (2.0*Ts)){
            double theta = M_PI/(4.0*a);
            h = (a/sqrt(2.0)) * ((1.0 + 2.0/M_PI)*sin(theta) + (1.0 - 2.0/M_PI)*cos(theta));
        }else{
            double num1 = cos((1.0 + a)*M_PI*t/T);
            double num2 = (T/(4.0*a*t)) * sin((1.0 - a)*M_PI*t/T);
            double den1 = M_PI*t/T;
            double den2 = 1.0 - x*x;
            h = (4.0*a) * (num1 + num2)/(den1 * den2); 
        }

        taps_out[n] = (float)h;
        sum += h*h; 
    }
    double scale = (sum > 0.0) ? (1.0/sqrt(sum)) : 1.0;
        for(int n = 0; n < num_taps; n++){
            taps_out[n] = (float)(taps_out[n]*scale);
        } 
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
    (void)sps;
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
    printf("IQ Data (unfiltered): \n");
    for(int i = 0; i < num_iq_samps; i++){
        int idx = 2 * i; 
        printf("%+.1f, %+.1f | ", out_iq[idx], out_iq[idx + 1]);
    }
    printf("\n");
    printf("Number of IQ Samples: %d\n", num_iq_samps);
    printf("\n");
    printf("Shaped IQ Data: \n");
    for(size_t i = 0; i < shaped_out_iq_length/2; i++){
        int idx = 2 * i; 
        printf("%+.1f, %+.1f | ", shaped_out_iq[idx], shaped_out_iq[idx + 1]);
    }
    printf("\n");
    printf("Shaped IQ size: %zu\n", shaped_out_iq_length);
}

void write_cf32(FILE *file, float *shaped_out_iq, size_t shaped_out_iq_length){
    size_t wrote = fwrite(shaped_out_iq, sizeof(float), shaped_out_iq_length, file);
    printf("Wrote cf32 IQ to file.\n");
}

int main(){
    printf("Hello World\n");
    const char* in_path = "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\hello.txt";
    const char* out_path = "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\goodbye.cfile";
    FILE* inputFile = OpenFile(in_path, true);
    FILE* outputFile = OpenFile(out_path, false); 
    Packet pkt;
    uint8_t preamble[4] = {0xAA, 0xAA, 0xAA, 0xAA};
    uint16_t sequence = 0;
    build_CRC16_table();
    if(crc_table_built){
        printf("CRC Table Built\n");
    }else{
        printf("CRC Table Error\n");
    }

    while ((fread(pkt.payload, 1, sizeof(pkt.payload), inputFile)) > 0){
        memcpy(pkt.preamble, preamble, sizeof(pkt.preamble));
        pkt.sequence = sequence;
        pkt.length = sizeof(pkt.payload);
        sequence++;
        pkt.crc = CRC16(pkt);
        bool packet_bits[2096];
        size_t num_bits = packet_to_bits(&pkt, packet_bits, 2096);
        int sps = 8;
        float out_iq[2*num_bits*sps];
        int num_iq_samps = bits_to_float_bpsk_iq(packet_bits, num_bits, out_iq, sps);

        float alpha = 0.22f;
        int span = 11;

        float *shaped_out_iq = NULL;
        size_t shaped_out_iq_length = 0;

        tx_pulse_shape(out_iq, 2*(size_t)num_iq_samps, sps, alpha, span, &shaped_out_iq, &shaped_out_iq_length);


        debugPrints(pkt, packet_bits, num_bits, out_iq, num_iq_samps, shaped_out_iq, shaped_out_iq_length);
        
        write_cf32(outputFile, shaped_out_iq, shaped_out_iq_length);

        free(shaped_out_iq);
    }

    fclose(inputFile);
    fclose(outputFile);
    return 0;
}
