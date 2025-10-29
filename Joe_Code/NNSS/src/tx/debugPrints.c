#include "api.h"
/*
 * debugPrints(pkt, packet_bits, num_bits, out_iq, num_iq_samps, shaped_out_iq, shaped_len)
 *  Input : packet + optional buffers (IQ can be NULL)
 *  Output: prints human-readable summary to stdout
 */
void debugPrints(Packet pkt, bool *packet_bits, size_t num_bits,
                 float *out_iq, int num_iq_samps, float *shaped_out_iq, size_t shaped_out_iq_length){
    (void)out_iq;
    (void)shaped_out_iq;
    (void)num_iq_samps;
    (void)shaped_out_iq_length;
    printf("\nPacket Contents\nPreamble:\n");
    for(size_t i = 0; i < sizeof(pkt.preamble); i++) printf("%02X ", pkt.preamble[i]);
    printf("\nSequence:\n%d\nLength:\n%d\nPayload:\n", pkt.sequence, pkt.length);
    for(size_t i = 0; i < sizeof(pkt.payload); i++) printf("%02X ", pkt.payload[i]);
    printf("\nCRC:\n%04X\n\n", pkt.crc);
    printf("Size: %d\n\n", (int)sizeof(pkt));

    printf("Packet Bits:\n");
    for(size_t i = 0; i < num_bits; i++) printf("%c", packet_bits[i] ? '1' : '0');
    printf("\nNumber of Bits: %d\n\n", (int)num_bits);

    //printf("Number of IQ Samples: %d\n", num_iq_samps);
    //printf("Shaped IQ size: %zu\n", shaped_out_iq_length);
}
