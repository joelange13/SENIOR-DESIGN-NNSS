#include "api.h"
/*
 * packet_to_bits(pkt, packet_bits, bits_cap)
 *  Input : pkt (preamble, sequence, length, payload, crc),
 *          packet_bits buffer + capacity in bits
 *  Output: number of bits written; 0 on error/overflow
 */
size_t packet_to_bits(const Packet *pkt, bool *packet_bits, size_t bits_cap ){
    size_t need_bits = 64u + 16u + 8u + (size_t)pkt->length * 8u + 16u;
    if(bits_cap < need_bits) return 0;

    size_t bit_pos = 0;
    int rc;
    for(int i = 0; i < 8; i++){
        rc = byte_to_bits(pkt->preamble[i], packet_bits, bits_cap, &bit_pos);
        if(rc != 0) return 0;
    }
    rc = byte_to_bits((uint8_t)(pkt->sequence & 0xFF), packet_bits, bits_cap, &bit_pos); if(rc) return 0;
    rc = byte_to_bits((uint8_t)((pkt->sequence >> 8) & 0xFF), packet_bits, bits_cap, &bit_pos); if(rc) return 0;
    rc = byte_to_bits(pkt->length, packet_bits, bits_cap, &bit_pos); if(rc) return 0;

    for(size_t i = 0; i < pkt->length; i++){
        rc = byte_to_bits(pkt->payload[i], packet_bits, bits_cap, &bit_pos);
        if(rc != 0) return 0;
    }
    rc = byte_to_bits((uint8_t)(pkt->crc & 0xFF), packet_bits, bits_cap, &bit_pos); if(rc) return 0;
    rc = byte_to_bits((uint8_t)((pkt->crc >> 8) & 0xFF), packet_bits, bits_cap, &bit_pos); if(rc) return 0;

    return bit_pos;
}
