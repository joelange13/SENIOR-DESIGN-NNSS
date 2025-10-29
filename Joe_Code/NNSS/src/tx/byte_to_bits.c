#include "api.h"
/*
 * byte_to_bits(byte, packet_bits, bits_cap, bit_pos)
 *  Input : byte (MSB-first), packet_bits buffer, capacity, in/out bit_pos
 *  Output: 0 on success, -1 if capacity exceeded (no write past end)
 */
int byte_to_bits(uint8_t byte, bool* packet_bits, size_t bits_cap, size_t* bit_pos){
    for (int b = 0; b < 8; ++b){
        if(*bit_pos >= bits_cap) return -1;
        packet_bits[(*bit_pos)++] = ((byte & (uint8_t)(0x80u >> b)) != 0);
    }
    return 0;
}
