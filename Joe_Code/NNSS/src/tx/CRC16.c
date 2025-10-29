#include "api.h"
/*
 * CRC16(pkt)
 *  Input : Packet pkt (uses pkt.payload[0..length-1])
 *  Output: 16-bit CRC (requires build_CRC16_table() to be called once)
 */
uint16_t CRC16(Packet pkt){
    uint16_t crc = 0xFFFF;
    for(int i = 0; i < pkt.length; i++){
        uint8_t idx = (uint8_t)(((crc >> 8) & 0xFF) ^ pkt.payload[i]);
        crc = (uint16_t)((crc << 8) ^ CRC16_table[idx]);
    }
    return crc;
}
