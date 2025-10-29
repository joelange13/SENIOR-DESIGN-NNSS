#include "api.h"
/*
 * build_CRC16_table()
 *  Input : none
 *  Output: fills global CRC16_table[256], sets crc_table_built=true
 */
void build_CRC16_table(void){
    uint16_t crc_polynomial = 0x1021;
    for(uint16_t i = 0; i < 256; ++i){
        uint16_t crc = (uint16_t)(i << 8);
        for(int b = 0; b < 8; ++b){
            if((crc & 0x8000) != 0) crc = (uint16_t)((crc << 1) ^ crc_polynomial);
            else                    crc = (uint16_t)(crc << 1);
        }
        CRC16_table[i] = crc;
    }
    crc_table_built = true;
}
