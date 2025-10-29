#include "api.h"
/* Global storage for CRC16 table + built flag */
uint16_t CRC16_table[256];
bool     crc_table_built = false;
