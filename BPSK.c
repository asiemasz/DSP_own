#include "BPSK.h"

static void generateModData(BPSK_parameters* param, uint8_t* data, uint16_t length, int8_t* output) {
    for(uint16_t i = 0; i < length; i++) {
        for(uint8_t j = 0; j < 8; j++) {
            *(output + i*8 + j) = (data[i] & (0x80 >> j)) ? 1 : -1;
        }
    }
}
