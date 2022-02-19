/**
 * It's int8_t version of CMSIS DSP's correlate function for float32 type.
 */

#ifndef XCORR_H
#define XCORR_H

#include <stdint.h>

void xcorrelate_int8(int8_t *data1, uint16_t data1_length, int8_t *data2,
                     uint16_t data2_length, int8_t *corr, uint16_t corr_length);

#endif