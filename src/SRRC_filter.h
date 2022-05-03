#ifndef SRRC_FILTER_H
#define SRRC_FILTER_H
#include "arm_math.h"
#include <stdint.h>

void SRRC_getFIRCoeffs(uint8_t symbolSpan, uint8_t SPB, float32_t beta,
                       float32_t *firCoeffs, uint16_t length,
                       float32_t linearGain);

#endif
