#ifndef BPSK_H
#define BPSK_H
#include <stdint.h>
#include <arm_math.h>

typedef struct BPSK_parameters {
    uint16_t Fc; // carrier frequency
    uint16_t Fs; // sampling frequency
    uint16_t Fb; // bit rate (bps)
    float32_t* firCoeffs;
    uint16_t firCoeffsLength;
} BPSK_parameters;

void BPSK_getModSamples(BPSK_parameters* params, uint8_t* data, uint16_t length, float32_t* outData, uint16_t outLength);

void BPSK_getOutputSignal(BPSK_parameters* params, uint8_t* data, uint16_t dataLength, float32_t* outSignal, uint16_t outLength);
#endif
