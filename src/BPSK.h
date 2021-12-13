#ifndef BPSK_H
#define BPSK_H
#include <arm_math.h>
#include <stdbool.h>
#include <stdint.h>


typedef struct BPSK_parameters {
  uint16_t Fc; // carrier frequency
  uint32_t Fs; // sampling frequency
  uint16_t Fb; // bit rate (bps)
  uint16_t samplesPerBit;
  uint16_t FSpan;        // filter span in symbols
  uint16_t prefixLength; // prefix length (for synchro purposes)
  uint16_t frameLength;  // one data frame length
  float32_t *firCoeffs;
  uint16_t firCoeffsLength;
  bool differential;
} BPSK_parameters;

void BPSK_getModSamples(BPSK_parameters *params, uint8_t *data, uint16_t length,
                        float32_t *outData, uint16_t outLength);

void BPSK_getOutputSignal(BPSK_parameters *params, uint8_t *data,
                          uint16_t dataLength, float32_t *outSignal,
                          uint16_t outLength);

void BPSK_demodulateSignal(BPSK_parameters *params, float32_t *signal,
                           uint16_t signalLength, uint8_t *outData,
                           uint16_t outLength);

void BPSK_syncInputSignal(BPSK_parameters *params, float32_t *signal,
                          uint16_t signalLength, uint32_t *startIdx,
                          uint16_t *foundIdx);

#endif
