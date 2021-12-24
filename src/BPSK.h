#ifndef BPSK_H
#define BPSK_H
#include "IIR_filter.h"
#include <arm_math.h>
#include <stdbool.h>
#include <stdint.h>

#define COSTAS_LOOP_MAX_ERROR 0.001f;

typedef struct {
  float32_t alpha;
  float32_t beta;
  IIR_filter *LP_filterI;
  IIR_filter *LP_filterQ;
  float32_t omega;
  float32_t error;
  float32_t phase;
} costasLoop_parameters;

typedef struct {
  uint16_t Fc; // carrier frequency
  uint32_t Fs; // sampling frequency
  uint16_t Fb; // bit rate (bps)
  uint16_t samplesPerBit;
  uint16_t FSpan;        // filter span in symbols
  uint16_t prefixLength; // prefix length (for synchro purposes)
  uint16_t frameLength;  // one data frame length
  float32_t *matchedFilterCoeffs;
  uint16_t matchedFilterCoeffsLength;
  costasLoop_parameters *costas;
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
                          uint16_t signalLength, uint16_t *startIdx,
                          uint16_t *foundIdx);

void BPSK_syncInputSignal_(BPSK_parameters *params, float32_t *signal,
                           uint16_t signalLength, uint16_t *startIdx,
                           uint16_t *foundIdx);

void BPSK_init(BPSK_parameters *params);

/** Costas Loop for carrier frequency and phase estimation */
void BPSK_syncSignalCarrier(BPSK_parameters *params, float32_t *signal,
                            uint16_t signalLength);

#endif
