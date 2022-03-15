#ifndef BPSK_H
#define BPSK_H
#include "IIR_filter.h"
#include "xcorr.h"
#include <arm_math.h>
#include <stdbool.h>
#include <stdint.h>

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
  float32_t loop_gain;
  uint8_t max_error;
  float32_t error;
  uint16_t curr_idx;
} gardnerTimingRecovery_parameters;

typedef struct {
  uint16_t Fc; // carrier frequency
  uint32_t Fs; // sampling frequency
  uint16_t Fb; // bit rate (bps)
  uint16_t samplesPerBit;
  uint16_t FSpan;       // matched filter span in symbols
  uint16_t frameLength; // one byte frame length
  float32_t *matchedFilterCoeffs;
  uint16_t matchedFilterCoeffsLength;
  costasLoop_parameters *costas;
  gardnerTimingRecovery_parameters *gardner;
  int8_t *preambleCode;
  uint16_t preambleCodeLength;
} BPSK_parameters;

void BPSK_getModSamples(BPSK_parameters *params, const uint8_t *data,
                        const uint16_t length, float32_t *outData,
                        const uint16_t outLength);

void BPSK_getOutputSignalWithPrefix(BPSK_parameters *params,
                                    const uint8_t *data,
                                    const uint16_t dataLength,
                                    float32_t *outSignal,
                                    const uint16_t outLength);

void BPSK_getOutputSignalWithPreamble(BPSK_parameters *params,
                                      const uint8_t *data,
                                      const uint16_t dataLength,
                                      float32_t *outSignal,
                                      const uint16_t outLength);

void BPSK_demodulateSignal(BPSK_parameters *params, const float32_t *signal,
                           const uint16_t signalLength, uint8_t *outData,
                           const uint16_t outLength);

void BPSK_demodulateSignal_decimated(BPSK_parameters *params,
                                     const int8_t *signal,
                                     const uint16_t signalLength,
                                     const uint16_t *startIdx,
                                     const uint16_t startNum, uint8_t *outData,
                                     const uint16_t outLength);

void BPSK_findSymbolsStarts_decimated(BPSK_parameters *params, int8_t *signal,
                                      const uint16_t signalLength,
                                      uint16_t *startIdx, uint16_t *foundIdx);

void BPSK_findSymbolsStarts(BPSK_parameters *params, const float32_t *signal,
                            const uint16_t signalLength,
                            const float32_t *preamble,
                            const uint16_t preambleLength, uint16_t *startIdx,
                            uint16_t *foundIdx);

void BPSK_reset(BPSK_parameters *params);

/** Costas Loop for carrier frequency and phase estimation */
void BPSK_carrierRecovery(BPSK_parameters *params, float32_t *signal,
                          const uint16_t signalLength);

void BPSK_timingRecovery(BPSK_parameters *params, float32_t *signal,
                         const uint16_t signalLength, int8_t *output,
                         const uint16_t outputLength);

#endif
