#include "BPSK.h"
#include "uart.h"

#define MICROSECONDS 1000000.0f
// turn bytes into -1 (corresponding to 0) and 1 (corresponding to 1) array
void BPSK_generateModData(const uint8_t *data, const uint16_t length,
                          int8_t *output) {
  for (uint16_t i = 0; i < length; i++) {
    for (uint8_t j = 0; j < 8U; j++) {
      *(output + i * 8 + j) = (data[i] & (0x80 >> j)) ? -1 : 1;
    }
  }
}

// generate modulation samples (with oversampling)
void BPSK_getModSamples(BPSK_parameters *params, const uint8_t *data,
                        const uint16_t length, float32_t *outData,
                        const uint16_t outLength) {

  assert(params->samplesPerBit > 0);
  assert(outLength == params->samplesPerBit * length * 8);

  int8_t modData[length * 8];

  BPSK_generateModData(data, length, modData);

  for (uint16_t i = 0; i < outLength; i += params->samplesPerBit) {
    *(outData + i) = modData[i];
  }

  if (params->matchedFilterCoeffsLength) {
    float32_t tempData[outLength + params->samplesPerBit * params->FSpan];

    arm_conv_f32(outData, outLength, params->matchedFilterCoeffs,
                 params->matchedFilterCoeffsLength, tempData);
    float32_t maxVal;
    void *x;
    arm_max_f32(tempData, outLength, &maxVal, x);
    uint32_t k = 0;
    for (uint16_t i = params->FSpan * params->samplesPerBit / 2;
         i < outLength + params->FSpan * params->samplesPerBit / 2; i++) {
      outData[k++] = tempData[i] / maxVal;
    }
  }
}

// generate output signal
void BPSK_getOutputSignalWithPrefix(BPSK_parameters *params,
                                    const uint8_t *data,
                                    const uint16_t dataLength,
                                    const uint16_t prefixLength,
                                    float32_t *outSignal,
                                    const uint16_t outLength) {

  assert(params->samplesPerBit > 0);
  assert(outLength == params->samplesPerBit * dataLength * 8 + prefixLength);

  BPSK_getModSamples(params, data, dataLength, outSignal + prefixLength,
                     outLength - prefixLength);

  for (uint16_t i = 0; i < prefixLength; i++) {
    outSignal[i] = outSignal[outLength - prefixLength + i];
  }
}

void BPSK_getOutputSignalWithPreamble(BPSK_parameters *params,
                                      const uint8_t *data,
                                      const uint16_t dataLength,
                                      float32_t *outSignal,
                                      const uint16_t outLength) {
  assert(params->samplesPerBit > 0);
  assert(outLength ==
         params->samplesPerBit * (dataLength * 8 + params->preambleCodeLength));

  int8_t modData[dataLength * 8];

  BPSK_generateModData(data, dataLength, modData);

  uint16_t i;
  uint16_t k = 0;
  for (i = 0; i < params->preambleCodeLength * params->samplesPerBit;
       i += params->samplesPerBit) {
    *(outSignal + i) = (float32_t)params->preambleCode[k++];
  }

  k = 0;
  for (; i < outLength; i += params->samplesPerBit) {
    *(outSignal + i) = (float32_t)modData[k++];
  }

  if (params->matchedFilterCoeffsLength) {
    float32_t tempData[outLength + params->samplesPerBit * params->FSpan];

    arm_conv_f32(outSignal, outLength, params->matchedFilterCoeffs,
                 params->matchedFilterCoeffsLength, tempData);
    float32_t maxVal;
    void *x;
    arm_max_f32(tempData, outLength, &maxVal, x);
    k = 0;
    for (i = params->FSpan * params->samplesPerBit / 2;
         i < outLength + params->FSpan * params->samplesPerBit / 2; i++) {
      outSignal[k++] = tempData[i] / maxVal;
    }
  }
}

void BPSK_demodulateSignal_decimated(BPSK_parameters *params,
                                     const int8_t *signal,
                                     const uint16_t signalLength,
                                     const uint16_t *startIdx,
                                     const uint16_t startNum, uint8_t *outData,
                                     const uint16_t outLength) {
  assert(startNum <= outLength);

  for (uint16_t i = 0; i < startNum; ++i) {
    *(outData + i) = 0;
    if ((*startIdx + i) > signalLength - params->frameLength)
      break;
    for (uint16_t j = 1; j <= params->frameLength; ++j) {
      uint8_t sym = *(signal + *(startIdx + i) + j) > 0 ? 0 : 1;
      *(outData + i) += sym << (params->frameLength - j);
    }
  }
}

void BPSK_demodulateSignal(BPSK_parameters *params, const float32_t *signal,
                           const uint16_t signalLength, uint8_t *outData,
                           uint16_t outLength) {

  uint16_t k = 0;

  for (uint16_t i = params->samplesPerBit / 2;
       i < signalLength - params->samplesPerBit / 2;
       i = i + params->samplesPerBit * 8) {
    outData[k] = 0;
    if (i == params->samplesPerBit / 2) {
      if (signal[i] < 0)
        outData[k] += 1;
      for (uint16_t j = params->samplesPerBit; j < 8 * params->samplesPerBit;
           j = j + params->samplesPerBit) {
        if (signal[i + j] * signal[i + j - params->samplesPerBit] < 0)
          outData[k] += (1 << (7 - j / params->samplesPerBit));
      }
    } else {
      for (uint16_t j = 0; j < 8 * params->samplesPerBit;
           j = j + params->samplesPerBit) {
        if (signal[i + j] * signal[i + j - params->samplesPerBit] < 0)
          outData[k] += (1 << (7 - j / params->samplesPerBit));
      }
    }
    ++k;
  }
}

void BPSK_findSymbolsStarts_decimated(BPSK_parameters *params, int8_t *signal,
                                      const uint16_t signalLength,
                                      uint16_t *startIdx, uint16_t *foundIdx) {
  uint16_t corrLength = signalLength > params->preambleCodeLength
                            ? (2 * signalLength - 1)
                            : (2 * params->preambleCodeLength - 1);

  int8_t corr[corrLength];
  xcorrelate_int8(signal, signalLength, params->preambleCode,
                  params->preambleCodeLength, corr, corrLength);

  uint16_t plusMax[signalLength / 13 * 2]; // TODO: It shouldn't be hardcoded
  uint16_t minusMax[signalLength / 13 * 2];
  uint16_t plusCount = 0;
  uint16_t minusCount = 0;

  for (uint16_t i = signalLength + 1; i < corrLength; ++i) {
    if (corr[i] == params->preambleCodeLength) {
      plusMax[plusCount++] = i - signalLength;
    } else if (corr[i] == -params->preambleCodeLength) {
      minusMax[minusCount++] = i - signalLength;
    }
  }

  uint16_t *maximas;
  uint16_t maximasCount;
  if (minusCount > plusCount) {
    for (uint16_t i = 0; i < signalLength; i++) {
      *(signal + i) = -(*(signal + i));
    }
    maximas = minusMax;
    maximasCount = minusCount;
  } else {
    maximas = plusMax;
    maximasCount = plusCount;
  }

  uint16_t prevStart, nextStart;
  prevStart = *(maximas);
  nextStart = prevStart;
  if (prevStart >= 13U - params->preambleCodeLength) {
    uint16_t num = prevStart / 13U;
    if (num == 0)
      ++num;
    while (num) {
      nextStart -= 13U;
      *(startIdx + *foundIdx + num - 1U) =
          nextStart + params->preambleCodeLength;
      --num;
    }
    *(foundIdx) = *(foundIdx) + prevStart / 13U;
  }

  for (uint16_t i = 1; i < maximasCount; i++) {
    nextStart = *(maximas + i);
    uint16_t distance = nextStart - prevStart;
    if (distance == 13U) {
      *(startIdx + (*foundIdx)) = prevStart + params->preambleCodeLength;
      *foundIdx = *foundIdx + 1U;
      prevStart = nextStart;
    } else if (distance >= 13U - 2U && distance <= 13U + 2U) {
      *(startIdx + (*foundIdx)) = nextStart - 13U + params->preambleCodeLength;
      prevStart = nextStart;
      *foundIdx = *foundIdx + 1U;
    } else if (distance > 13U + 2U &&
               ((distance % 13U) <= 1U || (distance % 13U) >= 12U)) {
      prevStart = nextStart;
      uint16_t num = distance / 13U;
      while (num) {
        nextStart -= 13U;
        *(startIdx + *foundIdx + num - 1U) =
            nextStart + params->preambleCodeLength;
        --num;
      }
      *(foundIdx) = *(foundIdx) + distance / 13U;
    } else {
      uint16_t distance_ = *(maximas + i + 1) - nextStart;
      if (distance_ % 13 <= 1U || distance_ % 13 >= 12U) {
        prevStart = nextStart;
      }
    }
  }
  nextStart = *(maximas + maximasCount - 1U);

  if ((nextStart - (*(startIdx + *foundIdx - 1U) -
                    params->preambleCodeLength)) >= (13U - 2U) &&
      (nextStart - (*(startIdx + *foundIdx - 1U) -
                    params->preambleCodeLength)) <= (13U + 2U)) {
    *(startIdx + *foundIdx) = nextStart + params->preambleCodeLength;
    *foundIdx = *foundIdx + 1U;
  }
}

void BPSK_findSymbolsStarts(BPSK_parameters *params, const float32_t *signal,
                            const uint16_t signalLength,
                            const float32_t *preamble,
                            const uint16_t preambleLength, uint16_t *startIdx,
                            uint16_t *foundIdx) {
  bool locked = false; // If any symbol detected
  float32_t corr[2 * signalLength - 1];
  arm_correlate_f32(signal, signalLength, preamble, preambleLength, corr);
  uint16_t i = signalLength;

  float32_t absMaxVal;
  uint8_t dumm;
  arm_max_f32(corr + i, signalLength, &absMaxVal, &dumm);

  float32_t maxVal;
  uint16_t idx;

  while (i < 2 * signalLength - 1 - params->frameLength) {
    if (!locked) {
      arm_max_f32(corr + i, params->frameLength + preambleLength, &maxVal,
                  &idx);
      if (maxVal > 0.6 * absMaxVal) {
        *(startIdx + *foundIdx) = i - signalLength + idx + preambleLength;
        ++(*foundIdx);
        i += idx + preambleLength + params->frameLength;
        locked = true;
      } else {
        i += params->frameLength + preambleLength;
      }
    } else {
      arm_max_f32(corr + i - params->samplesPerBit, 2 * params->samplesPerBit,
                  &maxVal, &idx);
      if (maxVal > 0.6 * absMaxVal) {
        *(startIdx + *foundIdx) =
            i - signalLength - params->samplesPerBit + idx + preambleLength;
        i += idx - params->samplesPerBit + preambleLength + params->frameLength;
        ++(*foundIdx);
      } else {
        locked = false;
        i += params->samplesPerBit;
      }
    }
  }
}

void BPSK_reset(BPSK_parameters *params) {
  // Initialize costas loop
  params->costas->error_int = 0.0f;
  params->costas->lock = 0.0f;
  params->costas->ask = 0.0f;
  params->costas->period = 1000000.0f / (float32_t)params->Fc;
  params->costas->f = -(MICROSECONDS / (float32_t)params->Fs) * 2.0f * M_PI /
                      (MICROSECONDS / (float32_t)params->Fc);
  params->costas->phase = 0;
  FIR_filter_reset(params->costas->LP_filterI);
  FIR_filter_reset(params->costas->LP_filterQ);

  params->gardner->curr_idx = 0;
}

void BPSK_carrierRecovery(BPSK_parameters *params, float32_t *signal,
                          const uint16_t signalLength) {
  float32_t REF_PERIOD = MICROSECONDS / (float32_t)params->Fc;
  float32_t SAMPLE_PERIOD = MICROSECONDS / (float32_t)params->Fs;

  float32_t a = params->costas->Kp;
  float32_t b = params->costas->Ki;
  float32_t intgralf = 1.2f * a / REF_PERIOD;
  float32_t I, Q, S;
  float32_t S_I, S_Q;
  float32_t error;
  float32_t locked, all;

    params->costas->omega += params->costas->beta * params->costas->error;

  for (uint16_t i = 0; i < signalLength; i++) {
    S = *(signal + i);

    params->costas->f += SAMPLE_PERIOD * 2 * M_PI / params->costas->period;

    if (params->costas->f > 2 * M_PI)
      params->costas->f -= 2 * M_PI;

    I = arm_cos_f32(params->costas->f + params->costas->phase);
    Q = -arm_sin_f32(params->costas->f + params->costas->phase);

    S_I = FIR_filter_step(params->costas->LP_filterI, S * I);
    S_Q = FIR_filter_step(params->costas->LP_filterQ, S * Q);

    error = (S_I > 0.0f) ? S_Q : -S_Q;
    *(signal + i) = S_I;
    params->costas->error_int += error * SAMPLE_PERIOD;
    error += params->costas->error_int * intgralf;

    params->costas->phase += a * error;
    params->costas->period -= b * error;

    params->costas->lock = (S_I * S_I) - (S_Q * S_Q);
    params->costas->ask = (S_I * S_I) + (S_Q * S_Q);
    ++all;
    if (params->costas->lock > 0.011f)
      ++locked;
  }

  // TODO: Checking if locked?
}

void BPSK_timingRecovery(BPSK_parameters *params, float32_t *signal,
                         const uint16_t signalLength, int8_t *output,
                         const uint16_t outputLength) {
  uint16_t idx = params->gardner->curr_idx;
  uint16_t si = 0;
  uint16_t left_idx, right_idx, mid_idx;

  float32_t gain = params->gardner->loop_gain;

  uint16_t spb = params->samplesPerBit;
  uint16_t half_spb = spb / 2;
  uint16_t quart_spb = spb / 4;

  while (idx < signalLength) {
    right_idx = idx + half_spb * 3;
    left_idx = idx + half_spb;
    mid_idx = idx + spb;
    *(output + si++) = *(signal + left_idx) > 0 ? -1 : 1;
    float32_t error =
        (*(signal + left_idx) - *(signal + right_idx)) * *(signal + mid_idx);

    if (error > params->gardner->max_error)
      error = params->gardner->max_error;
    else if (error < -params->gardner->max_error)
      error = -params->gardner->max_error;

    uint8_t correction_offset = half_spb * error * gain;
    idx += spb + correction_offset;
  }

  params->gardner->curr_idx = idx % signalLength;
}