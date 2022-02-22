#include "BPSK.h"
#include "uart.h"

// turn bytes into -1 (corresponding to 0) and 1 (corresponding to 1) array
static void BPSK_generateModData(const uint8_t *data, const uint16_t length,
                                 int8_t *output) {
  for (uint16_t i = 0; i < length; i++) {
    for (uint8_t j = 0; j < 8; j++) {
      *(output + i * 8 + j) = (data[i] & (0x80 >> j)) ? -1 : 1;
    }
  }
}

static void BPSK_generateDifferentialModData(const uint8_t *data,
                                             const uint16_t length,
                                             int8_t *output) {}

// generate modulation samples (with oversampling)
void BPSK_getModSamples(BPSK_parameters *params, const uint8_t *data,
                        const uint16_t length, float32_t *outData,
                        const uint16_t outLength) {

  assert(params->samplesPerBit > 0);
  assert(outLength == params->samplesPerBit * length * 8);

  int8_t modData[length * 8];

  if (params->differential)
    BPSK_generateDifferentialModData(data, length, modData);
  else
    BPSK_generateModData(data, length, modData);

  for (uint16_t i = 0; i < outLength; i = i + params->samplesPerBit) {
    for (uint16_t j = 0; j < params->samplesPerBit; j++) {
      *(outData + i + j) = modData[i / params->samplesPerBit];
    }
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

void BPSK_setPreamble(BPSK_parameters *params, const int16_t *code,
                      const uint16_t length, float32_t *preamble,
                      const uint16_t preambleLength) {

  assert(preambleLength == params->samplesPerBit * length);

  for (uint16_t i = 0; i < preambleLength; i = i + params->samplesPerBit) {
    for (uint16_t j = 0; j < params->samplesPerBit; j++) {
      *(preamble + i + j) =
          (code[i / params->samplesPerBit] > 0) ? 1.0f : -1.0f;
    }
  }

  params->preamble = preamble;
  params->preambleLength = preambleLength;
}

// generate output signal
void BPSK_getOutputSignalWithPrefix(BPSK_parameters *params,
                                    const uint8_t *data,
                                    const uint16_t dataLength,
                                    float32_t *outSignal,
                                    const uint16_t outLength) {

  assert(params->samplesPerBit > 0);
  assert(outLength ==
         params->samplesPerBit * dataLength * 8 + params->prefixLength);

  BPSK_getModSamples(params, data, dataLength, outSignal + params->prefixLength,
                     outLength - params->prefixLength);

  for (uint16_t i = 0; i < params->prefixLength; i++) {
    outSignal[i] = outSignal[outLength - params->prefixLength + i];
  }
}

void BPSK_getOutputSignalWithPreamble(BPSK_parameters *params,
                                      const uint8_t *data,
                                      const uint16_t dataLength,
                                      float32_t *outSignal,
                                      const uint16_t outLength) {
  assert(params->samplesPerBit > 0);
  assert(outSignal ==
         params->samplesPerBit * dataLength * 8 + params->preambleLength);

  BPSK_getModSamples(params, data, dataLength,
                     outSignal + params->preambleLength,
                     outLength - params->preambleLength);

  for (uint16_t i = 0; i < params->preambleLength; i++) {
    outSignal[i] = params->preamble[i];
  }
}

void BPSK_demodulateSignal(BPSK_parameters *params, const float32_t *signal,
                           const uint16_t signalLength, uint8_t *outData,
                           uint16_t outLength) {

  uint16_t k = 0;

  if (!params->differential) {
    for (uint16_t i = params->samplesPerBit / 2;
         i < signalLength - params->samplesPerBit / 2;
         i = i + params->samplesPerBit * 8) {
      outData[k] = 0;
      for (uint16_t j = 0; j < 8 * params->samplesPerBit;
           j = j + params->samplesPerBit) {
        if (signal[i + j] < 0)
          outData[k] += (1 << (7 - j / params->samplesPerBit));
      }
      ++k;
    }
  } else {
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
  prevStart = maximas[0];
  if (prevStart >= 13) {
    *(startIdx) = prevStart + params->preambleCodeLength - 13;
    *(foundIdx) = *(foundIdx) + 1;
  }

  for (uint16_t i = 1; i < maximasCount; i++) {
    nextStart = *(maximas + i);
    uint16_t distance = nextStart - prevStart;
    if (distance == 13) {
      *(startIdx + (*foundIdx)) = prevStart + params->preambleCodeLength;
      *foundIdx = *foundIdx + 1;
      prevStart = nextStart;
    } else if (distance >= 13 - 2 && distance <= 13 + 2) {
      *(startIdx + (*foundIdx)) = nextStart - 13 + params->preambleCodeLength;
      prevStart = nextStart;
      *foundIdx = *foundIdx + 1;
    } else if (distance > 13 + 2 && (distance % 13) <= 3) {
      prevStart = nextStart;
      uint8_t num = distance / 13U;
      if (num == 1) {
        *(startIdx + *foundIdx) = prevStart + params->preambleCodeLength;
        *foundIdx = *foundIdx + 1;
      }
      while (num) {
        *(startIdx + (*foundIdx) + --num) =
            nextStart - 13 + params->preambleCodeLength;
        nextStart -= 13;
      }
      *(foundIdx) = *(foundIdx) + num;
    }
  }
  nextStart = maximas[maximasCount - 1];

  if (nextStart - (*(startIdx + *foundIdx - 1) - params->preambleCodeLength) >=
          13 - 2 &&
      nextStart - (*(startIdx + *foundIdx - 1) - params->preambleCodeLength) <=
          13 + 2) {
    *(startIdx + *foundIdx) = nextStart + params->preambleCodeLength;
    *foundIdx = *foundIdx + 1;
  }
}

void BPSK_findSymbolsStarts(BPSK_parameters *params, const float32_t *signal,
                            const uint16_t signalLength, uint16_t *startIdx,
                            uint16_t *foundIdx) {
  bool locked = false; // If any symbol detected
  float32_t corr[2 * signalLength - 1];
  arm_correlate_f32(signal, signalLength, params->preamble,
                    params->preambleLength, corr);
  uint16_t i = signalLength;

  float32_t absMaxVal;
  uint8_t dumm;
  arm_max_f32(corr + i, signalLength, &absMaxVal, &dumm);

  float32_t maxVal;
  uint16_t idx;

  while (i < 2 * signalLength - 1 - params->frameLength) {
    if (!locked) {
      arm_max_f32(corr + i, params->frameLength + params->preambleLength,
                  &maxVal, &idx);
      if (maxVal > 0.7 * absMaxVal) {
        *(startIdx + *foundIdx) =
            i - signalLength + idx + params->preambleLength;
        ++(*foundIdx);
        i += idx + params->preambleLength + params->frameLength;
        locked = true;
      } else {
        i += params->frameLength + params->preambleLength;
      }
    } else {
      arm_max_f32(corr + i - params->samplesPerBit, 2 * params->samplesPerBit,
                  &maxVal, &idx);
      if (maxVal > 0.7 * absMaxVal) {
        *(startIdx + *foundIdx) = i - signalLength - params->samplesPerBit +
                                  idx + params->preambleLength;
        i += idx - params->samplesPerBit + params->preambleLength +
             params->frameLength;
        ++(*foundIdx);
      } else {
        locked = false;
        i += params->samplesPerBit;
      }
    }
  }
}

void BPSK_init(BPSK_parameters *params) {
  // Initialize costas loop
  params->costas->error = 0.0f;
  params->costas->omega = 2.0f * PI * params->Fc / params->Fs;
  params->costas->phase = 0.0f;
  params->gardner->error = 0.0f;
  params->gardner->curr_idx = 0;
}

void BPSK_syncSignalCarrier(BPSK_parameters *params, float32_t *signal,
                            const uint16_t signalLength) {
  float32_t errorTot = 0.0f;
  float32_t si, sq, sim, sqm;
  for (uint16_t i = 0; i < signalLength; ++i) {
    params->costas->phase += params->costas->omega;
    params->costas->phase += params->costas->alpha * params->costas->error;

    params->costas->omega += params->costas->beta * params->costas->error;

    float32_t freq = params->costas->omega * params->Fs / (2 * PI);

    if (params->costas->phase > 2 * PI) {
      params->costas->phase -= 2 * PI;
    }

    si = arm_cos_f32(params->costas->phase);
    sq = -arm_sin_f32(params->costas->phase);

    sim = si * signal[i];
    sqm = sq * signal[i];

    sim = IIR_filter_step(params->costas->LP_filterI, sim);
    sqm = IIR_filter_step(params->costas->LP_filterQ, sqm);

    params->costas->error = sim * sqm;

    signal[i] = sim;

    errorTot += params->costas->error;
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