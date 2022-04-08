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

  float32_t a = params->costas->alpha;
  float32_t b = params->costas->beta;
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

static float32_t interpolateLinear(float32_t *signal, uint16_t m_k,
                                   float32_t mu) {
  if (mu < 0) {
    m_k--;
    mu = mu + 1.0f;
  } else if (mu >= 1) {
    m_k++;
    mu = mu - 1.0f;
  }

  return mu * signal[m_k + 1] + (1 - mu) * signal[m_k];
}

void BPSK_timingRecovery(BPSK_parameters *params, float32_t *signal,
                         const uint16_t signalLength, int8_t *output,
                         const uint16_t outputLength) {
  uint16_t midpoint_offset = params->samplesPerBit / 2;
  float32_t Ki = params->gardner->Ki;
  float32_t Kp = params->gardner->Kp;

  uint16_t symbols_num = signalLength / params->samplesPerBit + 1;
  uint16_t start_idx = 0;

  float32_t min_val, max_val;
  uint16_t min_idx, max_idx;
  arm_min_f32(signal, params->samplesPerBit, &min_val, &min_idx);
  arm_max_f32(signal, params->samplesPerBit, &max_val, &max_idx);
  if (min_val < 0) {
    start_idx = min_idx;
  } else {
    start_idx = max_idx;
  }

  uint8_t strobe = 0;
  uint16_t k = 0;         // current symbol index
  uint16_t m_k = 0;       //
  float32_t mu = 0;       // fractional symbol timing offset estimate
  float32_t v_p, v_i;     // proportional and integral output
  float32_t v = 0.0f;     // PI output
  float32_t error = 0.0f; // Error from TED
  float32_t cnt = 1.0f;   // modulo - 1 counter
  float32_t W;

  uint16_t zc_idx;
  float32_t zc_mu;
  float32_t sample_zc;
  float32_t sample;

  float32_t last_sample = signal[start_idx];
  for (uint16_t i = start_idx; i < signalLength; i++) {
    if (strobe) {
      sample = interpolateLinear(signal, m_k, mu);
      zc_idx = m_k - midpoint_offset;
      zc_mu = mu;
      sample_zc = interpolateLinear(signal, zc_idx, zc_mu);
      error = sample_zc * (last_sample - sample);
      last_sample = sample;
      output[k] = sample > 0 ? 1 : -1;
    } else {
      error = 0.0f;
    }

    v_p = Kp * error;
    v_i += (Ki * error);
    v = v_p + v_i;

    W = 1.0f / params->samplesPerBit + v;
    strobe = cnt < W;
    if (strobe) {
      k++;
      m_k = i;
      mu = cnt / W;
    }

    cnt = cnt - W;
    if (cnt > 0) {
      cnt = cnt - (float32_t)((int16_t)cnt);
    } else {
      cnt = 1.0f + cnt - (float32_t)((int16_t)cnt);
    }
  }
}