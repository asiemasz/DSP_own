#include "BPSK.h"

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
  assert(outLength == dataLength * params->samplesPerBit *
                          (params->frameLength + params->preambleCodeLength));

  int8_t modData[dataLength * 8];

  BPSK_generateModData(data, dataLength, modData);

  uint16_t i;
  uint16_t k = 0;

  for (uint16_t data_num = 0; data_num < dataLength; data_num++) {
    uint16_t preamble_iter = 0;
    for (i = 0; i < params->preambleCodeLength * params->samplesPerBit;
         i += params->samplesPerBit) {
      *(outSignal +
        data_num * params->samplesPerBit *
            (params->frameLength + params->preambleCodeLength) +
        i) = (float32_t)params->preambleCode[preamble_iter++];
    }

    for (; i < (params->preambleCodeLength + params->frameLength) *
                   params->samplesPerBit;
         i += params->samplesPerBit) {
      *(outSignal +
        data_num * params->samplesPerBit *
            (params->frameLength + params->preambleCodeLength) +
        i) = (float32_t)modData[k++];
    }
  }

  if (params->matchedFilterCoeffsLength) {
    float32_t tempData[outLength + params->samplesPerBit * params->FSpan];

    arm_conv_f32(outSignal, outLength, params->matchedFilterCoeffs,
                 params->matchedFilterCoeffsLength, tempData);
    for (i = params->FSpan * params->samplesPerBit / 2;
         i < outLength + params->FSpan * params->samplesPerBit / 2; i++) {
      outSignal[i - params->FSpan * params->samplesPerBit / 2] = tempData[i];
    }
  }
}

void BPSK_demodulateSignal(BPSK_parameters *params, const int8_t *signal,
                           const uint16_t signalLength,
                           const uint16_t *startIdx, const uint16_t startNum,
                           uint8_t *outData, uint16_t *outLength) {
  assert(startNum <= *outLength);
  uint8_t k = 0;

  if (params->cachePrev && params->cachePrev->symbols_left) {
    (*outLength)++;
    k++;
    for (uint16_t j = 0; j < params->cachePrev->symbols_left; j++) {
      params->cachePrev->symbols_buffer[params->frameLength -
                                        params->cachePrev->symbols_left + j] =
          signal[j];
    }

    *(outData) = 0;
    for (uint16_t j = 1; j <= params->frameLength; j++) {
      *(outData) += params->cachePrev->symbols_buffer[j - 1] > 0
                        ? 0
                        : 1 << (params->frameLength - j);
    }
  }

  for (uint16_t i = 0; i < startNum; ++i) {
    *(outData + i + k) = 0;
    for (uint16_t j = 1; j <= params->frameLength; ++j) {
      uint8_t sym = *(signal + *(startIdx + i) + j) > 0 ? 0 : 1;
      *(outData + i + k) += sym << (params->frameLength - j);
    }
  }
}

void BPSK_findSymbolsStarts(BPSK_parameters *params, int8_t *signal,
                            const uint16_t signalLength, uint16_t *startIdx,
                            uint16_t *foundIdx) {

  uint16_t frame_length = params->frameLength + params->preambleCodeLength;
  if (params->cacheNext && params->cacheNext->symbols_left) {
    for (uint8_t i = 0; i < params->frameLength; i++) {
      params->cachePrev->symbols_buffer[i] =
          params->cacheNext->symbols_buffer[i];
    }
  }
  params->cachePrev->symbols_left = params->cacheNext->symbols_left;

  uint16_t corrLength = signalLength > params->preambleCodeLength
                            ? (2 * signalLength - 1)
                            : (2 * params->preambleCodeLength - 1);

  int8_t corr[corrLength];
  xcorrelate_int8(signal, signalLength, params->preambleCode,
                  params->preambleCodeLength, corr, corrLength);

  uint16_t plusMax[signalLength / frame_length * 2];
  uint16_t minusMax[signalLength / frame_length * 2];
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
  if (!maximasCount) {
    return;
  }
  uint16_t prevStart, nextStart;
  prevStart = *(maximas);
  nextStart = prevStart;

  for (uint16_t i = 1; i < maximasCount; i++) {
    nextStart = *(maximas + i);
    uint16_t distance = nextStart - prevStart;
    if (distance == frame_length) {
      *(startIdx + (*foundIdx)) = prevStart + params->preambleCodeLength;
      *foundIdx = *foundIdx + 1U;
      prevStart = nextStart;
    } /*else if (distance >= frame_length - 1U && distance <= frame_length + 1U)
    {
      *(startIdx + (*foundIdx)) = nextStart - frame_length +
    params->preambleCodeLength; prevStart = nextStart; *foundIdx = *foundIdx +
    1U;
    } */
    else if (distance > frame_length + 2U && (distance % frame_length) == 0) {
      prevStart = nextStart;
      uint16_t num = distance / frame_length;
      while (num) {
        nextStart -= frame_length;
        *(startIdx + *foundIdx + num - 1U) =
            nextStart + params->preambleCodeLength;
        --num;
      }
      *(foundIdx) = *(foundIdx) + distance / frame_length;
    } else {
      uint16_t distance_ = *(maximas + i + 1) - nextStart;
      if (distance_ % frame_length == 0) {
        prevStart = nextStart;
      }
    }
    if (*foundIdx == 1) {
      uint16_t first = *startIdx;
      if (first > frame_length) {
        uint16_t num = first / frame_length > 0 ? first / frame_length : 1;
        if (first - num * frame_length < params->nextStart)
          num--;
        *(startIdx + num) = *startIdx;
        while (num--) {
          first -= frame_length;
          *(startIdx + num) = first;
          *foundIdx = *foundIdx + 1;
        }
      } else {
        if (first < params->nextStart)
          *startIdx = *startIdx - 1;
      }
    }
  }
  uint16_t last = *(startIdx + *foundIdx - 1U);

  if (last + (params->frameLength + params->preambleCodeLength) <
      signalLength) {
    uint16_t num = ((signalLength - last) / frame_length > 0)
                       ? (signalLength - last) / frame_length
                       : 1;
    for (uint16_t i = 0; i < num; i++) {
      *(startIdx + *foundIdx) = last + (i + 1) * frame_length;
      *foundIdx = *foundIdx + 1;
    }
  }

  if (params->cacheNext && *foundIdx > 1) {
    uint16_t lastIdx = *(startIdx + *foundIdx - 1);
    if (lastIdx + params->frameLength >= signalLength) {
      *foundIdx = *foundIdx - 1;
      params->nextStart = (lastIdx + params->frameLength) % signalLength;
      params->cacheNext->symbols_left =
          (lastIdx + params->frameLength) % (signalLength - 1);
      uint8_t len = params->frameLength - params->cacheNext->symbols_left;
      for (uint8_t i = 0; i < len; i++) {
        params->cacheNext->symbols_buffer[i] = signal[signalLength - len + i];
      }
    } else {
      params->nextStart = 0;
      params->cacheNext->symbols_left = 0;
    }
  }
}

void BPSK_reset(BPSK_parameters *params) {
  // Initialize costas loop
  params->costas->error_int = 0.0f;
  params->costas->period = 1000000.0f / (float32_t)params->Fc;
  params->costas->f = -(MICROSECONDS / (float32_t)params->Fs) * 2.0f * M_PI /
                      (MICROSECONDS / (float32_t)params->Fc);
  params->costas->phase = 0;
  FIR_filter_reset(params->costas->LP_filterI);
  FIR_filter_reset(params->costas->LP_filterQ);

  params->gardner->v_i = 0.0f;
  params->gardner->last_sample = 0.0f;
  params->gardner->sample = 0.0f;
  params->gardner->sample_zc = 0.0f;
  params->gardner->error = 0;
  params->gardner->mu = 0.0f;
  params->gardner->cnt = 1.0f;
  params->gardner->strobe = 0;

  params->cacheNext->symbols_left = 0;
  params->cachePrev->symbols_left = 0;
  params->nextStart = 0;
}

float BPSK_carrierRecovery(BPSK_parameters *params, float32_t *signal,
                           const uint16_t signalLength) {
  float32_t REF_PERIOD = MICROSECONDS / (float32_t)params->Fc;
  float32_t SAMPLE_PERIOD = MICROSECONDS / (float32_t)params->Fs;

  float32_t a = params->costas->alpha;
  float32_t b = params->costas->beta;
  float32_t intgralf = 0.6f * a / REF_PERIOD;
  float32_t I, Q, S;
  float32_t S_I, S_Q;
  float32_t error;
  float32_t lock;
  float32_t locked = 0;
  float32_t all = 0;

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

    lock = (S_I * S_I) - (S_Q * S_Q);
    ++all;
    if (lock > 0.011f)
      ++locked;
  }
  return locked / all;
}

static float32_t interpolateLinear(const float32_t *signal, uint16_t m_k,
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
                         uint16_t *outputLength) {
  gardnerTimingRecovery_parameters *gardner = params->gardner;
  float32_t Ki = params->gardner->Ki;
  float32_t Kp = params->gardner->Kp;
  uint16_t midpoint_offset = params->samplesPerBit / 2;

  float32_t sps_invert = 1.0f / params->samplesPerBit;

  assert(*outputLength == signalLength / params->samplesPerBit + 1);

  uint16_t k = 0;  // current symbol index
  float32_t v;     // PI output
  float32_t error; // Error from TED
  float32_t W = 0.0f;
  uint16_t zc_idx = 0;
  uint16_t m_k = 0;

  for (uint16_t i = 0; i < signalLength; i++) {
    if (params->gardner->strobe) {
      error = gardner->error;
      output[k++] = gardner->sample > 0 ? 1 : -1;
    } else {
      error = 0.0f;
    }

    gardner->v_i += (Ki * error);
    v = Kp * error + gardner->v_i;

    W = sps_invert + v;
    params->gardner->strobe = params->gardner->cnt < W;
    if (params->gardner->strobe) {
      m_k = i;
      gardner->mu = params->gardner->cnt / W;
      gardner->sample = signal[m_k] =
          interpolateLinear(signal, m_k, gardner->mu);
      zc_idx = m_k - midpoint_offset;
      if (m_k - midpoint_offset > 0)
        gardner->sample_zc = interpolateLinear(signal, zc_idx, gardner->mu);
      gardner->error =
          gardner->sample_zc * (gardner->last_sample - gardner->sample);
      gardner->last_sample = gardner->sample;
    }

    gardner->cnt = gardner->cnt - W;
    if (gardner->cnt > 0) {
      gardner->cnt = gardner->cnt - (float32_t)((int16_t)gardner->cnt);
    } else {
      gardner->cnt = 1.0f + gardner->cnt - (float32_t)((int16_t)gardner->cnt);
    }
  }
  if (m_k + midpoint_offset < signalLength) {
    gardner->sample_zc =
        interpolateLinear(signal, m_k + midpoint_offset, gardner->mu);
  }
  *outputLength = k;
}
