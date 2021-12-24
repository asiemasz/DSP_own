#include "BPSK.h"
#include "uart.h"

static UART_initStruct uart2;

// turn bytes into -1 (corresponding to 0) and 1 (corresponding to 1) array
static void BPSK_generateModData(uint8_t *data, uint16_t length,
                                 int8_t *output) {
  for (uint16_t i = 0; i < length; i++) {
    for (uint8_t j = 0; j < 8; j++) {
      *(output + i * 8 + j) = (data[i] & (0x80 >> j)) ? -1 : 1;
    }
  }
}

static void BPSK_generateDifferentialModData(uint8_t *data, uint16_t length,
                                             int8_t *output) {}

// generate modulation samples (with oversampling)
void BPSK_getModSamples(BPSK_parameters *params, uint8_t *data, uint16_t length,
                        float32_t *outData, uint16_t outLength) {

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

// generate output signal
void BPSK_getOutputSignal(BPSK_parameters *params, uint8_t *data,
                          uint16_t dataLength, float32_t *outSignal,
                          uint16_t outLength) {

  assert(params->samplesPerBit > 0);
  assert(outLength ==
         params->samplesPerBit * dataLength * 8 + params->prefixLength);

  BPSK_getModSamples(params, data, dataLength, outSignal + params->prefixLength,
                     outLength - params->prefixLength);

  for (uint16_t i = 0; i < params->prefixLength; i++) {
    outSignal[i] = outSignal[outLength - params->prefixLength + i];
  }
}

void BPSK_demodulateSignal(BPSK_parameters *params, float32_t *signal,
                           uint16_t signalLength, uint8_t *outData,
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

void BPSK_syncInputSignal(BPSK_parameters *params, float32_t *signal,
                          uint16_t signalLength, uint16_t *startIdx,
                          uint16_t *foundIdx) {

  /*uart2.baudRate = 115200;
  uart2.mode = UART_TRANSMITTER_ONLY;
  uart2.oversampling = UART_OVERSAMPLING_BY_16;
  uart2.parityControl = UART_PARITY_CONTROL_DISABLED;
  uart2.stopBits = UART_STOP_BITS_1;
  uart2.wordLength = UART_WORD_LENGTH_8;
  uart2.uart = USART2;
*/
  *foundIdx = 0;

  uint16_t syncDataLength =
      signalLength - params->prefixLength - params->frameLength;

  float32_t sync[syncDataLength];
  sync[0] = 0;

  for (uint16_t i = 1; i < syncDataLength; i++) {
    sync[i] = sync[i - 1] -
              signal[i - 1] * signal[i - 1 + params->frameLength] +
              signal[i - 1 + params->prefixLength] *
                  signal[i - 1 + params->prefixLength + params->frameLength];
  }

  /* char buf[40];
   sprintf(buf, "\r\n Sync: \r\n");
   uart_sendString(&uart2, buf);
   for (uint16_t i = 0; i < syncDataLength; i++) {
     sprintf(buf, "%f \r\n", sync[i]);
     uart_sendString(&uart2, buf);
   }*/

  uint16_t idx = 0;
  uint16_t start = 0;

  float32_t absMax;
  uint16_t absMaxIdx;
  arm_max_f32(sync, syncDataLength, &absMax, &absMaxIdx);

  float32_t maxVal;
  uint16_t maxIdx;

  while (start < syncDataLength) {
    if (start + params->frameLength < syncDataLength) {
      arm_max_f32(sync + start, params->frameLength, &maxVal, &maxIdx);
    } else {
      arm_max_f32(sync + start, syncDataLength - start - 1, &maxVal, &maxIdx);
    }
    maxIdx += start;
    if ((absMax - maxVal) < 0.2 * absMax) {
      start = maxIdx + params->frameLength;
      *(startIdx + *foundIdx) = maxIdx;
      ++(*foundIdx);
    } else {
      start = start + params->frameLength;
    }
  }
  /*sprintf(buf, "Found %d starts \r\n", *(foundIdx));
  uart_sendString(&uart2, buf);

  for (uint16_t i = 0; i < *foundIdx; i++) {
    sprintf(buf, "\r\n %ld \r\n", *(startIdx + i));
    uart_sendString(&uart2, buf);
  }

  sprintf(buf, "\r\n Koniec \r\n");
  uart_sendString(&uart2, buf);*/
}

void BPSK_syncInputSignal_(BPSK_parameters *params, float32_t *signal,
                           uint16_t signalLength, uint16_t *startIdx,
                           uint16_t *foundIdx) {
  /*uart2.baudRate = 115200;
   uart2.mode = UART_TRANSMITTER_ONLY;
   uart2.oversampling = UART_OVERSAMPLING_BY_16;
   uart2.parityControl = UART_PARITY_CONTROL_DISABLED;
   uart2.stopBits = UART_STOP_BITS_1;
   uart2.wordLength = UART_WORD_LENGTH_8;
   uart2.uart = USART2;
   */
  float32_t B[params->prefixLength];
  arm_fill_f32(1.0f / params->prefixLength, B, params->prefixLength);
  float32_t temp[signalLength - params->frameLength];
  for (uint16_t i = 0; i < signalLength - params->frameLength; ++i) {
    temp[i] = signal[i] * signal[i + params->frameLength];
  }
  uint16_t syncDataLength =
      signalLength - params->frameLength + params->prefixLength - 1;

  float32_t sync[syncDataLength];
  arm_conv_f32(temp, signalLength - params->frameLength, B,
               params->prefixLength, sync);

  /* char buf[40];
   sprintf(buf, "\r\n Sync: \r\n");
   uart_sendString(&uart2, buf);
   for (uint16_t i = 0; i < syncDataLength; i++) {
     sprintf(buf, "%f \r\n", sync[i]);
     uart_sendString(&uart2, buf);
   }*/

  uint16_t idx = 0;
  uint16_t start = 0;

  float32_t absMax;
  uint16_t absMaxIdx;
  arm_max_f32(sync, syncDataLength, &absMax, &absMaxIdx);

  float32_t maxVal;
  uint16_t maxIdx;

  while (start < syncDataLength) {
    if (start + params->frameLength < syncDataLength) {
      arm_max_f32(sync + start, params->frameLength, &maxVal, &maxIdx);
    } else {
      arm_max_f32(sync + start, syncDataLength - start - 1, &maxVal, &maxIdx);
    }
    maxIdx += start;
    if ((absMax - maxVal) < 0.2 * absMax) {
      start = maxIdx + params->frameLength;
      *(startIdx + *foundIdx) = maxIdx;
      ++(*foundIdx);
    } else {
      start = start + params->frameLength;
    }
  }
  /* sprintf(buf, "Found %d starts \r\n", *(foundIdx));
   uart_sendString(&uart2, buf);

   for (uint16_t i = 0; i < *foundIdx; i++) {
     sprintf(buf, "\r\n %ld \r\n", *(startIdx + i));
     uart_sendString(&uart2, buf);
   }

   sprintf(buf, "\r\n Koniec \r\n");
   uart_sendString(&uart2, buf);*/
}

void BPSK_init(BPSK_parameters *params) {
  params->costas->error = 0;
  params->costas->omega = 2.0f * PI * params->Fc / params->Fs;
  params->costas->phase = 0;
}

void BPSK_syncSignalCarrier(BPSK_parameters *params, float32_t *signal,
                            uint16_t signalLength) {
  float32_t errorTot = 0.0f;
  for (uint16_t i = 0; i < signalLength; ++i) {
    float32_t si, sq, sample, sim, sqm;
    sample = signal[i];
    params->costas->phase += params->costas->omega;
    params->costas->phase += params->costas->alpha * params->costas->error;

    params->costas->omega += params->costas->beta * params->costas->error;

    float32_t freq = params->costas->omega * params->Fs / (2 * PI);

    if (params->costas->phase > 2 * PI) {
      params->costas->phase -= 2 * PI;
    }

    si = arm_cos_f32(params->costas->phase);
    sq = -arm_sin_f32(params->costas->phase);
    sim = si * sample;
    sqm = sq * sample;

    sim = IIR_filter_step(params->costas->LP_filterI, sim);
    sqm = IIR_filter_step(params->costas->LP_filterQ, sqm);

    params->costas->error = sim * sqm;

    errorTot += params->costas->error;
  }
}
