#include "FIR_filter.h"

FIR_filter FIR_filter_init(float32_t *coeffs, uint16_t coeffs_length,
                           float32_t *buf) {
  FIR_filter ret;
  ret.coeffs = coeffs;
  ret.length = coeffs_length;
  ret.buf = buf;
  ret.idx = 0;
  return ret;
}

FIR_filter FIR_filter_reset(FIR_filter *filter) {
  for (uint16_t i = 0; i < filter->length; ++i) {
    filter->buf[i] = 0.0f;
  }
  filter->idx = 0;
}

float32_t FIR_filter_step(FIR_filter *filter, float32_t input_value) {
  float32_t ret;
  int16_t i, index;

  filter->buf[filter->idx++] = input_value;
  filter->idx = filter->idx % filter->length;
  ret = 0.0f;

  index = filter->idx;
  for (i = 0; i < filter->length; i++) {
    --index;
    if (index < 0)
      index = filter->length - 1;
    ret += filter->buf[index] * filter->coeffs[i];
  }
  return ret;
}

float32_t FIR_filter_run(FIR_filter *filter, float32_t *input,
                         uint16_t input_length, float32_t *output) {
  for (uint16_t i = 0; i < input_length; i++) {
    output[i] = FIR_filter_step(filter, input[i]);
  }

  return output[input_length];
}