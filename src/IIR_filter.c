#include "IIR_filter.h"

IIR_filter IIR_filter_init(float32_t *coeffs_A, uint16_t length_A,
                           float32_t *coeffs_B, uint16_t length_B,
                           float32_t *buf_A, float32_t *buf_B) {
  IIR_filter ret;
  ret.coeffs_A = coeffs_A;
  ret.coeffs_B = coeffs_B;
  ret.n_A = length_A - 1;
  ret.n_B = length_B - 1;

  ret.buf_B = buf_B;
  ret.buf_A = buf_A;

  ret.pos_B = 0;
  ret.pos_A = 0;

  return ret;
}

float32_t IIR_filter_step(IIR_filter *filter, float32_t input_value) {
  float32_t acc = filter->coeffs_B[0] * input_value;

  for (uint16_t i = 1; i <= filter->n_B; ++i) {
    uint16_t p = (filter->pos_B + filter->n_B - i) % filter->n_B;
    acc += filter->coeffs_B[i] * filter->buf_B[p];
  }
  for (uint16_t i = 1; i <= filter->n_A; ++i) {
    uint16_t p = (filter->pos_A + filter->n_A - i) % filter->n_A;
    acc -= filter->coeffs_A[i] * filter->buf_A[p];
  }

  if (filter->n_B > 0) {
    filter->buf_B[filter->pos_B] = input_value;
    filter->pos_B = (filter->pos_B + 1) % filter->n_B;
  }

  if (filter->n_A > 0) {
    filter->buf_A[filter->pos_A] = acc;
    filter->pos_A = (filter->pos_A + 1) % filter->n_A;
  }

  return acc;
}

float32_t IIR_filter_run(IIR_filter *filter, float32_t *input,
                         uint16_t input_length, float32_t *output) {
  for (uint16_t i = 0; i < input_length; i++) {
    output[i] = IIR_filter_step(filter, input[i]);
  }
  return output[input_length - 1];
}