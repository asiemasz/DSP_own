#ifndef FIR_FILTER_H
#define FIR_FILTER_H

#include <arm_math.h>
#include <stdint.h>

typedef struct {
  float32_t *coeffs;
  uint16_t length;
  float32_t *buf;
  uint16_t idx;
} FIR_filter;

FIR_filter FIR_filter_init(float32_t *coeffs, uint16_t coeffs_length,
                           float32_t *buf);

FIR_filter FIR_filter_reset(FIR_filter *filter);

float32_t FIR_filter_step(FIR_filter *filter, float32_t input_value);

float32_t FIR_filter_run(FIR_filter *filter, float32_t *input,
                         uint16_t input_length, float32_t *output);

#endif