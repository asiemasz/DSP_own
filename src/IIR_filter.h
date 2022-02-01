/** Created as a port of biz library for DSP in Java (www.source-code.biz,
www.inventec.ch/chdh) for C ARM */
#ifndef IIR_FILTER_H
#define IIR_FILTER_H

#include <arm_math.h>
#include <stdint.h>

typedef struct {
  float32_t *coeffs_A; // Applied to output values
  float32_t *coeffs_B; // Applied to input values
  uint16_t n_A;        // Output signal delay line length
  uint16_t n_B;        // Input signal delay line length
  float32_t *buf_B;    // Input signal delay line
  float32_t *buf_A;    // Output signal delay line
  uint16_t pos_B;      // Position in input signal delay line
  uint16_t pos_A;      // Position in output signal delay line
} IIR_filter;

IIR_filter IIR_filter_init(float32_t *coeffs_A, uint16_t length_A,
                           float32_t *coeffs_B, uint16_t length_B,
                           float32_t *buf_A, float32_t *buf_B);

float32_t IIR_filter_step(IIR_filter *filter, float32_t input_value);

#endif
