#ifndef BPSK_H
#define BPSK_H
#include <stdint.h>

typedef struct BPSK_parameters {
    uint16_t Fc; // carrier frequency
    uint16_t Fs; // sampling frequency
    uint16_t Fb; // bit rate (bps)
} BPSK_parameters;


#endif
