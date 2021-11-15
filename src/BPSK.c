#include "BPSK.h"

//turn bytes into -1 (corresponding to 0) and 1 (corresponding to 1) array
static void BPSK_generateModData(uint8_t* data, uint16_t length, int8_t* output) {
    for(uint16_t i = 0; i < length; i++) {
        for(uint8_t j = 0; j < 8; j++) {
            *(output + i*8 + j) = (data[i] & (0x80 >> j)) ? -1 : 1;
        }
    }
}

//generate modulation samples (with oversampling)
void BPSK_getModSamples(BPSK_parameters* params, uint8_t* data, uint16_t length, float32_t* outData, uint16_t outLength) {
    uint8_t samplesPerBit = params->Fs / params->Fb; //calculate samples per bit parameter (Fs should be multiply of Fb)

    assert(samplesPerBit > 0);
    assert(outLength == samplesPerBit * length * 8);

    int8_t modData[length * 8];
    BPSK_generateModData(data, length, modData);

    for(uint16_t i = 0; i < outLength; i = i + samplesPerBit) {
        for(uint16_t j = 0; j < samplesPerBit; j++) {
            *(outData + i +j) = modData[i/samplesPerBit];
        }
    }

    if(params->firCoeffsLength) {
        float32_t tempData[outLength + params->firCoeffsLength - 1];

        arm_conv_f32(outData, outLength, params->firCoeffs, params->firCoeffsLength, tempData);
        float32_t maxVal;
        void* x;
        arm_max_f32(tempData, outLength, &maxVal, x);
        for(uint16_t i = 0; i < outLength; i++ ) {
           outData[i] = tempData[i]/maxVal;
        }
    }
}

//generate output signal
void BPSK_getOutputSignal(BPSK_parameters* params, uint8_t* data, uint16_t dataLength, float32_t* outSignal, uint16_t outLength) {
    uint8_t samplesPerBit = params->Fs / params->Fb; //calculate samples per bit parameter (Fs should be multiply of Fb)
    
    assert(samplesPerBit > 0);
    assert(outLength == samplesPerBit * dataLength * 8);

    BPSK_getModSamples(params, data, dataLength, outSignal, outLength);

    float32_t fn = (float32_t)params->Fc / params-> Fs;

    for(uint16_t i = 0; i < outLength; i++) {
        outSignal[i] = outSignal[i] * arm_sin_f32(fn*i*2*PI);
    }
}
