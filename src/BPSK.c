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
        float32_t tempData[outLength + samplesPerBit * params->FSpan];

        arm_conv_f32(outData, outLength, params->firCoeffs, params->firCoeffsLength, tempData);
        float32_t maxVal;
        void* x;
        arm_max_f32(tempData, outLength, &maxVal, x);
        uint32_t k = 0;
        for(uint16_t i = samplesPerBit/2 - 1; i < outLength + samplesPerBit/2; i++ ) {
           outData[k++] = tempData[i] / maxVal;
        }
    }
}

//generate output signal
void BPSK_getOutputSignal(BPSK_parameters* params, uint8_t* data, uint16_t dataLength, float32_t* outSignal, uint16_t outLength) {
    uint8_t samplesPerBit = params->Fs / params->Fb; //calculate samples per bit parameter (Fs should be multiply of Fb)
    
    assert(samplesPerBit > 0);
    assert(outLength == samplesPerBit * dataLength * 8 + params->prefixLength);

    BPSK_getModSamples(params, data, dataLength, outSignal + params->prefixLength, outLength - params->prefixLength);

    float32_t fn = (float32_t)params->Fc / params-> Fs;

    for(uint16_t i = params->prefixLength; i < outLength; i++) {
        outSignal[i] = outSignal[i] * arm_sin_f32(fn*i*2*PI);
    }
    for(uint16_t i = 0; i < params->prefixLength; i++) {
        outSignal[i] = outSignal[outLength - params->prefixLength + i];
    }
}

void BPSK_demodulateSignal(BPSK_parameters* params, float32_t* signal, uint16_t signalLength, uint8_t* outData, uint16_t outLength) {
    uint8_t samplesPerBit = params->Fs / params->Fb; //calculate samples per bit parameter (Fs should be multiply of Fb)
    float32_t fn = (float32_t)params->Fc / params-> Fs;

    for(uint16_t i = 0; i < signalLength; i++) {
        signal[i] = signal[i] * arm_sin_f32(fn*i*2*PI);
    }

    float32_t outSignal[signalLength + samplesPerBit*params->FSpan];

    //arm_conv_partial_f32(signal, signalLength, params->firCoeffs, params->firCoeffsLength, outSignal, samplesPerBit * params->FSpan , signalLength);
    arm_conv_f32(signal, signalLength, params->firCoeffs, params->firCoeffsLength, outSignal);

    uint16_t k = 0;

    for(uint16_t i = samplesPerBit * params->FSpan; i < signalLength + samplesPerBit*params->FSpan; i = i + samplesPerBit*8) { 
        outData[k] = 0;
        for(uint16_t j = 0; j < 8*samplesPerBit; j = j + samplesPerBit) {
            if(outSignal[i+j] > 0)
                outData[k] += (1 << (7 - j/samplesPerBit));
        }
        ++k;
    }
}

void BPSK_syncInputSignal(BPSK_parameters* params, float32_t* signal, uint16_t signalLength, uint32_t* startIdx) {
    uint16_t syncDataLength = signalLength - params->prefixLength - params->frameLength;

    float32_t sync[syncDataLength];
    sync[0] = 0;

    for(uint16_t i = 1; i < syncDataLength; i++) {
        sync[i] = sync[i-1] - signal[i-1]*signal[i-1 + params->frameLength] + signal[i - 1 + params->prefixLength] * signal[i - 1 + params->prefixLength + params->frameLength];
    }

    float32_t dumm;

    arm_max_f32(sync, syncDataLength, &dumm, startIdx);
    *startIdx += params->prefixLength + 1;
}

