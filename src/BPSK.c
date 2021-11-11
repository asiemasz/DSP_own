#include "BPSK.h"

//turn bytes into -1 (corresponding to 0) and 1 (corresponding to 1) array
static void generateModData(uint8_t* data, uint16_t length, int8_t* output) {
    for(uint16_t i = 0; i < length; i++) {
        for(uint8_t j = 0; j < 8; j++) {
            *(output + i*8 + j) = (data[i] & (0x80 >> j)) ? -1 : 1;
        }
    }
}

//generate modulation samples (with oversampling)
void getModSamples(BPSK_parameters* params, uint8_t* data, uint16_t length, int8_t* outData, uint16_t outLength) {
    uint8_t samplesPerBit = params->Fs / params->Fb; //calculate samples per bit parameter (Fs should be multiply of Fb)

    assert(samplesPerBit > 0);
    assert(outLength == samplesPerBit * length * 8);

    int8_t modData[length * 8];
    generateModData(data, length, modData);

    for(uint16_t i = 0; i < outLength; i = i + samplesPerBit) {
        for(uint16_t j = 0; j < samplesPerBit; j++) {
            *(outData + i +j) = modData[i/samplesPerBit];
        }
    }
}
