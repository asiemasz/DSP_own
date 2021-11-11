#include "SRRC_filter.h"
#include <arm_math.h>

/* Generate FIR coeffs of Square Root Raised Cosine filter
    symbolSpan - filter span in symbols,
    SPB - samples per bit
    beta - rolloff factor (usually used 0.2)
*/
void SRRC_getFIRCoeffs(uint8_t symbolSpan, uint8_t SPB, float32_t beta, float32_t* firCoeffs, uint16_t length) {
    assert(length == symbolSpan*SPB + 1);
    assert(symbolSpan * SPB % 2 == 0);

    float32_t presc;
    arm_sqrt_f32(symbolSpan, &presc);

    for( uint8_t i = 0; i < length; i++) {
        int16_t n = i - length/2;
        if(n == 0) {
            *(firCoeffs + i) = (1 / presc) * (1 - beta + 4*beta/PI);
        }
        else if(n == symbolSpan/4/beta || n == -symbolSpan/4/beta) {
            *(firCoeffs + i) = beta/presc/1.4142f * ((1 + 2/PI) * sin(PI/4/beta) + (1 - 2/PI) * cos(PI/4/beta));
        }
        else 
            *(firCoeffs + i) = (sin(PI*n*(1 - beta)/symbolSpan) + 4*beta*n/symbolSpan*cos(PI*n*(1+beta)/symbolSpan)) / (PI*n/symbolSpan * (1 - 16*beta*beta*n*n/(symbolSpan*symbolSpan))) / presc;
    }
}
