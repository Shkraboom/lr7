#ifndef CALCFFT_ORIGINAL_H
#define CALCFFT_ORIGINAL_H

#include "param.h"
#include "Sample.hpp"
#include "SampleComplex.hpp"

class CalcFFTOriginal {
public:
    CalcFFTOriginal() = default;
    ~CalcFFTOriginal() = default;

    void Calc(const Param& oParam_p, Sample<double>& oFunc_p, Sample<double>& oRes_p);
    void CreateFunction(const Param& oParam_p, Sample<double>& oFunction_p);
    void CalcFourier(SampleComplex& oSample_p);
    void ShiftSample(SampleComplex& oSample_p);
};

#endif
