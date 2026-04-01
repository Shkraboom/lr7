#ifndef CALCFFT_OPTIMIZED_H
#define CALCFFT_OPTIMIZED_H

#include "param.h"
#include "Sample.hpp"
#include "SampleComplex.hpp"

class CalcFFTOptimized {
public:
    CalcFFTOptimized() = default;
    ~CalcFFTOptimized() = default;

    void Calc(const Param& oParam_p, Sample<double>& oFunc_p, Sample<double>& oRes_p);

private:
    void CreateFunction(const Param& oParam_p, Sample<double>& oFunction_p);
    void CalcFourier(SampleComplex& oSample_p);
    void ShiftSample(SampleComplex& oSample_p);
};

#endif
