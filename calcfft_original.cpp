#include "calcfft_original.h"
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

void CalcFFTOriginal::CreateFunction(const Param& oParam_p, Sample<double>& oFunction_p) {
    const int n = oParam_p.sampleSize();
    const double dx = oParam_p.stepFunction();
    if (n <= 0 || dx <= 0) {
        throw std::invalid_argument("CalcFFTOriginal::CreateFunction: invalid N or step.");
    }
    const double half = 0.5 * oParam_p.m_dFunctionSize;
    const double cx = 0.5 * static_cast<double>(n - 1);
    const double cy = 0.5 * static_cast<double>(n - 1);

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            const double x_mm = (static_cast<double>(i) - cx) * dx;
            const double y_mm = (static_cast<double>(j) - cy) * dx;
            double v = 0.0;
            if (oParam_p.m_bIsCircle) {
                if (x_mm * x_mm + y_mm * y_mm <= half * half) {
                    v = 1.0;
                }
            } else {
                if (std::abs(x_mm) <= half && std::abs(y_mm) <= half) {
                    v = 1.0;
                }
            }
            oFunction_p(i, j) = v;
        }
    }
}

void CalcFFTOriginal::ShiftSample(SampleComplex& oSample_p) {
    const int nx = oSample_p.GetSizeX();
    const int ny = oSample_p.GetSizeY();
    const int sx = nx / 2;
    const int sy = ny / 2;
    SampleComplex tmp(nx, ny);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int ti = (i + sx) % nx;
            const int tj = (j + sy) % ny;
            tmp(ti, tj) = oSample_p(i, j);
        }
    }
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            oSample_p(i, j) = tmp(i, j);
        }
    }
}

void CalcFFTOriginal::CalcFourier(SampleComplex& oSample_p) {
    const int nx = oSample_p.GetSizeX();
    const int ny = oSample_p.GetSizeY();
    ShiftSample(oSample_p);

    fftw_complex* in = fftw_alloc_complex(static_cast<size_t>(nx) * static_cast<size_t>(ny));
    fftw_complex* out = fftw_alloc_complex(static_cast<size_t>(nx) * static_cast<size_t>(ny));
    if (!in || !out) {
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        throw std::runtime_error("CalcFFTOriginal::CalcFourier: fftw_alloc_complex failed.");
    }

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const size_t idx = static_cast<size_t>(j) * static_cast<size_t>(nx) + static_cast<size_t>(i);
            const Complex z = oSample_p(i, j);
            in[idx][0] = z.real();
            in[idx][1] = z.imag();
        }
    }

    fftw_plan plan = fftw_plan_dft_2d(ny, nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    const double scale = 1.0 / (static_cast<double>(nx) * static_cast<double>(ny));
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const size_t idx = static_cast<size_t>(j) * static_cast<size_t>(nx) + static_cast<size_t>(i);
            oSample_p(i, j) = Complex(out[idx][0] * scale, out[idx][1] * scale);
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    ShiftSample(oSample_p);
}

void CalcFFTOriginal::Calc(const Param& oParam_p, Sample<double>& oFunc_p, Sample<double>& oRes_p) {
    const int n = oParam_p.sampleSize();
    if (n <= 0) {
        throw std::invalid_argument("CalcFFTOriginal::Calc: sample size must be positive.");
    }
    oFunc_p = Sample<double>(n, n);
    oRes_p = Sample<double>(n, n);
    CreateFunction(oParam_p, oFunc_p);
    SampleComplex sc = SampleComplex::FromSampleDouble(oFunc_p);
    CalcFourier(sc);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            oRes_p(i, j) = sc.Modulus(i, j);
        }
    }
}
