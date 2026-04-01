#include "calcfft_optimized.h"
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

/*
 * Оптимизация 3: Вынос инвариантного условия (isCircle) из внутреннего цикла.
 * Оптимизация 2: Прямой доступ к памяти через GetPointer() без проверки границ.
 *                Предвычисление half^2 и y_mm^2 вне внутреннего цикла.
 */
void CalcFFTOptimized::CreateFunction(const Param& oParam_p, Sample<double>& oFunction_p) {
    const int n = oParam_p.sampleSize();
    const double dx = oParam_p.stepFunction();
    if (n <= 0 || dx <= 0) {
        throw std::invalid_argument("CalcFFTOptimized::CreateFunction: invalid N or step.");
    }
    const double half = 0.5 * oParam_p.m_dFunctionSize;
    const double cx = 0.5 * static_cast<double>(n - 1);
    const double cy = 0.5 * static_cast<double>(n - 1);

    double* data = oFunction_p.GetPointer();
    const size_t stride = static_cast<size_t>(n);

    if (oParam_p.m_bIsCircle) {
        const double r2 = half * half;
        for (int j = 0; j < n; ++j) {
            const double y_mm = (static_cast<double>(j) - cy) * dx;
            const double y2 = y_mm * y_mm;
            double* row = data + static_cast<size_t>(j) * stride;
            for (int i = 0; i < n; ++i) {
                const double x_mm = (static_cast<double>(i) - cx) * dx;
                row[i] = (x_mm * x_mm + y2 <= r2) ? 1.0 : 0.0;
            }
        }
    } else {
        for (int j = 0; j < n; ++j) {
            const double y_mm = (static_cast<double>(j) - cy) * dx;
            double* row = data + static_cast<size_t>(j) * stride;
            if (std::abs(y_mm) > half) {
                for (int i = 0; i < n; ++i) {
                    row[i] = 0.0;
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    const double x_mm = (static_cast<double>(i) - cx) * dx;
                    row[i] = (std::abs(x_mm) <= half) ? 1.0 : 0.0;
                }
            }
        }
    }
}

/*
 * Оптимизация 1: Сдвиг квадрантов на месте (in-place) без выделения
 * временного SampleComplex. Вместо полного копирования через tmp
 * используем std::swap для обмена элементов четырёх квадрантов.
 */
void CalcFFTOptimized::ShiftSample(SampleComplex& oSample_p) {
    const int nx = oSample_p.GetSizeX();
    const int ny = oSample_p.GetSizeY();
    const int sx = nx / 2;
    const int sy = ny / 2;

    Complex* data = oSample_p.GetPointer();
    const size_t stride = static_cast<size_t>(nx);

    for (int j = 0; j < sy; ++j) {
        Complex* row_top = data + static_cast<size_t>(j) * stride;
        Complex* row_bot = data + static_cast<size_t>(j + sy) * stride;
        for (int i = 0; i < sx; ++i) {
            std::swap(row_top[i], row_bot[i + sx]);
            std::swap(row_top[i + sx], row_bot[i]);
        }
    }
}

/*
 * Оптимизация 4: In-place FFT — используем reinterpret_cast для прямого
 * преобразования std::complex<double>* в fftw_complex* (гарантировано
 * совместимы по раскладке памяти). Исключаем выделение двух массивов
 * fftw_complex и два цикла копирования N^2 элементов.
 * Оптимизация 2: Прямой доступ к данным без bounds-checking.
 */
void CalcFFTOptimized::CalcFourier(SampleComplex& oSample_p) {
    const int nx = oSample_p.GetSizeX();
    const int ny = oSample_p.GetSizeY();
    ShiftSample(oSample_p);

    Complex* data = oSample_p.GetPointer();
    fftw_complex* fdata = reinterpret_cast<fftw_complex*>(data);

    fftw_plan plan = fftw_plan_dft_2d(ny, nx, fdata, fdata, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    const double scale = 1.0 / (static_cast<double>(nx) * static_cast<double>(ny));
    const size_t total = static_cast<size_t>(nx) * static_cast<size_t>(ny);
    for (size_t k = 0; k < total; ++k) {
        data[k] *= scale;
    }

    ShiftSample(oSample_p);
}

void CalcFFTOptimized::Calc(const Param& oParam_p, Sample<double>& oFunc_p, Sample<double>& oRes_p) {
    const int n = oParam_p.sampleSize();
    if (n <= 0) {
        throw std::invalid_argument("CalcFFTOptimized::Calc: sample size must be positive.");
    }
    oFunc_p = Sample<double>(n, n);
    oRes_p = Sample<double>(n, n);
    CreateFunction(oParam_p, oFunc_p);

    SampleComplex sc(n, n);
    const double* srcData = oFunc_p.GetPointer();
    Complex* dstData = sc.GetPointer();
    const size_t total = static_cast<size_t>(n) * static_cast<size_t>(n);
    for (size_t k = 0; k < total; ++k) {
        dstData[k] = Complex(srcData[k], 0.0);
    }

    CalcFourier(sc);

    double* resPtr = oRes_p.GetPointer();
    const Complex* scData = sc.GetPointer();
    for (size_t k = 0; k < total; ++k) {
        resPtr[k] = std::abs(scData[k]);
    }
}
