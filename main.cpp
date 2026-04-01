#include "calcfft_original.h"
#include "calcfft_optimized.h"
#include "param.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <cmath>

using Clock = std::chrono::high_resolution_clock;

template <typename CalcType>
double benchmarkCalc(CalcType& calc, const Param& param, int iterations) {
    std::vector<double> times;
    times.reserve(iterations);

    for (int iter = 0; iter < iterations; ++iter) {
        Sample<double> func, res;
        auto t0 = Clock::now();
        calc.Calc(param, func, res);
        auto t1 = Clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        times.push_back(ms);
    }

    return std::accumulate(times.begin(), times.end(), 0.0) / times.size();
}

void runBenchmark(int sampleSize, int iterations) {
    Param param;
    param.setSampleSize(sampleSize);
    param.m_dFunctionSize = 1.0;
    param.m_bIsCircle = true;
    param.RecalculateSteps();

    std::cout << "=== N = " << sampleSize
              << ", iterations = " << iterations << " ===" << std::endl;

    CalcFFTOriginal original;
    double tOrig = benchmarkCalc(original, param, iterations);

    CalcFFTOptimized optimized;
    double tOpt = benchmarkCalc(optimized, param, iterations);

    double speedup = tOrig / tOpt;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Original:  " << tOrig << " ms (avg)" << std::endl;
    std::cout << "  Optimized: " << tOpt  << " ms (avg)" << std::endl;
    std::cout << "  Speedup:   " << speedup << "x" << std::endl;

    Sample<double> funcOrig, resOrig, funcOpt, resOpt;
    original.Calc(param, funcOrig, resOrig);
    optimized.Calc(param, funcOpt, resOpt);

    double maxDiff = 0.0;
    for (int j = 0; j < sampleSize; ++j) {
        for (int i = 0; i < sampleSize; ++i) {
            maxDiff = std::max(maxDiff, std::abs(resOrig(i, j) - resOpt(i, j)));
        }
    }
    std::cout << "  Max diff:  " << std::scientific << std::setprecision(6)
              << maxDiff << " (correctness check)" << std::endl;
    std::cout << std::endl;
}

int main() {
    std::cout << "Лабораторная работа №7: Профилирование и оптимизация" << std::endl;
    std::cout << "Сравнение оригинальной и оптимизированной версий CalcFFT" << std::endl;
    std::cout << std::endl;
    std::cout << "Оптимизации:" << std::endl;
    std::cout << "  1. In-place ShiftSample (без временного массива)" << std::endl;
    std::cout << "  2. Прямой доступ к памяти (без проверки границ в operator())" << std::endl;
    std::cout << "  3. Вынос условия isCircle из внутреннего цикла CreateFunction" << std::endl;
    std::cout << "  4. In-place FFT (без копирования в/из fftw_complex массивов)" << std::endl;
    std::cout << std::endl;

    runBenchmark(256, 20);
    runBenchmark(512, 10);
    runBenchmark(1024, 5);

    return 0;
}
