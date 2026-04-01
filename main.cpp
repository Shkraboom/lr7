#include "calcfft_original.h"
#include "calcfft_optimized.h"
#include "param.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <cmath>
#include <functional>

using Clock = std::chrono::high_resolution_clock;

double measure(int iterations, const std::function<void()>& fn) {
    std::vector<double> times;
    times.reserve(iterations);
    for (int i = 0; i < iterations; ++i) {
        auto t0 = Clock::now();
        fn();
        auto t1 = Clock::now();
        times.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    return std::accumulate(times.begin(), times.end(), 0.0) / times.size();
}

int main() {
    const int N = 1024;
    const int iters = 20;

    Param param;
    param.setSampleSize(N);
    param.m_dFunctionSize = 1.0;
    param.m_bIsCircle = true;
    param.RecalculateSteps();

    CalcFFTOriginal orig;
    CalcFFTOptimized opt;

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "N = " << N << ", iterations = " << iters << "\n\n";

    {
        double tOrig = measure(iters, [&]() {
            Sample<double> f(N, N);
            orig.CreateFunction(param, f);
        });
        double tOpt = measure(iters, [&]() {
            Sample<double> f(N, N);
            opt.CreateFunction(param, f);
        });
        std::cout << "CreateFunction (opt 2+3):\n";
        std::cout << "  Before: " << tOrig << " us\n";
        std::cout << "  After:  " << tOpt  << " us\n";
        std::cout << "  Saved:  " << (tOrig - tOpt) << " us ("
                  << ((tOrig - tOpt) / tOrig * 100.0) << "%)\n\n";
    }

    {
        double tOrig = measure(iters, [&]() {
            SampleComplex sc(N, N);
            for (int j = 0; j < N; ++j)
                for (int i = 0; i < N; ++i)
                    sc(i, j) = Complex(i + j * 0.01, 0);
            orig.ShiftSample(sc);
        });
        double tOpt = measure(iters, [&]() {
            SampleComplex sc(N, N);
            for (int j = 0; j < N; ++j)
                for (int i = 0; i < N; ++i)
                    sc(i, j) = Complex(i + j * 0.01, 0);
            opt.ShiftSample(sc);
        });
        std::cout << "ShiftSample (opt 1):\n";
        std::cout << "  Before: " << tOrig << " us\n";
        std::cout << "  After:  " << tOpt  << " us\n";
        std::cout << "  Saved:  " << (tOrig - tOpt) << " us ("
                  << ((tOrig - tOpt) / tOrig * 100.0) << "%)\n\n";
    }

    {
        Sample<double> funcTmp(N, N);
        orig.CreateFunction(param, funcTmp);

        double tOrig = measure(iters, [&]() {
            SampleComplex sc = SampleComplex::FromSampleDouble(funcTmp);
            orig.CalcFourier(sc);
        });
        double tOpt = measure(iters, [&]() {
            SampleComplex sc = SampleComplex::FromSampleDouble(funcTmp);
            opt.CalcFourier(sc);
        });
        std::cout << "CalcFourier (opt 1+2+4):\n";
        std::cout << "  Before: " << tOrig << " us\n";
        std::cout << "  After:  " << tOpt  << " us\n";
        std::cout << "  Saved:  " << (tOrig - tOpt) << " us ("
                  << ((tOrig - tOpt) / tOrig * 100.0) << "%)\n\n";
    }

    {
        Sample<double> funcTmp(N, N);
        orig.CreateFunction(param, funcTmp);
        SampleComplex sc = SampleComplex::FromSampleDouble(funcTmp);
        orig.CalcFourier(sc);

        double tOrig = measure(iters, [&]() {
            Sample<double> res(N, N);
            for (int j = 0; j < N; ++j)
                for (int i = 0; i < N; ++i)
                    res(i, j) = sc.Modulus(i, j);
        });
        double tOpt = measure(iters, [&]() {
            Sample<double> res(N, N);
            double* resPtr = res.GetPointer();
            const Complex* scData = sc.GetPointer();
            const size_t total = static_cast<size_t>(N) * N;
            for (size_t k = 0; k < total; ++k)
                resPtr[k] = std::abs(scData[k]);
        });
        std::cout << "Modulus extraction (opt 2):\n";
        std::cout << "  Before: " << tOrig << " us\n";
        std::cout << "  After:  " << tOpt  << " us\n";
        std::cout << "  Saved:  " << (tOrig - tOpt) << " us ("
                  << ((tOrig - tOpt) / tOrig * 100.0) << "%)\n\n";
    }

    {
        double tOrig = measure(iters, [&]() {
            Sample<double> f, r;
            orig.Calc(param, f, r);
        });
        double tOpt = measure(iters, [&]() {
            Sample<double> f, r;
            opt.Calc(param, f, r);
        });
        std::cout << "=== OVERALL Calc ===\n";
        std::cout << "  Before: " << tOrig << " us  (" << tOrig / 1000.0 << " ms)\n";
        std::cout << "  After:  " << tOpt  << " us  (" << tOpt / 1000.0 << " ms)\n";
        std::cout << "  Saved:  " << (tOrig - tOpt) << " us ("
                  << ((tOrig - tOpt) / tOrig * 100.0) << "%)\n";
        std::cout << "  Speedup: " << std::setprecision(2) << (tOrig / tOpt) << "x\n";
    }

    return 0;
}
