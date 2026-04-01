#ifndef SAMPLECOMPLEX_HPP
#define SAMPLECOMPLEX_HPP

#include "Sample.hpp"
#include <complex>
#include <fstream>
#include <string>
#include <cmath>

using Complex = std::complex<double>;

class SampleComplex : public Sample<Complex> {
public:
    SampleComplex() : Sample<Complex>() {}

    SampleComplex(int size_x, int size_y) : Sample<Complex>(size_x, size_y) {}

    explicit SampleComplex(int n) : Sample<Complex>(n, n) {}

    double Modulus(int i, int j) const {
        return std::abs((*this)(i, j));
    }

    double Phase(int i, int j) const {
        return std::arg((*this)(i, j));
    }

    double Intensity(int i, int j) const {
        const Complex& z = (*this)(i, j);
        return z.real() * z.real() + z.imag() * z.imag();
    }

    Sample<double> ToSampleDouble() const {
        Sample<double> result(m_size_x, m_size_y);
        for (int j = 0; j < m_size_y; ++j) {
            for (int i = 0; i < m_size_x; ++i) {
                result(i, j) = (*this)(i, j).real();
            }
        }
        return result;
    }

    static SampleComplex FromSampleDouble(const Sample<double>& s) {
        SampleComplex result(s.GetSizeX(), s.GetSizeY());
        for (int j = 0; j < s.GetSizeY(); ++j) {
            for (int i = 0; i < s.GetSizeX(); ++i) {
                result(i, j) = Complex(s(i, j), 0.0);
            }
        }
        return result;
    }

    Sample<double> GetIntensitySample() const {
        Sample<double> result(m_size_x, m_size_y);
        for (int j = 0; j < m_size_y; ++j) {
            for (int i = 0; i < m_size_x; ++i) {
                result(i, j) = Intensity(i, j);
            }
        }
        return result;
    }

    void WriteRealPartToFile(const std::string& path) const {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("Cannot open file: " + path);
        for (int j = 0; j < m_size_y; ++j) {
            for (int i = 0; i < m_size_x; ++i) {
                if (i > 0) f << " ";
                f << (*this)(i, j).real();
            }
            f << "\n";
        }
    }

    void WriteImagPartToFile(const std::string& path) const {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("Cannot open file: " + path);
        for (int j = 0; j < m_size_y; ++j) {
            for (int i = 0; i < m_size_x; ++i) {
                if (i > 0) f << " ";
                f << (*this)(i, j).imag();
            }
            f << "\n";
        }
    }

    void WriteModulusToFile(const std::string& path) const {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("Cannot open file: " + path);
        for (int j = 0; j < m_size_y; ++j) {
            for (int i = 0; i < m_size_x; ++i) {
                if (i > 0) f << " ";
                f << Modulus(i, j);
            }
            f << "\n";
        }
    }

    void WritePhaseToFile(const std::string& path) const {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("Cannot open file: " + path);
        for (int j = 0; j < m_size_y; ++j) {
            for (int i = 0; i < m_size_x; ++i) {
                if (i > 0) f << " ";
                f << Phase(i, j);
            }
            f << "\n";
        }
    }

    void WriteIntensityToFile(const std::string& path) const {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("Cannot open file: " + path);
        for (int j = 0; j < m_size_y; ++j) {
            for (int i = 0; i < m_size_x; ++i) {
                if (i > 0) f << " ";
                f << Intensity(i, j);
            }
            f << "\n";
        }
    }
};

#endif // SAMPLECOMPLEX_HPP
