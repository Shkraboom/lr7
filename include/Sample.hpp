#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <vector>
#include <stdexcept>
#include <string>

template <typename PAR>
class Sample {
protected:
    std::vector<PAR> m_data;
    int m_size_x;
    int m_size_y;

public:
    Sample() : m_size_x(0), m_size_y(0) {}

    Sample(int size_x, int size_y) : m_size_x(size_x), m_size_y(size_y) {
        if (size_x <= 0 || size_y <= 0) {
            throw std::invalid_argument("Sample: dimensions must be positive.");
        }
        m_data.resize(static_cast<size_t>(size_x) * size_y);
    }

    virtual ~Sample() = default;

    int GetSizeX() const { return m_size_x; }
    int GetSizeY() const { return m_size_y; }

    PAR* GetPointer() {
        return m_data.empty() ? nullptr : m_data.data();
    }

    const PAR* GetPointer() const {
        return m_data.empty() ? nullptr : m_data.data();
    }

    PAR& operator()(int i, int j) {
        if ((i >= m_size_x) || (j >= m_size_y) || (i < 0) || (j < 0)) {
            throw std::out_of_range("Sample: index is out of bounds.");
        }
        return m_data[static_cast<size_t>(m_size_x) * j + i];
    }

    const PAR& operator()(int i, int j) const {
        if ((i >= m_size_x) || (j >= m_size_y) || (i < 0) || (j < 0)) {
            throw std::out_of_range("Sample: index is out of bounds.");
        }
        return m_data[static_cast<size_t>(m_size_x) * j + i];
    }

    virtual Sample& operator/=(const PAR& value) {
        if (value == PAR(0)) {
            throw std::invalid_argument("Division by zero.");
        }
        for (auto& elem : m_data) {
            elem /= value;
        }
        return *this;
    }
};

#endif // SAMPLE_HPP
