#include "param.h"

namespace {
constexpr int kDefaultSize = 512;
constexpr double kDefaultFunctionSizeMm = 1.0;
}

Param::Param() {
    Reset();
}

void Param::Reset() {
    m_iSize = kDefaultSize;
    m_dFunctionSize = kDefaultFunctionSizeMm;
    m_bIsCircle = true;
    RecalculateSteps();
}

void Param::RecalculateSteps() {
    if (m_iSize <= 0 || m_dFunctionSize <= 0) {
        m_dStepFunction = 0;
        m_dStepSpectr = 0;
        return;
    }
    m_dStepFunction = m_dFunctionSize / static_cast<double>(m_iSize);
    m_dStepSpectr = 1.0 / m_dFunctionSize;
}
