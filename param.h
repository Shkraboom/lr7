#ifndef PARAM_H
#define PARAM_H

/**
 * Параметры вычислений
 * Шаги по функции и спектру связаны с размером выборки и апертурой.
 */
class Param {
private:
    /// шаг по функции [мм]
    double m_dStepFunction = 0;
    /// шаг по спектру [1/мм] (пространственная частота)
    double m_dStepSpectr = 0;
    /// размер выборки (число отсчётов по стороне)
    int m_iSize = 512;

public:
    /// размер стороны квадрата или диаметр круга [мм]
    double m_dFunctionSize = 0;
    /// true — круг, false — квадрат
    bool m_bIsCircle = true;

    Param();
    ~Param() = default;

    /// Сброс к умолчаниям и пересчёт шагов из m_dFunctionSize и m_iSize
    void Reset();

    /// Пересчитать шаги из текущих m_dFunctionSize и m_iSize (связь ДПФ)
    void RecalculateSteps();

    double stepFunction() const { return m_dStepFunction; }
    double stepSpectr() const { return m_dStepSpectr; }
    int sampleSize() const { return m_iSize; }

    void setStepFunction(double v) { m_dStepFunction = v; }
    void setStepSpectr(double v) { m_dStepSpectr = v; }
    void setSampleSize(int n) { m_iSize = n; }
};

#endif
