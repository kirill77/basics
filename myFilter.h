#pragma once

template <NvU32 LOGN>
struct MyFilter
{
    MyFilter() { memset(this, 0, sizeof(*this)); }

    void addValue(double f)
    {
        nvAssert(!isnan(f));
        NvU32 u = (m_nValues++) & MASK;
        m_fSum -= m_fValues[u];
        m_fValues[u] = f;
        m_fSum += f;
        if (u % 1024 == 0)
        {
            resetSum();
        }
    }
    double getAverage() const
    {
        return m_nValues == 0 ? 0 : m_fSum / std::min(m_nValues, N);
    }

private:
    void resetSum()
    {
        m_fSum = m_fValues[0];
        for (NvU32 u = 1; u < N; ++u) m_fSum += m_fValues[u];
    }
    static const NvU32 N = (1 << LOGN);
    static const NvU32 MASK = N - 1;
    double m_fSum;
    double m_fValues[N];
    NvU32 m_nValues;
};
