#pragma once

#include "vectors.h"

template <class T>
struct BBox3
{
    BBox3() { }
    BBox3(const rtvector<T, 3>& vMin, const rtvector<T, 3>& vMax) : m_vMin(vMin), m_vMax(vMax) { }

    void include(const rtvector<T, 3>& v)
    {
        m_vMin = vmin(v, m_vMin);
        m_vMax = vmax(v, m_vMax);
    }
    void include(const BBox3<T>& b)
    {
        m_vMin = vmin(b.m_vMin, m_vMin);
        m_vMax = vmax(b.m_vMax, m_vMax);
    }
    inline bool operator ==(const BBox3<T> &other) const { return all(m_vMin == other.m_vMin) && all(m_vMax == other.m_vMax); }
    bool includes(const rtvector<T, 3>& v) const { return all(v >= m_vMin) && all(v < m_vMax); }
    const rtvector<T, 3>& operator[](NvU32 u) const { nvAssert(u < 2); return (&m_vMin)[u]; }
    rtvector<T, 3>& operator[](NvU32 u) { nvAssert(u < 2); return (&m_vMin)[u]; }
    rtvector<T, 3> computeCenter() const { return (m_vMin + m_vMax) / (T)2; }
    T evalVolume() const
    {
        T f = m_vMax[0] - m_vMin[0];
        for (NvU32 u = 1; u < 2; ++u) f *= m_vMax[u] - m_vMin[0];
        return f;
    }

    rtvector<T, 3> m_vMin, m_vMax;
};

typedef BBox3<float> BBox3f;
typedef BBox3<double> BBox3d;
