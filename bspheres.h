#pragma once

#include "vectors.h"

template <class T>
struct BSphere
{
    bool intersects(const rtvector<T, 3>& p)
    {
        return lengthSquared(p - m_vCenter) <= sqr(m_fRad);
    }
    rtvector<T, 3> m_vCenter;
    T m_fRad;
};

typedef BSphere<float> BSphereF;
typedef BSphere<double> BSphereD;