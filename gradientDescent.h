#pragma once

#include <array>
#include "mybasics.h"

template <class Functor>
double2 computeFirstAndSecondDerivatives(Functor f, double fX)
{
    double fDx = fabs(fX) / 1024;
    std::array<double, 3> fValues = { f(fX - fDx), f(fX), f(fX + fDx) };
    std::array<double, 2> fFirstDerivs = { (fValues[1] - fValues[0]) / fDx, (fValues[2] - fValues[1]) / fDx };
    double2 r;
    r[0] = (fValues[2] - fValues[0]) / (2 * fDx);
    r[1] = (fFirstDerivs[1] - fFirstDerivs[0]) / fDx;
    return r;
}

// finds point with zero derivative
template <class Functor>
double searchStationaryPoint(Functor f, double fX0, NvU32 nMaxIterations)
{
    double fX = fX0;
    // use newton method to find the point where derivative is zero
    for (NvU32 u = 0; u < nMaxIterations; ++u)
    {
        double2 derivs = computeFirstAndSecondDerivatives(f, fX);
        if (derivs[0] == 0)
            return fX;
        double fSign = derivs[0] * derivs[1] < 0 ? -1 : 1;
        derivs = double2({ fabs(derivs[0]), fabs(derivs[1]) });
        derivs[1] = std::max(derivs[0] / (abs(fX) * 0.1), derivs[1]);
        fX = fX - fSign * derivs[0] / derivs[1];
    }
    return fX;
}
