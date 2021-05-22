#pragma once

#include "bboxes.h"
#include "vectors.h"

struct MyUnitsTest
{
    static void test();
    static inline bool wasTested() { return s_bTested; }
private:
    static bool s_bTested;
};

#define UNIT_MEMBERS(X) \
X() { m_value = 0; } \
explicit X(T value) : m_value(value) { } \
X operator -() const { return X(-MyUnits<T>::m_value); } \
X operator +=(X s) { MyUnits<T>::m_value += s.m_value; return *this; } \
X operator *=(X s) { MyUnits<T>::m_value *= s.m_value; return *this; } \
X operator -=(X s) { MyUnits<T>::m_value -= s.m_value; return *this; } \
X operator /=(T s) { MyUnits<T>::m_value /= s; return *this; } \
X sqrt() const { return X((T)::sqrt(this->m_value)); } \
X operator -(X s) const { return X(this->m_value - s.m_value); } \
X operator +(X s) const { return X(this->m_value + s.m_value); } \
X operator /(X s) const { return X(this->m_value / s.m_value); } \
X operator *(X s) const { return X(this->m_value * s.m_value); } \
X operator *(T s) const { return X(this->m_value * s); } \
bool operator == (X s) const { return this->m_value == s.m_value; } \
bool operator == (T s) const { return this->m_value == s; } \
bool operator != (X s) const { return this->m_value != s.m_value; } \
bool operator != (T s) const { return this->m_value != s; } \
bool operator >= (X s) const { return this->m_value >= s.m_value; } \
bool operator < (X s) const { return this->m_value < s.m_value; } \
bool operator > (X s) const { return this->m_value > s.m_value; }

#define COLOUMB_DEC(X) 8.9875517923##X
#define COLOUMB_CONSTANT COLOUMB_DEC(e9)
#define ELECTRON_DEC(X) 1.602176634##X
#define ELECTRON_CHARGE ELECTRON_DEC(e-19)
#define DALTON_DEC(X) 1.660539066605##X
#define DALTON DALTON_DEC(e-27)

// those constants define base units used at atomic scales
#define TO_KILOGRAMS 1e-20
#define TO_METERS 1e-5
#define TO_COLOUMBS (ELECTRON_CHARGE)
#define TO_SECONDS 1e-6
template <class _T>
struct MyUnits
{
    typedef _T T;
    T m_value;

    UNIT_MEMBERS(MyUnits)

    //*** constructors
    static MyUnits<T> milliSecond() { return MyUnits<T>((T)(1e-3 / TO_SECONDS)); }
    static MyUnits<T> nanoSecond() { return MyUnits<T>((T)(1e-6 / TO_SECONDS)); }
    static MyUnits<T> microMeter() { return MyUnits<T>((T)(1e-6 / TO_METERS)); }
    static MyUnits<T> electron() { return MyUnits<T>((T)(ELECTRON_CHARGE / TO_COLOUMBS)); }
    static MyUnits<T> dalton() { return MyUnits<T>((T)(DALTON / TO_KILOGRAMS)); }

    //*** base units
    T toKilograms() { return (T)(this->m_value * TO_KILOGRAMS); }
    T toMeters() { return (T)(this->m_value * TO_METERS); }
    T toSeconds() { return (T)this->m_value * TO_SECONDS; }
    // we have double instead of T here because charge is declared as MyUnits<int> and we don't want to convert charge to int
    double toColoumbs() { return this->m_value * TO_COLOUMBS; }

    //*** derived units
    T toNewtons() { return (T)(this->m_value * (TO_KILOGRAMS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toSIAcceleration() { return (T)(this->m_value * (TO_METERS / TO_SECONDS / TO_SECONDS)); }
};

template <class T> MyUnits<T> sqrt(MyUnits<T> s) { return MyUnits<T>(::sqrt(s.m_value)); }

template <class UNITS>
inline const typename UNITS::T& removeUnits(const UNITS& s) { return (const typename UNITS::T&)s; }
template <class UNITS, NvU32 N>
inline const rtvector<typename UNITS::T, N> &removeUnits(const rtvector<UNITS, N>& v) { return (const rtvector<typename UNITS::T, N> &)v; }
template <class UNITS>
inline const BBox3<typename UNITS::T>& removeUnits(const BBox3<UNITS>& box) { return (const BBox3<typename UNITS::T> &)box; }

template <class UNITS>
inline const rtvector<UNITS, 3>& setUnits(const rtvector<typename UNITS::T, 3>& v) { return (const rtvector<UNITS, 3> &)v; }
template <class UNITS>
inline const BBox3<UNITS>& setUnits(const BBox3<typename UNITS::T>& box) { return (const BBox3<UNITS> &)box; }

template <class T>
rtvector<MyUnits<T>, 3> coloumbLaw(const rtvector<MyUnits<T>, 3> &dstPos, MyUnits<int> dstCharge, const rtvector<MyUnits<T>, 3>& srcPos, MyUnits<int> srcCharge)
{
    nvAssert(MyUnitsTest::wasTested());

    rtvector<T, 3> vDir = removeUnits(dstPos) - removeUnits(srcPos);
    T fDistSqr = lengthSquared(vDir);
    nvAssert(fDistSqr > 1e-15);
    T fDist = sqrt(fDistSqr);
    nvAssert(dstCharge.m_value * srcCharge.m_value != 0); // performance warning
    T fForce = dstCharge.m_value * srcCharge.m_value / fDistSqr * ((ELECTRON_CHARGE / (TO_METERS * TO_METERS * TO_METERS)) * (COLOUMB_CONSTANT * TO_SECONDS) * TO_SECONDS * (ELECTRON_CHARGE / TO_KILOGRAMS));
    rtvector<T, 3> vForce = vDir * (fForce / fDist);

    return setUnits<MyUnits<T>>(vForce);
}

template <class T>
T coloumbLaw(T fCharge1, T fCharge2, T fDist, T fDistSQR) { return fCharge1 * fCharge2 * COLOUMB_CONSTANT / fDistSQR; }

template <class T>
rtvector<MyUnits<T>, 3> newtonLaw(MyUnits<T> fMass, const rtvector<MyUnits<T>, 3>& vForce)
{
    nvAssert(MyUnitsTest::wasTested());
    nvAssert(fMass.m_value > 0);
    return vForce / fMass;
}

template <class T>
T newtonLaw(T fMass, T fForce)
{
    return fForce / fMass;
}