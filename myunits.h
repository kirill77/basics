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
X operator *=(X s) { MyUnits<T>::m_value *= s.m_value; return *this; } \
X operator *=(T s) { MyUnits<T>::m_value *= s; return *this; } \
X operator -=(X s) { MyUnits<T>::m_value -= s.m_value; return *this; } \
X operator /=(T s) { MyUnits<T>::m_value /= s; return *this; } \
X sqrt() const { return X((T)::sqrt(this->m_value)); } \
X operator -(X s) const { return X(this->m_value - s.m_value); } \
X operator +(X s) const { return X(this->m_value + s.m_value); } \
X operator /(X s) const { return X(this->m_value / s.m_value); } \
X operator /(T s) const { return X(this->m_value / s); } \
X operator *(X s) const { return X(this->m_value * s.m_value); } \
X operator *(T s) const { return X(this->m_value * s); } \
bool operator == (X s) const { return this->m_value == s.m_value; } \
bool operator == (T s) const { return this->m_value == s; } \
bool operator != (X s) const { return this->m_value != s.m_value; } \
bool operator != (T s) const { return this->m_value != s; } \
bool operator >= (X s) const { return this->m_value >= s.m_value; } \
bool operator <= (X s) const { return this->m_value <= s.m_value; } \
bool operator <= (T s) const { return this->m_value <= s; } \
bool operator < (X s) const { return this->m_value < s.m_value; } \
bool operator < (T s) const { return this->m_value < s; } \
bool operator > (X s) const { return this->m_value > s.m_value; } \
bool operator > (T s) const { return this->m_value > s; }

#define COLOUMB_DEC(X) 8.9875517923##X
#define COLOUMB_CONSTANT COLOUMB_DEC(e9)
#define ELECTRON_DEC(X) 1.602176634##X
#define ELECTRON_CHARGE ELECTRON_DEC(e-19)
#define DALTON_DEC(X) 1.660539066605##X
#define DALTON DALTON_DEC(e-27)
#define AVOGADRO 6.02214076e23

// those constants define base units used at atomic scales
#define TO_GRAMS 1e-17
#define TO_KILOGRAMS 1e-20
#define TO_METERS 1e-8
#define TO_COLOUMBS (ELECTRON_CHARGE)
constexpr double TO_SECONDS = 1e-9;
constexpr double TO_NANOSECONDS = TO_SECONDS * 1e9;
constexpr double TO_FEMTOSECONDS = TO_SECONDS * 1e15;
constexpr double oaue1 = TO_NANOSECONDS / TO_SECONDS;
static_assert(oaue1 > 0.999999999e9 && oaue1 < 1.00000001e9, "Contradictory time defines");
constexpr double oaue2 = TO_FEMTOSECONDS / TO_SECONDS;
static_assert(oaue2 > 0.99999999e15 && oaue2 < 1.0000001e15, "Contradictory time defines");

template <class _T>
struct MyUnits
{
    typedef _T T;
    T m_value;

    UNIT_MEMBERS(MyUnits);
    template <class T1> MyUnits<T> operator +=(MyUnits<T1> s) { m_value += (T)s.m_value; return *this; }
    void clear() { m_value = 0; }

    //*** constructors
    static MyUnits<T> second() { return MyUnits<T>((T)(1. / TO_SECONDS)); }
    static MyUnits<T> milliSecond() { return MyUnits<T>((T)(1e-3 / TO_SECONDS)); }
    static MyUnits<T> nanoSecond() { return MyUnits<T>((T)(1e-6 / TO_SECONDS)); }
    static MyUnits<T> meter() { return MyUnits<T>((T)(1. / TO_METERS)); }
    static MyUnits<T> milliLiter() { return MyUnits<T>((T)(1e-6 / TO_METERS / TO_METERS / TO_METERS)); }
    static MyUnits<T> microMeter() { return MyUnits<T>((T)(1e-6 / TO_METERS)); }
    static MyUnits<T> nanoMeter() { return MyUnits<T>((T)(1e-9 / TO_METERS)); }
    static MyUnits<T> angstrom() { return MyUnits<T>((T)(1e-10 / TO_METERS)); }
    static MyUnits<T> picometer() { return MyUnits<T>((T)(1e-12 / TO_METERS)); }
    static MyUnits<T> electron() { return MyUnits<T>((T)(ELECTRON_CHARGE / TO_COLOUMBS)); }
    static MyUnits<T> dalton() { return MyUnits<T>((T)(DALTON / TO_KILOGRAMS)); }
    static MyUnits<T> gram() { return MyUnits<T>((T)(1e-3 / TO_KILOGRAMS)); }
    static MyUnits<T> kJperMole() { return MyUnits<T>((T)(1e3 / AVOGADRO / TO_KILOGRAMS / TO_METERS / TO_METERS * TO_SECONDS * TO_SECONDS)); }
    static MyUnits<T> joule() { return MyUnits<T>((T)(1. / TO_KILOGRAMS / TO_METERS / TO_METERS * TO_SECONDS * TO_SECONDS)); }
    static MyUnits<T> evalTemperature(MyUnits<T> fAvgKinEnergy) { nvAssert(MyUnitsTest::wasTested()); return fAvgKinEnergy; }
    static MyUnits<T> evalPressure(MyUnits<T> fTtlKinEnergy, MyUnits<T> fVolume, NvU32 nParticles) { nvAssert(MyUnitsTest::wasTested()); return fTtlKinEnergy / fVolume; }

    // convertion between kinetic energy and temperature
    T toCelcius() const { return (T)(this->m_value * (57206.3436653589e-18 / TO_SECONDS / TO_SECONDS) - 273.15); }
    static MyUnits<T> fromCelcius(T fValue) { return MyUnits<T>((fValue + 273.15) / (57206.3436653589e-18 / TO_SECONDS / TO_SECONDS)); }

    //*** base units
    T toKilograms() { return (T)(this->m_value * TO_KILOGRAMS); }
    T toMeters() { return (T)(this->m_value * TO_METERS); }
    T toAngstroms() { return (T)(this->m_value * (TO_METERS * 1e10)); }
    T toSeconds() { return (T)this->m_value * TO_SECONDS; }
    T toColoumbs() { return (T)(this->m_value * TO_COLOUMBS); }

    //*** derived units
    T toNewtons() const { return (T)(this->m_value * (TO_KILOGRAMS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toNanoseconds() const { return (T)(this->m_value * TO_NANOSECONDS); }
    T toFemtoseconds() const { return (T)(this->m_value * TO_FEMTOSECONDS); }
    T toMetersPerSecond2() const { return (T)(this->m_value * (TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toJoules() const { return (T)(this->m_value * (TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toKJperMole() const { return (T)(this->m_value * (1e-3 * AVOGADRO * TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toAtmospheres() const { return (T)this->m_value * 1e-10; }
};

#undef UNIT_MEMBERS

namespace std
{
    template <class T>
    inline MyUnits<T> abs(const MyUnits<T>& s) { return MyUnits<T>(std::abs(s.m_value)); }
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

template <class T1, class T2, class T3>
rtvector<MyUnits<T1>, 3> coloumbLaw(const rtvector<MyUnits<T1>, 3> &dstPos, MyUnits<T2> dstCharge, const rtvector<MyUnits<T1>, 3>& srcPos, MyUnits<T3> srcCharge)
{
    nvAssert(MyUnitsTest::wasTested());

    rtvector<T1, 3> vDir = removeUnits(dstPos) - removeUnits(srcPos);
    T1 fDistSqr = lengthSquared(vDir);
    nvAssert(fDistSqr > 1e-15);
    T1 fDist = sqrt(fDistSqr);
    nvAssert(dstCharge.m_value * srcCharge.m_value != 0); // performance warning
    T1 fForce = dstCharge.m_value * srcCharge.m_value / fDistSqr * ((ELECTRON_CHARGE / (TO_METERS * TO_METERS * TO_METERS)) * (COLOUMB_CONSTANT * TO_SECONDS) * TO_SECONDS * (ELECTRON_CHARGE / TO_KILOGRAMS));
    rtvector<T1, 3> vForce = vDir * (fForce / fDist);

    return setUnits<MyUnits<T1>>(vForce);
}

template <class T>
T coloumbLaw(T fCharge1, T fCharge2, T fDist, T fDistSQR) { return fCharge1 * fCharge2 * COLOUMB_CONSTANT / fDistSQR; }

template <class T1, class T2>
rtvector<MyUnits<T2>, 3> newtonLaw(MyUnits<T1> fMass, const rtvector<MyUnits<T2>, 3>& vForce)
{
    nvAssert(MyUnitsTest::wasTested());
    nvAssert(fMass.m_value > 0);
    return vForce / MyUnits<T2>(fMass.m_value);
}

template <class T>
T newtonLaw(T fMass, T fForce)
{
    return fForce / fMass;
}
