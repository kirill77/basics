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

static constexpr double COLOUMB_CONSTANT = 8.987551792314e9; // kg * m^3 / C / s^2
static constexpr double ELECTRON_CHARGE = 1.602176634e-19; // C
static constexpr double DALTON = 1.660539066605e-27;
static constexpr double AVOGADRO = 6.02214076e23;
static constexpr double PLANCK_CONSTANT = 6.626070041e-34; // m^2 * kg / s

// those constants define base units used at atomic scales
static constexpr double TO_KILOGRAMS = 1e-20;
static constexpr double TO_METERS = 1e-8;
static constexpr double TO_COLOUMB = ELECTRON_CHARGE;
constexpr double TO_SECONDS = 1e-9;
constexpr double TO_NANOSECONDS = TO_SECONDS * 1e9;
constexpr double TO_FEMTOSECONDS = TO_SECONDS * 1e15;

constexpr double oaue1 = TO_NANOSECONDS / TO_SECONDS;
static_assert(oaue1 > 0.999999999e9 && oaue1 < 1.00000001e9, "Contradictory time defines");
constexpr double oaue2 = TO_FEMTOSECONDS / TO_SECONDS;
static_assert(oaue2 > 0.99999999e15 && oaue2 < 1.0000001e15, "Contradictory time defines");

static constexpr double MY_PLANCK_CONSTANT = PLANCK_CONSTANT / TO_METERS / TO_METERS / TO_KILOGRAMS * TO_SECONDS;

template <class _T>
struct MyUnits
{
    typedef _T T;
    T m_value;

    MyUnits() { m_value = 0; }
    explicit MyUnits(T value) : m_value(value) { }
    MyUnits<_T> operator -() const { return MyUnits<_T>(-m_value); }
    MyUnits<_T> operator *=(MyUnits<_T> s) { m_value *= s.m_value; return *this; }
    MyUnits<_T> operator *=(T s) { m_value *= s; return *this; }
    MyUnits<_T> operator -=(MyUnits<_T> s) { m_value -= s.m_value; return *this; }
    MyUnits<_T> operator /=(T s) { m_value /= s; return *this; }
    MyUnits<_T> operator -(MyUnits<_T> s) const { return MyUnits<_T>(m_value - s.m_value); }
    MyUnits<_T> operator +(MyUnits<_T> s) const { return MyUnits<_T>(m_value + s.m_value); }
    MyUnits<_T> operator /(MyUnits<_T> s) const { return MyUnits<_T>(m_value / s.m_value); }
    MyUnits<_T> operator /=(MyUnits<_T> s) { return MyUnits<_T>(m_value /= s.m_value); }
    MyUnits<_T> operator /(T s) const { return MyUnits<_T>(m_value / s); }
    MyUnits<_T> operator *(MyUnits<_T> s) const { return MyUnits<_T>(m_value * s.m_value); }
    MyUnits<_T> operator *(T s) const { return MyUnits<_T>(m_value * s); }
    bool operator == (MyUnits<_T> s) const { return m_value == s.m_value; }
    bool operator == (T s) const { return m_value == s; }
    bool operator != (MyUnits<_T> s) const { return m_value != s.m_value; }
    bool operator != (T s) const { return m_value != s; }
    bool operator >= (MyUnits<_T> s) const { return m_value >= s.m_value; }
    bool operator <= (MyUnits<_T> s) const { return m_value <= s.m_value; }
    bool operator <= (T s) const { return m_value <= s; }
    bool operator < (MyUnits<_T> s) const { return m_value < s.m_value; }
    bool operator < (T s) const { return m_value < s; }
    bool operator > (MyUnits<_T> s) const { return m_value > s.m_value; }
    bool operator > (T s) const { return m_value > s; }

    template <class T1> MyUnits<T> operator +=(MyUnits<T1> s) { m_value += (T)s.m_value; return *this; }

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
    static MyUnits<T> electronCharge() { return MyUnits<T>((T)(ELECTRON_CHARGE / TO_COLOUMB)); }
    static MyUnits<T> dalton() { return MyUnits<T>((T)(DALTON / TO_KILOGRAMS)); }
    static MyUnits<T> electronMass() { return MyUnits<T>((T)(9.1093837015e-31 / TO_KILOGRAMS)); }
    static MyUnits<T> gram() { return MyUnits<T>((T)(1e-3 / TO_KILOGRAMS)); }
    static MyUnits<T> kJperMole() { return MyUnits<T>((T)(1e3 / AVOGADRO / TO_KILOGRAMS / TO_METERS / TO_METERS * TO_SECONDS * TO_SECONDS)); }
    static MyUnits<T> joule() { return MyUnits<T>((T)(1. / TO_KILOGRAMS / TO_METERS / TO_METERS * TO_SECONDS * TO_SECONDS)); }
    static MyUnits<T> evalTemperature(MyUnits<T> fAvgKinEnergy) { nvAssert(MyUnitsTest::wasTested()); return fAvgKinEnergy; }
    static MyUnits<T> evalPressure(MyUnits<T> fTtlKinEnergy, MyUnits<T> fVolume, NvU32 nParticles) { nvAssert(MyUnitsTest::wasTested()); return fTtlKinEnergy * (2./3) / fVolume; }

    // convertion between kinetic energy and temperature
    T toCelcius() const { return (T)(m_value * (57206.3436653589e-2 * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS) - 273.15); }
    static MyUnits<T> fromCelcius(T fValue) { return MyUnits<T>((fValue + 273.15) / (57206.3436653589e-2 * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }

    //*** base units
    T toKilograms() { return (T)(m_value * TO_KILOGRAMS); }
    T toMeters() { return (T)(m_value * TO_METERS); }
    T toAngstroms() { return (T)(m_value * (TO_METERS * 1e10)); }
    T toSeconds() { return (T)m_value * TO_SECONDS; }
    T toColoumbs() { return (T)(m_value * TO_COLOUMB); }

    //*** derived units
    T toNewtons() const { return (T)(m_value * (TO_KILOGRAMS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toNanoseconds() const { return (T)(m_value * TO_NANOSECONDS); }
    T toFemtoseconds() const { return (T)(m_value * TO_FEMTOSECONDS); }
    T toMetersPerSecond2() const { return (T)(m_value * (TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toJoules() const { return (T)(m_value * (TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toKJperMole() const { return (T)(m_value * (1e-3 * AVOGADRO * TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    T toAtmospheres() const { return (T)(m_value * (TO_KILOGRAMS / TO_SECONDS / TO_SECONDS / TO_METERS / 101325)); }
};

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

template <class T>
MyUnits<T> coloumbLaw(MyUnits<T> fCharge1, MyUnits<T> fCharge2, MyUnits<T> fDist, MyUnits<T> fDistSQR)
{
    nvAssert(MyUnitsTest::wasTested());
    return fCharge1 * fCharge2 / fDistSQR * ((TO_COLOUMB / (TO_METERS * TO_METERS * TO_METERS)) * (COLOUMB_CONSTANT * TO_SECONDS) * TO_SECONDS * (TO_COLOUMB / TO_KILOGRAMS));
}

template <class T>
rtvector<MyUnits<T>, 3> coloumbLaw(const rtvector<MyUnits<T>, 3> &dstPos, MyUnits<T> dstCharge, const rtvector<MyUnits<T>, 3>& srcPos, MyUnits<T> srcCharge)
{
    nvAssert(MyUnitsTest::wasTested());

    auto vDir = dstPos - srcPos;
    MyUnits<T> fDistSqr = lengthSquared(vDir);
    nvAssert(fDistSqr > 1e-15);
    MyUnits<T> fDist = sqrt(fDistSqr);
    nvAssert(dstCharge.m_value * srcCharge.m_value != 0); // performance warning
    MyUnits<T> fForce = coloumbLaw<T>(dstCharge, srcCharge, fDist, fDistSqr);
    auto vForce = vDir * (fForce / fDist);

    return vForce;
}

template <class T>
T coloumbLawSI(T fCharge1, T fCharge2, T fDist, T fDistSQR)
{
    return fCharge1 * fCharge2 * COLOUMB_CONSTANT / fDistSQR;
}

template <class T>
T chargePotentialEnergySI(T fCharge1, T fCharge2, T fDist)
{
    nvAssert(MyUnitsTest::wasTested());
    return -fCharge1 * fCharge2 / fDist * COLOUMB_CONSTANT;
}
template <class T>
MyUnits<T> chargePotentialEnergy(MyUnits<T> fCharge1, MyUnits<T> fCharge2, MyUnits<T> fDist)
{
    nvAssert(MyUnitsTest::wasTested());
    return -fCharge1 * fCharge2 / fDist * ((TO_COLOUMB / (TO_METERS * TO_METERS * TO_METERS)) * (COLOUMB_CONSTANT * TO_SECONDS) * TO_SECONDS * (TO_COLOUMB / TO_KILOGRAMS));
}

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
