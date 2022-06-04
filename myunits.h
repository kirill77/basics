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
static constexpr double TO_KILOGRAMS = 1e-26;
static constexpr double TO_METERS = 1e-10;
static constexpr double TO_COLOUMB = ELECTRON_CHARGE;
constexpr double TO_SECONDS = 1e-14;
constexpr double TO_NANOSECONDS = TO_SECONDS * 1e9;
constexpr double TO_FEMTOSECONDS = TO_SECONDS * 1e15;

constexpr double oaue1 = TO_NANOSECONDS / TO_SECONDS;
static_assert(oaue1 > 0.999999999e9 && oaue1 < 1.00000001e9, "Contradictory time defines");
constexpr double oaue2 = TO_FEMTOSECONDS / TO_SECONDS;
static_assert(oaue2 > 0.99999999e15 && oaue2 < 1.0000001e15, "Contradictory time defines");

static constexpr double MY_PLANCK_CONSTANT = PLANCK_CONSTANT / TO_METERS / TO_METERS / TO_KILOGRAMS * TO_SECONDS;

// alias MyUnits<T> to just T
template <class T> using MyUnits = T;

template <class _T>
struct MyUnits1
{
    typedef _T T;

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
    static double evalPressure(MyUnits<T> fTtlKinEnergy, MyUnits<T> fVolume) { nvAssert(MyUnitsTest::wasTested()); return fTtlKinEnergy * (2./3) / fVolume; }

    // convertion between kinetic energy and temperature
    static T toCelcius(MyUnits<T> fValue) { return (T)(fValue * (57206.3436653589e18 * TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS) - 273.15); }
    static MyUnits<T> fromCelcius(T fValue) { return MyUnits<T>((fValue + 273.15) / (57206.3436653589e18 * TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }

    //*** base units
    static T toKilograms(MyUnits<T> fValue) { return (T)(fValue * TO_KILOGRAMS); }
    static T toMeters(MyUnits<T> fValue) { return (T)(fValue * TO_METERS); }
    static T toAngstroms(MyUnits<T> fValue) { return (T)(fValue * (TO_METERS * 1e10)); }
    static T toSeconds(MyUnits<T> fValue) { return (T)fValue * TO_SECONDS; }
    static T toColoumbs(MyUnits<T> fValue) { return (T)(fValue * TO_COLOUMB); }

    //*** derived units
    static T toNewtons(MyUnits<T> fValue) { return (T)(fValue * (TO_KILOGRAMS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    static T toNanoseconds(MyUnits<T> fValue) { return (T)(fValue * TO_NANOSECONDS); }
    static T toFemtoseconds(MyUnits<T> fValue) { return (T)(fValue * TO_FEMTOSECONDS); }
    static T toMetersPerSecond2(MyUnits<T> fValue) { return (T)(fValue * (TO_METERS / TO_SECONDS / TO_SECONDS)); }
    static T toJoules(MyUnits<T> fValue) { return (T)(fValue * (TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    static T toKJperMole(MyUnits<T> fValue) { return (T)(fValue * (1e-3 * AVOGADRO * TO_KILOGRAMS * TO_METERS * TO_METERS / TO_SECONDS / TO_SECONDS)); }
    static T toAtmospheres(MyUnits<T> fValue) { return (T)(fValue * (TO_KILOGRAMS / TO_SECONDS / TO_SECONDS / TO_METERS / 101325)); }
};

template <class T> using MyUnits3 = rtvector<MyUnits<T>, 3>;

namespace std
{
    template <class T>
    inline MyUnits<T> abs(const MyUnits<T>& s) { return MyUnits<T>(std::abs(s)); }
};

template <class T> MyUnits<T> sqrt(MyUnits<T> s) { return MyUnits<T>(::sqrt(s)); }
template <class T> MyUnits<T> pow(MyUnits<T> s, T p) { return MyUnits<T>(::pow(s, p)); }

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
    nvAssert(dstCharge * srcCharge != 0); // performance warning
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
    nvAssert(fMass > 0);
    return vForce / MyUnits<T2>(fMass);
}

template <class T>
T newtonLaw(T fMass, T fForce)
{
    return fForce / fMass;
}
