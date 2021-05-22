#include "myunits.h"

void MyUnitsTest::test()
{
    s_bTested = true;

    //**** test aboutEqual() function
    s_bTested = (1 != 1.000000000000001);
    nvAssert(s_bTested);
    s_bTested = s_bTested && aboutEqual<double>(1, 1.000000000000001);
    nvAssert(s_bTested);
    s_bTested = s_bTested && !aboutEqual<double>(1, 1.001);
    nvAssert(s_bTested);

    //**** test coloumb law
    // compute force in MyUnits
    MyUnits<double> deltaS = MyUnits<double>::microMeter() * 10;
    for (NvU32 uPass = 0; uPass < 3; ++uPass)
    {
        rtvector<MyUnits<double>, 3> p1, p2;
        p2[uPass] = p1[uPass] + deltaS;
        // since we usually represent electrons with MyUnits<int>, check that charge of one electron is == 1
        auto oneElectron = MyUnits<double>::electron();
        s_bTested = s_bTested && aboutEqual(oneElectron.m_value, 1.);
        nvAssert(s_bTested);
        rtvector<MyUnits<double>, 3> vMyForce = coloumbLaw(p1, MyUnits<int>::electron() * -1, p2, MyUnits<int>::electron() * -1);
        // compute force in SI units
        double fSIForce = coloumbLaw<double>(1.60217662e-19, 1.60217662e-19, deltaS.toMeters(), sqr(deltaS.toMeters()));
        // convert and compare
        double fMySIForce = -vMyForce[uPass].toNewtons();
        double fPercentDifference = std::abs(fSIForce - fMySIForce) / fSIForce * 100;
        s_bTested = s_bTested && (fPercentDifference < 0.001);
        nvAssert(s_bTested);

        //**** test newton law
        MyUnits<double> fMyMass = MyUnits<double>::dalton() * 22.989769282; // that must be sodium mass
        rtvector<MyUnits<double>, 3> vMyAcceleration = newtonLaw(fMyMass, vMyForce);
        double fSIAcceleration = newtonLaw(fMyMass.toKilograms(), vMyForce[uPass].toNewtons());
        double fMySIAcceleration = vMyAcceleration[uPass].toSIAcceleration();
        fPercentDifference = std::abs(fSIAcceleration - fMySIAcceleration) / fSIAcceleration * 100;
        s_bTested = s_bTested && (fPercentDifference < 0.001);
        nvAssert(s_bTested);

        deltaS /= 3;
    }
}
bool MyUnitsTest::s_bTested = false;
