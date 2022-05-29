#include <time.h>
#include "MonteCarlo/RNGUniform.h"
#include "bonds.h"

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
    MyUnits<double> deltaS = MyUnits1<double>::microMeter() * 10;
    for (NvU32 uPass = 0; uPass < 3; ++uPass)
    {
        rtvector<MyUnits<double>, 3> p1, p2;
        p2[uPass] = p1[uPass] + deltaS;
        // since we usually represent electrons with MyUnits<int>, check that charge of one electron is == 1
        auto oneElectron = MyUnits1<double>::electronCharge();
        s_bTested = s_bTested && aboutEqual(oneElectron, 1.);
        nvAssert(s_bTested);
        rtvector<MyUnits<double>, 3> vMyForce = coloumbLaw(p1, MyUnits1<double>::electronCharge() * -1, p2, MyUnits1<double>::electronCharge() * -1);
        // compute force in SI units
        double fSIForce = coloumbLawSI<double>(1.60217662e-19, 1.60217662e-19, MyUnits1<double>::toMeters(deltaS), sqr(MyUnits1<double>::toMeters(deltaS)));
        // convert and compare
        double fMySIForce = MyUnits1<double>::toNewtons(-vMyForce[uPass]);
        double fPercentDifference = std::abs(fSIForce - fMySIForce) / fSIForce * 100;
        s_bTested = s_bTested && (fPercentDifference < 0.001);
        nvAssert(s_bTested);

        //**** test newton law
        MyUnits<double> fMyMass = MyUnits1<double>::dalton() * 22.989769282; // that must be sodium mass
        rtvector<MyUnits<double>, 3> vMyAcceleration = newtonLaw(fMyMass, vMyForce);
        double fSIAcceleration = newtonLaw(MyUnits1<double>::toKilograms(fMyMass), MyUnits1<double>::toNewtons(vMyForce[uPass]));
        double fMySIAcceleration = MyUnits1<double>::toMetersPerSecond2(vMyAcceleration[uPass]);
        fPercentDifference = std::abs(fSIAcceleration - fMySIAcceleration) / fSIAcceleration * 100;
        s_bTested = s_bTested && (fPercentDifference < 0.001);
        nvAssert(s_bTested);

        deltaS /= 3;
    }

    BondsDataBase<double>::init();

    // check that we can translate between our energy and SI energy without problems
    {
        MyUnits<double> fMass = MyUnits1<double>::dalton() * 100;
        MyUnits<double> fDeltaTime = MyUnits1<double>::nanoSecond() * 100;
        MyUnits<double> fDeltaDistance = MyUnits1<double>::nanoMeter() * 100;
        MyUnits<double> fExternalPotential = MyUnits1<double>::joule() * 1e-10;
        MyUnits<double> fSpeed = fDeltaDistance / fDeltaTime;

        MyUnits<double> fTotalEnergy = fMass * sqr(fSpeed) / 2 + fExternalPotential;
        double fSITotalEnergy = MyUnits1<double>::toKilograms(fMass) * sqr(MyUnits1<double>::toMeters(fDeltaDistance) / MyUnits1<double>::toSeconds(fDeltaTime)) / 2 + MyUnits1<double>::toJoules(fExternalPotential);
        double fMySITotalEnergy = MyUnits1<double>::toJoules(fTotalEnergy);
        double fPercentDifference = std::abs(fMySITotalEnergy - fSITotalEnergy) / fSITotalEnergy * 100;
        s_bTested = s_bTested && fPercentDifference < 1;
        nvAssert(s_bTested);
    }

    {
        // must have data for H-H
        auto& eBond = BondsDataBase<double>::getEBond(1, 1, 1);
        s_bTested = s_bTested && eBond.isValid();
        nvAssert(s_bTested);
    }
    {
        // must have data for H-O
        auto& eBond = BondsDataBase<double>::getEBond(1, 8, 1);
        s_bTested = s_bTested && eBond.isValid();
        nvAssert(s_bTested);
    }

    //*** use lennard-Jones force to figure out bond lengths and bond energies
    {
        RNGUniform rng((NvU32)time(nullptr));
        const auto& aBonds = BondsDataBase<double>::getABonds();
        for (auto it = aBonds.begin(); it != aBonds.end(); ++it)
        {
            const auto& aBond = it->second;
            for (NvU32 nElectrons = 1; nElectrons < MAX_ELECTRONS_PER_BOND; ++nElectrons)
            {
                const auto& eBond = aBond[nElectrons];
                if (!eBond.isValid())
                    continue;
                NvU32 uDim = rng.generateUnsigned(0, 3);
                rtvector<MyUnits<double>, 3> vPos[2];
                MyUnits<double> fPrevForce(-1);
                auto deltaX = MyUnits1<double>::angstrom() * 0.005;
                vPos[1][uDim] = deltaX;
                double fEnergy = 0, fPotential = 0;
                for ( ; vPos[1][uDim] < MyUnits1<double>::angstrom() * 10; vPos[1][uDim] += deltaX)
                {
                    BondsDataBase<double>::LJ_Out ljOut;
                    rtvector<MyUnits<double>, 3> vDir = vPos[0] - vPos[1];
                    MyUnits<double> fDistSqr = dot(vDir, vDir);
                    bool bNonZero = eBond.lennardJones(fDistSqr, ljOut);
                    if (bNonZero)
                    {
                        MyUnits<double> fForce = vDir[uDim] * (ljOut.fForceTimesR / fDistSqr);
                        if (fForce * fPrevForce <= 0.) // we go until the force changes sign
                        {
                            MyUnits<double> bondLength = vPos[1][uDim];
                            double fPercentDifference = std::abs(bondLength - eBond.getLength()) / eBond.getLength() * 100;
                            s_bTested = s_bTested && fPercentDifference < 1;
                            nvAssert(s_bTested);
                            fEnergy = 0;
                            fPotential = ljOut.fPotential;
                        }
                        else
                        {
                            fEnergy += (fPrevForce + fForce) / 2;
                        }
                        fPrevForce = fForce;
                    }
                    else
                    {
                        fPrevForce = 0;
                    }
                }
                fEnergy *= deltaX;
                double fPercentDifference = std::abs(fEnergy - eBond.getEnergy()) / eBond.getEnergy() * 100;
                s_bTested = s_bTested && fPercentDifference < 1;
                nvAssert(s_bTested);
                // potential is negative - that's why plus instead of minus here
                fPercentDifference = std::abs(fPotential + eBond.getEnergy()) / eBond.getEnergy() * 100;
                s_bTested = s_bTested && fPercentDifference < 1;
                nvAssert(s_bTested);
            }
        }
    }

    // check that our temperature and pressure formulas match reality
    {
        MyUnits<double> fMass = BondsDataBase<double>::getElement(NPROTONS_O).m_fMass + BondsDataBase<double>::getElement(NPROTONS_H).m_fMass * 2;
        // the mass of water molecule must be about 18.01528 g/mol - check that
        MyUnits<double> fMass1 = MyUnits1<double>::gram() * 18.01528 / AVOGADRO;
        s_bTested = s_bTested && aboutEqual<double>(fMass, fMass1, 0.1);
        nvAssert(s_bTested);
        // https://www.verticallearning.org/curriculum/science/gr7/student/unit01/page05.html
        // we know "approximate" speed of water molecules must be at different temperatures
        {
            MyUnits<double> fSpeedC0 = MyUnits1<double>::meter() * 565 / MyUnits1<double>::second();
            MyUnits<double> fKinEnergy = fMass1 * fSpeedC0 * fSpeedC0 / 2;
            MyUnits<double> fTempC0 = MyUnits1<double>::evalTemperature(fKinEnergy);
            double f = MyUnits1<double>::toCelcius(fTempC0);
            s_bTested = s_bTested && aboutEqual<double>(f, 0, 0.1);
            nvAssert(s_bTested);
        }
        {
            MyUnits<double> fSpeedC20 = MyUnits1<double>::meter() * 590 / MyUnits1<double>::second();
            auto fTempC20 = MyUnits1<double>::evalTemperature(fMass1 * fSpeedC20 * fSpeedC20 / 2);
            double f = MyUnits1<double>::toCelcius(fTempC20);
            s_bTested = s_bTested && aboutEqual<double>(f, 24.707357663090249, 0.1);
            nvAssert(s_bTested);
        }
        {
            MyUnits<double> fSpeedC100 = MyUnits1<double>::meter() * 660 / MyUnits1<double>::second();
            auto fTempC100 = MyUnits1<double>::evalTemperature(fMass1 * fSpeedC100 * fSpeedC100 / 2);
            double f = MyUnits1<double>::toCelcius(fTempC100);
            s_bTested = s_bTested && aboutEqual<double>(f, 99.578138460333548, 0.1);
            nvAssert(s_bTested);
        }
    }

    typedef double T;

    // test for correctness of constants for computing kinetic and potential energy of an electron - check that energy is conserved
    // as particle is accelerating in electric field
    // first compute in SI units
    double fDisparity = 0;
    {
        double fInitialDistance = 1;
        double fCurDistance = fInitialDistance;
        double fInitialPotentialEnergy = chargePotentialEnergySI(1., 1., fCurDistance);
        double fCurSpeed = 0;
        while (fCurDistance * 2 > fInitialDistance)
        {
            double fCurForce = coloumbLawSI(1., 1., fCurDistance, fCurDistance * fCurDistance);
            double fCurAcceleration = fCurForce;
            double fWantedDeltaDistance = fCurDistance * 0.0001;
            double fDeltaTime = sqrt(fWantedDeltaDistance * 2 / fCurAcceleration);
            double fNextSpeed = fCurSpeed + fCurAcceleration * fDeltaTime;
            double fAverageSpeed = (fNextSpeed + fCurSpeed) / 2;
            double fDeltaDistance = fAverageSpeed * fDeltaTime;
            double fNextDistance = fCurDistance - fDeltaDistance;
            nvAssert(fNextDistance < fCurDistance);
            double fNextPotentialEnergy = chargePotentialEnergySI(1., 1., fNextDistance);
            double fNextKineticEnergy = fNextSpeed * fNextSpeed / 2;
            double fNextTotalEnergy = fNextPotentialEnergy + fNextKineticEnergy;

            fDisparity = ((fNextTotalEnergy - fInitialPotentialEnergy) / (fNextTotalEnergy + fInitialPotentialEnergy));
            s_bTested = (s_bTested && fDisparity < 0.01);
            nvAssert(s_bTested);

            fCurDistance = fNextDistance;
            fCurSpeed = fNextSpeed;
        }
        fDisparity = fDisparity;
    }
    // now compute in my units
    {
        T fInitialDistance = MyUnits1<T>::angstrom();
        T fCurDistance = fInitialDistance;
        T fInitialPotentialEnergy = chargePotentialEnergy(MyUnits1<T>::electronCharge(), MyUnits1<T>::electronCharge(), fCurDistance);
        T fCurSpeed = 0;
        while (fCurDistance * 2 > fInitialDistance)
        {
            T fCurForce = coloumbLaw(MyUnits1<T>::electronCharge(), MyUnits1<T>::electronCharge(), fCurDistance, fCurDistance * fCurDistance);
            T fCurAcceleration = fCurForce / MyUnits1<T>::electronMass();
            MyUnits<T> fWantedDeltaDistance = fCurDistance * 0.0001;
            MyUnits<T> fDeltaTime = sqrt(fWantedDeltaDistance * 2 / fCurAcceleration);
            MyUnits<T> fNextSpeed = fCurSpeed + fCurAcceleration * fDeltaTime;
            MyUnits<T> fAverageSpeed = (fNextSpeed + fCurSpeed) / 2;
            MyUnits<T> fDeltaDistance = fAverageSpeed * fDeltaTime;
            MyUnits<T> fNextDistance = fCurDistance - fDeltaDistance;
            nvAssert(fNextDistance < fCurDistance);
            T fNextPotentialEnergy = chargePotentialEnergy(MyUnits1<T>::electronCharge(), MyUnits1<T>::electronCharge(), fNextDistance);
            T fNextKineticEnergy = MyUnits1<T>::electronMass() * fNextSpeed * fNextSpeed / 2;
            MyUnits<T> fNextTotalEnergy = fNextPotentialEnergy + fNextKineticEnergy;

            fDisparity = ((fNextTotalEnergy - fInitialPotentialEnergy) / (fNextTotalEnergy + fInitialPotentialEnergy));
            s_bTested = (s_bTested && fDisparity < 0.01);
            nvAssert(s_bTested);

            fCurDistance = fNextDistance;
            fCurSpeed = fNextSpeed;
        }
        fDisparity = fDisparity;
    }
}
bool MyUnitsTest::s_bTested = false;

std::unordered_map<ATOM_KEY, BondsDataBase<double>::ABond> BondsDataBase<double>::m_aBonds;
std::unordered_map<NvU32, BondsDataBase<double>::Element> BondsDataBase<double>::m_elements;
MyUnits<double> BondsDataBase<double>::s_zeroForceDist = MyUnits1<double>::angstrom() * 6;
MyUnits<double> BondsDataBase<double>::s_zeroForceDistSqr = BondsDataBase<double>::s_zeroForceDist * BondsDataBase<double>::s_zeroForceDist;

template <class T>
void BondsDataBase<T>::init()
{
    setAtom(NPROTONS_O, MyUnits1<T>::dalton() * 15.999, MyUnits1<T>::picometer() * 152, 3.44, 2 /*valence*/);
    setAtom(NPROTONS_H, MyUnits1<T>::dalton() * 1.008, MyUnits1<T>::picometer() * 120, 2.2, 1 /*valence*/);

    setBond(NPROTONS_O, NPROTONS_O, 1, MyUnits1<T>::angstrom() * 1.278, MyUnits1<T>::kJperMole() * 140);
    setBond(NPROTONS_O, NPROTONS_O, 2, MyUnits1<T>::angstrom() * 1.2, MyUnits1<T>::kJperMole() * 498);
    setBond(NPROTONS_O, NPROTONS_H, 1, MyUnits1<T>::angstrom() * 1.2, MyUnits1<T>::kJperMole() * 464);
    setBond(NPROTONS_H, NPROTONS_H, 1, MyUnits1<T>::angstrom() * 0.74, MyUnits1<T>::kJperMole() * 436);
}
template <class T>
void BondsDataBase<T>::setBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons, MyUnits<T> fBondLength, MyUnits<T> fBondEnergy)
{
    auto& eBond = accessEBond(nProtons1, nProtons2, nElectrons);
    eBond = EBond(fBondLength, fBondEnergy);
    accessEBond(nProtons2, nProtons1, nElectrons) = eBond;
}
template <class T>
void BondsDataBase<T>::setAtom(NvU32 nProtons, MyUnits<T> fMass, MyUnits<T> fRadius, T fElectroNegativity, NvU32 uValence)
{
    auto& element = m_elements[nProtons];
    element.m_fMass = fMass;
    element.m_fRadius = fRadius;
    element.m_fElectroNegativity = fElectroNegativity;
    element.m_uValence = uValence;
}
