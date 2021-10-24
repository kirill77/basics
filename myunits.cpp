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
    MyUnits<double> deltaS = MyUnits<double>::microMeter() * 10;
    for (NvU32 uPass = 0; uPass < 3; ++uPass)
    {
        rtvector<MyUnits<double>, 3> p1, p2;
        p2[uPass] = p1[uPass] + deltaS;
        // since we usually represent electrons with MyUnits<int>, check that charge of one electron is == 1
        auto oneElectron = MyUnits<double>::electronCharge();
        s_bTested = s_bTested && aboutEqual(oneElectron.m_value, 1.);
        nvAssert(s_bTested);
        rtvector<MyUnits<double>, 3> vMyForce = coloumbLaw(p1, MyUnits<double>::electronCharge() * -1, p2, MyUnits<double>::electronCharge() * -1);
        // compute force in SI units
        double fSIForce = coloumbLawSI<double>(1.60217662e-19, 1.60217662e-19, deltaS.toMeters(), sqr(deltaS.toMeters()));
        // convert and compare
        double fMySIForce = -vMyForce[uPass].toNewtons();
        double fPercentDifference = std::abs(fSIForce - fMySIForce) / fSIForce * 100;
        s_bTested = s_bTested && (fPercentDifference < 0.001);
        nvAssert(s_bTested);

        //**** test newton law
        MyUnits<double> fMyMass = MyUnits<double>::dalton() * 22.989769282; // that must be sodium mass
        rtvector<MyUnits<double>, 3> vMyAcceleration = newtonLaw(fMyMass, vMyForce);
        double fSIAcceleration = newtonLaw(fMyMass.toKilograms(), vMyForce[uPass].toNewtons());
        double fMySIAcceleration = vMyAcceleration[uPass].toMetersPerSecond2();
        fPercentDifference = std::abs(fSIAcceleration - fMySIAcceleration) / fSIAcceleration * 100;
        s_bTested = s_bTested && (fPercentDifference < 0.001);
        nvAssert(s_bTested);

        deltaS /= 3;
    }

    BondsDataBase<double>::init();

    // check that we can translate between our energy and SI energy without problems
    {
        MyUnits<double> fMass = MyUnits<double>::dalton() * 100;
        MyUnits<double> fDeltaTime = MyUnits<double>::nanoSecond() * 100;
        MyUnits<double> fDeltaDistance = MyUnits<double>::nanoMeter() * 100;
        MyUnits<double> fExternalPotential = MyUnits<double>::joule() * 1e-10;
        MyUnits<double> fSpeed = fDeltaDistance / fDeltaTime;

        MyUnits<double> fTotalEnergy = fMass * sqr(fSpeed) / 2 + fExternalPotential;
        double fSITotalEnergy = fMass.toKilograms() * sqr(fDeltaDistance.toMeters() / fDeltaTime.toSeconds()) / 2 + fExternalPotential.toJoules();
        double fMySITotalEnergy = fTotalEnergy.toJoules();
        double fPercentDifference = std::abs(fMySITotalEnergy - fSITotalEnergy) / fSITotalEnergy * 100;
        s_bTested = s_bTested && fPercentDifference < 1;
        nvAssert(s_bTested);
    }

    {
        // must have data for H-H
        auto& eBond = BondsDataBase<double>::getEBond(1, 1, 1);
        s_bTested = s_bTested && eBond.m_fEnergy > 0;
        nvAssert(s_bTested);
    }
    {
        // must have data for H-O
        auto& eBond = BondsDataBase<double>::getEBond(1, 8, 1);
        s_bTested = s_bTested && eBond.m_fEnergy > 0;
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
                if (eBond.m_fEnergy == 0) // invalid eBond?
                    continue;
                NvU32 uDim = rng.generateUnsigned(0, 3);
                rtvector<MyUnits<double>, 3> vPos[2];
                MyUnits<double> fPrevForce(-1);
                auto deltaX = MyUnits<double>::angstrom() * 0.005;
                vPos[1][uDim] = deltaX;
                MyUnits<double> fEnergy, fPotential;
                for ( ; vPos[1][uDim] < MyUnits<double>::angstrom() * 10; vPos[1][uDim] += deltaX)
                {
                    BondsDataBase<double>::LJ_Out ljOut;
                    ljOut.vForce = vPos[0] - vPos[1];
                    bool bNonZero = eBond.lennardJones(ljOut.vForce, ljOut);
                    if (bNonZero)
                    {
                        if (ljOut.vForce[uDim] * fPrevForce <= 0.)
                        {
                            MyUnits<double> bondLength = vPos[1][uDim];
                            double fPercentDifference = std::abs(bondLength.m_value - eBond.m_fLength.m_value) / eBond.m_fLength.m_value * 100;
                            s_bTested = s_bTested && fPercentDifference < 1;
                            nvAssert(s_bTested);
                            fEnergy.clear();
                            fPotential = ljOut.fPotential;
                        }
                        else
                        {
                            fEnergy += (fPrevForce + ljOut.vForce[uDim]) / 2;
                        }
                        fPrevForce = ljOut.vForce[uDim];
                    }
                    else
                    {
                        fPrevForce.clear();
                    }
                }
                fEnergy *= deltaX;
                double fPercentDifference = std::abs(fEnergy.m_value - eBond.m_fEnergy.m_value) / eBond.m_fEnergy.m_value * 100;
                s_bTested = s_bTested && fPercentDifference < 1;
                nvAssert(s_bTested);
                // potential is negative - that's why plus instead of minus here
                fPercentDifference = std::abs(fPotential.m_value + eBond.m_fEnergy.m_value) / eBond.m_fEnergy.m_value * 100;
                s_bTested = s_bTested && fPercentDifference < 1;
                nvAssert(s_bTested);
            }
        }
    }

    // check that our temperature and pressure formulas match reality
    {
        MyUnits<double> fMass = BondsDataBase<double>::getAtom(NPROTONS_O).m_fMass + BondsDataBase<double>::getAtom(NPROTONS_H).m_fMass * 2;
        // the mass of water molecule must be about 18.01528 g/mol - check that
        MyUnits<double> fMass1 = MyUnits<double>::gram() * 18.01528 / AVOGADRO;
        s_bTested = s_bTested && aboutEqual<double>(fMass.m_value, fMass1.m_value, 0.1);
        nvAssert(s_bTested);
        // https://www.verticallearning.org/curriculum/science/gr7/student/unit01/page05.html
        // we know "approximate" speed of water molecules must be at different temperatures
        {
            MyUnits<double> fSpeedC0 = MyUnits<double>::meter() * 565 / MyUnits<double>::second();
            MyUnits<double> fKinEnergy = fMass1 * fSpeedC0 * fSpeedC0 / 2;
            MyUnits<double> fTempC0 = MyUnits<double>::evalTemperature(fKinEnergy);
            double f = fTempC0.toCelcius();
            s_bTested = s_bTested && aboutEqual<double>(f, 0, 0.1);
            nvAssert(s_bTested);
        }
        {
            MyUnits<double> fSpeedC20 = MyUnits<double>::meter() * 590 / MyUnits<double>::second();
            auto fTempC20 = MyUnits<double>::evalTemperature(fMass1 * fSpeedC20 * fSpeedC20 / 2);
            double f = fTempC20.toCelcius();
            s_bTested = s_bTested && aboutEqual<double>(f, 24.707357663090249, 0.1);
            nvAssert(s_bTested);
        }
        {
            MyUnits<double> fSpeedC100 = MyUnits<double>::meter() * 660 / MyUnits<double>::second();
            auto fTempC100 = MyUnits<double>::evalTemperature(fMass1 * fSpeedC100 * fSpeedC100 / 2);
            double f = fTempC100.toCelcius();
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
        MyUnits<T> fInitialDistance = MyUnits<T>::angstrom();
        MyUnits<T> fCurDistance = fInitialDistance;
        MyUnits<T> fInitialPotentialEnergy = chargePotentialEnergy(MyUnits<T>::electronCharge(), MyUnits<T>::electronCharge(), fCurDistance);
        MyUnits<T> fCurSpeed;
        while (fCurDistance * 2 > fInitialDistance)
        {
            MyUnits<T> fCurForce = coloumbLaw(MyUnits<T>::electronCharge(), MyUnits<T>::electronCharge(), fCurDistance, fCurDistance * fCurDistance);
            MyUnits<T> fCurAcceleration = fCurForce / MyUnits<T>::electronMass();
            MyUnits<T> fWantedDeltaDistance = fCurDistance * 0.0001;
            MyUnits<T> fDeltaTime = sqrt(fWantedDeltaDistance * 2 / fCurAcceleration);
            MyUnits<T> fNextSpeed = fCurSpeed + fCurAcceleration * fDeltaTime;
            MyUnits<T> fAverageSpeed = (fNextSpeed + fCurSpeed) / 2;
            MyUnits<T> fDeltaDistance = fAverageSpeed * fDeltaTime;
            MyUnits<T> fNextDistance = fCurDistance - fDeltaDistance;
            nvAssert(fNextDistance < fCurDistance);
            MyUnits<T> fNextPotentialEnergy = chargePotentialEnergy(MyUnits<T>::electronCharge(), MyUnits<T>::electronCharge(), fNextDistance);
            MyUnits<T> fNextKineticEnergy = MyUnits<T>::electronMass() * fNextSpeed * fNextSpeed / 2;
            MyUnits<T> fNextTotalEnergy = fNextPotentialEnergy + fNextKineticEnergy;

            fDisparity = ((fNextTotalEnergy - fInitialPotentialEnergy) / (fNextTotalEnergy + fInitialPotentialEnergy)).m_value;
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
std::unordered_map<NvU32, BondsDataBase<double>::Atom> BondsDataBase<double>::m_atoms;
MyUnits<double> BondsDataBase<double>::s_zeroForceDist = MyUnits<double>::angstrom() * 4;
MyUnits<double> BondsDataBase<double>::s_zeroForceDistSqr = BondsDataBase<double>::s_zeroForceDist * BondsDataBase<double>::s_zeroForceDist;

template <class T>
void BondsDataBase<T>::init()
{
    setAtom(NPROTONS_O, MyUnits<T>::dalton() * 15.999, MyUnits<T>::picometer() * 152, 3.44);
    setAtom(NPROTONS_H, MyUnits<T>::dalton() * 1.008, MyUnits<T>::picometer() * 120, 2.2);

    setBond(NPROTONS_O, NPROTONS_O, 1, MyUnits<T>::angstrom() * 1.278, MyUnits<T>::kJperMole() * 140);
    setBond(NPROTONS_O, NPROTONS_O, 2, MyUnits<T>::angstrom() * 1.2, MyUnits<T>::kJperMole() * 498);
    setBond(NPROTONS_O, NPROTONS_H, 1, MyUnits<T>::angstrom() * 1.2, MyUnits<T>::kJperMole() * 464);
    setBond(NPROTONS_H, NPROTONS_H, 1, MyUnits<T>::angstrom() * 0.74, MyUnits<T>::kJperMole() * 436);
}
template <class T>
void BondsDataBase<T>::setBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons, MyUnits<T> fBondLength, MyUnits<T> fBondEnergy)
{
    auto& eBond = accessEBond(nProtons1, nProtons2, nElectrons);
    eBond.m_fLength = fBondLength;
    eBond.m_fEnergy = fBondEnergy;
    // we need to compute sigma and epsilon to match fBondLength and fBondEnergy
    eBond.m_fSigma = fBondLength * pow(2, -1. / 6);
    eBond.m_fEpsilon = fBondEnergy;
    accessEBond(nProtons2, nProtons1, nElectrons) = eBond;
}
template <class T>
void BondsDataBase<T>::setAtom(NvU32 nProtons, MyUnits<T> fMass, MyUnits<T> fRadius, T fElectroNegativity)
{
    auto& atom = m_atoms[nProtons];
    atom.m_fMass = fMass;
    atom.m_fRadius = fRadius;
    atom.m_fElectroNegativity = fElectroNegativity;
}
