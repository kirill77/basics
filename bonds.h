#pragma once

// The goal of those objects is to compute interaction between 2 atoms. we consider several kinds of forces:
// 1. Covalent bond. Covalent bond forms between two adjacent atoms when they both have free valence electrons.
//    1.1 It is modeled by lennard-jones potential;
//    1.2 It may store some of the atoms kinetic energy - so the atoms stick together instead of executing an elastic collision;
//    1.3 May change atom fractional charge - that way atoms may start experiencing Coloumb force with other atoms;
// 2. Coloumb force. Atoms are assigned this force if they are not covalently bonded and charge1*charge2 != 0.
// 3. Van-der-Waals force. Atoms are assigned this force if they are not covalently bonded and charge1*charge2 == 0.
// 4. Bond angle force. This force is acting on 3 covalently-bound atoms and is trying to support angles between them that are known
//    from experiments.

#include "myunits.h"
#include <unordered_map>

const NvU32 NPROTONS_MAX = 32;
const NvU32 NPROTONS_O = 8;
const NvU32 NPROTONS_H = 1;
const NvU32 MAX_ELECTRONS_PER_BOND = 3;
typedef NvU32 ATOM_KEY;

// bond length and energy, and corresponding lennard-jones parameters
template <class T>
struct LJ_Out
{
    T fForceTimesR, fPotential;
    T fNormalizedForceTimesR; // at extremal force point this equals to just R (because 1 * R)
};

template <class T>
struct EBond
{
    static const T s_zeroForceDist, s_zeroForceDistSqr;

    EBond() { }
    static constexpr NvU32 MAX_ALLOWED_ENERGY_KJ_PER_MOLE = 600; // O-O bond energy is 140 KJ per mole
    EBond(double fBondLength, double fBondEnergy)
    {
#if ASSERT_ONLY_CODE
        m_fDbgBondLength = (T)fBondLength;
#endif
        m_fLength = (T)fBondLength;
        m_fLengthSqr = (T)sqr(fBondLength);
        m_fDissocLengthSqr = (T)sqr(fBondLength * 2); // TODO: this is ad-hoc - figure out better way
        m_fAssocLengthSqr = (T)sqr(fBondLength * 2);
        // we need to compute sigma and epsilon to match fBondLength and fBondEnergy
        m_fSigma = (T)(fBondLength * pow(2, -1. / 6));
        m_fEpsilon = (T)fBondEnergy;
        // solved in wolfram:
        // fDistSqr = fDist*fDist
        // fSigma = fBondLength * 2 ^ (-1. / 6)
        // fPow2 = fSigma * fSigma / fDistSqr
        // fPow6 = fPow2 * fPow2 * fPow2
        // fPow12 = fPow6 * fPow6
        // fPotential = fEpsilon * 4 * (fPow12 - fPow6)
        // fForceTimesR = fEpsilon * 24 * (fPow12 * 2 - fPow6)
        // dForce=D[fForceTimesR/fDist, fDist]
        // Solve[dForce == 0, fDist]
        T fExtremalForceDist = (T)(fBondLength * 1.10868); // that value is from wolfram
        LJ_Out<T> out;
        m_fExtremalForce = 1; // to avoid division by zero inside lennardJones
        bool bHasForce = lennardJones(fExtremalForceDist * fExtremalForceDist, out);
        nvAssert(bHasForce);
        m_fExtremalForce = out.fForceTimesR / fExtremalForceDist;
    }
    bool isValid() const { return m_fEpsilon > 0; }
    MyUnits<T> getEnergy() const { return m_fEpsilon; }
    MyUnits<T> getLength() const { return m_fLength; }
    T getAssocLengthSqr() const { return m_fAssocLengthSqr; }
    T getDissocLengthSqr() const { return m_fDissocLengthSqr; }
    bool lennardJones(T fDistSqr, LJ_Out<T>& out, bool isCovalentBond = true) const
    {
        nvAssert(m_fEpsilon > 0 && m_fSigma > 0);
        if (fDistSqr >= s_zeroForceDistSqr)
            return false;
        MyUnits<T> fPow2 = m_fSigma * m_fSigma / fDistSqr;
        MyUnits<T> fPow6 = fPow2 * fPow2 * fPow2;
        MyUnits<T> fPow12 = fPow6 * fPow6;
        if (!isCovalentBond) fPow6 = MyUnits<T>(); // remove attractive term for non-covalent bonds
        out.fPotential = m_fEpsilon * 4 * (fPow12 - fPow6);
        out.fForceTimesR = m_fEpsilon * 24 * (fPow12 * 2 - fPow6);
        out.fNormalizedForceTimesR = (out.fForceTimesR / m_fExtremalForce);
        nvAssert(!isnan(out.fForceTimesR));
        return true;
    }
private:
    T m_fLengthSqr = 0, m_fSigma = 0, m_fEpsilon = 0, m_fLength = 0, m_fAssocLengthSqr = 0, m_fDissocLengthSqr = 0, m_fExtremalForce = 0;
#if ASSERT_ONLY_CODE
    T m_fDbgBondLength;
#endif
};

template <class T>
struct BondsDataBase
{
    // describes different bonds that may happen between particular types of atoms (for instance O-O and O=O would be in the same ABond, but O-H would be in different ABond)
    struct ABond
    {
        const EBond<T>& operator[](NvU32 nElectrons) const { nvAssert(nElectrons < MAX_ELECTRONS_PER_BOND); return m_eBonds[nElectrons]; }
        EBond<T>& operator[](NvU32 nElectrons) { nvAssert(nElectrons < MAX_ELECTRONS_PER_BOND); return m_eBonds[nElectrons]; }
    private:
        // there can be up to 3 electrons creating the bond. more electrons = shorter length and larger energy
        EBond<T> m_eBonds[MAX_ELECTRONS_PER_BOND];
    };
    static inline const EBond<T>& getEBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons) { return accessEBond(nProtons1, nProtons2, nElectrons); }
    static void init();
    static const std::unordered_map<NvU32, ABond>& getABonds() { return m_aBonds; }

    struct Element
    {
        double m_fMass = 0, m_fRadius = 0;
        T m_fElectroNegativity = 0;
        NvU32 m_uValence = -1;
    };
    static const Element& getElement(NvU32 nProtons)
    {
        const auto& element = m_elements[nProtons];
        nvAssert(element.m_fMass > 0);
        return element;
    }

private:
    static inline ATOM_KEY computeAtomDatasKey(NvU32 nProtons1, NvU32 nProtons2)
    {
        nvAssert(nProtons1 <= NPROTONS_MAX && nProtons2 <= NPROTONS_MAX);
        return nProtons1 * NPROTONS_MAX + nProtons2;
    }
    static inline EBond<T>& accessEBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons)
    {
        NvU32 key = computeAtomDatasKey(nProtons1, nProtons2);
        return m_aBonds[key][nElectrons];
    }
    static void setBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons, double fBondLength, double fBondEnergy);
    static void setAtom(NvU32 nProtons, double fMass, double fRadius, double fElectroNegativity, NvU32 uValence);
    static std::unordered_map<ATOM_KEY, ABond> m_aBonds;
    static std::unordered_map<NvU32, Element> m_elements;
};

template <class T>
struct Atom
{
    explicit Atom(NvU32 nProtons = 1) : m_nProtons(nProtons), m_uValence(BondsDataBase<T>::getElement(m_nProtons).m_uValence), m_nCovBonds(0)
    {
        nvAssert(m_uValence != 0); // we don't work with noble gasses
    }

    NvU32 getNProtons() const { return m_nProtons; }
    T getMass() const { return (T)BondsDataBase<T>::getElement(m_nProtons).m_fMass; }
    NvU32 getValence() const { return m_uValence; }
    NvU32 getNCovBonds() const { return m_nCovBonds; }
    void setNCovBonds(NvU32 nCovBonds)
    {
        m_nCovBonds = nCovBonds;
        nvAssert(m_nCovBonds == nCovBonds);
    }

    rtvector<MyUnits<T>, 3> m_vPos, m_vSpeed;

private:
    union
    {
        NvU32 flags;
        struct
        {
            NvU32 m_nProtons : 8;
            NvU32 m_uValence : 3;
            NvU32 m_nCovBonds : 3;
        };
    };
};

template <class T>
struct ForceData
{
    rtvector<MyUnits<T>, 3> vForce;
    T fNormalizedForce;
};
template <class T>
struct Force
{
    Force() : m_isCovalentBond(0) { }

    Force(NvU32 uAtom1, NvU32 uAtom2) : m_uAtom1(uAtom1), m_uAtom2(uAtom2), m_isCovalentBond(0), m_prevCovalentState(0)
    {
        nvAssert(m_uAtom1 != m_uAtom2 || !isValid());
    }
    bool isValid() const { return m_uAtom1 != INVALID_UINT32; }
    NvU32 getAtom1Index() const { return m_uAtom1; }
    NvU32 getAtom2Index() const { return m_uAtom2; }
#if ASSERT_ONLY_CODE
    template <class T>
    bool dbgAreIndicesSane(const Atom<T>& atom1, const Atom<T>& atom2)
    {
        return &atom2 - &atom1 == (int)m_uAtom2 - (int)m_uAtom1;
    }
#endif

    // returns true if there is a force, false if the force is 0
    template <class WRAPPER>
    bool computeForce(const Atom<T>& atom1, const Atom<T>& atom2, const WRAPPER &w, ForceData<T> &outForceData) const
    {
        rtvector<MyUnits<T>, 3> vDir = w.computeDir(atom1, atom2);
        MyUnits<T> fDistSqr = dot(vDir, vDir);
        LJ_Out<T> out;
        const EBond<T>& eBond = BondsDataBase<T>::getEBond(atom1.getNProtons(), atom2.getNProtons(), 1);
        bool hasForce = eBond.lennardJones(fDistSqr, out, m_isCovalentBond);
        if (hasForce)
        {
            outForceData.vForce = vDir * (out.fForceTimesR / fDistSqr);
            outForceData.fNormalizedForce = out.fNormalizedForceTimesR / sqrt(fDistSqr);
            nvAssert(outForceData.fNormalizedForce < 1.01);
        }
        return hasForce;
    }
    bool createCovalentBondIfNeeded(Atom<T>& atom1, Atom<T>& atom2, T fDistSqr)
    {
        if (isCovalentBond() || atom1.getNCovBonds() >= atom1.getValence() || atom2.getNCovBonds() >= atom2.getValence())
            return false;
        const EBond<T>& eBond = BondsDataBase<T>::getEBond(atom1.getNProtons(), atom2.getNProtons(), 1);
        if (fDistSqr > eBond.getAssocLengthSqr())
            return false;
        setCovalentBond(atom1, atom2);
        return true;
    }
    // returns true if force is now zero, returns false otherwise
    template <class WRAPPER>
    bool dissociateWeakBond(Atom<T> &atom1, Atom<T> &atom2, const WRAPPER &w)
    {
        nvAssert(dbgAreIndicesSane(atom1, atom2));
        rtvector<T, 3> vDir = w.computeDir(atom1, atom2);
        auto fDistSqr = dot(vDir, vDir);

        const EBond<T>& bond = BondsDataBase<T>::getEBond(atom1.getNProtons(), atom2.getNProtons(), 1);
        // if atoms are too far apart - erase the force
        if (fDistSqr >= EBond<T>::s_zeroForceDistSqr)
        {
            if (isCovalentBond())
            {
                dropCovalentBond(atom1, atom2);
            }
            return true;
        }
        if (isCovalentBond() && fDistSqr >= bond.getDissocLengthSqr())
        {
            dropCovalentBond(atom1, atom2);
        }

        return false;
    }

    bool shouldDraw() const { return m_isCovalentBond; } // should we draw this?
    bool isCovalentBond() const { return m_isCovalentBond; }
    void setPrevCovalentState(bool value) { m_prevCovalentState = value; }
    bool getPrevCovalentState() const { return m_prevCovalentState; }

private:
    void setCovalentBond(Atom<T> &atom1, Atom<T> &atom2)
    {
        nvAssert(!isCovalentBond());
        m_isCovalentBond = 1;
        atom1.setNCovBonds(atom1.getNCovBonds() + 1);
        atom2.setNCovBonds(atom2.getNCovBonds() + 1);
    }
    void dropCovalentBond(Atom<T> &atom1, Atom<T> &atom2)
    {
        nvAssert(isCovalentBond());
        m_isCovalentBond = 0;
        atom1.setNCovBonds(atom1.getNCovBonds() - 1);
        atom2.setNCovBonds(atom2.getNCovBonds() - 1);
    }
    NvU32 m_isCovalentBond : 1;
    NvU32 m_prevCovalentState : 1;
    NvU32 m_uAtom1 = INVALID_UINT32, m_uAtom2 = INVALID_UINT32;
};
