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

template <class T>
struct BondsDataBase
{
    static MyUnits<T> s_zeroForceDist, s_zeroForceDistSqr;

    // bond length and energy, and corresponding lennard-jones parameters
    struct LJ_Out
    {
        rtvector<MyUnits<T>, 3> vForce;
        MyUnits<T> fDistSqr, fPotential;
    };
    struct EBond
    {
        bool lennardJones(const rtvector<MyUnits<T>, 3>& vSrcToDstDir, LJ_Out &out) const
        {
            nvAssert(m_fEpsilon > 0 && m_fSigma > 0);
            out.fDistSqr = lengthSquared(vSrcToDstDir);
            if (out.fDistSqr >= s_zeroForceDistSqr)
                return false;
            MyUnits<T> fPow2 = m_fSigma * m_fSigma / out.fDistSqr;
            MyUnits<T> fPow6 = fPow2 * fPow2 * fPow2;
            MyUnits<T> fPow12 = fPow6 * fPow6;
            out.fPotential = m_fEpsilon * 4 * (fPow12 - fPow6);
            auto fForceTimesR = m_fEpsilon * 24 * (fPow12 * 2 - fPow6);
            nvAssert(!isnan(fForceTimesR.m_value));
            out.vForce = vSrcToDstDir * (fForceTimesR / out.fDistSqr);
            return true;
        }
        MyUnits<T> m_fLength, m_fLengthSqr, m_fDissocLengthSqr, m_fEnergy, m_fEpsilon, m_fSigma;
    };
    // describes different bonds that may happen between particular types of atoms (for instance O-O and O=O would be in the same ABond, but O-H would be in different ABond)
    struct ABond
    {
        const EBond& operator[](NvU32 nElectrons) const { nvAssert(nElectrons < MAX_ELECTRONS_PER_BOND); return m_eBonds[nElectrons]; }
        EBond& operator[](NvU32 nElectrons) { nvAssert(nElectrons < MAX_ELECTRONS_PER_BOND); return m_eBonds[nElectrons]; }
    private:
        // there can be up to 3 electrons creating the bond. more electrons = shorter length and larger energy
        EBond m_eBonds[MAX_ELECTRONS_PER_BOND];
    };
    static inline const EBond& getEBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons) { return accessEBond(nProtons1, nProtons2, nElectrons); }
    static void init();
    static const std::unordered_map<NvU32, ABond>& getABonds() { return m_aBonds; }

    struct Element
    {
        MyUnits<double> m_fMass, m_fRadius;
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
    static inline EBond& accessEBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons)
    {
        NvU32 key = computeAtomDatasKey(nProtons1, nProtons2);
        return m_aBonds[key][nElectrons];
    }
    static void setBond(NvU32 nProtons1, NvU32 nProtons2, NvU32 nElectrons, MyUnits<T> fBondLength, MyUnits<T> fBondEnergy);
    static void setAtom(NvU32 nProtons, MyUnits<T> fMass, MyUnits<T> fRadius, T fElectroNegativity, NvU32 uValence);
    static std::unordered_map<ATOM_KEY, ABond> m_aBonds;
    static std::unordered_map<NvU32, Element> m_elements;
};

template <class T>
struct Atom
{
    Atom(NvU32 nProtons = 1) : m_nBondedAtoms(0), m_nProtons(nProtons), m_uValence(BondsDataBase<T>::getElement(m_nProtons).m_uValence)
    {
        for (NvU32 u = 0; u < m_bondedAtoms.size(); ++u) m_bondedAtoms[u] = -1;
        nvAssert(m_uValence != 0); // we don't work with noble gasses
    }

    NvU32 getNProtons() const { return m_nProtons; }
    MyUnits<T> getMass() const { return BondsDataBase<T>::getElement(m_nProtons).m_fMass; }
    NvU32 getValence() const { return m_uValence; }
    NvU32 getNBonds() const { return m_nBondedAtoms; }
    NvU32 getBond(NvU32 uBond) const { nvAssert(uBond < m_nBondedAtoms); return m_bondedAtoms[uBond]; }

    void addBond(NvU32 uAtom)
    {
        nvAssert(m_nBondedAtoms < m_uValence);
        m_bondedAtoms[m_nBondedAtoms] = uAtom;
        nvAssert(m_bondedAtoms[m_nBondedAtoms] == uAtom); // check that type conversion didn't loose information
        ++m_nBondedAtoms;
    }
    void removeBond(NvU32 uAtom)
    {
        for (NvU32 u = 0; ; ++u)
        {
            nvAssert(u < m_nBondedAtoms);
            if (m_bondedAtoms[u] == uAtom)
            {
                m_bondedAtoms[u] = m_bondedAtoms[--m_nBondedAtoms];
                return;
            }
        }
    }

    rtvector<MyUnits<T>, 3> m_vPos[2], m_vSpeed[2], m_vForce;

private:
    union
    {
        NvU32 flags;
        struct
        {
            NvU32 m_nProtons : 8;
            NvU32 m_nBondedAtoms : 3;
            NvU32 m_uValence : 3;
        };
    };
    std::array<unsigned short, 4> m_bondedAtoms;
};

// force acting between two atoms (when we have force acting between 3 atoms we'll call it Force3)
struct ForceKey
{
    ForceKey(NvU32 uAtom1, NvU32 uAtom2)
    {
        if (uAtom1 < uAtom2) // sort indices to avoid duplicate forces 1<->2 and 2<->1
        {
            m_uAtom1 = uAtom1;
            m_uAtom2 = uAtom2;
            nvAssert(m_uAtom1 == uAtom1 && m_uAtom2 == uAtom2 && m_uAtom1 != m_uAtom2);
            return;
        }
        m_uAtom1 = uAtom2;
        m_uAtom2 = uAtom1;
        nvAssert(m_uAtom1 == uAtom2 && m_uAtom2 == uAtom1 && m_uAtom1 != m_uAtom2);
    }
    bool operator ==(const ForceKey& other) const { return m_uAtom1 == other.m_uAtom1 && m_uAtom2 == other.m_uAtom2; }
    NvU32 getAtom1Index() const { return m_uAtom1; }
    NvU32 getAtom2Index() const { return m_uAtom2; }
private:
    NvU32 m_uAtom1 : 16;
    NvU32 m_uAtom2 : 16;
};
// custom specialization of std::hash can be injected in namespace std
template<>
struct std::hash<ForceKey>
{
    std::size_t operator()(ForceKey const& s) const noexcept
    {
        static_assert(sizeof(s) == sizeof(NvU32), "error: wrong ForceKey size");
        return std::hash<NvU32>{}((NvU32&)s);
    }
};
template <class T>
struct Force
{
    Force() : m_collisionDetected(0), m_isCovalentBond(0) { }

    bool hadCollision() const { return m_collisionDetected; }
    bool isCovalentBond() const { return m_isCovalentBond; }

    void notifyCollision() { m_collisionDetected = 1; }
    void dropCovalentBond() { nvAssert(m_isCovalentBond == 1); m_isCovalentBond = 0; }
    void setCovalentBond() { nvAssert(m_isCovalentBond == 0); m_isCovalentBond = 1; }

    MyUnits<T> m_fPotential[2]; // potentials corresponding prev and next state of the system
    MyUnits<T> m_fDistSqr[2]; // distances between atoms corresponding to prev and next state of the system

    template <NvU32 index>
    bool computeForce(NvU32 nProtons1, NvU32 nProtons2, rtvector<MyUnits<T>, 3> vInR, rtvector<MyUnits<T>, 3>& vOutForce)
    {
        typename BondsDataBase<T>::LJ_Out out;
        auto& eBond = BondsDataBase<T>::getEBond(nProtons1, nProtons2, 1);
        bool hasForce = eBond.lennardJones(vInR, out);
        m_fDistSqr[index] = out.fDistSqr; // this is needed even if force is 0
        if (hasForce)
        {
            vOutForce = out.vForce;
            m_fPotential[index] = out.fPotential;
        }
        return hasForce;
    }

private:
    NvU32 m_collisionDetected : 1; // collision detected during time step
    NvU32 m_isCovalentBond : 1;
};
