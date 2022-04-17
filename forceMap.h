#pragma once

#include <vector>
#include <intrin.h>
#include <emmintrin.h>
#include "mybasics.h"
#include "bonds.h"

// this struct is per atom. it stores 4 connections to other atoms
struct SSEForceBinding
{
    SSEForceBinding() : m_uOtherAtoms(_mm_set1_epi32(INVALID_UINT32)), m_uForces(_mm_set1_epi32(INVALID_UINT32)) { }

    // find a slot with specified other atom in it (INVALID_UINT32 corresponds to an empty slots)
    bool findBoundAtom(NvU32 uOtherAtom, NvU32& uSlot)
    {
        NvU32 equalMask = _mm_movemask_epi8(_mm_cmpeq_epi32(m_uOtherAtoms, _mm_set1_epi32(uOtherAtom)));
        bool isFound = (equalMask != 0);
        if (isFound)
        {
            if (equalMask & 0x11)
            {
                uSlot = (equalMask & 0x001) ? 0 : 1;
            }
            else
            {
                uSlot = (equalMask & 0x100) ? 2 : 3;
            }
        }
#if 1 // can be used to check that the slot logic above is correct
        NvU32 dbgSlotIndex = 0;
        nvAssert(dbgFindSlot(uOtherAtom, dbgSlotIndex) == isFound && (!isFound || dbgSlotIndex == uSlot));
#endif
        return isFound;
    }
    NvU32 getForceIndex(NvU32 uSlot) const { nvAssert(uSlot < 4); return m_uForces.m128i_u32[uSlot]; }
    void setSlot(NvU32 uSlot, NvU32 uOtherAtom, NvU32 uForce)
    {
        nvAssert(uSlot < 4);
        m_uOtherAtoms.m128i_u32[uSlot] = uOtherAtom;
        m_uForces.m128i_u32[uSlot] = uForce;
    }

private:
    bool dbgFindSlot(NvU32 uOtherAtom, NvU32& uOutSlot)
    {
        for (uOutSlot = 0; uOutSlot < 4; ++uOutSlot)
        {
            if (m_uOtherAtoms.m128i_u32[uOutSlot] == uOtherAtom) return true;
        }
        return false;
    }
    __m128i m_uOtherAtoms; // connected atoms (INVALID_UINT32 means no connection)
    __m128i m_uForces; // force indices corresponding to those connected atoms
};
template <NvU32 N>
struct ForceBinding
{
    // find a slot with specified other atom in it (INVALID_UINT32 corresponds to an empty slots)
    bool findBoundAtom(NvU32 uOtherAtom, NvU32& uOutSlot)
    {
        for (NvU32 u = 0; u < N; ++u)
        {
            if (m_bindings[u].findBoundAtom(uOtherAtom, uOutSlot))
            {
                uOutSlot += u * 4;
                return true;
            }
        }
        return false;
    }
    NvU32 getForceIndex(NvU32 uSlot) const { nvAssert(uSlot < 4 * N); return m_bindings[uSlot / 4].getForceIndex(uSlot & 3); }
    void setSlot(NvU32 uSlot, NvU32 uOtherAtom, NvU32 uForce)
    {
        nvAssert(uSlot < 4 * N);
        m_bindings[uSlot / 4].setSlot(uSlot & 3, uOtherAtom, uForce);
    }

private:
    SSEForceBinding m_bindings[N];
};

template <class T>
struct ForceMap
{
    void init(NvU32 nAtoms)
    {
        m_atomForceIndices.resize(nAtoms);
    }
    NvU32 getNForces() const { return (NvU32)m_forces.size(); }
    NvU32 createForce(NvU32 uAtom1, NvU32 uAtom2)
    {
        // if such force already exists - just return its index
        NvU32 uSlot1 = 0;
        if (m_atomForceIndices[uAtom1].findBoundAtom(uAtom2, uSlot1))
        {
            NvU32 uSlot2 = 0;
            nvAssert(m_atomForceIndices[uAtom2].findBoundAtom(uAtom1, uSlot2) &&
                     m_atomForceIndices[uAtom1].getForceIndex(uSlot1) == m_atomForceIndices[uAtom2].getForceIndex(uSlot2)); // must be symmetric
            return m_atomForceIndices[uAtom1].getForceIndex(uSlot1);
        }
        NvU32 uNewForce = allocateNewForce(uAtom1, uAtom2);
        bindAtomsInternal(uAtom1, uAtom2, uNewForce);
        bindAtomsInternal(uAtom2, uAtom1, uNewForce);
        return uNewForce;
    }
    NvU32 findFirstForceIndex() const
    {
        for (NvU32 u = 0; u < m_forces.size(); ++u)
        {
            if (m_forces[u].isValid()) return u;
        }
        return INVALID_UINT32;
    }
    NvU32 findNextForceIndex(NvU32 u) const
    {
        for (++u; u < m_forces.size(); ++u)
        {
            if (m_forces[u].isValid()) return u;
        }
        return INVALID_UINT32;
    }
    Force<T>& accessForceByIndex(NvU32 u) { nvAssert(m_forces[u].isValid()); return m_forces[u]; }
    const Force<T>& accessForceByIndex(NvU32 u) const { nvAssert(m_forces[u].isValid()); return m_forces[u]; }
    NvU32 dissociateForce(NvU32 uDissociatedForce)
    {
        nvAssert(m_forces[uDissociatedForce].isValid());
        m_forces[uDissociatedForce] = Force<T>(INVALID_UINT32, m_firstUnusedForce);
        m_firstUnusedForce = uDissociatedForce;
        nvAssert(!m_forces[uDissociatedForce].isValid());
        return findNextForceIndex(uDissociatedForce);
    }

    struct ConstIt
    {
        ConstIt(NvU32 u, const ForceMap<T>& forces) : m_u(u), m_forces(forces) { }
        bool operator != (const ConstIt& other) { return m_u != other.m_u; }
        void operator ++() { m_u = m_forces.findNextForceIndex(m_u); }
        const Force<T> &operator *() const { return m_forces.accessForceByIndex(m_u); }
    private:
        NvU32 m_u;
        const ForceMap<T>& m_forces;
    };
    ConstIt begin() const
    {
        return ConstIt(findFirstForceIndex(), *this);
    }
    ConstIt end() const
    {
        return ConstIt(INVALID_UINT32, *this);
    }

private:
    NvU32 allocateNewForce(NvU32 uAtom1, NvU32 uAtom2)
    {
        NvU32 uForce = 0;
        if (m_firstUnusedForce != INVALID_UINT32)
        {
            uForce = m_firstUnusedForce;
            m_firstUnusedForce = m_forces[m_firstUnusedForce].getAtom2Index();
        }
        else
        {
            uForce = (NvU32)m_forces.size();
            m_forces.resize(m_forces.size() + 1);
        }
        m_forces[uForce] = Force<T>(uAtom1, uAtom2);
        nvAssert(m_forces[uForce].isValid());
        return uForce;
    }
    void bindAtomsInternal(NvU32 uAtom1, NvU32 uAtom2, NvU32 uForce)
    {
        NvU32 uSlot = 0;
        if (!m_atomForceIndices[uAtom1].findBoundAtom(INVALID_UINT32, uSlot))
        {
            __debugbreak(); // this means we don't have enough index slots
        }
        m_atomForceIndices[uAtom1].setSlot(uSlot, uAtom2, uForce);
        nvAssert(m_forces[uForce].isValid());
    }
    NvU32 m_firstUnusedForce = INVALID_UINT32;
    std::vector<Force<T>> m_forces;
    std::vector<ForceBinding<16>> m_atomForceIndices;
};
