#pragma once

#include <array>
#include <vector>
#include <intrin.h>
#include <emmintrin.h>
#include "mybasics.h"
#include "bonds.h"

extern NvU32 g_debugCount;

// this struct is per atom. it stores 8 connections to other atoms
struct SSEForceBinding
{
    SSEForceBinding() :
        m_uOtherAtomsHi(_mm_set1_epi32(INVALID_UINT32)), m_uOtherAtomsLo(_mm_set1_epi32(INVALID_UINT32))
    {
        for (NvU32 u = 0; u < 8; ++u)
        {
            m_uForces[u] = INVALID_UINT32;
        }
    }
    // find a slot with the specified other atom in it (INVALID_UINT32 corresponds to empty slots)
    bool findBoundAtom(NvU32 uOtherAtom, NvU32& uSlot) const
    {
        NvU32 maskLo = _mm_movemask_epi8(_mm_cmpeq_epi16(m_uOtherAtomsLo, _mm_set1_epi16(uOtherAtom)));
        if (maskLo == 0) return false;
        NvU32 maskHi = _mm_movemask_epi8(_mm_cmpeq_epi16(m_uOtherAtomsHi, _mm_set1_epi16(uOtherAtom >> 16)));
        NvU32 equalMask = maskLo & maskHi;
        bool isFound = (equalMask != 0);
        if (isFound)
        {
            // skip first 4 slots if they're empty (perf optimization)
            uSlot = 0;
            if ((equalMask & 0xff) == 0)
            {
                uSlot += 4;
                equalMask >>= 8;
            }
            for ( ;(equalMask & 1) == 0; ++uSlot, equalMask >>= 2);
        }
#if 1 // can be used to check that the slot logic above is correct
        NvU32 dbgSlotIndex = 0;
        nvAssert(dbgFindSlot(uOtherAtom, dbgSlotIndex) == isFound && (!isFound || dbgSlotIndex == uSlot));
#endif
        return isFound;
    }
    bool isEmpty() const
    {
        NvU32 maskLo = _mm_movemask_epi8(_mm_cmpeq_epi16(m_uOtherAtomsLo, _mm_set1_epi16(-1)));
        if (maskLo != 0xffff) return false;
        NvU32 maskHi = _mm_movemask_epi8(_mm_cmpeq_epi16(m_uOtherAtomsHi, _mm_set1_epi16(-1)));
        if (maskHi != 0xffff) return false;
        return true;
    }
    NvU32 getForceIndex(NvU32 uSlot) const
    {
        return m_uForces[uSlot];
    }
    void setSlot(NvU32 uSlot, NvU32 uOtherAtom, NvU32 uForce)
    {
        m_uForces[uSlot] = uForce;
        m_uOtherAtomsHi.m128i_u16[uSlot] = uOtherAtom >> 16;
        m_uOtherAtomsLo.m128i_u16[uSlot] = uOtherAtom;
    }

private:
    bool dbgFindSlot(NvU32 uOtherAtom, NvU32& uOutSlot) const
    {
        for (uOutSlot = 0; uOutSlot < 8; ++uOutSlot)
        {
            NvU32 uTmp = (m_uOtherAtomsLo.m128i_u16[uOutSlot] & 0xffff) | (m_uOtherAtomsHi.m128i_u16[uOutSlot] << 16);
            if (uTmp == uOtherAtom) return true;
        }
        return false;
    }
    __m128i m_uOtherAtomsHi;
    __m128i m_uOtherAtomsLo;
    std::array<NvU32, 8> m_uForces;
};

struct ForceBindings
{
    void init(NvU32 nAtoms)
    {
        nvAssert(m_nAtoms == 0);
        m_nAtoms = nAtoms;
        m_atomLists.resize(m_nAtoms);
        m_uNext.resize(m_nAtoms);
        std::fill(m_uNext.begin(), m_uNext.end(), INVALID_UINT32);
    }
    NvU32 getExistingBindingSlot(NvU32 uAtom1, NvU32 uAtom2, NvU32& uOutSlot) const
    {
        NvU32 uPrev = 0;
        NvU32 u = findSlotForForce(uAtom1, uAtom2, uOutSlot, uPrev);
        nvAssert(u != INVALID_UINT32); // must exist
        return u;
    }
    NvU32 findOrCreateBindingSlot(NvU32 uAtom1, NvU32 uAtom2, NvU32& uOutSlot)
    {
        NvU32 u = 0, uPrev = 0;
        u = findSlotForForce(uAtom1, uAtom2, uOutSlot, uPrev);
        if (u == INVALID_UINT32)
        {
            u = findSlotForForce(uAtom1, INVALID_UINT32, uOutSlot, uPrev);
            if (u == INVALID_UINT32)
            {
                u = allocateNewBinding(uAtom1);
                uOutSlot = 0;
            }
        }
        return u;
    }
    void notifyForceDissociated(NvU32 uAtom1, NvU32 uAtom2)
    {
        NvU32 uSlot = INVALID_UINT32, uPrev = INVALID_UINT32;
        NvU32 u = findSlotForForce(uAtom1, uAtom2, uSlot, uPrev);
        nvAssert(u != INVALID_UINT32); // if they want to dissociated a force - it must be there

        SSEForceBinding& p = m_atomLists[u];
        p.setSlot(uSlot, INVALID_UINT32, INVALID_UINT32);
        if (p.isEmpty() && u != uAtom1) // if u == uAtom1 - can't unlink because that slot always belongs to uAtom1
        {
            nvAssert(uPrev != INVALID_UINT32);
            m_uNext[uPrev] = m_uNext[u]; // remove it from the list where it has been before...
            m_uNext[u] = m_firstFreeBinding; //... and add it to the free list
            m_firstFreeBinding = u;
        }
    }
    SSEForceBinding& operator[](NvU32 u) { return m_atomLists[u]; }
    const SSEForceBinding& operator[](NvU32 u) const { return m_atomLists[u]; }

private:
    NvU32 findSlotForForce(NvU32 uAtom1, NvU32 uAtom2, NvU32& uOutSlot, NvU32 &_uPrev) const
    {
        for (NvU32 u = uAtom1, uPrev = INVALID_UINT32; ; )
        {
            const SSEForceBinding& p1 = m_atomLists[u];
            if (p1.findBoundAtom(uAtom2, uOutSlot))
            {
                _uPrev = uPrev;
                return u;
            }
            uPrev = u;
            u = m_uNext[u];
            if (u == INVALID_UINT32) return INVALID_UINT32;
        }
    }
    NvU32 allocateNewBinding(NvU32 uAtom)
    {
        NvU32 u = 0;
        if (m_firstFreeBinding != INVALID_UINT32) // if free list isn't empty - get it from there...
        {
            NvU32 uTmp = m_uNext[m_firstFreeBinding];
            m_uNext[m_firstFreeBinding] = m_uNext[uAtom];
            m_uNext[uAtom] = u = m_firstFreeBinding;
            m_firstFreeBinding = uTmp;
        }
        else //... otherwise allocate a new one
        {
            u = (NvU32)m_atomLists.size();
            m_atomLists.resize(u + 1);
            // each binding is 8 forces. assume an atom can't have more than that many forces acting on it
            nvAssert(m_atomLists.size() < m_nAtoms * 128 / 8);
            m_uNext.resize(u + 1);
            m_uNext[u] = m_uNext[uAtom];
            m_uNext[uAtom] = u;
        }
        nvAssert(m_atomLists[u].isEmpty());
        return u;
    }
    NvU32 m_nAtoms = 0, m_firstFreeBinding = INVALID_UINT32;
    std::vector<SSEForceBinding> m_atomLists;
    std::vector<NvU32> m_uNext;
};

template <class T>
struct ForceMap
{
    void init(NvU32 nAtoms)
    {
        m_nAtoms = nAtoms;
        m_bindings.init(nAtoms);
    }
    bool isValid(NvU32 uForce) const { return m_forces[uForce].isValid(); }
    size_t size() const { return m_forces.size(); }
    NvU32 findExistingForceIndex(NvU32 uAtom1, NvU32 uAtom2) const
    {
        NvU32 uSlot = 0;
        NvU32 u = m_bindings.getExistingBindingSlot(uAtom1, uAtom2, uSlot);
        NvU32 uForce = m_bindings[u].getForceIndex(uSlot);
        nvAssert(uForce != INVALID_UINT32);
        return uForce;
    }
    NvU32 createForce(NvU32 uAtom1, NvU32 uAtom2)
    {
        // if such force already exists - just return its index
        NvU32 uSlot1 = 0, uSlot2 = 0;
        NvU32 u1 = m_bindings.findOrCreateBindingSlot(uAtom1, uAtom2, uSlot1);
        NvU32 uForce = m_bindings[u1].getForceIndex(uSlot1);
        if (uForce != INVALID_UINT32)
        {
#if ASSERT_ONLY_CODE
            NvU32 u2 = m_bindings.findOrCreateBindingSlot(uAtom2, uAtom1, uSlot2);
            nvAssert(uForce == m_bindings[u2].getForceIndex(uSlot2)); // must be symmetric
#endif
            return uForce;
        }
        uForce = allocateNewForce(uAtom1, uAtom2);
        NvU32 u2 = m_bindings.findOrCreateBindingSlot(uAtom2, uAtom1, uSlot2);
        nvAssert(m_bindings[u2].getForceIndex(uSlot2) == INVALID_UINT32); // must be symmetric
        m_bindings[u1].setSlot(uSlot1, uAtom2, uForce);
        m_bindings[u2].setSlot(uSlot2, uAtom1, uForce);
        return uForce;
    }
    NvU32 getFirstInvalidIndex() const
    {
        return m_firstUnusedForce;
    }
    NvU32 getNextInvalidIndex(NvU32 u) const
    {
        nvAssert(!isValid(u));
        return m_forces[u].getAtom2Index();
    }
    Force<T>& accessForceByIndex(NvU32 u) { nvAssert(m_forces[u].isValid()); return m_forces[u]; }
    const Force<T>& accessForceByIndex(NvU32 u) const { nvAssert(m_forces[u].isValid()); return m_forces[u]; }
    void notifyForceDissociated(NvU32 uDissociatedForce)
    {
        Force<T>& force = m_forces[uDissociatedForce];
        nvAssert(force.isValid());
        m_bindings.notifyForceDissociated(force.getAtom1Index(), force.getAtom2Index());
        m_bindings.notifyForceDissociated(force.getAtom2Index(), force.getAtom1Index());
        force = Force<T>(INVALID_UINT32, m_firstUnusedForce);
        m_firstUnusedForce = uDissociatedForce;
        nvAssert(!force.isValid());
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
            // assume an atom can't have more than that many forces acting on it
            nvAssert(m_forces.size() < m_nAtoms * 128);
        }
        m_forces[uForce] = Force<T>(uAtom1, uAtom2);
        nvAssert(m_forces[uForce].isValid());
        return uForce;
    }
    NvU32 m_nAtoms = 0, m_firstUnusedForce = INVALID_UINT32;
    std::vector<Force<T>> m_forces;
    ForceBindings m_bindings;
};
