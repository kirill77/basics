#pragma once

#include <vector>
#include <intrin.h>
#include <emmintrin.h>
#include "mybasics.h"
#include "bonds.h"

// this struct is per atom. it stores 8 connections to other atoms
struct SSEForceBinding
{
    SSEForceBinding() :
        m_uOtherAtomsHi(_mm_set1_epi32(INVALID_UINT32)), m_uOtherAtomsLo(_mm_set1_epi32(INVALID_UINT32)),
        m_uForcesHi(_mm_set1_epi32(INVALID_UINT32)), m_uForcesLo(_mm_set1_epi32(INVALID_UINT32))
    { }

    // find a slot with the specified other atom in it (INVALID_UINT32 corresponds to empty slots)
    bool findBoundAtom(NvU32 uOtherAtom, NvU32& uSlot)
    {
        NvU32 equalMask = _mm_movemask_epi8(_mm_and_si128(
            _mm_cmpeq_epi16(m_uOtherAtomsHi, _mm_set1_epi16(uOtherAtom >> 16)),
            _mm_cmpeq_epi16(m_uOtherAtomsLo, _mm_set1_epi16(uOtherAtom & 0xffff))));
        bool isFound = (equalMask != 0);
        if (isFound)
        {
            for (uSlot = 0; (equalMask & 1) == 0; ++uSlot, equalMask >>= 2);
        }
#if 1 // can be used to check that the slot logic above is correct
        NvU32 dbgSlotIndex = 0;
        nvAssert(dbgFindSlot(uOtherAtom, dbgSlotIndex) == isFound && (!isFound || dbgSlotIndex == uSlot));
#endif
        return isFound;
    }
    NvU32 getForceIndex(NvU32 uSlot) const
    {
        return (m_uForcesLo.m128i_u16[uSlot] & 0xffff) | (m_uForcesHi.m128i_u16[uSlot] << 16);
    }
    void setSlot(NvU32 uSlot, NvU32 uOtherAtom, NvU32 uForce)
    {
        nvAssert(uSlot < 8);
        m_uOtherAtomsHi.m128i_u16[uSlot] = uOtherAtom >> 16;
        m_uForcesHi.m128i_u16[uSlot] = uForce >> 16;
        m_uOtherAtomsLo.m128i_u16[uSlot] = uOtherAtom;
        m_uForcesLo.m128i_u16[uSlot] = uForce;
    }

private:
    bool dbgFindSlot(NvU32 uOtherAtom, NvU32& uOutSlot)
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
    __m128i m_uForcesHi;
    __m128i m_uForcesLo;
};
template <NvU32 N>
struct ForceBinding
{
    // find a slot with specified other atom in it (INVALID_UINT32 corresponds to an empty slots)
    SSEForceBinding *findBoundAtom(NvU32 uOtherAtom, NvU32& uOutSlot)
    {
        for (NvU32 u = 0; u < N; ++u)
        {
            if (m_bindings[u].findBoundAtom(uOtherAtom, uOutSlot))
            {
                return &m_bindings[u];
            }
        }
        return nullptr;
    }
    NvU32 getForceIndex(NvU32 uSlot) const { nvAssert(uSlot < 8 * N); return m_bindings[uSlot / 8].getForceIndex(uSlot & 7); }
    void setSlot(NvU32 uSlot, NvU32 uOtherAtom, NvU32 uForce)
    {
        nvAssert(uSlot < 8 * N);
        m_bindings[uSlot / 8].setSlot(uSlot & 7, uOtherAtom, uForce);
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
    bool isValid(NvU32 uForce) const { return m_forces[uForce].isValid(); }
    NvU32 size() const { return (NvU32)m_forces.size(); }
    NvU32 createForce(NvU32 uAtom1, NvU32 uAtom2)
    {
        // if such force already exists - just return its index
        NvU32 uSlot1 = 0;
        SSEForceBinding* p1 = m_atomForceIndices[uAtom1].findBoundAtom(uAtom2, uSlot1);
        if (p1)
        {
            NvU32 uSlot2 = 0;
            SSEForceBinding* p2 = m_atomForceIndices[uAtom2].findBoundAtom(uAtom1, uSlot2);
            nvAssert(p2 && p1->getForceIndex(uSlot1) == p2->getForceIndex(uSlot2)); // must be symmetric
            return p1->getForceIndex(uSlot1);
        }
        NvU32 uNewForce = allocateNewForce(uAtom1, uAtom2);
        bindAtomsInternal(uAtom1, uAtom2, uNewForce);
        bindAtomsInternal(uAtom2, uAtom1, uNewForce);
        return uNewForce;
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
        nvAssert(m_forces[uDissociatedForce].isValid());
        m_forces[uDissociatedForce] = Force<T>(INVALID_UINT32, m_firstUnusedForce);
        m_firstUnusedForce = uDissociatedForce;
        nvAssert(!m_forces[uDissociatedForce].isValid());
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
        SSEForceBinding* p = m_atomForceIndices[uAtom1].findBoundAtom(INVALID_UINT32, uSlot);
        if (!p)
        {
            __debugbreak(); // this means we don't have enough index slots
        }
        p->setSlot(uSlot, uAtom2, uForce);
        nvAssert(m_forces[uForce].isValid());
    }
    NvU32 m_firstUnusedForce = INVALID_UINT32;
    std::vector<Force<T>> m_forces;
    std::vector<ForceBinding<16>> m_atomForceIndices;
};
