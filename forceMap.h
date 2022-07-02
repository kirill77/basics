#pragma once

#include <vector>
#include "bonds.h"

struct ForceIndex
{
    NvU32 m_index = INVALID_UINT32;
};

template <class T>
struct ForceMap
{
    void init(NvU32 nAtoms)
    {
        m_nAtoms = nAtoms;
    }
    bool isValid(NvU32 uForce) const { return m_forces[uForce].isValid(); }
    size_t size() const { return m_forces.size(); }
    NvU32 findExistingForce(NvU32 uAtom1, NvU32 uAtom2) const
    {
        auto it = m_forceMap.find(makeKey(uAtom1, uAtom2));
        nvAssert(isValid(it->second.m_index));
        return it->second.m_index;
    }
    NvU32 findExistingOrCreateNewForce(NvU32 uAtom1, NvU32 uAtom2)
    {
        ForceIndex &i = m_forceMap[makeKey(uAtom1, uAtom2)];
        if (i.m_index == INVALID_UINT32)
        {
            i.m_index = allocateNewForce(uAtom1, uAtom2);
        }
        nvAssert(isValid(i.m_index));
        return i.m_index;
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
        m_forceMap.erase(makeKey(force.getAtom1Index(), force.getAtom2Index()));
        force = Force<T>(INVALID_UINT32, m_firstUnusedForce);
        m_firstUnusedForce = uDissociatedForce;
        nvAssert(!force.isValid());
    }

private:
    static NvU64 makeKey(NvU32 uAtom1, NvU32 uAtom2)
    {
        NvU64 uKey;
        if (uAtom1 < uAtom2)
        {
            ((NvU32*)&uKey)[0] = uAtom1;
            ((NvU32*)&uKey)[1] = uAtom2;
        }
        else
        {
            ((NvU32*)&uKey)[1] = uAtom1;
            ((NvU32*)&uKey)[0] = uAtom2;
        }
        return uKey;
    }
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
            // this is for sanity - we check that our force array doesn't grow uncontrollably due to some bug
            nvAssert(m_forces.size() < m_nAtoms * 128);
            nvAssert(m_forceMap.size() < m_nAtoms * 128);
        }
        m_forces[uForce] = Force<T>(uAtom1, uAtom2);
        nvAssert(m_forces[uForce].isValid());
        return uForce;
    }
    NvU32 m_nAtoms = 0, m_firstUnusedForce = INVALID_UINT32;
    std::vector<Force<T>> m_forces;
    std::unordered_map<NvU64, ForceIndex> m_forceMap;
};
