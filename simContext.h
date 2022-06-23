#pragma once

#include "bonds.h"
#include "forceMap.h"
#include "boxWrapper.h"

template <class T>
struct SimContext
{
    std::vector<Atom<T>> m_atoms;
    ForceMap<T> m_forces;
    BoxWrapper<T> m_bBox;
};
