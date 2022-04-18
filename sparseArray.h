#pragma once
// those classes provide SparseArray wrappers for commonly used array types

#include <vector>

// wrapper for std::vector<>
template <class T>
struct DenseStdArray
{
    void init(const std::vector<T>& p)
    {
    }
    bool isIndexValid(NvU32 u) const { return true; }
    NvU32 findFirstValidIndex(const std::vector<T>& p)
    {
        return 0;
    }
    NvU32 findNextValidIndex(const std::vector<T>& p, NvU32 u)
    {
        return ++u < p.size() ? u : INVALID_UINT32;
    }
    NvU32 getNumValidIndices(const std::vector<T>& p) const
    {
        return (NvU32)p.size();
    }
};

// wrapper for the arrays that already have findFirst/NextValidIndex
template <class BaseArray>
struct DenseArray
{
    void init(const BaseArray& p)
    {
    }
    NvU32 findFirstValidIndex(const BaseArray& p)
    {
        return p.findFirstValidIndex();
    }
    NvU32 findNextValidIndex(const BaseArray& p, NvU32 u)
    {
        return p.findNextValidIndex(u);
    }
};

// wrapper for any kind of array having a size() method
template <class BaseArray>
struct SparseArray
{
    void init(const BaseArray& p)
    {
#if ASSERT_ONLY_CODE
        m_nDbgBaseSize = (NvU32)p.size();
#endif
        m_validIndicesBitfield.assign(((NvU32)p.size() + 31) / 32, 0);
        m_validIndices.resize(0);
    }
    bool isIndexValid(NvU32 u) const { return m_validIndicesBitfield[u / 32] & (1 << (u & 31)); }
    void makeIndexValid(NvU32 u)
    {
        nvAssert(u < m_nDbgBaseSize);
        if (isIndexValid(u)) // that's to avoid doing push_back() twice on the same index
            return;
        m_validIndicesBitfield[u / 32] |= (1 << (u & 31));
        m_validIndices.push_back(u);
        nvAssert(m_validIndices.size() <= m_nDbgBaseSize);
    }
    NvU32 findFirstValidIndex(const BaseArray& p)
    {
        nvAssert(m_nDbgBaseSize == p.size());
        return m_validIndices.size() != 0 ? m_validIndices[m_curIndex = 0] : INVALID_UINT32;
    }
    NvU32 findNextValidIndex(const BaseArray& p, NvU32 u)
    {
        nvAssert(m_nDbgBaseSize == p.size());
        return ++m_curIndex < m_validIndices.size() ? m_validIndices[m_curIndex] : INVALID_UINT32;
    }
    NvU32 getNumValidIndices(const BaseArray& p) const
    {
        nvAssert(m_nDbgBaseSize == p.size());
        return (NvU32)m_validIndices.size();
    }
private:
    NvU32 m_curIndex = 0;
#if ASSERT_ONLY_CODE
    NvU32 m_nDbgBaseSize = 0;
#endif
    std::vector<NvU32> m_validIndices;
    std::vector<NvU32> m_validIndicesBitfield;
};