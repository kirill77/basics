#pragma once
// those classes provide SparseArray wrappers for commonly used array types

#include <vector>

struct DenseStdArray
{
    template <class T>
    void init(const std::vector<T>& p)
    {
    }
    template <class T>
    bool isValidIndex(const std::vector<T>& p, NvU32 u) const { return true; }
    void setValidIndex(NvU32 u)
    {
        // DenseStdArray assumes that all indices are valid
        nvAssert(false);
    }
    template <class T>
    NvU32 getFirstValidIndex(const std::vector<T>& p)
    {
        return 0;
    }
    template <class T>
    NvU32 getNextValidIndex(const std::vector<T>& p, NvU32 u)
    {
        return ++u < p.size() ? u : INVALID_UINT32;
    }
    template <class T>
    NvU32 getNumValidIndices(const std::vector<T>& p) const
    {
        return (NvU32)p.size();
    }
};

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
    bool isValidIndex(NvU32 u) const { return m_validIndicesBitfield[u / 32] & (1 << (u & 31)); }
    bool isValidIndex(const BaseArray& p, NvU32 u) const { return isValidIndex(u); }
    void makeIndexValid(NvU32 u)
    {
        nvAssert(u < m_nDbgBaseSize);
        if (isValidIndex(u)) // that's to avoid doing push_back() twice on the same index
            return;
        m_validIndicesBitfield[u / 32] |= (1 << (u & 31));
        m_validIndices.push_back(u);
        nvAssert(m_validIndices.size() <= m_nDbgBaseSize);
    }
    NvU32 getFirstValidIndex(const BaseArray& p)
    {
        nvAssert(m_nDbgBaseSize == p.size());
        return m_validIndices.size() != 0 ? m_validIndices[m_curIndex = 0] : INVALID_UINT32;
    }
    NvU32 getNextValidIndex(const BaseArray& p, NvU32 u)
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