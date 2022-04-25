#pragma once
// those classes provide SparseArray wrappers for commonly used array types
#include <array>
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

// SparseHierarchy is used for multi-layer hierarchical simulation. The very first layer contains all elements.
// As the simulation progresses - some elements need to be simulated in more detail than the others. This class
// allows efficient movement of elements between layers.
struct SparseHierarchy
{
    void init(NvU32 nElements)
    {
        nvAssert(nElements > 0);
        // initialize layers
        m_nLayers = 1;
        // initialize elements
        if (m_pElements.size() != nElements)
        {
            m_pLayers[0].m_uFirstEl = 0;
            m_pLayers[0].m_uLastEl = nElements - 1;
            m_pElements.resize(nElements);
            m_pElements.begin()->m_uPrev = INVALID_UINT32;
            m_pElements.rbegin()->m_uNext = INVALID_UINT32;
            for (NvU32 u = 1; u < nElements; ++u)
            {
                m_pElements[u - 1].m_uNext = u;
                m_pElements[u].m_uPrev = u - 1;
            }
        }
        nvAssert(dbgIsLayerSane(0));
    }
    bool hasElements(NvU32 uLayer) const
    {
        nvAssert(uLayer < m_nLayers);
        return m_pLayers[uLayer].m_uFirstEl != INVALID_UINT32;
    }
    bool hasElement(NvU32 uLayer, NvU32 uElement) const
    {
        nvAssert(uLayer < m_nLayers);
        return m_pElements[uElement].m_maxLayerId >= m_pLayers[uLayer].m_id;
    }
    void moveToMostDetailedLayer(NvU32 uMostDetailedLayer, NvU32 uElement)
    {
        nvAssert(uMostDetailedLayer == m_nLayers - 1);
        if (hasElement(uMostDetailedLayer, uElement)) return;
        m_pElements[uElement].m_maxLayerId = m_pLayers[uMostDetailedLayer].m_id;
        NvU32 uPrevLayer = uMostDetailedLayer - 1;
        nvAssert(dbgIsLayerSane(uMostDetailedLayer));
        nvAssert(dbgIsLayerSane(uPrevLayer));
        NvU32 uPrevEl = m_pElements[uElement].m_uPrev;
        NvU32 uNextEl = m_pElements[uElement].m_uNext;
        // unlinking
        if (uPrevEl == INVALID_UINT32)
        {
            m_pLayers[uPrevLayer].m_uFirstEl = uNextEl;
        }
        else
        {
            m_pElements[uPrevEl].m_uNext = uNextEl;
            m_pElements[uElement].m_uPrev = INVALID_UINT32;
        }
        if (uNextEl == INVALID_UINT32)
        {
            m_pLayers[uPrevLayer].m_uLastEl = uPrevEl;
        }
        else
        {
            m_pElements[uNextEl].m_uPrev = uPrevEl;
        }
        // linking
        NvU32 uFirstEl = m_pLayers[uMostDetailedLayer].m_uFirstEl;
        m_pLayers[uMostDetailedLayer].m_uFirstEl = uElement;
        m_pElements[uElement].m_uNext = uFirstEl;
        if (uFirstEl != INVALID_UINT32)
        {
            m_pElements[uFirstEl].m_uPrev = uElement;
        }
        if (m_pLayers[uMostDetailedLayer].m_uLastEl == INVALID_UINT32)
        {
            m_pLayers[uMostDetailedLayer].m_uLastEl = uElement;
        }
        nvAssert(dbgIsLayerSane(uMostDetailedLayer));
        nvAssert(dbgIsLayerSane(uPrevLayer));
    }
    void notifyLayerCreated(NvU32 uLayer)
    {
        nvAssert(uLayer == m_nLayers && m_nLayers < m_pLayers.size()); // each time must increase layer index by one
        m_nLayers = uLayer + 1;
        m_pLayers[uLayer].m_id = ++m_lastLayerId;
        m_pLayers[uLayer].m_uFirstEl = INVALID_UINT32;
        m_pLayers[uLayer].m_uLastEl = INVALID_UINT32;
    }
    void notifyLayerDestroyed(NvU32 uLayer)
    {
        nvAssert(uLayer == m_nLayers - 1); // each time must decrease layer index by one
        if (!hasElements(uLayer))
        {
            // if an empty layer is being destroyed - we can re-use its id in future
            --m_lastLayerId;
        }
        else
        {
            // since this layer is destroyed - all its elements are moved to the previous layer
            NvU32 uPrevLayer = uLayer - 1;
            if (hasElements(uPrevLayer))
            {
                // here we're linking two linked lists together (constant-time)
                m_pElements[m_pLayers[uLayer].m_uFirstEl].m_uPrev = m_pLayers[uPrevLayer].m_uLastEl;
                m_pElements[m_pLayers[uPrevLayer].m_uLastEl].m_uNext = m_pLayers[uLayer].m_uFirstEl;
                m_pLayers[uPrevLayer].m_uLastEl = m_pLayers[uLayer].m_uLastEl;
            }
            else
            {
                m_pLayers[uPrevLayer].m_uFirstEl = m_pLayers[uLayer].m_uFirstEl;
                m_pLayers[uPrevLayer].m_uLastEl = m_pLayers[uLayer].m_uLastEl;
            }
        }
        m_nLayers = uLayer;
        nvAssert(dbgIsLayerSane(uLayer - 1));
    }
    NvU32 getFirstLayerElement(NvU32 uLayer) const
    {
        nvAssert(uLayer < m_nLayers);
        return m_pLayers[uLayer].m_uFirstEl;
    }
    NvU32 getNextElement(NvU32 uElement) const
    {
        return m_pElements[uElement].m_uNext;
    }
private:
#if ASSERT_ONLY_CODE
    bool dbgIsLayerSane(NvU32 uLayer)
    {
        if ((m_pLayers[uLayer].m_uFirstEl == INVALID_UINT32) != (m_pLayers[uLayer].m_uLastEl == INVALID_UINT32))
        {
            nvAssert(false);
            return false;
        }
        NvU32 uFirst = m_pLayers[uLayer].m_uFirstEl;
        if (uFirst != INVALID_UINT32 && m_pElements[uFirst].m_uPrev != INVALID_UINT32)
        {
            nvAssert(false);
            return false;
        }
        NvU32 uLast = m_pLayers[uLayer].m_uLastEl;
        if (uLast != INVALID_UINT32 && m_pElements[uLast].m_uNext != INVALID_UINT32)
        {
            nvAssert(false);
            return false;
        }
        return true;
    }
#endif
    NvU32 m_nLayers = 0, m_lastLayerId = 0;
    struct Element
    {
        NvU32 m_maxLayerId = 0; // used to quickly determine if given element is present at the given layer
        NvU32 m_uPrev, m_uNext; // elements valid for the layer are linked
    };
    std::vector<Element> m_pElements;
    struct Layer
    {
        NvU32 m_id = 0, m_uFirstEl, m_uLastEl;
    };
    std::array<Layer, 32> m_pLayers;
};