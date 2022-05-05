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
// IMPORTANT: internally the first and the last elements of each layer are sentinels. It greatly simplifies the
// logic because otherwise layer being empty is a corner case that has to be taken into account in many places
struct SparseHierarchy
{
    void init(NvU32 nElements, NvU32 uGroundLayer)
    {
        nvAssert(nElements > 0);
        if (m_nOrigElements != nElements) // size has changed? reinitialize
        {
            m_nOrigElements = nElements;
            m_pElements.resize(nElements + m_pLayers.size() * 2);

            m_nLayers = 0;
            for (NvU32 u = 0; u <= uGroundLayer; ++u)
            {
                notifyLayerCreated(u);
            }

            // link all elements into list
            for (NvU32 u = 1; u < m_nOrigElements; ++u)
            {
                m_pElements[u - 1].m_uNext = u;
                m_pElements[u].m_uPrev = u - 1;
                m_pElements[u].m_maxLayerId = m_pLayers[uGroundLayer].m_id;
            }

            // put all elements onto ground layer
            m_pElements[_firstLayerEl(uGroundLayer)].m_uNext = 0;
            m_pElements[0].m_uPrev = _firstLayerEl(uGroundLayer);
            m_pElements[_lastLayerEl(uGroundLayer)].m_uPrev = m_nOrigElements - 1;
            m_pElements[m_nOrigElements - 1].m_uNext = _lastLayerEl(uGroundLayer);
        }
    }
    bool hasElements(NvU32 uLayer) const
    {
        return m_pElements[_firstLayerEl(uLayer)].m_uNext < m_nOrigElements;
    }
    bool hasElement(NvU32 uLayer, NvU32 uElement) const
    {
        nvAssert(uLayer < m_nLayers);
        return m_pElements[uElement].m_maxLayerId >= m_pLayers[uLayer].m_id;
    }
    void moveToLayer(NvU32 uLayer, NvU32 uElement)
    {
        if (hasElement(uLayer, uElement)) return;
        auto& el = m_pElements[uElement];
        el.m_maxLayerId = m_pLayers[uLayer].m_id;
        // unlinking from whatever layer the element is currently on
        m_pElements[el.m_uPrev].m_uNext = el.m_uNext;
        m_pElements[el.m_uNext].m_uPrev = el.m_uPrev;
        // linking to uLayer
        NvU32 uEl0 = _firstLayerEl(uLayer);
        auto& el0 = m_pElements[uEl0];
        NvU32 uEl2 = el0.m_uNext;
        el0.m_uNext = uElement;
        el.m_uPrev = uEl0;
        el.m_uNext = uEl2;
        m_pElements[uEl2].m_uPrev = uElement;
    }
    void notifyLayerCreated(NvU32 uLayer)
    {
        nvAssert(uLayer == m_nLayers && m_nLayers < m_pLayers.size()); // each time must increase layer index by one
        m_nLayers = uLayer + 1;
        m_pLayers[uLayer].m_id = ++m_lastLayerId;
        NvU32 uEl0 = _firstLayerEl(uLayer);
        NvU32 uEl1 = _lastLayerEl(uLayer);
        m_pElements[uEl0].m_uNext = uEl1;
        m_pElements[uEl1].m_uPrev = uEl0;
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
            if (hasElements(uLayer))
            {
                NvU32 uEl0 = _firstLayerEl(uPrevLayer);
                NvU32 uEl1 = m_pElements[_firstLayerEl(uLayer)].m_uNext;
                NvU32 uEl2 = m_pElements[_lastLayerEl(uLayer)].m_uPrev;
                NvU32 uEl3 = m_pElements[uEl0].m_uNext;

                m_pElements[uEl0].m_uNext = uEl1;
                m_pElements[uEl1].m_uPrev = uEl0;
                m_pElements[uEl2].m_uNext = uEl3;
                m_pElements[uEl3].m_uPrev = uEl2;
            }
        }
        m_nLayers = uLayer;
    }
    NvU32 getFirstLayerElement(NvU32 uLayer) const
    {
        // first element is a sentinel - so return the one after that
        return m_pElements[_firstLayerEl(uLayer)].m_uNext;
    }
    NvU32 getNextElement(NvU32 uElement) const
    {
        return m_pElements[uElement].m_uNext;
    }
private:
    NvU32 _firstLayerEl(NvU32 uLayer) const
    {
        nvAssert(uLayer < m_nLayers);
        return m_nOrigElements + uLayer * 2;
    }
    NvU32 _lastLayerEl(NvU32 uLayer) const
    {
        nvAssert(uLayer < m_nLayers);
        return m_nOrigElements + uLayer * 2 + 1;
    }
    NvU32 m_nLayers = 0, m_lastLayerId = 0;
    struct Element
    {
        NvU32 m_maxLayerId = 0; // used to quickly determine if given element is present at the given layer
        NvU32 m_uPrev = INVALID_UINT32, m_uNext = INVALID_UINT32; // elements valid for the layer are linked
    };
    std::vector<Element> m_pElements;
    NvU32 m_nOrigElements = 0;
    struct Layer
    {
        NvU32 m_id = 0;
    };
    std::array<Layer, 32> m_pLayers;
};