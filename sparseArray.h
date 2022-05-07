#pragma once

#include <array>
#include <vector>

// SparseHierarchy is used for multi-layer hierarchical simulation. The very first layer contains all elements.
// As the simulation progresses - some elements need to be simulated in more detail than the others. This class
// allows efficient movement of elements between layers.
struct SparseHierarchy
{
    void init(NvU32 nElements, NvU32 uGroundLayer) // uGroundLayer is the layer where all elements are initially
    {
        nvAssert(nElements > 0); // don't need this corner case at the moment
        if (m_nOrigElements != nElements) // size has changed? reinitialize
        {
            m_nOrigElements = nElements;
            // resize to 0 and then back to full size for constructor to be called on each element
            m_pElements.resize(0);
            m_pElements.resize(nElements + m_pLayers.size() * 2);

            for (NvU32 u = uGroundLayer; u <= uGroundLayer; --u)
            {
                m_pLayers[u] = Layer(); // re-initialize the layers completely
                createLayer(u);
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
        nvAssert(m_pLayers[uLayer].m_isValid);
        return m_pElements[uElement].m_maxLayerId >= m_pLayers[uLayer].m_id;
    }
    void moveToLayer(NvU32 uLayer, NvU32 uElement)
    {
        nvAssert(m_pLayers[uLayer].m_isValid);
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
    void createLayer(NvU32 uLayer)
    {
        nvAssert(!m_pLayers[uLayer].m_isValid);
        m_pLayers[uLayer].m_isValid = true;
        m_pLayers[uLayer].m_id = ++m_lastLayerId;
        NvU32 uEl0 = _firstLayerEl(uLayer);
        NvU32 uEl1 = _lastLayerEl(uLayer);
        m_pElements[uEl0].m_uNext = uEl1;
        m_pElements[uEl1].m_uPrev = uEl0;
    }
    void destroyLayer(NvU32 uLayer)
    {
        nvAssert(m_pLayers[uLayer].m_isValid);
        if (hasElements(uLayer))
        {
            // since this layer is destroyed - all its elements are moved to the previous layer
            moveAllElements(uLayer - 1, uLayer);
        }
        m_pLayers[uLayer].m_isValid = false;
    }
    void moveAllElements(NvU32 uDstLayer, NvU32 uSrcLayer)
    {
        nvAssert(m_pLayers[uDstLayer].m_isValid && m_pLayers[uSrcLayer].m_isValid);
        if (!hasElements(uSrcLayer))
            return;
        // this condition is required for hasElement() to work correctly
        nvAssert(m_pLayers[uDstLayer].m_id < m_pLayers[uSrcLayer].m_id);
        NvU32 uSrcFirst = _firstLayerEl(uSrcLayer);
        NvU32 uSrcLast = _lastLayerEl(uSrcLayer);

        NvU32 uEl0 = _firstLayerEl(uDstLayer);
        NvU32 uEl1 = m_pElements[uSrcFirst].m_uNext;
        NvU32 uEl2 = m_pElements[uSrcLast].m_uPrev;
        NvU32 uEl3 = m_pElements[uEl0].m_uNext;
        m_pElements[uEl0].m_uNext = uEl1;
        m_pElements[uEl1].m_uPrev = uEl0;
        m_pElements[uEl2].m_uNext = uEl3;
        m_pElements[uEl3].m_uPrev = uEl2;

        m_pElements[uSrcFirst].m_uNext = uSrcLast;
        m_pElements[uSrcLast].m_uPrev = uSrcFirst;
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
        nvAssert(m_pLayers[uLayer].m_isValid);
        return m_nOrigElements + uLayer * 2;
    }
    NvU32 _lastLayerEl(NvU32 uLayer) const
    {
        nvAssert(m_pLayers[uLayer].m_isValid);
        return m_nOrigElements + uLayer * 2 + 1;
    }
    NvU32 m_lastLayerId = INVALID_UINT32;
    struct Element
    {
        NvU32 m_maxLayerId = 0; // used to quickly determine if given element is present at the given layer
        NvU32 m_uPrev = INVALID_UINT32, m_uNext = INVALID_UINT32; // elements valid for the layer are linked
    };
    std::vector<Element> m_pElements;
    NvU32 m_nOrigElements = 0;
    struct Layer
    {
        bool m_isValid = false;
        NvU32 m_id = INVALID_UINT32;
    };
    std::array<Layer, 32> m_pLayers;
};