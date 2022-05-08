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
        nvAssert(nElements >= m_nOrigElements); // I only implemented adding new elements - can't shrink the array
        if (nElements <= m_nOrigElements)
            return;
#if ASSERT_ONLY_CODE
        for (NvU32 uLayer = 0; uLayer < m_pLayers.size(); ++uLayer)
        {
            // below I assume that all layers except the ground one are invalid
            nvAssert(uLayer == uGroundLayer || !m_pLayers[uLayer].m_isValid);
        }
#endif
        NvU32 uGrSentinel0 = INVALID_UINT32, uGrSentinel1 = INVALID_UINT32;
        if (m_pLayers[uGroundLayer].m_isValid && hasElements(uGroundLayer))
        {
            uGrSentinel0 = _firstLayerEl(uGroundLayer);
            uGrSentinel1 = _lastLayerEl(uGroundLayer);
        }
        NvU32 uFirstNewElement = m_nOrigElements;
        m_nOrigElements = nElements;
        m_pElements.resize(nElements + m_pLayers.size() * 2);
        if (!m_pLayers[uGroundLayer].m_isValid)
        {
            createLayer(uGroundLayer);
        }
        else
        {
            // copy sentinels to new place
            m_pElements[_lastLayerEl(uGroundLayer)] = m_pElements[uGrSentinel1];
            m_pElements[_firstLayerEl(uGroundLayer)] = m_pElements[uGrSentinel0];
            // fix links from the existing elements to the new sentinel positions
            linkTwoElements(_firstLayerEl(uGroundLayer), m_pElements[_firstLayerEl(uGroundLayer)].m_uNext);
            linkTwoElements(m_pElements[_lastLayerEl(uGroundLayer)].m_uPrev, _lastLayerEl(uGroundLayer));
        }

        // link all new elements together
        m_pElements[uFirstNewElement].m_maxLayerId = m_pLayers[uGroundLayer].m_id;
        for (NvU32 u = uFirstNewElement + 1; u < m_nOrigElements; ++u)
        {
            linkTwoElements(u - 1, u);
            m_pElements[u].m_maxLayerId = m_pLayers[uGroundLayer].m_id;
        }

        // add list of new elements to the ground layer list
        linkTwoElements(m_pElements[_lastLayerEl(uGroundLayer)].m_uPrev, uFirstNewElement);
        linkTwoElements(m_nOrigElements - 1, _lastLayerEl(uGroundLayer));
    }
    bool hasElements(NvU32 uLayer) const
    {
        NvU32 uSentinel0 = _firstLayerEl(uLayer);
        return m_pElements[uSentinel0].m_uNext != uSentinel0 + 1;
    }
    bool hasEverBeenAtLayer(NvU32 uLayer, NvU32 uElement) const
    {
        nvAssert(m_pLayers[uLayer].m_isValid);
        return m_pElements[uElement].m_maxLayerId - m_pLayers[uLayer].m_id < 0x80000000;
    }
    void moveToLayer(NvU32 uLayer, NvU32 uElement)
    {
        nvAssert(m_pLayers[uLayer].m_isValid);
        Element &el = m_pElements[uElement];
        el.m_maxLayerId = m_pLayers[uLayer].m_id;
        // unlinking from whatever layer the element is currently on
        linkTwoElements(el.m_uPrev, el.m_uNext);
        // linking to uLayer
        linkTwoElements(uElement, m_pElements[_firstLayerEl(uLayer)].m_uNext);
        linkTwoElements(_firstLayerEl(uLayer), uElement);
    }
    void createLayer(NvU32 uLayer)
    {
        nvAssert(!m_pLayers[uLayer].m_isValid);
        m_pLayers[uLayer].m_isValid = true;
        m_pLayers[uLayer].m_id = ++m_lastLayerId;
        linkTwoElements(_firstLayerEl(uLayer), _lastLayerEl(uLayer));
    }
    void destroyLayer(NvU32 uLayer)
    {
        nvAssert(m_pLayers[uLayer].m_isValid && !hasElements(uLayer));
        m_pLayers[uLayer].m_isValid = false;
    }
    void moveAllElements(NvU32 uDstLayer, NvU32 uSrcLayer)
    {
        nvAssert(m_pLayers[uDstLayer].m_isValid && m_pLayers[uSrcLayer].m_isValid);
        if (!hasElements(uSrcLayer))
            return;
        // this condition is required for hasEverBeenAtLayer() to work correctly
        nvAssert(m_pLayers[uSrcLayer].m_id - m_pLayers[uDstLayer].m_id < 0x80000000);
        // link the end of src chain to the beginning of dest chain
        linkTwoElements(m_pElements[_lastLayerEl(uSrcLayer)].m_uPrev, m_pElements[_firstLayerEl(uDstLayer)].m_uNext);
        // link first dst sentinel with the beginning of src chain
        linkTwoElements(_firstLayerEl(uDstLayer), m_pElements[_firstLayerEl(uSrcLayer)].m_uNext);
        // link two src sentinels together
        linkTwoElements(_firstLayerEl(uSrcLayer), _lastLayerEl(uSrcLayer));
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
    void linkTwoElements(NvU32 uEl0, NvU32 uEl1)
    {
        m_pElements[uEl0].m_uNext = uEl1;
        m_pElements[uEl1].m_uPrev = uEl0;
    }
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