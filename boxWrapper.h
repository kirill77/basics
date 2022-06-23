#pragma once

// class used to wrap coordinates and directions so that everything stays inside the simulation boundind box
template <class T>
struct BoxWrapper : public BBox3<MyUnits<T>>
{
    BoxWrapper(MyUnits<T> fBoxSide = MyUnits<T>(1.))
    {
        m_fBoxSize = fBoxSide;
        m_fHalfBoxSize = m_fBoxSize / 2;
        this->m_vMin = makeVector<MyUnits<T>, 3>(MyUnits<T>(-m_fHalfBoxSize));
        this->m_vMax = makeVector<MyUnits<T>, 3>(MyUnits<T>( m_fHalfBoxSize));
    }
    // if the atom exits bounding box, it enters from the other side
    MyUnits3<T> wrapThePos(const MyUnits3<T>& vPos) const
    {
        MyUnits3<T> vNewPos = vPos;
        for (NvU32 uDim = 0; uDim < 3; ++uDim)
        {
            if (vNewPos[uDim] < this->m_vMin[uDim])
            {
                auto fOvershoot = (this->m_vMin[uDim] - vNewPos[uDim]);
                int nBoxSizes = 1 + (int)(fOvershoot / m_fBoxSize);
                vNewPos[uDim] += m_fBoxSize * nBoxSizes;
                nvAssert(this->m_vMin[uDim] <= vNewPos[uDim] && vNewPos[uDim] <= this->m_vMax[uDim]);
                continue;
            }
            if (vNewPos[uDim] > this->m_vMax[uDim])
            {
                auto fOvershoot = (vNewPos[uDim] - this->m_vMax[uDim]);
                int nBoxSizes = 1 + (int)(fOvershoot / m_fBoxSize);
                vNewPos[uDim] -= m_fBoxSize * nBoxSizes;
                nvAssert(this->m_vMin[uDim] <= vNewPos[uDim] && vNewPos[uDim] <= this->m_vMax[uDim]);
            }
        }
        nvAssert(this->includes(vNewPos)); // atom must be inside the bounding box
        return vNewPos;
    }
    rtvector<MyUnits<T>, 3> computeDir(const rtvector<MyUnits<T>, 3> &vPos1, const rtvector<MyUnits<T>, 3>& vPos2) const
    {
        rtvector<MyUnits<T>, 3> vOutDir = vPos1 - vPos2;
        for (NvU32 uDim = 0; uDim < 3; ++uDim) // particles positions must wrap around the boundary of bounding box
        {
            if (vOutDir[uDim] < -m_fHalfBoxSize) vOutDir[uDim] += m_fBoxSize;
            else if (vOutDir[uDim] > m_fHalfBoxSize) vOutDir[uDim] -= m_fBoxSize;
        }
        return vOutDir;
    }
    bool areBoxesFartherThan(const BBox3<MyUnits<T>>& b1, const BBox3<MyUnits<T>>& b2, MyUnits<T> _fDistSqr)
    {
        T fDistSqr = 0;
        for (NvU32 uDim = 0; uDim < 3; ++uDim)
        {
            // compute distance between centers
            T fCenterDist = abs((b1.m_vMax[uDim] + b1.m_vMin[uDim]) - (b2.m_vMax[uDim] + b2.m_vMin[uDim])) / 2;
            // wrap the distance around
            fCenterDist -= m_fBoxSize * (int)(fCenterDist / m_fBoxSize);
            if (fCenterDist > m_fHalfBoxSize) fCenterDist = m_fBoxSize - fCenterDist;
            // subtract half box sizes
            fCenterDist -= ((b1.m_vMax[uDim] - b1.m_vMin[uDim]) + (b2.m_vMax[uDim] - b2.m_vMin[uDim])) / 2;
            fCenterDist = std::max((T)0, fCenterDist);
            fDistSqr += sqr(fCenterDist);
            if (fDistSqr >= _fDistSqr)
                return true;
        }
        return false;
    }

private:
    T m_fBoxSize, m_fHalfBoxSize;
};