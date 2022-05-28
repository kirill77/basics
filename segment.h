#pragma once

#include "vectors.h"

template <class T>
T distSqrToSegment(const rtvector<T, 3>& v, const rtvector<T, 3>& a, const rtvector<T, 3>& b)
{
    rtvector<T, 3> ab = b - a;
    rtvector<T, 3> av = v - a;

    if (dot(av, ab) <= 0.0)  // Point is lagging behind start of the segment, so perpendicular distance is not viable.
        return dot(av, av);  // Use distance to start of segment instead.

    rtvector<T, 3> bv = v - b;

    if (dot(bv, ab) >= 0.0)           // Point is advanced past the end of the segment, so perpendicular distance is not viable.
        return dot(bv, bv);         // Use distance to end of the segment instead.

    rtvector<T, 3> c = cross(ab, av);
    return dot(c, c) / dot(ab, ab);       // Perpendicular distance of point to segment.
}