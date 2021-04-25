#pragma once

typedef unsigned int NvU32;
typedef unsigned long long NvU64;

#ifdef NDEBUG
#define ASSERT_ONLY_CODE 0
#define nvAssert(x)
#else
#define ASSERT_ONLY_CODE 1
#define nvAssert(x) if (!(x)) { __debugbreak(); }
#endif

#define ARRAY_ELEMENT_COUNT(x) (sizeof(x) / sizeof(x[0]))

template <class T> inline T sqr(const T & a) { return a * a; }

template <class T>
class MyNumericLimits
{
};
template <>
struct MyNumericLimits<float>
{
    typedef NvU32 SubstituteUintType;
    static NvU32 topBit() { return 0x80000000U; }
};
template <>
struct MyNumericLimits<double>
{
    typedef NvU64 SubstituteUintType;
    static NvU64 topBit() { return 0x8000000000000000ULL; }
};
template <class floatType>
inline bool aboutEqual(floatType f1, floatType f2)
{
    typedef MyNumericLimits<floatType>::SubstituteUintType uintType;
    const uintType topBit = MyNumericLimits<floatType>::topBit();
    uintType u1 = (uintType &)f1, u2 = (uintType &)f2;
    bool haveDifferentSign = (u1 ^ u2) & topBit;
    u1 &= ~topBit;
    u2 &= ~topBit;
    uintType dist = haveDifferentSign ? (u1 + u2) : (u1 > u2 ? u1 - u2 : u2 - u1);
    return dist < 20;
}