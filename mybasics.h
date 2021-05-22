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
class MyNumeric
{
};
template <>
struct MyNumeric<float>
{
    typedef NvU32 SubstituteUintType;
    static NvU32 topBit() { return 0x80000000U; }
};
template <>
struct MyNumeric<double>
{
    typedef NvU64 SubstituteUintType;
    static NvU64 topBit() { return 0x8000000000000000ULL; }
};
template <class floatType>
inline bool aboutEqual(floatType f1, floatType f2)
{
    typedef MyNumeric<floatType>::SubstituteUintType uintType;
    const uintType topBit = MyNumeric<floatType>::topBit();
    uintType u1 = (uintType &)f1, u2 = (uintType &)f2;
    bool haveDifferentSign = (u1 ^ u2) & topBit;
    u1 &= ~topBit;
    u2 &= ~topBit;
    uintType dist = haveDifferentSign ? (u1 + u2) : (u1 > u2 ? u1 - u2 : u2 - u1);
    return dist < 20;
}

template <class T>
void nvSwap(T& a, T& b)
{
    T c;
    c = a;
    a = b;
    b = c;
}