#pragma once
#include <cmath>
#include <initializer_list>
#include "mybasics.h"

// Macro to define conversion and subscript operators
#define VECTOR_MEMBERS(T, n) \
			/* Conversions to C arrays of fixed size */ \
			typedef T (&array_t)[n]; \
			operator array_t () \
				{ return m_data; } \
			typedef const T (&const_array_t)[n]; \
			operator const_array_t () const \
				{ return m_data; } \
			/* Subscript operators - built-in subscripts are ambiguous without these */ \
			T & operator [] (int i) \
				{ nvAssert(i < n); return m_data[i]; } \
			const T & operator [] (int i) const \
				{ nvAssert(i < n); return m_data[i]; } \
			/* Conversion to bool is not allowed (otherwise would \
			   happen implicitly through array conversions) */ \
            rtvector() { } \
            rtvector(std::initializer_list<T> inputs) { nvAssert(inputs.size() == n); for (NvU32 u = 0; u < n; ++u) m_data[u] = inputs.begin()[u]; } \
            void set(T s) { for (int i = 0; i < n; ++i) m_data[i] = s; } \
            private: operator bool();

// Generic rtvector struct, providing storage, using partial
// specialization to get names (xyzw) for n <= 4

template <typename T, int n>
struct rtvector
{
    T m_data[n];
    VECTOR_MEMBERS(T, n)
};

#pragma warning(push)
#pragma warning(disable: 4201)	// Nameless struct/union

template <typename T>
struct rtvector<T, 2>
{
    T m_data[2];
    VECTOR_MEMBERS(T, 2)
};

template <typename T>
struct rtvector<T, 3>
{
    T m_data[3];
    VECTOR_MEMBERS(T, 3)
};

template <typename T>
struct rtvector<T, 4>
{
    T m_data[4];
    VECTOR_MEMBERS(T, 4)
};

#pragma warning(pop)
#undef VECTOR_MEMBERS

// Generic maker functions

template <typename T, int n>
rtvector<T, n> makeVector(T a)
{
    rtvector<T, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = a;
    return result;
}

template <typename T, int n, typename U>
rtvector<T, n> makeVector(const U* a)
{
    rtvector<T, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = T(a[i]);
    return result;
}

template <typename T, int n, typename U, int n_from>
rtvector<T, n> makeVector(rtvector<U, n_from> const& a)
{
    auto result = makeVector<T, n>(T(0));
    for (int i = 0; i < min(n, n_from); ++i)
        result[i] = T(a[i]);
    return result;
}

// Concrete vectors, and their maker functions,
// for the most common types and dimensions

#define DEFINE_CONCRETE_VECTORS(type) \
			typedef rtvector<type, 2> type##2; \
			typedef rtvector<type, 3> type##3; \
			typedef rtvector<type, 4> type##4;

DEFINE_CONCRETE_VECTORS(float);
DEFINE_CONCRETE_VECTORS(double);
//DEFINE_CONCRETE_VECTORS(half);	// !!!UNDONE: need to de-constructorize half
DEFINE_CONCRETE_VECTORS(int);
typedef NvU32 uint;
DEFINE_CONCRETE_VECTORS(uint);
//DEFINE_CONCRETE_VECTORS(byte);
//DEFINE_CONCRETE_VECTORS(bool);

#undef DEFINE_CONCRETE_VECTORS

// Overloaded math operators

#define DEFINE_UNARY_OPERATOR(op) \
			template <typename T, int n> \
			rtvector<T, n> operator op (rtvector<T, n> const & a) \
			{ \
				rtvector<T, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = op a[i]; \
				return result; \
			}

#define DEFINE_BINARY_OPERATORS(op) \
			/* Vector-rtvector op */ \
			template <typename T, int n> \
			rtvector<T, n> operator op (rtvector<T, n> const & a, rtvector<T, n> const & b) \
			{ \
				rtvector<T, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = a[i] op b[i]; \
				return result; \
			} \
			/* Scalar-rtvector op */ \
			template <typename T, int n> \
			rtvector<T, n> operator op (T a, rtvector<T, n> const & b) \
			{ \
				rtvector<T, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = a op b[i]; \
				return result; \
			} \
			/* Vector-scalar op */ \
			template <typename T, int n> \
			rtvector<T, n> operator op (rtvector<T, n> const & a, T b) \
			{ \
				rtvector<T, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = a[i] op b; \
				return result; \
			} \
			/* for rtvector<MyUnits<>, n> to work */ \
			template <typename T, int n> \
			rtvector<T, n> operator op (rtvector<T, n> const & a, typename T::T b) \
			{ \
				rtvector<T, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = a[i] op b; \
				return result; \
			}

#define DEFINE_INPLACE_OPERATORS(op) \
			/* Vector-rtvector op */ \
			template <typename T, int n> \
			rtvector<T, n> & operator op (rtvector<T, n> & a, rtvector<T, n> const & b) \
			{ \
				for (int i = 0; i < n; ++i) \
					a[i] op b[i]; \
				return a; \
			} \
			/* Vector-scalar op */ \
			template <typename T, int n> \
			rtvector<T, n> & operator op (rtvector<T, n> & a, T b) \
			{ \
				for (int i = 0; i < n; ++i) \
					a[i] op b; \
				return a; \
			} \
			/* for rtvector<MyUnits<>, n> to work */ \
			template <typename T, int n> \
			rtvector<T, n> & operator op (rtvector<T, n> & a, typename T::T b) \
			{ \
				for (int i = 0; i < n; ++i) \
					a[i] op b; \
				return a; \
			}

#define DEFINE_RELATIONAL_OPERATORS(op) \
			/* Vector-rtvector op */ \
			template <typename T, int n> \
			rtvector<bool, n> operator op (rtvector<T, n> const & a, rtvector<T, n> const & b) \
			{ \
				rtvector<bool, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = a[i] op b[i]; \
				return result; \
			} \
			/* Scalar-rtvector op */ \
			template <typename T, int n> \
			rtvector<bool, n> operator op (T a, rtvector<T, n> const & b) \
			{ \
				rtvector<bool, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = a op b[i]; \
				return result; \
			} \
			/* Vector-scalar op */ \
			template <typename T, int n> \
			rtvector<bool, n> operator op (rtvector<T, n> const & a, T b) \
			{ \
				rtvector<bool, n> result; \
				for (int i = 0; i < n; ++i) \
					result[i] = a[i] op b; \
				return result; \
			}

DEFINE_BINARY_OPERATORS(+);
DEFINE_BINARY_OPERATORS(-);
DEFINE_UNARY_OPERATOR(-);
DEFINE_BINARY_OPERATORS(*);
DEFINE_BINARY_OPERATORS(/ );
DEFINE_BINARY_OPERATORS(&);
DEFINE_BINARY_OPERATORS(| );
DEFINE_BINARY_OPERATORS(^);
DEFINE_UNARY_OPERATOR(!);
DEFINE_UNARY_OPERATOR(~);

DEFINE_INPLACE_OPERATORS(+= );
DEFINE_INPLACE_OPERATORS(-= );
DEFINE_INPLACE_OPERATORS(*= );
DEFINE_INPLACE_OPERATORS(/= );
DEFINE_INPLACE_OPERATORS(&= );
DEFINE_INPLACE_OPERATORS(|= );
DEFINE_INPLACE_OPERATORS(^= );

DEFINE_RELATIONAL_OPERATORS(== );
DEFINE_RELATIONAL_OPERATORS(!= );
DEFINE_RELATIONAL_OPERATORS(< );
DEFINE_RELATIONAL_OPERATORS(> );
DEFINE_RELATIONAL_OPERATORS(<= );
DEFINE_RELATIONAL_OPERATORS(>= );

#undef DEFINE_UNARY_OPERATOR
#undef DEFINE_BINARY_OPERATORS
#undef DEFINE_INPLACE_OPERATORS
#undef DEFINE_RELATIONAL_OPERATORS

// Other math functions    

template <typename T, int n>
T dot(rtvector<T, n> const& a, rtvector<T, n> const& b)
{
    T result(0);
    for (int i = 0; i < n; ++i)
        result += a[i] * b[i];
    return result;
}

template <typename T, int n>
T lengthSquared(rtvector<T, n> const& a)
{
    return dot(a, a);
}

template <typename T, int n>
T length(rtvector<T, n> const& a)
{
    return sqrt(lengthSquared(a));
}

template <typename T, int n>
rtvector<T, n> normalize(rtvector<T, n> const& a)
{
    return a / length(a);
}

template <typename T, int n>
rtvector<T, n> pow(rtvector<T, n> const& a, float p)
{
    rtvector<T, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = ::pow(a[i], p);
    return result;
}

template <typename T, int n>
rtvector<bool, n> isnear(rtvector<T, n> const& a, rtvector<T, n> const& b, float epsilon = util::epsilon)
{
    rtvector<bool, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = isnear(a[i], b[i], epsilon);
    return result;
}

template <typename T, int n>
rtvector<bool, n> isnear(rtvector<T, n> const& a, T b, float epsilon = util::epsilon)
{
    rtvector<bool, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = isnear(a[i], b, epsilon);
    return result;
}

template <typename T, int n>
rtvector<bool, n> isnear(T a, rtvector<T, n> const& b, float epsilon = util::epsilon)
{
    rtvector<bool, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = isnear(a, b[i], epsilon);
    return result;
}

template <typename T, int n>
rtvector<bool, n> isfinite(rtvector<T, n> const& a)
{
    rtvector<bool, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = isfinite(a[i]);
    return result;
}

template <typename T, int n>
rtvector<int, n> round(rtvector<T, n> const& a)
{
    rtvector<int, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = round(a[i]);
    return result;
}

template <typename toT, typename T, int n>
rtvector<toT, n> convTo(rtvector<T, n> const& a)
{
    rtvector<toT, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = (toT)(a[i]);
    return result;
}

template <typename T>
rtvector<T, 3> cross(rtvector<T, 3> const& a, rtvector<T, 3> const& b)
{
    rtvector<T, 3> result =
    {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
    return result;
}

template <typename T>
rtvector<T, 2> orthogonal(rtvector<T, 2> const& a)
{
    rtvector<T, 2> result = { -a.y, a.x };
    return result;
}

template <typename T>
rtvector<T, 3> orthogonal(rtvector<T, 3> const& a)
{
    // Implementation due to Sam Hocevar - see blog post:
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-rtvector-combing-coconuts
    if (abs(a.x) > abs(a.z))
    {
        rtvector<T, 3> result = { -a.y, a.x, T(0) };
        return result;
    }
    else
    {
        rtvector<T, 3> result = { T(0), -a.z, a.y };
        return result;
    }
}

// Utilities for bool vectors

template <int n>
bool any(rtvector<bool, n> const& a)
{
    bool result = false;
    for (int i = 0; i < n; ++i)
        result = result || a[i];
    return result;
}

template <int n>
bool all(rtvector<bool, n> const& a)
{
    bool result = true;
    for (int i = 0; i < n; ++i)
        result = result && a[i];
    return result;
}

template <typename T, int n>
rtvector<T, n> select(rtvector<bool, n> const& cond, rtvector<T, n> const& a, rtvector<T, n> const& b)
{
    rtvector<T, n> result;
    for (int i = 0; i < n; ++i)
        result[i] = cond[i] ? a[i] : b[i];
    return result;
}

template <typename T, int n>
rtvector<T, n> vmin(rtvector<T, n> const& a, rtvector<T, n> const& b)
{
    return select(a < b, a, b);
}

template <typename T, int n>
rtvector<T, n> vmax(rtvector<T, n> const& a, rtvector<T, n> const& b)
{
    return select(a < b, b, a);
}

template <typename T, int n>
rtvector<T, n> abs(rtvector<T, n> const& a)
{
    return select(a < T(0), -a, a);
}

template <typename T, int n>
rtvector<T, n> saturate(rtvector<T, n> const& value)
{
    return clamp(value, makeVector<T, n>(0), makeVector<T, n>(1));
}

template <typename T, int n>
T minComponent(rtvector<T, n> const& a)
{
    T result = a[0];
    for (int i = 1; i < n; ++i)
        result = min(result, a[i]);
    return result;
}

template <typename T, int n>
T maxComponent(rtvector<T, n> const& a)
{
    T result = a[0];
    for (int i = 1; i < n; ++i)
        result = max(result, a[i]);
    return result;
}
