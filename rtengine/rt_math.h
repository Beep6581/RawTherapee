#ifndef _MYMATH_
#define _MYMATH_
#include <cmath>
#include <algorithm>
#include <type_traits>


namespace rtengine
{
static const int MAXVAL = 0xffff;
static const float MAXVALF = float(MAXVAL);  // float version of MAXVAL
static const double MAXVALD = double(MAXVAL); // double version of MAXVAL

template <typename _Tp>
inline const _Tp SQR (_Tp x)
{
//      return std::pow(x,2); Slower than:
    return (x * x);
}

template<typename _Tp>
inline const _Tp& min(const _Tp& a, const _Tp& b)
{
    return std::min(a, b);
}

template<typename _Tp>
inline const _Tp& max(const _Tp& a, const _Tp& b)
{
    return std::max(a, b);
}


template<typename _Tp>
inline const _Tp LIM(const _Tp& a, const _Tp& b, const _Tp& c)
{
    return std::max(b, std::min(a, c));
}

template<typename _Tp>
inline const _Tp LIM01(const _Tp& a)
{
    return std::max(_Tp(0), std::min(a, _Tp(1)));
}

template<typename _Tp>
inline const _Tp ULIM(const _Tp& a, const _Tp& b, const _Tp& c)
{
    return ((b < c) ? LIM(a, b, c) : LIM(a, c, b));
}

template<typename _Tp>
inline const _Tp CLIP(const _Tp& a)
{
    return LIM(a, static_cast<_Tp>(0), static_cast<_Tp>(MAXVAL));
}


template<typename _Tp>
inline const _Tp& min(const _Tp& a, const _Tp& b, const _Tp& c)
{
    return std::min(c, std::min(a, b));
}

template<typename _Tp>
inline const _Tp& max(const _Tp& a, const _Tp& b, const _Tp& c)
{
    return std::max(c, std::max(a, b));
}

template<typename _Tp>
inline const _Tp& min(const _Tp& a, const _Tp& b, const _Tp& c, const _Tp& d)
{
    return std::min(d, std::min(c, std::min(a, b)));
}

template<typename _Tp>
inline const _Tp& max(const _Tp& a, const _Tp& b, const _Tp& c, const _Tp& d)
{
    return std::max(d, std::max(c, std::max(a, b)));
}

template<typename _Tp>
inline const _Tp intp(const _Tp a, const _Tp b, const _Tp c)
{
    // calculate a * b + (1 - a) * c
    // following is valid:
    // intp(a, b+x, c+x) = intp(a, b, c) + x
    // intp(a, b*x, c*x) = intp(a, b, c) * x
    return a * (b-c) + c;
}

template<typename T>
T norm1(const T& x, const T& y)
{
    return std::abs(x) + std::abs(y);
}

template<typename T>
T norm2(const T& x, const T& y)
{
    return std::sqrt(x * x + y * y);
}

template<typename T>
T norminf(const T& x, const T& y)
{
    return std::max(std::abs(x), std::abs(y));
}

template<
    typename T,
    typename = typename std::enable_if<std::is_integral<T>::value>::type,
    typename = typename std::enable_if<std::is_unsigned<T>::value>::type
>
T adds(T a, T b)
{
    T res = a + b;
    res |= -(res < a);
    return res;
}

}
#endif
