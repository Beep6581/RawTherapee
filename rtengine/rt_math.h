#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace rtengine
{

constexpr int MAXVAL = 0xffff;
constexpr float MAXVALF = static_cast<float>(MAXVAL);  // float version of MAXVAL
constexpr double MAXVALD = static_cast<double>(MAXVAL); // double version of MAXVAL

template <typename _Tp>
inline _Tp SQR (_Tp x)
{
//      return std::pow(x,2); Slower than:
    return x * x;
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
inline const _Tp& LIM(const _Tp& a, const _Tp& b, const _Tp& c)
{
    return std::max(b, std::min(a, c));
}

template<typename _Tp>
inline _Tp LIM01(const _Tp& a)
{
    return std::max(_Tp(0), std::min(a, _Tp(1)));
}

template<typename _Tp>
inline _Tp CLIP(const _Tp& a)
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
inline const _Tp& max(const _Tp& a, const _Tp& b, const _Tp& c, const _Tp& d, const _Tp& e)
{
    return max(max(a,b,c),std::max(d,e));
}

template<typename _Tp>
inline const _Tp& max(const _Tp& a, const _Tp& b, const _Tp& c, const _Tp& d, const _Tp& e, const _Tp& f, const _Tp& g)
{
    return max(max(a,b,c,d),max(e,f,g));
}


template<typename _Tp>
inline _Tp intp(_Tp a, _Tp b, _Tp c)
{
    // calculate a * b + (1 - a) * c
    // following is valid:
    // intp(a, b+x, c+x) = intp(a, b, c) + x
    // intp(a, b*x, c*x) = intp(a, b, c) * x
    return a * (b - c) + c;
}

template<typename T>
inline T norm1(const T& x, const T& y)
{
    return std::abs(x) + std::abs(y);
}

template<typename T>
inline T norm2(const T& x, const T& y)
{
    return std::sqrt(x * x + y * y);
}

template< typename T >
inline T norminf(const T& x, const T& y)
{
    return std::max(std::abs(x), std::abs(y));
}

inline int float2uint16range(float d) // clips input to [0;65535] and rounds
{
    d = CLIP(d); // clip to [0;65535]
    return d + 0.5f;
}

constexpr std::uint8_t uint16ToUint8Rounded(std::uint16_t i)
{
    return ((i + 128) - ((i + 128) >> 8)) >> 8;
}

}
