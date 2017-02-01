#pragma once

#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdint>

namespace rtengine
{

constexpr int MAXVAL = 0xffff;
constexpr float MAXVALF = static_cast<float>(MAXVAL);  // float version of MAXVAL
constexpr double MAXVALD = static_cast<double>(MAXVAL); // double version of MAXVAL

constexpr double RT_PI = 3.14159265358979323846; // pi
constexpr double RT_PI_2 = 1.57079632679489661923; // pi/2
constexpr double RT_1_PI = 0.31830988618379067154; // 1/pi
constexpr double RT_2_PI = 0.63661977236758134308; // 2/pi
constexpr double RT_SQRT1_2 = 0.70710678118654752440; // 1/sqrt(2)

constexpr double RT_INFINITY = std::numeric_limits<double>::infinity();
constexpr double RT_NAN = std::numeric_limits<double>::quiet_NaN();

constexpr float RT_PI_F = RT_PI;
constexpr float RT_PI_F_2 = RT_PI_2;

constexpr float RT_INFINITY_F = std::numeric_limits<float>::infinity();
constexpr float RT_NAN_F = std::numeric_limits<float>::quiet_NaN();

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
