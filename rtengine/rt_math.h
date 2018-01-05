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
constexpr double RT_PI_180 = 0.017453292519943295769; // pi/180
constexpr double RT_1_PI = 0.31830988618379067154; // 1/pi
constexpr double RT_2_PI = 0.63661977236758134308; // 2/pi
constexpr double RT_SQRT1_2 = 0.70710678118654752440; // 1/sqrt(2)

constexpr double RT_INFINITY = std::numeric_limits<double>::infinity();
constexpr double RT_NAN = std::numeric_limits<double>::quiet_NaN();

constexpr float RT_PI_F = RT_PI;
constexpr float RT_PI_F_2 = RT_PI_2;
constexpr float RT_PI_F_180 = RT_PI_180;

constexpr float RT_INFINITY_F = std::numeric_limits<float>::infinity();
constexpr float RT_NAN_F = std::numeric_limits<float>::quiet_NaN();

template<typename T>
constexpr T SQR(T x)
{
    return x * x;
}

template<typename T>
constexpr const T& min(const T& a)
{
    return a;
}

template<typename T>
constexpr const T& min(const T& a, const T& b)
{
    return b < a ? b : a;
}

template<typename T, typename... ARGS>
constexpr const T& min(const T& a, const T& b, const ARGS&... args)
{
    return min(min(a, b), min(args...));
}

template<typename T>
constexpr const T& max(const T& a)
{
    return a;
}

template<typename T>
constexpr const T& max(const T& a, const T& b)
{
    return a < b ? b : a;
}

template<typename T, typename... ARGS>
constexpr const T& max(const T& a, const T& b, const ARGS&... args)
{
    return max(max(a, b), max(args...));
}

template<typename T>
constexpr const T& LIM(const T& a, const T& b, const T& c)
{
    return max(b, min(a, c));
}

template<typename T>
constexpr T LIM01(const T& a)
{
    return max(T(0), min(a, T(1)));
}

template<typename T>
constexpr T CLIP(const T& a)
{
    return LIM(a, static_cast<T>(0), static_cast<T>(MAXVAL));
}

template <typename T>
constexpr T SGN(const T& a)
{
    // returns -1 for a < 0, 0 for a = 0 and +1 for a > 0
    return (T(0) < a) - (a < T(0));
}

template<typename T>
constexpr T intp(T a, T b, T c)
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
    return max(std::abs(x), std::abs(y));
}

constexpr int float2uint16range(float d)
{
    // clips input to [0;65535] and rounds
    return CLIP(d) + 0.5f;
}

constexpr std::uint8_t uint16ToUint8Rounded(std::uint16_t i)
{
    return ((i + 128) - ((i + 128) >> 8)) >> 8;
}

}
