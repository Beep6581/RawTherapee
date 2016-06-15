/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <array>
#include <algorithm>

#include "opthelper.h"

template<typename T, std::size_t N>
inline T median(std::array<T, N> array)
{
    const typename std::array<T, N>::iterator middle = array.begin() + N / 2;
    std::nth_element(array.begin(), middle, array.end());

    return
        N % 2
            ? *middle
            : ((*middle + *std::max_element(array.begin(), middle)) / static_cast<T>(2));
}

template<typename T, typename... ARGS>
inline T median(T arg, ARGS... args)
{
    return median(std::array<T, sizeof...(args) + 1>{std::move(arg), std::move(args)...});
}

template<typename T>
inline T median(T a, T b, T c)
{
    return std::max(std::min(a, b), std::min(c, std::max(a, b)));
}

template<>
inline vfloat median(vfloat a, vfloat b, vfloat c)
{
    return vmaxf(vminf(a, b), vminf(c, vmaxf(a, b)));
}

// See http://stackoverflow.com/questions/480960/code-to-calculate-median-of-five-in-c-sharp
template<typename T>
inline T median(T a, T b, T c, T d, T e)
{
    if (b < a) {
        std::swap(a, b);
    }
    if (d < c) {
        std::swap(c, d);
    }

    if (c < a) {
        std::swap(b, d);
        c = a;
    }

    a = e;

    if (b < a) {
        std::swap(a, b);
    }

    if (a < c) {
        std::swap(b, d);
        a = c;
    }

    return std::min(a, d);
}

template<>
inline vfloat median(vfloat a, vfloat b, vfloat c, vfloat d, vfloat e)
{
    const vfloat f = vmaxf(vminf(a, b), vminf(c, d));
    const vfloat g = vminf(vmaxf(a, b), vmaxf(c, d));
    return median(e, f, g);
}

// See http://www.cs.hut.fi/~cessu/selection/V_7_4
// Hand unrolled algorithm by Flössie ;)
template<typename T>
inline T median(T a, T b, T c, T d, T e, T f, T g)
{
    if (b < a) {
        std::swap(a, b);
    }

    if (d < c) {
        std::swap(b, c);
        std::swap(b, d);
    } else {
        std::swap(b, c);
    }

    if (b < a) {
        std::swap(a, b);
    } else {
        std::swap(c, d);
    }

    if (e < d) {
        std::swap(c, d);
        std::swap(c, e);

        if (c < b) {
            std::swap(b, c);

            if (g < f) {
                std::swap(d, g);
                std::swap(e, g);
                std::swap(f, g);
            } else {
                std::swap(d, f);
                std::swap(e, f);
            }

            if (g < c) {
                std::swap(a, d);

                if (d < b) {
                    std::swap(a, b);
                } else {
                    std::swap(a, d);
                    b = d;
                }

                if (g < f) {
                    return std::max(a, g);
                }

                return std::max(b, f);
            }

            if (f < d) {
                std::swap(d, f);
            }

            if (d < c) {
                return std::min(c, f);
            }

            return std::min(d, e);
        }

        if (d < c) {
            std::swap(c, e);

            if (f < b) {
                if (g < b) {
                    return b;
                }

                return std::min(d, g);
            }

            if (g < f) {
                std::swap(e, g);
                std::swap(f, g);
            } else {
                std::swap(e, f);
            }

            if (e < d) {
                return std::min(d, g);
            }

            return std::min(e, f);
        }

        if (g < f) {
            std::swap(f, g);
        }

        std::swap(e, f);

        if (e < c) {
            if (g < b) {
                return b;
            }

            return std::min(c, g);
        }

        if (d < e) {
            return std::min(d, f);
        }

        return std::min(e, f);
    }

    if (d < b) {
        std::swap(b, d);
    } else {
        std::swap(c, e);
    }

    if (e < d) {
        std::swap(d, e);

        if (f < b) {
            if (g < b) {
                return b;
            }

            return std::min(d, g);
        }

        if (g < f) {
            std::swap(e, g);
            std::swap(f, g);
        } else {
            std::swap(e, f);
        }

        if (e < d) {
            return std::min(d, g);
        }

        return std::min(e, f);
    }

    if (g < f) {
        std::swap(c, g);
        std::swap(f, g);
    } else {
        std::swap(c, f);
    }

    if (c < d) {
        if (g < b) {
            return b;
        }

        return std::min(d, g);
    }

    if (e < c) {
        return std::min(e, f);
    }

    return std::min(c, f);
}

// middle 4 of 6 elements,
#define MIDDLE4OF6(s0,s1,s2,s3,s4,s5,d0,d1,d2,d3,d4,d5,temp) \
{\
d1 = std::min(s1,s2);\
d2 = std::max(s1,s2);\
d0 = std::min(s0,d2);\
d2 = std::max(s0,d2);\
temp = std::min(d0,d1);\
d1 = std::max(d0,d1);\
d0 = temp;\
d4 = std::min(s4,s5);\
d5 = std::max(s4,s5);\
d3 = std::min(s3,d5);\
d5 = std::max(s3,d5);\
temp = std::min(d3,d4);\
d4 = std::max(d3,d4);\
d3 = std::max(d0,temp);\
d2 = std::min(d2,d5);\
}

// middle 4 of 6 elements, vectorized
#define VMIDDLE4OF6(s0,s1,s2,s3,s4,s5,d0,d1,d2,d3,d4,d5,temp) \
{\
d1 = vminf(s1,s2);\
d2 = vmaxf(s1,s2);\
d0 = vminf(s0,d2);\
d2 = vmaxf(s0,d2);\
temp = vminf(d0,d1);\
d1 = vmaxf(d0,d1);\
d0 = temp;\
d4 = vminf(s4,s5);\
d5 = vmaxf(s4,s5);\
d3 = vminf(s3,d5);\
d5 = vmaxf(s3,d5);\
temp = vminf(d3,d4);\
d4 = vmaxf(d3,d4);\
d3 = vmaxf(d0,temp);\
d2 = vminf(d2,d5);\
}


#define MEDIAN7(s0,s1,s2,s3,s4,s5,s6,t0,t1,t2,t3,t4,t5,t6,median) \
{\
t0 = std::min(s0,s5);\
t5 = std::max(s0,s5);\
t3 = std::max(t0,s3);\
t0 = std::min(t0,s3);\
t1 = std::min(s1,s6);\
t6 = std::max(s1,s6);\
t2 = std::min(s2,s4);\
t4 = std::max(s2,s4);\
t1 = std::max(t0,t1);\
median = std::min(t3,t5);\
t5 = std::max(t3,t5);\
t3 = median;\
median = std::min(t2,t6);\
t6 = std::max(t2,t6);\
t3 = std::max(median,t3);\
t3 = std::min(t3,t6);\
t4 = std::min(t4,t5);\
median = std::min(t1,t4);\
t4 = std::max(t1,t4);\
t3 = std::max(median,t3);\
median = std::min(t3,t4);\
}

#define VMEDIAN7(s0,s1,s2,s3,s4,s5,s6,t0,t1,t2,t3,t4,t5,t6,median) \
{\
t0 = vminf(s0,s5);\
t5 = vmaxf(s0,s5);\
t3 = vmaxf(t0,s3);\
t0 = vminf(t0,s3);\
t1 = vminf(s1,s6);\
t6 = vmaxf(s1,s6);\
t2 = vminf(s2,s4);\
t4 = vmaxf(s2,s4);\
t1 = vmaxf(t0,t1);\
median = vminf(t3,t5);\
t5 = vmaxf(t3,t5);\
t3 = median;\
median = vminf(t2,t6);\
t6 = vmaxf(t2,t6);\
t3 = vmaxf(median,t3);\
t3 = vminf(t3,t6);\
t4 = vminf(t4,t5);\
median = vminf(t1,t4);\
t4 = vmaxf(t1,t4);\
t3 = vmaxf(median,t3);\
median = vminf(t3,t4);\
}
