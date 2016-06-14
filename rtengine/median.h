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
T median(std::array<T, N> array)
{
    const typename std::array<T, N>::iterator middle = array.begin() + array.size() / 2;
    std::nth_element(array.begin(), middle, array.end());

    return
        array.size() % 2
            ? ((*middle + *std::min_element(middle + 1, array.end())) / static_cast<T>(2))
            : *middle;
}

template<typename T, typename... ARGS>
T median(T arg, ARGS... args)
{
    return median(std::array<T, sizeof...(args) + 1>{std::move(arg), std::move(args)...});
}

template<typename T>
inline T median3(T a, T b, T c)
{
    return std::max(std::min(a, b), std::min(c, std::max(a, b)));
}

template<>
inline vfloat median3(vfloat a, vfloat b, vfloat c)
{
    return vmaxf(vminf(a, b), vminf(c, vmaxf(a, b)));
}

// See http://stackoverflow.com/questions/480960/code-to-calculate-median-of-five-in-c-sharp
template<typename T>
inline T median5(T a, T b, T c, T d, T e)
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
inline vfloat median5(vfloat a, vfloat b, vfloat c, vfloat d, vfloat e)
{
    const vfloat f = vmaxf(vminf(a, b), vminf(c, d));
    const vfloat g = vminf(vmaxf(a, b), vmaxf(c, d));
    return median3(e, f, g);
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
