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

template<typename T>
inline T median(std::array<T, 3> array)
{
    return std::max(std::min(array[0], array[1]), std::min(array[2], std::max(array[0], array[1])));
}

template<>
inline vfloat median(std::array<vfloat, 3> array)
{
    return vmaxf(vminf(array[0], array[1]), vminf(array[2], vmaxf(array[0], array[1])));
}

template<typename T>
inline T median(std::array<T, 5> array)
{
    T tmp = std::min(array[0], array[1]);
    array[1] = std::max(array[0], array[1]);
    array[0] = tmp;
    tmp = std::min(array[3], array[4]);
    array[4] = std::max(array[3], array[4]);
    array[3] = std::max(array[0], tmp);
    array[1] = std::min(array[1], array[4]);
    tmp = std::min(array[1], array[2]);
    array[2] = std::max(array[1], array[2]);
    array[1] = tmp;
    tmp = std::min(array[2], array[3]);
    return std::max(array[1], tmp);
}

template<>
inline vfloat median(std::array<vfloat, 5> array)
{
    vfloat tmp = vminf(array[0], array[1]);
    array[1] = vmaxf(array[0], array[1]);
    array[0] = tmp;
    tmp = vminf(array[3], array[4]);
    array[4] = vmaxf(array[3], array[4]);
    array[3] = vmaxf(array[0], tmp);
    array[1] = vminf(array[1], array[4]);
    tmp = vminf(array[1], array[2]);
    array[2] = vmaxf(array[1], array[2]);
    array[1] = tmp;
    tmp = vminf(array[2], array[3]);
    return vmaxf(array[1], tmp);
}

template<typename T>
inline T median(std::array<T, 7> array)
{
    T tmp = std::min(array[0], array[5]);
    array[5] = std::max(array[0], array[5]);
    array[0] = tmp;
    tmp = std::min(array[0], array[3]);
    array[3] = std::max(array[0], array[3]);
    array[0] = tmp;
    tmp = std::min(array[1], array[6]);
    array[6] = std::max(array[1], array[6]);
    array[1] = tmp;
    tmp = std::min(array[2], array[4]);
    array[4] = std::max(array[2], array[4]);
    array[2] = tmp;
    array[1] = std::max(array[0], array[1]);
    tmp = std::min(array[3], array[5]);
    array[5] = std::max(array[3], array[5]);
    array[3] = tmp;
    tmp = std::min(array[2], array[6]);
    array[6] = std::max(array[2], array[6]);
    array[3] = std::max(tmp, array[3]);
    array[3] = std::min(array[3], array[6]);
    tmp = std::min(array[4], array[5]);
    array[4] = std::max(array[1], tmp);
    tmp = std::min(array[1], tmp);
    array[3] = std::max(tmp, array[3]);
    return std::min(array[3], array[4]);
}

template<>
inline vfloat median(std::array<vfloat, 7> array)
{
    vfloat tmp = vminf(array[0], array[5]);
    array[5] = vmaxf(array[0], array[5]);
    array[0] = tmp;
    tmp = vminf(array[0], array[3]);
    array[3] = vmaxf(array[0], array[3]);
    array[0] = tmp;
    tmp = vminf(array[1], array[6]);
    array[6] = vmaxf(array[1], array[6]);
    array[1] = tmp;
    tmp = vminf(array[2], array[4]);
    array[4] = vmaxf(array[2], array[4]);
    array[2] = tmp;
    array[1] = vmaxf(array[0], array[1]);
    tmp = vminf(array[3], array[5]);
    array[5] = vmaxf(array[3], array[5]);
    array[3] = tmp;
    tmp = vminf(array[2], array[6]);
    array[6] = vmaxf(array[2], array[6]);
    array[3] = vmaxf(tmp, array[3]);
    array[3] = vminf(array[3], array[6]);
    tmp = vminf(array[4], array[5]);
    array[4] = vmaxf(array[1], tmp);
    tmp = vminf(array[1], tmp);
    array[3] = vmaxf(tmp, array[3]);
    return vminf(array[3], array[4]);
}

template<typename T>
inline T median(std::array<T, 9> array)
{
    T tmp = std::min(array[1], array[2]);
    array[2] = std::max(array[1], array[2]);
    array[1] = tmp;
    tmp = std::min(array[4], array[5]);
    array[5] = std::max(array[4], array[5]);
    array[4] = tmp;
    tmp = std::min(array[7], array[8]);
    array[8] = std::max(array[7], array[8]);
    array[7] = tmp;
    tmp = std::min(array[0], array[1]);
    array[1] = std::max(array[0], array[1]);
    array[0] = tmp;
    tmp = std::min(array[3], array[4]);
    array[4] = std::max(array[3], array[4]);
    array[3] = tmp;
    tmp = std::min(array[6], array[7]);
    array[7] = std::max(array[6], array[7]);
    array[6] = tmp;
    tmp = std::min(array[1], array[2]);
    array[2] = std::max(array[1], array[2]);
    array[1] = tmp;
    tmp = std::min(array[4], array[5]);
    array[5] = std::max(array[4], array[5]);
    array[4] = tmp;
    tmp = std::min(array[7], array[8]);
    array[8] = std::max(array[7], array[8]);
    array[3] = std::max(array[0], array[3]);
    array[5] = std::min(array[5], array[8]);
    array[7] = std::max(array[4], tmp);
    tmp = std::min(array[4], tmp);
    array[6] = std::max(array[3], array[6]);
    array[4] = std::max(array[1], tmp);
    array[2] = std::min(array[2], array[5]);
    array[4] = std::min(array[4], array[7]);
    tmp = std::min(array[4], array[2]);
    array[2] = std::max(array[4], array[2]);
    array[4] = std::max(array[6], tmp);
    return std::min(array[4], array[2]);
}

template<>
inline vfloat median(std::array<vfloat, 9> array)
{
    vfloat tmp = vminf(array[1], array[2]);
    array[2] = vmaxf(array[1], array[2]);
    array[1] = tmp;
    tmp = vminf(array[4], array[5]);
    array[5] = vmaxf(array[4], array[5]);
    array[4] = tmp;
    tmp = vminf(array[7], array[8]);
    array[8] = vmaxf(array[7], array[8]);
    array[7] = tmp;
    tmp = vminf(array[0], array[1]);
    array[1] = vmaxf(array[0], array[1]);
    array[0] = tmp;
    tmp = vminf(array[3], array[4]);
    array[4] = vmaxf(array[3], array[4]);
    array[3] = tmp;
    tmp = vminf(array[6], array[7]);
    array[7] = vmaxf(array[6], array[7]);
    array[6] = tmp;
    tmp = vminf(array[1], array[2]);
    array[2] = vmaxf(array[1], array[2]);
    array[1] = tmp;
    tmp = vminf(array[4], array[5]);
    array[5] = vmaxf(array[4], array[5]);
    array[4] = tmp;
    tmp = vminf(array[7], array[8]);
    array[8] = vmaxf(array[7], array[8]);
    array[3] = vmaxf(array[0], array[3]);
    array[5] = vminf(array[5], array[8]);
    array[7] = vmaxf(array[4], tmp);
    tmp = vminf(array[4], tmp);
    array[6] = vmaxf(array[3], array[6]);
    array[4] = vmaxf(array[1], tmp);
    array[2] = vminf(array[2], array[5]);
    array[4] = vminf(array[4], array[7]);
    tmp = vminf(array[4], array[2]);
    array[2] = vmaxf(array[4], array[2]);
    array[4] = vmaxf(array[6], tmp);
    return vminf(array[4], array[2]);
}

template<typename T>
inline T median(std::array<T, 25> array)
{
    T tmp = std::min(array[0], array[1]);
    array[1] = std::max(array[0], array[1]);
    array[0] = tmp;
    tmp = std::min(array[3], array[4]);
    array[4] = std::max(array[3], array[4]);
    array[3] = tmp;
    tmp = std::min(array[2], array[4]);
    array[4] = std::max(array[2], array[4]);
    array[2] = std::min(tmp, array[3]);
    array[3] = std::max(tmp, array[3]);
    tmp = std::min(array[6], array[7]);
    array[7] = std::max(array[6], array[7]);
    array[6] = tmp;
    tmp = std::min(array[5], array[7]);
    array[7] = std::max(array[5], array[7]);
    array[5] = std::min(tmp, array[6]);
    array[6] = std::max(tmp, array[6]);
    tmp = std::min(array[9], array[10]);
    array[10] = std::max(array[9], array[10]);
    array[9] = tmp;
    tmp = std::min(array[8], array[10]);
    array[10] = std::max(array[8], array[10]);
    array[8] = std::min(tmp, array[9]);
    array[9] = std::max(tmp, array[9]);
    tmp = std::min(array[12], array[13]);
    array[13] = std::max(array[12], array[13]);
    array[12] = tmp;
    tmp = std::min(array[11], array[13]);
    array[13] = std::max(array[11], array[13]);
    array[11] = std::min(tmp, array[12]);
    array[12] = std::max(tmp, array[12]);
    tmp = std::min(array[15], array[16]);
    array[16] = std::max(array[15], array[16]);
    array[15] = tmp;
    tmp = std::min(array[14], array[16]);
    array[16] = std::max(array[14], array[16]);
    array[14] = std::min(tmp, array[15]);
    array[15] = std::max(tmp, array[15]);
    tmp = std::min(array[18], array[19]);
    array[19] = std::max(array[18], array[19]);
    array[18] = tmp;
    tmp = std::min(array[17], array[19]);
    array[19] = std::max(array[17], array[19]);
    array[17] = std::min(tmp, array[18]);
    array[18] = std::max(tmp, array[18]);
    tmp = std::min(array[21], array[22]);
    array[22] = std::max(array[21], array[22]);
    array[21] = tmp;
    tmp = std::min(array[20], array[22]);
    array[22] = std::max(array[20], array[22]);
    array[20] = std::min(tmp, array[21]);
    array[21] = std::max(tmp, array[21]);
    tmp = std::min(array[23], array[24]);
    array[24] = std::max(array[23], array[24]);
    array[23] = tmp;
    tmp = std::min(array[2], array[5]);
    array[5] = std::max(array[2], array[5]);
    array[2] = tmp;
    tmp = std::min(array[3], array[6]);
    array[6] = std::max(array[3], array[6]);
    array[3] = tmp;
    tmp = std::min(array[0], array[6]);
    array[6] = std::max(array[0], array[6]);
    array[0] = std::min(tmp, array[3]);
    array[3] = std::max(tmp, array[3]);
    tmp = std::min(array[4], array[7]);
    array[7] = std::max(array[4], array[7]);
    array[4] = tmp;
    tmp = std::min(array[1], array[7]);
    array[7] = std::max(array[1], array[7]);
    array[1] = std::min(tmp, array[4]);
    array[4] = std::max(tmp, array[4]);
    tmp = std::min(array[11], array[14]);
    array[14] = std::max(array[11], array[14]);
    array[11] = tmp;
    tmp = std::min(array[8], array[14]);
    array[14] = std::max(array[8], array[14]);
    array[8] = std::min(tmp, array[11]);
    array[11] = std::max(tmp, array[11]);
    tmp = std::min(array[12], array[15]);
    array[15] = std::max(array[12], array[15]);
    array[12] = tmp;
    tmp = std::min(array[9], array[15]);
    array[15] = std::max(array[9], array[15]);
    array[9] = std::min(tmp, array[12]);
    array[12] = std::max(tmp, array[12]);
    tmp = std::min(array[13], array[16]);
    array[16] = std::max(array[13], array[16]);
    array[13] = tmp;
    tmp = std::min(array[10], array[16]);
    array[16] = std::max(array[10], array[16]);
    array[10] = std::min(tmp, array[13]);
    array[13] = std::max(tmp, array[13]);
    tmp = std::min(array[20], array[23]);
    array[23] = std::max(array[20], array[23]);
    array[20] = tmp;
    tmp = std::min(array[17], array[23]);
    array[23] = std::max(array[17], array[23]);
    array[17] = std::min(tmp, array[20]);
    array[20] = std::max(tmp, array[20]);
    tmp = std::min(array[21], array[24]);
    array[24] = std::max(array[21], array[24]);
    array[21] = tmp;
    tmp = std::min(array[18], array[24]);
    array[24] = std::max(array[18], array[24]);
    array[18] = std::min(tmp, array[21]);
    array[21] = std::max(tmp, array[21]);
    tmp = std::min(array[19], array[22]);
    array[22] = std::max(array[19], array[22]);
    array[19] = tmp;
    array[17] = std::max(array[8], array[17]);
    tmp = std::min(array[9], array[18]);
    array[18] = std::max(array[9], array[18]);
    array[9] = tmp;
    tmp = std::min(array[0], array[18]);
    array[18] = std::max(array[0], array[18]);
    array[9] = std::max(tmp, array[9]);
    tmp = std::min(array[10], array[19]);
    array[19] = std::max(array[10], array[19]);
    array[10] = tmp;
    tmp = std::min(array[1], array[19]);
    array[19] = std::max(array[1], array[19]);
    array[1] = std::min(tmp, array[10]);
    array[10] = std::max(tmp, array[10]);
    tmp = std::min(array[11], array[20]);
    array[20] = std::max(array[11], array[20]);
    array[11] = tmp;
    tmp = std::min(array[2], array[20]);
    array[20] = std::max(array[2], array[20]);
    array[11] = std::max(tmp, array[11]);
    tmp = std::min(array[12], array[21]);
    array[21] = std::max(array[12], array[21]);
    array[12] = tmp;
    tmp = std::min(array[3], array[21]);
    array[21] = std::max(array[3], array[21]);
    array[3] = std::min(tmp, array[12]);
    array[12] = std::max(tmp, array[12]);
    tmp = std::min(array[13], array[22]);
    array[22] = std::max(array[13], array[22]);
    array[4] = std::min(array[4], array[22]);
    array[13] = std::max(array[4], tmp);
    tmp = std::min(array[4], tmp);
    array[4] = tmp;
    tmp = std::min(array[14], array[23]);
    array[23] = std::max(array[14], array[23]);
    array[14] = tmp;
    tmp = std::min(array[5], array[23]);
    array[23] = std::max(array[5], array[23]);
    array[5] = std::min(tmp, array[14]);
    array[14] = std::max(tmp, array[14]);
    tmp = std::min(array[15], array[24]);
    array[24] = std::max(array[15], array[24]);
    array[15] = tmp;
    array[6] = std::min(array[6], array[24]);
    tmp = std::min(array[6], array[15]);
    array[15] = std::max(array[6], array[15]);
    array[6] = tmp;
    tmp = std::min(array[7], array[16]);
    array[7] = std::min(tmp, array[19]);
    tmp = std::min(array[13], array[21]);
    array[15] = std::min(array[15], array[23]);
    tmp = std::min(array[7], tmp);
    array[7] = std::min(tmp, array[15]);
    array[9] = std::max(array[1], array[9]);
    array[11] = std::max(array[3], array[11]);
    array[17] = std::max(array[5], array[17]);
    array[17] = std::max(array[11], array[17]);
    array[17] = std::max(array[9], array[17]);
    tmp = std::min(array[4], array[10]);
    array[10] = std::max(array[4], array[10]);
    array[4] = tmp;
    tmp = std::min(array[6], array[12]);
    array[12] = std::max(array[6], array[12]);
    array[6] = tmp;
    tmp = std::min(array[7], array[14]);
    array[14] = std::max(array[7], array[14]);
    array[7] = tmp;
    tmp = std::min(array[4], array[6]);
    array[6] = std::max(array[4], array[6]);
    array[7] = std::max(tmp, array[7]);
    tmp = std::min(array[12], array[14]);
    array[14] = std::max(array[12], array[14]);
    array[12] = tmp;
    array[10] = std::min(array[10], array[14]);
    tmp = std::min(array[6], array[7]);
    array[7] = std::max(array[6], array[7]);
    array[6] = tmp;
    tmp = std::min(array[10], array[12]);
    array[12] = std::max(array[10], array[12]);
    array[10] = std::max(array[6], tmp);
    tmp = std::min(array[6], tmp);
    array[17] = std::max(tmp, array[17]);
    tmp = std::min(array[12], array[17]);
    array[17] = std::max(array[12], array[17]);
    array[12] = tmp;
    array[7] = std::min(array[7], array[17]);
    tmp = std::min(array[7], array[10]);
    array[10] = std::max(array[7], array[10]);
    array[7] = tmp;
    tmp = std::min(array[12], array[18]);
    array[18] = std::max(array[12], array[18]);
    array[12] = std::max(array[7], tmp);
    array[10] = std::min(array[10], array[18]);
    tmp = std::min(array[12], array[20]);
    array[20] = std::max(array[12], array[20]);
    array[12] = tmp;
    tmp = std::min(array[10], array[20]);
    return std::max(tmp, array[12]);
}

template<>
inline vfloat median(std::array<vfloat, 25> array)
{
    vfloat tmp = vminf(array[0], array[1]);
    array[1] = vmaxf(array[0], array[1]);
    array[0] = tmp;
    tmp = vminf(array[3], array[4]);
    array[4] = vmaxf(array[3], array[4]);
    array[3] = tmp;
    tmp = vminf(array[2], array[4]);
    array[4] = vmaxf(array[2], array[4]);
    array[2] = vminf(tmp, array[3]);
    array[3] = vmaxf(tmp, array[3]);
    tmp = vminf(array[6], array[7]);
    array[7] = vmaxf(array[6], array[7]);
    array[6] = tmp;
    tmp = vminf(array[5], array[7]);
    array[7] = vmaxf(array[5], array[7]);
    array[5] = vminf(tmp, array[6]);
    array[6] = vmaxf(tmp, array[6]);
    tmp = vminf(array[9], array[10]);
    array[10] = vmaxf(array[9], array[10]);
    array[9] = tmp;
    tmp = vminf(array[8], array[10]);
    array[10] = vmaxf(array[8], array[10]);
    array[8] = vminf(tmp, array[9]);
    array[9] = vmaxf(tmp, array[9]);
    tmp = vminf(array[12], array[13]);
    array[13] = vmaxf(array[12], array[13]);
    array[12] = tmp;
    tmp = vminf(array[11], array[13]);
    array[13] = vmaxf(array[11], array[13]);
    array[11] = vminf(tmp, array[12]);
    array[12] = vmaxf(tmp, array[12]);
    tmp = vminf(array[15], array[16]);
    array[16] = vmaxf(array[15], array[16]);
    array[15] = tmp;
    tmp = vminf(array[14], array[16]);
    array[16] = vmaxf(array[14], array[16]);
    array[14] = vminf(tmp, array[15]);
    array[15] = vmaxf(tmp, array[15]);
    tmp = vminf(array[18], array[19]);
    array[19] = vmaxf(array[18], array[19]);
    array[18] = tmp;
    tmp = vminf(array[17], array[19]);
    array[19] = vmaxf(array[17], array[19]);
    array[17] = vminf(tmp, array[18]);
    array[18] = vmaxf(tmp, array[18]);
    tmp = vminf(array[21], array[22]);
    array[22] = vmaxf(array[21], array[22]);
    array[21] = tmp;
    tmp = vminf(array[20], array[22]);
    array[22] = vmaxf(array[20], array[22]);
    array[20] = vminf(tmp, array[21]);
    array[21] = vmaxf(tmp, array[21]);
    tmp = vminf(array[23], array[24]);
    array[24] = vmaxf(array[23], array[24]);
    array[23] = tmp;
    tmp = vminf(array[2], array[5]);
    array[5] = vmaxf(array[2], array[5]);
    array[2] = tmp;
    tmp = vminf(array[3], array[6]);
    array[6] = vmaxf(array[3], array[6]);
    array[3] = tmp;
    tmp = vminf(array[0], array[6]);
    array[6] = vmaxf(array[0], array[6]);
    array[0] = vminf(tmp, array[3]);
    array[3] = vmaxf(tmp, array[3]);
    tmp = vminf(array[4], array[7]);
    array[7] = vmaxf(array[4], array[7]);
    array[4] = tmp;
    tmp = vminf(array[1], array[7]);
    array[7] = vmaxf(array[1], array[7]);
    array[1] = vminf(tmp, array[4]);
    array[4] = vmaxf(tmp, array[4]);
    tmp = vminf(array[11], array[14]);
    array[14] = vmaxf(array[11], array[14]);
    array[11] = tmp;
    tmp = vminf(array[8], array[14]);
    array[14] = vmaxf(array[8], array[14]);
    array[8] = vminf(tmp, array[11]);
    array[11] = vmaxf(tmp, array[11]);
    tmp = vminf(array[12], array[15]);
    array[15] = vmaxf(array[12], array[15]);
    array[12] = tmp;
    tmp = vminf(array[9], array[15]);
    array[15] = vmaxf(array[9], array[15]);
    array[9] = vminf(tmp, array[12]);
    array[12] = vmaxf(tmp, array[12]);
    tmp = vminf(array[13], array[16]);
    array[16] = vmaxf(array[13], array[16]);
    array[13] = tmp;
    tmp = vminf(array[10], array[16]);
    array[16] = vmaxf(array[10], array[16]);
    array[10] = vminf(tmp, array[13]);
    array[13] = vmaxf(tmp, array[13]);
    tmp = vminf(array[20], array[23]);
    array[23] = vmaxf(array[20], array[23]);
    array[20] = tmp;
    tmp = vminf(array[17], array[23]);
    array[23] = vmaxf(array[17], array[23]);
    array[17] = vminf(tmp, array[20]);
    array[20] = vmaxf(tmp, array[20]);
    tmp = vminf(array[21], array[24]);
    array[24] = vmaxf(array[21], array[24]);
    array[21] = tmp;
    tmp = vminf(array[18], array[24]);
    array[24] = vmaxf(array[18], array[24]);
    array[18] = vminf(tmp, array[21]);
    array[21] = vmaxf(tmp, array[21]);
    tmp = vminf(array[19], array[22]);
    array[22] = vmaxf(array[19], array[22]);
    array[19] = tmp;
    array[17] = vmaxf(array[8], array[17]);
    tmp = vminf(array[9], array[18]);
    array[18] = vmaxf(array[9], array[18]);
    array[9] = tmp;
    tmp = vminf(array[0], array[18]);
    array[18] = vmaxf(array[0], array[18]);
    array[9] = vmaxf(tmp, array[9]);
    tmp = vminf(array[10], array[19]);
    array[19] = vmaxf(array[10], array[19]);
    array[10] = tmp;
    tmp = vminf(array[1], array[19]);
    array[19] = vmaxf(array[1], array[19]);
    array[1] = vminf(tmp, array[10]);
    array[10] = vmaxf(tmp, array[10]);
    tmp = vminf(array[11], array[20]);
    array[20] = vmaxf(array[11], array[20]);
    array[11] = tmp;
    tmp = vminf(array[2], array[20]);
    array[20] = vmaxf(array[2], array[20]);
    array[11] = vmaxf(tmp, array[11]);
    tmp = vminf(array[12], array[21]);
    array[21] = vmaxf(array[12], array[21]);
    array[12] = tmp;
    tmp = vminf(array[3], array[21]);
    array[21] = vmaxf(array[3], array[21]);
    array[3] = vminf(tmp, array[12]);
    array[12] = vmaxf(tmp, array[12]);
    tmp = vminf(array[13], array[22]);
    array[22] = vmaxf(array[13], array[22]);
    array[4] = vminf(array[4], array[22]);
    array[13] = vmaxf(array[4], tmp);
    tmp = vminf(array[4], tmp);
    array[4] = tmp;
    tmp = vminf(array[14], array[23]);
    array[23] = vmaxf(array[14], array[23]);
    array[14] = tmp;
    tmp = vminf(array[5], array[23]);
    array[23] = vmaxf(array[5], array[23]);
    array[5] = vminf(tmp, array[14]);
    array[14] = vmaxf(tmp, array[14]);
    tmp = vminf(array[15], array[24]);
    array[24] = vmaxf(array[15], array[24]);
    array[15] = tmp;
    array[6] = vminf(array[6], array[24]);
    tmp = vminf(array[6], array[15]);
    array[15] = vmaxf(array[6], array[15]);
    array[6] = tmp;
    tmp = vminf(array[7], array[16]);
    array[7] = vminf(tmp, array[19]);
    tmp = vminf(array[13], array[21]);
    array[15] = vminf(array[15], array[23]);
    tmp = vminf(array[7], tmp);
    array[7] = vminf(tmp, array[15]);
    array[9] = vmaxf(array[1], array[9]);
    array[11] = vmaxf(array[3], array[11]);
    array[17] = vmaxf(array[5], array[17]);
    array[17] = vmaxf(array[11], array[17]);
    array[17] = vmaxf(array[9], array[17]);
    tmp = vminf(array[4], array[10]);
    array[10] = vmaxf(array[4], array[10]);
    array[4] = tmp;
    tmp = vminf(array[6], array[12]);
    array[12] = vmaxf(array[6], array[12]);
    array[6] = tmp;
    tmp = vminf(array[7], array[14]);
    array[14] = vmaxf(array[7], array[14]);
    array[7] = tmp;
    tmp = vminf(array[4], array[6]);
    array[6] = vmaxf(array[4], array[6]);
    array[7] = vmaxf(tmp, array[7]);
    tmp = vminf(array[12], array[14]);
    array[14] = vmaxf(array[12], array[14]);
    array[12] = tmp;
    array[10] = vminf(array[10], array[14]);
    tmp = vminf(array[6], array[7]);
    array[7] = vmaxf(array[6], array[7]);
    array[6] = tmp;
    tmp = vminf(array[10], array[12]);
    array[12] = vmaxf(array[10], array[12]);
    array[10] = vmaxf(array[6], tmp);
    tmp = vminf(array[6], tmp);
    array[17] = vmaxf(tmp, array[17]);
    tmp = vminf(array[12], array[17]);
    array[17] = vmaxf(array[12], array[17]);
    array[12] = tmp;
    array[7] = vminf(array[7], array[17]);
    tmp = vminf(array[7], array[10]);
    array[10] = vmaxf(array[7], array[10]);
    array[7] = tmp;
    tmp = vminf(array[12], array[18]);
    array[18] = vmaxf(array[12], array[18]);
    array[12] = vmaxf(array[7], tmp);
    array[10] = vminf(array[10], array[18]);
    tmp = vminf(array[12], array[20]);
    array[20] = vmaxf(array[12], array[20]);
    array[12] = tmp;
    tmp = vminf(array[10], array[20]);
    return vmaxf(tmp, array[12]);
}

template<typename T, typename... ARGS>
inline T median(T arg, ARGS... args)
{
    return median(std::array<T, sizeof...(args) + 1>{std::move(arg), std::move(args)...});
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
