/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2021 Ingo Weyrich <heckflosse67@gmx.de>
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cstdint>
#include <vector>

struct GainMap
{
    uint32_t Top;
    uint32_t Left;
    uint32_t Bottom;
    uint32_t Right;
    uint32_t Plane;
    uint32_t Planes;
    uint32_t RowPitch;
    uint32_t ColPitch;
    uint32_t MapPointsV;
    uint32_t MapPointsH;
    double MapSpacingV;
    double MapSpacingH;
    double MapOriginV;
    double MapOriginH;
    uint32_t MapPlanes;
    std::vector<float> MapGain;
};
