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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "coord.h"

#include "rt_math.h"

namespace rtengine
{

Coord& Coord::operator= (const PolarCoord& other)
{
    const auto radius = other.radius;
    const auto angle = other.angle / 180.0 * rtengine::RT_PI;

    x = radius * std::cos (angle);
    y = radius * std::sin (angle);

    return *this;
}

PolarCoord& PolarCoord::operator= (const Coord& other)
{
    const double x = other.x;
    const double y = other.y;

    radius = rtengine::norm2 (x, y);
    angle = std::atan2 (y, x) * 180.0 / rtengine::RT_PI;

    return *this;
}

/// @brief Clip the coord to stay in the width x height bounds
/// @return true if the x or y coordinate has changed
bool Coord::clip (const int width, const int height)
{
    const auto newX = rtengine::LIM<int> (x, 0, width);
    const auto newY = rtengine::LIM<int> (y, 0, height);

    if (x != newX || y != newY) {

        x = newX;
        y = newY;

        return true;
    } else {
        return false;
    }
}

}
