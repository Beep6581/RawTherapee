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

#include "coord.h"

namespace rtengine
{

void Coord::setFromPolar(PolarCoord polar)
{
    while (polar.angle <   0.f) {
        polar.angle += 360.f;
    }

    while (polar.angle > 360.f) {
        polar.angle -= 360.f;
    }

    x = polar.radius * cos(polar.angle / 180.f * M_PI);
    y = polar.radius * sin(polar.angle / 180.f * M_PI);
}

}
