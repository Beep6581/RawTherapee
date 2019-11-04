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
#pragma once

#include "color.h"
#include "curves.h"

namespace rtengine
{

class SatAndValueBlendingToneCurve : public ToneCurve
{
public:
    void Apply(float& r, float& g, float& b) const;
};

// Tone curve modifying the value channel only, preserving hue and saturation
// values in 0xffff space
inline void SatAndValueBlendingToneCurve::Apply (float& ir, float& ig, float& ib) const
{

    assert (lutToneCurve);

    float r = CLIP(ir);
    float g = CLIP(ig);
    float b = CLIP(ib);

    const float lum = (r + g + b) / 3.f;
    const float newLum = lutToneCurve[lum];

    if (newLum == lum) {
        return;
    }

    float h, s, v;
    Color::rgb2hsvtc(r, g, b, h, s, v);

    float dV;
    if (newLum > lum) {
        // Linearly targeting Value = 1 and Saturation = 0
        const float coef = (newLum - lum) / (65535.f - lum);
        dV = (1.f - v) * coef;
        s *= 1.f - coef;
    } else {
        // Linearly targeting Value = 0
        const float coef = (newLum - lum) / lum ;
        dV = v * coef;
    }
    Color::hsv2rgbdcp(h, s, v + dV, r, g, b);

    setUnlessOOG(ir, ig, ib, r, g, b);
}

}
