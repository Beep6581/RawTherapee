/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "improcfun.h"

namespace rtengine {

void ImProcFunctions::softLight(float *red, float *green, float *blue, int istart, int jstart, int tW, int tH, int TS)
{
    if (!params->softlight.enabled || !params->softlight.strength) {
        return;
    }

    const float blend = params->softlight.strength / 100.f;
    const float orig = 1.f - blend;

    const auto apply =
        [=](float x) -> float
        {
            if (!OOG(x)) {
                float v = Color::gamma_srgb(x) / MAXVALF;
                // Pegtop's formula from
                // https://en.wikipedia.org/wiki/Blend_modes#Soft_Light
                float v2 = v * v;
                float v22 = v2 * 2.f;
                v = v2 + v22 - v22 * v;
                x = blend * Color::igamma_srgb(v * MAXVALF) + orig * x;
            }
            return x;
        };

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        int ti = 0;
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int i = istart; i < tH; i++) {
            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                const int idx = ti * TS + tj;
                red[idx] = apply(red[idx]);
                green[idx] = apply(green[idx]);
                blue[idx] = apply(blue[idx]);
            }
            ++ti;
        }
    }
}

} // namespace rtengine
