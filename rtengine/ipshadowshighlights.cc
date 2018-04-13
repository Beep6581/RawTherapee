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
#include "gauss.h"

namespace rtengine {

void ImProcFunctions::shadowsHighlights(LabImage *lab)
{
    if (!params->sh.enabled || (!params->sh.highlights && !params->sh.shadows)){
        return;
    }

    const int width = lab->W;
    const int height = lab->H;

    array2D<float> mask(width, height);
    const float sigma = params->sh.radius * 5.f / scale;

    const auto apply =
        [&](int amount, int tonalwidth, bool hl) -> void
        {
            // first highlights
            const float thresh = tonalwidth * 327.68f;
            const float scale = hl ? (thresh > 0.f ? 0.9f / thresh : 1.f) : thresh * 0.9f;
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float l = lab->L[y][x];
                    if (hl) {
                        mask[y][x] = (l > thresh) ? 1.f : std::pow(l * scale, 4);
                    } else {
                        mask[y][x] = l <= thresh ? 1.f : std::pow(scale / l, 4);
                    }
                }
            }
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                gaussianBlur(mask, mask, width, height, sigma);
            }

            const float base = std::pow(4.f, float(amount)/100.f);
            const float gamma = hl ? base : 1.f / base;

            LUTf f(32768);
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int l = 0; l < 32768; ++l) {
                f[l] = std::pow(l / 32768.f, gamma) * 32768.f;
            }
            
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float l = lab->L[y][x];
                    float blend = mask[y][x];
                    float orig = 1.f - blend;
                    if (l >= 0.f && l < 32768.f) {
                        lab->L[y][x] = f[l] * blend + l * orig;
                        if (!hl && l > 1.f) {
                            // when pushing shadows, scale also the chromaticity
                            float s = max(lab->L[y][x] / l * 0.5f, 1.f) * blend;
                            float a = lab->a[y][x];
                            float b = lab->b[y][x];
                            lab->a[y][x] = a * s + a * orig;
                            lab->b[y][x] = b * s + b * orig;
                        }
                    }
                }
            }
        };

    if (params->sh.highlights > 0) {
        apply(params->sh.highlights, params->sh.htonalwidth, true);
    }

    if (params->sh.shadows > 0) {
        apply(params->sh.shadows, params->sh.stonalwidth, false);
    }
}

} // namespace rtengine
