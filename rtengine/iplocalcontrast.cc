/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Ported from G'MIC by Alberto Griggio <alberto.griggio@gmail.com>
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
#include "settings.h"
#include "gauss.h"
#include "boxblur.h"
#include "array2D.h"
#define BENCHMARK
#include "StopWatch.h"

namespace rtengine {

void ImProcFunctions::localContrast(LabImage *lab)
{
    BENCHFUN
    if (!params->localContrast.enabled) {
        return;
    }

    const int width = lab->W;
    const int height = lab->H;
    const float a = -params->localContrast.amount;
    const float half = 0.5f;
    const float dark = params->localContrast.darkness;
    const float light = params->localContrast.lightness;
    array2D<float> buf(width, height);
    const float sigma = params->localContrast.radius / scale;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    gaussianBlur(lab->L, buf, width, height, sigma);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float bufval = (buf[y][x] - lab->L[y][x]) * a;

            if (dark != 1 || light != 1) {
                bufval = max(bufval, half) * light + min(bufval, half) * dark;
            }

            lab->L[y][x] += bufval;
        }
    }
}

} // namespace rtengine
