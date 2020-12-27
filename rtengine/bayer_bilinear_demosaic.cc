////////////////////////////////////////////////////////////////
//
//  Bilinear bayer demosaic, optimized for speed, intended use is for flat regions of dual-demosaic
//
//  copyright (c) 2020  Ingo Weyrich <heckflosse67@gmx.de>
//
//
//  code dated: May 09, 2020
//
//  bayer_bilinear_demosaic.cc is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include "rawimagesource.h"
#include "rt_math.h"

using namespace rtengine;

void RawImageSource::bayer_bilinear_demosaic(const float* const * blend, const array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue)
{

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 1; i < H - 1; ++i) {
        float **nonGreen1 = red;
        float **nonGreen2 = blue;
        if (FC(i, 0) == 2 || FC(i, 1) == 2) { // blue row => swap pointers
            std::swap(nonGreen1, nonGreen2);
        }
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 2 - (FC(i, 1) & 1); j < W - 2; j += 2) { // always begin with a green pixel
            green[i][j] = intp(blend[i][j], green[i][j], rawData[i][j]);
            nonGreen1[i][j] = intp(blend[i][j], nonGreen1[i][j], (rawData[i][j - 1] + rawData[i][j + 1]) * 0.5f);
            nonGreen2[i][j] = intp(blend[i][j], nonGreen2[i][j], (rawData[i - 1][j] + rawData[i + 1][j]) * 0.5f);
            green[i][j + 1] = intp(blend[i][j + 1], green[i][j + 1], ((rawData[i - 1][j + 1] + rawData[i][j]) + (rawData[i][j + 2] + rawData[i + 1][j + 1])) * 0.25f);
            nonGreen1[i][j + 1] = intp(blend[i][j + 1], nonGreen1[i][j + 1], rawData[i][j + 1]);
            nonGreen2[i][j + 1] = intp(blend[i][j + 1], nonGreen2[i][j + 1], ((rawData[i - 1][j] + rawData[i - 1][j + 2]) + (rawData[i + 1][j] + rawData[i + 1][j + 2])) * 0.25f);
        }
    }
}
