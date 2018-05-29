////////////////////////////////////////////////////////////////
//
//          AMaZE demosaic algorithm
// (Aliasing Minimization and Zipper Elimination)
//
//  copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//  optimized for speed by Ingo Weyrich
//
// incorporating ideas of Luis Sanz Rodrigues and Paul Lee
//
// code dated: May 27, 2010
// latest modification: Ingo Weyrich, January 25, 2016
//
//  amaze_interpolate_RT.cc is free software: you can redistribute it and/or modify
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
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include "rtengine.h"
#include "rawimagesource.h"
#include "rt_math.h"
#define BENCHMARK
#include "StopWatch.h"
#include "rt_algo.h"

using namespace std;

namespace rtengine
{

void RawImageSource::amaze_vng4_demosaic_RT(int winw, int winh, array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue)
{
    BENCHFUN

            vng4_demosaic ();
            array2D<float> redTmp(winw, winh);
            array2D<float> greenTmp(winw, winh);
            array2D<float> blueTmp(winw, winh);
            array2D<float> L(winw, winh);
            amaze_demosaic_RT (0, 0, winw, winh, rawData, redTmp, greenTmp, blueTmp);
            const float xyz_rgb[3][3] = {          // XYZ from RGB
                                        { 0.412453, 0.357580, 0.180423 },
                                        { 0.212671, 0.715160, 0.072169 },
                                        { 0.019334, 0.119193, 0.950227 }
                                        };
            #pragma omp parallel
            {
                #pragma omp for
                for(int i = 0; i < winh; ++i) {
                    Color::RGB2L(redTmp[i], greenTmp[i], blueTmp[i], L[i], xyz_rgb, winw);
                }
            }
            // calculate contrast based blend factors to reduce sharpening in regions with low contrast
            JaggedArray<float> blend(winw, winh);
            buildBlendMask(L, blend, winw, winh, 20.f / 100.f);
            #pragma omp parallel for
            for(int i = 0; i < winh; ++i) {
                for(int j = 0; j < winw; ++j) {
                    red[i][j] = intp(blend[i][j], redTmp[i][j], red[i][j]);
                    green[i][j] = intp(blend[i][j], greenTmp[i][j], green[i][j]);
                    blue[i][j] = intp(blend[i][j], blueTmp[i][j], blue[i][j]);
                }
            }

}
}
