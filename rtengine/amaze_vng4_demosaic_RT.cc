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
#include "sleef.c"
#include "opthelper.h"
#include "jaggedarray.h"
#include "gauss.h"
#include "StopWatch.h"

using namespace std;
namespace {

float calcBlendFactor(float val, float threshold) {
    // sigmoid function
    // result is in ]0;1] range
    // inflexion point is at (x, y) (threshold, 0.5)
    return 1.f / (1.f + xexpf(16.f - 16.f * val / threshold));
}

#ifdef __SSE2__
vfloat calcBlendFactor(vfloat valv, vfloat thresholdv) {
    // sigmoid function
    // result is in ]0;1] range
    // inflexion point is at (x, y) (threshold, 0.5)
    const vfloat onev = F2V(1.f);
    const vfloat c16v = F2V(16.f);
    return onev / (onev + xexpf(c16v - c16v * valv / thresholdv));
}
#endif

void buildBlendMask(float** luminance, rtengine::JaggedArray<float> &blend, int W, int H, float contrastThreshold, float amount = 1.f) {
BENCHFUN

    if(contrastThreshold == 0.f) {
        for(int j = 0; j < H; ++j) {
            for(int i = 0; i < W; ++i) {
                blend[j][i] = 1.f;
            }
        }
    } else {
        constexpr float scale = 0.0625f / 327.68f;
#ifdef _OPENMP
        #pragma omp parallel
#endif
    {
#ifdef __SSE2__
        const vfloat contrastThresholdv = F2V(contrastThreshold);
        const vfloat scalev = F2V(scale);
        const vfloat amountv = F2V(amount);
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for(int j = 2; j < H - 2; ++j) {
            int i = 2;
#ifdef __SSE2__
            for(; i < W - 5; i += 4) {
                vfloat contrastv = vsqrtf(SQRV(LVFU(luminance[j][i+1]) - LVFU(luminance[j][i-1])) + SQRV(LVFU(luminance[j+1][i]) - LVFU(luminance[j-1][i])) +
                                          SQRV(LVFU(luminance[j][i+2]) - LVFU(luminance[j][i-2])) + SQRV(LVFU(luminance[j+2][i]) - LVFU(luminance[j-2][i]))) * scalev;

                STVFU(blend[j][i], amountv * calcBlendFactor(contrastv, contrastThresholdv));
            }
#endif
            for(; i < W - 2; ++i) {

                float contrast = sqrtf(rtengine::SQR(luminance[j][i+1] - luminance[j][i-1]) + rtengine::SQR(luminance[j+1][i] - luminance[j-1][i]) + 
                                       rtengine::SQR(luminance[j][i+2] - luminance[j][i-2]) + rtengine::SQR(luminance[j+2][i] - luminance[j-2][i])) * scale;

                blend[j][i] = amount * calcBlendFactor(contrast, contrastThreshold);
            }
        }
#ifdef _OPENMP
        #pragma omp single
#endif
        {
            // upper border
            for(int j = 0; j < 2; ++j) {
                for(int i = 2; i < W - 2; ++i) {
                    blend[j][i] = blend[2][i];
                }
            }
            // lower border
            for(int j = H - 2; j < H; ++j) {
                for(int i = 2; i < W - 2; ++i) {
                    blend[j][i] = blend[H-3][i];
                }
            }
            for(int j = 0; j < H; ++j) {
                // left border
                blend[j][0] = blend[j][1] = blend[j][2];
                // right border
                blend[j][W - 2] = blend[j][W - 1] = blend[j][W - 3];
            }
        }
        // blur blend mask to smooth transitions
        gaussianBlur(blend, blend, W, H, 2.0);
    }
    }
}
}

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
                float a[winw] ALIGNED16;
                float b[winw] ALIGNED16;
                #pragma omp for
                for(int i = 0; i < winh; ++i) {
                    Color::RGB2Lab(redTmp[i], greenTmp[i], blueTmp[i], L[i], a, b, xyz_rgb, winw);
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
