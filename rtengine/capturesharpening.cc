/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Ingo Weyrich (heckflosse67@gmx.de)
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
#include <cmath>
#include <iostream>

#include "jaggedarray.h"
#include "rtengine.h"
#include "rawimagesource.h"
#include "rt_math.h"
#include "improcfun.h"
#include "procparams.h"
#include "color.h"
#include "gauss.h"
#include "rt_algo.h"
#define BENCHMARK
#include "StopWatch.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "opthelper.h"
#include "../rtgui/multilangmgr.h"

namespace {
void CaptureDeconvSharpening (float** luminance, float** tmp, const float * const * blend, int W, int H, double sigma, int iterations, rtengine::ProgressListener* plistener, double start, double step)
{

    rtengine::JaggedArray<float> tmpI(W, H);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int i = 0; i < H; i++) {
            for(int j = 0; j < W; j++) {
                tmpI[i][j] = max(luminance[i][j], 0.f);
            }
        }
        
        for (int k = 0; k < iterations; k++) {
            // apply gaussian blur and divide luminance by result of gaussian blur
            gaussianBlur(tmpI, tmp, W, H, sigma, nullptr, GAUSS_DIV, luminance);
            gaussianBlur(tmp, tmpI, W, H, sigma, nullptr, GAUSS_MULT);
            if (plistener) {
#ifdef _OPENMP
                #pragma omp single
#endif
                start += step;
                plistener->setProgress(start);
            }
        } // end for

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                luminance[i][j] = rtengine::intp(blend[i][j], max(tmpI[i][j], 0.0f), luminance[i][j]);
            }
        }
    } // end parallel
}

}

namespace rtengine
{

void RawImageSource::captureSharpening(const procparams::CaptureSharpeningParams &sharpeningParams, bool showMask, double &conrastThreshold) {
BENCHFUN


    if (plistener) {
        plistener->setProgressStr(M("TP_PDSHARPENING_LABEL"));
        plistener->setProgress(0.0);
    }

    const float xyz_rgb[3][3] = {          // XYZ from RGB
                                    { 0.412453, 0.357580, 0.180423 },
                                    { 0.212671, 0.715160, 0.072169 },
                                    { 0.019334, 0.119193, 0.950227 }
                                };

    float contrast = conrastThreshold / 100.f;

    array2D<float>& redVals = redCache ? *redCache : red;
    array2D<float>& greenVals = greenCache ? *greenCache : green;
    array2D<float>& blueVals = blueCache ? *blueCache : blue;

    if (showMask) {
        StopWatch Stop1("Show mask");
        array2D<float>& L = blue; // blue will be overridden anyway => we can use its buffer to store L
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; ++i) {
            Color::RGB2L(redVals[i], greenVals[i], blueVals[i], L[i], xyz_rgb, W);
        }
        if (plistener) {
            plistener->setProgress(0.1);
        }
        array2D<float>& blend = red; // red will be overridden anyway => we can use its buffer to store the blend mask
        buildBlendMask(L, blend, W, H, contrast, 1.f, sharpeningParams.autoContrast);
        if (plistener) {
            plistener->setProgress(0.2);
        }
        conrastThreshold = contrast * 100.f;
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                red[i][j] = green[i][j] = blue[i][j] = blend[i][j] * 16384.f;
            }
        }
        if (plistener) {
            plistener->setProgress(1.0);
        }
        return;
    }

    array2D<float>* Lbuffer = nullptr;
    if (!redCache) {
        Lbuffer = new array2D<float>(W, H);
    }

    array2D<float>* YOldbuffer = nullptr;
    if (!greenCache) {
        YOldbuffer = new array2D<float>(W, H);
    }

    array2D<float>* YNewbuffer = nullptr;
    if (!blueCache) {
        YNewbuffer = new array2D<float>(W, H);
    }
    array2D<float>& L = Lbuffer ? *Lbuffer : red;
    array2D<float>& YOld = YOldbuffer ? * YOldbuffer : green;
    array2D<float>& YNew = YNewbuffer ? * YNewbuffer : blue;

    StopWatch Stop1("rgb2YL");
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < H; ++i) {
        Color::RGB2L(redVals[i], greenVals[i], blueVals[i], L[i], xyz_rgb, W);
        Color::RGB2Y(redVals[i], greenVals[i], blueVals[i], YOld[i], YNew[i], sharpeningParams.gamma, W);
    }
    if (plistener) {
        plistener->setProgress(0.1);
    }
    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    buildBlendMask(L, blend, W, H, contrast, 1.f, sharpeningParams.autoContrast);
    if (plistener) {
        plistener->setProgress(0.2);
    }
    conrastThreshold = contrast * 100.f;

    Stop1.stop();
    array2D<float>& tmp = L; // L is not used anymore now => we can use its buffer as the needed temporary buffer
    CaptureDeconvSharpening(YNew, tmp, blend, W, H, sharpeningParams.deconvradius, sharpeningParams.deconviter, plistener, 0.2, (0.9 - 0.2) / sharpeningParams.deconviter);
    if (plistener) {
        plistener->setProgress(0.9);
    }
    StopWatch Stop2("Y2RGB");
    const float gamma = sharpeningParams.gamma;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < H; ++i) {
        int j = 0;
#ifdef __SSE2__
        const vfloat gammav = F2V(gamma);
        for (; j < W - 3; j += 4) {
            const vfloat factor = pow_F(vmaxf(LVFU(YNew[i][j]), ZEROV) / vmaxf(LVFU(YOld[i][j]), F2V(0.00001f)), gammav);
            STVFU(red[i][j], LVFU(redVals[i][j]) * factor);
            STVFU(green[i][j], LVFU(greenVals[i][j]) * factor);
            STVFU(blue[i][j], LVFU(blueVals[i][j]) * factor);
        }

#endif
        for (; j < W; ++j) {
            const float factor = pow_F(std::max(YNew[i][j], 0.f) / std::max(YOld[i][j], 0.00001f), gamma);
            red[i][j] = redVals[i][j] * factor;
            green[i][j] = greenVals[i][j] * factor;
            blue[i][j] = blueVals[i][j] * factor;
        }
    }
    Stop2.stop();
    delete Lbuffer;
    delete YOldbuffer;
    delete YNewbuffer;
    if (plistener) {
        plistener->setProgress(1.0);
    }
}

} /* namespace */
