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

void buildClipMaskBayer(const float * const *rawData, int W, int H, float** clipMask, const float whites[2][2])
{

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 0; row < H; ++row) {
        for (int col = 0; col < W; ++col) {
            clipMask[row][col] = 1.f;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 2; row < H - 2; ++row) {
        float clip0 = whites[row & 1][0];
        float clip1 = whites[row & 1][1];
        for (int col = 2; col < W - 2; ++col) {
            if (rawData[row][col] >= clip0) {
                clipMask[row - 2][col - 1] = clipMask[row - 2][col] = clipMask[row - 2][col + 1] = 0.f;
                clipMask[row - 1][col - 2] = clipMask[row - 1][col - 1] = clipMask[row - 1][col] = clipMask[row - 1][col + 1] = clipMask[row - 1][col + 2] = 0.f;
                clipMask[row][col - 2] = clipMask[row][col - 1] = clipMask[row][col] = clipMask[row][col + 1] = clipMask[row][col + 2] = 0.f;
                clipMask[row + 1][col - 2] = clipMask[row + 1][col - 1] = clipMask[row + 1][col] = clipMask[row + 1][col + 1] = clipMask[row + 1][col + 2] = 0.f;
                clipMask[row + 2][col - 1] = clipMask[row + 2][col] = clipMask[row + 2][col + 1] = 0.f;
            }
            std::swap(clip0, clip1);
        }
    }
}

void buildClipMaskXtrans(const float * const *rawData, int W, int H, float** clipMask, const float whites[6][6])
{

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 0; row < H; ++row) {
        for (int col = 0; col < W; ++col) {
            clipMask[row][col] = 1.f;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 2; row < H - 2; ++row) {
        for (int col = 2; col < W - 2; ++col) {
            const float clip = whites[row % 6][col % 6];
            if (rawData[row][col] >= clip) {
                clipMask[row - 2][col - 1] = clipMask[row - 2][col] = clipMask[row - 2][col + 1] = 0.f;
                clipMask[row - 1][col - 2] = clipMask[row - 1][col - 1] = clipMask[row - 1][col] = clipMask[row - 1][col + 1] = clipMask[row - 1][col + 2] = 0.f;
                clipMask[row][col - 2] = clipMask[row][col - 1] = clipMask[row][col] = clipMask[row][col + 1] = clipMask[row][col + 2] = 0.f;
                clipMask[row + 1][col - 2] = clipMask[row + 1][col - 1] = clipMask[row + 1][col] = clipMask[row + 1][col + 1] = clipMask[row + 1][col + 2] = 0.f;
                clipMask[row + 2][col - 1] = clipMask[row + 2][col] = clipMask[row + 2][col + 1] = 0.f;
            }
        }
    }
}

void buildClipMaskMono(const float * const *rawData, int W, int H, float** clipMask, float white)
{

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 0; row < H; ++row) {
        for (int col = 0; col < W; ++col) {
            clipMask[row][col] = 1.f;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 2; row < H - 2; ++row) {
        for (int col = 2; col < W - 2; ++col) {
            if (rawData[row][col] >= white) {
                clipMask[row - 2][col - 1] = clipMask[row - 2][col] = clipMask[row - 2][col + 1] = 0.f;
                clipMask[row - 1][col - 2] = clipMask[row - 1][col - 1] = clipMask[row - 1][col] = clipMask[row - 1][col + 1] = clipMask[row - 1][col + 2] = 0.f;
                clipMask[row][col - 2] = clipMask[row][col - 1] = clipMask[row][col] = clipMask[row][col + 1] = clipMask[row][col + 2] = 0.f;
                clipMask[row + 1][col - 2] = clipMask[row + 1][col - 1] = clipMask[row + 1][col] = clipMask[row + 1][col + 1] = clipMask[row + 1][col + 2] = 0.f;
                clipMask[row + 2][col - 1] = clipMask[row + 2][col] = clipMask[row + 2][col + 1] = 0.f;
            }
        }
    }
}

float calcRadiusBayer(const float * const *rawData, int W, int H, float lowerLimit, float upperLimit, const unsigned int fc[2])
{

    float maxRatio = 1.f;
#ifdef _OPENMP
    #pragma omp parallel for reduction(max:maxRatio) schedule(dynamic, 16)
#endif
    for (int row = 4; row < H - 4; ++row) {
        for (int col = 5 + (fc[row & 1] & 1); col < W - 4; col += 2) {
            const float val00 = rawData[row][col];
            if (val00 > 0.f) {
                const float val1m1 = rawData[row + 1][col - 1];
                const float val1p1 = rawData[row + 1][col + 1];
                const float maxVal0 = std::max(val00, val1m1);
                if (val1m1 > 0.f && maxVal0 > lowerLimit) {
                    const float minVal = std::min(val00, val1m1);
                    if (UNLIKELY(maxVal0 > maxRatio * minVal)) {
                        bool clipped = false;
                        if (maxVal0 == val00) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col - 1], rawData[row - 1][col + 1], val1p1) >= upperLimit) {
                                clipped = true;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row][col - 2], val00, rawData[row + 2][col - 2], rawData[row + 2][col]) >= upperLimit) {
                                clipped = true;
                            }
                        }
                        if (!clipped) {
                            maxRatio = maxVal0 / minVal;
                        }
                    }
                }
                const float maxVal1 = std::max(val00, val1p1);
                if (val1p1 > 0.f && maxVal1 > lowerLimit) {
                    const float minVal = std::min(val00, val1p1);
                    if (UNLIKELY(maxVal1 > maxRatio * minVal)) {
                        if (maxVal1 == val00) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col - 1], rawData[row - 1][col + 1], val1p1) >= upperLimit) {
                                continue;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(val00, rawData[row][col + 2], rawData[row + 2][col], rawData[row + 2][col + 2]) >= upperLimit) {
                                continue;
                            }
                        }
                        maxRatio = maxVal1 / minVal;
                    }
                }
            }
        }
    }
    return std::sqrt((1.f / (std::log(1.f / maxRatio) /  2.f)) / -2.f);
}

float calcRadiusXtrans(const float * const *rawData, int W, int H, float lowerLimit, float upperLimit, unsigned int starty, unsigned int startx)
{

    float maxRatio = 1.f;
#ifdef _OPENMP
    #pragma omp parallel for reduction(max:maxRatio) schedule(dynamic, 16)
#endif
    for (int row = starty + 3; row < H - 4; row += 3) {
        for (int col = startx + 3; col < W - 4; col += 3) {
            const float valtl = rawData[row][col];
            const float valtr = rawData[row][col + 1];
            const float valbl = rawData[row + 1][col];
            const float valbr = rawData[row + 1][col + 1];
            if (valtl > 1.f) {
                const float maxValtltr = std::max(valtl, valtr);
                if (valtr > 1.f && maxValtltr > lowerLimit) {
                    const float minVal = std::min(valtl, valtr);
                    if (UNLIKELY(maxValtltr > maxRatio * minVal)) {
                        bool clipped = false;
                        if (maxValtltr == valtl) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col - 1], valtr, valbl, valbr) >= upperLimit) {
                                clipped = true;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col + 2], valtl, valbl, valbr) >= upperLimit) {
                                clipped = true;
                            }
                        }
                        if (!clipped) {
                            maxRatio = maxValtltr / minVal;
                        }
                    }
                }
                const float maxValtlbl = std::max(valtl, valbl);
                if (valbl > 1.f && maxValtlbl > lowerLimit) {
                    const float minVal = std::min(valtl, valbl);
                    if (UNLIKELY(maxValtlbl > maxRatio * minVal)) {
                        bool clipped = false;
                        if (maxValtlbl == valtl) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col - 1], valtr, valbl, valbr) >= upperLimit) {
                                clipped = true;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(valtl, valtr, rawData[row + 2][col - 1], valbr) >= upperLimit) {
                                clipped = true;
                            }
                        }
                        if (!clipped) {
                            maxRatio = maxValtlbl / minVal;
                        }
                    }
                }
            }
            if (valbr > 1.f) {
                const float maxValblbr = std::max(valbl, valbr);
                if (valbl > 1.f && maxValblbr > lowerLimit) {
                    const float minVal = std::min(valbl, valbr);
                    if (UNLIKELY(maxValblbr > maxRatio * minVal)) {
                        bool clipped = false;
                        if (maxValblbr == valbr) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(valtl, valtr, valbl, rawData[row + 2][col + 2]) >= upperLimit) {
                                clipped = true;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(valtl, valtr, rawData[row + 2][col - 1], valbr) >= upperLimit) {
                                clipped = true;
                            }
                        }
                        if (!clipped) {
                            maxRatio = maxValblbr / minVal;
                        }
                    }
                }
                const float maxValtrbr = std::max(valtr, valbr);
                if (valtr > 1.f && maxValtrbr > lowerLimit) {
                    const float minVal = std::min(valtr, valbr);
                    if (UNLIKELY(maxValtrbr > maxRatio * minVal)) {
                        if (maxValtrbr == valbr) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(valtl, valtr, valbl, rawData[row + 2][col + 2]) >= upperLimit) {
                                continue;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col + 2], valtl, valbl, valbr) >= upperLimit) {
                                continue;
                            }
                        }
                        maxRatio = maxValtrbr / minVal;
                    }
                }
            }
        }
    }
    return std::sqrt((1.f / (std::log(1.f / maxRatio))) / -2.f);
}
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

void RawImageSource::captureSharpening(const procparams::CaptureSharpeningParams &sharpeningParams, bool showMask, double &conrastThreshold, double &radius) {

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

    const float clipVal = (ri->get_white(1) - ri->get_cblack(1)) * scale_mul[1];

    array2D<float>& redVals = redCache ? *redCache : red;
    array2D<float>& greenVals = greenCache ? *greenCache : green;
    array2D<float>& blueVals = blueCache ? *blueCache : blue;

    array2D<float> clipMask(W, H);
    constexpr float clipLimit = 0.95f;
    if (ri->getSensorType() == ST_BAYER) {
        const float whites[2][2] = {
                                    {(ri->get_white(FC(0,0)) - c_black[FC(0,0)]) * scale_mul[FC(0,0)] * clipLimit, (ri->get_white(FC(0,1)) - c_black[FC(0,1)]) * scale_mul[FC(0,1)] * clipLimit},
                                    {(ri->get_white(FC(1,0)) - c_black[FC(1,0)]) * scale_mul[FC(1,0)] * clipLimit, (ri->get_white(FC(1,1)) - c_black[FC(1,1)]) * scale_mul[FC(1,1)] * clipLimit}
                                   };
        buildClipMaskBayer(rawData, W, H, clipMask, whites);
        const unsigned int fc[2] = {FC(0,0), FC(1,0)};
        if (sharpeningParams.autoRadius) {
            radius = calcRadiusBayer(rawData, W, H, 1000.f, clipVal, fc);
        }
    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
        float whites[6][6];
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                const auto color = ri->XTRANSFC(i, j);
                whites[i][j] = (ri->get_white(color) - c_black[color]) * scale_mul[color] * clipLimit;
            }
        }
        buildClipMaskXtrans(rawData, W, H, clipMask, whites);
        bool found = false;
        int i, j;
        for (i = 6; i < 12 && !found; ++i) {
            for (j = 6; j < 12 && !found; ++j) {
                if (ri->XTRANSFC(i, j) == 1) {
                    if (ri->XTRANSFC(i, j - 1) != ri->XTRANSFC(i, j + 1)) {
                        if (ri->XTRANSFC(i - 1, j) != 1) {
                            if (ri->XTRANSFC(i, j - 1) != 1) {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (sharpeningParams.autoRadius) {
            radius = calcRadiusXtrans(rawData, W, H, 1000.f, clipVal, i, j);
        }

    } else if (ri->get_colors() == 1) {
        buildClipMaskMono(rawData, W, H, clipMask, (ri->get_white(0) - c_black[0]) * scale_mul[0] * clipLimit);
        if (sharpeningParams.autoRadius) {
            const unsigned int fc[2] = {0, 0};
            radius = calcRadiusBayer(rawData, W, H, 1000.f, clipVal, fc);
        }
    }

    if (showMask) {
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
        buildBlendMask(L, blend, W, H, contrast, 1.f, sharpeningParams.autoContrast, clipMask);
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
    buildBlendMask(L, blend, W, H, contrast, 1.f, sharpeningParams.autoContrast, clipMask);
    if (plistener) {
        plistener->setProgress(0.2);
    }
    conrastThreshold = contrast * 100.f;

    array2D<float>& tmp = L; // L is not used anymore now => we can use its buffer as the needed temporary buffer
    CaptureDeconvSharpening(YNew, tmp, blend, W, H, radius, sharpeningParams.deconviter, plistener, 0.2, (0.9 - 0.2) / sharpeningParams.deconviter);
    if (plistener) {
        plistener->setProgress(0.9);
    }
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

    delete Lbuffer;
    delete YOldbuffer;
    delete YNewbuffer;
    if (plistener) {
        plistener->setProgress(1.0);
    }
}

} /* namespace */
