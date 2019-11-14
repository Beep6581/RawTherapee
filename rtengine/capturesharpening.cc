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

#include "rtengine.h"
#include "rawimage.h"
#include "rawimagesource.h"
#include "rt_math.h"
#include "procparams.h"
#include "color.h"
#include "rt_algo.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "opthelper.h"
#include "../rtgui/multilangmgr.h"

namespace {

void compute7x7kernel(float sigma, float kernel[7][7]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -3; i <= 3; ++i) {
        for (int j = -3; j <= 3; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 1.15)) {
                kernel[i + 3][j + 3] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 3][j + 3];
            } else {
                kernel[i + 3][j + 3] = 0.f;
            }
        }
    }

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute5x5kernel(float sigma, float kernel[5][5]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 0.84)) {
                kernel[i + 2][j + 2] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 2][j + 2];
            } else {
                kernel[i + 2][j + 2] = 0.f;
            }
        }
    }

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute3x3kernel(float sigma, float kernel[3][3]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 0.84)) {
                kernel[i + 1][j + 1] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 1][j + 1];
            } else {
                kernel[i + 1][j + 1] = 0.f;
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

inline void initTile(float** dst, const int tileSize)
{

    // first rows
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < tileSize; ++j) {
            dst[i][j] = 1.f;
        }
    }

    // left and right border
    for (int i = 3; i < tileSize - 3; ++i) {
        dst[i][0] = dst[i][1] = dst[i][2] = 1.f;
        dst[i][tileSize - 3] = dst[i][tileSize - 2] = dst[i][tileSize - 1] = 1.f;
    }

    // last rows
    for (int i = tileSize - 3 ; i < tileSize; ++i) {
        for (int j = 0; j < tileSize; ++j) {
            dst[i][j] = 1.f;
        }
    }
}

inline void gauss3x3div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[3][3])
{

    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

inline void gauss5x5div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * (src[i - 2][j - 1] + src[i - 2][j + 1] + src[i - 1][j - 2] + src[i - 1][j + 2] + src[i + 1][j - 2] + src[i + 1][j + 2] + src[i + 2][j - 1] + src[i + 2][j + 1]) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

inline void gauss7x7div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[7][7])
{

    const float c31 = kernel[0][2];
    const float c30 = kernel[0][3];
    const float c22 = kernel[1][1];
    const float c21 = kernel[1][2];
    const float c20 = kernel[1][3];
    const float c11 = kernel[2][2];
    const float c10 = kernel[2][3];
    const float c00 = kernel[3][3];

    for (int i = 3; i < tileSize - 3; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * (src[i - 3][j - 1] + src[i - 3][j + 1] + src[i - 1][j - 3] + src[i - 1][j + 3] + src[i + 1][j - 3] + src[i + 1][j + 3] + src[i + 3][j - 1] + src[i + 3][j + 1]) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * (src[i - 2][j - 1] + src[i - 2][j + 1] * c21 + src[i - 1][j - 2] + src[i - 1][j + 2] + src[i + 1][j - 2] + src[i + 1][j + 2] + src[i + 2][j - 1] + src[i + 2][j + 1]) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

inline void gauss3x3mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[3][3])
{
    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] *= val;
        }
    }

}

inline void gauss5x5mult (float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * (src[i - 2][j - 1] + src[i - 2][j + 1] + src[i - 1][j - 2] + src[i - 1][j + 2] + src[i + 1][j - 2] + src[i + 1][j + 2] + src[i + 2][j - 1] + src[i + 2][j + 1]) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

inline void gauss7x7mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[7][7])
{

    const float c31 = kernel[0][2];
    const float c30 = kernel[0][3];
    const float c22 = kernel[1][1];
    const float c21 = kernel[1][2];
    const float c20 = kernel[1][3];
    const float c11 = kernel[2][2];
    const float c10 = kernel[2][3];
    const float c00 = kernel[3][3];

    for (int i = 3; i < tileSize - 3; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * (src[i - 3][j - 1] + src[i - 3][j + 1] + src[i - 1][j - 3] + src[i - 1][j + 3] + src[i + 1][j - 3] + src[i + 1][j + 3] + src[i + 3][j - 1] + src[i + 3][j + 1]) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * (src[i - 2][j - 1] + src[i - 2][j + 1] * c21 + src[i - 1][j - 2] + src[i - 1][j + 2] + src[i + 1][j - 2] + src[i + 1][j + 2] + src[i + 2][j - 1] + src[i + 2][j + 1]) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

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

bool checkForStop(float** tmpIThr, float** iterCheck, int fullTileSize, int border)
{
    bool stopped = false;
    for (int ii = border; !stopped && ii < fullTileSize - border; ++ii) {
#ifdef __SSE2__
        for (int jj = border; jj < fullTileSize - border; jj += 4) {
            if (_mm_movemask_ps((vfloat)vmaskf_lt(LVFU(tmpIThr[ii][jj]), LVFU(iterCheck[ii - border][jj - border])))) {
                stopped = true;
                break;
            }
        }
#else
        for (int jj = border; jj < fullTileSize - border; ++jj) {
            if (tmpIThr[ii][jj] < iterCheck[ii - border][jj - border]) {
                stopped = true;
                break;
            }
        }
#endif
    }
    return stopped;
}

void CaptureDeconvSharpening (float ** clipmask, float** luminance, float** oldLuminance, const float * const * blend, int W, int H, double sigma, double sigmaCornerOffset, int iterations, bool checkIterStop, rtengine::ProgressListener* plistener, double startVal, double endVal)
{
BENCHFUN
    const bool is5x5 = (sigma <= 0.84 && sigmaCornerOffset == 0.0);
    const bool is3x3 = (sigma < 0.6 && sigmaCornerOffset == 0.0);
    float kernel7[7][7];
    float kernel5[5][5];
    float kernel3[3][3];
    if (is3x3) {
        compute3x3kernel(sigma, kernel3);
    } else if (is5x5) {
        compute5x5kernel(sigma, kernel5);
    } else {
        compute7x7kernel(sigma, kernel7);
    }

    constexpr int tileSize = 32;
    constexpr int border = 5;
    constexpr int fullTileSize = tileSize + 2 * border;
    const float cornerRadius = std::min<float>(1.15f, sigma + sigmaCornerOffset);
    const float cornerDistance = sqrt(rtengine::SQR(W * 0.5f) + rtengine::SQR(H * 0.5f));
    const float distanceFactor = (cornerRadius - sigma) / cornerDistance;

    double progress = startVal;
    const double progressStep = (endVal - startVal) * rtengine::SQR(tileSize) / (W * H);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int progresscounter = 0;
        array2D<float> tmpIThr(fullTileSize, fullTileSize);
        array2D<float> tmpThr(fullTileSize, fullTileSize);
        array2D<float> lumThr(fullTileSize, fullTileSize);
        array2D<float> iterCheck(tileSize, tileSize);
        initTile(tmpThr, fullTileSize);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16) collapse(2)
#endif
        for (int i = border; i < H - border; i+= tileSize) {
            for(int j = border; j < W - border; j+= tileSize) {
                const bool endOfCol = (i + tileSize + border) >= H;
                const bool endOfRow = (j + tileSize + border) >= W;
                // fill tiles
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    if (checkIterStop) {
                        for (int k = 0, ii = endOfCol ? H - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? W - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                iterCheck[k][l] = oldLuminance[ii][jj] * clipmask[ii][jj] * 0.5f;
                            }
                        }
                    }
                    for (int k = 0, ii = endOfCol ? H - fullTileSize : i; k < fullTileSize; ++k, ++ii) {
                        for (int l = 0, jj = endOfRow ? W - fullTileSize : j; l < fullTileSize; ++l, ++jj) {
                            tmpIThr[k][l] = oldLuminance[ii - border][jj - border];
                            lumThr[k][l] = oldLuminance[ii - border][jj - border];
                        }
                    }
                } else {
                    if (checkIterStop) {
                        for (int ii = 0; ii < tileSize; ++ii) {
                            for (int jj = 0; jj < tileSize; ++jj) {
                                iterCheck[ii][jj] = oldLuminance[i + ii][j + jj] * clipmask[i + ii][j + jj] * 0.5f;
                            }
                        }
                    }
                    for (int ii = i; ii < i + fullTileSize; ++ii) {
                        for (int jj = j; jj < j + fullTileSize; ++jj) {
                            tmpIThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                            lumThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                        }
                    }
                }
                bool stopped = false;
                if (is3x3) {
                    for (int k = 0; k < iterations && !stopped; ++k) {
                        // apply 3x3 gaussian blur and divide luminance by result of gaussian blur
                        gauss3x3div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel3);
                        gauss3x3mult(tmpThr, tmpIThr, fullTileSize, kernel3);
                        if (checkIterStop) {
                            stopped = checkForStop(tmpIThr, iterCheck, fullTileSize, border);
                        }
                    }
                } else if (is5x5) {
                    for (int k = 0; k < iterations && !stopped; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss5x5div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel5);
                        gauss5x5mult(tmpThr, tmpIThr, fullTileSize, kernel5);
                        if (checkIterStop) {
                            stopped = checkForStop(tmpIThr, iterCheck, fullTileSize, border);
                        }
                    }
                } else {
                    if (sigmaCornerOffset != 0.0) {
                        const float distance = sqrt(rtengine::SQR(i + tileSize / 2 - H / 2) + rtengine::SQR(j + tileSize / 2 - W / 2));
                        const float sigmaTile = static_cast<float>(sigma) + distanceFactor * distance;
                        if (sigmaTile >= 0.4f) {
                            float lkernel7[7][7];
                            compute7x7kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel7);
                            for (int k = 0; k < iterations && !stopped; ++k) {
                                // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel7);
                                gauss7x7mult(tmpThr, tmpIThr, fullTileSize, lkernel7);
                                if (checkIterStop) {
                                    stopped = checkForStop(tmpIThr, iterCheck, fullTileSize, border);
                                }
                            }
                        }
                    } else {
                        for (int k = 0; k < iterations && !stopped; ++k) {
                            // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                            gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel7);
                            gauss7x7mult(tmpThr, tmpIThr, fullTileSize, kernel7);
                            if (checkIterStop) {
                                stopped = checkForStop(tmpIThr, iterCheck, fullTileSize, border);
                            }
                        }
                    }
                }
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    for (int k = border, ii = endOfCol ? H - fullTileSize - border : i - border; k < fullTileSize - border; ++k) {
                        for (int l = border, jj = endOfRow ? W - fullTileSize - border : j - border; l < fullTileSize - border; ++l) {
                            luminance[ii + k][jj + l] = rtengine::intp(blend[ii + k][jj + l], std::max(tmpIThr[k][l], 0.0f), luminance[ii + k][jj + l]);
                        }
                    }
                } else {
                    for (int ii = border; ii < fullTileSize - border; ++ii) {
                        for (int jj = border; jj < fullTileSize - border; ++jj) {
                            luminance[i + ii - border][j + jj - border] = rtengine::intp(blend[i + ii - border][j + jj - border], std::max(tmpIThr[ii][jj], 0.0f), luminance[i + ii - border][j + jj - border]);
                        }
                    }
                }
                if (plistener) {
                    if (++progresscounter % 16 == 0) {
#ifdef _OPENMP
                        #pragma omp critical(csprogress)
#endif
                        {
                            progress += 16.0 * progressStep;
                            progress = rtengine::min(progress, endVal);
                            plistener->setProgress(progress);
                        }
                    }
                }
            }
        }
    }
}

}

namespace rtengine
{

void RawImageSource::captureSharpening(const procparams::CaptureSharpeningParams &sharpeningParams, bool showMask, double &conrastThreshold, double &radius) {

    if (plistener) {
        plistener->setProgressStr(M("TP_PDSHARPENING_LABEL"));
        plistener->setProgress(0.0);
    }
BENCHFUN
    constexpr float xyz_rgb[3][3] = {          // XYZ from RGB
                                    { 0.412453, 0.357580, 0.180423 },
                                    { 0.212671, 0.715160, 0.072169 },
                                    { 0.019334, 0.119193, 0.950227 }
                                };

    float contrast = conrastThreshold / 100.f;

    const float clipVal = (ri->get_white(1) - ri->get_cblack(1)) * scale_mul[1];

    const array2D<float>& redVals = redCache ? *redCache : red;
    const array2D<float>& greenVals = greenCache ? *greenCache : green;
    const array2D<float>& blueVals = blueCache ? *blueCache : blue;

    array2D<float> clipMask(W, H);
    constexpr float clipLimit = 0.95f;
    constexpr float maxSigma = 1.15f;

    if (getSensorType() == ST_BAYER) {
        const float whites[2][2] = {
                                    {(ri->get_white(FC(0,0)) - c_black[FC(0,0)]) * scale_mul[FC(0,0)] * clipLimit, (ri->get_white(FC(0,1)) - c_black[FC(0,1)]) * scale_mul[FC(0,1)] * clipLimit},
                                    {(ri->get_white(FC(1,0)) - c_black[FC(1,0)]) * scale_mul[FC(1,0)] * clipLimit, (ri->get_white(FC(1,1)) - c_black[FC(1,1)]) * scale_mul[FC(1,1)] * clipLimit}
                                   };
        buildClipMaskBayer(rawData, W, H, clipMask, whites);
        const unsigned int fc[2] = {FC(0,0), FC(1,0)};
        if (sharpeningParams.autoRadius) {
            radius = std::min(calcRadiusBayer(rawData, W, H, 1000.f, clipVal, fc), maxSigma);
        }
    } else if (getSensorType() == ST_FUJI_XTRANS) {
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
            radius = std::min(calcRadiusXtrans(rawData, W, H, 1000.f, clipVal, i, j), maxSigma);
        }

    } else if (ri->get_colors() == 1) {
        buildClipMaskMono(rawData, W, H, clipMask, (ri->get_white(0) - c_black[0]) * scale_mul[0] * clipLimit);
        if (sharpeningParams.autoRadius) {
            const unsigned int fc[2] = {0, 0};
            radius = std::min(calcRadiusBayer(rawData, W, H, 1000.f, clipVal, fc), maxSigma);
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
        Color::RGB2Y(redVals[i], greenVals[i], blueVals[i], YOld[i], YNew[i], W);
    }
    if (plistener) {
        plistener->setProgress(0.1);
    }
    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    array2D<float>& blend = clipMask; // we can share blend and clipMask buffer here
    buildBlendMask(L, blend, W, H, contrast, 1.f, sharpeningParams.autoContrast, clipMask);
    if (plistener) {
        plistener->setProgress(0.2);
    }
    conrastThreshold = contrast * 100.f;

    CaptureDeconvSharpening(clipMask, YNew, YOld, blend, W, H, radius, sharpeningParams.deconvradiusOffset, sharpeningParams.deconviter, sharpeningParams.deconvitercheck, plistener, 0.2, 0.9);
    if (plistener) {
        plistener->setProgress(0.9);
    }
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < H; ++i) {
        int j = 0;
#ifdef __SSE2__
        for (; j < W - 3; j += 4) {
            const vfloat factor = vmaxf(LVFU(YNew[i][j]), ZEROV) / vmaxf(LVFU(YOld[i][j]), F2V(0.00001f));
            STVFU(red[i][j], LVFU(redVals[i][j]) * factor);
            STVFU(green[i][j], LVFU(greenVals[i][j]) * factor);
            STVFU(blue[i][j], LVFU(blueVals[i][j]) * factor);
        }

#endif
        for (; j < W; ++j) {
            const float factor = std::max(YNew[i][j], 0.f) / std::max(YOld[i][j], 0.00001f);
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
    rgbSourceModified = false;
}

} /* namespace */
