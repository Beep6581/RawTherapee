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
#include "rtgui/multilangmgr.h"
#include "improcfun.h"
#include "boxblur.h"
#include "median.h"

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
    return std::sqrt(1.f / std::log(maxRatio));
}

float calcRadiusXtrans(const float * const *rawData, int W, int H, float lowerLimit, float upperLimit, unsigned int starty, unsigned int startx)
{

    float maxRatio = 1.f;
#ifdef _OPENMP
    #pragma omp parallel for reduction(max:maxRatio) schedule(dynamic, 16)
#endif
    for (int row = starty + 2; row < H - 4; row += 3) {
        for (int col = startx + 2; col < W - 4; col += 3) {
            const float valp1p1 = rawData[row + 1][col + 1];
            const bool squareClipped = rtengine::max(valp1p1, rawData[row + 1][col + 2], rawData[row + 2][col + 1], rawData[row + 2][col + 2]) >= upperLimit;
            const float greenSolitary = rawData[row][col];
            if (greenSolitary > 1.f && std::max(rawData[row - 1][col - 1], rawData[row - 1][col + 1]) < upperLimit) {
                if (greenSolitary < upperLimit) {
                    const float valp1m1 = rawData[row + 1][col - 1];
                    if (valp1m1 > 1.f && rtengine::max(rawData[row + 1][col - 2], valp1m1, rawData[row + 2][col - 2], rawData[row + 1][col - 1]) < upperLimit) {
                        const float maxVal = std::max(greenSolitary, valp1m1);
                        if (maxVal > lowerLimit) {
                            const float minVal = std::min(greenSolitary, valp1m1);
                            if (UNLIKELY(maxVal > maxRatio * minVal)) {
                                maxRatio = maxVal / minVal;
                            }
                        }
                    }
                    if (valp1p1 > 1.f && !squareClipped) {
                        const float maxVal = std::max(greenSolitary, valp1p1);
                        if (maxVal > lowerLimit) {
                            const float minVal = std::min(greenSolitary, valp1p1);
                            if (UNLIKELY(maxVal > maxRatio * minVal)) {
                                maxRatio = maxVal / minVal;
                            }
                        }
                    }
                }
            }
            if (!squareClipped) {
                const float valp2p2 = rawData[row + 2][col + 2];
                if (valp2p2 > 1.f) {
                    if (valp1p1 > 1.f) {
                        const float maxVal = std::max(valp1p1, valp2p2);
                        if (maxVal > lowerLimit) {
                            const float minVal = std::min(valp1p1, valp2p2);
                            if (UNLIKELY(maxVal > maxRatio * minVal)) {
                                maxRatio = maxVal / minVal;
                            }
                        }
                    }
                    const float greenSolitaryRight = rawData[row + 3][col + 3];
                    if (rtengine::max(greenSolitaryRight, rawData[row + 4][col + 2], rawData[row + 4][col + 4]) < upperLimit) {
                        if (greenSolitaryRight > 1.f) {
                            const float maxVal = std::max(greenSolitaryRight, valp2p2);
                            if (maxVal > lowerLimit) {
                                const float minVal = std::min(greenSolitaryRight, valp2p2);
                                if (UNLIKELY(maxVal > maxRatio * minVal)) {
                                    maxRatio = maxVal / minVal;
                                }
                            }
                        }
                    }
                }
                const float valp1p2 = rawData[row + 1][col + 2];
                const float valp2p1 = rawData[row + 2][col + 1];
                if (valp2p1 > 1.f) {
                    if (valp1p2 > 1.f) {
                        const float maxVal = std::max(valp1p2, valp2p1);
                        if (maxVal > lowerLimit) {
                            const float minVal = std::min(valp1p2, valp2p1);
                            if (UNLIKELY(maxVal > maxRatio * minVal)) {
                                maxRatio = maxVal / minVal;
                            }
                        }
                    }
                    const float greenSolitaryLeft = rawData[row + 3][col];
                    if (rtengine::max(greenSolitaryLeft, rawData[row + 4][col - 1], rawData[row + 4][col + 1]) < upperLimit) {
                        if (greenSolitaryLeft > 1.f) {
                            const float maxVal = std::max(greenSolitaryLeft, valp2p1);
                            if (maxVal > lowerLimit) {
                                const float minVal = std::min(greenSolitaryLeft, valp2p1);
                                if (UNLIKELY(maxVal > maxRatio * minVal)) {
                                    maxRatio = maxVal / minVal;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return std::sqrt(1.f / std::log(maxRatio));
}

bool checkForStop(float** tmpIThr, float** iterCheck, int fullTileSize, int border)
{
    for (int ii = border; ii < fullTileSize - border; ++ii) {
#ifdef __SSE2__
        for (int jj = border; jj < fullTileSize - border; jj += 4) {
            if (UNLIKELY(_mm_movemask_ps((vfloat)vmaskf_lt(LVFU(tmpIThr[ii][jj]), LVFU(iterCheck[ii - border][jj - border]))))) {
                return true;
            }
        }
#else
        for (int jj = border; jj < fullTileSize - border; ++jj) {
            if (tmpIThr[ii][jj] < iterCheck[ii - border][jj - border]) {
                return true;
            }
        }
#endif
    }
    return false;
}

void CaptureDeconvSharpening (float** luminance, const float* const * oldLuminance, const float * const * blend, int W, int H, float sigma, float sigmaCornerOffset, int iterations, bool checkIterStop, rtengine::ProgressListener* plistener, double startVal, double endVal)
{
BENCHFUN
    const bool is9x9 = (sigma <= 1.5f && sigmaCornerOffset == 0.f);
    const bool is7x7 = (sigma <= 1.15f && sigmaCornerOffset == 0.f);
    const bool is5x5 = (sigma <= 0.84f && sigmaCornerOffset == 0.f);
    const bool is3x3 = (sigma < 0.6f && sigmaCornerOffset == 0.f);
    float kernel13[13][13];
    float kernel9[9][9];
    float kernel7[7][7];
    float kernel5[5][5];
    float kernel3[3][3];
    if (is3x3) {
        rtengine::compute3x3kernel2(sigma, kernel3);
    } else if (is5x5) {
        rtengine::compute5x5kernel2(sigma, kernel5);
    } else if (is7x7) {
        rtengine::compute7x7kernel2(sigma, kernel7);
    } else if (is9x9) {
        rtengine::compute9x9kernel2(sigma, kernel9);
    } else {
        rtengine::compute13x13kernel2(sigma, kernel13);
    }

    constexpr int tileSize = 32;
    const int border = (is3x3 || is5x5 || is7x7) ? iterations <= 30 ? 5 : 7 : 8;
    const int fullTileSize = tileSize + 2 * border;
    const float cornerRadius = std::min<float>(2.f, sigma + sigmaCornerOffset);
    const float cornerDistance = sqrt(rtengine::SQR(W * 0.5f) + rtengine::SQR(H * 0.5f));
    const float distanceFactor = (cornerRadius - sigma) / cornerDistance;

    double progress = startVal;
    const double progressStep = (endVal - startVal) * rtengine::SQR(tileSize) / (W * H);

    constexpr float minBlend = 0.01f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int progresscounter = 0;
        array2D<float> tmpIThr(fullTileSize, fullTileSize);
        array2D<float> tmpThr(fullTileSize, fullTileSize);
        tmpThr.fill(1.f);
        array2D<float> lumThr(fullTileSize, fullTileSize);
        array2D<float> iterCheck(tileSize, tileSize);
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
                    float maxVal = 0.f;
                    if (checkIterStop) {
                        for (int k = 0, ii = endOfCol ? H - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? W - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                iterCheck[k][l] = oldLuminance[ii][jj] * blend[ii][jj] * 0.5f;
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    } else {
                        for (int k = 0, ii = endOfCol ? H - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? W - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    }
                    if (maxVal < minBlend) {
                        // no pixel of the tile has a blend factor >= minBlend => skip the tile
                        continue;
                    }
                    for (int k = 0, ii = endOfCol ? H - fullTileSize : i - border; k < fullTileSize; ++k, ++ii) {
                        for (int l = 0, jj = endOfRow ? W - fullTileSize : j - border; l < fullTileSize; ++l, ++jj) {
                            tmpIThr[k][l] = oldLuminance[ii][jj];
                            lumThr[k][l] = oldLuminance[ii][jj];
                        }
                    }
                } else {
                    float maxVal = 0.f;
                    if (checkIterStop) {
                        for (int ii = 0; ii < tileSize; ++ii) {
                            for (int jj = 0; jj < tileSize; ++jj) {
                                iterCheck[ii][jj] = oldLuminance[i + ii][j + jj] * blend[i + ii][j + jj] * 0.5f;
                                maxVal = std::max(maxVal, blend[i + ii][j + jj]);
                            }
                        }
                    } else {
                        for (int ii = 0; ii < tileSize; ++ii) {
                            for (int jj = 0; jj < tileSize; ++jj) {
                                maxVal = std::max(maxVal, blend[i + ii][j + jj]);
                            }
                        }
                    }
                    if (maxVal < minBlend) {
                        // no pixel of the tile has a blend factor >= minBlend => skip the tile
                        continue;
                    }
                    for (int ii = i; ii < i + fullTileSize; ++ii) {
                        for (int jj = j; jj < j + fullTileSize; ++jj) {
                            tmpIThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                            lumThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                        }
                    }
                }
                if (is3x3) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 3x3 gaussian blur and divide luminance by result of gaussian blur
                        rtengine::gauss3x3div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel3);
                        rtengine::gauss3x3mult2(tmpThr, tmpIThr, fullTileSize, kernel3);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is5x5) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        rtengine::gauss5x5div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel5);
                        rtengine::gauss5x5mult2(tmpThr, tmpIThr, fullTileSize, kernel5);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is7x7) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        rtengine::gauss7x7div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel7);
                        rtengine::gauss7x7mult2(tmpThr, tmpIThr, fullTileSize, kernel7);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is9x9) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        rtengine::gauss9x9div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel9);
                        rtengine::gauss9x9mult2(tmpThr, tmpIThr, fullTileSize, kernel9);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else {
                    if (sigmaCornerOffset != 0.f) {
                        const float distance = sqrt(rtengine::SQR(i + tileSize / 2 - H / 2) + rtengine::SQR(j + tileSize / 2 - W / 2));
                        const float sigmaTile = static_cast<float>(sigma) + distanceFactor * distance;
                        if (sigmaTile >= 0.4f) {
                            if (sigmaTile > 1.5f) { // have to use 13x13 kernel
                                float lkernel13[13][13];
                                rtengine::compute13x13kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel13);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss13x13div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel13);
                                    rtengine::gauss13x13mult2(tmpThr, tmpIThr, fullTileSize, lkernel13);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 1.15f) { // have to use 9x9 kernel
                                float lkernel9[9][9];
                                rtengine::compute9x9kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel9);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 9x9 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss9x9div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel9);
                                    rtengine::gauss9x9mult2(tmpThr, tmpIThr, fullTileSize, lkernel9);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 0.84f) { // have to use 7x7 kernel
                                float lkernel7[7][7];
                                rtengine::compute7x7kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel7);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss7x7div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel7);
                                    rtengine::gauss7x7mult2(tmpThr, tmpIThr, fullTileSize, lkernel7);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else { // can use 5x5 kernel
                                float lkernel5[5][5];
                                rtengine::compute5x5kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel5);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss5x5div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel5);
                                    rtengine::gauss5x5mult2(tmpThr, tmpIThr, fullTileSize, lkernel5);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        for (int k = 0; k < iterations; ++k) {
                            // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                            rtengine::gauss13x13div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel13);
                            rtengine::gauss13x13mult2(tmpThr, tmpIThr, fullTileSize, kernel13);
                            if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                break;
                            }
                        }
                    }
                }
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    for (int k = border, ii = endOfCol ? H - fullTileSize : i - border; k < fullTileSize - border; ++k) {
                        for (int l = border, jj = endOfRow ? W - fullTileSize : j - border; l < fullTileSize - border; ++l) {
                            luminance[ii + k][jj + l] = rtengine::intp(blend[ii + k][jj + l], tmpIThr[k][l], luminance[ii + k][jj + l]);
                        }
                    }
                } else {
                    for (int ii = border; ii < fullTileSize - border; ++ii) {
                        for (int jj = border; jj < fullTileSize - border; ++jj) {
                            luminance[i + ii - border][j + jj - border] = rtengine::intp(blend[i + ii - border][j + jj - border], tmpIThr[ii][jj], luminance[i + ii - border][j + jj - border]);
                        }
                    }
                }
                if (plistener) {
                    if (++progresscounter % 32 == 0) {
#ifdef _OPENMP
                        #pragma omp critical(csprogress)
#endif
                        {
                            progress += 32.0 * progressStep;
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

    if (!(ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS || ri->get_colors() == 1)) {
        return;
    }

    if (plistener) {
        plistener->setProgressStr(M("TP_PDSHARPENING_LABEL88"));
        plistener->setProgress(0.0);
    }
BENCHFUN

    constexpr float xyz_rgb[3][3] = {          // XYZ from RGB
                                    { 0.412453, 0.357580, 0.180423 },
                                    { 0.212671, 0.715160, 0.072169 },
                                    { 0.019334, 0.119193, 0.950227 }
                                };

    float contrast = conrastThreshold / 100.0;

    const float clipVal = (ri->get_white(1) - ri->get_cblack(1)) * scale_mul[1];
    array2D<float> redVals (W, H);
    array2D<float> greenVals(W, H);
    array2D<float> blueVals(W, H);
     
    redVals = redCache ? *redCache : red;
    greenVals = greenCache ? *greenCache : green;
    blueVals = blueCache ? *blueCache : blue;

    typedef ImProcFunctions::Median Median;

    //predoise : small median to denoise before capture sharpening : allow CS to work correctly and reduce a little the noise
    //J.Desmis October 2024
    if(sharpeningParams.noisecap > 0.f) {
        //I have choose median due to its low aggressiveness and for a 3x3 its speed
        float denstr = 0.01 * sharpeningParams.noisecap;

        float** tmL;
        float** mR;
        float** mG;
        float** mB;
        int wid = W;
        int hei = H;
        tmL = new float*[hei];
        mR = new float*[hei];
        mG = new float*[hei];
        mB = new float*[hei];

        for (int i = 0; i < hei; ++i) {
            tmL[i] = new float[wid];
            mR[i] = new float[wid];
            mG[i] = new float[wid];
            mB[i] = new float[wid];
        }
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif                  
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                mR[i][j] = redVals[i][j];
                mG[i][j] = greenVals[i][j];
                mB[i][j] = blueVals[i][j];
            }
        }
        ImProcFunctions::Median medianTypeL = Median::TYPE_3X3_SOFT;

        int itera = 1;
        if(denstr < 0.2f) {
            medianTypeL = Median::TYPE_3X3_SOFT;
            itera = 1;
        } else if (denstr < 0.3f) {
            medianTypeL = Median::TYPE_3X3_SOFT;
            itera = 2;
        } else if (denstr < 0.45f) {
            medianTypeL = Median::TYPE_3X3_STRONG;
            itera = 2;
        } else if (denstr < 0.6f) {
            medianTypeL = Median::TYPE_3X3_STRONG;
            itera = 3;
        } else if (denstr < 0.8f) {
            medianTypeL = Median::TYPE_3X3_STRONG;
            itera = 4;
        } else if (denstr < 0.9f) {
            medianTypeL = Median::TYPE_5X5_STRONG;
            itera = 2;
        } else {
            medianTypeL = Median::TYPE_5X5_STRONG;
            itera = 3;            
        }
        ImProcFunctions::Median_Denoise(mR, mR, W, H, medianTypeL , itera, false, tmL);
        ImProcFunctions::Median_Denoise(mG, mG, W, H, medianTypeL , itera, false, tmL);
        ImProcFunctions::Median_Denoise(mB, mB, W, H, medianTypeL , itera, false, tmL);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                redVals[i][j] = intp(denstr, mR[i][j], redVals[i][j]); 
                greenVals[i][j] = intp(denstr, mG[i][j], greenVals[i][j]); 
                blueVals[i][j] = intp(denstr, mB[i][j], blueVals[i][j]); 
           }
        }


        for (int i = 0; i < hei; ++i) {
            delete[] tmL[i];
            delete[] mR[i];
            delete[] mG[i];
            delete[] mB[i];
        }

        delete[] tmL;
        delete[] mR;
        delete[] mG;
        delete[] mB;
    }
//end predenoise sharpening

    array2D<float> clipMask(W, H);
    constexpr float clipLimit = 0.95f;
    constexpr float maxSigma = 2.f;

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

    if (std::isnan(radius)) {
        return;
    }

    bool showMaskperso = sharpeningParams.showcap;


    if (showMask  ||  showMaskperso) {
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

        buildBlendMask(L, clipMask, W, H, contrast, sharpeningParams.autoContrast, clipMask);
        if (plistener) {
            plistener->setProgress(0.2);
        }
        conrastThreshold = contrast * 100.f;
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                red[i][j] = green[i][j] = blue[i][j] = clipMask[i][j] * 16384.f;
            }
        }
        if (plistener) {
            plistener->setProgress(1.0);
        }
        return;
    }

    std::unique_ptr<array2D<float>> Lbuffer;
    if (!redCache) {
        Lbuffer.reset(new array2D<float>(W, H));
    }

    std::unique_ptr<array2D<float>> YOldbuffer;
    if (!greenCache) {
        YOldbuffer.reset(new array2D<float>(W, H));
    }

    std::unique_ptr<array2D<float>> YNewbuffer;
    if (!blueCache) {
        YNewbuffer.reset(new array2D<float>(W, H));
    }
    array2D<float>& L = Lbuffer.get() ? *Lbuffer.get() : red;
    array2D<float>& YOld = YOldbuffer.get() ? *YOldbuffer.get() : green;
    array2D<float>& YNew = YNewbuffer.get() ? *YNewbuffer.get() : blue;

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
    buildBlendMask(L, clipMask, W, H, contrast, sharpeningParams.autoContrast, clipMask);
    if (plistener) {
        plistener->setProgress(0.2);
    }
    conrastThreshold = contrast * 100.f;
    CaptureDeconvSharpening(YNew, YOld, clipMask, W, H, radius, sharpeningParams.deconvradiusOffset, sharpeningParams.deconviter, sharpeningParams.deconvitercheck, plistener, 0.2, 0.9);
    if (plistener) {
        plistener->setProgress(0.9);
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < H; ++i) {
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 0; j < W; ++j) {
            const float factor = YNew[i][j] / std::max(YOld[i][j], 0.00001f);
            red[i][j] = redVals[i][j] * factor;
            green[i][j] = greenVals[i][j] * factor;
            blue[i][j] = blueVals[i][j] * factor;
        }
    }

    if (plistener) {
        plistener->setProgress(1.0);
    }
    rgbSourceModified = false;
}

// To call Capture Sharpening from Improcoordinator for Selective Editing
// call to 2 RAW functions of CaptureSharpening Raw - calcRadiusBayer and  calcRadiusXtrans
bool RawImageSource::getDeconvAutoRadius_capturesharpening_SE(float *out)
{
    const float clipVal = (ri->get_white(1) - ri->get_cblack(1)) * scale_mul[1];
    if (ri->getSensorType() == ST_BAYER) {//Bayer
        if (!out) {
            return true; // only check whether this is supported
        }
        const unsigned int fc[2] = {FC(0,0), FC(1,0)};
        *out = calcRadiusBayer(rawData, W, H, 1000.f, clipVal, fc);
        return true;
    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {//X trans
        if (!out) {
            return true; // only check whether this is supported
        }
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
        i-=7;
        j-=6;
        // std::cout << "found : " << found << std::endl;
        // std::cout << "i : " << i << std::endl;
        // std::cout << "j : " << j << std::endl;
        
        *out = calcRadiusXtrans(rawData, W, H, 1000.f, clipVal, i, j);
        return true;
    } else if (ri->get_colors() == 1) {
        if (out) {
            const unsigned int fc[2] = {0, 0};
            *out = calcRadiusBayer(rawData, W, H, 1000.f, clipVal, fc);
        }
        return true;
    }
    return false;
}


} /* namespace */
