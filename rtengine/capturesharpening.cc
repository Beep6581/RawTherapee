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

namespace {

void compute13x13kernel(float sigma, float kernel[13][13]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -6; i <= 6; ++i) {
        for (int j = -6; j <= 6; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 2.0)) {
                kernel[i + 6][j + 6] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 6][j + 6];
            } else {
                kernel[i + 6][j + 6] = 0.f;
            }
        }
    }

    for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 13; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute9x9kernel(float sigma, float kernel[9][9]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -4; i <= 4; ++i) {
        for (int j = -4; j <= 4; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 1.5)) {
                kernel[i + 4][j + 4] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 4][j + 4];
            } else {
                kernel[i + 4][j + 4] = 0.f;
            }
        }
    }

    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

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
            kernel[i + 1][j + 1] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
            sum += kernel[i + 1][j + 1];
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void gauss3x3div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[3][3])
{

    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss5x5div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss7x7div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[7][7])
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
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss9x9div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[9][9])
{

    const float c42 = kernel[0][2];
    const float c41 = kernel[0][3];
    const float c40 = kernel[0][4];
    const float c33 = kernel[1][1];
    const float c32 = kernel[1][2];
    const float c31 = kernel[1][3];
    const float c30 = kernel[1][4];
    const float c22 = kernel[2][2];
    const float c21 = kernel[2][3];
    const float c20 = kernel[2][4];
    const float c11 = kernel[3][3];
    const float c10 = kernel[3][4];
    const float c00 = kernel[4][4];

    for (int i = 4; i < tileSize - 4; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 4; j < tileSize - 4; ++j) {
            const float val = c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss13x13div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[13][13])
{
    const float c60 = kernel[0][6];
    const float c53 = kernel[1][3];
    const float c52 = kernel[1][4];
    const float c51 = kernel[1][5];
    const float c50 = kernel[1][6];
    const float c44 = kernel[2][2];
    const float c42 = kernel[2][4];
    const float c41 = kernel[2][5];
    const float c40 = kernel[2][6];
    const float c33 = kernel[3][3];
    const float c32 = kernel[3][4];
    const float c31 = kernel[3][5];
    const float c30 = kernel[3][6];
    const float c22 = kernel[4][4];
    const float c21 = kernel[4][5];
    const float c20 = kernel[4][6];
    const float c11 = kernel[5][5];
    const float c10 = kernel[5][6];
    const float c00 = kernel[6][6];

    for (int i = 6; i < tileSize - 6; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 6; j < tileSize - 6; ++j) {
            const float val = c60 * (src[i - 6][j] + src[i][j - 6] + src[i][j + 6] + src[i + 6][j]) +
                              c53 * ((src[i - 5][j - 3] + src[i - 5][j + 3]) + (src[i - 3][j - 5] + src[i - 3][j + 5]) + (src[i + 3][j - 5] + src[i + 3][j + 5]) + (src[i + 5][j - 3] + src[i + 5][j + 3])) +
                              c52 * ((src[i - 5][j - 2] + src[i - 5][j + 2]) + (src[i - 2][j - 5] + src[i - 2][j + 5]) + (src[i + 2][j - 5] + src[i + 2][j + 5]) + (src[i + 5][j - 2] + src[i + 5][j + 2])) +
                              c51 * ((src[i - 5][j - 1] + src[i - 5][j + 1]) + (src[i - 1][j - 5] + src[i - 1][j + 5]) + (src[i + 1][j - 5] + src[i + 1][j + 5]) + (src[i + 5][j - 1] + src[i + 5][j + 1])) +
                              c50 * ((src[i - 5][j] + src[i][j - 5] + src[i][j + 5] + src[i + 5][j]) + ((src[i - 4][j - 3] + src[i - 4][j + 3]) + (src[i - 3][j - 4] + src[i - 3][j + 4]) + (src[i + 3][j - 4] + src[i + 3][j + 4]) + (src[i + 4][j - 3] + src[i + 4][j + 3]))) +
                              c44 * (src[i - 4][j - 4] + src[i - 4][j + 4] + src[i + 4][j - 4] + src[i + 4][j + 4]) +
                              c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss3x3mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[3][3])
{
    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] *= val;
        }
    }

}

void gauss5x5mult (float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

void gauss7x7mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[7][7])
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
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

void gauss9x9mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[9][9])
{

    const float c42 = kernel[0][2];
    const float c41 = kernel[0][3];
    const float c40 = kernel[0][4];
    const float c33 = kernel[1][1];
    const float c32 = kernel[1][2];
    const float c31 = kernel[1][3];
    const float c30 = kernel[1][4];
    const float c22 = kernel[2][2];
    const float c21 = kernel[2][3];
    const float c20 = kernel[2][4];
    const float c11 = kernel[3][3];
    const float c10 = kernel[3][4];
    const float c00 = kernel[4][4];

    for (int i = 4; i < tileSize - 4; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 4; j < tileSize - 4; ++j) {
            const float val = c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];
            dst[i][j] *= val;
        }
    }
}

void gauss13x13mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[13][13])
{

    const float c60 = kernel[0][6];
    const float c53 = kernel[1][3];
    const float c52 = kernel[1][4];
    const float c51 = kernel[1][5];
    const float c50 = kernel[1][6];
    const float c44 = kernel[2][2];
    const float c42 = kernel[2][4];
    const float c41 = kernel[2][5];
    const float c40 = kernel[2][6];
    const float c33 = kernel[3][3];
    const float c32 = kernel[3][4];
    const float c31 = kernel[3][5];
    const float c30 = kernel[3][6];
    const float c22 = kernel[4][4];
    const float c21 = kernel[4][5];
    const float c20 = kernel[4][6];
    const float c11 = kernel[5][5];
    const float c10 = kernel[5][6];
    const float c00 = kernel[6][6];

    for (int i = 6; i < tileSize - 6; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 6; j < tileSize - 6; ++j) {
            const float val = c60 * (src[i - 6][j] + src[i][j - 6] + src[i][j + 6] + src[i + 6][j]) +
                              c53 * ((src[i - 5][j - 3] + src[i - 5][j + 3]) + (src[i - 3][j - 5] + src[i - 3][j + 5]) + (src[i + 3][j - 5] + src[i + 3][j + 5]) + (src[i + 5][j - 3] + src[i + 5][j + 3])) +
                              c52 * ((src[i - 5][j - 2] + src[i - 5][j + 2]) + (src[i - 2][j - 5] + src[i - 2][j + 5]) + (src[i + 2][j - 5] + src[i + 2][j + 5]) + (src[i + 5][j - 2] + src[i + 5][j + 2])) +
                              c51 * ((src[i - 5][j - 1] + src[i - 5][j + 1]) + (src[i - 1][j - 5] + src[i - 1][j + 5]) + (src[i + 1][j - 5] + src[i + 1][j + 5]) + (src[i + 5][j - 1] + src[i + 5][j + 1])) +
                              c50 * ((src[i - 5][j] + src[i][j - 5] + src[i][j + 5] + src[i + 5][j]) + ((src[i - 4][j - 3] + src[i - 4][j + 3]) + (src[i - 3][j - 4] + src[i - 3][j + 4]) + (src[i + 3][j - 4] + src[i + 3][j + 4]) + (src[i + 4][j - 3] + src[i + 4][j + 3]))) +
                              c44 * (src[i - 4][j - 4] + src[i - 4][j + 4] + src[i + 4][j - 4] + src[i + 4][j + 4]) +
                              c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
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
        compute3x3kernel(sigma, kernel3);
    } else if (is5x5) {
        compute5x5kernel(sigma, kernel5);
    } else if (is7x7) {
        compute7x7kernel(sigma, kernel7);
    } else if (is9x9) {
        compute9x9kernel(sigma, kernel9);
    } else {
        compute13x13kernel(sigma, kernel13);
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
                        gauss3x3div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel3);
                        gauss3x3mult(tmpThr, tmpIThr, fullTileSize, kernel3);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is5x5) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss5x5div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel5);
                        gauss5x5mult(tmpThr, tmpIThr, fullTileSize, kernel5);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is7x7) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel7);
                        gauss7x7mult(tmpThr, tmpIThr, fullTileSize, kernel7);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is9x9) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss9x9div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel9);
                        gauss9x9mult(tmpThr, tmpIThr, fullTileSize, kernel9);
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
                                compute13x13kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel13);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                                    gauss13x13div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel13);
                                    gauss13x13mult(tmpThr, tmpIThr, fullTileSize, lkernel13);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 1.15f) { // have to use 9x9 kernel
                                float lkernel9[9][9];
                                compute9x9kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel9);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 9x9 gaussian blur and divide luminance by result of gaussian blur
                                    gauss9x9div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel9);
                                    gauss9x9mult(tmpThr, tmpIThr, fullTileSize, lkernel9);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 0.84f) { // have to use 7x7 kernel
                                float lkernel7[7][7];
                                compute7x7kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel7);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel7);
                                    gauss7x7mult(tmpThr, tmpIThr, fullTileSize, lkernel7);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else { // can use 5x5 kernel
                                float lkernel5[5][5];
                                compute5x5kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel5);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    gauss5x5div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel5);
                                    gauss5x5mult(tmpThr, tmpIThr, fullTileSize, lkernel5);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        for (int k = 0; k < iterations; ++k) {
                            // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                            gauss13x13div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel13);
                            gauss13x13mult(tmpThr, tmpIThr, fullTileSize, kernel13);
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
        plistener->setProgressStr(M("TP_PDSHARPENING_LABEL"));
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

    const array2D<float>& redVals = redCache ? *redCache : red;
    const array2D<float>& greenVals = greenCache ? *greenCache : green;
    const array2D<float>& blueVals = blueCache ? *blueCache : blue;

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

} /* namespace */
