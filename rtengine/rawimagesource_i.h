/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "rawimagesource.h"

namespace rtengine
{

inline void RawImageSource::convert_row_to_YIQ (const float* const r, const float* const g, const float* const b, float* Y, float* I, float* Q, const int W)
{
#ifdef _OPENMP
    #pragma omp simd
#endif

    for (int j = 0; j < W; j++) {
        Y[j] = .299f * r[j] + .587f * g[j] + .114f * b[j];
        I[j] = .596f * r[j] - .275f * g[j] - .321f * b[j];
        Q[j] = .212f * r[j] - .523f * g[j] + .311f * b[j];
    }
}

inline void RawImageSource::convert_row_to_RGB (float* r, float* g, float* b, const float* const Y, const float* const I, const float* const Q, const int W)
{
#ifdef _OPENMP
    #pragma omp simd
#endif

    for (int j = 1; j < W - 1; j++) {
        r[j] = Y[j] + 0.956f * I[j] + 0.621f * Q[j];
        g[j] = Y[j] - 0.272f * I[j] - 0.647f * Q[j];
        b[j] = Y[j] - 1.105f * I[j] + 1.702f * Q[j];
    }
}

inline void RawImageSource::convert_to_RGB (float &r, float &g, float &b, const float Y, const float I, const float Q)
{
    r = Y + 0.956f * I + 0.621f * Q;
    g = Y - 0.272f * I - 0.647f * Q;
    b = Y - 1.105f * I + 1.702f * Q;
}

inline void RawImageSource::interpolate_row_rb_mul_pp (const array2D<float> &rawData, float* ar, float* ab, float* pg, float* cg, float* ng, int i, float r_mul, float g_mul, float b_mul, int x1, int width, int skip)
{

    if ((ri->ISRED(i, 0) || ri->ISRED(i, 1)) && pg && ng) {
        // RGRGR or GRGRGR line
        for (int j = x1, jx = 0; jx < width; j += skip, jx++) {
            if (ri->ISRED(i, j)) {
                // red is simple
                ar[jx] = r_mul * rawData[i][j];
                // blue: cross interpolation
                float b = 0;
                int n = 0;

                if (i > 0 && j > 0) {
                    b += b_mul * rawData[i - 1][j - 1] - g_mul * pg[j - 1];
                    n++;
                }

                if (i > 0 && j < W - 1) {
                    b += b_mul * rawData[i - 1][j + 1] - g_mul * pg[j + 1];
                    n++;
                }

                if (i < H - 1 && j > 0) {
                    b += b_mul * rawData[i + 1][j - 1] - g_mul * ng[j - 1];
                    n++;
                }

                if (i < H - 1 && j < W - 1) {
                    b += b_mul * rawData[i + 1][j + 1] - g_mul * ng[j + 1];
                    n++;
                }

                b = g_mul * cg[j] + b / std::max(1, n);
                ab[jx] = std::max(0.f, b);
            } else {
                // linear R-G interp. horizontally
                float r;

                if (j == 0) {
                    r = g_mul * cg[0] + r_mul * rawData[i][1] - g_mul * cg[1];
                } else if (j == W - 1) {
                    r = g_mul * cg[W - 1] + r_mul * rawData[i][W - 2] - g_mul * cg[W - 2];
                } else {
                    r = g_mul * cg[j] + (r_mul * rawData[i][j - 1] - g_mul * cg[j - 1] + r_mul * rawData[i][j + 1] - g_mul * cg[j + 1]) / 2;
                }

                ar[jx] = std::max(0.f, r);
                // linear B-G interp. vertically
                float b;

                if (i == 0) {
                    b = g_mul * ng[j] + b_mul * rawData[1][j] - g_mul * cg[j];
                } else if (i == H - 1) {
                    b = g_mul * pg[j] + b_mul * rawData[H - 2][j] - g_mul * cg[j];
                } else {
                    b = g_mul * cg[j] + (b_mul * rawData[i - 1][j] - g_mul * pg[j] + b_mul * rawData[i + 1][j] - g_mul * ng[j]) / 2;
                }

                ab[jx] = std::max(0.f, b);
            }
        }
    } else if(pg && ng) {
        // BGBGB or GBGBGB line
        for (int j = x1, jx = 0; jx < width; j += skip, jx++) {
            if (ri->ISBLUE(i, j)) {
                // red is simple
                ab[jx] = b_mul * rawData[i][j];
                // blue: cross interpolation
                float r = 0;
                int n = 0;

                if (i > 0 && j > 0) {
                    r += r_mul * rawData[i - 1][j - 1] - g_mul * pg[j - 1];
                    n++;
                }

                if (i > 0 && j < W - 1) {
                    r += r_mul * rawData[i - 1][j + 1] - g_mul * pg[j + 1];
                    n++;
                }

                if (i < H - 1 && j > 0) {
                    r += r_mul * rawData[i + 1][j - 1] - g_mul * ng[j - 1];
                    n++;
                }

                if (i < H - 1 && j < W - 1) {
                    r += r_mul * rawData[i + 1][j + 1] - g_mul * ng[j + 1];
                    n++;
                }

                r = g_mul * cg[j] + r / std::max(n, 1);

                ar[jx] = std::max(0.f, r);
            } else {
                // linear B-G interp. horizontally
                float b;

                if (j == 0) {
                    b = g_mul * cg[0] + b_mul * rawData[i][1] - g_mul * cg[1];
                } else if (j == W - 1) {
                    b = g_mul * cg[W - 1] + b_mul * rawData[i][W - 2] - g_mul * cg[W - 2];
                } else {
                    b = g_mul * cg[j] + (b_mul * rawData[i][j - 1] - g_mul * cg[j - 1] + b_mul * rawData[i][j + 1] - g_mul * cg[j + 1]) / 2;
                }

                ab[jx] = std::max(0.f, b);
                // linear R-G interp. vertically
                float r;

                if (i == 0) {
                    r = g_mul * ng[j] + r_mul * rawData[1][j] - g_mul * cg[j];
                } else if (i == H - 1) {
                    r = g_mul * pg[j] + r_mul * rawData[H - 2][j] - g_mul * cg[j];
                } else {
                    r = g_mul * cg[j] + (r_mul * rawData[i - 1][j] - g_mul * pg[j] + r_mul * rawData[i + 1][j] - g_mul * ng[j]) / 2;
                }

                ar[jx] = std::max(0.f, r);
            }
        }
    }
}

}
