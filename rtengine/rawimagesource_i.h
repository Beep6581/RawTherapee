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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RAWIMAGESOURCE_I_H_INCLUDED
#define RAWIMAGESOURCE_I_H_INCLUDED

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

inline void RawImageSource::convert_to_cielab_row (float* ar, float* ag, float* ab, float* oL, float* oa, float* ob)
{

    for (int j = 0; j < W; j++) {
        float r = ar[j];
        float g = ag[j];
        float b = ab[j];

        float x = lc00 * r + lc01 * g + lc02 * b;
        float y = lc10 * r + lc11 * g + lc12 * b;
        float z = lc20 * r + lc21 * g + lc22 * b;

        if (y > threshold) {
            oL[j] = cache[(int)y];
        } else {
            oL[j] = float(903.3f * y / MAXVALF);
        }

        oa[j] = 500.f * ((x > threshold ? cache[(int)x] : 7.787f * x / MAXVALF + 16.f / 116.f) - (y > threshold ? cache[(int)y] : 7.787f * y / MAXVALF + 16.f / 116.f));
        ob[j] = 200.f * ((y > threshold ? cache[(int)y] : 7.787f * y / MAXVALF + 16.f / 116.f) - (z > threshold ? cache[(int)z] : 7.787f * z / MAXVALF + 16.f / 116.f));
    }
}

inline void RawImageSource::interpolate_row_g (float* agh, float* agv, int i)
{

    for (int j = 0; j < W; j++) {
        if (ri->ISGREEN(i, j)) {
            agh[j] = rawData[i][j];
            agv[j] = rawData[i][j];
        } else {
            int gh = 0;
            int gv = 0;

            if (j > 1 && j < W - 2) {
                gh = (-rawData[i][j - 2] + 2 * rawData[i][j - 1] + 2 * rawData[i][j] + 2 * rawData[i][j + 1] - rawData[i][j + 2]) / 4;
                int maxgh = max(rawData[i][j - 1], rawData[i][j + 1]);
                int mingh = min(rawData[i][j - 1], rawData[i][j + 1]);

                if (gh > maxgh) {
                    gh = maxgh;
                } else if (gh < mingh) {
                    gh = mingh;
                }
            } else if (j == 0) {
                gh = rawData[i][1];
            } else if (j == 1) {
                gh = (rawData[i][0] + rawData[i][2]) / 2;
            } else if (j == W - 1) {
                gh = rawData[i][W - 2];
            } else if (j == W - 2) {
                gh = (rawData[i][W - 1] + rawData[i][W - 3]) / 2;
            }

            if (i > 1 && i < H - 2) {
                gv = (-rawData[i - 2][j] + 2 * rawData[i - 1][j] + 2 * rawData[i][j] + 2 * rawData[i + 1][j] - rawData[i + 2][j]) / 4;
                int maxgv = max(rawData[i - 1][j], rawData[i + 1][j]);
                int mingv = min(rawData[i - 1][j], rawData[i + 1][j]);

                if (gv > maxgv) {
                    gv = maxgv;
                } else if (gv < mingv) {
                    gv = mingv;
                }
            } else if (i == 0) {
                gv = rawData[1][j];
            } else if (i == 1) {
                gv = (rawData[0][j] + rawData[2][j]) / 2;
            } else if (i == H - 1) {
                gv = rawData[H - 2][j];
            } else if (i == H - 2) {
                gv = (rawData[H - 1][j] + rawData[H - 3][j]) / 2;
            }

            agh[j] = gh;
            agv[j] = gv;
        }
    }
}

inline void RawImageSource::interpolate_row_rb (float* ar, float* ab, float* pg, float* cg, float* ng, int i)
{
    const auto getPg = [pg](int index) {
        return
            pg
                ? pg[index]
                : 0.f;
    };

    const auto getNg = [ng](int index) {
        return
            ng
                ? ng[index]
                : 0.f;
    };

    float *nonGreen1 = ar;
    float *nonGreen2 = ab;

    if ((ri->ISBLUE(i, 0) || ri->ISBLUE(i, 1))) {
        std::swap(nonGreen1, nonGreen2);
    }
    int j = FC(i, 0) & 1;
    if (j) {
        // linear R-G interp. horizontally
        float val1;

        val1 = cg[0] + rawData[i][1] - cg[1];

        nonGreen1[0] = CLIP(val1);
        // linear B-G interp. vertically
        float val2;

        if (i == 0) {
            val2 = getNg(0) + rawData[1][0] - cg[0];
        } else if (i == H - 1) {
            val2 = getPg(0) + rawData[H - 2][0] - cg[0];
        } else {
            val2 = cg[0] + (rawData[i - 1][0] - getPg(0) + rawData[i + 1][0] - getNg(0)) / 2;
        }

        nonGreen2[0] = val2;
    }
    // RGRGR or GRGRGR line
    for (; j < W - 1; j += 2) {
        // nonGreen of row is simple
        nonGreen1[j] = rawData[i][j];
        // non green of next row: cross interpolation
        float nonGreen = 0.f;
        int n = 0;

        if (i > 0) {
            if (j > 0) {
                nonGreen += rawData[i - 1][j - 1] - getPg(j - 1);
                n++;
            }
            nonGreen += rawData[i - 1][j + 1] - getPg(j + 1);
            n++;
        }

        if (i < H - 1) {
            if (j > 0) {
                nonGreen += rawData[i + 1][j - 1] - getNg(j - 1);
                n++;
            }
            nonGreen += rawData[i + 1][j + 1] - getNg(j + 1);
            n++;
        }

        nonGreen2[j] = cg[j] + nonGreen / n;

        // linear R-G interp. horizontally
        float val1;

        if (j == W - 2) {
            val1 = cg[W - 1] + rawData[i][W - 2] - cg[W - 2];
        } else {
            val1 = cg[j + 1] + (rawData[i][j] - cg[j] + rawData[i][j + 2] - cg[j + 2]) / 2;
        }

        nonGreen1[j + 1] = CLIP(val1);
        // linear B-G interp. vertically
        float val2;

        if (i == 0) {
            val2 = getNg(j + 1) + rawData[1][j + 1] - cg[j + 1];
        } else if (i == H - 1) {
            val2 = getPg(j + 1) + rawData[H - 2][j + 1] - cg[j + 1];
        } else {
            val2 = cg[j + 1] + (rawData[i - 1][j + 1] - getPg(j + 1) + rawData[i + 1][j + 1] - getNg(j + 1)) / 2;
        }

        nonGreen2[j + 1] = val2;
    }

    if(j == W - 1) {
        // nonGreen of row is simple
        nonGreen1[j] = rawData[i][j];
        // non green of next row: cross interpolation
        float nonGreen = 0.f;
        int n = 0;

        if (i > 0) {
            nonGreen += rawData[i - 1][j - 1] - getPg(j - 1);
            n++;
        }

        if (i < H - 1) {
            nonGreen += rawData[i + 1][j - 1] - getNg(j - 1);
            n++;
        }

        nonGreen2[j] = cg[j] + nonGreen / n;
    }
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
                ab[jx] = b;
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

                ar[jx] = r;
                // linear B-G interp. vertically
                float b;

                if (i == 0) {
                    b = g_mul * ng[j] + b_mul * rawData[1][j] - g_mul * cg[j];
                } else if (i == H - 1) {
                    b = g_mul * pg[j] + b_mul * rawData[H - 2][j] - g_mul * cg[j];
                } else {
                    b = g_mul * cg[j] + (b_mul * rawData[i - 1][j] - g_mul * pg[j] + b_mul * rawData[i + 1][j] - g_mul * ng[j]) / 2;
                }

                ab[jx] = b;
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

                ar[jx] = r;
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

                ab[jx] = b;
                // linear R-G interp. vertically
                float r;

                if (i == 0) {
                    r = g_mul * ng[j] + r_mul * rawData[1][j] - g_mul * cg[j];
                } else if (i == H - 1) {
                    r = g_mul * pg[j] + r_mul * rawData[H - 2][j] - g_mul * cg[j];
                } else {
                    r = g_mul * cg[j] + (r_mul * rawData[i - 1][j] - g_mul * pg[j] + r_mul * rawData[i + 1][j] - g_mul * ng[j]) / 2;
                }

                ar[jx] = r;
            }
        }
    }
}

}

#endif
