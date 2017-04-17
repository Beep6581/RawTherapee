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

#include "improcfun.h"
#include "rt_math.h"
#include "sleef.c"
#include "opthelper.h"
//#define PROFILE

#ifdef PROFILE
#   include <iostream>
#endif

namespace rtengine
{

static inline float Lanc(float x, float a)
{
    if (x * x < 1e-6f) {
        return 1.0f;
    } else if (x * x > a * a) {
        return 0.0f;
    } else {
        x = static_cast<float>(rtengine::RT_PI) * x;
        return a * xsinf(x) * xsinf(x / a) / (x * x);
    }
}

void ImProcFunctions::Lanczos(const Image16* src, Image16* dst, float scale)
{

    const float delta = 1.0f / scale;
    const float a = 3.0f;
    const float sc = min(scale, 1.0f);
    const int support = static_cast<int>(2.0f * a / sc) + 1;

    #pragma omp parallel
    {
        // storage for precomputed parameters for horisontal interpolation
        float * wwh = new float[support * dst->getWidth()];
        int * jj0 = new int[dst->getWidth()];
        int * jj1 = new int[dst->getWidth()];

        // temporal storage for vertically-interpolated row of pixels
        float * lr = new float[src->getWidth()];
        float * lg = new float[src->getWidth()];
        float * lb = new float[src->getWidth()];

        // Phase 1: precompute coefficients for horisontal interpolation

        for (int j = 0; j < dst->getWidth(); j++) {

            // x coord of the center of pixel on src image
            float x0 = (static_cast<float>(j) + 0.5f) * delta - 0.5f;

            // weights for interpolation in horisontal direction
            float * w = wwh + j * support;

            // sum of weights used for normalization
            float ws = 0.0f;

            jj0[j] = max(0, static_cast<int>(floorf(x0 - a / sc)) + 1);
            jj1[j] = min(src->getWidth(), static_cast<int>(floorf(x0 + a / sc)) + 1);

            // calculate weights
            for (int jj = jj0[j]; jj < jj1[j]; jj++) {
                int k = jj - jj0[j];
                float z = sc * (x0 - static_cast<float>(jj));
                w[k] = Lanc(z, a);
                ws += w[k];
            }

            // normalize weights
            for (int k = 0; k < support; k++) {
                w[k] /= ws;
            }
        }

        // Phase 2: do actual interpolation
        #pragma omp for

        for (int i = 0; i < dst->getHeight(); i++) {

            // y coord of the center of pixel on src image
            float y0 = (static_cast<float>(i) + 0.5f) * delta - 0.5f;

            // weights for interpolation in y direction
            float w[support];

            // sum of weights used for normalization
            float ws = 0.0f;

            int ii0 = max(0, static_cast<int>(floorf(y0 - a / sc)) + 1);
            int ii1 = min(src->getHeight(), static_cast<int>(floorf(y0 + a / sc)) + 1);

            // calculate weights for vertical interpolation
            for (int ii = ii0; ii < ii1; ii++) {
                int k = ii - ii0;
                float z = sc * (y0 - static_cast<float>(ii));
                w[k] = Lanc(z, a);
                ws += w[k];
            }

            // normalize weights
            for (int k = 0; k < support; k++) {
                w[k] /= ws;
            }

            // Do vertical interpolation. Store results.
            for (int j = 0; j < src->getWidth(); j++) {

                float r = 0.0f, g = 0.0f, b = 0.0f;

                for (int ii = ii0; ii < ii1; ii++) {
                    int k = ii - ii0;

                    r += w[k] * src->r(ii, j);
                    g += w[k] * src->g(ii, j);
                    b += w[k] * src->b(ii, j);
                }

                lr[j] = r;
                lg[j] = g;
                lb[j] = b;
            }

            // Do horizontal interpolation
            for(int j = 0; j < dst->getWidth(); j++) {

                float * wh = wwh + support * j;

                float r = 0.0f, g = 0.0f, b = 0.0f;

                for (int jj = jj0[j]; jj < jj1[j]; jj++) {
                    int k = jj - jj0[j];

                    r += wh[k] * lr[jj];
                    g += wh[k] * lg[jj];
                    b += wh[k] * lb[jj];
                }

                dst->r(i, j) = CLIP(static_cast<int>(r));
                dst->g(i, j) = CLIP(static_cast<int>(g));
                dst->b(i, j) = CLIP(static_cast<int>(b));
            }
        }

        delete[] wwh;
        delete[] jj0;
        delete[] jj1;
        delete[] lr;
        delete[] lg;
        delete[] lb;
    }
}


SSEFUNCTION void ImProcFunctions::Lanczos(const LabImage* src, LabImage* dst, float scale)
{
    const float delta = 1.0f / scale;
    const float a = 3.0f;
    const float sc = min(scale, 1.0f);
    const int support = static_cast<int>(2.0f * a / sc) + 1;

    // storage for precomputed parameters for horizontal interpolation
    float * wwh = new float[support * dst->W];
    int * jj0 = new int[dst->W];
    int * jj1 = new int[dst->W];

    // Phase 1: precompute coefficients for horizontal interpolation
    for (int j = 0; j < dst->W; j++) {

        // x coord of the center of pixel on src image
        float x0 = (static_cast<float>(j) + 0.5f) * delta - 0.5f;

        // weights for interpolation in horizontal direction
        float * w = wwh + j * support;

        // sum of weights used for normalization
        float ws = 0.0f;

        jj0[j] = max(0, static_cast<int>(floorf(x0 - a / sc)) + 1);
        jj1[j] = min(src->W, static_cast<int>(floorf(x0 + a / sc)) + 1);

        // calculate weights
        for (int jj = jj0[j]; jj < jj1[j]; jj++) {
            int k = jj - jj0[j];
            float z = sc * (x0 - static_cast<float>(jj));
            w[k] = Lanc(z, a);
            ws += w[k];
        }

        // normalize weights
        for (int k = 0; k < support; k++) {
            w[k] /= ws;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        // temporal storage for vertically-interpolated row of pixels
        float * lL = new float[src->W];
        float * la = new float[src->W];
        float * lb = new float[src->W];
        // weights for interpolation in y direction
        float w[support] ALIGNED64;

        // Phase 2: do actual interpolation
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < dst->H; i++) {
            // y coord of the center of pixel on src image
            float y0 = (static_cast<float>(i) + 0.5f) * delta - 0.5f;

            // sum of weights used for normalization
            float ws = 0.0f;

            int ii0 = max(0, static_cast<int>(floorf(y0 - a / sc)) + 1);
            int ii1 = min(src->H, static_cast<int>(floorf(y0 + a / sc)) + 1);

            // calculate weights for vertical interpolation
            for (int ii = ii0; ii < ii1; ii++) {
                int k = ii - ii0;
                float z = sc * (y0 - static_cast<float>(ii));
                w[k] = Lanc(z, a);
                ws += w[k];
            }

            // normalize weights
            for (int k = 0; k < support; k++) {
                w[k] /= ws;
            }

            // Do vertical interpolation. Store results.
#ifdef __SSE2__
            int j;
            __m128 Lv, av, bv, wkv;

            for (j = 0; j < src->W - 3; j += 4) {
                Lv = _mm_setzero_ps();
                av = _mm_setzero_ps();
                bv = _mm_setzero_ps();

                for (int ii = ii0; ii < ii1; ii++) {
                    int k = ii - ii0;
                    wkv = _mm_set1_ps(w[k]);
                    Lv += wkv * LVFU(src->L[ii][j]);
                    av += wkv * LVFU(src->a[ii][j]);
                    bv += wkv * LVFU(src->b[ii][j]);
                }

                STVF(lL[j], Lv);
                STVF(la[j], av);
                STVF(lb[j], bv);
            }

#else
            int j = 0;
#endif

            for (; j < src->W; j++) {
                float L = 0.0f, a = 0.0f, b = 0.0f;

                for (int ii = ii0; ii < ii1; ii++) {
                    int k = ii - ii0;

                    L += w[k] * src->L[ii][j];
                    a += w[k] * src->a[ii][j];
                    b += w[k] * src->b[ii][j];
                }

                lL[j] = L;
                la[j] = a;
                lb[j] = b;
            }

            // Do horizontal interpolation
            for(int j = 0; j < dst->W; j++) {

                float * wh = wwh + support * j;

                float L = 0.0f, a = 0.0f, b = 0.0f;

                for (int jj = jj0[j]; jj < jj1[j]; jj++) {
                    int k = jj - jj0[j];

                    L += wh[k] * lL[jj];
                    a += wh[k] * la[jj];
                    b += wh[k] * lb[jj];
                }

                dst->L[i][j] = L;
                dst->a[i][j] = a;
                dst->b[i][j] = b;
            }
        }

        delete[] lL;
        delete[] la;
        delete[] lb;
    }
    delete[] jj0;
    delete[] jj1;
    delete[] wwh;
}

float ImProcFunctions::resizeScale (const ProcParams* params, int fw, int fh, int &imw, int &imh)
{
    imw = fw;
    imh = fh;

    if (!params || !params->resize.enabled) {
        return 1.0;
    }

    // get the resize parameters
    int refw, refh;
    double dScale;

    if (params->crop.enabled && params->resize.appliesTo == "Cropped area") {
        // the resize values applies to the crop dimensions
        refw = params->crop.w;
        refh = params->crop.h;
    } else {
        // the resize values applies to the image dimensions
        // if a crop exists, it will be resized to the calculated scale
        refw = fw;
        refh = fh;
    }

    switch(params->resize.dataspec) {
    case (1):
        // Width
        dScale = (double)params->resize.width / (double)refw;
        break;

    case (2):
        // Height
        dScale = (double)params->resize.height / (double)refh;
        break;

    case (3):

        // FitBox
        if ((double)refw / (double)refh > (double)params->resize.width / (double)params->resize.height) {
            dScale = (double)params->resize.width / (double)refw;
        } else {
            dScale = (double)params->resize.height / (double)refh;
        }

        break;

    default:
        // Scale
        dScale = params->resize.scale;
        break;
    }

    if (fabs(dScale - 1.0) <= 1e-5) {
        return 1.0;
    }

    if (params->crop.enabled && params->resize.appliesTo == "Full image") {
        imw = params->crop.w;
        imh = params->crop.h;
    } else {
        imw = refw;
        imh = refh;
    }

    imw = (int)( (double)imw * dScale + 0.5 );
    imh = (int)( (double)imh * dScale + 0.5 );
    return (float)dScale;
}

void ImProcFunctions::resize (Image16* src, Image16* dst, float dScale)
{
#ifdef PROFILE
    time_t t1 = clock();
#endif

    if(params->resize.method != "Nearest" ) {
        Lanczos(src, dst, dScale);
    } else {
        // Nearest neighbour algorithm
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif

        for (int i = 0; i < dst->getHeight(); i++) {
            int sy = i / dScale;
            sy = LIM(sy, 0, src->getHeight() - 1);

            for (int j = 0; j < dst->getWidth(); j++) {
                int sx = j / dScale;
                sx = LIM(sx, 0, src->getWidth() - 1);
                dst->r(i, j) = src->r(sy, sx);
                dst->g(i, j) = src->g(sy, sx);
                dst->b(i, j) = src->b(sy, sx);
            }
        }
    }

#ifdef PROFILE
    time_t t2 = clock();
    std::cout << "Resize: " << params->resize.method << ": "
              << (float)(t2 - t1) / CLOCKS_PER_SEC << std::endl;
#endif
}

}
