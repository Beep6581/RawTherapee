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

#include "improcfun.h"

#include "alignedbuffer.h"
#include "imagefloat.h"
#include "labimage.h"
#include "opthelper.h"
#include "rt_math.h"
#include "procparams.h"
#include "sleef.h"

//#define PROFILE

#ifdef PROFILE
#   include <iostream>
#endif


namespace rtengine
{

static inline float Lanc (float x, float a)
{
    if (x * x < 1e-6f) {
        return 1.0f;
    } else if (x * x > a * a) {
        return 0.0f;
    } else {
        x = static_cast<float> (rtengine::RT_PI) * x;
        return a * xsinf (x) * xsinf (x / a) / (x * x);
    }
}

void ImProcFunctions::Lanczos (const Imagefloat* src, Imagefloat* dst, float scale)
{

    const float delta = 1.0f / scale;
    const float a = 3.0f;
    const float sc = min (scale, 1.0f);
    const int support = static_cast<int> (2.0f * a / sc) + 1;

#ifdef _OPENMP
    #pragma omp parallel
#endif
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
            float x0 = (static_cast<float> (j) + 0.5f) * delta - 0.5f;

            // weights for interpolation in horisontal direction
            float * w = wwh + j * support;

            // sum of weights used for normalization
            float ws = 0.0f;

            jj0[j] = max (0, static_cast<int> (floorf (x0 - a / sc)) + 1);
            jj1[j] = min (src->getWidth(), static_cast<int> (floorf (x0 + a / sc)) + 1);

            // calculate weights
            for (int jj = jj0[j]; jj < jj1[j]; jj++) {
                int k = jj - jj0[j];
                float z = sc * (x0 - static_cast<float> (jj));
                w[k] = Lanc (z, a);
                ws += w[k];
            }

            // normalize weights
            for (int k = 0; k < support; k++) {
                w[k] /= ws;
            }
        }

        // Phase 2: do actual interpolation
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < dst->getHeight(); i++) {

            // y coord of the center of pixel on src image
            float y0 = (static_cast<float> (i) + 0.5f) * delta - 0.5f;

            // weights for interpolation in y direction
            float w[support];
            for (auto& f : w) {
                f = 0.f;
            }

            // sum of weights used for normalization
            float ws = 0.0f;

            int ii0 = max (0, static_cast<int> (floorf (y0 - a / sc)) + 1);
            int ii1 = min (src->getHeight(), static_cast<int> (floorf (y0 + a / sc)) + 1);

            // calculate weights for vertical interpolation
            for (int ii = ii0; ii < ii1; ii++) {
                int k = ii - ii0;
                float z = sc * (y0 - static_cast<float> (ii));
                w[k] = Lanc (z, a);
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

                    r += w[k] * src->r (ii, j);
                    g += w[k] * src->g (ii, j);
                    b += w[k] * src->b (ii, j);
                }

                lr[j] = r;
                lg[j] = g;
                lb[j] = b;
            }

            // Do horizontal interpolation
            for (int j = 0; j < dst->getWidth(); j++) {

                float * wh = wwh + support * j;

                float r = 0.0f, g = 0.0f, b = 0.0f;

                for (int jj = jj0[j]; jj < jj1[j]; jj++) {
                    int k = jj - jj0[j];

                    r += wh[k] * lr[jj];
                    g += wh[k] * lg[jj];
                    b += wh[k] * lb[jj];
                }

                dst->r (i, j) = /*CLIP*/ (r);//static_cast<int> (r));
                dst->g (i, j) = /*CLIP*/ (g);//static_cast<int> (g));
                dst->b (i, j) = /*CLIP*/ (b);//static_cast<int> (b));
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


void ImProcFunctions::Lanczos (const LabImage* src, LabImage* dst, float scale)
{
    const float delta = 1.0f / scale;
    constexpr float a = 3.0f;
    const float sc = min(scale, 1.0f);
    const int support = static_cast<int> (2.0f * a / sc) + 1;

    // storage for precomputed parameters for horizontal interpolation
    float* const wwh = new float[support * dst->W];
    int* const jj0 = new int[dst->W];
    int* const jj1 = new int[dst->W];

    // Phase 1: precompute coefficients for horizontal interpolation
    for (int j = 0; j < dst->W; j++) {

        // x coord of the center of pixel on src image
        float x0 = (static_cast<float> (j) + 0.5f) * delta - 0.5f;

        // weights for interpolation in horizontal direction
        float * w = wwh + j * support;

        // sum of weights used for normalization
        float ws = 0.0f;

        jj0[j] = max (0, static_cast<int> (floorf (x0 - a / sc)) + 1);
        jj1[j] = min (src->W, static_cast<int> (floorf (x0 + a / sc)) + 1);

        // calculate weights
        for (int jj = jj0[j]; jj < jj1[j]; jj++) {
            int k = jj - jj0[j];
            float z = sc * (x0 - static_cast<float> (jj));
            w[k] = Lanc (z, a);
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
        AlignedBuffer<float> aligned_buffer_ll(src->W);
        AlignedBuffer<float> aligned_buffer_la(src->W);
        AlignedBuffer<float> aligned_buffer_lb(src->W);
        float* const lL = aligned_buffer_ll.data;
        float* const la = aligned_buffer_la.data;
        float* const lb = aligned_buffer_lb.data;
        // weights for interpolation in y direction
        float w[support] ALIGNED64;
        memset(w, 0, sizeof(w));

        // Phase 2: do actual interpolation
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < dst->H; i++) {
            // y coord of the center of pixel on src image
            float y0 = (static_cast<float> (i) + 0.5f) * delta - 0.5f;

            // sum of weights used for normalization
            float ws = 0.0f;

            int ii0 = max (0, static_cast<int> (floorf (y0 - a / sc)) + 1);
            int ii1 = min (src->H, static_cast<int> (floorf (y0 + a / sc)) + 1);

            // calculate weights for vertical interpolation
            for (int ii = ii0; ii < ii1; ii++) {
                int k = ii - ii0;
                float z = sc * (y0 - static_cast<float> (ii));
                w[k] = Lanc (z, a);
                ws += w[k];
            }

            // normalize weights
            for (int k = 0; k < support; k++) {
                w[k] /= ws;
            }

            // Do vertical interpolation. Store results.
            int j = 0;
#ifdef __SSE2__
            __m128 Lv, av, bv, wkv;

            for (j = 0; j < src->W - 3; j += 4) {
                Lv = ZEROV;
                av = ZEROV;
                bv = ZEROV;

                for (int ii = ii0; ii < ii1; ii++) {
                    int k = ii - ii0;
                    wkv = F2V(w[k]);
                    Lv += wkv * LVFU(src->L[ii][j]);
                    av += wkv * LVFU(src->a[ii][j]);
                    bv += wkv * LVFU(src->b[ii][j]);
                }

                STVF(lL[j], Lv);
                STVF(la[j], av);
                STVF(lb[j], bv);
            }
#endif

            for (; j < src->W; ++j) {
                float Ll = 0.0f, La = 0.0f, Lb = 0.0f;

                for (int ii = ii0; ii < ii1; ++ii) {
                    int k = ii - ii0;

                    Ll += w[k] * src->L[ii][j];
                    La += w[k] * src->a[ii][j];
                    Lb += w[k] * src->b[ii][j];
                }

                lL[j] = Ll;
                la[j] = La;
                lb[j] = Lb;
            }

            // Do horizontal interpolation
            for (int x = 0; x < dst->W; ++x) {
                float * wh = wwh + support * x;
                float Ll = 0.0f, La = 0.0f, Lb = 0.0f;

                for (int jj = jj0[x]; jj < jj1[x]; ++jj) {
                    int k = jj - jj0[x];

                    Ll += wh[k] * lL[jj];
                    La += wh[k] * la[jj];
                    Lb += wh[k] * lb[jj];
                }

                dst->L[i][x] = Ll;
                dst->a[i][x] = La;
                dst->b[i][x] = Lb;
            }
        }
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

    switch (params->resize.dataspec) {
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
            dScale = (dScale > 1.0 && !params->resize.allowUpscaling) ? 1.0 : dScale;

            break;
            
        case (4):
        
            // Long Edge
            if (refw > refh) {
                dScale = (double)params->resize.longedge / (double)refw;
            } else {
                dScale = (double)params->resize.longedge / (double)refh;
            }
            break;
            
        case (5):
        
            // Short Edge
            if (refw > refh) {
                dScale = (double)params->resize.shortedge / (double)refh;
            } else {
                dScale = (double)params->resize.shortedge / (double)refw;
            }
            break;

        default:
            // Scale
            dScale = params->resize.scale;
            break;
    }

    if (fabs (dScale - 1.0) <= 1e-5) {
        return 1.0;
    }

    if (params->crop.enabled && params->resize.appliesTo == "Full image") {
        imw = params->crop.w;
        imh = params->crop.h;
    } else {
        imw = refw;
        imh = refh;
    }

    imw = (int) ( (double)imw * dScale + 0.5 );
    imh = (int) ( (double)imh * dScale + 0.5 );
    return (float)dScale;
}

void ImProcFunctions::resize (Imagefloat* src, Imagefloat* dst, float dScale)
{
#ifdef PROFILE
    time_t t1 = clock();
#endif

    if (params->resize.method != "Nearest" ) {
        Lanczos (src, dst, dScale);
    } else {
        // Nearest neighbour algorithm
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif

        for (int i = 0; i < dst->getHeight(); i++) {
            int sy = i / dScale;
            sy = LIM (sy, 0, src->getHeight() - 1);

            for (int j = 0; j < dst->getWidth(); j++) {
                int sx = j / dScale;
                sx = LIM (sx, 0, src->getWidth() - 1);
                dst->r (i, j) = src->r (sy, sx);
                dst->g (i, j) = src->g (sy, sx);
                dst->b (i, j) = src->b (sy, sx);
            }
        }
    }

#ifdef PROFILE
    time_t t2 = clock();
    std::cout << "Resize: " << params->resize.method << ": "
              << (float) (t2 - t1) / CLOCKS_PER_SEC << std::endl;
#endif
}

}
