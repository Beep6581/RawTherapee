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

#include "bilateral2.h"
#include "cieimage.h"
#include "gauss.h"
#include "improcfun.h"
#include "jaggedarray.h"
#include "labimage.h"
#include "opthelper.h"
#include "procparams.h"
#include "rt_algo.h"
#include "rt_math.h"
#include "settings.h"
#include "sleef.h"

//#define BENCHMARK
#include "StopWatch.h"

using namespace std;

namespace {

void sharpenHaloCtrl (float** luminance, float** blurmap, float** base, float** blend, int W, int H, const procparams::SharpeningParams &sharpenParam)
{

    const float scale = (100.f - sharpenParam.halocontrol_amount) * 0.01f;
    const float sharpFac = sharpenParam.amount * 0.01f;
    float** nL = base;

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 2; i < H - 2; i++) {
        float max1 = 0, max2 = 0, min1 = 0, min2 = 0;

        for (int j = 2; j < W - 2; j++) {
            // compute 3 iterations, only forward
            float np1 = 2.f * (nL[i - 2][j] + nL[i - 2][j + 1] + nL[i - 2][j + 2] + nL[i - 1][j] + nL[i - 1][j + 1] + nL[i - 1][j + 2] + nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2]) / 27.f + nL[i - 1][j + 1] / 3.f;
            float np2 = 2.f * (nL[i - 1][j] + nL[i - 1][j + 1] + nL[i - 1][j + 2] + nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2] + nL[i + 1][j] + nL[i + 1][j + 1] + nL[i + 1][j + 2]) / 27.f + nL[i]  [j + 1] / 3.f;
            float np3 = 2.f * (nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2] + nL[i + 1][j] + nL[i + 1][j + 1] + nL[i + 1][j + 2] + nL[i + 2][j] + nL[i + 2][j + 1] + nL[i + 2][j + 2]) / 27.f + nL[i + 1][j + 1] / 3.f;

            // Max/Min of all these deltas and the last two max/min
            float maxn = rtengine::max(np1, np2, np3);
            float minn = rtengine::min(np1, np2, np3);
            float max_ = rtengine::max(max1, max2, maxn);
            float min_ = rtengine::min(min1, min2, minn);

            // Shift the queue
            max1 = max2;
            max2 = maxn;
            min1 = min2;
            min2 = minn;
            float labL = luminance[i][j];

            if (max_ < labL) {
                max_ = labL;
            }

            if (min_ > labL) {
                min_ = labL;
            }

            // deviation from the environment as measurement
            float diff = nL[i][j] - blurmap[i][j];

            constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
            float delta = sharpenParam.threshold.multiply<float, float, float>(
                              rtengine::min(fabsf(diff), upperBound),   // X axis value = absolute value of the difference
                              sharpFac * diff               // Y axis max value = sharpening.amount * signed difference
                          );
            float newL = labL + delta;

            // applying halo control
            if (newL > max_) {
                newL = max_ + (newL - max_) * scale;
            } else if (newL < min_) {
                newL = min_ - (min_ - newL) * scale;
            }

            luminance[i][j] = intp(blend[i][j], newL, luminance[i][j]);
        }
    }
}

void dcdamping (float** aI, float** aO, float damping, int W, int H)
{

    const float dampingFac = -2.0 / (damping * damping);

#ifdef __SSE2__
    vfloat Iv, Ov, Uv, zerov, onev, fourv, fivev, dampingFacv, Tv, Wv, Lv;
    zerov = _mm_setzero_ps();
    onev = F2V(1.f);
    fourv = F2V(4.f);
    fivev = F2V(5.f);
    dampingFacv = F2V(dampingFac);
#endif
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < H; i++) {
        int j = 0;
#ifdef __SSE2__

        for (; j < W - 3; j += 4) {
            Iv = LVFU(aI[i][j]);
            Ov = LVFU(aO[i][j]);
            Lv = xlogf(Iv / Ov);
            Wv = Ov - Iv;
            Uv = (Ov * Lv + Wv) * dampingFacv;
            Uv = vminf(Uv, onev);
            Tv = Uv * Uv;
            Tv = Tv * Tv;
            Uv = Tv * (fivev - Uv * fourv);
            Uv = (Wv / Iv) * Uv + onev;
            Uv = vselfzero(vmaskf_gt(Iv, zerov), Uv);
            Uv = vselfzero(vmaskf_gt(Ov, zerov), Uv);
            STVFU(aI[i][j], Uv);
        }

#endif

        for(; j < W; j++) {
            float I = aI[i][j];
            float O = aO[i][j];

            if (O <= 0.f || I <= 0.f) {
                aI[i][j] = 0.f;
                continue;
            }

            float U = (O * xlogf(I / O) - I + O) * dampingFac;
            U = rtengine::min(U, 1.0f);
            U = U * U * U * U * (5.f - U * 4.f);
            aI[i][j] = (O - I) / I * U + 1.f;
        }
    }
}

}

namespace rtengine
{

void ImProcFunctions::deconvsharpening (float** luminance, float** tmp, const float * const * blend, int W, int H, const procparams::SharpeningParams &sharpenParam, double Scale)
{
    if (sharpenParam.deconvamount == 0 && sharpenParam.blurradius < 0.25f) {
        return;
    }
BENCHFUN
    JaggedArray<float> tmpI(W, H);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < H; i++) {
        for(int j = 0; j < W; j++) {
            tmpI[i][j] = max(luminance[i][j], 0.f);
        }
    }

    JaggedArray<float>* blurbuffer = nullptr;

    if (sharpenParam.blurradius >= 0.25f) {
        blurbuffer = new JaggedArray<float>(W, H);
        JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(tmpI, blur, W, H, sharpenParam.blurradius);
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    blur[i][j] = intp(blend[i][j], luminance[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    }
    const float damping = sharpenParam.deconvdamping / 5.0;
    const bool needdamp = sharpenParam.deconvdamping > 0;
    const double sigma = sharpenParam.deconvradius / Scale;
    const float amount = sharpenParam.deconvamount / 100.f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        for (int k = 0; k < sharpenParam.deconviter; k++) {
            if (!needdamp) {
                // apply gaussian blur and divide luminance by result of gaussian blur
                gaussianBlur(tmpI, tmp, W, H, sigma, false, GAUSS_DIV, luminance);
            } else {
                // apply gaussian blur + damping
                gaussianBlur(tmpI, tmp, W, H, sigma);
                dcdamping(tmp, luminance, damping, W, H);
            }
            gaussianBlur(tmp, tmpI, W, H, sigma, false, GAUSS_MULT);
        } // end for

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                luminance[i][j] = intp(blend[i][j] * amount, max(tmpI[i][j], 0.0f), luminance[i][j]);
            }
        }

        if (sharpenParam.blurradius >= 0.25f) {
            JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    luminance[i][j] = intp(blend[i][j], luminance[i][j], max(blur[i][j], 0.0f));
                }
            }
        }
    } // end parallel
    delete blurbuffer;
}

void ImProcFunctions::sharpening (LabImage* lab, const procparams::SharpeningParams &sharpenParam, bool showMask)
{

    if ((!sharpenParam.enabled) || sharpenParam.amount < 1 || lab->W < 8 || lab->H < 8) {
        return;
    }

    int W = lab->W, H = lab->H;

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    float contrast = sharpenParam.contrast / 100.f;
    buildBlendMask(lab->L, blend, W, H, contrast);

    if(showMask) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                lab->L[i][j] = blend[i][j] * 32768.f;
            }
        }
        return;
    }

    JaggedArray<float> b2(W, H);

    if (sharpenParam.method == "rld") {
        deconvsharpening (lab->L, b2, blend, lab->W, lab->H, sharpenParam, scale);
        return;
    }
BENCHFUN

    // Rest is UNSHARP MASK
    float** b3 = nullptr;

    if (sharpenParam.edgesonly) {
        b3 = new float*[H];

        for (int i = 0; i < H; i++) {
            b3[i] = new float[W];
        }
    }

    JaggedArray<float> blur(W, H);

    if (sharpenParam.blurradius >= 0.25f) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(lab->L, blur, W, H, sharpenParam.blurradius);
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    blur[i][j] = intp(blend[i][j], lab->L[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        if (!sharpenParam.edgesonly) {
            gaussianBlur (lab->L, b2, W, H, sharpenParam.radius / scale);
        } else {
            bilateral<float, float> (lab->L, (float**)b3, b2, W, H, sharpenParam.edges_radius / scale, sharpenParam.edges_tolerance, multiThread);
            gaussianBlur (b3, b2, W, H, sharpenParam.radius / scale);
        }
    }
    float** base = lab->L;

    if (sharpenParam.edgesonly) {
        base = b3;
    }

    if (!sharpenParam.halocontrol) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
                float diff = base[i][j] - b2[i][j];
                float delta = sharpenParam.threshold.multiply<float, float, float>(
                                  min(fabsf(diff), upperBound),                   // X axis value = absolute value of the difference, truncated to the max value of this field
                                  sharpenParam.amount * diff * 0.01f        // Y axis max value
                              );
                lab->L[i][j] = intp(blend[i][j], lab->L[i][j] + delta, lab->L[i][j]);
            }
    } else {
        if (!sharpenParam.edgesonly) {
            // make a deep copy of lab->L
            JaggedArray<float> labCopy(W, H);

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for( int i = 0; i < H; i++ )
                for( int j = 0; j < W; j++ ) {
                    labCopy[i][j] = lab->L[i][j];
                }

            sharpenHaloCtrl (lab->L, b2, labCopy, blend, W, H, sharpenParam);
        } else {
            sharpenHaloCtrl (lab->L, b2, base, blend, W, H, sharpenParam);
        }

    }

    if (sharpenParam.edgesonly) {
        for (int i = 0; i < H; i++) {
            delete [] b3[i];
        }

        delete [] b3;
    }

    if (sharpenParam.blurradius >= 0.25f) {
#ifdef _OPENMP
    #pragma omp parallel for
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                lab->L[i][j] = intp(blend[i][j], lab->L[i][j], max(blur[i][j], 0.0f));
            }
        }
    }

}

// To the extent possible under law, Manuel Llorens <manuelllorens@gmail.com>
// has waived all copyright and related or neighboring rights to this work.
// This work is published from: Spain.

// Thanks to Manuel for this excellent job (Jacques Desmis JDC or frej83)
void ImProcFunctions::MLsharpen (LabImage* lab)
{
    // JD: this algorithm maximize clarity of images; it does not play on accutance. It can remove (partially) the effects of the AA filter)
    // I think we can use this algorithm alone in most cases, or first to clarify image and if you want a very little USM (unsharp mask sharpening) after...
    if (!params->sharpenEdge.enabled) {
        return;
    }

    MyTime t1e, t2e;
    t1e.set();

    int offset, c, i, j, p, width2;
    int width = lab->W, height = lab->H;
    float *L, lumH, lumV, lumD1, lumD2, v, contrast, s;
    float difL, difR, difT, difB, difLT, difRB, difLB, difRT, wH, wV, wD1, wD2, chmax[3];
    float f1, f2, f3, f4;
    float templab;
    int iii, kkk;
    width2 = 2 * width;
    const float epsil = 0.01f; //prevent divide by zero
    const float eps2 = 0.001f; //prevent divide by zero
    float amount;
    amount = params->sharpenEdge.amount / 100.0f;

    if (amount < 0.00001f) {
        return;
    }

    if (settings->verbose) {
        printf ("SharpenEdge amount %f\n", amount);
    }

    L = new float[width * height];

    chmax[0] = 8.0f;
    chmax[1] = 3.0f;
    chmax[2] = 3.0f;

    int channels;

    if (params->sharpenEdge.threechannels) {
        channels = 0;
    } else {
        channels = 2;
    }

    if (settings->verbose) {
        printf ("SharpenEdge channels %d\n", channels);
    }

    int passes = params->sharpenEdge.passes;

    if (settings->verbose) {
        printf ("SharpenEdge passes %d\n", passes);
    }

    for (p = 0; p < passes; p++)
        for (c = 0; c <= channels; c++) { // c=0 Luminance only

#ifdef _OPENMP
            #pragma omp parallel for private(offset) shared(L)
#endif

            for (offset = 0; offset < width * height; offset++) {
                int ii = offset / width;
                int kk = offset - ii * width;

                if      (c == 0) {
                    L[offset] = lab->L[ii][kk] / 327.68f;    // adjust to RT and to 0..100
                } else if (c == 1) {
                    L[offset] = lab->a[ii][kk] / 327.68f;
                } else { /*if (c==2) */
                    L[offset] = lab->b[ii][kk] / 327.68f;
                }
            }

#ifdef _OPENMP
            #pragma omp parallel for private(j,i,iii,kkk, templab,offset,wH,wV,wD1,wD2,s,lumH,lumV,lumD1,lumD2,v,contrast,f1,f2,f3,f4,difT,difB,difL,difR,difLT,difLB,difRT,difRB) shared(lab,L,amount)
#endif

            for(j = 2; j < height - 2; j++)
                for(i = 2, offset = j * width + i; i < width - 2; i++, offset++) {
                    // weight functions
                    wH = eps2 + fabs(L[offset + 1] - L[offset - 1]);
                    wV = eps2 + fabs(L[offset + width] - L[offset - width]);

                    s = 1.0f + fabs(wH - wV) / 2.0f;
                    wD1 = eps2 + fabs(L[offset + width + 1] - L[offset - width - 1]) / s;
                    wD2 = eps2 + fabs(L[offset + width - 1] - L[offset - width + 1]) / s;
                    s = wD1;
                    wD1 /= wD2;
                    wD2 /= s;

                    // initial values
                    int ii = offset / width;
                    int kk = offset - ii * width;

                    if      (c == 0) {
                        lumH = lumV = lumD1 = lumD2 = v = lab->L[ii][kk] / 327.68f;
                    } else if (c == 1) {
                        lumH = lumV = lumD1 = lumD2 = v = lab->a[ii][kk] / 327.68f;
                    } else { /* if (c==2) */
                        lumH = lumV = lumD1 = lumD2 = v = lab->b[ii][kk] / 327.68f;
                    }


                    // contrast detection
                    contrast = sqrt(fabs(L[offset + 1] - L[offset - 1]) * fabs(L[offset + 1] - L[offset - 1]) + fabs(L[offset + width] - L[offset - width]) * fabs(L[offset + width] - L[offset - width])) / chmax[c];

                    if (contrast > 1.0f) {
                        contrast = 1.0f;
                    }

                    // new possible values
                    if (((L[offset] < L[offset - 1]) && (L[offset] > L[offset + 1])) || ((L[offset] > L[offset - 1]) && (L[offset] < L[offset + 1]))) {
                        f1 = fabs(L[offset - 2] - L[offset - 1]);
                        f2 = fabs(L[offset - 1] - L[offset]);
                        f3 = fabs(L[offset - 1] - L[offset - width]) * fabs(L[offset - 1] - L[offset + width]);
                        f4 = sqrt(fabs(L[offset - 1] - L[offset - width2]) * fabs(L[offset - 1] - L[offset + width2]));
                        difL = f1 * f2 * f2 * f3 * f3 * f4;
                        f1 = fabs(L[offset + 2] - L[offset + 1]);
                        f2 = fabs(L[offset + 1] - L[offset]);
                        f3 = fabs(L[offset + 1] - L[offset - width]) * fabs(L[offset + 1] - L[offset + width]);
                        f4 = sqrt(fabs(L[offset + 1] - L[offset - width2]) * fabs(L[offset + 1] - L[offset + width2]));
                        difR = f1 * f2 * f2 * f3 * f3 * f4;

                        if ((difR > epsil) && (difL > epsil)) {
                            lumH = (L[offset - 1] * difR + L[offset + 1] * difL) / (difL + difR);
                            lumH = v * (1.f - contrast) + lumH * contrast;
                        }
                    }

                    if (((L[offset] < L[offset - width]) && (L[offset] > L[offset + width])) || ((L[offset] > L[offset - width]) && (L[offset] < L[offset + width]))) {
                        f1 = fabs(L[offset - width2] - L[offset - width]);
                        f2 = fabs(L[offset - width] - L[offset]);
                        f3 = fabs(L[offset - width] - L[offset - 1]) * fabs(L[offset - width] - L[offset + 1]);
                        f4 = sqrt(fabs(L[offset - width] - L[offset - 2]) * fabs(L[offset - width] - L[offset + 2]));
                        difT = f1 * f2 * f2 * f3 * f3 * f4;
                        f1 = fabs(L[offset + width2] - L[offset + width]);
                        f2 = fabs(L[offset + width] - L[offset]);
                        f3 = fabs(L[offset + width] - L[offset - 1]) * fabs(L[offset + width] - L[offset + 1]);
                        f4 = sqrt(fabs(L[offset + width] - L[offset - 2]) * fabs(L[offset + width] - L[offset + 2]));
                        difB = f1 * f2 * f2 * f3 * f3 * f4;

                        if ((difB > epsil) && (difT > epsil)) {
                            lumV = (L[offset - width] * difB + L[offset + width] * difT) / (difT + difB);
                            lumV = v * (1.f - contrast) + lumV * contrast;
                        }
                    }

                    if (((L[offset] < L[offset - 1 - width]) && (L[offset] > L[offset + 1 + width])) || ((L[offset] > L[offset - 1 - width]) && (L[offset] < L[offset + 1 + width]))) {
                        f1 = fabs(L[offset - 2 - width2] - L[offset - 1 - width]);
                        f2 = fabs(L[offset - 1 - width] - L[offset]);
                        f3 = fabs(L[offset - 1 - width] - L[offset - width + 1]) * fabs(L[offset - 1 - width] - L[offset + width - 1]);
                        f4 = sqrt(fabs(L[offset - 1 - width] - L[offset - width2 + 2]) * fabs(L[offset - 1 - width] - L[offset + width2 - 2]));
                        difLT = f1 * f2 * f2 * f3 * f3 * f4;
                        f1 = fabs(L[offset + 2 + width2] - L[offset + 1 + width]);
                        f2 = fabs(L[offset + 1 + width] - L[offset]);
                        f3 = fabs(L[offset + 1 + width] - L[offset - width + 1]) * fabs(L[offset + 1 + width] - L[offset + width - 1]);
                        f4 = sqrt(fabs(L[offset + 1 + width] - L[offset - width2 + 2]) * fabs(L[offset + 1 + width] - L[offset + width2 - 2]));
                        difRB = f1 * f2 * f2 * f3 * f3 * f4;

                        if ((difLT > epsil) && (difRB > epsil)) {
                            lumD1 = (L[offset - 1 - width] * difRB + L[offset + 1 + width] * difLT) / (difLT + difRB);
                            lumD1 = v * (1.f - contrast) + lumD1 * contrast;
                        }
                    }

                    if (((L[offset] < L[offset + 1 - width]) && (L[offset] > L[offset - 1 + width])) || ((L[offset] > L[offset + 1 - width]) && (L[offset] < L[offset - 1 + width]))) {
                        f1 = fabs(L[offset - 2 + width2] - L[offset - 1 + width]);
                        f2 = fabs(L[offset - 1 + width] - L[offset]);
                        f3 = fabs(L[offset - 1 + width] - L[offset - width - 1]) * fabs(L[offset - 1 + width] - L[offset + width + 1]);
                        f4 = sqrt(fabs(L[offset - 1 + width] - L[offset - width2 - 2]) * fabs(L[offset - 1 + width] - L[offset + width2 + 2]));
                        difLB = f1 * f2 * f2 * f3 * f3 * f4;
                        f1 = fabs(L[offset + 2 - width2] - L[offset + 1 - width]);
                        f2 = fabs(L[offset + 1 - width] - L[offset]) * fabs(L[offset + 1 - width] - L[offset]);
                        f3 = fabs(L[offset + 1 - width] - L[offset + width + 1]) * fabs(L[offset + 1 - width] - L[offset - width - 1]);
                        f4 = sqrt(fabs(L[offset + 1 - width] - L[offset + width2 + 2]) * fabs(L[offset + 1 - width] - L[offset - width2 - 2]));
                        difRT = f1 * f2 * f2 * f3 * f3 * f4;

                        if ((difLB > epsil) && (difRT > epsil)) {
                            lumD2 = (L[offset + 1 - width] * difLB + L[offset - 1 + width] * difRT) / (difLB + difRT);
                            lumD2 = v * (1.f - contrast) + lumD2 * contrast;
                        }
                    }

                    s = amount;

                    // avoid sharpening diagonals too much
                    if (((fabs(wH / wV) < 0.45f) && (fabs(wH / wV) > 0.05f)) || ((fabs(wV / wH) < 0.45f) && (fabs(wV / wH) > 0.05f))) {
                        s = amount / 3.0f;
                    }

                    // final mix
                    if ((wH != 0.0f) && (wV != 0.0f) && (wD1 != 0.0f) && (wD2 != 0.0f)) {
                        iii = offset / width;
                        kkk = offset - iii * width;
                        float provL = lab->L[iii][kkk] / 327.68f;

                        if(c == 0) {
                            if(provL < 92.f) {
                                templab = v * (1.f - s) + (lumH * wH + lumV * wV + lumD1 * wD1 + lumD2 * wD2) / (wH + wV + wD1 + wD2) * s;
                            } else {
                                templab = provL;
                            }
                        } else {
                            templab = v * (1.f - s) + (lumH * wH + lumV * wV + lumD1 * wD1 + lumD2 * wD2) / (wH + wV + wD1 + wD2) * s;
                        }

                        if      (c == 0) {
                            lab->L[iii][kkk] = fabs(327.68f * templab);    // fabs because lab->L always >0
                        } else if (c == 1) {
                            lab->a[iii][kkk] =      327.68f * templab ;
                        } else if (c == 2) {
                            lab->b[iii][kkk] =      327.68f * templab ;
                        }
                    }

                }
        }

    delete [] L;

    t2e.set();

    if (settings->verbose) {
        printf("SharpenEdge gradient  %d usec\n", t2e.etime(t1e));
    }
}

// To the extent possible under law, Manuel Llorens <manuelllorens@gmail.com>
// has waived all copyright and related or neighboring rights to this work.
// This code is licensed under CC0 v1.0, see license information at
// http://creativecommons.org/publicdomain/zero/1.0/

//! MicroContrast is a sharpening method developed by Manuel Llorens and documented here: http://www.rawness.es/sharpening/?lang=en
//! <BR>The purpose is maximize clarity of the image without creating halo's.
//! <BR>Addition from JD : pyramid  + pondered contrast with matrix 5x5
//! <BR>2017 Ingo Weyrich : reduced processing time
//! \param luminance : Luminance channel of image
void ImProcFunctions::MLmicrocontrast(float** luminance, int W, int H)
{
    if (!params->sharpenMicro.enabled || params->sharpenMicro.contrast == 100 || params->sharpenMicro.amount < 1.0) {
        return;
    }
BENCHFUN
    const int k = params->sharpenMicro.matrix ? 1 : 2;
    // k=2 matrix 5x5  k=1 matrix 3x3
    const int width = W, height = H;
    const int unif = params->sharpenMicro.uniformity;
    const float amount = (k == 1 ? 2.7f : 1.f) * params->sharpenMicro.amount / 1500.0f; //amount 2000.0 quasi no artifacts ==> 1500 = maximum, after artifacts, 25/9 if 3x3

    if (settings->verbose) {
        printf ("Micro-contrast amount %f\n", amount);
        printf ("Micro-contrast uniformity %i\n", unif);
    }

    //modulation uniformity in function of luminance
    const float L98[11] = {0.001f, 0.0015f, 0.002f, 0.004f, 0.006f, 0.008f, 0.01f, 0.03f, 0.05f, 0.1f, 0.1f};
    const float L95[11] = {0.0012f, 0.002f, 0.005f, 0.01f, 0.02f, 0.05f, 0.1f, 0.12f, 0.15f, 0.2f, 0.25f};
    const float L92[11] = {0.01f, 0.015f, 0.02f, 0.06f, 0.10f, 0.13f, 0.17f, 0.25f, 0.3f, 0.32f, 0.35f};
    const float L90[11] = {0.015f, 0.02f, 0.04f, 0.08f, 0.12f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f};
    const float L87[11] = {0.025f, 0.03f, 0.05f, 0.1f, 0.15f, 0.25f, 0.3f, 0.4f, 0.5f, 0.63f, 0.75f};
    const float L83[11] = {0.055f, 0.08f, 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.75f, 0.85f};
    const float L80[11] = {0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f};
    const float L75[11] = {0.22f, 0.25f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.85f, 0.9f, 0.95f};
    const float L70[11] = {0.35f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.97f, 1.0f, 1.0f, 1.0f, 1.0f};
    const float L63[11] = {0.55f, 0.6f, 0.7f, 0.8f, 0.85f, 0.9f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    const float L58[11] = {0.75f, 0.77f, 0.8f, 0.9f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    //default 5
    //modulation contrast
    const float Cont0[11] = {0.05f, 0.1f, 0.2f, 0.25f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f};
    const float Cont1[11] = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 0.95f, 1.0f};
    const float Cont2[11] = {0.2f, 0.40f, 0.6f, 0.7f, 0.8f, 0.85f, 0.90f, 0.95f, 1.0f, 1.05f, 1.10f};
    const float Cont3[11] = {0.5f, 0.6f, 0.7f, 0.8f, 0.85f, 0.9f, 1.0f, 1.0f, 1.05f, 1.10f, 1.20f};
    const float Cont4[11] = {0.8f, 0.85f, 0.9f, 0.95f, 1.0f, 1.05f, 1.10f, 1.150f, 1.2f, 1.25f, 1.40f};
    const float Cont5[11] = {1.0f, 1.1f, 1.2f, 1.25f, 1.3f, 1.4f, 1.45f, 1.50f, 1.6f, 1.65f, 1.80f};

    const float sqrt2 = sqrt(2.0);
    const float sqrt1d25 = sqrt(1.25);
    float *LM = new float[width * height]; //allocation for Luminance

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    float contrast = params->sharpenMicro.contrast / 100.f;
    buildBlendMask(luminance, blend, W, H, contrast);

#ifdef _OPENMP
    #pragma omp parallel
#endif
{

#ifdef _OPENMP
    #pragma omp for schedule(dynamic,16)
#endif

    for(int j = 0; j < height; j++)
        for(int i = 0, offset = j * width + i; i < width; i++, offset++) {
            LM[offset] = luminance[j][i] / 327.68f; // adjust to [0;100] and to RT variables
        }

#ifdef _OPENMP
    #pragma omp for schedule(dynamic,16)
#endif

    for(int j = k; j < height - k; j++)
        for(int i = k, offset = j * width + i; i < width - k; i++, offset++) {
            float v = LM[offset];

            float contrast;
            if (k == 1) {
                contrast = sqrtf(SQR(LM[offset + 1] - LM[offset - 1]) + SQR(LM[offset + width] - LM[offset - width])) * 0.125f;    //for 3x3
            } else /* if (k==2) */ contrast = sqrtf(SQR(LM[offset + 1] - LM[offset - 1]) + SQR(LM[offset + width] - LM[offset - width])
                                                       + SQR(LM[offset + 2] - LM[offset - 2]) + SQR(LM[offset + 2 * width] - LM[offset - 2 * width])) * 0.0625f; //for 5x5

            contrast = std::min(contrast, 1.f);

            //matrix 5x5
            float temp = v + 4.f *( v * (amount + sqrt2 * amount)); //begin 3x3
            float temp1 = sqrt2 * amount *(LM[offset - width - 1] + LM[offset - width + 1] + LM[offset + width - 1] + LM[offset + width + 1]);
            temp1 += amount * (LM[offset - width] + LM[offset - 1] + LM[offset + 1] + LM[offset + width]);

            temp -= temp1;

            // add JD continue 5x5
            if (k == 2) {
                float temp2 = -(LM[offset + 2 * width] + LM[offset - 2 * width] + LM[offset - 2] + LM[offset + 2]);

                temp2 -= sqrt1d25 * (LM[offset + 2 * width - 1] + LM[offset + 2 * width + 1] + LM[offset +  width + 2] + LM[offset +  width - 2] + 
                                     LM[offset - 2 * width - 1] + LM[offset - 2 * width + 1] + LM[offset -  width + 2] + LM[offset -  width - 2]);

                temp2 -= sqrt2 * (LM[offset + 2 * width - 2] + LM[offset + 2 * width + 2] + LM[offset - 2 * width - 2] + LM[offset - 2 * width + 2]);
                temp2 += 18.601126159f * v ; // 18.601126159 = 4 + 4 * sqrt(2) + 8 * sqrt(1.25)
                temp2 *= 2.f * amount;
                temp += temp2;
            }

            temp = std::max(temp, 0.f);

            for(int row = j - k; row <= j + k; ++row) {
                for(int offset2 = row * width + i - k; offset2 <= row * width + i + k; ++offset2) {
                    if((LM[offset2] - temp) * (v - LM[offset2]) > 0.f) {
                        temp = intp(0.75f, temp, LM[offset2]);
                        goto breakout;
                    }
                }
            }
            breakout:

            if (LM[offset] > 95.0f || LM[offset] < 5.0f) {
                contrast *= Cont0[unif];    //+ JD : luminance  pyramid to adjust contrast by evaluation of LM[offset]
            } else if (LM[offset] > 90.0f || LM[offset] < 10.0f) {
                contrast *= Cont1[unif];
            } else if (LM[offset] > 80.0f || LM[offset] < 20.0f) {
                contrast *= Cont2[unif];
            } else if (LM[offset] > 70.0f || LM[offset] < 30.0f) {
                contrast *= Cont3[unif];
            } else if (LM[offset] > 60.0f || LM[offset] < 40.0f) {
                contrast *= Cont4[unif];
            } else {
                contrast *= Cont5[unif];    //(2.0f/k)*Cont5[unif];
            }

            contrast = std::min(contrast, 1.f);

            float tempL = intp(contrast, LM[offset], temp);
            // JD: modulation of microcontrast in function of original Luminance and modulation of luminance
            if (tempL > LM[offset]) {
                float temp2 = tempL / LM[offset]; //for highlights
                temp2 = std::min(temp2, 1.7f); //limit action
                temp2 -= 1.f;
                if (LM[offset] > 98.0f) {
                    temp = 0.f;
                } else if (LM[offset] > 95.0f) {
                    temp = L95[unif];
                } else if (LM[offset] > 92.0f) {
                    temp = L92[unif];
                } else if (LM[offset] > 90.0f) {
                    temp = L90[unif];
                } else if (LM[offset] > 87.0f) {
                    temp = L87[unif];
                } else if (LM[offset] > 83.0f) {
                    temp = L83[unif];
                } else if (LM[offset] > 80.0f) {
                    temp = L80[unif];
                } else if (LM[offset] > 75.0f) {
                    temp = L75[unif];
                } else if (LM[offset] > 70.0f) {
                    temp = L70[unif];
                } else if (LM[offset] > 63.0f) {
                    temp = L63[unif];
                } else if (LM[offset] > 58.0f) {
                    temp = L58[unif];
                } else if (LM[offset] > 42.0f) {
                    temp = L58[unif];
                } else if (LM[offset] > 37.0f) {
                    temp = L63[unif];
                } else if (LM[offset] > 30.0f) {
                    temp = L70[unif];
                } else if (LM[offset] > 25.0f) {
                    temp = L75[unif];
                } else if (LM[offset] > 20.0f) {
                    temp = L80[unif];
                } else if (LM[offset] > 17.0f) {
                    temp = L83[unif];
                } else if (LM[offset] > 13.0f) {
                    temp = L87[unif];
                } else if (LM[offset] > 10.0f) {
                    temp = L90[unif];
                } else if (LM[offset] > 5.0f) {
                    temp = L95[unif];
                } else {
                    temp = 0.f;
                }
                luminance[j][i] = intp(blend[j][i], luminance[j][i] * (temp * temp2 + 1.f), luminance[j][i]);
            } else {

                float temp4 = LM[offset] / tempL; //

                if (temp4 > 1.0f) {
                    temp4 = std::min(temp4, 1.7f); //limit action
                    temp4 -= 1.f;
                    if (LM[offset] < 2.0f) {
                        temp = L98[unif];
                    } else if (LM[offset] < 5.0f) {
                        temp = L95[unif];
                    } else if (LM[offset] < 8.0f) {
                        temp = L92[unif];
                    } else if (LM[offset] < 10.0f) {
                        temp = L90[unif];
                    } else if (LM[offset] < 13.0f) {
                        temp = L87[unif];
                    } else if (LM[offset] < 17.0f) {
                        temp = L83[unif];
                    } else if (LM[offset] < 20.0f) {
                        temp = L80[unif];
                    } else if (LM[offset] < 25.0f) {
                        temp = L75[unif];
                    } else if (LM[offset] < 30.0f) {
                        temp = L70[unif];
                    } else if (LM[offset] < 37.0f) {
                        temp = L63[unif];
                    } else if (LM[offset] < 42.0f) {
                        temp = L58[unif];
                    } else if (LM[offset] < 58.0f) {
                        temp = L58[unif];
                    } else if (LM[offset] < 63.0f) {
                        temp = L63[unif];
                    } else if (LM[offset] < 70.0f) {
                        temp = L70[unif];
                    } else if (LM[offset] < 75.0f) {
                        temp = L75[unif];
                    } else if (LM[offset] < 80.0f) {
                        temp = L80[unif];
                    } else if (LM[offset] < 83.0f) {
                        temp = L83[unif];
                    } else if (LM[offset] < 87.0f) {
                        temp = L87[unif];
                    } else if (LM[offset] < 90.0f) {
                        temp = L90[unif];
                    } else if (LM[offset] < 95.0f) {
                        temp = L95[unif];
                    } else {
                        temp = 0.f;
                    }
                    luminance[j][i] = intp(blend[j][i], luminance[j][i] / (temp * temp4 + 1.f), luminance[j][i]);
                }
            }
        }
}
    delete [] LM;
}

void ImProcFunctions::MLmicrocontrast(LabImage* lab)
{
    MLmicrocontrast(lab->L, lab->W, lab->H);
}

void ImProcFunctions::MLmicrocontrastcam(CieImage* ncie)
{
    MLmicrocontrast(ncie->sh_p, ncie->W, ncie->H);
}

void ImProcFunctions::sharpeningcam (CieImage* ncie, float** b2, bool showMask)
{
    if ((!params->sharpening.enabled) || params->sharpening.amount < 1 || ncie->W < 8 || ncie->H < 8) {
        return;
    }

    int W = ncie->W, H = ncie->H;

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    float contrast = params->sharpening.contrast / 100.f;
    buildBlendMask(ncie->sh_p, blend, W, H, contrast);
    if(showMask) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                ncie->sh_p[i][j] = blend[i][j] * 32768.f;
            }
        }
        return;
    }

    if (params->sharpening.method == "rld") {
        deconvsharpening (ncie->sh_p, b2, blend, ncie->W, ncie->H, params->sharpening, scale);
        return;
    }

    // Rest is UNSHARP MASK

    float** b3 = nullptr;

    if (params->sharpening.edgesonly) {
        b3 = new float*[H];

        for (int i = 0; i < H; i++) {
            b3[i] = new float[W];
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        if (!params->sharpening.edgesonly) {
            gaussianBlur (ncie->sh_p, b2, W, H, params->sharpening.radius / scale);
        } else {
            bilateral<float, float> (ncie->sh_p, (float**)b3, b2, W, H, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance, multiThread);
            gaussianBlur (b3, b2, W, H, params->sharpening.radius / scale);
        }
    }

    float** base = ncie->sh_p;

    if (params->sharpening.edgesonly) {
        base = b3;
    }

    if (!params->sharpening.halocontrol) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
                float diff = base[i][j] - b2[i][j];
                float delta = params->sharpening.threshold.multiply<float, float, float>(
                                  min(fabsf(diff), upperBound),                   // X axis value = absolute value of the difference, truncated to the max value of this field
                                  params->sharpening.amount * diff * 0.01f      // Y axis max value
                              );

                if(ncie->J_p[i][j] > 8.0f && ncie->J_p[i][j] < 92.0f) {
                    ncie->sh_p[i][j] = intp(blend[i][j], ncie->sh_p[i][j] + delta, ncie->sh_p[i][j]);
                }
            }
    } else {
        float** ncieCopy = nullptr;

        if (!params->sharpening.edgesonly) {
            // make deep copy of ncie->sh_p
            ncieCopy = new float*[H];

            for( int i = 0; i < H; i++ ) {
                ncieCopy[i] = new float[W];
            }

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for( int i = 0; i < H; i++ )
                for( int j = 0; j < W; j++ ) {
                    ncieCopy[i][j] = ncie->sh_p[i][j];
                }

            base = ncieCopy;
        }

        sharpenHaloCtrl (ncie->sh_p, b2, base, blend, W, H, params->sharpening);

        if(ncieCopy) {
            for( int i = 0; i < H; i++ ) {
                delete[] ncieCopy[i];
            }

            delete[] ncieCopy;
        }
    }

    if (params->sharpening.edgesonly) {
        for (int i = 0; i < H; i++) {
            delete [] b3[i];
        }

        delete [] b3;
    }
}

}
