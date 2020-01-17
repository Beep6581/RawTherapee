/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2020 Gabor Horvath <hgabor@rawtherapee.com>
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

#include <memory>
#include "improcfun.h"
#include "labimage.h"
#include "procparams.h"
#include "rt_math.h"

#define BENCHMARK
#include "StopWatch.h"

using namespace std;
namespace {
#ifdef __SSE2__
bool inintervall(float a, float b, float c)
{
    return a < std::max(b, c) && a > std::min(b, c);
}

float selectweight(float a, float b, float low, float high)
{
    const float minVal = std::min(a,b);
    const float maxVal = std::max(a,b);
    const float res = (minVal < 0.45f * maxVal) ? low : high;
    return (minVal > 0.05f * maxVal) ? res : high;
}

#else
bool inintervall(float a, float b, float c)
{
    return (a < b && a > c) || (a < c && a > b);
}

float selectweight(float a, float b, float low, float high)
{
    if ((a < 0.45f * b && a > 0.05f * b) || (b < 0.45f * a && b > 0.05f * a)) {
        return low;
    } else {
        return high;
    }
}

#endif
}
namespace rtengine
{

// To the extent possible under law, Manuel Llorens <manuelllorens@gmail.com>
// has waived all copyright and related or neighboring rights to this work.
// This work is published from: Spain.

// Thanks to Manuel for this excellent job (Jacques Desmis JDC or frej83)
void ImProcFunctions::MLsharpen (LabImage* lab)
{
    // JD: this algorithm maximize clarity of images; it does not play on accutance. It can remove (partially) the effects of the AA filter)
    // I think we can use this algorithm alone in most cases, or first to clarify image and if you want a very little USM (unsharp mask sharpening) after...
    if (!params->sharpenEdge.enabled || params->sharpenEdge.amount == 0) {
        return;
    }
    
BENCHFUN

    const int width = lab->W, height = lab->H;
    constexpr float chmax[3] = {1.f / 8.f, 1.f / 3.f, 1.f / 3.f};
    const int width2 = 2 * width;
    constexpr float eps2 = 0.001f; //prevent divide by zero
    const float amount = params->sharpenEdge.amount / 100.0f;
    const float amountby3 = params->sharpenEdge.amount / 300.0f;

    std::unique_ptr<float[]> L(new float[width * height]);

    const int channels = params->sharpenEdge.threechannels ? 1 : 3;
    const int passes = params->sharpenEdge.passes;

    for (int c = 0; c < channels; ++c) { // c=0 Luminance only
        float** channel = c == 0 ? lab->L : c == 1 ? lab->a : lab->b;
        for (int p = 0; p < passes; ++p) {
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int i = 0; i < height; ++i) {
                for (int j = 0; j < width; ++j) {
                    L[i * width + j] = channel[i][j] / 327.68f;
                }
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int j = 2; j < height - 2; j++) {
                for (int i = 2, offset = j * width + i; i < width - 2; i++, offset++) {
                    // weight functions
                    const float wH = eps2 + fabs(L[offset + 1] - L[offset - 1]);
                    const float wV = eps2 + fabs(L[offset + width] - L[offset - width]);

                    float s = 2.f / (2.f + fabs(wH - wV));
                    float wD1 = eps2 + fabs(L[offset + width + 1] - L[offset - width - 1]) * s;
                    float wD2 = eps2 + fabs(L[offset + width - 1] - L[offset - width + 1]) * s;
                    s = wD1;
                    wD1 /= wD2;
                    wD2 /= s;

                    const float v = L[offset];
                    float lumH, lumV, lumD1, lumD2;
                    lumH = lumV = lumD1 = lumD2 = v;

                    // contrast detection
                    const float contrast = std::min(std::sqrt(SQR(L[offset + 1] - L[offset - 1]) + SQR(L[offset + width] - L[offset - width])) * chmax[c], 1.f);

                    // new possible values
                    if (inintervall(L[offset], L[offset - 1], L[offset + 1])) {
                        float f1 = fabs(L[offset - 2] - L[offset - 1]);
                        float f2 = L[offset - 1] - L[offset];
                        float f3 = (L[offset - 1] - L[offset - width]) * (L[offset - 1] - L[offset + width]);
                        float f4 = std::sqrt(fabs(L[offset - 1] - L[offset - width2]) * fabs(L[offset - 1] - L[offset + width2]));
                        const float difL = f1 * SQR(f2) * SQR(f3) * f4;
                        if (difL > 0.f) {
                            f1 = fabs(L[offset + 2] - L[offset + 1]);
                            f2 = L[offset + 1] - L[offset];
                            f3 = (L[offset + 1] - L[offset - width]) * (L[offset + 1] - L[offset + width]);
                            f4 = std::sqrt(fabs(L[offset + 1] - L[offset - width2]) * fabs(L[offset + 1] - L[offset + width2]));
                            const float difR = f1 * SQR(f2) * SQR(f3) * f4;
                            if (difR > 0.f) {
                                lumH = (L[offset - 1] * difR + L[offset + 1] * difL) / (difL + difR);
                                lumH = intp(contrast, lumH, v);
                            }
                        }
                    }

                    if (inintervall(L[offset], L[offset - width], L[offset + width])) {
                        float f1 = fabs(L[offset - width2] - L[offset - width]);
                        float f2 = L[offset - width] - L[offset];
                        float f3 = (L[offset - width] - L[offset - 1]) * (L[offset - width] - L[offset + 1]);
                        float f4 = std::sqrt(fabs(L[offset - width] - L[offset - 2]) * fabs(L[offset - width] - L[offset + 2]));
                        const float difT = f1 * SQR(f2) * SQR(f3) * f4;
                        if (difT > 0.f) {
                            f1 = fabs(L[offset + width2] - L[offset + width]);
                            f2 = L[offset + width] - L[offset];
                            f3 = (L[offset + width] - L[offset - 1]) * (L[offset + width] - L[offset + 1]);
                            f4 = std::sqrt(fabs(L[offset + width] - L[offset - 2]) * fabs(L[offset + width] - L[offset + 2]));
                            const float difB = f1 * SQR(f2) * SQR(f3) * f4;
                            if (difB > 0.f) {
                                lumV = (L[offset - width] * difB + L[offset + width] * difT) / (difT + difB);
                                lumV = intp(contrast, lumV, v);
                            }
                        }
                    }

                    if (inintervall(L[offset], L[offset - 1 - width], L[offset + 1 + width])) {
                        float f1 = fabs(L[offset - 2 - width2] - L[offset - 1 - width]);
                        float f2 = L[offset - 1 - width] - L[offset];
                        float f3 = (L[offset - 1 - width] - L[offset - width + 1]) * (L[offset - 1 - width] - L[offset + width - 1]);
                        float f4 = std::sqrt(fabs(L[offset - 1 - width] - L[offset - width2 + 2]) * fabs(L[offset - 1 - width] - L[offset + width2 - 2]));
                        const float difLT = f1 * SQR(f2) * SQR(f3) * f4;
                        if (difLT > 0.f) {
                            f1 = fabs(L[offset + 2 + width2] - L[offset + 1 + width]);
                            f2 = L[offset + 1 + width] - L[offset];
                            f3 = (L[offset + 1 + width] - L[offset - width + 1]) * (L[offset + 1 + width] - L[offset + width - 1]);
                            f4 = std::sqrt(fabs(L[offset + 1 + width] - L[offset - width2 + 2]) * fabs(L[offset + 1 + width] - L[offset + width2 - 2]));
                            const float difRB = f1 * SQR(f2) * SQR(f3) * f4;
                            if (difRB > 0.f) {
                                lumD1 = (L[offset - 1 - width] * difRB + L[offset + 1 + width] * difLT) / (difLT + difRB);
                                lumD1 = intp(contrast, lumD1, v);
                            }
                        }
                    }

                    if (inintervall(L[offset], L[offset + 1 - width], L[offset] > L[offset - 1 + width])) {
                        float f1 = fabs(L[offset - 2 + width2] - L[offset - 1 + width]);
                        float f2 = L[offset - 1 + width] - L[offset];
                        float f3 = (L[offset - 1 + width] - L[offset - width - 1]) * (L[offset - 1 + width] - L[offset + width + 1]);
                        float f4 = std::sqrt(fabs(L[offset - 1 + width] - L[offset - width2 - 2]) * fabs(L[offset - 1 + width] - L[offset + width2 + 2]));
                        const float difLB = f1 * SQR(f2) * SQR(f3) * f4;
                        if (difLB > 0.f) {
                            f1 = fabs(L[offset + 2 - width2] - L[offset + 1 - width]);
                            f2 = L[offset + 1 - width] - L[offset];
                            f3 = (L[offset + 1 - width] - L[offset + width + 1]) * (L[offset + 1 - width] - L[offset - width - 1]);
                            f4 = std::sqrt(fabs(L[offset + 1 - width] - L[offset + width2 + 2]) * fabs(L[offset + 1 - width] - L[offset - width2 - 2]));
                            const float difRT = f1 * SQR(f2) * SQR(f3) * f4;
                            if (difRT > 0.f) {
                                lumD2 = (L[offset + 1 - width] * difLB + L[offset - 1 + width] * difRT) / (difLB + difRT);
                                lumD2 = intp(contrast, lumD2, v);
                            }
                        }
                    }

                    // final mix
                    // avoid sharpening diagonals too much
                    const float weight = selectweight(wH, wV, amountby3, amount);

                    if (c == 0) {
                        if (v < 92.f) {
                            channel[j][i] = fabs(327.68f * intp(weight, (lumH * wH + lumV * wV + lumD1 * wD1 + lumD2 * wD2) / (wH + wV + wD1 + wD2), v)); // fabs because lab->L always > 0
                        }
                    } else {
                        channel[j][i] = 327.68f * intp(weight, (lumH * wH + lumV * wV + lumD1 * wD1 + lumD2 * wD2) / (wH + wV + wD1 + wD2), v);
                    }
                }
            }
        }
    }
}

}
