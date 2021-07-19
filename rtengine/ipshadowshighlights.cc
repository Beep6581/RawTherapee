/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "array2D.h"
#include "color.h"
#include "curves.h"
#include "gauss.h"
#include "guidedfilter.h"
#include "iccstore.h"
#include "labimage.h"
#include "opthelper.h"
#include "procparams.h"
#include "sleef.h"

namespace rtengine {
//modifications to pass parameters needs by locallab, to avoid 2 functions - no change in process - J.Desmis march 2019
void ImProcFunctions::shadowsHighlights(LabImage *lab, bool ena, int labmode, int hightli, int shado, int rad, int scal, int hltonal, int shtonal)
{
    if (!ena || (!hightli && !shado)){
        return;
    }
    const int width = lab->W;
    const int height = lab->H;
    const bool lab_mode = labmode;

    array2D<float> mask(width, height);
    array2D<float> L(width, height);
    const float radius = float(rad) * 10 / scal;
    LUTf f(lab_mode ? 32768 : 65536);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);

    const auto rgb2lab =
        [&](float R, float G, float B, float &l, float &a, float &b) -> void
        {
            float x, y, z;
            Color::rgbxyz(R, G, B, x, y, z, ws);
            Color::XYZ2Lab(x, y, z, l, a, b);
        };

    const auto lab2rgb =
        [&](float l, float a, float b, float &R, float &G, float &B) -> void
        {
            float x, y, z;
            Color::Lab2XYZ(l, a, b, x, y, z);
            Color::xyz2rgb(x, y, z, R, G, B, iws);
        };
    
    const auto apply =
        [&](int amount, int tonalwidth, bool hl) -> void
        {
            const float thresh = tonalwidth * 327.68f;
            const float scale = hl ? (thresh > 0.f ? 0.9f / thresh : 1.f) : thresh * 0.9f;

#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float l = lab->L[y][x];
                    float l1 = l / 32768.f;
                    if (hl) {
                        mask[y][x] = (l > thresh) ? 1.f : pow4(l * scale);
                        L[y][x] = 1.f - l1;
                    } else {
                        mask[y][x] = l <= thresh ? 1.f : pow4(scale / l);
                        L[y][x] = l1;
                    }
                }
            }

            guidedFilter(L, mask, mask, radius, 0.075, multiThread, 4);

            const float base = std::pow(4.f, float(amount)/100.f);
            const float gamma = hl ? base : 1.f / base;

            const float contrast = std::pow(2.f, float(amount)/100.f);
            DiagonalCurve sh_contrast({
                    DCT_NURBS,
                    0, 0,
                    0.125, std::pow(0.125 / 0.25, contrast) * 0.25, 
                    0.25, 0.25,
                    0.375, std::pow(0.375 / 0.25, contrast) * 0.25,
                    1, 1
                });

            if(!hl) {
                if (lab_mode) {
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int l = 0; l < 32768; ++l) {
                        auto val = pow_F(l / 32768.f, gamma);
                        // get a bit more contrast in the shadows
                        val = sh_contrast.getVal(val);
                        f[l] = val * 32768.f;
                    }
                } else {
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int c = 0; c < 65536; ++c) {
                        float l, a, b;
                        float R = c, G = c, B = c;
                        rgb2lab(R, G, B, l, a, b);
                        auto val = pow_F(l / 32768.f, gamma);
                        // get a bit more contrast in the shadows
                        val = sh_contrast.getVal(val);
                        l = val * 32768.f;
                        lab2rgb(l, a, b, R, G, B);
                        f[c] = G;
                    }
                }
            } else {
                if (lab_mode) {
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int l = 0; l < 32768; ++l) {
                        auto val = pow_F(l / 32768.f, gamma);
                        f[l] = val * 32768.f;
                    }
                } else {
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int c = 0; c < 65536; ++c) {
                        float l, a, b;
                        float R = c, G = c, B = c;
                        rgb2lab(R, G, B, l, a, b);
                        auto val = pow_F(l / 32768.f, gamma);
                        l = val * 32768.f;
                        lab2rgb(l, a, b, R, G, B);
                        f[c] = G;
                    }
                }
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float l = lab->L[y][x];
                    float blend = LIM01(mask[y][x]);
                    float orig = 1.f - blend;
                    if (l >= 0.f && l < 32768.f) {
                        if (lab_mode) {
                            lab->L[y][x] = intp(blend, f[l], l);
                            if (!hl && l > 1.f) {
                                // when pushing shadows, scale also the chromaticity
                                float s = max(lab->L[y][x] / l * 0.5f, 1.f) * blend;
                                float a = lab->a[y][x];
                                float b = lab->b[y][x];
                                lab->a[y][x] = a * s + a * orig;
                                lab->b[y][x] = b * s + b * orig;
                            }
                        } else {
                            float rgb[3];
                            lab2rgb(l, lab->a[y][x], lab->b[y][x], rgb[0], rgb[1], rgb[2]);
                            for (int i = 0; i < 3; ++i) {
                                float c = rgb[i];
                                if (!OOG(c)) {
                                    rgb[i] = intp(blend, f[c], c);
                                }
                            }
                            rgb2lab(rgb[0], rgb[1], rgb[2], lab->L[y][x], lab->a[y][x], lab->b[y][x]);
                        }
                    }
                }
            }
        };

    if (hightli > 0) {
        apply(hightli * 0.7, hltonal, true);
    }

    if (shado > 0) {
        apply(shado * 0.6, shtonal, false);
    }
}

} // namespace rtengine
