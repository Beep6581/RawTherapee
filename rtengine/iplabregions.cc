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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "improcfun.h"
#include "guidedfilter.h"

namespace rtengine {

void ImProcFunctions::labColorCorrectionRegions(LabImage *lab)
{
    if (!params->colorToning.enabled || params->colorToning.method != "LabRegions") {
        return;
    }

    const float factor = 3.f;
    const float scaling = 1.f;

    int n = params->colorToning.labregions.size();
    int show_mask_idx = params->colorToning.labregionsShowMask;
    if (show_mask_idx >= n) {
        show_mask_idx = -1;
    }
    std::vector<std::unique_ptr<FlatCurve>> hmask(n);
    std::vector<std::unique_ptr<FlatCurve>> cmask(n);
    std::vector<std::unique_ptr<FlatCurve>> lmask(n);

    const int begin_idx = max(show_mask_idx, 0);
    const int end_idx = (show_mask_idx < 0 ? n : show_mask_idx+1);

    for (int i = begin_idx; i < end_idx; ++i) {
        auto &r = params->colorToning.labregions[i];
        if (!r.hueMask.empty() && r.hueMask[0] != FCT_Linear) {
            hmask[i].reset(new FlatCurve(r.hueMask, true));
        }
        if (!r.chromaticityMask.empty() && r.chromaticityMask[0] != FCT_Linear) {
            cmask[i].reset(new FlatCurve(r.chromaticityMask, false));
        }
        if (!r.lightnessMask.empty() && r.lightnessMask[0] != FCT_Linear) {
            lmask[i].reset(new FlatCurve(r.lightnessMask, false));
        }
    }

    std::vector<array2D<float>> abmask(n);
    std::vector<array2D<float>> Lmask(n);
    for (int i = begin_idx; i < end_idx; ++i) {
        abmask[i](lab->W, lab->H);
        Lmask[i](lab->W, lab->H);
    }
    
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < lab->H; ++y) {
        for (int x = 0; x < lab->W; ++x) {
            float l = lab->L[y][x];
            float a = lab->a[y][x];
            float b = lab->b[y][x];
            float c, h;
            Color::Lab2Lch(a, b, c, h);
            float c1 = lin2log(c * (327.68f / 48000.f), 10.f);
            float h1 = Color::huelab_to_huehsv2(h);
            h1 = h1 + 1.f/6.f;
            if (h1 > 1.f) {
                h1 -= 1.f;
            }
            h1 = lin2log(h1, 3.f);
            float l1 = l / 32768.f;

            for (int i = begin_idx; i < end_idx; ++i) {
                auto &hm = hmask[i];
                auto &cm = cmask[i];
                auto &lm = lmask[i];
                float blend = LIM01((hm ? hm->getVal(h1) : 1.f) * (cm ? cm->getVal(c1) : 1.f) * (lm ? lm->getVal(l1) : 1.f));
                Lmask[i][y][x] = abmask[i][y][x] = blend;
            }
        }
    }

    {
        array2D<float> guide(lab->W, lab->H, lab->L);
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < lab->H; ++y) {
            for (int x = 0; x < lab->W; ++x) {
                guide[y][x] = LIM01(lab->L[y][x] / 32768.f);
            }
        }
        
        for (int i = begin_idx; i < end_idx; ++i) {
            rtengine::guidedFilter(guide, abmask[i], abmask[i], max(int(4 / scale + 0.5), 1), 0.001, multiThread);
            rtengine::guidedFilter(guide, Lmask[i], Lmask[i], max(int(25 / scale + 0.5), 1), 0.0001, multiThread);
        }
    }

    if (show_mask_idx >= 0) {
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < lab->H; ++y) {
            for (int x = 0; x < lab->W; ++x) {
                auto blend = abmask[show_mask_idx][y][x];
                lab->a[y][x] = 0.f;
                lab->b[y][x] = blend * 42000.f;
                lab->L[y][x] = LIM(lab->L[y][x] + 32768.f * blend, 0.f, 32768.f);
            }
        }

        return;
    }

    const auto abcoord =
        [](float x) -> float
        {
            const float m = ColorToningParams::LABGRID_CORR_MAX;
            return SGN(x) * log2lin(std::abs(x) / m, 4.f) * m;
        };

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < lab->H; ++y) {
        for (int x = 0; x < lab->W; ++x) {
            float l = lab->L[y][x];
            float a = lab->a[y][x];
            float b = lab->b[y][x];

            for (int i = 0; i < n; ++i) {
                auto &r = params->colorToning.labregions[i];
                float blend = abmask[i][y][x];
                float s = 1.f + r.saturation / 100.f;
                float a_new = LIM(s * (a + 32768.f * abcoord(r.a) / factor / scaling), -42000.f, 42000.f);
                float b_new = LIM(s * (b + 32768.f * abcoord(r.b) / factor / scaling), -42000.f, 42000.f);
                float l_new = LIM(l * (1.f + float(r.lightness) / 1000.f), 0.f, 32768.f);
                l = intp(Lmask[i][y][x], l_new, l);
                a = intp(blend, a_new, a);
                b = intp(blend, b_new, b);
            }
            lab->L[y][x] = l;
            lab->a[y][x] = a;
            lab->b[y][x] = b;
        }
    }
}

} // namespace rtengine
