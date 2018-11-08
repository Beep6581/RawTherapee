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
//#define BENCHMARK
#include "StopWatch.h"
#include "sleef.c"

namespace {

#ifdef __SSE2__
void fastlin2log(float *x, float factor, float base, int w)
{
    float baseLog = 1.f / xlogf(base);
    vfloat baseLogv = F2V(baseLog);
    factor = factor * (base - 1.f);
    vfloat factorv = F2V(factor);
    vfloat onev = F2V(1.f);
    int i = 0;
    for (; i < w - 3; i += 4) {
        STVFU(x[i], xlogf(LVFU(x[i]) * factorv + onev) * baseLogv);
    }
    for (; i < w; ++i) {
        x[i] = xlogf(x[i] * factor + 1.f) * baseLog;
    }
}
#endif

}

namespace rtengine {

void ImProcFunctions::labColorCorrectionRegions(LabImage *lab)
{
    if (!params->colorToning.enabled || params->colorToning.method != "LabRegions") {
        return;
    }
BENCHFUN
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

    array2D<float> guide(lab->W, lab->H);

    // magic constant c_factor: normally chromaticity is in [0; 42000] (see color.h), but here we use the constant to match how the chromaticity pipette works (see improcfun.cc lines 4705-4706 and color.cc line 1930
    constexpr float c_factor = 327.68f / 48000.f;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float cBuffer[lab->W];
        float hBuffer[lab->W];
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16)
#endif
        for (int y = 0; y < lab->H; ++y) {
#ifdef __SSE2__
            // vectorized precalculation
            Color::Lab2Lch(lab->a[y], lab->b[y], cBuffer, hBuffer, lab->W);
            fastlin2log(cBuffer, c_factor, 10.f, lab->W);
#endif
            for (int x = 0; x < lab->W; ++x) {
                const float l = lab->L[y][x] / 32768.f;
                guide[y][x] = LIM01(l);
#ifdef __SSE2__
                // use precalculated values
                const float c = cBuffer[x];
                float h = hBuffer[x];
#else
                float c, h;
                Color::Lab2Lch(lab->a[y][x], lab->b[y][x], c, h);
                c = xlin2log(c * c_factor, 10.f);
#endif
                h = Color::huelab_to_huehsv2(h);
                h += 1.f/6.f; // offset the hue because we start from purple instead of red
                if (h > 1.f) {
                    h -= 1.f;
                }
                h = xlin2log(h, 3.f);

                for (int i = begin_idx; i < end_idx; ++i) {
                    auto &hm = hmask[i];
                    auto &cm = cmask[i];
                    auto &lm = lmask[i];
                    float blend = LIM01((hm ? hm->getVal(h) : 1.f) * (cm ? cm->getVal(c) : 1.f) * (lm ? lm->getVal(l) : 1.f));
                    Lmask[i][y][x] = abmask[i][y][x] = blend;
                }
            }
        }
    }

    for (int i = begin_idx; i < end_idx; ++i) {
        rtengine::guidedFilter(guide, abmask[i], abmask[i], max(int(4 / scale + 0.5), 1), 0.001, multiThread);
        rtengine::guidedFilter(guide, Lmask[i], Lmask[i], max(int(25 / scale + 0.5), 1), 0.0001, multiThread);
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
            return 12000.f * SGN(x) * xlog2lin(std::abs(x), 4.f);
        };

    float abca[n];
    float abcb[n];
    float rs[n];
    float rl[n];
    for (int i = 0; i < n; ++i) {
        auto &r = params->colorToning.labregions[i];
        abca[i] = abcoord(r.a);
        abcb[i] = abcoord(r.b);
        rs[i] = 1.f + r.saturation / 100.f;
        rl[i] = 1.f + r.lightness / 500.f;
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        vfloat c42000v = F2V(42000.f);
        vfloat cm42000v = F2V(-42000.f);
        vfloat c32768v = F2V(32768.f);
#endif
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int y = 0; y < lab->H; ++y) {
            int x = 0;
#ifdef __SSE2__
            for (; x < lab->W - 3; x += 4) {
                vfloat lv = LVFU(lab->L[y][x]);
                vfloat av = LVFU(lab->a[y][x]);
                vfloat bv = LVFU(lab->b[y][x]);

                for (int i = 0; i < n; ++i) {
                    vfloat blendv = LVFU(abmask[i][y][x]);
                    vfloat sv = F2V(rs[i]);
                    vfloat a_newv = vclampf(sv * (av + F2V(abca[i])), cm42000v, c42000v);
                    vfloat b_newv = vclampf(sv * (bv + F2V(abcb[i])), cm42000v, c42000v);
                    vfloat l_newv = vclampf(lv * F2V(rl[i]), ZEROV, c32768v);
                    lv = vintpf(LVFU(Lmask[i][y][x]), l_newv, lv);
                    av = vintpf(blendv, a_newv, av);
                    bv = vintpf(blendv, b_newv, bv);
                }
                STVFU(lab->L[y][x], lv);
                STVFU(lab->a[y][x], av);
                STVFU(lab->b[y][x], bv);
            }
#endif
            for (; x < lab->W; ++x) {
                float l = lab->L[y][x];
                float a = lab->a[y][x];
                float b = lab->b[y][x];

                for (int i = 0; i < n; ++i) {
                    float blend = abmask[i][y][x];
                    float s = rs[i];
                    float a_new = LIM(s * (a + abca[i]), -42000.f, 42000.f);
                    float b_new = LIM(s * (b + abcb[i]), -42000.f, 42000.f);
                    float l_new = LIM(l * rl[i], 0.f, 32768.f);
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
}

} // namespace rtengine
