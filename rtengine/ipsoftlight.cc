/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
  * Optimized 2019 Ingo Weyrich <heckflosse67@gmx.de>
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

#include "color.h"
#include "iccstore.h"
#include "improcfun.h"
#include "labimage.h"

#include "procparams.h"

namespace rtengine
{

namespace
{

inline float sl(float blend, float x)
{
    if (!rtengine::OOG(x)) {
        float v = rtengine::Color::gamma_srgb(x) / rtengine::MAXVALF;
        // using Pegtop's formula from
        // https://en.wikipedia.org/wiki/Blend_modes#Soft_Light
        // const float orig = 1.f - blend;
        // float v2 = v * v;
        // float v22 = v2 * 2.f;
        // v = v2 + v22 - v22 * v;
        // return blend * Color::igamma_srgb(v * MAXVALF) + orig * x;

        // using optimized formula (heckflosse67@gmx.de)
        return rtengine::intp(blend, rtengine::Color::igamma_srgb(v * v * (3.f - 2.f * v) * rtengine::MAXVALF), x);
    }

    return x;
}
} // namespace

#ifdef __SSE2__
inline vfloat sl(vfloat blend, vfloat x)
{
    const vfloat v = rtengine::Color::gammatab_srgb[x] / F2V(rtengine::MAXVALF);
    return vself(vmaskf_gt(x, F2V(rtengine::MAXVALF)), x, vself(vmaskf_lt(x, ZEROV), x, vintpf(blend, rtengine::Color::igammatab_srgb[v * v * (F2V(3.f) - (v + v)) * rtengine::MAXVALF], x)));
}
#endif

//} // namespace

void ImProcFunctions::softLight(LabImage *lab, const rtengine::procparams::SoftLightParams &softLightParams)
{
    if (!softLightParams.enabled || !softLightParams.strength) {
        return;
    }

    const TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    const float wp[3][3] = {
        {static_cast<float>(wprof[0][0]), static_cast<float>(wprof[0][1]), static_cast<float>(wprof[0][2])},
        {static_cast<float>(wprof[1][0]), static_cast<float>(wprof[1][1]), static_cast<float>(wprof[1][2])},
        {static_cast<float>(wprof[2][0]), static_cast<float>(wprof[2][1]), static_cast<float>(wprof[2][2])}
    };

    const TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);
    const float wip[3][3] = {
        {static_cast<float>(wiprof[0][0]), static_cast<float>(wiprof[0][1]), static_cast<float>(wiprof[0][2])},
        {static_cast<float>(wiprof[1][0]), static_cast<float>(wiprof[1][1]), static_cast<float>(wiprof[1][2])},
        {static_cast<float>(wiprof[2][0]), static_cast<float>(wiprof[2][1]), static_cast<float>(wiprof[2][2])}
    };

#ifdef __SSE2__
    const vfloat wpv[3][3] = {
        {F2V(wprof[0][0]), F2V(wprof[0][1]), F2V(wprof[0][2])},
        {F2V(wprof[1][0]), F2V(wprof[1][1]), F2V(wprof[1][2])},
        {F2V(wprof[2][0]), F2V(wprof[2][1]), F2V(wprof[2][2])}
    };

    const vfloat wipv[3][3] = {
        {F2V(wiprof[0][0]), F2V(wiprof[0][1]), F2V(wiprof[0][2])},
        {F2V(wiprof[1][0]), F2V(wiprof[1][1]), F2V(wiprof[1][2])},
        {F2V(wiprof[2][0]), F2V(wiprof[2][1]), F2V(wiprof[2][2])}
    };
#endif

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        const float blend = softLightParams.strength / 100.f;
#ifdef __SSE2__
        const vfloat blendv = F2V(blend);
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int i = 0; i < lab->H; ++i) {
            int j = 0;
#ifdef __SSE2__

            for (; j < lab->W - 3; j += 4) {
                vfloat Xv, Yv, Zv;
                vfloat Rv, Gv, Bv;
                Color::Lab2XYZ(LVFU(lab->L[i][j]), LVFU(lab->a[i][j]), LVFU(lab->b[i][j]), Xv, Yv, Zv);
                Color::xyz2rgb(Xv, Yv, Zv, Rv, Gv, Bv, wipv);
                Rv = sl(blendv, Rv);
                Gv = sl(blendv, Gv);
                Bv = sl(blendv, Bv);
                Color::rgbxyz(Rv, Gv, Bv, Xv, Yv, Zv, wpv);

                for (int k = 0; k < 4; ++k) {
                    Color::XYZ2Lab(Xv[k], Yv[k], Zv[k], lab->L[i][j + k], lab->a[i][j + k], lab->b[i][j + k]);
                }
            }

#endif

            for (; j < lab->W; j++) {
                float X, Y, Z;
                float R, G, B;
                Color::Lab2XYZ(lab->L[i][j], lab->a[i][j], lab->b[i][j], X, Y, Z);
                Color::xyz2rgb(X, Y, Z, R, G, B, wip);
                R = sl(blend, R);
                G = sl(blend, G);
                B = sl(blend, B);
                Color::rgbxyz(R, G, B, X, Y, Z, wp);
                Color::XYZ2Lab(X, Y, Z, lab->L[i][j], lab->a[i][j], lab->b[i][j]);
            }
        }
    }
}
}