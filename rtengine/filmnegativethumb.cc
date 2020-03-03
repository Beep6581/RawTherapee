/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Romei <aldrop8@gmail.com>
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
#include <cmath>

#include "colortemp.h"
#include "LUT.h"
#include "rtengine.h"
#include "rtthumbnail.h"
#include "opthelper.h"
#include "sleef.h"
#include "rt_algo.h"
#include "settings.h"
#include "procparams.h"
#define BENCHMARK
#include "StopWatch.h"

namespace
{

using rtengine::Imagefloat;
using rtengine::findMinMaxPercentile;


void calcMedians(
    const Imagefloat* baseImg,
    const int x1, const int y1,
    const int x2, const int y2,
    float &rmed, float &gmed, float &bmed
)
{
    // Channel vectors to calculate medians
    std::vector<float> rv, gv, bv;

    const int sz = std::max(0, (y2 - y1) * (x2 - x1));
    rv.reserve(sz);
    gv.reserve(sz);
    bv.reserve(sz);

    for (int i = y1; i < y2; i++) {
        for (int j = x1; j < x2; j++) {
            rv.push_back(baseImg->r(i, j));
            gv.push_back(baseImg->g(i, j));
            bv.push_back(baseImg->b(i, j));
        }
    }

    // Calculate channel medians from whole image
    findMinMaxPercentile(rv.data(), rv.size(), 0.5f, rmed, 0.5f, rmed, true);
    findMinMaxPercentile(gv.data(), gv.size(), 0.5f, gmed, 0.5f, gmed, true);
    findMinMaxPercentile(bv.data(), bv.size(), 0.5f, bmed, 0.5f, bmed, true);
}

}

void rtengine::Thumbnail::processFilmNegative(
    const procparams::ProcParams &params,
    const Imagefloat* baseImg,
    const int rwidth, const int rheight
)
{

    // Channel exponents
    const float rexp = -params.filmNegative.redRatio * params.filmNegative.greenExp;
    const float gexp = -params.filmNegative.greenExp;
    const float bexp = -params.filmNegative.blueRatio * params.filmNegative.greenExp;

    // Calculate output multipliers
    float rmult, gmult, bmult;

    const float MAX_OUT_VALUE = 65000.f;

    // For backwards compatibility with profiles saved by RT 5.7
    const bool oldChannelScaling = params.filmNegative.redBase == -1.f;

    // If the film base values are not set in params, estimate multipliers from each channel's median value.
    if (params.filmNegative.redBase <= 0.f) {

        // Channel medians
        float rmed, gmed, bmed;

        if (oldChannelScaling) {
            // If using the old method, calculate nedians on the whole image
            calcMedians(baseImg, 0, 0, rwidth, rheight, rmed, gmed, bmed);
        } else {
            // The new method cuts out a 20% border from medians calculation.
            const int bW = rwidth * 20 / 100;
            const int bH = rheight * 20 / 100;
            calcMedians(baseImg, bW, bH, rwidth - bW, rheight - bH, rmed, gmed, bmed);
        }

        if (settings->verbose) {
            printf("Thumbnail input channel medians: %g %g %g\n", rmed, gmed, bmed);
        }

        // Calculate output medians
        rmed = powf(rmed, rexp);
        gmed = powf(gmed, gexp);
        bmed = powf(bmed, bexp);

        // Calculate output multipliers so that the median value is 1/24 of the output range.
        rmult = (MAX_OUT_VALUE / 24.f) / rmed;
        gmult = (MAX_OUT_VALUE / 24.f) / gmed;
        bmult = (MAX_OUT_VALUE / 24.f) / bmed;

    } else {

        // Read film-base values from params
        float rbase = params.filmNegative.redBase;
        float gbase = params.filmNegative.greenBase;
        float bbase = params.filmNegative.blueBase;

        // Reconstruct scale_mul coefficients from thumbnail metadata:
        //   redMultiplier / camwbRed is pre_mul[0]
        //   pre_mul[0] * scaleGain is scale_mul[0]
        // Apply channel scaling to raw (unscaled) input base values, to
        // match with actual (scaled) data in baseImg
        rbase *= (redMultiplier / camwbRed)     * scaleGain;
        gbase *= (greenMultiplier / camwbGreen) * scaleGain;
        bbase *= (blueMultiplier / camwbBlue)   * scaleGain;

        if (settings->verbose) {
            printf("Thumbnail input film base values: %g %g %g\n", rbase, gbase, bbase);
        }

        // Apply exponents to get output film base values
        rbase = powf(rbase, rexp);
        gbase = powf(gbase, gexp);
        bbase = powf(bbase, bexp);

        if (settings->verbose) {
            printf("Thumbnail output film base values: %g %g %g\n", rbase, gbase, bbase);
        }

        // Calculate multipliers so that film base value is 1/512th of the output range.
        rmult = (MAX_OUT_VALUE / 512.f) / rbase;
        gmult = (MAX_OUT_VALUE / 512.f) / gbase;
        bmult = (MAX_OUT_VALUE / 512.f) / bbase;

    }


    if (oldChannelScaling) {
        // Need to calculate channel averages, to fake the same conditions
        // found in rawimagesource, where get_ColorsCoeff is called with
        // forceAutoWB=true.
        float rsum = 0.f, gsum = 0.f, bsum = 0.f;

        for (int i = 0; i < rheight; i++) {
            for (int j = 0; j < rwidth; j++) {
                rsum += baseImg->r(i, j);
                gsum += baseImg->g(i, j);
                bsum += baseImg->b(i, j);
            }
        }

        const float ravg = rsum / (rheight * rwidth);
        const float gavg = gsum / (rheight * rwidth);
        const float bavg = bsum / (rheight * rwidth);

        // Shifting current WB multipliers, based on channel averages.
        rmult /= gavg / ravg;
        // gmult /= gAvg / gAvg;  green chosen as reference channel
        bmult /= gavg / bavg;

    } else {

        // Get and un-apply multipliers to adapt the thumbnail to a known fixed WB setting,
        // as in the main image processing.

        double r, g, b;
        ColorTemp(3500., 1., 1., "Custom").getMultipliers(r, g, b);
        //iColorMatrix is cam_rgb
        const double rm = camwbRed   / (iColorMatrix[0][0] * r + iColorMatrix[0][1] * g + iColorMatrix[0][2] * b);
        const double gm = camwbGreen / (iColorMatrix[1][0] * r + iColorMatrix[1][1] * g + iColorMatrix[1][2] * b);
        const double bm = camwbBlue  / (iColorMatrix[2][0] * r + iColorMatrix[2][1] * g + iColorMatrix[2][2] * b);

        // Normalize max WB multiplier to 1.f
        const double m = max(rm, gm, bm);
        rmult /= rm / m;
        gmult /= gm / m;
        bmult /= bm / m;
    }


    if (settings->verbose) {
        printf("Thumbnail computed multipliers: %g %g %g\n", static_cast<double>(rmult), static_cast<double>(gmult), static_cast<double>(bmult));
    }


#ifdef __SSE2__
    const vfloat clipv = F2V(MAXVALF);
    const vfloat rexpv = F2V(rexp);
    const vfloat gexpv = F2V(gexp);
    const vfloat bexpv = F2V(bexp);
    const vfloat rmultv = F2V(rmult);
    const vfloat gmultv = F2V(gmult);
    const vfloat bmultv = F2V(bmult);
#endif

    for (int i = 0; i < rheight; i++) {
        float *rline = baseImg->r(i);
        float *gline = baseImg->g(i);
        float *bline = baseImg->b(i);
        int j = 0;
#ifdef __SSE2__

        for (; j < rwidth - 3; j += 4) {
            STVFU(rline[j], vminf(rmultv * pow_F(LVFU(rline[j]), rexpv), clipv));
            STVFU(gline[j], vminf(gmultv * pow_F(LVFU(gline[j]), gexpv), clipv));
            STVFU(bline[j], vminf(bmultv * pow_F(LVFU(bline[j]), bexpv), clipv));
        }

#endif

        for (; j < rwidth; ++j) {
            rline[j] = CLIP(rmult * pow_F(rline[j], rexp));
            gline[j] = CLIP(gmult * pow_F(gline[j], gexp));
            bline[j] = CLIP(bmult * pow_F(bline[j], bexp));
        }
    }
}
