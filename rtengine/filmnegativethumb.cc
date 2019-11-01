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

#include "rtengine.h"
#include "rtthumbnail.h"
#include "opthelper.h"
#include "rt_algo.h"
#include "rtengine.h"
#include "settings.h"
#include "procparams.h"
#define BENCHMARK
#include "StopWatch.h"

namespace rtengine
{

extern const Settings* settings;

}

void rtengine::Thumbnail::processFilmNegative(
    const procparams::ProcParams &params,
    const Imagefloat* baseImg,
    const int rwidth, const int rheight,
    float &rmi, float &gmi, float &bmi
)
{

    // Channel exponents
    const float rexp = -params.filmNegative.redRatio * params.filmNegative.greenExp;
    const float gexp = -params.filmNegative.greenExp;
    const float bexp = -params.filmNegative.blueRatio * params.filmNegative.greenExp;

    // If channel medians are not set in params, use previous method.
    // For backwards compatibility with profiles saved by RT 5.7
    const bool oldChannelScaling = params.filmNegative.redMedian == 0.f;

    // Reference channel medians
    float rmed, gmed, bmed;

    if (oldChannelScaling) {

        // Need to calculate channel averages, to fake the same conditions
        // found in rawimagesource, where get_ColorsCoeff is called with
        // forceAutoWB=true.
        float rsum = 0.f, gsum = 0.f, bsum = 0.f;

        // Channel vectors to calculate medians
        std::vector<float> rv, gv, bv;

        for (int i = 0; i < rheight; i++) {
            for (int j = 0; j < rwidth; j++) {
                const float r = baseImg->r(i, j);
                const float g = baseImg->g(i, j);
                const float b = baseImg->b(i, j);

                rsum += r;
                gsum += g;
                bsum += b;

                rv.push_back(r);
                gv.push_back(g);
                bv.push_back(b);
            }
        }

        const float ravg = rsum / (rheight * rwidth);
        const float gavg = gsum / (rheight * rwidth);
        const float bavg = bsum / (rheight * rwidth);

        // Shifting current WB multipliers, based on channel averages.
        rmi /= (gavg / ravg);
        // gmi /= (gAvg/gAvg);  green chosen as reference channel
        bmi /= (gavg / bavg);


        // Calculate channel medians from whole image
        findMinMaxPercentile(rv.data(), rv.size(), 0.5f, rmed, 0.5f, rmed, true);
        findMinMaxPercentile(gv.data(), gv.size(), 0.5f, gmed, 0.5f, gmed, true);
        findMinMaxPercentile(bv.data(), bv.size(), 0.5f, bmed, 0.5f, bmed, true);

        if (settings->verbose) {
            printf("Thumbnail input channel medians: %g %g %g\n", rmed, gmed, bmed);
        }

        // Calculate output medians
        rmed = powf(rmed, rexp);
        gmed = powf(gmed, gexp);
        bmed = powf(bmed, bexp);

    } else {

        // Read reference channel medians from params
        rmed = params.filmNegative.redMedian;
        gmed = params.filmNegative.greenMedian;
        bmed = params.filmNegative.blueMedian;

        // Reconstruct scale_mul coefficients from thumbnail metadata:
        //   redMultiplier / camwbRed is pre_mul[0]
        //   pre_mul[0] * scaleGain is scale_mul[0]
        // Apply channel scaling to raw (unscaled) reference input medians, to
        // match with actual (scaled) data in baseImg
        rmed *= (redMultiplier / camwbRed)     * scaleGain;
        gmed *= (greenMultiplier / camwbGreen) * scaleGain;
        bmed *= (blueMultiplier / camwbBlue)   * scaleGain;

        if (settings->verbose) {
            printf("Thumbnail input channel medians: %g %g %g\n", rmed, gmed, bmed);
        }

        // Calculate output medians
        rmed = powf(rmed, rexp);
        gmed = powf(gmed, gexp);
        bmed = powf(bmed, bexp);

    }

    const float MAX_OUT_VALUE = 65000.f;

    // Calculate output multipliers
    const float rmult = MAX_OUT_VALUE / (rmed * 24);
    const float gmult = MAX_OUT_VALUE / (gmed * 24);
    const float bmult = MAX_OUT_VALUE / (bmed * 24);

    if (settings->verbose) {
        printf("Thumbnail output channel medians: %g %g %g\n", rmed, gmed, bmed);
        printf("Thumbnail computed multipliers: %g %g %g\n", rmult, gmult, bmult);
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
