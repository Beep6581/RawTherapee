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
#include <cmath>
#include <iostream>

#include "rtengine.h"
#include "rawimagesource.h"
#include "mytime.h"
#include "procparams.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "opthelper.h"
#include "rt_algo.h"

namespace rtengine
{

extern const Settings* settings;

bool RawImageSource::channelsAvg(Coord spotPos, int spotSize, float avgs[3], const FilmNegativeParams &params)
{
    avgs[0] = avgs[1] = avgs[2] = 0.f;  // Channel averages

    if(ri->getSensorType() != ST_BAYER && ri->getSensorType() != ST_FUJI_XTRANS)
        return false;

    if (settings->verbose)
        printf("Spot coord:  x=%d y=%d\n", spotPos.x, spotPos.y);

    int x1 = spotPos.x - spotSize / 2;
    int x2 = spotPos.x + spotSize / 2;
    int y1 = spotPos.y - spotSize / 2;
    int y2 = spotPos.y + spotSize / 2;

    if(x1<0 || x2>W || y1<0 || y2>H)
        return false; // Spot goes outside bounds, bail out.

    int pxCount[3] = {0}; // Per-channel sample counts
    for(int c=spotPos.x-spotSize; c<spotPos.x+spotSize; c++) {
        for(int r=spotPos.y-spotSize; r<spotPos.y+spotSize; r++) {

            int ch = (ri->getSensorType() == ST_BAYER) ? FC(r,c) : ri->XTRANSFC(r,c);

            pxCount[ch]++;
            // If film negative is currently enabled, undo the effect by elevating to 1/exp,
            // in order to sample the original, linear value
            if(params.enabled)
                avgs[ch] += powf(rawData[r][c], -1 / (ch==0 ? params.redExp : ch==1 ? params.greenExp : params.blueExp));
            else
                avgs[ch] += rawData[r][c];
        }
    }

    for(int ch=0; ch<3; ch++)
        avgs[ch] = avgs[ch] / (pxCount[ch]);

    return true;
}

// Calculate logarithms in arbitrary base
float logBase(float base, float num) {
    return log(num) / log(base);
}

bool RawImageSource::getFilmNegativeExponents (Coord2D spotA, Coord2D spotB, int tran, const FilmNegativeParams &currentParams, float newExps[3])
{
    float clearVals[3], denseVals[3];

    newExps[0] = currentParams.redExp;
    newExps[1] = currentParams.greenExp;
    newExps[2] = currentParams.blueExp;

    int spotSize = 32; // TODO : make this confugurable ?
    Coord spot;
    // Sample first spot
    transformPosition (spotA.x, spotA.y, tran, spot.x, spot.y);
    if(!channelsAvg(spot, spotSize, clearVals, currentParams))
      return false;

    // Sample second spot
    transformPosition (spotB.x, spotB.y, tran, spot.x, spot.y);
    if(!channelsAvg(spot, spotSize, denseVals, currentParams))
      return false;

    // Detect which one is the dense spot, based on green channel
    if(clearVals[1] < denseVals[1])
      std::swap(clearVals, denseVals);

    if (settings->verbose) {
        printf("Clear film values:  R=%f G=%f B=%f\n", clearVals[0], clearVals[1], clearVals[2]);
        printf("Dense film values:  R=%f G=%f B=%f\n", denseVals[0], denseVals[1], denseVals[2]);
    }

    float denseGreenRatio = clearVals[1] / denseVals[1];

    // Calculate exponents for each channel, based on the ratio between the bright and dark values,
    // compared to the ratio in the reference channel (green)
    for(int ch=0; ch<3; ch++)
        if(ch==1)
            newExps[ch] = 1.f;  // Green is the reference channel
        else
            newExps[ch] = CLAMP(logBase(clearVals[ch] / denseVals[ch], denseGreenRatio), 0.3f, 3.f);

    if (settings->verbose)
        printf("New exponents:  R=%f G=%f B=%f\n", newExps[0], newExps[1], newExps[2]);

    return true;
}

void RawImageSource::filmNegativeProcess(const procparams::FilmNegativeParams &params)
{
    if(!params.enabled)
        return;

    float exps[3] = { (float)params.redExp, (float)params.greenExp, (float)params.blueExp };
    
    MyTime t1, t2, t3,t4, t5, t6;
    t1.set();

    if(ri->getSensorType() == ST_BAYER) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int row = 0; row < H; row ++) {
            int col = 0;
            // Exponents are expressed as positive in the parameters, so negate them in order
            // to get the reciprocals. Avoid trouble with zeroes, minimum pixel value is 1.
            const float exps0 = -exps[FC(row, col)];
            const float exps1 = -exps[FC(row, col + 1)];
#ifdef __SSE2__
            const vfloat expsv = _mm_setr_ps(exps0, exps1, exps0, exps1);
            const vfloat onev = F2V(1.f);
            for (; col < W - 3; col+=4) {
                STVFU(rawData[row][col], pow_F(vmaxf(LVFU(rawData[row][col]), onev), expsv));
            }
#endif // __SSE2__
            for (; col < W - 1; col+=2) {
                rawData[row][col] = pow_F(max(rawData[row][col], 1.f), exps0);
                rawData[row][col + 1] = pow_F(max(rawData[row][col + 1], 1.f), exps1);
            }
            if (col < W) {
                rawData[row][col] = pow_F(max(rawData[row][col], 1.f), exps0);
            }
        }
    } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int row = 0; row < H; row ++) {
            int col = 0;
            // Exponents are expressed as positive in the parameters, so negate them in order
            // to get the reciprocals. Avoid trouble with zeroes, minimum pixel value is 1.
            const float expsc[6] = {-exps[ri->XTRANSFC(row, 0)], -exps[ri->XTRANSFC(row, 1)], -exps[ri->XTRANSFC(row, 2)], -exps[ri->XTRANSFC(row, 3)], -exps[ri->XTRANSFC(row, 4)], -exps[ri->XTRANSFC(row, 5)]};
#ifdef __SSE2__
            const vfloat expsv0 = _mm_setr_ps(expsc[0], expsc[1], expsc[2], expsc[3]);
            const vfloat expsv1 = _mm_setr_ps(expsc[4], expsc[5], expsc[0], expsc[1]);
            const vfloat expsv2 = _mm_setr_ps(expsc[2], expsc[3], expsc[4], expsc[5]);
            const vfloat onev = F2V(1.f);
            for (; col < W - 11; col+=12) {
                STVFU(rawData[row][col], pow_F(vmaxf(LVFU(rawData[row][col]), onev), expsv0));
                STVFU(rawData[row][col + 4], pow_F(vmaxf(LVFU(rawData[row][col + 4]), onev), expsv1));
                STVFU(rawData[row][col + 8], pow_F(vmaxf(LVFU(rawData[row][col + 8]), onev), expsv2));
            }
#endif // __SSE2__
            for (; col < W - 5; col+=6) {
                for (int c = 0; c < 6; ++c) {
                    rawData[row][col + c] = pow_F(max(rawData[row][col + c], 1.f), expsc[c]);
                }
            }
            for (int c = 0; col < W; col++, c++) {
                rawData[row][col + c] = pow_F(max(rawData[row][col + c], 1.f), expsc[c]);
            }
        }
    }


    t2.set();
    if (settings->verbose)
        printf("Pow loop time us: %d\n", t2.etime(t1));

    // Channel vectors to calculate medians
    std::vector<float> cvs[3] = {
        std::vector<float>(),
        std::vector<float>(),
        std::vector<float>()
    };

    // Sample one every 5 pixels, and push the value in the appropriate channel vector.
    // Chose an odd step, not multiple of the CFA size, to get a chance to visit each channel.
    if(ri->getSensorType() == ST_BAYER) {
        for (int row = 0; row < H; row+=5) {
            for (int col = 0; col < W; col+=5) {
                int c  = FC(row, col);                        // three colors,  0=R, 1=G,  2=B
                cvs[c].push_back(rawData[row][col]);
            }
        }
    } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
        for (int row = 0; row < H; row+=5) {
            for (int col = 0; col < W; col+=5) {
                int c  = ri->XTRANSFC(row, col);                        // three colors,  0=R, 1=G,  2=B
                cvs[c].push_back(rawData[row][col]);
            }
        }
    }

    const float MAX_OUT_VALUE = 65000.f;

    t3.set();
    if (settings->verbose)
        printf("Median vector fill loop time us: %d\n", t3.etime(t2));

    float medians[3];  // Channel median values
    float mults[3] = { 1.f };  // Channel normalization multipliers

    for (int c=0; c<3; c++) {
        // Find median values for each channel using a histogram search function
        findMinMaxPercentile(&cvs[c][0], cvs[c].size(), 0.5f, medians[c], 0.5f, medians[c], true);
        // Determine the channel multipler so that N times the median becomes 65k. This clips away
        // the values in the dark border surrounding the negative (due to the film holder, for example),
        // the reciprocal of which have blown up to stellar values.
        mults[c] = MAX_OUT_VALUE / (medians[c] * 24);
    }


    t4.set();
    if (settings->verbose) {
        printf("Sample count : %lu, %lu, %lu\n", cvs[0].size(), cvs[1].size(), cvs[2].size());
        printf("Medians : %f %f %f\n", medians[0], medians[1], medians[2] );
        printf("Computed multipliers : %f %f %f\n", mults[0], mults[1], mults[2] );
        printf("Median calc time us: %d\n", t4.etime(t3));
    }


    if(ri->getSensorType() == ST_BAYER) {

#ifdef _OPENMP
        #pragma omp for nowait
#endif
        for (int row = 0; row < H; row ++) {
            for (int col = 0; col < W; col++) {
                int c  = FC(row, col);                        // three colors,  0=R, 1=G,  2=B
                // Apply the multipliers
                rawData[row][col] *= mults[c];
            }
        }
        
    } else if(ri->getSensorType() == ST_FUJI_XTRANS) {

#ifdef _OPENMP
        #pragma omp for nowait
#endif
        for (int row = 0; row < H; row ++) {
            for (int col = 0; col < W; col++) {
                int c  = ri->XTRANSFC(row, col);                        // three colors,  0=R, 1=G,  2=B
                // Apply the multipliers
                rawData[row][col] *= mults[c];
            }
        }

    }


    t5.set();
    if (settings->verbose)
        printf("Mult loop time us: %d\n", t5.etime(t4));


    PixelsMap bitmapBads(W, H);

    int totBP = 0; // Hold count of bad pixels to correct

    if(ri->getSensorType() == ST_BAYER) {


#ifdef _OPENMP
        #pragma omp parallel for reduction(+:totBP) schedule(dynamic,16)
#endif

        for(int i = 0; i < H; i++)
            for(int j = 0; j < W; j++) {
                if (rawData[i][j] >= MAX_OUT_VALUE) {
                    bitmapBads.set(j, i);
                    totBP++;
                }
            }

        if (totBP > 0) {
            interpolateBadPixelsBayer( bitmapBads, rawData );
        }

    } else if(ri->getSensorType() == ST_FUJI_XTRANS) {

        // TODO

    }

    t6.set();
    if (settings->verbose) {
        printf("Bad pixels count: %d\n", totBP);
        printf("Bad pixels interpolation time us: %d\n", t6.etime(t5));
    }
}

}