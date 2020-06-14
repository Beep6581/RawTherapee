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
#include <iostream>

#include "rawimage.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "improccoordinator.h"
#include "imagefloat.h"

#include "coord.h"
#include "mytime.h"
#include "opthelper.h"
#include "pixelsmap.h"
#include "procparams.h"
#include "rt_algo.h"
#include "rtengine.h"
#include "rtthumbnail.h"
#include "sleef.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "iccstore.h"
#include "rt_math.h"
#include "color.h"


namespace
{
using rtengine::settings;
using rtengine::Coord2D;
Coord2D translateCoord(const rtengine::ImProcFunctions& ipf, int fw, int fh, int x, int y) {

    const std::vector<Coord2D> points = {Coord2D(x, y)};

    std::vector<Coord2D> red;
    std::vector<Coord2D> green;
    std::vector<Coord2D> blue;
    ipf.transCoord(fw, fh, points, red, green, blue);

    return green[0];
}


constexpr double FILMNEG_MATRIX[3][3] = {
    // FilmNeg
    {0.47734,	0.34436,	0.1425 },
    {0.18016,	0.79044,	0.0294 },
    {0	    ,   0.00696,	0.81795}
};

constexpr double FILMNEG_INV_MATRIX[3][3] = {
    // FilmNeg Inv
    {  2.505614 , -1.088087 , -0.397408 },
    { -0.571270 ,  1.513598 ,  0.045120 },
    {  0.004861 , -0.012879 ,  1.222185 }
};


void getSpotAvgMax(rtengine::ImageSource *imgsrc, rtengine::ColorTemp currWB, const std::unique_ptr<rtengine::procparams::ProcParams> &params,
                 Coord2D p, int tr, int spotSize, std::array<float, 3> &rawValues, std::array<float, 3> &maxValues)
{
    using rtengine::Color;

    int x1 = MAX(0, (int)p.x - spotSize/2);
    int y1 = MAX(0, (int)p.y - spotSize/2);
    PreviewProps pp(x1, y1, spotSize, spotSize, 1);

    if (settings->verbose) {
        printf("Spot: %d,%d   %d,%d\n", x1, y1, x1+spotSize/2, y1+spotSize/2);
    }

    rtengine::Imagefloat spotImg(spotSize,spotSize);
    imgsrc->getImage(currWB, tr, &spotImg, pp, params->toneCurve, params->raw);
    imgsrc->convertColorSpace(&spotImg, params->icm, currWB);

    rtengine::TMatrix wprof = rtengine::ICCStore::getInstance()->workingSpaceMatrix (params->icm.workingProfile);

    float ravg = 0.f, gavg = 0.f, bavg = 0.f;
    float rmax = 0.f, gmax = 0.f, bmax = 0.f;
    for (int i = 0; i < spotSize; ++i) {
        for (int j = 0; j < spotSize; ++j) {

            float x, y, z;
            Color::rgbxyz (spotImg.r(i,j), spotImg.g(i,j), spotImg.b(i,j), x, y, z, wprof);
            float r, g, b;
            Color::xyz2rgb (x, y, z, r, g, b, FILMNEG_INV_MATRIX );

            ravg += r;
            gavg += g;
            bavg += b;
            rmax = MAX(rmax, r);
            gmax = MAX(gmax, g);
            bmax = MAX(bmax, b);
        }
    }

    rawValues[0] = ravg / (spotSize*spotSize);
    rawValues[1] = gavg / (spotSize*spotSize);
    rawValues[2] = bavg / (spotSize*spotSize);

    if (settings->verbose) {
        printf("Spot RGB: %f,%f,%f\n", rawValues[0], rawValues[1], rawValues[2]);
    }


    maxValues[0] = rmax;
    maxValues[1] = gmax;
    maxValues[2] = bmax;
}

void getSpotAvg(rtengine::ImageSource *imgsrc, rtengine::ColorTemp currWB, const std::unique_ptr<rtengine::procparams::ProcParams> &params,
                 Coord2D p, int tr, int spotSize, std::array<float, 3> &rawValues)
{
    std::array<float, 3> dummy = {};
    getSpotAvgMax(imgsrc, currWB, params, p, tr, spotSize, rawValues, dummy);
}


void calcMedians(
    const rtengine::Imagefloat* input,
    int x1, int y1, int x2, int y2,
    float &rmed, float &gmed, float &bmed
)
{
    using rtengine::findMinMaxPercentile;

    // Channel vectors to calculate medians
    std::vector<float> rv, gv, bv;

    const int sz = (x2 - x1) * (y2 - y1);
    rv.reserve(sz);
    gv.reserve(sz);
    bv.reserve(sz);


    for (int ii = y1; ii < y2; ii ++) {
        for (int jj = x1; jj < x2; jj ++) {
            rv.push_back( input->r (ii, jj) );
            gv.push_back( input->g (ii, jj) );
            bv.push_back( input->b (ii, jj) );
        }
    }

    // Calculate channel medians from whole image
    findMinMaxPercentile(rv.data(), rv.size(), 0.5f, rmed, 0.5f, rmed, true);
    findMinMaxPercentile(gv.data(), gv.size(), 0.5f, gmed, 0.5f, gmed, true);
    findMinMaxPercentile(bv.data(), bv.size(), 0.5f, bmed, 0.5f, bmed, true);
}


/**
 * Calculate red and blue input channel balance, so that the input r,g,b values will
 * result in neutral gray after exponentiation
 */
void calcBalance(float baseVal, float rexp, float gexp, float bexp,
                 float r, float g, float b,
                 float &rBal, float &bBal)
{
    float rout = pow_F(r, rexp) / pow_F(baseVal, rexp);
    float gout = pow_F(g, gexp) / pow_F(baseVal, gexp);
    float bout = pow_F(b, bexp) / pow_F(baseVal, bexp);

    rBal = pow_F(rout / gout, (1 / rexp));
    bBal = pow_F(bout / gout, (1 / bexp));
}

}



void rtengine::ImProcFunctions::filmNegativeProcess(rtengine::Imagefloat *input, rtengine::Imagefloat *output,
        const procparams::FilmNegativeParams &params, std::array<float, 3>& filmBaseValues)
{

    if (!params.enabled) {
        return;
    }

    float rexp = -(params.greenExp * params.redRatio);
    float gexp = -params.greenExp;
    float bexp = -(params.greenExp * params.blueRatio);


    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (this->params->icm.workingProfile);
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (this->params->icm.workingProfile);

    
#ifdef _OPENMP
        #pragma omp parallel for
#endif

    for (int i = 0; i < input->getHeight(); i++) {
        for (int j = 0; j < input->getWidth(); j++) {

            float newr = input->r(i, j);
            float newg = input->g(i, j);
            float newb = input->b(i, j);
            float x, y, z;
            Color::rgbxyz (newr, newg, newb, x, y, z, wprof);
            Color::xyz2rgb ( x, y, z, newr, newg, newb, FILMNEG_INV_MATRIX );

            input->r(i, j) = newr;
            input->g(i, j) = newg;
            input->b(i, j) = newb;

        }
    }

    float rmult, gmult, bmult;

    // If film base values are set in params, use those
    if (filmBaseValues[1] <= 0.f) {
        // ...otherwise, the film negative tool might have just been enabled on this image,
        // whithout any previous setting. So, estimate film base values from channel medians

        float rmed, gmed, bmed;
        // Cut 20% border from medians calculation. It will probably contain outlier values
        // from the film holder, which will bias the median result.
        const int bW = input->getWidth() * 20 / 100;
        const int bH = input->getHeight() * 20 / 100;
        calcMedians(input, bW, bH,
            input->getWidth() - bW, input->getHeight() - bH,
            rmed, gmed, bmed);

        // // Apply exponents to get output film base values
        // rmed = powf(rmed, rexp);
        // gmed = powf(gmed, gexp);
        // bmed = powf(bmed, bexp);

        if (settings->verbose) {
            printf("Channel medians: R=%g, G=%g, B=%g\n", rmed, gmed, bmed);
        }

        // Estimate film base values, so that in the output data, each channel
        // median will correspond to 1/24th of MAX.
        filmBaseValues[0] = pow_F(24.f / 512.f, 1.f / rexp) * rmed;
        filmBaseValues[1] = pow_F(24.f / 512.f, 1.f / gexp) * gmed;
        filmBaseValues[2] = pow_F(24.f / 512.f, 1.f / bexp) * bmed;

    }

    // Apply channel exponents to base values, and scale the output to 1/512 of max.
    rmult = (MAXVALF / 512.f) / pow_F(rtengine::max(filmBaseValues[0], 1.f), rexp);
    gmult = (MAXVALF / 512.f) / pow_F(rtengine::max(filmBaseValues[1], 1.f), gexp);
    bmult = (MAXVALF / 512.f) / pow_F(rtengine::max(filmBaseValues[2], 1.f), bexp);


#ifdef __SSE2__
    const vfloat clipv = F2V(MAXVALF);
    const vfloat rexpv = F2V(rexp);
    const vfloat gexpv = F2V(gexp);
    const vfloat bexpv = F2V(bexp);
    const vfloat rmultv = F2V(rmult);
    const vfloat gmultv = F2V(gmult);
    const vfloat bmultv = F2V(bmult);
#endif

    const int rheight = input->getHeight();
    const int rwidth = input->getWidth();

    for (int i = 0; i < rheight; i++) {
        float *rlinein = input->r(i);
        float *glinein = input->g(i);
        float *blinein = input->b(i);
        float *rlineout = output->r(i);
        float *glineout = output->g(i);
        float *blineout = output->b(i);
        int j = 0;
#ifdef __SSE2__

        for (; j < rwidth - 3; j += 4) {
            STVFU(rlineout[j], vminf(rmultv * pow_F(LVFU(rlinein[j]), rexpv), clipv));
            STVFU(glineout[j], vminf(gmultv * pow_F(LVFU(glinein[j]), gexpv), clipv));
            STVFU(blineout[j], vminf(bmultv * pow_F(LVFU(blinein[j]), bexpv), clipv));
        }

#endif

        for (; j < rwidth; ++j) {
            rlineout[j] = CLIP(rmult * pow_F(rlinein[j], rexp));
            glineout[j] = CLIP(gmult * pow_F(glinein[j], gexp));
            blineout[j] = CLIP(bmult * pow_F(blinein[j], bexp));
        }
    }


    // Convert back to working profile
#ifdef _OPENMP
        #pragma omp parallel for
#endif
    for (int i = 0; i < output->getHeight(); i++) {
        for (int j = 0; j < output->getWidth(); j++) {

            float newr = output->r(i, j);
            float newg = output->g(i, j);
            float newb = output->b(i, j);
            float x, y, z;

            Color::rgbxyz (newr, newg, newb, x, y, z, FILMNEG_MATRIX);
            Color::xyz2rgb ( x, y, z, newr, newg, newb, wiprof );

            output->r(i, j) = newr;
            output->g(i, j) = newg;
            output->b(i, j) = newb;

        }
    }

}


bool rtengine::ImProcCoordinator::getFilmNegativeExponents(int xA, int yA, int xB, int yB, std::array<float, 3>& newExps, float &rBal, float &bBal)
{
    MyMutex::MyLock lock(mProcessing);

    const int tr = getCoarseBitMask(params->coarse);

    const Coord2D p1 = translateCoord(ipf, fw, fh, xA, yA);
    const Coord2D p2 = translateCoord(ipf, fw, fh, xB, yB);

    // TEMP : Legacy mode, for backwards compatibility with RT 5.8 processing profiles
    if (imgsrc->isRAW() && params->filmNegative.greenBase == -1) {
        return imgsrc->getFilmNegativeExponents(p1, p2, tr, params->filmNegative, newExps);
    }

    FilmNegativeParams currentParams = params->filmNegative;

    newExps = {
        static_cast<float>(currentParams.redRatio * currentParams.greenExp),
        static_cast<float>(currentParams.greenExp),
        static_cast<float>(currentParams.blueRatio * currentParams.greenExp)
    };

    constexpr int spotSize = 32; // TODO: Make this configurable?

    std::array<float, 3> clearVals;
    getSpotAvg(imgsrc, currWB, params, p1, tr, spotSize, clearVals);

    std::array<float, 3> denseVals;
    getSpotAvg(imgsrc, currWB, params, p2, tr, spotSize, denseVals);

    // Detect which one is the dense spot, based on green channel
    if (clearVals[1] < denseVals[1]) {
        std::swap(clearVals, denseVals);
    }

    if (settings->verbose) {
        printf("Clear film values: R=%g G=%g B=%g\n", static_cast<double>(clearVals[0]), static_cast<double>(clearVals[1]), static_cast<double>(clearVals[2]));
        printf("Dense film values: R=%g G=%g B=%g\n", static_cast<double>(denseVals[0]), static_cast<double>(denseVals[1]), static_cast<double>(denseVals[2]));
    }

    const float denseGreenRatio = clearVals[1] / denseVals[1];

    // Calculate logarithms in arbitrary base
    const auto logBase =
        [](float base, float num) -> float
        {
            return std::log(num) / std::log(base);
        };

    // Calculate exponents for each channel, based on the ratio between the bright and dark values,
    // compared to the ratio in the reference channel (green)
    for (int ch = 0; ch < 3; ++ch) {
        if (ch == 1) {
            newExps[ch] = 1.f;  // Green is the reference channel
        } else {
            newExps[ch] = rtengine::LIM(logBase(clearVals[ch] / denseVals[ch], denseGreenRatio), 0.3f, 4.f);
        }
    }

    if (settings->verbose) {
        printf("New exponents:  R=%g G=%g B=%g\n", static_cast<double>(newExps[0]), static_cast<double>(newExps[1]), static_cast<double>(newExps[2]));
    }

    // Re-adjust color balance based on dense spot values and new exponents
    calcBalance(rtengine::max(static_cast<float>(params->filmNegative.greenBase), 1.f),
        -newExps[0], -newExps[1], -newExps[2],
        denseVals[0], denseVals[1], denseVals[2],
        rBal, bBal);

    return true;

}


float rtengine::ImProcCoordinator::getFilmBaseGreen(int x, int y, const int spotSize)
{
    MyMutex::MyLock lock(mProcessing);

    const int tr = getCoarseBitMask(params->coarse);

    const Coord2D p = translateCoord(ipf, fw, fh, x, y);
    
    // Return the maximum green value from the sampled spot
    std::array<float,3> avg = {}, max = {};
    getSpotAvgMax(imgsrc, currWB, params, p, tr, spotSize, avg, max);

    return max[1];
}


bool rtengine::ImProcCoordinator::getFilmNegativeBalance(int x, int y, int spotSize, float &rBal, float &bBal)
{
    MyMutex::MyLock lock(mProcessing);

    const int tr = getCoarseBitMask(params->coarse);

    const Coord2D p = translateCoord(ipf, fw, fh, x, y);

    std::array<float, 3> vals = { };
    getSpotAvg(imgsrc, currWB, params, p, tr, spotSize, vals);

    if (vals[0] <= 0.f || vals[1] <= 0.f || vals[2] <= 0.f) {
        return false;
    }

    calcBalance(rtengine::max(this->filmBaseValues[1], 1.f),
        - params->filmNegative.greenExp * params->filmNegative.redRatio,
        - params->filmNegative.greenExp,
        - params->filmNegative.greenExp * params->filmNegative.blueRatio,
        vals[0], vals[1], vals[2],
        rBal, bBal);

    return true;
}







// ---------- >>> legacy mode >>> ---------------

namespace
{
using rtengine::ST_BAYER;
using rtengine::ST_FUJI_XTRANS;
using rtengine::settings;

bool channelsAvg(
    const rtengine::RawImage* ri,
    int width,
    int height,
    const float* cblacksom,
    rtengine::Coord spotPos,
    int spotSize,
    std::array<float, 3>& avgs
)
{
    avgs = {}; // Channel averages

    if (ri->getSensorType() != ST_BAYER && ri->getSensorType() != ST_FUJI_XTRANS) {
        return false;
    }

    if (settings->verbose) {
        printf("Spot coord:  x=%d y=%d\n", spotPos.x, spotPos.y);
    }

    const int half_spot_size = spotSize / 2;

    const int& x1 = spotPos.x - half_spot_size;
    const int& x2 = spotPos.x + half_spot_size;
    const int& y1 = spotPos.y - half_spot_size;
    const int& y2 = spotPos.y + half_spot_size;

    if (x1 < 0 || x2 > width || y1 < 0 || y2 > height) {
        return false; // Spot goes outside bounds, bail out.
    }

    std::array<int, 3> pxCount = {}; // Per-channel sample counts

    for (int c = x1; c < x2; ++c) {
        for (int r = y1; r < y2; ++r) {
            const int ch = ri->getSensorType() == ST_BAYER ? ri->FC(r, c) : ri->XTRANSFC(r, c);

            ++pxCount[ch];

            // Sample the original unprocessed values from RawImage, subtracting black levels.
            // Scaling is irrelevant, as we are only interested in the ratio between two spots.
            avgs[ch] += ri->data[r][c] - cblacksom[ch];
        }
    }

    for (int ch = 0; ch < 3; ++ch) {
        avgs[ch] /= pxCount[ch];
    }

    return true;
}


}


bool rtengine::RawImageSource::getFilmNegativeExponents(Coord2D spotA, Coord2D spotB, int tran, const procparams::FilmNegativeParams &currentParams, std::array<float, 3>& newExps)
{
    newExps = {
        static_cast<float>(currentParams.redRatio * currentParams.greenExp),
        static_cast<float>(currentParams.greenExp),
        static_cast<float>(currentParams.blueRatio * currentParams.greenExp)
    };

    constexpr int spotSize = 32; // TODO: Make this configurable?

    Coord spot;
    std::array<float, 3> clearVals;
    std::array<float, 3> denseVals;

    // Get channel averages in the two spots, sampling from the original ri->data buffer.
    // NOTE: rawData values might be affected by CA corection, FlatField, etc, so:
    //   rawData[y][x] == (ri->data[y][x] - cblacksom[c]) * scale_mul[c]
    // is not always true. To calculate exponents on the exact values, we should keep
    // a copy of the rawData buffer after preprocessing. Worth the memory waste?

    // Sample first spot
    transformPosition(spotA.x, spotA.y, tran, spot.x, spot.y);

    if (!channelsAvg(ri, W, H, cblacksom, spot, spotSize, clearVals)) {
        return false;
    }

    // Sample second spot
    transformPosition(spotB.x, spotB.y, tran, spot.x, spot.y);

    if (!channelsAvg(ri, W, H, cblacksom, spot, spotSize, denseVals)) {
        return false;
    }

    // Detect which one is the dense spot, based on green channel
    if (clearVals[1] < denseVals[1]) {
        std::swap(clearVals, denseVals);
    }

    if (settings->verbose) {
        printf("Clear film values: R=%g G=%g B=%g\n", static_cast<double>(clearVals[0]), static_cast<double>(clearVals[1]), static_cast<double>(clearVals[2]));
        printf("Dense film values: R=%g G=%g B=%g\n", static_cast<double>(denseVals[0]), static_cast<double>(denseVals[1]), static_cast<double>(denseVals[2]));
    }

    const float denseGreenRatio = clearVals[1] / denseVals[1];

    // Calculate logarithms in arbitrary base
    const auto logBase =
        [](float base, float num) -> float
        {
            return std::log(num) / std::log(base);
        };

    // Calculate exponents for each channel, based on the ratio between the bright and dark values,
    // compared to the ratio in the reference channel (green)
    for (int ch = 0; ch < 3; ++ch) {
        if (ch == 1) {
            newExps[ch] = 1.f;  // Green is the reference channel
        } else {
            newExps[ch] = rtengine::LIM(logBase(clearVals[ch] / denseVals[ch], denseGreenRatio), 0.3f, 4.f);
        }
    }

    if (settings->verbose) {
        printf("New exponents:  R=%g G=%g B=%g\n", static_cast<double>(newExps[0]), static_cast<double>(newExps[1]), static_cast<double>(newExps[2]));
    }

    return true;
}



void rtengine::RawImageSource::filmNegativeProcess(const procparams::FilmNegativeParams &params)
{
//    BENCHFUNMICRO

    if (!params.enabled) {
        return;
    }

    // Exponents are expressed as positive in the parameters, so negate them in order
    // to get the reciprocals.
    const std::array<float, 3> exps = {
        static_cast<float>(-params.redRatio * params.greenExp),
        static_cast<float>(-params.greenExp),
        static_cast<float>(-params.blueRatio * params.greenExp)
    };

    MyTime t1, t2, t3,t4, t5;

    t1.set();

    // Channel vectors to calculate medians
    std::array<std::vector<float>, 3> cvs;

    // Sample one every 5 pixels, and push the value in the appropriate channel vector.
    // Choose an odd step, not a multiple of the CFA size, to get a chance to visit each channel.
    if (ri->getSensorType() == ST_BAYER) {
        for (int row = 0; row < H; row += 5) {
            const int c0 = ri->FC(row, 0);
            const int c1 = ri->FC(row, 5);
            int col = 0;
            for (; col < W - 5; col += 10) {
                cvs[c0].push_back(rawData[row][col]);
                cvs[c1].push_back(rawData[row][col + 5]);
            }
            if (col < W) {
                cvs[c0].push_back(rawData[row][col]);
            }
        }
    }
    else if (ri->getSensorType() == ST_FUJI_XTRANS) {
        for (int row = 0; row < H; row += 5) {
            const std::array<unsigned int, 6> cs = {
                ri->XTRANSFC(row, 0),
                ri->XTRANSFC(row, 5),
                ri->XTRANSFC(row, 10),
                ri->XTRANSFC(row, 15),
                ri->XTRANSFC(row, 20),
                ri->XTRANSFC(row, 25)
            };
            int col = 0;
            for (; col < W - 25; col += 30) {
                for (int c = 0; c < 6; ++c) {
                    cvs[cs[c]].push_back(rawData[row][col + c * 5]);
                }
            }
            for (int c = 0; col < W; col += 5, ++c) {
                cvs[cs[c]].push_back(rawData[row][col]);
            }
        }
    }

    constexpr float MAX_OUT_VALUE = 65000.f;

    t2.set();

    if (settings->verbose) {
        printf("Median vector fill loop time us: %d\n", t2.etime(t1));
    }

    t2.set();

    std::array<float, 3> medians; // Channel median values
    std::array<float, 3> mults = {
        1.f,
        1.f,
        1.f
    };  // Channel normalization multipliers

    for (int c = 0; c < 3; ++c) {
        // Find median values for each channel
        if (!cvs[c].empty()) {
            findMinMaxPercentile(cvs[c].data(), cvs[c].size(), 0.5f, medians[c], 0.5f, medians[c], true);
            medians[c] = pow_F(rtengine::max(medians[c], 1.f), exps[c]);
            // Determine the channel multiplier so that N times the median becomes 65k. This clips away
            // the values in the dark border surrounding the negative (due to the film holder, for example),
            // the reciprocal of which have blown up to stellar values.
            mults[c] = MAX_OUT_VALUE / (medians[c] * 24.f);
        }
    }

    t3.set();

    if (settings->verbose) {
        printf("Sample count: %zu, %zu, %zu\n", cvs[0].size(), cvs[1].size(), cvs[2].size());
        printf("Medians: %g %g %g\n", medians[0], medians[1], medians[2] );
        printf("Computed multipliers: %g %g %g\n", mults[0], mults[1], mults[2] );
        printf("Median calc time us: %d\n", t3.etime(t2));
    }

    constexpr float CLIP_VAL = 65535.f;

    t3.set();

    if (ri->getSensorType() == ST_BAYER) {
#ifdef __SSE2__
        const vfloat onev = F2V(1.f);
        const vfloat clipv = F2V(CLIP_VAL);
#endif

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int row = 0; row < H; ++row) {
            int col = 0;
            // Avoid trouble with zeroes, minimum pixel value is 1.
            const float exps0 = exps[FC(row, col)];
            const float exps1 = exps[FC(row, col + 1)];
            const float mult0 = mults[FC(row, col)];
            const float mult1 = mults[FC(row, col + 1)];
#ifdef __SSE2__
            const vfloat expsv = _mm_setr_ps(exps0, exps1, exps0, exps1);
            const vfloat multsv = _mm_setr_ps(mult0, mult1, mult0, mult1);
            for (; col < W - 3; col += 4) {
                STVFU(rawData[row][col], vminf(multsv * pow_F(vmaxf(LVFU(rawData[row][col]), onev), expsv), clipv));
            }
#endif // __SSE2__
            for (; col < W - 1; col += 2) {
                rawData[row][col] = rtengine::min(mult0 * pow_F(rtengine::max(rawData[row][col], 1.f), exps0), CLIP_VAL);
                rawData[row][col + 1] = rtengine::min(mult1 * pow_F(rtengine::max(rawData[row][col + 1], 1.f), exps1), CLIP_VAL);
            }
            if (col < W) {
                rawData[row][col] = rtengine::min(mult0 * pow_F(rtengine::max(rawData[row][col], 1.f), exps0), CLIP_VAL);
            }
        }
    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef __SSE2__
        const vfloat onev = F2V(1.f);
        const vfloat clipv = F2V(CLIP_VAL);
#endif

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int row = 0; row < H; row ++) {
            int col = 0;
            // Avoid trouble with zeroes, minimum pixel value is 1.
            const std::array<float, 6> expsc = {
                exps[ri->XTRANSFC(row, 0)],
                exps[ri->XTRANSFC(row, 1)],
                exps[ri->XTRANSFC(row, 2)],
                exps[ri->XTRANSFC(row, 3)],
                exps[ri->XTRANSFC(row, 4)],
                exps[ri->XTRANSFC(row, 5)]
            };
            const std::array<float, 6> multsc = {
                mults[ri->XTRANSFC(row, 0)],
                mults[ri->XTRANSFC(row, 1)],
                mults[ri->XTRANSFC(row, 2)],
                mults[ri->XTRANSFC(row, 3)],
                mults[ri->XTRANSFC(row, 4)],
                mults[ri->XTRANSFC(row, 5)]
            };
#ifdef __SSE2__
            const vfloat expsv0 = _mm_setr_ps(expsc[0], expsc[1], expsc[2], expsc[3]);
            const vfloat expsv1 = _mm_setr_ps(expsc[4], expsc[5], expsc[0], expsc[1]);
            const vfloat expsv2 = _mm_setr_ps(expsc[2], expsc[3], expsc[4], expsc[5]);
            const vfloat multsv0 = _mm_setr_ps(multsc[0], multsc[1], multsc[2], multsc[3]);
            const vfloat multsv1 = _mm_setr_ps(multsc[4], multsc[5], multsc[0], multsc[1]);
            const vfloat multsv2 = _mm_setr_ps(multsc[2], multsc[3], multsc[4], multsc[5]);
            for (; col < W - 11; col += 12) {
                STVFU(rawData[row][col], vminf(multsv0 * pow_F(vmaxf(LVFU(rawData[row][col]), onev), expsv0), clipv));
                STVFU(rawData[row][col + 4], vminf(multsv1 * pow_F(vmaxf(LVFU(rawData[row][col + 4]), onev), expsv1), clipv));
                STVFU(rawData[row][col + 8], vminf(multsv2 * pow_F(vmaxf(LVFU(rawData[row][col + 8]), onev), expsv2), clipv));
            }
#endif // __SSE2__
            for (; col < W - 5; col += 6) {
                for (int c = 0; c < 6; ++c) {
                    rawData[row][col + c] = rtengine::min(multsc[c] * pow_F(rtengine::max(rawData[row][col + c], 1.f), expsc[c]), CLIP_VAL);
                }
            }
            for (int c = 0; col < W; col++, c++) {
                rawData[row][col + c] = rtengine::min(multsc[c] * pow_F(rtengine::max(rawData[row][col + c], 1.f), expsc[c]), CLIP_VAL);
            }
        }
    }

    t4.set();

    if (settings->verbose) {
        printf("Pow loop time us: %d\n", t4.etime(t3));
    }

    t4.set();

    PixelsMap bitmapBads(W, H);

    int totBP = 0; // Hold count of bad pixels to correct

    if (ri->getSensorType() == ST_BAYER) {
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:totBP) schedule(dynamic,16)
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                if (rawData[i][j] >= MAX_OUT_VALUE) {
                    bitmapBads.set(j, i);
                    ++totBP;
                }
            }
        }

        if (totBP > 0) {
            interpolateBadPixelsBayer(bitmapBads, rawData);
        }

    }
    else if (ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:totBP) schedule(dynamic,16)
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                if (rawData[i][j] >= MAX_OUT_VALUE) {
                    bitmapBads.set(j, i);
                    totBP++;
                }
            }
        }

        if (totBP > 0) {
            interpolateBadPixelsXtrans(bitmapBads);
        }
    }

    t5.set();

    if (settings->verbose) {
        printf("Bad pixels count: %d\n", totBP);
        printf("Bad pixels interpolation time us: %d\n", t5.etime(t4));
    }
}



// For backwards compatibility with profiles saved by RT 5.7 - 5.8
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

    const float MAX_OUT_VALUE = 65000.f;

    // Channel medians
    float rmed, gmed, bmed;

    // If using the old method, calculate medians on the whole image
    calcMedians(baseImg, 0, 0, rwidth, rheight, rmed, gmed, bmed);

    if (settings->verbose) {
        printf("Thumbnail input channel medians: %g %g %g\n", rmed, gmed, bmed);
    }

    // Calculate output medians
    rmed = powf(rmed, rexp);
    gmed = powf(gmed, gexp);
    bmed = powf(bmed, bexp);

    // Calculate output multipliers so that the median value is 1/24 of the output range.
    float rmult, gmult, bmult;
    rmult = (MAX_OUT_VALUE / 24.f) / rmed;
    gmult = (MAX_OUT_VALUE / 24.f) / gmed;
    bmult = (MAX_OUT_VALUE / 24.f) / bmed;


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



// ----------------- <<< legacy mode <<< ------------


