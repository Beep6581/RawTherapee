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
#include <cassert>

#include "../rawimagesource.h"
#include "../rawimagesource_i.h"
#include "../jaggedarray.h"
#include "../rt_math.h"
#include "../../rtgui/multilangmgr.h"
#include "../procparams.h"
#include "../opthelper.h"
//#define BENCHMARK
#include "../StopWatch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace rtengine
{

#undef ABS

#define ABS(a) ((a)<0?-(a):(a))

void RawImageSource::hphd_vertical (float** hpmap, int col_from, int col_to)
{
    float* temp = new float[max(W, H)];
    float* avg = new float[max(W, H)];
    float* dev = new float[max(W, H)];

    memset (temp, 0, max(W, H)*sizeof(float));
    memset (avg, 0, max(W, H)*sizeof(float));
    memset (dev, 0, max(W, H)*sizeof(float));

    for (int k = col_from; k < col_to; k++) {
        for (int i = 5; i < H - 5; i++) {
            temp[i] = (rawData[i - 5][k] - 8 * rawData[i - 4][k] + 27 * rawData[i - 3][k] - 48 * rawData[i - 2][k] + 42 * rawData[i - 1][k] -
                       (rawData[i + 5][k] - 8 * rawData[i + 4][k] + 27 * rawData[i + 3][k] - 48 * rawData[i + 2][k] + 42 * rawData[i + 1][k])) / 100.0;
            temp[i] = ABS(temp[i]);
        }

        for (int j = 4; j < H - 4; j++) {
            float avgL = (temp[j - 4] + temp[j - 3] + temp[j - 2] + temp[j - 1] + temp[j] + temp[j + 1] + temp[j + 2] + temp[j + 3] + temp[j + 4]) / 9.0;
            avg[j] = avgL;
            float devL = ((temp[j - 4] - avgL) * (temp[j - 4] - avgL) + (temp[j - 3] - avgL) * (temp[j - 3] - avgL) + (temp[j - 2] - avgL) * (temp[j - 2] - avgL) + (temp[j - 1] - avgL) * (temp[j - 1] - avgL) + (temp[j] - avgL) * (temp[j] - avgL) + (temp[j + 1] - avgL) * (temp[j + 1] - avgL) + (temp[j + 2] - avgL) * (temp[j + 2] - avgL) + (temp[j + 3] - avgL) * (temp[j + 3] - avgL) + (temp[j + 4] - avgL) * (temp[j + 4] - avgL)) / 9.0;

            if (devL < 0.001) {
                devL = 0.001;
            }

            dev[j] = devL;
        }

        for (int j = 5; j < H - 5; j++) {
            float avgL = avg[j - 1];
            float avgR = avg[j + 1];
            float devL = dev[j - 1];
            float devR = dev[j + 1];
            hpmap[j][k] = avgL + (avgR - avgL) * devL / (devL + devR);
        }
    }

    delete [] temp;
    delete [] avg;
    delete [] dev;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hphd_horizontal (float** hpmap, int row_from, int row_to)
{
    float* temp = new float[max(W, H)];
    float* avg = new float[max(W, H)];
    float* dev = new float[max(W, H)];

    memset (temp, 0, max(W, H)*sizeof(float));
    memset (avg, 0, max(W, H)*sizeof(float));
    memset (dev, 0, max(W, H)*sizeof(float));

    for (int i = row_from; i < row_to; i++) {
        for (int j = 5; j < W - 5; j++) {
            temp[j] = (rawData[i][j - 5] - 8 * rawData[i][j - 4] + 27 * rawData[i][j - 3] - 48 * rawData[i][j - 2] + 42 * rawData[i][j - 1] -
                       (rawData[i][j + 5] - 8 * rawData[i][j + 4] + 27 * rawData[i][j + 3] - 48 * rawData[i][j + 2] + 42 * rawData[i][j + 1])) / 100;
            temp[j] = ABS(temp[j]);
        }

        for (int j = 4; j < W - 4; j++) {
            float avgL = (temp[j - 4] + temp[j - 3] + temp[j - 2] + temp[j - 1] + temp[j] + temp[j + 1] + temp[j + 2] + temp[j + 3] + temp[j + 4]) / 9.0;
            avg[j] = avgL;
            float devL = ((temp[j - 4] - avgL) * (temp[j - 4] - avgL) + (temp[j - 3] - avgL) * (temp[j - 3] - avgL) + (temp[j - 2] - avgL) * (temp[j - 2] - avgL) + (temp[j - 1] - avgL) * (temp[j - 1] - avgL) + (temp[j] - avgL) * (temp[j] - avgL) + (temp[j + 1] - avgL) * (temp[j + 1] - avgL) + (temp[j + 2] - avgL) * (temp[j + 2] - avgL) + (temp[j + 3] - avgL) * (temp[j + 3] - avgL) + (temp[j + 4] - avgL) * (temp[j + 4] - avgL)) / 9.0;

            if (devL < 0.001) {
                devL = 0.001;
            }

            dev[j] = devL;
        }

        for (int j = 5; j < W - 5; j++) {
            float avgL = avg[j - 1];
            float avgR = avg[j + 1];
            float devL = dev[j - 1];
            float devR = dev[j + 1];
            float hpv = avgL + (avgR - avgL) * devL / (devL + devR);

            if (hpmap[i][j] < 0.8 * hpv) {
                hpmap[i][j] = 2;
            } else if (hpv < 0.8 * hpmap[i][j]) {
                hpmap[i][j] = 1;
            } else {
                hpmap[i][j] = 0;
            }
        }
    }

    delete [] temp;
    delete [] avg;
    delete [] dev;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hphd_green (float** hpmap)
{
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 3; i < H - 3; i++) {
        for (int j = 3; j < W - 3; j++) {
            if (ri->ISGREEN(i, j)) {
                green[i][j] = rawData[i][j];
            } else {
                if (hpmap[i][j] == 1) {
                    int g2 = rawData[i][j + 1] + ((rawData[i][j] - rawData[i][j + 2]) / 2);
                    int g4 = rawData[i][j - 1] + ((rawData[i][j] - rawData[i][j - 2]) / 2);

                    int dx = rawData[i][j + 1] - rawData[i][j - 1];
                    int d1 = rawData[i][j + 3] - rawData[i][j + 1];
                    int d2 = rawData[i][j + 2] - rawData[i][j];
                    int d3 = (rawData[i - 1][j + 2] - rawData[i - 1][j]) / 2;
                    int d4 = (rawData[i + 1][j + 2] - rawData[i + 1][j]) / 2;

                    double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = rawData[i][j - 3] - rawData[i][j - 1];
                    d2 = rawData[i][j - 2] - rawData[i][j];
                    d3 = (rawData[i - 1][j - 2] - rawData[i - 1][j]) / 2;
                    d4 = (rawData[i + 1][j - 2] - rawData[i + 1][j]) / 2;

                    double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    green[i][j] = (e2 * g2 + e4 * g4) / (e2 + e4);
                } else if (hpmap[i][j] == 2) {
                    int g1 = rawData[i - 1][j] + ((rawData[i][j] - rawData[i - 2][j]) / 2);
                    int g3 = rawData[i + 1][j] + ((rawData[i][j] - rawData[i + 2][j]) / 2);

                    int dy = rawData[i + 1][j] - rawData[i - 1][j];
                    int d1 = rawData[i - 1][j] - rawData[i - 3][j];
                    int d2 = rawData[i][j] - rawData[i - 2][j];
                    int d3 = (rawData[i][j - 1] - rawData[i - 2][j - 1]) / 2;
                    int d4 = (rawData[i][j + 1] - rawData[i - 2][j + 1]) / 2;

                    double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = rawData[i + 1][j] - rawData[i + 3][j];
                    d2 = rawData[i][j] - rawData[i + 2][j];
                    d3 = (rawData[i][j - 1] - rawData[i + 2][j - 1]) / 2;
                    d4 = (rawData[i][j + 1] - rawData[i + 2][j + 1]) / 2;

                    double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    green[i][j] = (e1 * g1 + e3 * g3) / (e1 + e3);
                } else {
                    int g1 = rawData[i - 1][j] + ((rawData[i][j] - rawData[i - 2][j]) / 2);
                    int g2 = rawData[i][j + 1] + ((rawData[i][j] - rawData[i][j + 2]) / 2);
                    int g3 = rawData[i + 1][j] + ((rawData[i][j] - rawData[i + 2][j]) / 2);
                    int g4 = rawData[i][j - 1] + ((rawData[i][j] - rawData[i][j - 2]) / 2);

                    int dx = rawData[i][j + 1] - rawData[i][j - 1];
                    int dy = rawData[i + 1][j] - rawData[i - 1][j];

                    int d1 = rawData[i - 1][j] - rawData[i - 3][j];
                    int d2 = rawData[i][j] - rawData[i - 2][j];
                    int d3 = (rawData[i][j - 1] - rawData[i - 2][j - 1]) / 2;
                    int d4 = (rawData[i][j + 1] - rawData[i - 2][j + 1]) / 2;

                    double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = rawData[i][j + 3] - rawData[i][j + 1];
                    d2 = rawData[i][j + 2] - rawData[i][j];
                    d3 = (rawData[i - 1][j + 2] - rawData[i - 1][j]) / 2;
                    d4 = (rawData[i + 1][j + 2] - rawData[i + 1][j]) / 2;

                    double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = rawData[i + 1][j] - rawData[i + 3][j];
                    d2 = rawData[i][j] - rawData[i + 2][j];
                    d3 = (rawData[i][j - 1] - rawData[i + 2][j - 1]) / 2;
                    d4 = (rawData[i][j + 1] - rawData[i + 2][j + 1]) / 2;

                    double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = rawData[i][j - 3] - rawData[i][j - 1];
                    d2 = rawData[i][j - 2] - rawData[i][j];
                    d3 = (rawData[i - 1][j - 2] - rawData[i - 1][j]) / 2;
                    d4 = (rawData[i + 1][j - 2] - rawData[i + 1][j]) / 2;

                    double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    green[i][j] = (e1 * g1 + e2 * g2 + e3 * g3 + e4 * g4) / (e1 + e2 + e3 + e4);
                }
            }
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hphd_demosaic ()
{
    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::HPHD)));
        plistener->setProgress (0.0);
    }

    JaggedArray<float> hpmap (W, H, true);

#ifdef _OPENMP
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int blk = W / nthreads;

        if (tid < nthreads - 1) {
            hphd_vertical (hpmap, tid * blk, (tid + 1)*blk);
        } else {
            hphd_vertical (hpmap, tid * blk, W);
        }
    }
#else
    hphd_vertical (hpmap, 0, W);
#endif

    if (plistener) {
        plistener->setProgress (0.33);
    }

#ifdef _OPENMP
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int blk = H / nthreads;

        if (tid < nthreads - 1) {
            hphd_horizontal (hpmap, tid * blk, (tid + 1)*blk);
        } else {
            hphd_horizontal (hpmap, tid * blk, H);
        }
    }
#else
    hphd_horizontal (hpmap, 0, H);
#endif

    hphd_green (hpmap);

    if (plistener) {
        plistener->setProgress (0.66);
    }

    for (int i = 0; i < H; i++) {
        if (i == 0) {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], nullptr, green[i], green[i + 1], i, 1.0, 1.0, 1.0, 0, W, 1);
        } else if (i == H - 1) {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], green[i - 1], green[i], nullptr, i, 1.0, 1.0, 1.0, 0, W, 1);
        } else {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], green[i - 1], green[i], green[i + 1], i, 1.0, 1.0, 1.0, 0, W, 1);
        }
    }

    bayerborder_demosaic(W, H, 4, rawData, red, green, blue);

    if (plistener) {
        plistener->setProgress (1.0);
    }
}

} /* namespace */
