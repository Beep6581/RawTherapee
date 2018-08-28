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
#include "../rawimage.h"
#include "../mytime.h"
#include "../iccmatrices.h"
#include "../iccstore.h"
#include "../image8.h"
#include "../curves.h"
#include "../dfmanager.h"
#include "../slicer.h"
#include "../rt_math.h"
#include "../color.h"
#include "../../rtgui/multilangmgr.h"
#include "../procparams.h"
#include "../sleef.c"
#include "../opthelper.h"
#include "../median.h"
#define BENCHMARK
#include "../StopWatch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace rtengine
{

extern const Settings* settings;

#undef ABS

#define ABS(a) ((a)<0?-(a):(a))
#define CLIREF(x) LIM(x,-200000.0f,200000.0f) // avoid overflow : do not act directly on image[] or pix[]
#define x1125(a) (a + xdivf(a, 3))
#define x0875(a) (a - xdivf(a, 3))
#define x0250(a) xdivf(a, 2)
#define x00625(a) xdivf(a, 4)
#define x0125(a) xdivf(a, 3)

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

    border_interpolate2(W, H, 4, rawData, red, green, blue);

    if (plistener) {
        plistener->setProgress (1.0);
    }
}

#undef fc
#define fc(row,col) \
    (ri->get_filters() >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)

#define FORCC for (unsigned int c=0; c < colors; c++)

void RawImageSource::border_interpolate(unsigned int border, float (*image)[4], unsigned int start, unsigned int end)
{
    unsigned row, col, y, x, f;
    float sum[8];
    unsigned int width = W, height = H;
    unsigned int colors = 3;

    if (end == 0 ) {
        end = H;
    }

    for (row = start; row < end; row++)
        for (col = 0; col < width; col++) {
            if (col == border && row >= border && row < height - border) {
                col = width - border;
            }

            memset (sum, 0, sizeof sum);

            for (y = row - 1; y != row + 2; y++)
                for (x = col - 1; x != col + 2; x++)
                    if (y < height && x < width) {
                        f = fc(y, x);
                        sum[f] += image[y * width + x][f];
                        sum[f + 4]++;
                    }

            f = fc(row, col);

            FORCC if (c != f && sum[c + 4]) {
                image[row * width + col][c] = sum[c] / sum[c + 4];
            }
        }
}

void RawImageSource::border_interpolate2( int winw, int winh, int lborders, const array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue)
{
    int bord = lborders;
    int width = winw;
    int height = winh;

    for (int i = 0; i < height; i++) {

        float sum[6];

        for (int j = 0; j < bord; j++) { //first few columns
            for (int c = 0; c < 6; c++) {
                sum[c] = 0;
            }

            for (int i1 = i - 1; i1 < i + 2; i1++)
                for (int j1 = j - 1; j1 < j + 2; j1++) {
                    if ((i1 > -1) && (i1 < height) && (j1 > -1)) {
                        int c = FC(i1, j1);
                        sum[c] += rawData[i1][j1];
                        sum[c + 3]++;
                    }
                }

            int c = FC(i, j);

            if (c == 1) {
                red[i][j] = sum[0] / sum[3];
                green[i][j] = rawData[i][j];
                blue[i][j] = sum[2] / sum[5];
            } else {
                green[i][j] = sum[1] / sum[4];

                if (c == 0) {
                    red[i][j] = rawData[i][j];
                    blue[i][j] = sum[2] / sum[5];
                } else {
                    red[i][j] = sum[0] / sum[3];
                    blue[i][j] = rawData[i][j];
                }
            }
        }//j

        for (int j = width - bord; j < width; j++) { //last few columns
            for (int c = 0; c < 6; c++) {
                sum[c] = 0;
            }

            for (int i1 = i - 1; i1 < i + 2; i1++)
                for (int j1 = j - 1; j1 < j + 2; j1++) {
                    if ((i1 > -1) && (i1 < height ) && (j1 < width)) {
                        int c = FC(i1, j1);
                        sum[c] += rawData[i1][j1];
                        sum[c + 3]++;
                    }
                }

            int c = FC(i, j);

            if (c == 1) {
                red[i][j] = sum[0] / sum[3];
                green[i][j] = rawData[i][j];
                blue[i][j] = sum[2] / sum[5];
            } else {
                green[i][j] = sum[1] / sum[4];

                if (c == 0) {
                    red[i][j] = rawData[i][j];
                    blue[i][j] = sum[2] / sum[5];
                } else {
                    red[i][j] = sum[0] / sum[3];
                    blue[i][j] = rawData[i][j];
                }
            }
        }//j
    }//i

    for (int i = 0; i < bord; i++) {

        float sum[6];

        for (int j = bord; j < width - bord; j++) { //first few rows
            for (int c = 0; c < 6; c++) {
                sum[c] = 0;
            }

            for (int i1 = i - 1; i1 < i + 2; i1++)
                for (int j1 = j - 1; j1 < j + 2; j1++) {
                    if ((i1 > -1) && (i1 < height) && (j1 > -1)) {
                        int c = FC(i1, j1);
                        sum[c] += rawData[i1][j1];
                        sum[c + 3]++;
                    }
                }

            int c = FC(i, j);

            if (c == 1) {
                red[i][j] = sum[0] / sum[3];
                green[i][j] = rawData[i][j];
                blue[i][j] = sum[2] / sum[5];
            } else {
                green[i][j] = sum[1] / sum[4];

                if (c == 0) {
                    red[i][j] = rawData[i][j];
                    blue[i][j] = sum[2] / sum[5];
                } else {
                    red[i][j] = sum[0] / sum[3];
                    blue[i][j] = rawData[i][j];
                }
            }
        }//j
    }

    for (int i = height - bord; i < height; i++) {

        float sum[6];

        for (int j = bord; j < width - bord; j++) { //last few rows
            for (int c = 0; c < 6; c++) {
                sum[c] = 0;
            }

            for (int i1 = i - 1; i1 < i + 2; i1++)
                for (int j1 = j - 1; j1 < j + 2; j1++) {
                    if ((i1 > -1) && (i1 < height) && (j1 < width)) {
                        int c = FC(i1, j1);
                        sum[c] += rawData[i1][j1];
                        sum[c + 3]++;
                    }
                }

            int c = FC(i, j);

            if (c == 1) {
                red[i][j] = sum[0] / sum[3];
                green[i][j] = rawData[i][j];
                blue[i][j] = sum[2] / sum[5];
            } else {
                green[i][j] = sum[1] / sum[4];

                if (c == 0) {
                    red[i][j] = rawData[i][j];
                    blue[i][j] = sum[2] / sum[5];
                } else {
                    red[i][j] = sum[0] / sum[3];
                    blue[i][j] = rawData[i][j];
                }
            }
        }//j
    }

}

// LSMME demosaicing algorithm
// L. Zhang and X. Wu,
// Color demozaicing via directional Linear Minimum Mean Square-error Estimation,
// IEEE Trans. on Image Processing, vol. 14, pp. 2167-2178,
// Dec. 2005.
// Adapted to RawTherapee by Jacques Desmis 3/2013
// Improved speed and reduced memory consumption by Ingo Weyrich 2/2015
//TODO Tiles to reduce memory consumption
void RawImageSource::lmmse_interpolate_omp(int winw, int winh, array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue, int iterations)
{
    const int width = winw, height = winh;
    const int ba = 10;
    const int rr1 = height + 2 * ba;
    const int cc1 = width + 2 * ba;
    const int w1 = cc1;
    const int w2 = 2 * w1;
    const int w3 = 3 * w1;
    const int w4 = 4 * w1;
    float h0, h1, h2, h3, h4, hs;
    h0 = 1.0f;
    h1 = exp( -1.0f / 8.0f);
    h2 = exp( -4.0f / 8.0f);
    h3 = exp( -9.0f / 8.0f);
    h4 = exp(-16.0f / 8.0f);
    hs = h0 + 2.0f * (h1 + h2 + h3 + h4);
    h0 /= hs;
    h1 /= hs;
    h2 /= hs;
    h3 /= hs;
    h4 /= hs;
    int passref = 0;
    int iter = 0;

    if(iterations <= 4) {
        iter = iterations - 1;
        passref = 0;
    } else if (iterations <= 6) {
        iter = 3;
        passref = iterations - 4;
    } else if (iterations <= 8) {
        iter = 3;
        passref = iterations - 6;
    }

    bool applyGamma = true;

    if(iterations == 0) {
        applyGamma = false;
        iter = 0;
    } else {
        applyGamma = true;
    }

    float *rix[5];
    float *qix[5];
    float *buffer = (float *)calloc(rr1 * cc1 * 5 * sizeof(float), 1);

    if(buffer == nullptr) { // allocation of big block of memory failed, try to get 5 smaller ones
        printf("lmmse_interpolate_omp: allocation of big memory block failed, try to get 5 smaller ones now...\n");
        bool allocationFailed = false;

        for(int i = 0; i < 5; i++) {
            qix[i] = (float *)calloc(rr1 * cc1 * sizeof(float), 1);

            if(!qix[i]) { // allocation of at least one small block failed
                allocationFailed = true;
            }
        }

        if(allocationFailed) { // fall back to igv_interpolate
            printf("lmmse_interpolate_omp: allocation of 5 small memory blocks failed, falling back to igv_interpolate...\n");

            for(int i = 0; i < 5; i++) { // free the already allocated buffers
                if(qix[i]) {
                    free(qix[i]);
                }
            }

            igv_interpolate(winw, winh);
            return;
        }
    } else {
        qix[0] = buffer;

        for(int i = 1; i < 5; i++) {
            qix[i] = qix[i - 1] + rr1 * cc1;
        }
    }

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::LMMSE)));
        plistener->setProgress (0.0);
    }


    LUTf *gamtab;

    if(applyGamma) {
        gamtab = &(Color::gammatab_24_17a);
    } else {
        gamtab = new LUTf(65536, LUT_CLIP_ABOVE | LUT_CLIP_BELOW);
        gamtab->makeIdentity(65535.f);
    }


#ifdef _OPENMP
    #pragma omp parallel private(rix)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rrr = ba; rrr < rr1 - ba; rrr++) {
            for (int ccc = ba, row = rrr - ba; ccc < cc1 - ba; ccc++) {
                int col = ccc - ba;
                float *rix = qix[4] + rrr * cc1 + ccc;
                rix[0] = (*gamtab)[rawData[row][col]];
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.1);
            }
        }

        // G-R(B)
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int rr = 2; rr < rr1 - 2; rr++) {
            // G-R(B) at R(B) location
            for (int cc = 2 + (FC(rr, 2) & 1); cc < cc1 - 2; cc += 2) {
                rix[4] = qix[4] + rr * cc1 + cc;
                float v0 = x00625(rix[4][-w1 - 1] + rix[4][-w1 + 1] + rix[4][w1 - 1] + rix[4][w1 + 1]) + x0250(rix[4][0]);
                // horizontal
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[0][0] = - x0250(rix[4][ -2] + rix[4][ 2]) + xdiv2f(rix[4][ -1] + rix[4][0] + rix[4][ 1]);
                float Y = v0 + xdiv2f(rix[0][0]);

                if (rix[4][0] > 1.75f * Y) {
                    rix[0][0] = median(rix[0][0], rix[4][ -1], rix[4][ 1]);
                } else {
                    rix[0][0] = LIM(rix[0][0], 0.0f, 1.0f);
                }

                rix[0][0] -= rix[4][0];
                // vertical
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[1][0] = -x0250(rix[4][-w2] + rix[4][w2]) + xdiv2f(rix[4][-w1] + rix[4][0] + rix[4][w1]);
                Y = v0 + xdiv2f(rix[1][0]);

                if (rix[4][0] > 1.75f * Y) {
                    rix[1][0] = median(rix[1][0], rix[4][-w1], rix[4][w1]);
                } else {
                    rix[1][0] = LIM(rix[1][0], 0.0f, 1.0f);
                }

                rix[1][0] -= rix[4][0];
            }

            // G-R(B) at G location
            for (int ccc = 2 + (FC(rr, 3) & 1); ccc < cc1 - 2; ccc += 2) {
                rix[0] = qix[0] + rr * cc1 + ccc;
                rix[1] = qix[1] + rr * cc1 + ccc;
                rix[4] = qix[4] + rr * cc1 + ccc;
                rix[0][0] = x0250(rix[4][ -2] + rix[4][ 2]) - xdiv2f(rix[4][ -1] + rix[4][0] + rix[4][ 1]);
                rix[1][0] = x0250(rix[4][-w2] + rix[4][w2]) - xdiv2f(rix[4][-w1] + rix[4][0] + rix[4][w1]);
                rix[0][0] = LIM(rix[0][0], -1.0f, 0.0f) + rix[4][0];
                rix[1][0] = LIM(rix[1][0], -1.0f, 0.0f) + rix[4][0];
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.2);
            }
        }


        // apply low pass filter on differential colors
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 4; rr < rr1 - 4; rr++)
            for (int cc = 4; cc < cc1 - 4; cc++) {
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[2] = qix[2] + rr * cc1 + cc;
                rix[2][0] = h0 * rix[0][0] + h1 * (rix[0][ -1] + rix[0][ 1]) + h2 * (rix[0][ -2] + rix[0][ 2]) + h3 * (rix[0][ -3] + rix[0][ 3]) + h4 * (rix[0][ -4] + rix[0][ 4]);
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[3] = qix[3] + rr * cc1 + cc;
                rix[3][0] = h0 * rix[1][0] + h1 * (rix[1][-w1] + rix[1][w1]) + h2 * (rix[1][-w2] + rix[1][w2]) + h3 * (rix[1][-w3] + rix[1][w3]) + h4 * (rix[1][-w4] + rix[1][w4]);
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.3);
            }
        }

        // interpolate G-R(B) at R(B)
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 4; rr < rr1 - 4; rr++) {
            int cc = 4 + (FC(rr, 4) & 1);
#ifdef __SSE2__
            __m128 p1v, p2v, p3v, p4v, p5v, p6v, p7v, p8v, p9v, muv, vxv, vnv, xhv, vhv, xvv, vvv;
            __m128 epsv = _mm_set1_ps(1e-7);
            __m128 ninev = _mm_set1_ps(9.f);

            for (; cc < cc1 - 10; cc += 8) {
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[2] = qix[2] + rr * cc1 + cc;
                rix[3] = qix[3] + rr * cc1 + cc;
                rix[4] = qix[4] + rr * cc1 + cc;
                // horizontal
                p1v = LC2VFU(rix[2][-4]);
                p2v = LC2VFU(rix[2][-3]);
                p3v = LC2VFU(rix[2][-2]);
                p4v = LC2VFU(rix[2][-1]);
                p5v = LC2VFU(rix[2][ 0]);
                p6v = LC2VFU(rix[2][ 1]);
                p7v = LC2VFU(rix[2][ 2]);
                p8v = LC2VFU(rix[2][ 3]);
                p9v = LC2VFU(rix[2][ 4]);
                muv = (p1v + p2v + p3v + p4v + p5v + p6v + p7v + p8v + p9v) / ninev;
                vxv = epsv + SQRV(p1v - muv) + SQRV(p2v - muv) + SQRV(p3v - muv) + SQRV(p4v - muv) + SQRV(p5v - muv) + SQRV(p6v - muv) + SQRV(p7v - muv) + SQRV(p8v - muv) + SQRV(p9v - muv);
                p1v -= LC2VFU(rix[0][-4]);
                p2v -= LC2VFU(rix[0][-3]);
                p3v -= LC2VFU(rix[0][-2]);
                p4v -= LC2VFU(rix[0][-1]);
                p5v -= LC2VFU(rix[0][ 0]);
                p6v -= LC2VFU(rix[0][ 1]);
                p7v -= LC2VFU(rix[0][ 2]);
                p8v -= LC2VFU(rix[0][ 3]);
                p9v -= LC2VFU(rix[0][ 4]);
                vnv = epsv + SQRV(p1v) + SQRV(p2v) + SQRV(p3v) + SQRV(p4v) + SQRV(p5v) + SQRV(p6v) + SQRV(p7v) + SQRV(p8v) + SQRV(p9v);
                xhv = (LC2VFU(rix[0][0]) * vxv + LC2VFU(rix[2][0]) * vnv) / (vxv + vnv);
                vhv = vxv * vnv / (vxv + vnv);

                // vertical
                p1v = LC2VFU(rix[3][-w4]);
                p2v = LC2VFU(rix[3][-w3]);
                p3v = LC2VFU(rix[3][-w2]);
                p4v = LC2VFU(rix[3][-w1]);
                p5v = LC2VFU(rix[3][  0]);
                p6v = LC2VFU(rix[3][ w1]);
                p7v = LC2VFU(rix[3][ w2]);
                p8v = LC2VFU(rix[3][ w3]);
                p9v = LC2VFU(rix[3][ w4]);
                muv = (p1v + p2v + p3v + p4v + p5v + p6v + p7v + p8v + p9v) / ninev;
                vxv = epsv + SQRV(p1v - muv) + SQRV(p2v - muv) + SQRV(p3v - muv) + SQRV(p4v - muv) + SQRV(p5v - muv) + SQRV(p6v - muv) + SQRV(p7v - muv) + SQRV(p8v - muv) + SQRV(p9v - muv);
                p1v -= LC2VFU(rix[1][-w4]);
                p2v -= LC2VFU(rix[1][-w3]);
                p3v -= LC2VFU(rix[1][-w2]);
                p4v -= LC2VFU(rix[1][-w1]);
                p5v -= LC2VFU(rix[1][  0]);
                p6v -= LC2VFU(rix[1][ w1]);
                p7v -= LC2VFU(rix[1][ w2]);
                p8v -= LC2VFU(rix[1][ w3]);
                p9v -= LC2VFU(rix[1][ w4]);
                vnv = epsv + SQRV(p1v) + SQRV(p2v) + SQRV(p3v) + SQRV(p4v) + SQRV(p5v) + SQRV(p6v) + SQRV(p7v) + SQRV(p8v) + SQRV(p9v);
                xvv = (LC2VFU(rix[1][0]) * vxv + LC2VFU(rix[3][0]) * vnv) / (vxv + vnv);
                vvv = vxv * vnv / (vxv + vnv);
                // interpolated G-R(B)
                muv = (xhv * vvv + xvv * vhv) / (vhv + vvv);
                STC2VFU(rix[4][0], muv);
            }

#endif

            for (; cc < cc1 - 4; cc += 2) {
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[2] = qix[2] + rr * cc1 + cc;
                rix[3] = qix[3] + rr * cc1 + cc;
                rix[4] = qix[4] + rr * cc1 + cc;
                // horizontal
                float p1 = rix[2][-4];
                float p2 = rix[2][-3];
                float p3 = rix[2][-2];
                float p4 = rix[2][-1];
                float p5 = rix[2][ 0];
                float p6 = rix[2][ 1];
                float p7 = rix[2][ 2];
                float p8 = rix[2][ 3];
                float p9 = rix[2][ 4];
                float mu = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) / 9.f;
                float vx = 1e-7 + SQR(p1 - mu) + SQR(p2 - mu) + SQR(p3 - mu) + SQR(p4 - mu) + SQR(p5 - mu) + SQR(p6 - mu) + SQR(p7 - mu) + SQR(p8 - mu) + SQR(p9 - mu);
                p1 -= rix[0][-4];
                p2 -= rix[0][-3];
                p3 -= rix[0][-2];
                p4 -= rix[0][-1];
                p5 -= rix[0][ 0];
                p6 -= rix[0][ 1];
                p7 -= rix[0][ 2];
                p8 -= rix[0][ 3];
                p9 -= rix[0][ 4];
                float vn = 1e-7 + SQR(p1) + SQR(p2) + SQR(p3) + SQR(p4) + SQR(p5) + SQR(p6) + SQR(p7) + SQR(p8) + SQR(p9);
                float xh = (rix[0][0] * vx + rix[2][0] * vn) / (vx + vn);
                float vh = vx * vn / (vx + vn);

                // vertical
                p1 = rix[3][-w4];
                p2 = rix[3][-w3];
                p3 = rix[3][-w2];
                p4 = rix[3][-w1];
                p5 = rix[3][  0];
                p6 = rix[3][ w1];
                p7 = rix[3][ w2];
                p8 = rix[3][ w3];
                p9 = rix[3][ w4];
                mu = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) / 9.f;
                vx = 1e-7 + SQR(p1 - mu) + SQR(p2 - mu) + SQR(p3 - mu) + SQR(p4 - mu) + SQR(p5 - mu) + SQR(p6 - mu) + SQR(p7 - mu) + SQR(p8 - mu) + SQR(p9 - mu);
                p1 -= rix[1][-w4];
                p2 -= rix[1][-w3];
                p3 -= rix[1][-w2];
                p4 -= rix[1][-w1];
                p5 -= rix[1][  0];
                p6 -= rix[1][ w1];
                p7 -= rix[1][ w2];
                p8 -= rix[1][ w3];
                p9 -= rix[1][ w4];
                vn = 1e-7 + SQR(p1) + SQR(p2) + SQR(p3) + SQR(p4) + SQR(p5) + SQR(p6) + SQR(p7) + SQR(p8) + SQR(p9);
                float xv = (rix[1][0] * vx + rix[3][0] * vn) / (vx + vn);
                float vv = vx * vn / (vx + vn);
                // interpolated G-R(B)
                rix[4][0] = (xh * vv + xv * vh) / (vh + vv);
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.4);
            }
        }

        // copy CFA values
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 0; rr < rr1; rr++)
            for (int cc = 0, row = rr - ba; cc < cc1; cc++) {
                int col = cc - ba;
                int c = FC(rr, cc);
                rix[c] = qix[c] + rr * cc1 + cc;

                if ((row >= 0) & (row < height) & (col >= 0) & (col < width)) {
                    rix[c][0] = (*gamtab)[rawData[row][col]];
                } else {
                    rix[c][0] = 0.f;
                }

                if (c != 1) {
                    rix[1] = qix[1] + rr * cc1 + cc;
                    rix[4] = qix[4] + rr * cc1 + cc;
                    rix[1][0] = rix[c][0] + rix[4][0];
                }
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.5);
            }
        }

        // bilinear interpolation for R/B
        // interpolate R/B at G location
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 1; rr < rr1 - 1; rr++)
            for (int cc = 1 + (FC(rr, 2) & 1), c = FC(rr, cc + 1); cc < cc1 - 1; cc += 2) {
                rix[c] = qix[c] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[c][0] = rix[1][0] + xdiv2f(rix[c][ -1] - rix[1][ -1] + rix[c][ 1] - rix[1][ 1]);
                c = 2 - c;
                rix[c] = qix[c] + rr * cc1 + cc;
                rix[c][0] = rix[1][0] + xdiv2f(rix[c][-w1] - rix[1][-w1] + rix[c][w1] - rix[1][w1]);
                c = 2 - c;
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.6);
            }
        }

        // interpolate R/B at B/R location
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 1; rr < rr1 - 1; rr++)
            for (int cc = 1 + (FC(rr, 1) & 1), c = 2 - FC(rr, cc); cc < cc1 - 1; cc += 2) {
                rix[c] = qix[c] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[c][0] = rix[1][0] + x0250(rix[c][-w1] - rix[1][-w1] + rix[c][ -1] - rix[1][ -1] + rix[c][  1] - rix[1][  1] + rix[c][ w1] - rix[1][ w1]);
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.7);
            }
        }

    }// End of parallelization 1

    // median filter/
    for (int pass = 0; pass < iter; pass++) {
        // Apply 3x3 median filter
        // Compute median(R-G) and median(B-G)

#ifdef _OPENMP
        #pragma omp parallel for private(rix)
#endif

        for (int rr = 1; rr < rr1 - 1; rr++) {
            for (int c = 0; c < 3; c += 2) {
                int d = c + 3 - (c == 0 ? 0 : 1);
                int cc = 1;
#ifdef __SSE2__

                for (; cc < cc1 - 4; cc += 4) {
                    rix[d] = qix[d] + rr * cc1 + cc;
                    rix[c] = qix[c] + rr * cc1 + cc;
                    rix[1] = qix[1] + rr * cc1 + cc;
                    // Assign 3x3 differential color values
                    const std::array<vfloat, 9> p = {
                        LVFU(rix[c][-w1 - 1]) - LVFU(rix[1][-w1 - 1]),
                        LVFU(rix[c][-w1]) - LVFU(rix[1][-w1]),
                        LVFU(rix[c][-w1 + 1]) - LVFU(rix[1][-w1 + 1]),
                        LVFU(rix[c][   -1]) - LVFU(rix[1][   -1]),
                        LVFU(rix[c][  0]) - LVFU(rix[1][  0]),
                        LVFU(rix[c][    1]) - LVFU(rix[1][    1]),
                        LVFU(rix[c][ w1 - 1]) - LVFU(rix[1][ w1 - 1]),
                        LVFU(rix[c][ w1]) - LVFU(rix[1][ w1]),
                        LVFU(rix[c][ w1 + 1]) - LVFU(rix[1][ w1 + 1])
                    };
                    _mm_storeu_ps(&rix[d][0], median(p));
                }

#endif

                for (; cc < cc1 - 1; cc++) {
                    rix[d] = qix[d] + rr * cc1 + cc;
                    rix[c] = qix[c] + rr * cc1 + cc;
                    rix[1] = qix[1] + rr * cc1 + cc;
                    // Assign 3x3 differential color values
                    const std::array<float, 9> p = {
                        rix[c][-w1 - 1] - rix[1][-w1 - 1],
                        rix[c][-w1] - rix[1][-w1],
                        rix[c][-w1 + 1] - rix[1][-w1 + 1],
                        rix[c][   -1] - rix[1][   -1],
                        rix[c][  0] - rix[1][  0],
                        rix[c][    1] - rix[1][    1],
                        rix[c][ w1 - 1] - rix[1][ w1 - 1],
                        rix[c][ w1] - rix[1][ w1],
                        rix[c][ w1 + 1] - rix[1][ w1 + 1]
                    };
                    rix[d][0] = median(p);
                }
            }
        }

        // red/blue at GREEN pixel locations & red/blue and green at BLUE/RED pixel locations
#ifdef _OPENMP
        #pragma omp parallel for private (rix)
#endif

        for (int rr = 0; rr < rr1; rr++) {
            rix[0] = qix[0] + rr * cc1;
            rix[1] = qix[1] + rr * cc1;
            rix[2] = qix[2] + rr * cc1;
            rix[3] = qix[3] + rr * cc1;
            rix[4] = qix[4] + rr * cc1;
            int c0 = FC(rr, 0);
            int c1 = FC(rr, 1);

            if(c0 == 1) {
                c1 = 2 - c1;
                int d = c1 + 3 - (c1 == 0 ? 0 : 1);
                int cc;

                for (cc = 0; cc < cc1 - 1; cc += 2) {
                    rix[0][0] = rix[1][0] + rix[3][0];
                    rix[2][0] = rix[1][0] + rix[4][0];
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                    rix[c1][0] = rix[1][0] + rix[d][0];
                    rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                }

                if(cc < cc1) { // remaining pixel, only if width is odd
                    rix[0][0] = rix[1][0] + rix[3][0];
                    rix[2][0] = rix[1][0] + rix[4][0];
                }
            } else {
                c0 = 2 - c0;
                int d = c0 + 3 - (c0 == 0 ? 0 : 1);
                int cc;

                for (cc = 0; cc < cc1 - 1; cc += 2) {
                    rix[c0][0] = rix[1][0] + rix[d][0];
                    rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                    rix[0][0] = rix[1][0] + rix[3][0];
                    rix[2][0] = rix[1][0] + rix[4][0];
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                }

                if(cc < cc1) { // remaining pixel, only if width is odd
                    rix[c0][0] = rix[1][0] + rix[d][0];
                    rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
                }
            }
        }
    }

    if (plistener) {
        plistener->setProgress (0.8);
    }

    if(applyGamma) {
        gamtab = &(Color::igammatab_24_17);
    } else {
        gamtab->makeIdentity();
    }

    array2D<float>* rgb[3];
    rgb[0] = &red;
    rgb[1] = &green;
    rgb[2] = &blue;

    // copy result back to image matrix
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int row = 0; row < height; row++) {
        for (int col = 0, rr = row + ba; col < width; col++) {
            int cc = col + ba;
            int c = FC(row, col);

            for (int ii = 0; ii < 3; ii++)
                if (ii != c) {
                    float *rix = qix[ii] + rr * cc1 + cc;
                    (*(rgb[ii]))[row][col] = (*gamtab)[65535.f * rix[0]];
                } else {
                    (*(rgb[ii]))[row][col] = CLIP(rawData[row][col]);
                }
        }
    }

    if (plistener) {
        plistener->setProgress (1.0);
    }

    if(buffer) {
        free(buffer);
    } else
        for(int i = 0; i < 5; i++) {
            free(qix[i]);
        }

    if(!applyGamma) {
        delete gamtab;
    }

    if(iterations > 4 && iterations <= 6) {
        refinement(passref);
    } else if(iterations > 6) {
        refinement_lassus(passref);
    }

}

/***
*
*   Bayer CFA Demosaicing using Integrated Gaussian Vector on Color Differences
*   Revision 1.0 - 2013/02/28
*
*   Copyright (c) 2007-2013 Luis Sanz Rodriguez
*   Using High Order Interpolation technique by Jim S, Jimmy Li, and Sharmil Randhawa
*
*   Contact info: luis.sanz.rodriguez@gmail.com
*
*   This code is distributed under a GNU General Public License, version 3.
*   Visit <http://www.gnu.org/licenses/> for more information.
*
***/
// Adapted to RawTherapee by Jacques Desmis 3/2013
// SSE version by Ingo Weyrich 5/2013
#ifdef __SSE2__
#define CLIPV(a) LIMV(a,zerov,c65535v)
void RawImageSource::igv_interpolate(int winw, int winh)
{
    static const float eps = 1e-5f, epssq = 1e-5f; //mod epssq -10f =>-5f Jacques 3/2013 to prevent artifact (divide by zero)

    static const int h1 = 1, h2 = 2, h3 = 3, h5 = 5;
    const int width = winw, height = winh;
    const int v1 = 1 * width, v2 = 2 * width, v3 = 3 * width, v5 = 5 * width;
    float* rgb[2];
    float* chr[4];
    float *rgbarray, *vdif, *hdif, *chrarray;
    rgbarray    = (float (*)) malloc((width * height) * sizeof( float ) );
    rgb[0] = rgbarray;
    rgb[1] = rgbarray + (width * height) / 2;

    vdif  = (float (*)) calloc( width * height / 2, sizeof * vdif );
    hdif  = (float (*)) calloc( width * height / 2, sizeof * hdif );

    chrarray    = (float (*)) calloc( width * height, sizeof( float ) );
    chr[0] = chrarray;
    chr[1] = chrarray + (width * height) / 2;

    // mapped chr[2] and chr[3] to hdif and hdif, because these are out of use, when chr[2] and chr[3] are used
    chr[2] = hdif;
    chr[3] = vdif;

    border_interpolate2(winw, winh, 7, rawData, red, green, blue);

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::IGV)));
        plistener->setProgress (0.0);
    }

#ifdef _OPENMP
    #pragma omp parallel default(none) shared(rgb,vdif,hdif,chr)
#endif
    {
        __m128 ngv, egv, wgv, sgv, nvv, evv, wvv, svv, nwgv, negv, swgv, segv, nwvv, nevv, swvv, sevv, tempv, temp1v, temp2v, temp3v, temp4v, temp5v, temp6v, temp7v, temp8v;
        __m128 epsv = _mm_set1_ps( eps );
        __m128 epssqv = _mm_set1_ps( epssq );
        __m128 c65535v = _mm_set1_ps( 65535.f );
        __m128 c23v = _mm_set1_ps( 23.f );
        __m128 c40v = _mm_set1_ps( 40.f );
        __m128 c51v = _mm_set1_ps( 51.f );
        __m128 c32v = _mm_set1_ps( 32.f );
        __m128 c8v = _mm_set1_ps( 8.f );
        __m128 c7v = _mm_set1_ps( 7.f );
        __m128 c6v = _mm_set1_ps( 6.f );
        __m128 c10v = _mm_set1_ps( 10.f );
        __m128 c21v = _mm_set1_ps( 21.f );
        __m128 c78v = _mm_set1_ps( 78.f );
        __m128 c69v = _mm_set1_ps( 69.f );
        __m128 c3145680v = _mm_set1_ps( 3145680.f );
        __m128 onev = _mm_set1_ps ( 1.f );
        __m128 zerov = _mm_set1_ps ( 0.f );
        __m128 d725v = _mm_set1_ps ( 0.725f );
        __m128 d1375v = _mm_set1_ps ( 0.1375f );

        float *dest1, *dest2;
        float ng, eg, wg, sg, nv, ev, wv, sv, nwg, neg, swg, seg, nwv, nev, swv, sev;
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 0; row < height - 0; row++) {
            dest1 = rgb[FC(row, 0) & 1];
            dest2 = rgb[FC(row, 1) & 1];
            int col, indx;

            for (col = 0, indx = row * width + col; col < width - 7; col += 8, indx += 8) {
                temp1v = LVFU( rawData[row][col] );
                temp1v = CLIPV( temp1v );
                temp2v = LVFU( rawData[row][col + 4] );
                temp2v = CLIPV( temp2v );
                tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 2, 0, 2, 0 ) );
                _mm_storeu_ps( &dest1[indx >> 1], tempv );
                tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 3, 1, 3, 1 ) );
                _mm_storeu_ps( &dest2[indx >> 1], tempv );
            }

            for (; col < width; col++, indx += 2) {
                dest1[indx >> 1] = CLIP(rawData[row][col]); //rawData = RT datas
                col++;
                if(col < width)
                    dest2[indx >> 1] = CLIP(rawData[row][col]); //rawData = RT datas
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.13);
            }
        }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 5; row < height - 5; row++) {
            int col, indx, indx1;

            for (col = 5 + (FC(row, 1) & 1), indx = row * width + col, indx1 = indx >> 1; col < width - 12; col += 8, indx += 8, indx1 += 4) {
                //N,E,W,S Gradients
                ngv = (epsv + (vabsf(LVFU(rgb[1][(indx - v1) >> 1]) - LVFU(rgb[1][(indx - v3) >> 1])) + vabsf(LVFU(rgb[0][indx1]) - LVFU(rgb[0][(indx1 - v1)]))) / c65535v);
                egv = (epsv + (vabsf(LVFU(rgb[1][(indx + h1) >> 1]) - LVFU(rgb[1][(indx + h3) >> 1])) + vabsf(LVFU(rgb[0][indx1]) - LVFU(rgb[0][(indx1 + h1)]))) / c65535v);
                wgv = (epsv + (vabsf(LVFU(rgb[1][(indx - h1) >> 1]) - LVFU(rgb[1][(indx - h3) >> 1])) + vabsf(LVFU(rgb[0][indx1]) - LVFU(rgb[0][(indx1 - h1)]))) / c65535v);
                sgv = (epsv + (vabsf(LVFU(rgb[1][(indx + v1) >> 1]) - LVFU(rgb[1][(indx + v3) >> 1])) + vabsf(LVFU(rgb[0][indx1]) - LVFU(rgb[0][(indx1 + v1)]))) / c65535v);
                //N,E,W,S High Order Interpolation (Li & Randhawa)
                //N,E,W,S Hamilton Adams Interpolation
                // (48.f * 65535.f) = 3145680.f
                tempv = c40v * LVFU(rgb[0][indx1]);
                nvv = LIMV(((c23v * LVFU(rgb[1][(indx - v1) >> 1]) + c23v * LVFU(rgb[1][(indx - v3) >> 1]) + LVFU(rgb[1][(indx - v5) >> 1]) + LVFU(rgb[1][(indx + v1) >> 1]) + tempv - c32v * LVFU(rgb[0][(indx1 - v1)]) - c8v * LVFU(rgb[0][(indx1 - v2)]))) / c3145680v, zerov, onev);
                evv = LIMV(((c23v * LVFU(rgb[1][(indx + h1) >> 1]) + c23v * LVFU(rgb[1][(indx + h3) >> 1]) + LVFU(rgb[1][(indx + h5) >> 1]) + LVFU(rgb[1][(indx - h1) >> 1]) + tempv - c32v * LVFU(rgb[0][(indx1 + h1)]) - c8v * LVFU(rgb[0][(indx1 + h2)]))) / c3145680v, zerov, onev);
                wvv = LIMV(((c23v * LVFU(rgb[1][(indx - h1) >> 1]) + c23v * LVFU(rgb[1][(indx - h3) >> 1]) + LVFU(rgb[1][(indx - h5) >> 1]) + LVFU(rgb[1][(indx + h1) >> 1]) + tempv - c32v * LVFU(rgb[0][(indx1 - h1)]) - c8v * LVFU(rgb[0][(indx1 - h2)]))) / c3145680v, zerov, onev);
                svv = LIMV(((c23v * LVFU(rgb[1][(indx + v1) >> 1]) + c23v * LVFU(rgb[1][(indx + v3) >> 1]) + LVFU(rgb[1][(indx + v5) >> 1]) + LVFU(rgb[1][(indx - v1) >> 1]) + tempv - c32v * LVFU(rgb[0][(indx1 + v1)]) - c8v * LVFU(rgb[0][(indx1 + v2)]))) / c3145680v, zerov, onev);
                //Horizontal and vertical color differences
                tempv = LVFU( rgb[0][indx1] ) / c65535v;
                _mm_storeu_ps( &vdif[indx1], (sgv * nvv + ngv * svv) / (ngv + sgv) - tempv );
                _mm_storeu_ps( &hdif[indx1], (wgv * evv + egv * wvv) / (egv + wgv) - tempv );
            }

            // borders without SSE
            for (; col < width - 5; col += 2, indx += 2, indx1++) {
                //N,E,W,S Gradients
                ng = (eps + (fabsf(rgb[1][(indx - v1) >> 1] - rgb[1][(indx - v3) >> 1]) + fabsf(rgb[0][indx1] - rgb[0][(indx1 - v1)])) / 65535.f);;
                eg = (eps + (fabsf(rgb[1][(indx + h1) >> 1] - rgb[1][(indx + h3) >> 1]) + fabsf(rgb[0][indx1] - rgb[0][(indx1 + h1)])) / 65535.f);
                wg = (eps + (fabsf(rgb[1][(indx - h1) >> 1] - rgb[1][(indx - h3) >> 1]) + fabsf(rgb[0][indx1] - rgb[0][(indx1 - h1)])) / 65535.f);
                sg = (eps + (fabsf(rgb[1][(indx + v1) >> 1] - rgb[1][(indx + v3) >> 1]) + fabsf(rgb[0][indx1] - rgb[0][(indx1 + v1)])) / 65535.f);
                //N,E,W,S High Order Interpolation (Li & Randhawa)
                //N,E,W,S Hamilton Adams Interpolation
                // (48.f * 65535.f) = 3145680.f
                nv = LIM(((23.0f * rgb[1][(indx - v1) >> 1] + 23.0f * rgb[1][(indx - v3) >> 1] + rgb[1][(indx - v5) >> 1] + rgb[1][(indx + v1) >> 1] + 40.0f * rgb[0][indx1] - 32.0f * rgb[0][(indx1 - v1)] - 8.0f * rgb[0][(indx1 - v2)])) / 3145680.f, 0.0f, 1.0f);
                ev = LIM(((23.0f * rgb[1][(indx + h1) >> 1] + 23.0f * rgb[1][(indx + h3) >> 1] + rgb[1][(indx + h5) >> 1] + rgb[1][(indx - h1) >> 1] + 40.0f * rgb[0][indx1] - 32.0f * rgb[0][(indx1 + h1)] - 8.0f * rgb[0][(indx1 + h2)])) / 3145680.f, 0.0f, 1.0f);
                wv = LIM(((23.0f * rgb[1][(indx - h1) >> 1] + 23.0f * rgb[1][(indx - h3) >> 1] + rgb[1][(indx - h5) >> 1] + rgb[1][(indx + h1) >> 1] + 40.0f * rgb[0][indx1] - 32.0f * rgb[0][(indx1 - h1)] - 8.0f * rgb[0][(indx1 - h2)])) / 3145680.f, 0.0f, 1.0f);
                sv = LIM(((23.0f * rgb[1][(indx + v1) >> 1] + 23.0f * rgb[1][(indx + v3) >> 1] + rgb[1][(indx + v5) >> 1] + rgb[1][(indx - v1) >> 1] + 40.0f * rgb[0][indx1] - 32.0f * rgb[0][(indx1 + v1)] - 8.0f * rgb[0][(indx1 + v2)])) / 3145680.f, 0.0f, 1.0f);
                //Horizontal and vertical color differences
                vdif[indx1] = (sg * nv + ng * sv) / (ng + sg) - (rgb[0][indx1]) / 65535.f;
                hdif[indx1] = (wg * ev + eg * wv) / (eg + wg) - (rgb[0][indx1]) / 65535.f;
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.26);
            }
        }
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row++) {
            int col, d, indx1;

            for (col = 7 + (FC(row, 1) & 1), indx1 = (row * width + col) >> 1, d = FC(row, col) / 2; col < width - 14; col += 8, indx1 += 4) {
                //H&V integrated gaussian vector over variance on color differences
                //Mod Jacques 3/2013
                ngv = LIMV(epssqv + c78v * SQRV(LVFU(vdif[indx1])) + c69v * (SQRV(LVFU(vdif[indx1 - v1])) + SQRV(LVFU(vdif[indx1 + v1]))) + c51v * (SQRV(LVFU(vdif[indx1 - v2])) + SQRV(LVFU(vdif[indx1 + v2]))) + c21v * (SQRV(LVFU(vdif[indx1 - v3])) + SQRV(LVFU(vdif[indx1 + v3]))) - c6v * SQRV(LVFU(vdif[indx1 - v1]) + LVFU(vdif[indx1]) + LVFU(vdif[indx1 + v1]))
                           - c10v * (SQRV(LVFU(vdif[indx1 - v2]) + LVFU(vdif[indx1 - v1]) + LVFU(vdif[indx1])) + SQRV(LVFU(vdif[indx1]) + LVFU(vdif[indx1 + v1]) + LVFU(vdif[indx1 + v2]))) - c7v * (SQRV(LVFU(vdif[indx1 - v3]) + LVFU(vdif[indx1 - v2]) + LVFU(vdif[indx1 - v1])) + SQRV(LVFU(vdif[indx1 + v1]) + LVFU(vdif[indx1 + v2]) + LVFU(vdif[indx1 + v3]))), zerov, onev);
                egv = LIMV(epssqv + c78v * SQRV(LVFU(hdif[indx1])) + c69v * (SQRV(LVFU(hdif[indx1 - h1])) + SQRV(LVFU(hdif[indx1 + h1]))) + c51v * (SQRV(LVFU(hdif[indx1 - h2])) + SQRV(LVFU(hdif[indx1 + h2]))) + c21v * (SQRV(LVFU(hdif[indx1 - h3])) + SQRV(LVFU(hdif[indx1 + h3]))) - c6v * SQRV(LVFU(hdif[indx1 - h1]) + LVFU(hdif[indx1]) + LVFU(hdif[indx1 + h1]))
                           - c10v * (SQRV(LVFU(hdif[indx1 - h2]) + LVFU(hdif[indx1 - h1]) + LVFU(hdif[indx1])) + SQRV(LVFU(hdif[indx1]) + LVFU(hdif[indx1 + h1]) + LVFU(hdif[indx1 + h2]))) - c7v * (SQRV(LVFU(hdif[indx1 - h3]) + LVFU(hdif[indx1 - h2]) + LVFU(hdif[indx1 - h1])) + SQRV(LVFU(hdif[indx1 + h1]) + LVFU(hdif[indx1 + h2]) + LVFU(hdif[indx1 + h3]))), zerov, onev);
                //Limit chrominance using H/V neighbourhood
                nvv = median(d725v * LVFU(vdif[indx1]) + d1375v * LVFU(vdif[indx1 - v1]) + d1375v * LVFU(vdif[indx1 + v1]), LVFU(vdif[indx1 - v1]), LVFU(vdif[indx1 + v1]));
                evv = median(d725v * LVFU(hdif[indx1]) + d1375v * LVFU(hdif[indx1 - h1]) + d1375v * LVFU(hdif[indx1 + h1]), LVFU(hdif[indx1 - h1]), LVFU(hdif[indx1 + h1]));
                //Chrominance estimation
                tempv = (egv * nvv + ngv * evv) / (ngv + egv);
                _mm_storeu_ps(&(chr[d][indx1]), tempv);
                //Green channel population
                temp1v = c65535v * tempv + LVFU(rgb[0][indx1]);
                _mm_storeu_ps( &(rgb[0][indx1]), temp1v );
            }

            for (; col < width - 7; col += 2, indx1++) {
                //H&V integrated gaussian vector over variance on color differences
                //Mod Jacques 3/2013
                ng = LIM(epssq + 78.0f * SQR(vdif[indx1]) + 69.0f * (SQR(vdif[indx1 - v1]) + SQR(vdif[indx1 + v1])) + 51.0f * (SQR(vdif[indx1 - v2]) + SQR(vdif[indx1 + v2])) + 21.0f * (SQR(vdif[indx1 - v3]) + SQR(vdif[indx1 + v3])) - 6.0f * SQR(vdif[indx1 - v1] + vdif[indx1] + vdif[indx1 + v1])
                         - 10.0f * (SQR(vdif[indx1 - v2] + vdif[indx1 - v1] + vdif[indx1]) + SQR(vdif[indx1] + vdif[indx1 + v1] + vdif[indx1 + v2])) - 7.0f * (SQR(vdif[indx1 - v3] + vdif[indx1 - v2] + vdif[indx1 - v1]) + SQR(vdif[indx1 + v1] + vdif[indx1 + v2] + vdif[indx1 + v3])), 0.f, 1.f);
                eg = LIM(epssq + 78.0f * SQR(hdif[indx1]) + 69.0f * (SQR(hdif[indx1 - h1]) + SQR(hdif[indx1 + h1])) + 51.0f * (SQR(hdif[indx1 - h2]) + SQR(hdif[indx1 + h2])) + 21.0f * (SQR(hdif[indx1 - h3]) + SQR(hdif[indx1 + h3])) - 6.0f * SQR(hdif[indx1 - h1] + hdif[indx1] + hdif[indx1 + h1])
                         - 10.0f * (SQR(hdif[indx1 - h2] + hdif[indx1 - h1] + hdif[indx1]) + SQR(hdif[indx1] + hdif[indx1 + h1] + hdif[indx1 + h2])) - 7.0f * (SQR(hdif[indx1 - h3] + hdif[indx1 - h2] + hdif[indx1 - h1]) + SQR(hdif[indx1 + h1] + hdif[indx1 + h2] + hdif[indx1 + h3])), 0.f, 1.f);
                //Limit chrominance using H/V neighbourhood
                nv = median(0.725f * vdif[indx1] + 0.1375f * vdif[indx1 - v1] + 0.1375f * vdif[indx1 + v1], vdif[indx1 - v1], vdif[indx1 + v1]);
                ev = median(0.725f * hdif[indx1] + 0.1375f * hdif[indx1 - h1] + 0.1375f * hdif[indx1 + h1], hdif[indx1 - h1], hdif[indx1 + h1]);
                //Chrominance estimation
                chr[d][indx1] = (eg * nv + ng * ev) / (ng + eg);
                //Green channel population
                rgb[0][indx1] = rgb[0][indx1] + 65535.f * chr[d][indx1];
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.39);
            }
        }
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row++) {
            int col, indx, c;

            for (col = 7 + (FC(row, 1) & 1), indx = row * width + col, c = 1 - FC(row, col) / 2; col < width - 14; col += 8, indx += 8) {
                //NW,NE,SW,SE Gradients
                nwgv = onev / (epsv + vabsf(LVFU(chr[c][(indx - v1 - h1) >> 1]) - LVFU(chr[c][(indx - v3 - h3) >> 1])) + vabsf(LVFU(chr[c][(indx + v1 + h1) >> 1]) - LVFU(chr[c][(indx - v3 - h3) >> 1])));
                negv = onev / (epsv + vabsf(LVFU(chr[c][(indx - v1 + h1) >> 1]) - LVFU(chr[c][(indx - v3 + h3) >> 1])) + vabsf(LVFU(chr[c][(indx + v1 - h1) >> 1]) - LVFU(chr[c][(indx - v3 + h3) >> 1])));
                swgv = onev / (epsv + vabsf(LVFU(chr[c][(indx + v1 - h1) >> 1]) - LVFU(chr[c][(indx + v3 + h3) >> 1])) + vabsf(LVFU(chr[c][(indx - v1 + h1) >> 1]) - LVFU(chr[c][(indx + v3 - h3) >> 1])));
                segv = onev / (epsv + vabsf(LVFU(chr[c][(indx + v1 + h1) >> 1]) - LVFU(chr[c][(indx + v3 - h3) >> 1])) + vabsf(LVFU(chr[c][(indx - v1 - h1) >> 1]) - LVFU(chr[c][(indx + v3 + h3) >> 1])));
                //Limit NW,NE,SW,SE Color differences
                nwvv = median(LVFU(chr[c][(indx - v1 - h1) >> 1]), LVFU(chr[c][(indx - v3 - h1) >> 1]), LVFU(chr[c][(indx - v1 - h3) >> 1]));
                nevv = median(LVFU(chr[c][(indx - v1 + h1) >> 1]), LVFU(chr[c][(indx - v3 + h1) >> 1]), LVFU(chr[c][(indx - v1 + h3) >> 1]));
                swvv = median(LVFU(chr[c][(indx + v1 - h1) >> 1]), LVFU(chr[c][(indx + v3 - h1) >> 1]), LVFU(chr[c][(indx + v1 - h3) >> 1]));
                sevv = median(LVFU(chr[c][(indx + v1 + h1) >> 1]), LVFU(chr[c][(indx + v3 + h1) >> 1]), LVFU(chr[c][(indx + v1 + h3) >> 1]));
                //Interpolate chrominance: R@B and B@R
                tempv = (nwgv * nwvv + negv * nevv + swgv * swvv + segv * sevv) / (nwgv + negv + swgv + segv);
                _mm_storeu_ps( &(chr[c][indx >> 1]), tempv);
            }

            for (; col < width - 7; col += 2, indx += 2) {
                //NW,NE,SW,SE Gradients
                nwg = 1.0f / (eps + fabsf(chr[c][(indx - v1 - h1) >> 1] - chr[c][(indx - v3 - h3) >> 1]) + fabsf(chr[c][(indx + v1 + h1) >> 1] - chr[c][(indx - v3 - h3) >> 1]));
                neg = 1.0f / (eps + fabsf(chr[c][(indx - v1 + h1) >> 1] - chr[c][(indx - v3 + h3) >> 1]) + fabsf(chr[c][(indx + v1 - h1) >> 1] - chr[c][(indx - v3 + h3) >> 1]));
                swg = 1.0f / (eps + fabsf(chr[c][(indx + v1 - h1) >> 1] - chr[c][(indx + v3 + h3) >> 1]) + fabsf(chr[c][(indx - v1 + h1) >> 1] - chr[c][(indx + v3 - h3) >> 1]));
                seg = 1.0f / (eps + fabsf(chr[c][(indx + v1 + h1) >> 1] - chr[c][(indx + v3 - h3) >> 1]) + fabsf(chr[c][(indx - v1 - h1) >> 1] - chr[c][(indx + v3 + h3) >> 1]));
                //Limit NW,NE,SW,SE Color differences
                nwv = median(chr[c][(indx - v1 - h1) >> 1], chr[c][(indx - v3 - h1) >> 1], chr[c][(indx - v1 - h3) >> 1]);
                nev = median(chr[c][(indx - v1 + h1) >> 1], chr[c][(indx - v3 + h1) >> 1], chr[c][(indx - v1 + h3) >> 1]);
                swv = median(chr[c][(indx + v1 - h1) >> 1], chr[c][(indx + v3 - h1) >> 1], chr[c][(indx + v1 - h3) >> 1]);
                sev = median(chr[c][(indx + v1 + h1) >> 1], chr[c][(indx + v3 + h1) >> 1], chr[c][(indx + v1 + h3) >> 1]);
                //Interpolate chrominance: R@B and B@R
                chr[c][indx >> 1] = (nwg * nwv + neg * nev + swg * swv + seg * sev) / (nwg + neg + swg + seg);
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.65);
            }
        }
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row++) {
            int col, indx;

            for (col = 7 + (FC(row, 0) & 1), indx = row * width + col; col < width - 14; col += 8, indx += 8) {
                //N,E,W,S Gradients
                ngv = onev / (epsv + vabsf(LVFU(chr[0][(indx - v1) >> 1]) - LVFU(chr[0][(indx - v3) >> 1])) + vabsf(LVFU(chr[0][(indx + v1) >> 1]) - LVFU(chr[0][(indx - v3) >> 1])));
                egv = onev / (epsv + vabsf(LVFU(chr[0][(indx + h1) >> 1]) - LVFU(chr[0][(indx + h3) >> 1])) + vabsf(LVFU(chr[0][(indx - h1) >> 1]) - LVFU(chr[0][(indx + h3) >> 1])));
                wgv = onev / (epsv + vabsf(LVFU(chr[0][(indx - h1) >> 1]) - LVFU(chr[0][(indx - h3) >> 1])) + vabsf(LVFU(chr[0][(indx + h1) >> 1]) - LVFU(chr[0][(indx - h3) >> 1])));
                sgv = onev / (epsv + vabsf(LVFU(chr[0][(indx + v1) >> 1]) - LVFU(chr[0][(indx + v3) >> 1])) + vabsf(LVFU(chr[0][(indx - v1) >> 1]) - LVFU(chr[0][(indx + v3) >> 1])));
                //Interpolate chrominance: R@G and B@G
                tempv = ((ngv * LVFU(chr[0][(indx - v1) >> 1]) + egv * LVFU(chr[0][(indx + h1) >> 1]) + wgv * LVFU(chr[0][(indx - h1) >> 1]) + sgv * LVFU(chr[0][(indx + v1) >> 1])) / (ngv + egv + wgv + sgv));
                _mm_storeu_ps( &chr[0 + 2][indx >> 1], tempv);
            }

            for (; col < width - 7; col += 2, indx += 2) {
                //N,E,W,S Gradients
                ng = 1.0f / (eps + fabsf(chr[0][(indx - v1) >> 1] - chr[0][(indx - v3) >> 1]) + fabsf(chr[0][(indx + v1) >> 1] - chr[0][(indx - v3) >> 1]));
                eg = 1.0f / (eps + fabsf(chr[0][(indx + h1) >> 1] - chr[0][(indx + h3) >> 1]) + fabsf(chr[0][(indx - h1) >> 1] - chr[0][(indx + h3) >> 1]));
                wg = 1.0f / (eps + fabsf(chr[0][(indx - h1) >> 1] - chr[0][(indx - h3) >> 1]) + fabsf(chr[0][(indx + h1) >> 1] - chr[0][(indx - h3) >> 1]));
                sg = 1.0f / (eps + fabsf(chr[0][(indx + v1) >> 1] - chr[0][(indx + v3) >> 1]) + fabsf(chr[0][(indx - v1) >> 1] - chr[0][(indx + v3) >> 1]));
                //Interpolate chrominance: R@G and B@G
                chr[0 + 2][indx >> 1] = ((ng * chr[0][(indx - v1) >> 1] + eg * chr[0][(indx + h1) >> 1] + wg * chr[0][(indx - h1) >> 1] + sg * chr[0][(indx + v1) >> 1]) / (ng + eg + wg + sg));
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.78);
            }
        }
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row++) {
            int col, indx;

            for (col = 7 + (FC(row, 0) & 1), indx = row * width + col; col < width - 14; col += 8, indx += 8) {
                //N,E,W,S Gradients
                ngv = onev / (epsv + vabsf(LVFU(chr[1][(indx - v1) >> 1]) - LVFU(chr[1][(indx - v3) >> 1])) + vabsf(LVFU(chr[1][(indx + v1) >> 1]) - LVFU(chr[1][(indx - v3) >> 1])));
                egv = onev / (epsv + vabsf(LVFU(chr[1][(indx + h1) >> 1]) - LVFU(chr[1][(indx + h3) >> 1])) + vabsf(LVFU(chr[1][(indx - h1) >> 1]) - LVFU(chr[1][(indx + h3) >> 1])));
                wgv = onev / (epsv + vabsf(LVFU(chr[1][(indx - h1) >> 1]) - LVFU(chr[1][(indx - h3) >> 1])) + vabsf(LVFU(chr[1][(indx + h1) >> 1]) - LVFU(chr[1][(indx - h3) >> 1])));
                sgv = onev / (epsv + vabsf(LVFU(chr[1][(indx + v1) >> 1]) - LVFU(chr[1][(indx + v3) >> 1])) + vabsf(LVFU(chr[1][(indx - v1) >> 1]) - LVFU(chr[1][(indx + v3) >> 1])));
                //Interpolate chrominance: R@G and B@G
                tempv = ((ngv * LVFU(chr[1][(indx - v1) >> 1]) + egv * LVFU(chr[1][(indx + h1) >> 1]) + wgv * LVFU(chr[1][(indx - h1) >> 1]) + sgv * LVFU(chr[1][(indx + v1) >> 1])) / (ngv + egv + wgv + sgv));
                _mm_storeu_ps( &chr[1 + 2][indx >> 1], tempv);
            }

            for (; col < width - 7; col += 2, indx += 2) {
                //N,E,W,S Gradients
                ng = 1.0f / (eps + fabsf(chr[1][(indx - v1) >> 1] - chr[1][(indx - v3) >> 1]) + fabsf(chr[1][(indx + v1) >> 1] - chr[1][(indx - v3) >> 1]));
                eg = 1.0f / (eps + fabsf(chr[1][(indx + h1) >> 1] - chr[1][(indx + h3) >> 1]) + fabsf(chr[1][(indx - h1) >> 1] - chr[1][(indx + h3) >> 1]));
                wg = 1.0f / (eps + fabsf(chr[1][(indx - h1) >> 1] - chr[1][(indx - h3) >> 1]) + fabsf(chr[1][(indx + h1) >> 1] - chr[1][(indx - h3) >> 1]));
                sg = 1.0f / (eps + fabsf(chr[1][(indx + v1) >> 1] - chr[1][(indx + v3) >> 1]) + fabsf(chr[1][(indx - v1) >> 1] - chr[1][(indx + v3) >> 1]));
                //Interpolate chrominance: R@G and B@G
                chr[1 + 2][indx >> 1] = ((ng * chr[1][(indx - v1) >> 1] + eg * chr[1][(indx + h1) >> 1] + wg * chr[1][(indx - h1) >> 1] + sg * chr[1][(indx + v1) >> 1]) / (ng + eg + wg + sg));
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.91);
            }
        }
        float *src1, *src2, *redsrc0, *redsrc1, *bluesrc0, *bluesrc1;
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int row = 7; row < height - 7; row++) {
            int col, indx, fc;
            fc = FC(row, 7) & 1;
            src1 = rgb[fc];
            src2 = rgb[fc ^ 1];
            redsrc0 = chr[fc << 1];
            redsrc1 = chr[(fc ^ 1) << 1];
            bluesrc0 = chr[(fc << 1) + 1];
            bluesrc1 = chr[((fc ^ 1) << 1) + 1];

            for(col = 7, indx = row * width + col; col < width - 14; col += 8, indx += 8) {
                temp1v = LVFU( src1[indx >> 1] );
                temp2v = LVFU( src2[(indx + 1) >> 1] );
                tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 1, 0, 1, 0 ) );
                tempv = _mm_shuffle_ps( tempv, tempv, _MM_SHUFFLE( 3, 1, 2, 0 ) );
                _mm_storeu_ps( &green[row][col], CLIPV( tempv ));
                temp5v = LVFU(redsrc0[indx >> 1]);
                temp6v = LVFU(redsrc1[(indx + 1) >> 1]);
                temp3v = _mm_shuffle_ps( temp5v, temp6v, _MM_SHUFFLE( 1, 0, 1, 0 ) );
                temp3v = _mm_shuffle_ps( temp3v, temp3v, _MM_SHUFFLE( 3, 1, 2, 0 ) );
                temp3v = CLIPV( tempv - c65535v * temp3v );
                _mm_storeu_ps( &red[row][col], temp3v);
                temp7v = LVFU(bluesrc0[indx >> 1]);
                temp8v = LVFU(bluesrc1[(indx + 1) >> 1]);
                temp4v = _mm_shuffle_ps( temp7v, temp8v, _MM_SHUFFLE( 1, 0, 1, 0 ) );
                temp4v = _mm_shuffle_ps( temp4v, temp4v, _MM_SHUFFLE( 3, 1, 2, 0 ) );
                temp4v = CLIPV( tempv - c65535v * temp4v );
                _mm_storeu_ps( &blue[row][col], temp4v);

                tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 3, 2, 3, 2 ) );
                tempv = _mm_shuffle_ps( tempv, tempv, _MM_SHUFFLE( 3, 1, 2, 0 ) );
                _mm_storeu_ps( &green[row][col + 4], CLIPV( tempv ));

                temp3v = _mm_shuffle_ps( temp5v, temp6v, _MM_SHUFFLE( 3, 2, 3, 2 ) );
                temp3v = _mm_shuffle_ps( temp3v, temp3v, _MM_SHUFFLE( 3, 1, 2, 0 ) );
                temp3v = CLIPV( tempv - c65535v * temp3v );
                _mm_storeu_ps( &red[row][col + 4], temp3v);
                temp4v = _mm_shuffle_ps( temp7v, temp8v, _MM_SHUFFLE( 3, 2, 3, 2 ) );
                temp4v = _mm_shuffle_ps( temp4v, temp4v, _MM_SHUFFLE( 3, 1, 2, 0 ) );
                temp4v = CLIPV( tempv - c65535v * temp4v );
                _mm_storeu_ps( &blue[row][col + 4], temp4v);
            }

            for(; col < width - 7; col++, indx += 2) {
                red  [row][col] = CLIP(src1[indx >> 1] - 65535.f * redsrc0[indx >> 1]);
                green[row][col] = CLIP(src1[indx >> 1]);
                blue [row][col] = CLIP(src1[indx >> 1] - 65535.f * bluesrc0[indx >> 1]);
                col++;
                red  [row][col] = CLIP(src2[(indx + 1) >> 1] - 65535.f * redsrc1[(indx + 1) >> 1]);
                green[row][col] = CLIP(src2[(indx + 1) >> 1]);
                blue [row][col] = CLIP(src2[(indx + 1) >> 1] - 65535.f * bluesrc1[(indx + 1) >> 1]);
            }
        }
    }// End of parallelization

    if (plistener) {
        plistener->setProgress (1.0);
    }

    free(chrarray);
    free(rgbarray);
    free(vdif);
    free(hdif);
}
#undef CLIPV
#else
void RawImageSource::igv_interpolate(int winw, int winh)
{
    static const float eps = 1e-5f, epssq = 1e-5f; //mod epssq -10f =>-5f Jacques 3/2013 to prevent artifact (divide by zero)
    static const int h1 = 1, h2 = 2, h3 = 3, h4 = 4, h5 = 5, h6 = 6;
    const int width = winw, height = winh;
    const int v1 = 1 * width, v2 = 2 * width, v3 = 3 * width, v4 = 4 * width, v5 = 5 * width, v6 = 6 * width;
    float* rgb[3];
    float* chr[2];
    float (*rgbarray), *vdif, *hdif, (*chrarray);

    rgbarray    = (float (*)) calloc(width * height * 3, sizeof( float));
    rgb[0] = rgbarray;
    rgb[1] = rgbarray + (width * height);
    rgb[2] = rgbarray + 2 * (width * height);

    chrarray    = (float (*)) calloc(width * height * 2, sizeof( float));
    chr[0] = chrarray;
    chr[1] = chrarray + (width * height);

    vdif  = (float (*))    calloc(width * height / 2, sizeof * vdif);
    hdif  = (float (*))    calloc(width * height / 2, sizeof * hdif);

    border_interpolate2(winw, winh, 7, rawData, red, green, blue);

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::IGV)));
        plistener->setProgress (0.0);
    }

#ifdef _OPENMP
    #pragma omp parallel default(none) shared(rgb,vdif,hdif,chr)
#endif
    {

        float ng, eg, wg, sg, nv, ev, wv, sv, nwg, neg, swg, seg, nwv, nev, swv, sev;

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 0; row < height - 0; row++)
            for (int col = 0, indx = row * width + col; col < width - 0; col++, indx++) {
                int c = FC(row, col);
                rgb[c][indx] = CLIP(rawData[row][col]); //rawData = RT datas
            }

//  border_interpolate2(7, rgb);

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.13);
            }
        }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 5; row < height - 5; row++)
            for (int col = 5 + (FC(row, 1) & 1), indx = row * width + col, c = FC(row, col); col < width - 5; col += 2, indx += 2) {
                //N,E,W,S Gradients
                ng = (eps + (fabsf(rgb[1][indx - v1] - rgb[1][indx - v3]) + fabsf(rgb[c][indx] - rgb[c][indx - v2])) / 65535.f);;
                eg = (eps + (fabsf(rgb[1][indx + h1] - rgb[1][indx + h3]) + fabsf(rgb[c][indx] - rgb[c][indx + h2])) / 65535.f);
                wg = (eps + (fabsf(rgb[1][indx - h1] - rgb[1][indx - h3]) + fabsf(rgb[c][indx] - rgb[c][indx - h2])) / 65535.f);
                sg = (eps + (fabsf(rgb[1][indx + v1] - rgb[1][indx + v3]) + fabsf(rgb[c][indx] - rgb[c][indx + v2])) / 65535.f);
                //N,E,W,S High Order Interpolation (Li & Randhawa)
                //N,E,W,S Hamilton Adams Interpolation
                // (48.f * 65535.f) = 3145680.f
                nv = LIM(((23.0f * rgb[1][indx - v1] + 23.0f * rgb[1][indx - v3] + rgb[1][indx - v5] + rgb[1][indx + v1] + 40.0f * rgb[c][indx] - 32.0f * rgb[c][indx - v2] - 8.0f * rgb[c][indx - v4])) / 3145680.f, 0.0f, 1.0f);
                ev = LIM(((23.0f * rgb[1][indx + h1] + 23.0f * rgb[1][indx + h3] + rgb[1][indx + h5] + rgb[1][indx - h1] + 40.0f * rgb[c][indx] - 32.0f * rgb[c][indx + h2] - 8.0f * rgb[c][indx + h4])) / 3145680.f, 0.0f, 1.0f);
                wv = LIM(((23.0f * rgb[1][indx - h1] + 23.0f * rgb[1][indx - h3] + rgb[1][indx - h5] + rgb[1][indx + h1] + 40.0f * rgb[c][indx] - 32.0f * rgb[c][indx - h2] - 8.0f * rgb[c][indx - h4])) / 3145680.f, 0.0f, 1.0f);
                sv = LIM(((23.0f * rgb[1][indx + v1] + 23.0f * rgb[1][indx + v3] + rgb[1][indx + v5] + rgb[1][indx - v1] + 40.0f * rgb[c][indx] - 32.0f * rgb[c][indx + v2] - 8.0f * rgb[c][indx + v4])) / 3145680.f, 0.0f, 1.0f);
                //Horizontal and vertical color differences
                vdif[indx >> 1] = (sg * nv + ng * sv) / (ng + sg) - (rgb[c][indx]) / 65535.f;
                hdif[indx >> 1] = (wg * ev + eg * wv) / (eg + wg) - (rgb[c][indx]) / 65535.f;
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.26);
            }
        }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row++)
            for (int col = 7 + (FC(row, 1) & 1), indx = row * width + col, c = FC(row, col), d = c / 2; col < width - 7; col += 2, indx += 2) {
                //H&V integrated gaussian vector over variance on color differences
                //Mod Jacques 3/2013
                ng = LIM(epssq + 78.0f * SQR(vdif[indx >> 1]) + 69.0f * (SQR(vdif[(indx - v2) >> 1]) + SQR(vdif[(indx + v2) >> 1])) + 51.0f * (SQR(vdif[(indx - v4) >> 1]) + SQR(vdif[(indx + v4) >> 1])) + 21.0f * (SQR(vdif[(indx - v6) >> 1]) + SQR(vdif[(indx + v6) >> 1])) - 6.0f * SQR(vdif[(indx - v2) >> 1] + vdif[indx >> 1] + vdif[(indx + v2) >> 1])
                         - 10.0f * (SQR(vdif[(indx - v4) >> 1] + vdif[(indx - v2) >> 1] + vdif[indx >> 1]) + SQR(vdif[indx >> 1] + vdif[(indx + v2) >> 1] + vdif[(indx + v4) >> 1])) - 7.0f * (SQR(vdif[(indx - v6) >> 1] + vdif[(indx - v4) >> 1] + vdif[(indx - v2) >> 1]) + SQR(vdif[(indx + v2) >> 1] + vdif[(indx + v4) >> 1] + vdif[(indx + v6) >> 1])), 0.f, 1.f);
                eg = LIM(epssq + 78.0f * SQR(hdif[indx >> 1]) + 69.0f * (SQR(hdif[(indx - h2) >> 1]) + SQR(hdif[(indx + h2) >> 1])) + 51.0f * (SQR(hdif[(indx - h4) >> 1]) + SQR(hdif[(indx + h4) >> 1])) + 21.0f * (SQR(hdif[(indx - h6) >> 1]) + SQR(hdif[(indx + h6) >> 1])) - 6.0f * SQR(hdif[(indx - h2) >> 1] + hdif[indx >> 1] + hdif[(indx + h2) >> 1])
                         - 10.0f * (SQR(hdif[(indx - h4) >> 1] + hdif[(indx - h2) >> 1] + hdif[indx >> 1]) + SQR(hdif[indx >> 1] + hdif[(indx + h2) >> 1] + hdif[(indx + h4) >> 1])) - 7.0f * (SQR(hdif[(indx - h6) >> 1] + hdif[(indx - h4) >> 1] + hdif[(indx - h2) >> 1]) + SQR(hdif[(indx + h2) >> 1] + hdif[(indx + h4) >> 1] + hdif[(indx + h6) >> 1])), 0.f, 1.f);
                //Limit chrominance using H/V neighbourhood
                nv = median(0.725f * vdif[indx >> 1] + 0.1375f * vdif[(indx - v2) >> 1] + 0.1375f * vdif[(indx + v2) >> 1], vdif[(indx - v2) >> 1], vdif[(indx + v2) >> 1]);
                ev = median(0.725f * hdif[indx >> 1] + 0.1375f * hdif[(indx - h2) >> 1] + 0.1375f * hdif[(indx + h2) >> 1], hdif[(indx - h2) >> 1], hdif[(indx + h2) >> 1]);
                //Chrominance estimation
                chr[d][indx] = (eg * nv + ng * ev) / (ng + eg);
                //Green channel population
                rgb[1][indx] = rgb[c][indx] + 65535.f * chr[d][indx];
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.39);
            }
        }

//  free(vdif); free(hdif);
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row += 2)
            for (int col = 7 + (FC(row, 1) & 1), indx = row * width + col, c = 1 - FC(row, col) / 2; col < width - 7; col += 2, indx += 2) {
                //NW,NE,SW,SE Gradients
                nwg = 1.0f / (eps + fabsf(chr[c][indx - v1 - h1] - chr[c][indx - v3 - h3]) + fabsf(chr[c][indx + v1 + h1] - chr[c][indx - v3 - h3]));
                neg = 1.0f / (eps + fabsf(chr[c][indx - v1 + h1] - chr[c][indx - v3 + h3]) + fabsf(chr[c][indx + v1 - h1] - chr[c][indx - v3 + h3]));
                swg = 1.0f / (eps + fabsf(chr[c][indx + v1 - h1] - chr[c][indx + v3 + h3]) + fabsf(chr[c][indx - v1 + h1] - chr[c][indx + v3 - h3]));
                seg = 1.0f / (eps + fabsf(chr[c][indx + v1 + h1] - chr[c][indx + v3 - h3]) + fabsf(chr[c][indx - v1 - h1] - chr[c][indx + v3 + h3]));
                //Limit NW,NE,SW,SE Color differences
                nwv = median(chr[c][indx - v1 - h1], chr[c][indx - v3 - h1], chr[c][indx - v1 - h3]);
                nev = median(chr[c][indx - v1 + h1], chr[c][indx - v3 + h1], chr[c][indx - v1 + h3]);
                swv = median(chr[c][indx + v1 - h1], chr[c][indx + v3 - h1], chr[c][indx + v1 - h3]);
                sev = median(chr[c][indx + v1 + h1], chr[c][indx + v3 + h1], chr[c][indx + v1 + h3]);
                //Interpolate chrominance: R@B and B@R
                chr[c][indx] = (nwg * nwv + neg * nev + swg * swv + seg * sev) / (nwg + neg + swg + seg);
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.52);
            }
        }
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 8; row < height - 7; row += 2)
            for (int col = 7 + (FC(row, 1) & 1), indx = row * width + col, c = 1 - FC(row, col) / 2; col < width - 7; col += 2, indx += 2) {
                //NW,NE,SW,SE Gradients
                nwg = 1.0f / (eps + fabsf(chr[c][indx - v1 - h1] - chr[c][indx - v3 - h3]) + fabsf(chr[c][indx + v1 + h1] - chr[c][indx - v3 - h3]));
                neg = 1.0f / (eps + fabsf(chr[c][indx - v1 + h1] - chr[c][indx - v3 + h3]) + fabsf(chr[c][indx + v1 - h1] - chr[c][indx - v3 + h3]));
                swg = 1.0f / (eps + fabsf(chr[c][indx + v1 - h1] - chr[c][indx + v3 + h3]) + fabsf(chr[c][indx - v1 + h1] - chr[c][indx + v3 - h3]));
                seg = 1.0f / (eps + fabsf(chr[c][indx + v1 + h1] - chr[c][indx + v3 - h3]) + fabsf(chr[c][indx - v1 - h1] - chr[c][indx + v3 + h3]));
                //Limit NW,NE,SW,SE Color differences
                nwv = median(chr[c][indx - v1 - h1], chr[c][indx - v3 - h1], chr[c][indx - v1 - h3]);
                nev = median(chr[c][indx - v1 + h1], chr[c][indx - v3 + h1], chr[c][indx - v1 + h3]);
                swv = median(chr[c][indx + v1 - h1], chr[c][indx + v3 - h1], chr[c][indx + v1 - h3]);
                sev = median(chr[c][indx + v1 + h1], chr[c][indx + v3 + h1], chr[c][indx + v1 + h3]);
                //Interpolate chrominance: R@B and B@R
                chr[c][indx] = (nwg * nwv + neg * nev + swg * swv + seg * sev) / (nwg + neg + swg + seg);
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.65);
            }
        }
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row++)
            for (int col = 7 + (FC(row, 0) & 1), indx = row * width + col; col < width - 7; col += 2, indx += 2) {
                //N,E,W,S Gradients
                ng = 1.0f / (eps + fabsf(chr[0][indx - v1] - chr[0][indx - v3]) + fabsf(chr[0][indx + v1] - chr[0][indx - v3]));
                eg = 1.0f / (eps + fabsf(chr[0][indx + h1] - chr[0][indx + h3]) + fabsf(chr[0][indx - h1] - chr[0][indx + h3]));
                wg = 1.0f / (eps + fabsf(chr[0][indx - h1] - chr[0][indx - h3]) + fabsf(chr[0][indx + h1] - chr[0][indx - h3]));
                sg = 1.0f / (eps + fabsf(chr[0][indx + v1] - chr[0][indx + v3]) + fabsf(chr[0][indx - v1] - chr[0][indx + v3]));
                //Interpolate chrominance: R@G and B@G
                chr[0][indx] = ((ng * chr[0][indx - v1] + eg * chr[0][indx + h1] + wg * chr[0][indx - h1] + sg * chr[0][indx + v1]) / (ng + eg + wg + sg));
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.78);
            }
        }
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 7; row < height - 7; row++)
            for (int col = 7 + (FC(row, 0) & 1), indx = row * width + col; col < width - 7; col += 2, indx += 2) {

                //N,E,W,S Gradients
                ng = 1.0f / (eps + fabsf(chr[1][indx - v1] - chr[1][indx - v3]) + fabsf(chr[1][indx + v1] - chr[1][indx - v3]));
                eg = 1.0f / (eps + fabsf(chr[1][indx + h1] - chr[1][indx + h3]) + fabsf(chr[1][indx - h1] - chr[1][indx + h3]));
                wg = 1.0f / (eps + fabsf(chr[1][indx - h1] - chr[1][indx - h3]) + fabsf(chr[1][indx + h1] - chr[1][indx - h3]));
                sg = 1.0f / (eps + fabsf(chr[1][indx + v1] - chr[1][indx + v3]) + fabsf(chr[1][indx - v1] - chr[1][indx + v3]));
                //Interpolate chrominance: R@G and B@G
                chr[1][indx] = ((ng * chr[1][indx - v1] + eg * chr[1][indx + h1] + wg * chr[1][indx - h1] + sg * chr[1][indx + v1]) / (ng + eg + wg + sg));
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.91);
            }

            //Interpolate borders
//  border_interpolate2(7, rgb);
        }
        /*
#ifdef _OPENMP
        #pragma omp for
#endif
            for (int row=0; row < height; row++)  //borders
                for (int col=0; col < width; col++) {
                    if (col==7 && row >= 7 && row < height-7)
                        col = width-7;
                    int indxc=row*width+col;
                    red  [row][col] = rgb[indxc][0];
                    green[row][col] = rgb[indxc][1];
                    blue [row][col] = rgb[indxc][2];
                }
        */

#ifdef _OPENMP
        #pragma omp for
#endif

        for(int row = 7; row < height - 7; row++)
            for(int col = 7, indx = row * width + col; col < width - 7; col++, indx++) {
                red  [row][col] = CLIP(rgb[1][indx] - 65535.f * chr[0][indx]);
                green[row][col] = CLIP(rgb[1][indx]);
                blue [row][col] = CLIP(rgb[1][indx] - 65535.f * chr[1][indx]);
            }
    }// End of parallelization

    if (plistener) {
        plistener->setProgress (1.0);
    }

    free(chrarray);
    free(rgbarray);
    free(vdif);
    free(hdif);
}
#endif

void RawImageSource::nodemosaic(bool bw)
{
    red(W, H);
    green(W, H);
    blue(W, H);
    #pragma omp parallel for

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (bw) {
                red[i][j] = green[i][j] = blue[i][j] = rawData[i][j];
            } else if(ri->getSensorType() != ST_FUJI_XTRANS) {
                switch( FC(i, j)) {
                case 0:
                    red[i][j] = rawData[i][j];
                    green[i][j] = blue[i][j] = 0;
                    break;

                case 1:
                    green[i][j] = rawData[i][j];
                    red[i][j] = blue[i][j] = 0;
                    break;

                case 2:
                    blue[i][j] = rawData[i][j];
                    red[i][j] = green[i][j] = 0;
                    break;
                }
            } else {
                switch( ri->XTRANSFC(i, j)) {
                case 0:
                    red[i][j] = rawData[i][j];
                    green[i][j] = blue[i][j] = 0;
                    break;

                case 1:
                    green[i][j] = rawData[i][j];
                    red[i][j] = blue[i][j] = 0;
                    break;

                case 2:
                    blue[i][j] = rawData[i][j];
                    red[i][j] = green[i][j] = 0;
                    break;
                }
            }
        }
    }
}

/*
   Refinement based on EECI demosaicing algorithm by L. Chang and Y.P. Tan
   Paul Lee
   Adapted for RawTherapee - Jacques Desmis 04/2013
*/

#ifdef __SSE2__
#define CLIPV(a) LIMV(a,ZEROV,c65535v)
#endif
void RawImageSource::refinement(int PassCount)
{
    MyTime t1e, t2e;
    t1e.set();

    int width = W;
    int height = H;
    int w1 = width;
    int w2 = 2 * w1;

    if (plistener) {
        plistener->setProgressStr (M("TP_RAW_DMETHOD_PROGRESSBAR_REFINE"));
    }

    array2D<float> *rgb[3];
    rgb[0] = &red;
    rgb[1] = &green;
    rgb[2] = &blue;

    for (int b = 0; b < PassCount; b++) {
        if (plistener) {
            plistener->setProgress ((float)b / PassCount);
        }


#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float *pix[3];

            /* Reinforce interpolated green pixels on RED/BLUE pixel locations */
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 2; row < height - 2; row++) {
                int col = 2 + (FC(row, 2) & 1);
                int c = FC(row, col);
#ifdef __SSE2__
                __m128 dLv, dRv, dUv, dDv, v0v;
                __m128 onev = _mm_set1_ps(1.f);
                __m128 zd5v = _mm_set1_ps(0.5f);
                __m128 c65535v = _mm_set1_ps(65535.f);

                for (; col < width - 8; col += 8) {
                    int indx = row * width + col;
                    pix[c] = (float*)(*rgb[c]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    dLv = onev / (onev + vabsf(LC2VFU(pix[c][ -2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dRv = onev / (onev + vabsf(LC2VFU(pix[c][  2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dUv = onev / (onev + vabsf(LC2VFU(pix[c][-w2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    dDv = onev / (onev + vabsf(LC2VFU(pix[c][ w2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    v0v = CLIPV(LC2VFU(pix[c][0]) + zd5v + ((LC2VFU(pix[1][-1]) - LC2VFU(pix[c][-1])) * dLv + (LC2VFU(pix[1][1]) - LC2VFU(pix[c][1])) * dRv + (LC2VFU(pix[1][-w1]) - LC2VFU(pix[c][-w1])) * dUv + (LC2VFU(pix[1][w1]) - LC2VFU(pix[c][w1])) * dDv ) / (dLv + dRv + dUv + dDv));
                    STC2VFU(pix[1][0], v0v);
                }

#endif

                for (; col < width - 2; col += 2) {
                    int indx = row * width + col;
                    pix[c] = (float*)(*rgb[c]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    float dL = 1.f / (1.f + fabsf(pix[c][ -2] - pix[c][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dR = 1.f / (1.f + fabsf(pix[c][  2] - pix[c][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dU = 1.f / (1.f + fabsf(pix[c][-w2] - pix[c][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float dD = 1.f / (1.f + fabsf(pix[c][ w2] - pix[c][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float v0 = (pix[c][0] + 0.5f + ((pix[1][ -1] - pix[c][ -1]) * dL + (pix[1][  1] - pix[c][  1]) * dR + (pix[1][-w1] - pix[c][-w1]) * dU + (pix[1][ w1] - pix[c][ w1]) * dD ) / (dL + dR + dU + dD));
                    pix[1][0] = CLIP(v0);
                }
            }

            /* Reinforce interpolated red/blue pixels on GREEN pixel locations */
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 2; row < height - 2; row++) {
                int col = 2 + (FC(row, 3) & 1);
                int c = FC(row, col + 1);
#ifdef __SSE2__
                __m128 dLv, dRv, dUv, dDv, v0v;
                __m128 onev = _mm_set1_ps(1.f);
                __m128 zd5v = _mm_set1_ps(0.5f);
                __m128 c65535v = _mm_set1_ps(65535.f);

                for (; col < width - 8; col += 8) {
                    int indx = row * width + col;
                    pix[1] = (float*)(*rgb[1]) + indx;

                    for (int i = 0; i < 2; c = 2 - c, i++) {
                        pix[c] = (float*)(*rgb[c]) + indx;
                        dLv = onev / (onev + vabsf(LC2VFU(pix[1][ -2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][ 1]) - LC2VFU(pix[c][ -1])));
                        dRv = onev / (onev + vabsf(LC2VFU(pix[1][  2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][ 1]) - LC2VFU(pix[c][ -1])));
                        dUv = onev / (onev + vabsf(LC2VFU(pix[1][-w2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][w1]) - LC2VFU(pix[c][-w1])));
                        dDv = onev / (onev + vabsf(LC2VFU(pix[1][ w2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][w1]) - LC2VFU(pix[c][-w1])));
                        v0v = CLIPV(LC2VFU(pix[1][0]) + zd5v - ((LC2VFU(pix[1][-1]) - LC2VFU(pix[c][-1])) * dLv + (LC2VFU(pix[1][1]) - LC2VFU(pix[c][1])) * dRv + (LC2VFU(pix[1][-w1]) - LC2VFU(pix[c][-w1])) * dUv + (LC2VFU(pix[1][w1]) - LC2VFU(pix[c][w1])) * dDv ) / (dLv + dRv + dUv + dDv));
                        STC2VFU(pix[c][0], v0v);
                    }
                }

#endif

                for (; col < width - 2; col += 2) {
                    int indx = row * width + col;
                    pix[1] = (float*)(*rgb[1]) + indx;

                    for (int i = 0; i < 2; c = 2 - c, i++) {
                        pix[c] = (float*)(*rgb[c]) + indx;
                        float dL = 1.f / (1.f + fabsf(pix[1][ -2] - pix[1][0]) + fabsf(pix[c][ 1] - pix[c][ -1]));
                        float dR = 1.f / (1.f + fabsf(pix[1][  2] - pix[1][0]) + fabsf(pix[c][ 1] - pix[c][ -1]));
                        float dU = 1.f / (1.f + fabsf(pix[1][-w2] - pix[1][0]) + fabsf(pix[c][w1] - pix[c][-w1]));
                        float dD = 1.f / (1.f + fabsf(pix[1][ w2] - pix[1][0]) + fabsf(pix[c][w1] - pix[c][-w1]));
                        float v0 = (pix[1][0] + 0.5f - ((pix[1][ -1] - pix[c][ -1]) * dL + (pix[1][  1] - pix[c][  1]) * dR + (pix[1][-w1] - pix[c][-w1]) * dU + (pix[1][ w1] - pix[c][ w1]) * dD ) / (dL + dR + dU + dD));
                        pix[c][0] = CLIP(v0);
                    }
                }
            }

            /* Reinforce integrated red/blue pixels on BLUE/RED pixel locations */
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 2; row < height - 2; row++) {
                int col = 2 + (FC(row, 2) & 1);
                int c = 2 - FC(row, col);
#ifdef __SSE2__
                __m128 dLv, dRv, dUv, dDv, v0v;
                __m128 onev = _mm_set1_ps(1.f);
                __m128 zd5v = _mm_set1_ps(0.5f);
                __m128 c65535v = _mm_set1_ps(65535.f);

                for (; col < width - 8; col += 8) {
                    int indx = row * width + col;
                    pix[0] = (float*)(*rgb[0]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    pix[2] = (float*)(*rgb[2]) + indx;
                    int d = 2 - c;
                    dLv = onev / (onev + vabsf(LC2VFU(pix[d][ -2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dRv = onev / (onev + vabsf(LC2VFU(pix[d][  2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dUv = onev / (onev + vabsf(LC2VFU(pix[d][-w2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    dDv = onev / (onev + vabsf(LC2VFU(pix[d][ w2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    v0v = CLIPV(LC2VFU(pix[1][0]) + zd5v - ((LC2VFU(pix[1][-1]) - LC2VFU(pix[c][-1])) * dLv + (LC2VFU(pix[1][1]) - LC2VFU(pix[c][1])) * dRv + (LC2VFU(pix[1][-w1]) - LC2VFU(pix[c][-w1])) * dUv + (LC2VFU(pix[1][w1]) - LC2VFU(pix[c][w1])) * dDv ) / (dLv + dRv + dUv + dDv));
                    STC2VFU(pix[c][0], v0v);
                }

#endif

                for (; col < width - 2; col += 2) {
                    int indx = row * width + col;
                    pix[0] = (float*)(*rgb[0]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    pix[2] = (float*)(*rgb[2]) + indx;
                    int d = 2 - c;
                    float dL = 1.f / (1.f + fabsf(pix[d][ -2] - pix[d][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dR = 1.f / (1.f + fabsf(pix[d][  2] - pix[d][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dU = 1.f / (1.f + fabsf(pix[d][-w2] - pix[d][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float dD = 1.f / (1.f + fabsf(pix[d][ w2] - pix[d][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float v0 = (pix[1][0] + 0.5f - ((pix[1][ -1] - pix[c][ -1]) * dL + (pix[1][  1] - pix[c][  1]) * dR + (pix[1][-w1] - pix[c][-w1]) * dU + (pix[1][ w1] - pix[c][ w1]) * dD ) / (dL + dR + dU + dD));
                    pix[c][0] = CLIP(v0);
                }
            }
        } // end parallel
    }

    t2e.set();

    if (settings->verbose) {
        printf("Refinement Lee %d usec\n", t2e.etime(t1e));
    }
}
#ifdef __SSE2__
#undef CLIPV
#endif


// Refinement based on EECI demozaicing algorithm by L. Chang and Y.P. Tan
// from "Lassus" : Luis Sanz Rodriguez, adapted by Jacques Desmis - JDC - and Oliver Duis for RawTherapee
// increases the signal to noise ratio (PSNR) # +1 to +2 dB : tested with Dcraw : eg: Lighthouse + AMaZE : whitout refinement:39.96dB, with refinement:41.86 dB
// reduce color artifacts, improves the interpolation
// but it's relatively slow
//
// Should be DISABLED if it decreases image quality by increases some image noise and generates blocky edges
void RawImageSource::refinement_lassus(int PassCount)
{
    // const int PassCount=1;

    // if (settings->verbose) printf("Refinement\n");

    MyTime t1e, t2e;
    t1e.set();
    int u = W, v = 2 * u, w = 3 * u, x = 4 * u, y = 5 * u;
    float (*image)[3];
    image = (float(*)[3]) calloc(W * H, sizeof * image);
#ifdef _OPENMP
    #pragma omp parallel shared(image)
#endif
    {
        // convert red, blue, green to image
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                image[i * W + j][0] = red  [i][j];
                image[i * W + j][1] = green[i][j];
                image[i * W + j][2] = blue [i][j];
            }
        }

        for (int b = 0; b < PassCount; b++) {
            if (plistener) {
                plistener->setProgressStr (M("TP_RAW_DMETHOD_PROGRESSBAR_REFINE"));
                plistener->setProgress ((float)b / PassCount);
            }

            // Reinforce interpolated green pixels on RED/BLUE pixel locations
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 6; row < H - 6; row++) {
                for (int col = 6 + (FC(row, 2) & 1), c = FC(row, col); col < W - 6; col += 2) {
                    float (*pix)[3] = image + row * W + col;

                    // Cubic Spline Interpolation by Li and Randhawa, modified by Luis Sanz Rodriguez

                    float f[4];
                    f[0] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-v][c]) - x0875(pix[0][c]) - x0250(pix[-x][c]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[-w][1])) + fabs(x0875(pix[-w][1]) - x1125(pix[-u][1]) + x0250(pix[-y][1])));
                    f[1] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+2][c]) - x0875(pix[0][c]) - x0250(pix[+4][c]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[+3][1])) + fabs(x0875(pix[+3][1]) - x1125(pix[+1][1]) + x0250(pix[+5][1])));
                    f[2] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-2][c]) - x0875(pix[0][c]) - x0250(pix[-4][c]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[-3][1])) + fabs(x0875(pix[-3][1]) - x1125(pix[-1][1]) + x0250(pix[-5][1])));
                    f[3] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+v][c]) - x0875(pix[0][c]) - x0250(pix[+x][c]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[+w][1])) + fabs(x0875(pix[+w][1]) - x1125(pix[+u][1]) + x0250(pix[+y][1])));

                    float g[4];//CLIREF avoid overflow
                    g[0] = pix[0][c] + (x0875(CLIREF(pix[-u][1] - pix[-u][c])) + x0125(CLIREF(pix[+u][1] - pix[+u][c])));
                    g[1] = pix[0][c] + (x0875(CLIREF(pix[+1][1] - pix[+1][c])) + x0125(CLIREF(pix[-1][1] - pix[-1][c])));
                    g[2] = pix[0][c] + (x0875(CLIREF(pix[-1][1] - pix[-1][c])) + x0125(CLIREF(pix[+1][1] - pix[+1][c])));
                    g[3] = pix[0][c] + (x0875(CLIREF(pix[+u][1] - pix[+u][c])) + x0125(CLIREF(pix[-u][1] - pix[-u][c])));

                    pix[0][1] = (f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]);

                }
            }

            // Reinforce interpolated red/blue pixels on GREEN pixel locations
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 6; row < H - 6; row++) {
                for (int col = 6 + (FC(row, 3) & 1), c = FC(row, col + 1); col < W - 6; col += 2) {
                    float (*pix)[3] = image + row * W + col;

                    for (int i = 0; i < 2; c = 2 - c, i++) {
                        float f[4];
                        f[0] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[-v][1]) - x1125(pix[0][1]) + x0250(pix[-x][1]))) + fabs(pix[u] [c] - pix[-u][c]) + fabs(pix[-w][c] - pix[-u][c]));
                        f[1] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[+2][1]) - x1125(pix[0][1]) + x0250(pix[+4][1]))) + fabs(pix[+1][c] - pix[-1][c]) + fabs(pix[+3][c] - pix[+1][c]));
                        f[2] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[-2][1]) - x1125(pix[0][1]) + x0250(pix[-4][1]))) + fabs(pix[+1][c] - pix[-1][c]) + fabs(pix[-3][c] - pix[-1][c]));
                        f[3] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[+v][1]) - x1125(pix[0][1]) + x0250(pix[+x][1]))) + fabs(pix[u] [c] - pix[-u][c]) + fabs(pix[+w][c] - pix[+u][c]));

                        float g[5];//CLIREF avoid overflow
                        g[0] = CLIREF(pix[-u][1] - pix[-u][c]);
                        g[1] = CLIREF(pix[+1][1] - pix[+1][c]);
                        g[2] = CLIREF(pix[-1][1] - pix[-1][c]);
                        g[3] = CLIREF(pix[+u][1] - pix[+u][c]);
                        g[4] = ((f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]));
                        pix[0][c] = pix[0][1] - (0.65f * g[4] + 0.35f * CLIREF(pix[0][1] - pix[0][c]));
                    }
                }
            }

            // Reinforce integrated red/blue pixels on BLUE/RED pixel locations
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 6; row < H - 6; row++) {
                for (int col = 6 + (FC(row, 2) & 1), c = 2 - FC(row, col), d = 2 - c; col < W - 6; col += 2) {
                    float (*pix)[3] = image + row * W + col;

                    float f[4];
                    f[0] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-v][d]) - x0875(pix[0][d]) - x0250(pix[-x][d]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[-w][1])) + fabs(x0875(pix[-w][1]) - x1125(pix[-u][1]) + x0250(pix[-y][1])));
                    f[1] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+2][d]) - x0875(pix[0][d]) - x0250(pix[+4][d]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[+3][1])) + fabs(x0875(pix[+3][1]) - x1125(pix[+1][1]) + x0250(pix[+5][1])));
                    f[2] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-2][d]) - x0875(pix[0][d]) - x0250(pix[-4][d]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[-3][1])) + fabs(x0875(pix[-3][1]) - x1125(pix[-1][1]) + x0250(pix[-5][1])));
                    f[3] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+v][d]) - x0875(pix[0][d]) - x0250(pix[+x][d]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[+w][1])) + fabs(x0875(pix[+w][1]) - x1125(pix[+u][1]) + x0250(pix[+y][1])));

                    float g[5];
                    g[0] = (x0875((pix[-u][1] - pix[-u][c])) + x0125((pix[-v][1] - pix[-v][c])));
                    g[1] = (x0875((pix[+1][1] - pix[+1][c])) + x0125((pix[+2][1] - pix[+2][c])));
                    g[2] = (x0875((pix[-1][1] - pix[-1][c])) + x0125((pix[-2][1] - pix[-2][c])));
                    g[3] = (x0875((pix[+u][1] - pix[+u][c])) + x0125((pix[+v][1] - pix[+v][c])));

                    g[4] = (f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]);

                    const std::array<float, 9> p = {
                        pix[-u - 1][1] - pix[-u - 1][c],
                        pix[-u + 0][1] - pix[-u + 0][c],
                        pix[-u + 1][1] - pix[-u + 1][c],
                        pix[+0 - 1][1] - pix[+0 - 1][c],
                        pix[+0 + 0][1] - pix[+0 + 0][c],
                        pix[+0 + 1][1] - pix[+0 + 1][c],
                        pix[+u - 1][1] - pix[+u - 1][c],
                        pix[+u + 0][1] - pix[+u + 0][c],
                        pix[+u + 1][1] - pix[+u + 1][c]
                    };

                    const float med = median(p);

                    pix[0][c] = LIM(pix[0][1] - (1.30f * g[4] - 0.30f * (pix[0][1] - pix[0][c])), 0.99f * (pix[0][1] - med), 1.01f * (pix[0][1] - med));

                }
            }

        }

        // put modified values to red, green, blue
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                red  [i][j] = image[i * W + j][0];
                green[i][j] = image[i * W + j][1];
                blue [i][j] = image[i * W + j][2];
            }
        }
    }

    free(image);

    t2e.set();

    if (settings->verbose) {
        printf("Refinement Lassus %d usec\n", t2e.etime(t1e));
    }
}



} /* namespace */
