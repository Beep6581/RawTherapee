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

#include "rawimagesource.h"
#include "rawimagesource_i.h"
#include "jaggedarray.h"
#include "rawimage.h"
#include "mytime.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "image8.h"
#include "curves.h"
#include "dfmanager.h"
#include "slicer.h"
#include "rt_math.h"
#include "color.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#include "sleef.c"
#include "opthelper.h"
#include "median.h"
//#define BENCHMARK
#include "StopWatch.h"
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

/*
   Patterned Pixel Grouping Interpolation by Alain Desbiolles
*/
void RawImageSource::ppg_demosaic()
{
    int width = W, height = H;
    int dir[5] = { 1, width, -1, -width, 1 };
    int row, col, diff[2] = {}, guess[2], c, d, i;
    float (*pix)[4];

    float (*image)[4];

    if (plistener) {
        // looks like ppg isn't supported anymore
        //plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::ppg)));
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "xxx"));
        plistener->setProgress (0.0);
    }

    image = (float (*)[4]) calloc (H * W, sizeof * image);

    for (int ii = 0; ii < H; ii++)
        for (int jj = 0; jj < W; jj++) {
            image[ii * W + jj][fc(ii, jj)] = rawData[ii][jj];
        }

    border_interpolate(3, image);

    /*  Fill in the green layer with gradients and pattern recognition: */
    for (row = 3; row < height - 3; row++) {
        for (col = 3 + (FC(row, 3) & 1), c = FC(row, col); col < width - 3; col += 2) {
            pix = image + row * width + col;

            for (i = 0; (d = dir[i]) > 0; i++) {
                guess[i] = (pix[-d][1] + pix[0][c] + pix[d][1]) * 2
                           - pix[-2 * d][c] - pix[2 * d][c];
                diff[i] = ( ABS(pix[-2 * d][c] - pix[ 0][c]) +
                            ABS(pix[ 2 * d][c] - pix[ 0][c]) +
                            ABS(pix[  -d][1] - pix[ d][1]) ) * 3 +
                          ( ABS(pix[ 3 * d][1] - pix[ d][1]) +
                            ABS(pix[-3 * d][1] - pix[-d][1]) ) * 2;
            }

            d = dir[i = diff[0] > diff[1]];
            pix[0][1] = median(static_cast<float>(guess[i] >> 2), pix[d][1], pix[-d][1]);
        }

        if(plistener) {
            plistener->setProgress(0.33 * row / (height - 3));
        }
    }

    /*  Calculate red and blue for each green pixel:        */
    for (row = 1; row < height - 1; row++) {
        for (col = 1 + (FC(row, 2) & 1), c = FC(row, col + 1); col < width - 1; col += 2) {
            pix = image + row * width + col;

            for (i = 0; (d = dir[i]) > 0; c = 2 - c, i++)
                pix[0][c] = CLIP(0.5 * (pix[-d][c] + pix[d][c] + 2 * pix[0][1]
                                        - pix[-d][1] - pix[d][1]) );
        }

        if(plistener) {
            plistener->setProgress(0.33 + 0.33 * row / (height - 1));
        }
    }

    /*  Calculate blue for red pixels and vice versa:       */
    for (row = 1; row < height - 1; row++) {
        for (col = 1 + (FC(row, 1) & 1), c = 2 - FC(row, col); col < width - 1; col += 2) {
            pix = image + row * width + col;

            for (i = 0; (d = dir[i] + dir[i + 1]) > 0; i++) {
                diff[i] = ABS(pix[-d][c] - pix[d][c]) +
                          ABS(pix[-d][1] - pix[0][1]) +
                          ABS(pix[ d][1] - pix[0][1]);
                guess[i] = pix[-d][c] + pix[d][c] + 2 * pix[0][1]
                           - pix[-d][1] - pix[d][1];
            }

            if (diff[0] != diff[1]) {
                pix[0][c] = CLIP(guess[diff[0] > diff[1]] / 2);
            } else {
                pix[0][c] = CLIP((guess[0] + guess[1]) / 4);
            }
        }

        if(plistener) {
            plistener->setProgress(0.67 + 0.33 * row / (height - 1));
        }
    }

    red(W, H);

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            red[i][j] = image[i * W + j][0];
        }

    green(W, H);

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            green[i][j] = image[i * W + j][1];
        }

    blue(W, H);

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            blue[i][j] = image[i * W + j][2];
        }

    free (image);
}

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

// Joint Demosaicing and Denoising using High Order Interpolation Techniques
// Revision 0.9.1a - 09/02/2010 - Contact info: luis.sanz.rodriguez@gmail.com
// Copyright Luis Sanz Rodriguez 2010
// Adapted to RawTherapee by Jacques Desmis 3/2013

void RawImageSource::jdl_interpolate_omp()  // from "Lassus"
{
    int width = W, height = H;
    int row, col, c, d, i, u = width, v = 2 * u, w = 3 * u, x = 4 * u, y = 5 * u, z = 6 * u, indx, (*dif)[2], (*chr)[2];
    float f[4], g[4];
    float (*image)[4];
    image = (float (*)[4]) calloc (width * height, sizeof * image);
    dif = (int (*)[2]) calloc(width * height, sizeof * dif);
    chr = (int (*)[2]) calloc(width * height, sizeof * chr);

    if (plistener) {
        // this function seems to be unused
        //plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::jdl)));
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "xxx"));
        plistener->setProgress (0.0);
    }

#ifdef _OPENMP
    #pragma omp parallel default(none) shared(image,width,height,u,w,v,y,x,z,dif,chr) private(row,col,f,g,indx,c,d,i)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int ii = 0; ii < height; ii++)
            for (int jj = 0; jj < width; jj++) {
                image[ii * width + jj][fc(ii, jj)] = rawData[ii][jj];
            }

        border_interpolate(6, image);

#ifdef _OPENMP
        #pragma omp for
#endif

        for (row = 5; row < height - 5; row++)
            for (col = 5 + (FC(row, 1) & 1), indx = row * width + col, c = FC(row, col); col < u - 5; col += 2, indx += 2) {
                f[0] = 1.f + abs(image[indx - u][1] - image[indx - w][1]) + abs(image[indx - u][1] - image[indx + u][1]) + abs(image[indx][c] - image[indx - v][c]) + abs(image[indx - v][c] - image[indx - x][c]);
                f[1] = 1.f + abs(image[indx + 1][1] - image[indx + 3][1]) + abs(image[indx + 1][1] - image[indx - 1][1]) + abs(image[indx][c] - image[indx + 2][c]) + abs(image[indx + 2][c] - image[indx + 4][c]);
                f[2] = 1.f + abs(image[indx - 1][1] - image[indx - 3][1]) + abs(image[indx - 1][1] - image[indx + 1][1]) + abs(image[indx][c] - image[indx - 2][c]) + abs(image[indx - 2][c] - image[indx - 4][c]);
                f[3] = 1.f + abs(image[indx + u][1] - image[indx + w][1]) + abs(image[indx + u][1] - image[indx - u][1]) + abs(image[indx][c] - image[indx + v][c]) + abs(image[indx + v][c] - image[indx + x][c]);
                g[0] = CLIP((22.f * image[indx - u][1] + 22.f * image[indx - w][1] + 2.f * image[indx - y][1] + 2.f * image[indx + u][1] + 40.f * image[indx][c] - 32.f * image[indx - v][c] - 8.f * image[indx - x][c]) / 48.f);
                g[1] = CLIP((22.f * image[indx + 1][1] + 22.f * image[indx + 3][1] + 2.f * image[indx + 5][1] + 2.f * image[indx - 1][1] + 40.f * image[indx][c] - 32.f * image[indx + 2][c] - 8.f * image[indx + 4][c]) / 48.f);
                g[2] = CLIP((22.f * image[indx - 1][1] + 22.f * image[indx - 3][1] + 2.f * image[indx - 5][1] + 2.f * image[indx + 1][1] + 40.f * image[indx][c] - 32.f * image[indx - 2][c] - 8.f * image[indx - 4][c]) / 48.f);
                g[3] = CLIP((22.f * image[indx + u][1] + 22.f * image[indx + w][1] + 2.f * image[indx + y][1] + 2.f * image[indx - u][1] + 40.f * image[indx][c] - 32.f * image[indx + v][c] - 8.f * image[indx + x][c]) / 48.f);
                dif[indx][0] = CLIP((f[3] * g[0] + f[0] * g[3]) / (f[0] + f[3])) - image[indx][c];
                dif[indx][1] = CLIP((f[2] * g[1] + f[1] * g[2]) / (f[1] + f[2])) - image[indx][c];
            }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (row = 6; row < height - 6; row++)
            for (col = 6 + (FC(row, 2) & 1), indx = row * width + col, c = FC(row, col) / 2; col < u - 6; col += 2, indx += 2) {
                f[0] = 1.f + 78.f * SQR((float)dif[indx][0]) + 69.f * (SQR((float) dif[indx - v][0]) + SQR((float)dif[indx + v][0])) + 51.f * (SQR((float)dif[indx - x][0]) + SQR((float)dif[indx + x][0])) + 21.f * (SQR((float)dif[indx - z][0]) + SQR((float)dif[indx + z][0])) - 6.f * SQR((float)dif[indx - v][0] + dif[indx][0] + dif[indx + v][0]) - 10.f * (SQR((float)dif[indx - x][0] + dif[indx - v][0] + dif[indx][0]) + SQR((float)dif[indx][0] + dif[indx + v][0] + dif[indx + x][0])) - 7.f * (SQR((float)dif[indx - z][0] + dif[indx - x][0] + dif[indx - v][0]) + SQR((float)dif[indx + v][0] + dif[indx + x][0] + dif[indx + z][0]));
                f[1] = 1.f + 78.f * SQR((float)dif[indx][1]) + 69.f * (SQR((float)dif[indx - 2][1]) + SQR((float)dif[indx + 2][1])) + 51.f * (SQR((float)dif[indx - 4][1]) + SQR((float)dif[indx + 4][1])) + 21.f * (SQR((float)dif[indx - 6][1]) + SQR((float)dif[indx + 6][1])) - 6.f * SQR((float)dif[indx - 2][1] + dif[indx][1] + dif[indx + 2][1]) - 10.f * (SQR((float)dif[indx - 4][1] + dif[indx - 2][1] + dif[indx][1]) + SQR((float)dif[indx][1] + dif[indx + 2][1] + dif[indx + 4][1])) - 7.f * (SQR((float)dif[indx - 6][1] + dif[indx - 4][1] + dif[indx - 2][1]) + SQR((float)dif[indx + 2][1] + dif[indx + 4][1] + dif[indx + 6][1]));
                g[0] = median(0.725f * dif[indx][0] + 0.1375f * dif[indx - v][0] + 0.1375f * dif[indx + v][0], static_cast<float>(dif[indx - v][0]), static_cast<float>(dif[indx + v][0]));
                g[1] = median(0.725f * dif[indx][1] + 0.1375f * dif[indx - 2][1] + 0.1375f * dif[indx + 2][1], static_cast<float>(dif[indx - 2][1]), static_cast<float>(dif[indx + 2][1]));
                chr[indx][c] = (f[1] * g[0] + f[0] * g[1]) / (f[0] + f[1]);
            }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (row = 6; row < height - 6; row++)
            for (col = 6 + (FC(row, 2) & 1), indx = row * width + col, c = 1 - FC(row, col) / 2, d = 2 * c; col < u - 6; col += 2, indx += 2) {
                f[0] = 1.f / (float)(1.f + fabs((float)chr[indx - u - 1][c] - chr[indx + u + 1][c]) + fabs((float)chr[indx - u - 1][c] - chr[indx - w - 3][c]) + fabs((float)chr[indx + u + 1][c] - chr[indx - w - 3][c]));
                f[1] = 1.f / (float)(1.f + fabs((float)chr[indx - u + 1][c] - chr[indx + u - 1][c]) + fabs((float)chr[indx - u + 1][c] - chr[indx - w + 3][c]) + fabs((float)chr[indx + u - 1][c] - chr[indx - w + 3][c]));
                f[2] = 1.f / (float)(1.f + fabs((float)chr[indx + u - 1][c] - chr[indx - u + 1][c]) + fabs((float)chr[indx + u - 1][c] - chr[indx + w + 3][c]) + fabs((float)chr[indx - u + 1][c] - chr[indx + w - 3][c]));
                f[3] = 1.f / (float)(1.f + fabs((float)chr[indx + u + 1][c] - chr[indx - u - 1][c]) + fabs((float)chr[indx + u + 1][c] - chr[indx + w - 3][c]) + fabs((float)chr[indx - u - 1][c] - chr[indx + w + 3][c]));
                g[0] = median(chr[indx - u - 1][c], chr[indx - w - 1][c], chr[indx - u - 3][c]);
                g[1] = median(chr[indx - u + 1][c], chr[indx - w + 1][c], chr[indx - u + 3][c]);
                g[2] = median(chr[indx + u - 1][c], chr[indx + w - 1][c], chr[indx + u - 3][c]);
                g[3] = median(chr[indx + u + 1][c], chr[indx + w + 1][c], chr[indx + u + 3][c]);
                chr[indx][c] = (f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]);
                image[indx][1] = CLIP(image[indx][2 - d] + chr[indx][1 - c]);
                image[indx][d] = CLIP(image[indx][1] - chr[indx][c]);
            }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (row = 6; row < height - 6; row++)
            for (col = 6 + (FC(row, 1) & 1), indx = row * width + col, c = FC(row, col + 1) / 2, d = 2 * c; col < u - 6; col += 2, indx += 2)
                for(i = 0; i <= 1; c = 1 - c, d = 2 * c, i++) {
                    f[0] = 1.f / (float)(1.f + fabs((float)chr[indx - u][c] - chr[indx + u][c]) + fabs((float)chr[indx - u][c] - chr[indx - w][c]) + fabs((float)chr[indx + u][c] - chr[indx - w][c]));
                    f[1] = 1.f / (float)(1.f + fabs((float)chr[indx + 1][c] - chr[indx - 1][c]) + fabs((float)chr[indx + 1][c] - chr[indx + 3][c]) + fabs((float)chr[indx - 1][c] - chr[indx + 3][c]));
                    f[2] = 1.f / (float)(1.f + fabs((float)chr[indx - 1][c] - chr[indx + 1][c]) + fabs((float)chr[indx - 1][c] - chr[indx - 3][c]) + fabs((float)chr[indx + 1][c] - chr[indx - 3][c]));
                    f[3] = 1.f / (float)(1.f + fabs((float)chr[indx + u][c] - chr[indx - u][c]) + fabs((float)chr[indx + u][c] - chr[indx + w][c]) + fabs((float)chr[indx - u][c] - chr[indx + w][c]));
                    g[0] = 0.875f * chr[indx - u][c] + 0.125f * chr[indx - w][c];
                    g[1] = 0.875f * chr[indx + 1][c] + 0.125f * chr[indx + 3][c];
                    g[2] = 0.875f * chr[indx - 1][c] + 0.125f * chr[indx - 3][c];
                    g[3] = 0.875f * chr[indx + u][c] + 0.125f * chr[indx + w][c];
                    image[indx][d] = CLIP(image[indx][1] - (f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]));
                }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int ii = 0; ii < height; ii++) {
            for (int jj = 0; jj < width; jj++) {
                red[ii][jj]  = CLIP(image[ii * width + jj][0]);
                green[ii][jj]    = CLIP(image[ii * width + jj][1]);
                blue[ii][jj]     = CLIP(image[ii * width + jj][2]);
            }
        }
    } // End of parallelization
    free (image);
    free(dif);
    free(chr);
}

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

}