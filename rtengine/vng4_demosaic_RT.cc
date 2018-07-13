////////////////////////////////////////////////////////////////
//
//          VNG4 demosaic algorithm
// 
// optimized for speed by Ingo Weyrich
//
//
//  vng4_interpolate_RT.cc is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include "rtengine.h"
#include "rawimagesource.h"
#include "rawimagesource_i.h"
#include "../rtgui/multilangmgr.h"
//#define BENCHMARK
#include "StopWatch.h"

namespace rtengine
{
#define fc(row,col) (prefilters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)
typedef unsigned short ushort;
void RawImageSource::vng4_demosaic (const array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue, bool keepGreens)
{
    BENCHFUN
    const signed short int *cp, terms[] = {
        -2, -2, +0, -1, 0, 0x01, -2, -2, +0, +0, 1, 0x01, -2, -1, -1, +0, 0, 0x01,
        -2, -1, +0, -1, 0, 0x02, -2, -1, +0, +0, 0, 0x03, -2, -1, +0, +1, 1, 0x01,
        -2, +0, +0, -1, 0, 0x06, -2, +0, +0, +0, 1, 0x02, -2, +0, +0, +1, 0, 0x03,
        -2, +1, -1, +0, 0, 0x04, -2, +1, +0, -1, 1, 0x04, -2, +1, +0, +0, 0, 0x06,
        -2, +1, +0, +1, 0, 0x02, -2, +2, +0, +0, 1, 0x04, -2, +2, +0, +1, 0, 0x04,
        -1, -2, -1, +0, 0, 0x80, -1, -2, +0, -1, 0, 0x01, -1, -2, +1, -1, 0, 0x01,
        -1, -2, +1, +0, 1, 0x01, -1, -1, -1, +1, 0, 0x88, -1, -1, +1, -2, 0, 0x40,
        -1, -1, +1, -1, 0, 0x22, -1, -1, +1, +0, 0, 0x33, -1, -1, +1, +1, 1, 0x11,
        -1, +0, -1, +2, 0, 0x08, -1, +0, +0, -1, 0, 0x44, -1, +0, +0, +1, 0, 0x11,
        -1, +0, +1, -2, 1, 0x40, -1, +0, +1, -1, 0, 0x66, -1, +0, +1, +0, 1, 0x22,
        -1, +0, +1, +1, 0, 0x33, -1, +0, +1, +2, 1, 0x10, -1, +1, +1, -1, 1, 0x44,
        -1, +1, +1, +0, 0, 0x66, -1, +1, +1, +1, 0, 0x22, -1, +1, +1, +2, 0, 0x10,
        -1, +2, +0, +1, 0, 0x04, -1, +2, +1, +0, 1, 0x04, -1, +2, +1, +1, 0, 0x04,
        +0, -2, +0, +0, 1, 0x80, +0, -1, +0, +1, 1, 0x88, +0, -1, +1, -2, 0, 0x40,
        +0, -1, +1, +0, 0, 0x11, +0, -1, +2, -2, 0, 0x40, +0, -1, +2, -1, 0, 0x20,
        +0, -1, +2, +0, 0, 0x30, +0, -1, +2, +1, 1, 0x10, +0, +0, +0, +2, 1, 0x08,
        +0, +0, +2, -2, 1, 0x40, +0, +0, +2, -1, 0, 0x60, +0, +0, +2, +0, 1, 0x20,
        +0, +0, +2, +1, 0, 0x30, +0, +0, +2, +2, 1, 0x10, +0, +1, +1, +0, 0, 0x44,
        +0, +1, +1, +2, 0, 0x10, +0, +1, +2, -1, 1, 0x40, +0, +1, +2, +0, 0, 0x60,
        +0, +1, +2, +1, 0, 0x20, +0, +1, +2, +2, 0, 0x10, +1, -2, +1, +0, 0, 0x80,
        +1, -1, +1, +1, 0, 0x88, +1, +0, +1, +2, 0, 0x08, +1, +0, +2, -1, 0, 0x40,
        +1, +0, +2, +1, 0, 0x10
    },
    chood[] = { -1, -1, -1, 0, -1, +1, 0, +1, +1, +1, +1, 0, +1, -1, 0, -1 };

    double progress = 0.0;
    const bool plistenerActive = plistener;

    if (plistenerActive) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::VNG4)));
        plistener->setProgress (progress);
    }

    const unsigned prefilters = ri->prefilters;
    const int width = W, height = H;
    constexpr unsigned int colors = 4;
    float (*image)[4];

    image = (float (*)[4]) calloc (height * width, sizeof * image);

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int ii = 0; ii < H; ii++)
        for (int jj = 0; jj < W; jj++) {
            image[ii * W + jj][fc(ii, jj)] = rawData[ii][jj];
        }

    {
        int lcode[16][16][32];
        float mul[16][16][8];
        float csum[16][16][4];

// first linear interpolation
        for (int row = 0; row < 16; row++)
            for (int col = 0; col < 16; col++) {
                int * ip = lcode[row][col];
                int mulcount = 0;
                float sum[4];
                memset (sum, 0, sizeof sum);

                for (int y = -1; y <= 1; y++)
                    for (int x = -1; x <= 1; x++) {
                        int shift = (y == 0) + (x == 0);

                        if (shift == 2) {
                            continue;
                        }

                        int color = fc(row + y, col + x);
                        *ip++ = (width * y + x) * 4 + color;

                        mul[row][col][mulcount] = (1 << shift);
                        *ip++ = color;
                        sum[color] += (1 << shift);
                        mulcount++;
                    }

                int colcount = 0;

                for (unsigned int c = 0; c < colors; c++)
                    if (c != fc(row, col)) {
                        *ip++ = c;
                        csum[row][col][colcount] = sum[c];
                        colcount ++;
                    }
            }

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int row = 1; row < height - 1; row++) {
            for (int col = 1; col < width - 1; col++) {
                float * pix = image[row * width + col];
                int * ip = lcode[row & 15][col & 15];
                float sum[4];
                memset (sum, 0, sizeof sum);

                for (int i = 0; i < 8; i++, ip += 2) {
                    sum[ip[1]] += pix[ip[0]] * mul[row & 15][col & 15][i];
                }

                for (unsigned int i = 0; i < colors - 1; i++, ip++) {
                    pix[ip[0]] = sum[ip[0]] / csum[row & 15][col & 15][i];
                }
            }
        }
    }

    const int prow = 7, pcol = 1;
    int *code[8][2];
    int t, g;
    int * ip = (int *) calloc ((prow + 1) * (pcol + 1), 1280);

    for (int row = 0; row <= prow; row++)   /* Precalculate for VNG */
        for (int col = 0; col <= pcol; col++) {
            code[row][col] = ip;

            for (cp = terms, t = 0; t < 64; t++) {
                int y1 = *cp++;
                int x1 = *cp++;
                int y2 = *cp++;
                int x2 = *cp++;
                int weight = *cp++;
                int grads = *cp++;
                unsigned int color = fc(row + y1, col + x1);

                if (fc(row + y2, col + x2) != color) {
                    continue;
                }

                int diag = (fc(row, col + 1) == color && fc(row + 1, col) == color) ? 2 : 1;

                if (abs(y1 - y2) == diag && abs(x1 - x2) == diag) {
                    continue;
                }

                *ip++ = (y1 * width + x1) * 4 + color;
                *ip++ = (y2 * width + x2) * 4 + color;
                *ip++ = weight;

                for (g = 0; g < 8; g++)
                    if (grads & (1 << g)) {
                        *ip++ = g;
                    }

                *ip++ = -1;
            }

            *ip++ = INT_MAX;

            for (cp = chood, g = 0; g < 8; g++) {
                int y = *cp++;
                int x = *cp++;
                *ip++ = (y * width + x) * 4;
                unsigned int color = fc(row, col);

                if (fc(row + y, col + x) != color && fc(row + y * 2, col + x * 2) == color) {
                    *ip++ = (y * width + x) * 8 + color;
                } else {
                    *ip++ = 0;
                }
            }
        }

    if(plistenerActive) {
        progress = 0.1;
        plistener->setProgress (progress);
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float gval[8], thold, sum[3];
        int g;
        const int progressStep = 64;
        const double progressInc = (0.98 - progress) / ((height - 2) / progressStep);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16) nowait
#endif

        for (int row = 2; row < height - 2; row++) {    /* Do VNG interpolation */
            for (int col = 2; col < width - 2; col++) {
                float * pix = image[row * width + col];
                int color = fc(row, col);
                if (keepGreens && (color & 1)) {
                    green[row][col] = pix[color];
                } else {
                    int * ip = code[row & prow][col & pcol];
                    memset (gval, 0, sizeof gval);

                    while ((g = ip[0]) != INT_MAX) {        /* Calculate gradients */
                        float diff = fabsf(pix[g] - pix[ip[1]]) * (1 << ip[2]);
                        gval[ip[3]] += diff;
                        ip += 4;

                        while ((g = *ip++) != -1) {
                            gval[g] += diff;
                        }
                    }

                    ip++;
                    {
                        float gmin, gmax;
                        gmin = gmax = gval[0];          /* Choose a threshold */

                        for (g = 1; g < 8; g++) {
                            if (gmin > gval[g]) {
                                gmin = gval[g];
                            }

                            if (gmax < gval[g]) {
                                gmax = gval[g];
                            }
                        }

                        thold = gmin + (gmax / 2);
                    }
                    memset (sum, 0, sizeof sum);
                    float t1, t2;
                    t1 = t2 = pix[color];

                    if(color & 1) {
                        int num = 0;

                        for (g = 0; g < 8; g++, ip += 2) {  /* Average the neighbors */
                            if (gval[g] <= thold) {
                                if(ip[1]) {
                                    sum[0] += (t1 + pix[ip[1]]) * 0.5f;
                                }

                                sum[1] += pix[ip[0] + (color ^ 2)];
                                num++;
                            }
                        }

                        t1 += (sum[1] - sum[0]) / num;
                    } else {
                        int num = 0;

                        for (g = 0; g < 8; g++, ip += 2) {  /* Average the neighbors */
                            if (gval[g] <= thold) {
                                sum[1] += pix[ip[0] + 1];
                                sum[2] += pix[ip[0] + 3];

                                if(ip[1]) {
                                    sum[0] += (t1 + pix[ip[1]]) * 0.5f;
                                }

                                num++;
                            }
                        }

                        t1 += (sum[1] - sum[0]) / num;
                        t2 += (sum[2] - sum[0]) / num;
                    }

                    green[row][col] = 0.5f * (t1 + t2);
                }
            }

            if(plistenerActive) {
                if((row % progressStep) == 0)
#ifdef _OPENMP
                    #pragma omp critical (updateprogress)
#endif
                {
                    progress += progressInc;
                    plistener->setProgress (progress);
                }
            }
        }

    }
    free (code[0][0]);
    free (image);

    if(plistenerActive) {
        plistener->setProgress (0.98);
    }

    // Interpolate R and B
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < H; i++) {
        if (i == 0)
            // rm, gm, bm must be recovered
            //interpolate_row_rb_mul_pp (red, blue, NULL, green[i], green[i+1], i, rm, gm, bm, 0, W, 1);
        {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], nullptr, green[i], green[i + 1], i, 1.0, 1.0, 1.0, 0, W, 1);
        } else if (i == H - 1) {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], green[i - 1], green[i], nullptr, i, 1.0, 1.0, 1.0, 0, W, 1);
        } else {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], green[i - 1], green[i], green[i + 1], i, 1.0, 1.0, 1.0, 0, W, 1);
        }
    }
    border_interpolate2(W, H, 3, rawData, red, green, blue);

    if(plistenerActive) {
        plistener->setProgress (1.0);
    }
}
}
