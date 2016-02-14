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
#include <climits>
#include <cstdio>
#include "rawimagesource.h"
#include "rawimagesource_i.h"
#include "jaggedarray.h"
#include "../rtgui/options.h"

namespace rtengine
{

// computes highlight recovery multipliers. Needs a possibly downscaled image where
// the highlights are indicated by INT_MAX
void hlmultipliers (int** rec[3], int max_3[3], int dh, int dw)
{

    // STEP I. recover color with two-color information
    int phase = -1;
    int k = 0;

    for (k = 0; k < 1000; k++) {
        int changed = 0;

        for (int i = 1; i < dh - 1; i++)
            for (int j = 1; j < dw - 1; j++) {
                int co, c1, c2;

//                  if (phase==2)
//                    phase++;
//                if (phase>0)// && phase!=4)
//                    continue;
                if (phase == -1 || phase == 0 || phase == 2) {
                    if (rec[0][i][j] == INT_MAX && rec[1][i][j] != INT_MAX && rec[1][i][j] >= 0 && rec[2][i][j] != INT_MAX && rec[2][i][j] >= 0) {
                        co = 0;
                        c1 = 1;
                        c2 = 2;
                    } else if (rec[1][i][j] == INT_MAX && rec[0][i][j] != INT_MAX && rec[0][i][j] >= 0 && rec[2][i][j] != INT_MAX && rec[2][i][j] >= 0) {
                        co = 1;
                        c1 = 0;
                        c2 = 2;
                    } else if (rec[2][i][j] == INT_MAX && rec[1][i][j] != INT_MAX && rec[1][i][j] >= 0 && rec[0][i][j] != INT_MAX && rec[0][i][j] >= 0) {
                        co = 2;
                        c1 = 1;
                        c2 = 0;
                    } else {
                        continue;
                    }

                    double ratio[2] = {0.0, 0.0};
                    int count = 0;
                    double rato = (double)rec[c1][i][j] / rec[c2][i][j];
                    double arato = 0.0;

                    if (phase == 2) {
                        for (int x = -1; x <= 1; x++)
                            for (int y = -1; y <= 1; y++) {
                                // average m/c color ratios in the surrounding pixels
                                if (rec[co][i + x][j + y] >= 0 && rec[co][i + x][j + y] != INT_MAX && rec[c1][i + x][j + y] >= 0 && rec[c1][i + x][j + y] != INT_MAX && rec[c2][i + x][j + y] > 0 && rec[c2][i + x][j + y] != INT_MAX) {
                                    double ratt = (double)rec[c1][i + x][j + y] / rec[c2][i + x][j + y];

                                    if (ratt > rato * 1.2 || ratt < rato / 1.2 || rec[co][i + x][j + y] < max_3[co] * 1 / 2) {
                                        continue;
                                    }

                                    ratio[0] += (double)rec[c1][i + x][j + y] / rec[co][i + x][j + y];
                                    ratio[1] += (double)rec[c2][i + x][j + y] / rec[co][i + x][j + y];
                                    count++;
                                }
                            }
                    } else if (phase == -1) {
                        for (int x = -1; x <= 1; x++)
                            for (int y = -1; y <= 1; y++) {
                                // average m/c color ratios in the surrounding pixels
                                if (rec[co][i + x][j + y] >= 0 && rec[co][i + x][j + y] != INT_MAX && rec[c1][i + x][j + y] >= 0 && rec[c1][i + x][j + y] != INT_MAX && rec[c2][i + x][j + y] > 0 && rec[c2][i + x][j + y] != INT_MAX) {
                                    double ratt = (double)rec[c1][i + x][j + y] / rec[c2][i + x][j + y];

                                    if (ratt > rato * 1.05 || ratt < rato / 1.05 || rec[co][i + x][j + y] < max_3[co] * 4 / 5) {
                                        continue;
                                    }

                                    arato += ratt;
                                    ratio[0] += (double)rec[c1][i + x][j + y] / rec[co][i + x][j + y];
                                    ratio[1] += (double)rec[c2][i + x][j + y] / rec[co][i + x][j + y];

                                    count++;
                                }
                            }
                    } else {
                        for (int x = -1; x <= 1; x++)
                            for (int y = -1; y <= 1; y++) {
                                // average m/c color ratios in the surrounding pixels
                                if (rec[co][i + x][j + y] >= 0 && rec[co][i + x][j + y] != INT_MAX && rec[c1][i + x][j + y] >= 0 && rec[c1][i + x][j + y] != INT_MAX && rec[c2][i + x][j + y] > 0 && rec[c2][i + x][j + y] != INT_MAX) {
                                    double ratt = (double)rec[c1][i + x][j + y] / rec[c2][i + x][j + y];

                                    if (ratt > rato * 1.1 || ratt < rato / 1.1 || rec[co][i + x][j + y] < max_3[co] * 3 / 4) {
                                        continue;
                                    }

                                    arato += ratt;
                                    ratio[0] += (double)rec[c1][i + x][j + y] / rec[co][i + x][j + y];
                                    ratio[1] += (double)rec[c2][i + x][j + y] / rec[co][i + x][j + y];

                                    count++;
                                }
                            }
                    }

                    // compute new pixel values from the surrounding color ratios
                    if (count > 1) { //(phase==0 && count>1) || (phase==2 && count>1)) {
                        rec[co][i][j] = -(int)((rec[c1][i][j] / ratio[0] * count + rec[c2][i][j] / ratio[1] * count) / 2);
                        changed++;
                    }
                } else if (phase == 1 || phase == 3) {
                    if (rec[0][i][j] == INT_MAX && rec[1][i][j] == INT_MAX && rec[2][i][j] != INT_MAX && rec[2][i][j] >= 0) {
                        co = 2;
                        c1 = 0;
                        c2 = 1;
                    } else if (rec[0][i][j] == INT_MAX && rec[2][i][j] == INT_MAX && rec[1][i][j] != INT_MAX && rec[1][i][j] >= 0) {
                        co = 1;
                        c1 = 0;
                        c2 = 2;
                    } else if (rec[1][i][j] == INT_MAX && rec[2][i][j] == INT_MAX && rec[0][i][j] != INT_MAX && rec[0][i][j] >= 0) {
                        co = 0;
                        c1 = 1;
                        c2 = 2;
                    } else {
                        continue;
                    }

                    double ratio[2] = {0.0, 0.0};
                    int count[2] = {0, 0};

                    for (int x = -1; x <= 1; x++)
                        for (int y = -1; y <= 1; y++) {
                            // average m/c color ratios in the surrounding pixels
                            if (rec[co][i + x][j + y] >= 0 && rec[co][i + x][j + y] != INT_MAX && rec[c1][i + x][j + y] > 0 && rec[c1][i + x][j + y] != INT_MAX) {
                                if ((phase == 1 && rec[c1][i + x][j + y] < max_3[c1] * 3 / 4) || (phase == 3 && rec[c1][i + x][j + y] < max_3[c1] * 1 / 2)) {
                                    continue;
                                }

                                ratio[0] += (double)rec[co][i + x][j + y] / rec[c1][i + x][j + y];
                                count[0] ++;
                            }

                            if (rec[co][i + x][j + y] >= 0 && rec[co][i + x][j + y] != INT_MAX && rec[c2][i + x][j + y] > 0 && rec[c2][i + x][j + y] != INT_MAX) {
                                if ((phase == 1 && rec[c2][i + x][j + y] < max_3[c2] * 3 / 4) || (phase == 3 && rec[c2][i + x][j + y] < max_3[c2] * 1 / 2))
//                                if (/*phase!=3 && */rec[c2][i+x][j+y]<max[c2]*3/4)
                                {
                                    continue;
                                }

                                ratio[1] += (double)rec[co][i + x][j + y] / rec[c2][i + x][j + y];
                                count[1] ++;
                            }
                        }

                    // compute new pixel values from the surrounding color ratios
                    if ((phase == 1 && count[0] > 2) || (phase == 3 && count[0] > 1)) {
                        rec[c1][i][j] = - (int) ((double)rec[co][i][j] / ratio[0] * count[0]);
                        changed++;
                    }

                    if ((phase == 1 && count[1] > 2) || (phase == 3 && count[1] > 1)) {
                        rec[c2][i][j] = - (int) ((double)rec[co][i][j] / ratio[1] * count[1]);
                        changed++;
                    }
                } else {
                    int val = 0;
                    int num = 0;

                    for (int c = 0; c < 3; c++)
                        if (rec[c][i][j] != INT_MAX) {
                            val += rec[c][i][j];
                            num++;
                        }

                    if (num < 3 && num > 0) {
                        for (int c = 0; c < 3; c++) {
                            rec[c][i][j] = val / num;
                        }
                    }
                }
            }


        bool change = false;

        for (int i = 1; i < dh - 1; i++)
            for (int j = 1; j < dw - 1; j++)
                for (int c = 0; c < 3; c++) {
                    if (rec[c][i][j] < 0) {
                        rec[c][i][j] = -rec[c][i][j];
                        change = true;
                    }
                }

        if (!change && phase < 4) {
            phase++;

            if( options.rtSettings.verbose ) {
                printf ("phc %d: %d\n", phase, k);
            }
        } else if (!change) {
            break;
        }

        if (k % 20 == 0 && options.rtSettings.verbose ) {
            printf ("changed %d\n", changed);
        }
    }

    if( options.rtSettings.verbose ) {
        printf ("Highlight recovery ends in %d iterations\n", k);
    }

    int maxval = max(max_3[0], max_3[1], max_3[2]);

    for (int i = 0; i < dh; i++)
        for (int j = 0; j < dw; j++)
            if (rec[0][i][j] == INT_MAX || rec[1][i][j] == INT_MAX || rec[2][i][j] == INT_MAX) {
                rec[0][i][j] = maxval;
                rec[1][i][j] = maxval;
                rec[2][i][j] = maxval;

            }
}

void RawImageSource::HLRecovery_ColorPropagation (float* red, float* green, float* blue, int i, int sx1, int width, int skip)
{

    int blr = (i + HR_SCALE / 2) / HR_SCALE - 1;

    if (blr < 0 || blr >= H / HR_SCALE - 2) {
        return;
    }

    double mr1 = 1.0 - ((double)((i + HR_SCALE / 2) % HR_SCALE) / HR_SCALE + 0.5 / HR_SCALE);
    int maxcol = W / HR_SCALE - 2;

    for (int j = sx1, jx = 0; jx < width; j += skip, jx++) {
        if (needhr[i][j]) {
            int blc = (j + HR_SCALE / 2) / HR_SCALE - 1;

            if (blc < 0 || blc >= maxcol) {
                continue;
            }

            double mc1 = 1.0 - ((double)((j + HR_SCALE / 2) % HR_SCALE) / HR_SCALE + 0.5 / HR_SCALE);
            double mulr = mr1 * mc1 * hrmap[0][blr][blc] + mr1 * (1.0 - mc1) * hrmap[0][blr][blc + 1] + (1.0 - mr1) * mc1 * hrmap[0][blr + 1][blc] + (1.0 - mr1) * (1.0 - mc1) * hrmap[0][blr + 1][blc + 1];
            double mulg = mr1 * mc1 * hrmap[1][blr][blc] + mr1 * (1.0 - mc1) * hrmap[1][blr][blc + 1] + (1.0 - mr1) * mc1 * hrmap[1][blr + 1][blc] + (1.0 - mr1) * (1.0 - mc1) * hrmap[1][blr + 1][blc + 1];
            double mulb = mr1 * mc1 * hrmap[2][blr][blc] + mr1 * (1.0 - mc1) * hrmap[2][blr][blc + 1] + (1.0 - mr1) * mc1 * hrmap[2][blr + 1][blc] + (1.0 - mr1) * (1.0 - mc1) * hrmap[2][blr + 1][blc + 1];
            red[jx] = (red[jx] * mulr);
            green[jx] = (green[jx] * mulg);
            blue[jx] = (blue[jx] * mulb);
        } else {
            red[jx] = (red[jx]);
            green[jx] = (green[jx]);
            blue[jx] = (blue[jx]);
        }
    }
}

void RawImageSource::updateHLRecoveryMap_ColorPropagation ()
{

    // detect maximal pixel values
    float* red = new float[W];
    float* blue = new float[W];
    int maxr = 0, maxg = 0, maxb = 0;

    for (int i = 32; i < H - 32; i++) {
        interpolate_row_rb (red, blue, green[i - 1], green[i], green[i + 1], i);

        for (int j = 32; j < W - 32; j++) {
            if ((ri->ISRED(i, j)   || ri->getSensorType() != ST_BAYER) && red[j] > maxr) {
                maxr = red[j];
            }

            if ((ri->ISGREEN(i, j) || ri->getSensorType() != ST_BAYER) && green[i][j] > maxg) {
                maxg = green[i][j];
            }

            if ((ri->ISBLUE(i, j)  || ri->getSensorType() != ST_BAYER) && blue[j] > maxb) {
                maxb = blue[j];
            }
        }
    }

    delete [] red;
    delete [] blue;

    maxr = maxr * 19 / 20;
    maxg = maxg * 19 / 20;
    maxb = maxb * 19 / 20;
    max_3[0] = maxr;
    max_3[1] = maxg;
    max_3[2] = maxb;

    // downscale image
    int dw = W / HR_SCALE;
    int dh = H / HR_SCALE;
    Image16* ds = new Image16 (dw, dh);

    // overburnt areas
    int** rec[3];

    for (int i = 0; i < 3; i++) {
        rec[i] = allocJaggedArray<int> (dw, dh);
    }

    float* reds[HR_SCALE];
    float* blues[HR_SCALE];

    for (int i = 0; i < HR_SCALE; i++) {
        reds[i] = new float[W];
        blues[i] = new float[W];
    }

    if (needhr) {
        freeJaggedArray<char>(needhr);
    }

    needhr = allocJaggedArray<char> (W, H);

    for (int i = 0; i < dh; i++) {
        for (int j = 0; j < HR_SCALE; j++) {
            interpolate_row_rb (reds[j], blues[j], green[HR_SCALE * i + j - 1], green[HR_SCALE * i + j], green[HR_SCALE * i + j + 1], HR_SCALE * i + j);

            for (int k = 0; k < W; k++)
                if (reds[j][k] >= max_3[0] || green[HR_SCALE * i + j][k] >= max_3[1] || blues[j][k] >= max_3[2]) {
                    needhr[HR_SCALE * i + j][k] = 1;
                } else {
                    needhr[HR_SCALE * i + j][k] = 0;
                }
        }

        for (int j = 0; j < dw; j++) {
            int sumr = 0;
            int cr = 0;
            int sumg = 0;
            int cg = 0;
            int sumb = 0;
            int cb = 0;

            for (int x = 0; x < HR_SCALE; x++)
                for (int y = 0; y < HR_SCALE; y++) {
                    int ix = HR_SCALE * i + x;
                    int jy = HR_SCALE * j + y;
                    sumr += reds[x][jy];

                    if (reds[x][jy] < maxr) {
                        cr++;
                    }

                    sumg += green[ix][jy];

                    if (green[ix][jy] < maxg) {
                        cg++;
                    }

                    sumb += blues[x][jy];

                    if (blues[x][jy] < maxb) {
                        cb++;
                    }
                }

            if (cr < HR_SCALE * HR_SCALE) {
                rec[0][i][j] = INT_MAX;
            } else {
                rec[0][i][j] = sumr / HR_SCALE / HR_SCALE;
            }

            if (cg < HR_SCALE * HR_SCALE) {
                rec[1][i][j] = INT_MAX;
            } else {
                rec[1][i][j] = sumg / HR_SCALE / HR_SCALE;
            }

            if (cb < HR_SCALE * HR_SCALE) {
                rec[2][i][j] = INT_MAX;
            } else {
                rec[2][i][j] = sumb / HR_SCALE / HR_SCALE;
            }

            ds->r(i, j) = sumr / HR_SCALE / HR_SCALE;
            ds->g(i, j) = sumg / HR_SCALE / HR_SCALE;
            ds->b(i, j) = sumb / HR_SCALE / HR_SCALE;
        }
    }

    for (int i = 0; i < HR_SCALE; i++) {
        delete [] reds[i];
        delete [] blues[i];
    }

    hlmultipliers (rec, max_3, dh, dw);

    if (hrmap[0] != NULL) {
        freeJaggedArray<float> (hrmap[0]);
        freeJaggedArray<float> (hrmap[1]);
        freeJaggedArray<float> (hrmap[2]);
    }

    hrmap[0] = allocJaggedArray<float> (dw, dh);
    hrmap[1] = allocJaggedArray<float> (dw, dh);
    hrmap[2] = allocJaggedArray<float> (dw, dh);

    for (int i = 0; i < dh; i++)
        for (int j = 0; j < dw; j++) {
            hrmap[0][i][j] = ds->r(i, j) > 0 ? (double)rec[0][i][j] / ds->r(i, j) : 1.0;
            hrmap[1][i][j] = ds->g(i, j) > 0 ? (double)rec[1][i][j] / ds->g(i, j) : 1.0;
            hrmap[2][i][j] = ds->b(i, j) > 0 ? (double)rec[2][i][j] / ds->b(i, j) : 1.0;
        }

    delete ds;

    freeJaggedArray<int> (rec[0]);
    freeJaggedArray<int> (rec[1]);
    freeJaggedArray<int> (rec[2]);
}

}

