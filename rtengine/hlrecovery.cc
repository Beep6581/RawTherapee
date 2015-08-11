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

namespace rtengine
{

template<class T> T** allocArray (int W, int H)
{

    T** t = new T*[H];

    for (int i = 0; i < H; i++) {
        t[i] = new T[W];
    }

    return t;
}

void RawImageSource::updateHLRecoveryMap (bool needred, bool needgreen, bool needblue, bool full)
{

    // detect maximal pixel values
    unsigned short* red = new unsigned short[W];
    unsigned short* blue = new unsigned short[W];
    int maxr = 0, maxg = 0, maxb = 0;

    for (int i = 32; i < H - 32; i++) {
        interpolate_row_rb (red, blue, green[i - 1], green[i], green[i + 1], i);

        for (int j = 32; j < W - 32; j++) {
            if (red[j] > maxr) {
                maxr = red[j];
            }

            if (green[i][j] > maxg) {
                maxg = green[i][j];
            }

            if (blue[j] > maxb) {
                maxb = blue[j];
            }
        }
    }

    delete [] red;
    delete [] blue;

    maxr = maxr * 19 / 20;
    maxg = maxg * 19 / 20;
    maxb = maxb * 19 / 20;
    int max[3];
    max[0] = maxr;
    max[1] = maxg;
    max[2] = maxb;

    if( options.rtSettings.verbose ) {
        printf ("HLRecoveryMap Maximum: R: %d, G: %d, B: %d\n", maxr, maxg, maxb);
    }

    // downscale image
    int dw = W / SCALE;
    int dh = H / SCALE;
    Image16* ds = new Image16 (dw, dh);

    // overburnt areas
    int** rec[3];

    for (int i = 0; i < 3; i++) {
        rec[i] = allocArray<int> (dw, dh);
    }

    unsigned short* reds[SCALE];
    unsigned short* blues[SCALE];

    for (int i = 0; i < SCALE; i++) {
        reds[i] = new unsigned short[W];
        blues[i] = new unsigned short[W];
    }

    for (int i = 0; i < dh; i++) {
        for (int j = 0; j < SCALE; j++) {
            interpolate_row_rb (reds[j], blues[j], green[SCALE * i + j - 1], green[SCALE * i + j], green[SCALE * i + j + 1], SCALE * i + j);
        }

        for (int j = 0; j < dw; j++) {
            int sumr = 0;
            int cr = 0;
            int sumg = 0;
            int cg = 0;
            int sumb = 0;
            int cb = 0;

            for (int x = 0; x < SCALE; x++)
                for (int y = 0; y < SCALE; y++) {
                    int ix = SCALE * i + x;
                    int jy = SCALE * j + y;
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

            if (cr < SCALE * SCALE && needred) {
                rec[0][i][j] = INT_MAX;
            } else {
                rec[0][i][j] = sumr / SCALE / SCALE;
            }

            if (cg < SCALE * SCALE && needgreen) {
                rec[1][i][j] = INT_MAX;
            } else {
                rec[1][i][j] = sumg / SCALE / SCALE;
            }

            if (cb < SCALE * SCALE && needblue) {
                rec[2][i][j] = INT_MAX;
            } else {
                rec[2][i][j] = sumb / SCALE / SCALE;
            }

            ds->r(i, j) = sumr / SCALE / SCALE;
            ds->g(i, j) = sumg / SCALE / SCALE;
            ds->b(i, j) = sumb / SCALE / SCALE;
        }
    }

    for (int i = 0; i < SCALE; i++) {
        delete [] reds[i];
        delete [] blues[i];
    }


    // STEP I. recover color from the partially lost areas
    bool phase2 = false;

    for (int k = 0; k < 400; k++) {
        if (k > 200) {
            phase2 = true;
        }

        for (int i = 1; i < dh - 1; i++)
            for (int j = 1; j < dw - 1; j++) {
                for (int c = 0; c < 3; c++) {
                    // if channel c is lost
                    if (rec[c][i][j] == INT_MAX) {
                        double ratio[2] = {0.0, 0.0};
                        double w[2] = {0.0, 0.0};
                        int count[2] = {0, 0};
                        int ix = 0;

                        for (int m = 0; m < 3; m++) {
                            if (m == c) {
                                continue;
                            }

                            // if channel m is not lost at this point (or already recovered)
                            if (rec[m][i][j] != INT_MAX && rec[m][i][j] >= 0) {
                                for (int x = -1; x <= 1; x++)
                                    for (int y = -1; y <= 1; y++)

                                        // average m/c color ratios in the surrounding pixels
                                        if (rec[m][i + x][j + y] >= 0 && rec[m][i + x][j + y] != INT_MAX && rec[c][i + x][j + y] > 0 && rec[c][i + x][j + y] != INT_MAX) {
                                            double ww = 1.0;

                                            if (!phase2 && (/*(double)(rec[m][i+x][j+y] - rec[m][i][j])/max[m]*(rec[m][i+x][j+y] - rec[m][i][j])/max[m] > 1.0/2 || */rec[c][i + x][j + y] < max[c] * 3 / 4)) {
                                                continue;
                                            }

                                            w[ix] += ww;
                                            ratio[ix] += ww * rec[m][i + x][j + y] / rec[c][i + x][j + y];
                                            count[ix] ++;
                                        }
                            }

                            ix++;
                        }

                        // compute new pixel values from the surrounding color ratios
                        double newc = 0.0;
                        int nc = 0;
                        ix = 0;

                        for (int m = 0; m < 3; m++) {
                            if (c == m) {
                                continue;
                            }

                            if (count[ix]) {
                                newc += (double)rec[m][i][j] / ratio[ix] * w[ix];
                                nc++;
                            }

                            ix++;
                        }

                        if (nc) {
                            rec[c][i][j] = - (int) (newc / nc);
                        }
                    }
                }
            }

        bool change = false;

        for (int i = 0; i < dh; i++)
            for (int j = 0; j < dw; j++)
                for (int c = 0; c < 3; c++) {
                    if (rec[c][i][j] < 0) {
                        rec[c][i][j] = -rec[c][i][j];
                    }

                    change = true;
                }

        if (!change) {
            break;
        }
    }

    printf ("Phase1 vege\n");

    // STEP II. recover fully lost pixels
    if (full) {
        int maxY = (299 * max[0] + 587 * max[1] + 114 * max[2]) / 1000;
        phase2 = false;

        for (int k = 0; k < 600; k++) {
            if (k > 200) {
                phase2 = true;
            }

            for (int i = 1; i < dh - 1; i++)
                for (int j = 1; j < dw - 1; j++) {
                    if (rec[0][i][j] == INT_MAX || rec[1][i][j] == INT_MAX || rec[2][i][j] == INT_MAX) {
                        int count = 0;
                        double yavg = 0, iavg = 0, qavg = 0, weight = 0.0;

                        for (int x = -1; x <= 1; x++)
                            for (int y = -1; y <= 1; y++)
                                if (rec[0][i + x][j + y] > 0 && rec[0][i + x][j + y] != INT_MAX && rec[1][i + x][j + y] > 0 && rec[1][i + x][j + y] != INT_MAX && rec[2][i + x][j + y] > 0 && rec[2][i + x][j + y] != INT_MAX) {
                                    // convert to yiq
                                    double Y = 0.299 * rec[0][i + x][j + y] + 0.587 * rec[1][i + x][j + y] + 0.114 * rec[2][i + x][j + y];
                                    double I = 0.596 * rec[0][i + x][j + y] - 0.275 * rec[1][i + x][j + y] - 0.321 * rec[2][i + x][j + y];
                                    double Q = 0.212 * rec[0][i + x][j + y] - 0.523 * rec[1][i + x][j + y] + 0.311 * rec[2][i + x][j + y];

                                    if (Y > maxY * 7 / 10) {
                                        double w = 1.0;// / (I*I+Q*Q);
                                        yavg += Y * w;
                                        iavg += I * w;
                                        qavg += Q * w;
                                        weight += w;
                                        count++;
                                    }
                                }

                        if ((!phase2 && count > 5) || (phase2 && count > 3)) {
                            double Y = yavg / weight;
                            double I = iavg / weight;
                            double Q = qavg / weight;
                            rec[0][i][j] = - (Y + 0.956 * I + 0.621 * Q);
                            rec[1][i][j] = - (Y - 0.272 * I - 0.647 * Q);
                            rec[2][i][j] = - (Y - 1.105 * I + 1.702 * Q);
                        }
                    }

                }

            bool change = false;

            for (int i = 0; i < dh; i++)
                for (int j = 0; j < dw; j++)
                    for (int c = 0; c < 3; c++) {
                        if (rec[c][i][j] < 0) {
                            rec[c][i][j] = -rec[c][i][j];
                        }

                        change = true;
                    }

            if (!change) {
                break;
            }
        }
    }

    int maxval = 0;

    for (int i = 0; i < dh; i++)
        for (int j = 0; j < dw; j++)
            for (int c = 0; c < 3; c++)
                if (rec[c][i][j] != INT_MAX && rec[c][i][j] > maxval) {
                    maxval = rec[c][i][j];
                }

    for (int i = 0; i < dh; i++)
        for (int j = 0; j < dw; j++)
            if (rec[0][i][j] == INT_MAX || rec[1][i][j] == INT_MAX || rec[2][i][j] == INT_MAX) {
                rec[0][i][j] = maxval;
                rec[1][i][j] = maxval;
                rec[2][i][j] = maxval;
            }

    if (hrmap[0] != NULL) {
        freeArray<float> (hrmap[0], dh);
        freeArray<float> (hrmap[1], dh);
        freeArray<float> (hrmap[2], dh);
    }

    hrmap[0] = allocArray<float> (dw, dh);
    hrmap[1] = allocArray<float> (dw, dh);
    hrmap[2] = allocArray<float> (dw, dh);

    this->full = full;

    for (int i = 0; i < dh; i++)
        for (int j = 0; j < dw; j++) {
            hrmap[0][i][j] = (double)rec[0][i][j] / ds->r(i, j);
            hrmap[1][i][j] = (double)rec[1][i][j] / ds->g(i, j);
            hrmap[2][i][j] = (double)rec[2][i][j] / ds->b(i, j);
        }

    /*    for (int i=0; i<dh; i++)
            for (int j=0; j<dh; j++) {
                ds->r(i,j) = CLIP (rec[0][i][j]);
                ds->g(i,j) = CLIP (rec[1][i][j]);
                ds->b(i,j) = CLIP (rec[2][i][j]);
            }
        ds->save ("test.png");
    */
    delete ds;
    freeArray<int> (rec[0], dh);
    freeArray<int> (rec[1], dh);
    freeArray<int> (rec[2], dh);

    printf ("HLMap vege\n");
}

void RawImageSource::hlRecovery (unsigned short* red, unsigned short* green, unsigned short* blue, int i, int sx1, int sx2, int skip)
{

    int blr = (i + SCALE / 2) / SCALE - 1;

    if (blr < 0 || blr >= H / SCALE - 1) {
        return;
    }

    double mr1 = 1.0 - ((double)((i + SCALE / 2) % SCALE) / SCALE + 0.5 / SCALE);
    int jx = 0;
    int maxcol = W / SCALE;

    for (int j = sx1, jx = 0; j < sx2; j += skip, jx++) {
        int blc = (j + SCALE / 2) / SCALE - 1;

        if (blc < 0 || blc >= maxcol - 1) {
            continue;
        }

        double mc1 = 1.0 - ((double)((j + SCALE / 2) % SCALE) / SCALE + 0.5 / SCALE);
        double mulr = mr1 * mc1 * hrmap[0][blr][blc] + mr1 * (1.0 - mc1) * hrmap[0][blr][blc + 1] + (1.0 - mr1) * mc1 * hrmap[0][blr + 1][blc] + (1.0 - mr1) * (1.0 - mc1) * hrmap[0][blr + 1][blc + 1];
        double mulg = mr1 * mc1 * hrmap[1][blr][blc] + mr1 * (1.0 - mc1) * hrmap[1][blr][blc + 1] + (1.0 - mr1) * mc1 * hrmap[1][blr + 1][blc] + (1.0 - mr1) * (1.0 - mc1) * hrmap[1][blr + 1][blc + 1];
        double mulb = mr1 * mc1 * hrmap[2][blr][blc] + mr1 * (1.0 - mc1) * hrmap[2][blr][blc + 1] + (1.0 - mr1) * mc1 * hrmap[2][blr + 1][blc] + (1.0 - mr1) * (1.0 - mc1) * hrmap[2][blr + 1][blc + 1];
        red[jx] = CLIP(red[jx] * mulr);
        green[jx] = CLIP(green[jx] * mulg);
        blue[jx] = CLIP(blue[jx] * mulb);
    }
}
}

