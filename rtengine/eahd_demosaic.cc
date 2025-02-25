/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2018 Gabor Horvath <hgabor@rawtherapee.com> and other RawTherapee contributors
 *  Split out to own compilation unit and made speedup 2018 Ingo Weyrich (heckflosse67@gmx.de)
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

#include "color.h"
#include "rawimage.h"
#include "rawimagesource.h"
#include "rawimagesource_i.h"
#include "jaggedarray.h"
#include "iccmatrices.h"
#include "rt_math.h"
#include "rtgui/multilangmgr.h"
//#define BENCHMARK
#include "StopWatch.h"

using namespace std;

namespace rtengine
{


inline void RawImageSource::interpolate_row_g (float* agh, float* agv, int i)
{

    for (int j = 0; j < W; j++) {
        if (ri->ISGREEN(i, j)) {
            agh[j] = rawData[i][j];
            agv[j] = rawData[i][j];
        } else {
            int gh = 0;
            int gv = 0;

            if (j > 1 && j < W - 2) {
                gh = (-rawData[i][j - 2] + 2 * rawData[i][j - 1] + 2 * rawData[i][j] + 2 * rawData[i][j + 1] - rawData[i][j + 2]) / 4;
                int maxgh = max(rawData[i][j - 1], rawData[i][j + 1]);
                int mingh = min(rawData[i][j - 1], rawData[i][j + 1]);

                if (gh > maxgh) {
                    gh = maxgh;
                } else if (gh < mingh) {
                    gh = mingh;
                }
            } else if (j == 0) {
                gh = rawData[i][1];
            } else if (j == 1) {
                gh = (rawData[i][0] + rawData[i][2]) / 2;
            } else if (j == W - 1) {
                gh = rawData[i][W - 2];
            } else if (j == W - 2) {
                gh = (rawData[i][W - 1] + rawData[i][W - 3]) / 2;
            }

            if (i > 1 && i < H - 2) {
                gv = (-rawData[i - 2][j] + 2 * rawData[i - 1][j] + 2 * rawData[i][j] + 2 * rawData[i + 1][j] - rawData[i + 2][j]) / 4;
                int maxgv = max(rawData[i - 1][j], rawData[i + 1][j]);
                int mingv = min(rawData[i - 1][j], rawData[i + 1][j]);

                if (gv > maxgv) {
                    gv = maxgv;
                } else if (gv < mingv) {
                    gv = mingv;
                }
            } else if (i == 0) {
                gv = rawData[1][j];
            } else if (i == 1) {
                gv = (rawData[0][j] + rawData[2][j]) / 2;
            } else if (i == H - 1) {
                gv = rawData[H - 2][j];
            } else if (i == H - 2) {
                gv = (rawData[H - 1][j] + rawData[H - 3][j]) / 2;
            }

            agh[j] = gh;
            agv[j] = gv;
        }
    }
}

inline void RawImageSource::interpolate_row_rb (float* ar, float* ab, float* pg, float* cg, float* ng, int i)
{
    const auto getPg = [pg](int index) {
        return
            pg
                ? pg[index]
                : 0.f;
    };

    const auto getNg = [ng](int index) {
        return
            ng
                ? ng[index]
                : 0.f;
    };

    float *nonGreen1 = ar;
    float *nonGreen2 = ab;

    if ((ri->ISBLUE(i, 0) || ri->ISBLUE(i, 1))) {
        std::swap(nonGreen1, nonGreen2);
    }
    int j = FC(i, 0) & 1;
    if (j) {
        // linear R-G interp. horizontally
        float val1;

        val1 = cg[0] + rawData[i][1] - cg[1];

        nonGreen1[0] = CLIP(val1);
        // linear B-G interp. vertically
        float val2;

        if (i == 0) {
            val2 = getNg(0) + rawData[1][0] - cg[0];
        } else if (i == H - 1) {
            val2 = getPg(0) + rawData[H - 2][0] - cg[0];
        } else {
            val2 = cg[0] + (rawData[i - 1][0] - getPg(0) + rawData[i + 1][0] - getNg(0)) / 2;
        }

        nonGreen2[0] = val2;
    }
    // RGRGR or GRGRGR line
    for (; j < W - 1; j += 2) {
        // nonGreen of row is simple
        nonGreen1[j] = rawData[i][j];
        // non green of next row: cross interpolation
        float nonGreen = 0.f;
        int n = 0;

        if (i > 0) {
            if (j > 0) {
                nonGreen += rawData[i - 1][j - 1] - getPg(j - 1);
                n++;
            }
            nonGreen += rawData[i - 1][j + 1] - getPg(j + 1);
            n++;
        }

        if (i < H - 1) {
            if (j > 0) {
                nonGreen += rawData[i + 1][j - 1] - getNg(j - 1);
                n++;
            }
            nonGreen += rawData[i + 1][j + 1] - getNg(j + 1);
            n++;
        }

        nonGreen2[j] = cg[j] + nonGreen / std::max(n, 1);

        // linear R-G interp. horizontally
        float val1;

        if (j == W - 2) {
            val1 = cg[W - 1] + rawData[i][W - 2] - cg[W - 2];
        } else {
            val1 = cg[j + 1] + (rawData[i][j] - cg[j] + rawData[i][j + 2] - cg[j + 2]) / 2;
        }

        nonGreen1[j + 1] = CLIP(val1);
        // linear B-G interp. vertically
        float val2;

        if (i == 0) {
            val2 = getNg(j + 1) + rawData[1][j + 1] - cg[j + 1];
        } else if (i == H - 1) {
            val2 = getPg(j + 1) + rawData[H - 2][j + 1] - cg[j + 1];
        } else {
            val2 = cg[j + 1] + (rawData[i - 1][j + 1] - getPg(j + 1) + rawData[i + 1][j + 1] - getNg(j + 1)) / 2;
        }

        nonGreen2[j + 1] = val2;
    }

    if(j == W - 1) {
        // nonGreen of row is simple
        nonGreen1[j] = rawData[i][j];
        // non green of next row: cross interpolation
        float nonGreen = 0.f;
        int n = 0;

        if (i > 0) {
            nonGreen += rawData[i - 1][j - 1] - getPg(j - 1);
            n++;
        }

        if (i < H - 1) {
            nonGreen += rawData[i + 1][j - 1] - getNg(j - 1);
            n++;
        }

        nonGreen2[j] = cg[j] + nonGreen / std::max(n, 1);
    }
}

#define DIST(a,b) (std::fabs(a-b))

void RawImageSource::eahd_demosaic ()
{
    BENCHFUN
    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), M("TP_RAW_EAHD")));
        plistener->setProgress (0.0);
    }

    // prepare constants for cielab conversion
    //TODO: revisit after conversion to D50 illuminant
    const float lc00 = (0.412453 * imatrices.rgb_cam[0][0] + 0.357580 * imatrices.rgb_cam[0][1] + 0.180423 * imatrices.rgb_cam[0][2]) ;// / 0.950456;
    const float lc01 = (0.412453 * imatrices.rgb_cam[1][0] + 0.357580 * imatrices.rgb_cam[1][1] + 0.180423 * imatrices.rgb_cam[1][2]) ;// / 0.950456;
    const float lc02 = (0.412453 * imatrices.rgb_cam[2][0] + 0.357580 * imatrices.rgb_cam[2][1] + 0.180423 * imatrices.rgb_cam[2][2]) ;// / 0.950456;

    const float lc10 = 0.212671 * imatrices.rgb_cam[0][0] + 0.715160 * imatrices.rgb_cam[0][1] + 0.072169 * imatrices.rgb_cam[0][2];
    const float lc11 = 0.212671 * imatrices.rgb_cam[1][0] + 0.715160 * imatrices.rgb_cam[1][1] + 0.072169 * imatrices.rgb_cam[1][2];
    const float lc12 = 0.212671 * imatrices.rgb_cam[2][0] + 0.715160 * imatrices.rgb_cam[2][1] + 0.072169 * imatrices.rgb_cam[2][2];

    const float lc20 = (0.019334 * imatrices.rgb_cam[0][0] + 0.119193 * imatrices.rgb_cam[0][1] + 0.950227 * imatrices.rgb_cam[0][2]) ;// / 1.088754;
    const float lc21 = (0.019334 * imatrices.rgb_cam[1][0] + 0.119193 * imatrices.rgb_cam[1][1] + 0.950227 * imatrices.rgb_cam[1][2]) ;// / 1.088754;
    const float lc22 = (0.019334 * imatrices.rgb_cam[2][0] + 0.119193 * imatrices.rgb_cam[2][1] + 0.950227 * imatrices.rgb_cam[2][2]) ;// / 1.088754;

    const float wp[3][3] = {{lc00, lc01, lc02}, {lc10, lc11, lc12}, {lc20, lc21, lc22}};

    // end of cielab preparation

    JaggedArray<float>
    rh (W, 3), gh (W, 4), bh (W, 3),
    rv (W, 3), gv (W, 4), bv (W, 3),
    lLh (W, 3), lah (W, 3), lbh (W, 3),
    lLv (W, 3), lav (W, 3), lbv (W, 3);
    JaggedArray<uint16_t> homh (W, 3), homv (W, 3);

    // interpolate first two lines
    interpolate_row_g (gh[0], gv[0], 0);
    interpolate_row_g (gh[1], gv[1], 1);
    interpolate_row_g (gh[2], gv[2], 2);
    interpolate_row_rb (rh[0], bh[0], nullptr, gh[0], gh[1], 0);
    interpolate_row_rb (rv[0], bv[0], nullptr, gv[0], gv[1], 0);
    interpolate_row_rb (rh[1], bh[1], gh[0], gh[1], gh[2], 1);
    interpolate_row_rb (rv[1], bv[1], gv[0], gv[1], gv[2], 1);

    Color::RGB2Lab(rh[0], gh[0], bh[0], lLh[0], lah[0], lbh[0], wp, W);
    Color::RGB2Lab(rv[0], gv[0], bv[0], lLv[0], lav[0], lbv[0], wp, W);
    Color::RGB2Lab(rh[1], gh[1], bh[1], lLh[1], lah[1], lbh[1], wp, W);
    Color::RGB2Lab(rv[1], gv[1], bv[1], lLv[1], lav[1], lbv[1], wp, W);

    for (int j = 0; j < W; j++) {
        homh[0][j] = 0;
        homv[0][j] = 0;
        homh[1][j] = 0;
        homv[1][j] = 0;
    }

    float dLmaph[9];
    float dLmapv[9];
    float dCamaph[9];
    float dCamapv[9];
    float dCbmaph[9];
    float dCbmapv[9];

    for (int i = 1; i < H - 1; i++) {
        int mod[3] = {(i-1) % 3, i % 3, (i+1) %3};
        int ix = i % 3;
        int imx = (i - 1) % 3;
        int ipx = (i + 1) % 3;

        if (i < H - 2) {
            interpolate_row_g  (gh[(i + 2) % 4], gv[(i + 2) % 4], i + 2);
            interpolate_row_rb (rh[(i + 1) % 3], bh[(i + 1) % 3], gh[i % 4], gh[(i + 1) % 4], gh[(i + 2) % 4], i + 1);
            interpolate_row_rb (rv[(i + 1) % 3], bv[(i + 1) % 3], gv[i % 4], gv[(i + 1) % 4], gv[(i + 2) % 4], i + 1);
        } else {
            interpolate_row_rb (rh[(i + 1) % 3], bh[(i + 1) % 3], gh[i % 4], gh[(i + 1) % 4], nullptr, i + 1);
            interpolate_row_rb (rv[(i + 1) % 3], bv[(i + 1) % 3], gv[i % 4], gv[(i + 1) % 4], nullptr, i + 1);
        }

        Color::RGB2Lab(rh[(i + 1) % 3], gh[(i + 1) % 4], bh[(i + 1) % 3], lLh[(i + 1) % 3], lah[(i + 1) % 3], lbh[(i + 1) % 3], wp, W);
        Color::RGB2Lab(rv[(i + 1) % 3], gv[(i + 1) % 4], bv[(i + 1) % 3], lLv[(i + 1) % 3], lav[(i + 1) % 3], lbv[(i + 1) % 3], wp, W);

        for (int j = 0; j < W; j++) {
            homh[ipx][j] = 0;
            homv[ipx][j] = 0;
        }

        for (int j = 1; j < W - 1; j++) {
            int dmi = 0;
            for (int x = -1; x <= 0; x++) {
                int idx = mod[x + 1];

                for (int y = -1; y <= 1; y++) {
                    // compute distance in a, b, and L
                    if (dmi < 4) {
                        int sh = homh[idx][j + y];
                        int sv = homv[idx][j + y];

                        if (sh > sv) { // fixate horizontal pixel
                            dLmaph[dmi]  = DIST(lLh[ix][j], lLh[idx][j + y]);
                            dCamaph[dmi] = DIST(lah[ix][j], lah[idx][j + y]);
                            dCbmaph[dmi] = DIST(lbh[ix][j], lbh[idx][j + y]);
                            dLmapv[dmi]  = DIST(lLv[ix][j], lLh[idx][j + y]);
                            dCamapv[dmi] = DIST(lav[ix][j], lah[idx][j + y]);
                            dCbmapv[dmi] = DIST(lbv[ix][j], lbh[idx][j + y]);
                        } else if (sh < sv) {
                            dLmaph[dmi]  = DIST(lLh[ix][j], lLv[idx][j + y]);
                            dCamaph[dmi] = DIST(lah[ix][j], lav[idx][j + y]);
                            dCbmaph[dmi] = DIST(lbh[ix][j], lbv[idx][j + y]);
                            dLmapv[dmi]  = DIST(lLv[ix][j], lLv[idx][j + y]);
                            dCamapv[dmi] = DIST(lav[ix][j], lav[idx][j + y]);
                            dCbmapv[dmi] = DIST(lbv[ix][j], lbv[idx][j + y]);
                        } else {
                            dLmaph[dmi]  = DIST(lLh[ix][j], lLh[idx][j + y]);
                            dCamaph[dmi] = DIST(lah[ix][j], lah[idx][j + y]);
                            dCbmaph[dmi] = DIST(lbh[ix][j], lbh[idx][j + y]);
                            dLmapv[dmi]  = DIST(lLv[ix][j], lLv[idx][j + y]);
                            dCamapv[dmi] = DIST(lav[ix][j], lav[idx][j + y]);
                            dCbmapv[dmi] = DIST(lbv[ix][j], lbv[idx][j + y]);
                        }
                    } else {
                        dLmaph[dmi]  = DIST(lLh[ix][j], lLh[idx][j + y]);
                        dCamaph[dmi] = DIST(lah[ix][j], lah[idx][j + y]);
                        dCbmaph[dmi] = DIST(lbh[ix][j], lbh[idx][j + y]);
                        dLmapv[dmi]  = DIST(lLv[ix][j], lLv[idx][j + y]);
                        dCamapv[dmi] = DIST(lav[ix][j], lav[idx][j + y]);
                        dCbmapv[dmi] = DIST(lbv[ix][j], lbv[idx][j + y]);
                    }

                    dmi++;
                }
            }
            int idx = mod[2];

            for (int y = -1; y <= 1; y++) {
                // compute distance in a, b, and L
                dLmaph[dmi]  = DIST(lLh[ix][j], lLh[idx][j + y]);
                dCamaph[dmi] = DIST(lah[ix][j], lah[idx][j + y]);
                dCbmaph[dmi] = DIST(lbh[ix][j], lbh[idx][j + y]);
                dLmapv[dmi]  = DIST(lLv[ix][j], lLv[idx][j + y]);
                dCamapv[dmi] = DIST(lav[ix][j], lav[idx][j + y]);
                dCbmapv[dmi] = DIST(lbv[ix][j], lbv[idx][j + y]);
                dmi++;
            }

            // compute eL & eC
            float eL = min(max(dLmaph[3], dLmaph[5]), max(dLmapv[1], dLmapv[7]));
            float eCa = min(max(dCamaph[3], dCamaph[5]), max(dCamapv[1], dCamapv[7]));
            float eCb = min(max(dCbmaph[3], dCbmaph[5]), max(dCbmapv[1], dCbmapv[7]));

            int wh = 0;

            for (int d = 0; d < 9; ++d) {
                wh += (dLmaph[d] <= eL) * (dCamaph[d] <= eCa) * (dCbmaph[d] <= eCb);
            }

            homh[imx][j - 1] += wh;
            homh[imx][j]  += wh;
            homh[imx][j + 1] += wh;
            homh[ix][j - 1] += wh;
            homh[ix][j]   += wh;
            homh[ix][j + 1] += wh;
            homh[ipx][j - 1] += wh;
            homh[ipx][j]  += wh;
            homh[ipx][j + 1] += wh;

            int wv = 0;

            for (int d = 0; d < 9; ++d) {
                wv += (dLmapv[d] <= eL) * (dCamapv[d] <= eCa) * (dCbmapv[d] <= eCb);
            }

            homv[imx][j - 1] += wv;
            homv[imx][j]  += wv;
            homv[imx][j + 1] += wv;
            homv[ix][j - 1] += wv;
            homv[ix][j]   += wv;
            homv[ix][j + 1] += wv;
            homv[ipx][j - 1] += wv;
            homv[ipx][j]  += wv;
            homv[ipx][j + 1] += wv;
        }

        // finalize image
        for (int j = 0; j < W; j++) {
            if (ri->ISGREEN(i - 1, j)) {
                green[i - 1][j] = rawData[i - 1][j];
            } else {
                int hc = homh[imx][j];
                int vc = homv[imx][j];

                if (hc > vc) {
                    green[i - 1][j] = std::max(0.f, gh[(i - 1) % 4][j]);
                } else if (hc < vc) {
                    green[i - 1][j] = std::max(0.f, gv[(i - 1) % 4][j]);
                } else {
                    green[i - 1][j] = std::max(0.f, (gh[(i - 1) % 4][j] + gv[(i - 1) % 4][j]) / 2);
                }
            }
        }

        if (!(i % 20) && plistener) {
            plistener->setProgress ((double)i / (H - 2));
        }
    }

    // finish H-2th and H-1th row, homogeneity value is still available
    for (int i = H - 1; i < H + 1; i++)
        for (int j = 0; j < W; j++) {
            int hc = homh[(i - 1) % 3][j];
            int vc = homv[(i - 1) % 3][j];

            if (hc > vc) {
                green[i - 1][j] = std::max(0.f, gh[(i - 1) % 4][j]);
            } else if (hc < vc) {
                green[i - 1][j] = std::max(0.f, gv[(i - 1) % 4][j]);
            } else {
                green[i - 1][j] = std::max(0.f, (gh[(i - 1) % 4][j] + gv[(i - 1) % 4][j]) / 2);
            }
        }

    // Interpolate R and B
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < H; i++) {
        if (i == 0) {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], nullptr, green[i], green[i + 1], i, 1.0, 1.0, 1.0, 0, W, 1);
        } else if (i == H - 1) {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], green[i - 1], green[i], nullptr, i, 1.0, 1.0, 1.0, 0, W, 1);
        } else {
            interpolate_row_rb_mul_pp (rawData, red[i], blue[i], green[i - 1], green[i], green[i + 1], i, 1.0, 1.0, 1.0, 0, W, 1);
        }
    }
}

}
