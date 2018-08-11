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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <cmath>

#include "rawimagesource.h"
#include "rawimagesource_i.h"
#include "jaggedarray.h"
#include "rawimage.h"
#include "iccmatrices.h"
#include "rt_math.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#define BENCHMARK
#include "StopWatch.h"

using namespace std;

namespace rtengine
{

#define DIST(a,b) (std::fabs(a-b))

extern const Settings* settings;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::eahd_demosaic ()
{
    BENCHFUN
    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::EAHD)));
        plistener->setProgress (0.0);
    }

    // prepare cache and constants for cielab conversion
    //TODO: revisit after conversion to D50 illuminant
    lc00 = (0.412453 * imatrices.rgb_cam[0][0] + 0.357580 * imatrices.rgb_cam[0][1] + 0.180423 * imatrices.rgb_cam[0][2]) ;// / 0.950456;
    lc01 = (0.412453 * imatrices.rgb_cam[1][0] + 0.357580 * imatrices.rgb_cam[1][1] + 0.180423 * imatrices.rgb_cam[1][2]) ;// / 0.950456;
    lc02 = (0.412453 * imatrices.rgb_cam[2][0] + 0.357580 * imatrices.rgb_cam[2][1] + 0.180423 * imatrices.rgb_cam[2][2]) ;// / 0.950456;

    lc10 = 0.212671 * imatrices.rgb_cam[0][0] + 0.715160 * imatrices.rgb_cam[0][1] + 0.072169 * imatrices.rgb_cam[0][2];
    lc11 = 0.212671 * imatrices.rgb_cam[1][0] + 0.715160 * imatrices.rgb_cam[1][1] + 0.072169 * imatrices.rgb_cam[1][2];
    lc12 = 0.212671 * imatrices.rgb_cam[2][0] + 0.715160 * imatrices.rgb_cam[2][1] + 0.072169 * imatrices.rgb_cam[2][2];

    lc20 = (0.019334 * imatrices.rgb_cam[0][0] + 0.119193 * imatrices.rgb_cam[0][1] + 0.950227 * imatrices.rgb_cam[0][2]) ;// / 1.088754;
    lc21 = (0.019334 * imatrices.rgb_cam[1][0] + 0.119193 * imatrices.rgb_cam[1][1] + 0.950227 * imatrices.rgb_cam[1][2]) ;// / 1.088754;
    lc22 = (0.019334 * imatrices.rgb_cam[2][0] + 0.119193 * imatrices.rgb_cam[2][1] + 0.950227 * imatrices.rgb_cam[2][2]) ;// / 1.088754;

    int maxindex = 3 * 65536; //2*65536 3 = avoid crash 3/2013 J.Desmis
    cache = new float[maxindex];
    threshold = 0.008856 * MAXVALD;

    for (int i = 0; i < maxindex; i++) {
        cache[i] = std::cbrt(double(i) / MAXVALD);
    }

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

    convert_to_cielab_row (rh[0], gh[0], bh[0], lLh[0], lah[0], lbh[0]);
    convert_to_cielab_row (rv[0], gv[0], bv[0], lLv[0], lav[0], lbv[0]);
    convert_to_cielab_row (rh[1], gh[1], bh[1], lLh[1], lah[1], lbh[1]);
    convert_to_cielab_row (rv[1], gv[1], bv[1], lLv[1], lav[1], lbv[1]);

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

        convert_to_cielab_row (rh[(i + 1) % 3], gh[(i + 1) % 4], bh[(i + 1) % 3], lLh[(i + 1) % 3], lah[(i + 1) % 3], lbh[(i + 1) % 3]);
        convert_to_cielab_row (rv[(i + 1) % 3], gv[(i + 1) % 4], bv[(i + 1) % 3], lLv[(i + 1) % 3], lav[(i + 1) % 3], lbv[(i + 1) % 3]);

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

            for (int dmi = 0; dmi < 9; dmi++) {
                wh += (dLmaph[dmi] <= eL) * (dCamaph[dmi] <= eCa) * (dCbmaph[dmi] <= eCb);
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

            for (int dmi = 0; dmi < 9; dmi++) {
                wv += (dLmapv[dmi] <= eL) * (dCamapv[dmi] <= eCa) * (dCbmapv[dmi] <= eCb);
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
                    green[i - 1][j] = gh[(i - 1) % 4][j];
                } else if (hc < vc) {
                    green[i - 1][j] = gv[(i - 1) % 4][j];
                } else {
                    green[i - 1][j] = (gh[(i - 1) % 4][j] + gv[(i - 1) % 4][j]) / 2;
                }
            }
        }

        if (!(i % 20) && plistener) {
            plistener->setProgress ((double)i / (H - 2));
        }
    }

    // finish H-2th and H-1th row, homogenity value is still valailable
    for (int i = H - 1; i < H + 1; i++)
        for (int j = 0; j < W; j++) {
            int hc = homh[(i - 1) % 3][j];
            int vc = homv[(i - 1) % 3][j];

            if (hc > vc) {
                green[i - 1][j] = gh[(i - 1) % 4][j];
            } else if (hc < vc) {
                green[i - 1][j] = gv[(i - 1) % 4][j];
            } else {
                green[i - 1][j] = (gh[(i - 1) % 4][j] + gv[(i - 1) % 4][j]) / 2;
            }
        }

    // Interpolate R and B
    #pragma omp parallel for
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