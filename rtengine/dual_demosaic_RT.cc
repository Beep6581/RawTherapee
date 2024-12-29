////////////////////////////////////////////////////////////////
//
//          combine demosaic algorithms
//
//
//  copyright (c) 2018  Ingo Weyrich <heckflosse67@gmx.de>
//
//  blends output of two demosaicers based on contrast
//
//
//  dual_demosaic_RT.cc is free software: you can redistribute it and/or modify
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
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include "color.h"
#include "jaggedarray.h"
#include "procparams.h"
#include "rawimagesource.h"
#include "rt_algo.h"
#include "rt_math.h"
#include "rtengine.h"

#include "rtgui/options.h"

using namespace std;

namespace rtengine
{

void RawImageSource::dual_demosaic_RT(bool isBayer, const procparams::RAWParams &raw, int winw, int winh, const array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue, double &contrast, bool autoContrast)
{

    if (contrast == 0.0 && !autoContrast) {
        // contrast == 0.0 means only first demosaicer will be used
        if(isBayer) {
            if (raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEBILINEAR) ||
                raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEVNG4)) {
                amaze_demosaic_RT(0, 0, winw, winh, rawData, red, green, blue, options.chunkSizeAMAZE, options.measure);
            } else if (raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBBILINEAR) ||
                       raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBVNG4)) {
                dcb_demosaic(raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
            } else if (raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDBILINEAR) ||
                       raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDVNG4)) {
                rcd_demosaic(options.chunkSizeRCD, options.measure);
            }
        } else {
            if (raw.xtranssensor.method == procparams::RAWParams::XTransSensor::getMethodString(procparams::RAWParams::XTransSensor::Method::FOUR_PASS)) {
                xtrans_interpolate (3, true, options.chunkSizeXT, options.measure);
            } else {
                xtrans_interpolate (1, false, options.chunkSizeXT, options.measure);
            }
        }

        return;
    }

    array2D<float> L(winw, winh);

    if (isBayer) {
        if (raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEBILINEAR) ||
            raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEVNG4) ||
            raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::PIXELSHIFT)) {
            amaze_demosaic_RT(0, 0, winw, winh, rawData, red, green, blue, options.chunkSizeAMAZE, options.measure);
        } else if (raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBBILINEAR) ||
                   raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBVNG4)) {
            dcb_demosaic(raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
        } else if (raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDBILINEAR) ||
                   raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDVNG4)) {
            rcd_demosaic(options.chunkSizeRCD, options.measure);
        }
    } else {
        if (raw.xtranssensor.method == procparams::RAWParams::XTransSensor::getMethodString(procparams::RAWParams::XTransSensor::Method::FOUR_PASS)) {
            xtrans_interpolate (3, true, options.chunkSizeXT, options.measure);
        } else {
            xtrans_interpolate (1, false, options.chunkSizeXT, options.measure);
        }
    }

    const float xyz_rgb[3][3] = {          // XYZ from RGB
                                { 0.412453, 0.357580, 0.180423 },
                                { 0.212671, 0.715160, 0.072169 },
                                { 0.019334, 0.119193, 0.950227 }
                                };

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for(int i = 0; i < winh; ++i) {
        Color::RGB2L(red[i], green[i], blue[i], L[i], xyz_rgb, winw);
    }

    // calculate contrast based blend factors to use flat demosaicer in regions with low contrast
    JaggedArray<float> blend(winw, winh);
    float contrastf = contrast / 100.0;

    buildBlendMask(L, blend, winw, winh, contrastf, autoContrast);
    contrast = contrastf * 100.f;

    if (isBayer) {
        if (raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEBILINEAR) ||
            raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDBILINEAR) ||
            raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBBILINEAR)) {
            bayer_bilinear_demosaic(blend, rawData, red, green, blue);
        } else {
            array2D<float>& redTmp = L; // L is not needed anymore => reuse it
            array2D<float> greenTmp(winw, winh);
            array2D<float> blueTmp(winw, winh);
            vng4_demosaic(rawData, redTmp, greenTmp, blueTmp);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif
            for(int i = 0; i < winh; ++i) {
                // the following is split into 3 loops intentionally to avoid cache conflicts on CPUs with only 4-way cache
                for(int j = 0; j < winw; ++j) {
                    red[i][j] = intp(blend[i][j], red[i][j], redTmp[i][j]);
                }
                for(int j = 0; j < winw; ++j) {
                    green[i][j] = intp(blend[i][j], green[i][j], greenTmp[i][j]);
                }
                for(int j = 0; j < winw; ++j) {
                    blue[i][j] = intp(blend[i][j], blue[i][j], blueTmp[i][j]);
                }
            }
        }
    } else {
        fast_xtrans_interpolate_blend(blend, rawData, red, green, blue);
    }
}
}
