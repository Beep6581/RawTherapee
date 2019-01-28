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
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include "jaggedarray.h"
#include "rtengine.h"
#include "rawimagesource.h"
#include "rt_math.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "rt_algo.h"
#include "../rtgui/multilangmgr.h"

using namespace std;

namespace rtengine
{

void RawImageSource::dual_demosaic(bool isBayer, const RAWParams &raw, int winw, int winh, array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue, double &contrast, bool autoContrast)
{
    BENCHFUN
    const unsigned cfa4[2][2] = {{PFC(0,0), PFC(0,1)},{PFC(1,0),PFC(1,1)}};
    const unsigned cfa[2][2] = {{FC(0,0), FC(0,1)},{FC(1,0),FC(1,1)}};

    if (contrast == 0.0 && !autoContrast) {
        // contrast == 0.0 means only first demosaicer will be used
        if(isBayer) {
            std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                if (plistener)
                    plistener->setProgress(p);
                return false;
            };
            if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZEVNG4) ) {
                if (plistener) {
                    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZE)));
                }
                amaze_demosaic(winw, winh, 0, 0, winw, winh, rawData, red, green, blue, cfa, setProgCancel, initialGain, border, 65535.f, 65535.f);
            } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::DCBVNG4) ) {
                if (plistener) {
                    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::DCB)));
                }
                dcb_demosaic(winw, winh, rawData, red, green, blue, cfa, setProgCancel, raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
            } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::RCDVNG4) ) {
                if (plistener) {
                    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::RCD)));
                }
                rcd_demosaic(winw, winh, rawData, red, green, blue, cfa, setProgCancel);
            }
        } else {
            unsigned xtrans[6][6];
            ri->getXtransMatrix(xtrans);
            float rgb_cam[3][4];
            ri->getRgbCam(rgb_cam);
            if (plistener) {
                plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "Xtrans"));
            }
            std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                if (plistener)
                    plistener->setProgress(p);
                return false;
            };
            if (raw.xtranssensor.method == RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FOUR_PASS) ) {
                markesteijn_demosaic(winw, winh, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 3, true);
            } else {
                markesteijn_demosaic(winw, winh, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 1, false);
            }
        }

        return;
    }

    array2D<float> redTmp(winw, winh);
    array2D<float> greenTmp(winw, winh);
    array2D<float> blueTmp(winw, winh);
    array2D<float> L(winw, winh);

    if (isBayer) {
        if (plistener) {
            plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::VNG4)));
        }
        std::function<bool(double)> setProgCancel = [this](double p) -> bool {
            if(plistener)
                plistener->setProgress(p);
            return false;
        };
        vng4_demosaic(winw, winh, rawData, redTmp, greenTmp, blueTmp, cfa4, setProgCancel);

        if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZEVNG4) || raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::PIXELSHIFT)) {
            if (plistener) {
                plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZE)));
            }
            amaze_demosaic(winw, winh, 0, 0, winw, winh, rawData, red, green, blue, cfa, setProgCancel, initialGain, border, 65535.f, 65535.f);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::DCBVNG4) ) {
            if (plistener) {
                plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::DCB)));
            }
            dcb_demosaic(winw, winh, rawData, red, green, blue, cfa, setProgCancel, raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::RCDVNG4) ) {
            if (plistener) {
                plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::RCD)));
            }
            rcd_demosaic(winw, winh, rawData, red, green, blue, cfa, setProgCancel);
        }
    } else {
        if (plistener) {
            plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), raw.xtranssensor.method));
        }
        std::function<bool(double)> setProgCancel = [this](double p) -> bool {
            if (plistener)
                plistener->setProgress(p);
            return false;
        };
        unsigned xtrans[6][6];
        ri->getXtransMatrix(xtrans);
        float rgb_cam[3][4];
        ri->getRgbCam(rgb_cam);
        if (raw.xtranssensor.method == RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FOUR_PASS) ) {
            markesteijn_demosaic(winw, winh, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 3, true);
        } else {
            markesteijn_demosaic(winw, winh, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 1, false);
        }
        xtransfast_demosaic(winw, winh, rawData, redTmp, greenTmp, blueTmp, xtrans, setProgCancel);
    }

    const float xyz_rgb[3][3] = {          // XYZ from RGB
                                { 0.412453, 0.357580, 0.180423 },
                                { 0.212671, 0.715160, 0.072169 },
                                { 0.019334, 0.119193, 0.950227 }
                                };

    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < winh; ++i) {
            Color::RGB2L(red[i], green[i], blue[i], L[i], xyz_rgb, winw);
        }
    }
    // calculate contrast based blend factors to use vng4 in regions with low contrast
    JaggedArray<float> blend(winw, winh);
    float contrastf = contrast / 100.f;

    buildBlendMask(L, blend, winw, winh, contrastf, 1.f, autoContrast);
    contrast = contrastf * 100.f;

    // the following is split into 3 loops intentionally to avoid cache conflicts on CPUs with only 4-way cache
    #pragma omp parallel for
    for(int i = 0; i < winh; ++i) {
        for(int j = 0; j < winw; ++j) {
            red[i][j] = intp(blend[i][j], red[i][j], redTmp[i][j]);
        }
    }
    #pragma omp parallel for
    for(int i = 0; i < winh; ++i) {
        for(int j = 0; j < winw; ++j) {
            green[i][j] = intp(blend[i][j], green[i][j], greenTmp[i][j]);
        }
    }
    #pragma omp parallel for
    for(int i = 0; i < winh; ++i) {
        for(int j = 0; j < winw; ++j) {
            blue[i][j] = intp(blend[i][j], blue[i][j], blueTmp[i][j]);
        }
    }

}
}
