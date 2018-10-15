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
//  dual_demosaic.cc is free software: you can redistribute it and/or modify
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

#include "../jaggedarray.h"
#include "../rtengine.h"
#include "../rawimagesource.h"
#include "../rt_math.h"
#include "../StopWatch.h"
#include "../rt_algo.h"
#include "../../rtgui/multilangmgr.h"

using namespace std;

namespace rtengine
{

void RawImageSource::dual_demosaic(bool isBayer, const RAWParams &raw, int winw, int winh, array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue, double &contrast, bool autoContrast, int autoX, int autoY)
{
    BENCHFUN

    if (contrast == 0.0 && !autoContrast) {
        // contrast == 0.0 means only first demosaicer will be used
        if(isBayer) {
            if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZEVNG4) ) {
                if (plistener) {
                    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZE)));
                    plistener->setProgress(0.0);
                }
                std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                    if (plistener)
                        plistener->setProgress(p);
                    return false;
                };
                librtprocess::amaze_demosaic(winw, winh, 0, 0, winw, winh, rawData, red, green, blue, {{{FC(0,0), FC(0,1)},{FC(1,0),FC(1,1)}}}, setProgCancel, initialGain, border);
            } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::DCBVNG4) ) {
                dcb_demosaic(raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
            } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::RCDVNG4) ) {
                rcd_demosaic();
            }
        } else {
            if (raw.xtranssensor.method == RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FOUR_PASS) ) {
                if (plistener) {
                    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "Xtrans"));
                    plistener->setProgress(0.0);
                }
                std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                    if (plistener)
                        plistener->setProgress(p);
                    return false;
                };
                int xtrans[6][6];
                ri->getXtransMatrix(xtrans);
                float rgb_cam[3][4];
                ri->getRgbCam(rgb_cam);
                librtprocess::markesteijn_demosaic(W, H, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 3, true);
            } else {
                if (plistener) {
                    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "Xtrans"));
                    plistener->setProgress(0.0);
                }
                std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                    if (plistener)
                        plistener->setProgress(p);
                    return false;
                };
                int xtrans[6][6];
                ri->getXtransMatrix(xtrans);
                float rgb_cam[3][4];
                ri->getRgbCam(rgb_cam);
                librtprocess::markesteijn_demosaic(W, H, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 1, false);
            }
        }

        return;
    }

    array2D<float> redTmp(winw, winh);
    array2D<float> greenTmp(winw, winh);
    array2D<float> blueTmp(winw, winh);
    array2D<float> L(winw, winh);

    if (isBayer) {
        vng4_demosaic(rawData, redTmp, greenTmp, blueTmp);

        if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZEVNG4) || raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::PIXELSHIFT)) {
            if (plistener) {
                plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::AMAZE)));
                plistener->setProgress(0.0);
            }
            std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                if (plistener)
                    plistener->setProgress(p);
                return false;
            };
            librtprocess::amaze_demosaic(winw, winh, 0, 0, winw, winh, rawData, red, green, blue, {{{FC(0,0), FC(0,1)},{FC(1,0),FC(1,1)}}}, setProgCancel, initialGain, border);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::DCBVNG4) ) {
            dcb_demosaic(raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::RCDVNG4) ) {
            rcd_demosaic();
        }
    } else {
        if (raw.xtranssensor.method == RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FOUR_PASS) ) {
            if (plistener) {
                plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "Xtrans"));
                plistener->setProgress(0.0);
            }
            std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                if (plistener)
                    plistener->setProgress(p);
                return false;
            };
            int xtrans[6][6];
            ri->getXtransMatrix(xtrans);
            float rgb_cam[3][4];
            ri->getRgbCam(rgb_cam);
            librtprocess::markesteijn_demosaic(W, H, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 3, true);
        } else {
            if (plistener) {
                plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "Xtrans"));
                plistener->setProgress(0.0);
            }
            std::function<bool(double)> setProgCancel = [this](double p) -> bool {
                if (plistener)
                    plistener->setProgress(p);
                return false;
            };
            int xtrans[6][6];
            ri->getXtransMatrix(xtrans);
            float rgb_cam[3][4];
            ri->getRgbCam(rgb_cam);
            librtprocess::markesteijn_demosaic(W, H, rawData, red, green, blue, xtrans, rgb_cam, setProgCancel, 1, false);
        }
        int xtrans[6][6];
        ri->getXtransMatrix(xtrans);
        fast_xtransdemosaic(rawData, redTmp, greenTmp, blueTmp, xtrans);
    }

    const float xyz_rgb[3][3] = {          // XYZ from RGB
                                { 0.412453, 0.357580, 0.180423 },
                                { 0.212671, 0.715160, 0.072169 },
                                { 0.019334, 0.119193, 0.950227 }
                                };

    if (autoContrast && autoX >= 0 && autoY >= 0) {
        constexpr int rectSize = 40;
        const int autoWidth = min(rectSize, winw - autoX);
        const int autoHeight = min(rectSize, winh - autoY);
        if (std::min(autoWidth, autoHeight) > 20) {
            array2D<float> autoL(autoWidth, autoHeight);
            for(int i = 0; i < autoHeight; ++i) {
                Color::RGB2L(red[i + autoY] + autoX, green[i + autoY] + autoX, blue[i + autoY] + autoX, autoL[i], xyz_rgb, autoWidth);
            }
            // calculate contrast based blend factors to use vng4 in regions with low contrast
            JaggedArray<float> blend(autoWidth - 2, autoHeight - 2);
            int c = calcContrastThreshold(autoL, blend, autoWidth, autoHeight);
            if(c < 100) {
                contrast = c; // alternative : contrast = c - 1
            }
        }
    }

    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < winh; ++i) {
            Color::RGB2L(red[i], green[i], blue[i], L[i], xyz_rgb, winw);
        }
    }
    // calculate contrast based blend factors to use vng4 in regions with low contrast
    JaggedArray<float> blend(winw, winh);
    buildBlendMask(L, blend, winw, winh, contrast / 100.f);

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
