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
//#include <glib/gstdio.h>
#include "procparams.h"
#include "rt_math.h"
#include "safegtk.h"
#include "safekeyfile.h"
#include "dcp.h"
#include "../rtgui/multilangmgr.h"
#include "../rtgui/version.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/paramsedited.h"
#include "../rtgui/options.h"
#include <locale.h>
#define APPVERSION VERSION

using namespace std;
extern Options options;

namespace rtengine {
namespace procparams {
	const int tr=(int) options.rtSettings.top_right;
	const int br=(int) options.rtSettings.bot_right;
	const int tl=(int) options.rtSettings.top_left;
	const int bl=(int) options.rtSettings.bot_left;
	
    const char *RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::numMethods]={"amaze","igv","lmmse","eahd", "hphd", "vng4", "dcb", "ahd", "fast", "mono", "none" };
    const char *RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::numMethods]={"3-pass (best)", "1-pass (medium)", "fast", "mono", "none" };

const char *RAWParams::ff_BlurTypestring[RAWParams::numFlatFileBlurTypes]={/*"Parametric",*/ "Area Flatfield", "Vertical Flatfield", "Horizontal Flatfield", "V+H Flatfield"};
std::vector<WBEntry*> WBParams::wbEntries;

bool ToneCurveParams::HLReconstructionNecessary(LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw) {
    if (options.rtSettings.verbose)
        printf("histRedRaw[  0]=%07d, histGreenRaw[  0]=%07d, histBlueRaw[  0]=%07d\nhistRedRaw[255]=%07d, histGreenRaw[255]=%07d, histBlueRaw[255]=%07d\n",
                histRedRaw[0], histGreenRaw[0], histBlueRaw[0], histRedRaw[255], histGreenRaw[255], histBlueRaw[255]);

    return histRedRaw[255]>50 || histGreenRaw[255]>50 || histBlueRaw[255]>50 || histRedRaw[0]>50 || histGreenRaw[0]>50 || histBlueRaw[0]>50;
}

void WBParams::init() {
    // Creation of the different methods and its associated temperature value
    wbEntries.push_back(new WBEntry("Camera"              ,WBT_CAMERA,      M("TP_WBALANCE_CAMERA"),        0, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Auto"                ,WBT_AUTO,        M("TP_WBALANCE_AUTO"),          0, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Daylight"            ,WBT_DAYLIGHT,    M("TP_WBALANCE_DAYLIGHT"),   5300, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Cloudy"              ,WBT_CLOUDY,      M("TP_WBALANCE_CLOUDY"),     6200, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Shade"               ,WBT_SHADE,       M("TP_WBALANCE_SHADE"),      7600, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Water 1"             ,WBT_WATER,       M("TP_WBALANCE_WATER1"),    35000, 0.3f,   1.1f));
    wbEntries.push_back(new WBEntry("Water 2"             ,WBT_WATER,       M("TP_WBALANCE_WATER2"),    48000, 0.63f, 1.38f));
    wbEntries.push_back(new WBEntry("Tungsten"            ,WBT_TUNGSTEN,    M("TP_WBALANCE_TUNGSTEN"),   2856, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F1"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO1"),      6430, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F2"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO2"),      4230, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F3"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO3"),      3450, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F4"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO4"),      2940, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F5"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO5"),      6350, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F6"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO6"),      4150, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F7"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO7"),      6500, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F8"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO8"),      5020, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F9"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO9"),      4330, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F10"            ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO10"),     5300, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F11"            ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO11"),     4000, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Fluo F12"            ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO12"),     3000, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("HMI Lamp"            ,WBT_LAMP,        M("TP_WBALANCE_HMI"),        4800, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("GTI Lamp"            ,WBT_LAMP,        M("TP_WBALANCE_GTI"),        5000, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("JudgeIII Lamp"       ,WBT_LAMP,        M("TP_WBALANCE_JUDGEIII"),   5100, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Solux Lamp 3500K"    ,WBT_LAMP,        M("TP_WBALANCE_SOLUX35"),    3480, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Solux Lamp 4100K"    ,WBT_LAMP,        M("TP_WBALANCE_SOLUX41"),    3930, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Solux Lamp 4700K"    ,WBT_LAMP,        M("TP_WBALANCE_SOLUX47"),    4700, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("NG Solux Lamp 4700K" ,WBT_LAMP,        M("TP_WBALANCE_SOLUX47_NG"), 4480, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("LED LSI Lumelex 2040",WBT_LED,         M("TP_WBALANCE_LED_LSI"),    2970, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("LED CRS SP12 WWMR16" ,WBT_LED,         M("TP_WBALANCE_LED_CRS"),    3050, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Flash 5500K"         ,WBT_FLASH,       M("TP_WBALANCE_FLASH55"),    5500, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Flash 6000K"         ,WBT_FLASH,       M("TP_WBALANCE_FLASH60"),    6000, 1.f,     1.f));
    wbEntries.push_back(new WBEntry("Flash 6500K"         ,WBT_FLASH,       M("TP_WBALANCE_FLASH65"),    6500, 1.f,     1.f));
    // Should remain the last one
    wbEntries.push_back(new WBEntry("Custom"              ,WBT_CUSTOM,      M("TP_WBALANCE_CUSTOM"),        0, 1.f,     1.f));
}

void WBParams::cleanup() {
    for (unsigned int i=0; i<wbEntries.size(); i++) {
        delete wbEntries[i];
    }
}

// Maps crop to resized width (e.g. smaller previews)
void CropParams::mapToResized(int resizedWidth, int resizedHeight, int scale, int &x1, int &x2, int &y1, int &y2) const {
    x1 = 0, x2 = resizedWidth, y1 = 0, y2 = resizedHeight;
    if (enabled) {
        x1 = min(resizedWidth-1,  max(0, x / scale));
        y1 = min(resizedHeight-1, max(0, y / scale));
        x2 = min(resizedWidth,    max(0, (x+w) / scale));
        y2 = min(resizedHeight,   max(0, (y+h) / scale));
    }
}

ColorToningParams::ColorToningParams () : hlColSat(60, 80, false), shadowsColSat(80, 208, false) {
    setDefaults();
}

void ColorToningParams::getDefaultColorCurve(std::vector<double> &curve) {
    double v[8]= { 0.050, 0.62, 0.25, 0.25,
                   0.585, 0.11, 0.25, 0.25 };

    curve.resize(9);
    curve.at(0) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}


void ColorToningParams::getDefaultOpacityCurve(std::vector<double> &curve) {
    double v[16]={ 0.00, 0.3, 0.35, 0.00,
                   0.25, 0.8, 0.35, 0.35,
                   0.70, 0.8, 0.35, 0.35,
                   1.00, 0.3, 0.00, 0.00 };
    curve.resize(17);
    curve.at(0 ) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}

void ColorToningParams::getDefaultCLCurve(std::vector<double> &curve) {
    double v[6]= { 0.00, 0.00,
                   0.35, 0.65,
                   1.00, 1.00 };

    curve.resize(7);
    curve.at(0) = double(DCT_NURBS);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}

void ColorToningParams::getDefaultCL2Curve(std::vector<double> &curve) {
    double v[6]= { 0.00, 0.00,
                   0.35, 0.65,
                   1.00, 1.00 };

    curve.resize(7);
    curve.at(0) = double(DCT_NURBS);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}

void ColorToningParams::setDefaults() {
    enabled = false;
    autosat=true;
    method = "Lab";

    getDefaultColorCurve(colorCurve);
    getDefaultOpacityCurve(opacityCurve);
    getDefaultCLCurve(clcurve);
    getDefaultCL2Curve(cl2curve);

    hlColSat.setValues(60, 80);
    shadowsColSat.setValues(80, 208);
    balance = 0;
    satProtectionThreshold = 30;
    saturatedOpacity = 80;
    strength = 50;
    lumamode = true;
    twocolor = "Std";
    redlow = 0.0;
    greenlow = 0.0;
    bluelow = 0.0;
    satlow = 0.0;
    sathigh = 0.0;
    redmed = 0.0;
    greenmed = 0.0;
    bluemed = 0.0;
    redhigh = 0.0;
    greenhigh = 0.0;
    bluehigh = 0.0;
}

void ColorToningParams::mixerToCurve(std::vector<double> &colorCurve, std::vector<double> &opacityCurve) const {
    // check if non null first
    if (!redlow && !greenlow && !bluelow && !redmed && !greenmed && !bluemed && !redhigh && !greenhigh && !bluehigh) {
        colorCurve.resize(1);
        colorCurve.at(0) = FCT_Linear;
        opacityCurve.resize(1);
        opacityCurve.at(0) = FCT_Linear;
        return;
    }

    float low[3];  // RGB color for shadows
    float med[3];  // RGB color for mid-tones
    float high[3]; // RGB color for highlights
    float lowSat  = 0.f;
    float medSat  = 0.f;
    float highSat = 0.f;
    float minTmp, maxTmp;

    // Fill the shadow mixer values of the Color TOning tool
    low[0] = float(redlow  )/100.f; // [-1. ; +1.]
    low[1] = float(greenlow)/100.f; // [-1. ; +1.]
    low[2] = float(bluelow )/100.f; // [-1. ; +1.]
    minTmp = min<float>(low[0], low[1], low[2]);
    maxTmp = max<float>(low[0], low[1], low[2]);
    if (maxTmp-minTmp > 0.005f) {
        float v[3];
        lowSat = (maxTmp-minTmp)/2.f;
        if      (low[0] == minTmp) v[0] = 0.f;
        else if (low[1] == minTmp) v[1] = 0.f;
        else if (low[2] == minTmp) v[2] = 0.f;
        if      (low[0] == maxTmp) v[0] = 1.f;
        else if (low[1] == maxTmp) v[1] = 1.f;
        else if (low[2] == maxTmp) v[2] = 1.f;
        if      (low[0] != minTmp && low[0] != maxTmp) v[0] = (low[0]-minTmp)/(maxTmp-minTmp);
        else if (low[1] != minTmp && low[1] != maxTmp) v[1] = (low[1]-minTmp)/(maxTmp-minTmp);
        else if (low[2] != minTmp && low[2] != maxTmp) v[2] = (low[2]-minTmp)/(maxTmp-minTmp);
        low[0] = v[0];
        low[1] = v[1];
        low[2] = v[2];
    }
    else {
        low[0] = low[1] = low[2] = 1.f;
    }

    // Fill the mid-tones mixer values of the Color TOning tool
    med[0] = float(redmed  )/100.f; // [-1. ; +1.]
    med[1] = float(greenmed)/100.f; // [-1. ; +1.]
    med[2] = float(bluemed )/100.f; // [-1. ; +1.]
    minTmp = min<float>(med[0], med[1], med[2]);
    maxTmp = max<float>(med[0], med[1], med[2]);
    if (maxTmp-minTmp > 0.005f) {
        float v[3];
        medSat = (maxTmp-minTmp)/2.f;
        if      (med[0] == minTmp) v[0] = 0.f;
        else if (med[1] == minTmp) v[1] = 0.f;
        else if (med[2] == minTmp) v[2] = 0.f;
        if      (med[0] == maxTmp) v[0] = 1.f;
        else if (med[1] == maxTmp) v[1] = 1.f;
        else if (med[2] == maxTmp) v[2] = 1.f;
        if      (med[0] != minTmp && med[0] != maxTmp) v[0] = (med[0]-minTmp)/(maxTmp-minTmp);
        else if (med[1] != minTmp && med[1] != maxTmp) v[1] = (med[1]-minTmp)/(maxTmp-minTmp);
        else if (med[2] != minTmp && med[2] != maxTmp) v[2] = (med[2]-minTmp)/(maxTmp-minTmp);
        med[0] = v[0];
        med[1] = v[1];
        med[2] = v[2];
    }
    else {
        med[0] = med[1] = med[2] = 1.f;
    }

    // Fill the highlight mixer values of the Color TOning tool
    high[0] = float(redhigh  )/100.f; // [-1. ; +1.]
    high[1] = float(greenhigh)/100.f; // [-1. ; +1.]
    high[2] = float(bluehigh )/100.f; // [-1. ; +1.]
    minTmp = min<float>(high[0], high[1], high[2]);
    maxTmp = max<float>(high[0], high[1], high[2]);
    if (maxTmp-minTmp > 0.005f) {
        float v[3];
        highSat = (maxTmp-minTmp)/2.f;
        if      (high[0] == minTmp) v[0] = 0.f;
        else if (high[1] == minTmp) v[1] = 0.f;
        else if (high[2] == minTmp) v[2] = 0.f;
        if      (high[0] == maxTmp) v[0] = 1.f;
        else if (high[1] == maxTmp) v[1] = 1.f;
        else if (high[2] == maxTmp) v[2] = 1.f;
        if      (high[0] != minTmp && high[0] != maxTmp) v[0] = (high[0]-minTmp)/(maxTmp-minTmp);
        else if (high[1] != minTmp && high[1] != maxTmp) v[1] = (high[1]-minTmp)/(maxTmp-minTmp);
        else if (high[2] != minTmp && high[2] != maxTmp) v[2] = (high[2]-minTmp)/(maxTmp-minTmp);
        high[0] = v[0];
        high[1] = v[1];
        high[2] = v[2];
    }
    else {
        high[0] = high[1] = high[2] = 1.f;
    }




    const double xPosLow  = 0.1;
    const double xPosMed  = 0.4;
    const double xPosHigh = 0.7;



    colorCurve.resize( medSat!=0.f ? 13 : 9 );
    colorCurve.at(0) = FCT_MinMaxCPoints;
    opacityCurve.resize(13);
    opacityCurve.at(0) = FCT_MinMaxCPoints;

    float h, s, l;
    int idx = 1;

    if (lowSat == 0.f) {
        if (medSat != 0.f)
            Color::rgb2hsl(med[0], med[1], med[2], h, s, l);
        else// highSat can't be null if the 2 other ones are!
            Color::rgb2hsl(high[0], high[1], high[2], h, s, l);
    }
    else
        Color::rgb2hsl(low[0], low[1], low[2], h, s, l);
    colorCurve.at(idx++)  = xPosLow;
    colorCurve.at(idx++)  = h;
    colorCurve.at(idx++)  = 0.35;
    colorCurve.at(idx++)  = 0.35;

    if (medSat != 0.f) {
        Color::rgb2hsl(med[0], med[1], med[2], h, s, l);
        colorCurve.at(idx++)  = xPosMed;
        colorCurve.at(idx++)  = h;
        colorCurve.at(idx++)  = 0.35;
        colorCurve.at(idx++)  = 0.35;
    }

    if (highSat == 0.f) {
        if (medSat != 0.f)
            Color::rgb2hsl(med[0], med[1], med[2], h, s, l);
        else// lowSat can't be null if the 2 other ones are!
            Color::rgb2hsl(low[0], low[1], low[2], h, s, l);
    }
    else
        Color::rgb2hsl(high[0], high[1], high[2], h, s, l);
    colorCurve.at(idx++)  = xPosHigh;
    colorCurve.at(idx++)  = h;
    colorCurve.at(idx++)  = 0.35;
    colorCurve.at(idx)    = 0.35;

    opacityCurve.at(1)  = xPosLow;
    opacityCurve.at(2)  = double(lowSat);
    opacityCurve.at(3)  = 0.35;
    opacityCurve.at(4)  = 0.35;
    opacityCurve.at(5)  = xPosMed;
    opacityCurve.at(6)  = double(medSat);
    opacityCurve.at(7)  = 0.35;
    opacityCurve.at(8)  = 0.35;
    opacityCurve.at(9)  = xPosHigh;
    opacityCurve.at(10) = double(highSat);
    opacityCurve.at(11) = 0.35;
    opacityCurve.at(12) = 0.35;
}

void ColorToningParams::slidersToCurve(std::vector<double> &colorCurve, std::vector<double> &opacityCurve) const {
    if (hlColSat.value[0]==0 && shadowsColSat.value[0]==0) { // if both opacity are null, set both curves to Linear
        colorCurve.resize(1);
        colorCurve.at(0) = FCT_Linear;
        opacityCurve.resize(1);
        opacityCurve.at(0) = FCT_Linear;
        return;
    }

    colorCurve.resize(9);
    colorCurve.at(0) = FCT_MinMaxCPoints;
    colorCurve.at(1) = 0.26 + 0.12*double(balance)/100.;
    colorCurve.at(2) = double(shadowsColSat.value[1])/360.;
    colorCurve.at(3) = 0.35;
    colorCurve.at(4) = 0.35;
    colorCurve.at(5) = 0.64 + 0.12*double(balance)/100.;
    colorCurve.at(6) = double(hlColSat.value[1])/360.;
    colorCurve.at(7) = 0.35;
    colorCurve.at(8) = 0.35;

    opacityCurve.resize(9);
    opacityCurve.at(0) = FCT_MinMaxCPoints;
    opacityCurve.at(1) = colorCurve.at(1);
    opacityCurve.at(2) = double(shadowsColSat.value[0])/100.;
    opacityCurve.at(3) = 0.35;
    opacityCurve.at(4) = 0.35;
    opacityCurve.at(5) = colorCurve.at(5);
    opacityCurve.at(6) = double(hlColSat.value[0])/100.;
    opacityCurve.at(7) = 0.35;
    opacityCurve.at(8) = 0.35;
}

void ColorToningParams::getCurves(ColorGradientCurve &colorCurveLUT, OpacityCurve &opacityCurveLUT, const double xyz_rgb[3][3], const double rgb_xyz[3][3], bool &opautili) const {
    float satur=0.8f;
    float lumin=0.5f;//middle of luminance for optimization of gamut - no real importance...as we work in XYZ and gamut control

    // Transform slider values to control points
    std::vector<double> cCurve, oCurve;
    if (method=="RGBSliders" || method=="Splitlr")
        slidersToCurve(cCurve, oCurve);
    else if (method=="Splitco")
        mixerToCurve(cCurve, oCurve);
    else {
        cCurve = this->colorCurve;
        oCurve = this->opacityCurve;
    }

    if(method=="Lab") {
        if(twocolor=="Separ") satur=0.9f;
        if(twocolor=="All"  || twocolor=="Two") satur=0.9f;
        colorCurveLUT.SetXYZ(cCurve, xyz_rgb, rgb_xyz, satur, lumin);
        opacityCurveLUT.Set(oCurve, opautili);		
    }
    else if(method=="Splitlr" || method=="Splitco") {
        colorCurveLUT.SetXYZ(cCurve, xyz_rgb, rgb_xyz, satur, lumin);
        opacityCurveLUT.Set(oCurve, opautili);
    }
    else if(method.substr(0,3)=="RGB") {
        colorCurveLUT.SetRGB(cCurve, xyz_rgb, rgb_xyz);
        opacityCurveLUT.Set(oCurve, opautili);
     }
}
//WaveletParams::WaveletParams (): hueskin(-5, 25, 170, 120, false), hueskin2(-260, -250, -130, -140, false), hllev(50, 75, 100, 98, false), bllev(0, 2, 50, 25, false), pastlev(0, 2, 30, 20, false), satlev(30, 45, 130, 100, false), edgcont(0, 20, 100, 75, false){

WaveletParams::WaveletParams (): hueskin(-5, 25, 170, 120, false), hueskin2(-260, -250, -130, -140, false), hllev(50, 75, 100, 98, false), bllev(0, 2, 50, 25, false), 
	pastlev(0, 2, 30, 20, false), satlev(30, 45, 130, 100, false),  edgcont(bl, tl, br, tr, false), /* edgcont(0, 10, 75, 40, false),*/level0noise(0, 0, false),level1noise(0, 0, false), level2noise(0, 0, false){
    setDefaults (); 
}

void WaveletParams::getDefaultOpacityCurveRG(std::vector<double> &curve) {
    double v[8]= { 0.0, 0.50,0.35,0.35,
                   1.00, 0.50,0.35,0.35};			   

    curve.resize(9);
    curve.at(0) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}
void WaveletParams::getDefaultOpacityCurveBY(std::vector<double> &curve) {
    double v[8]= { 0.0, 0.50,0.35,0.35,
                   1.00, 0.50,0.35,0.35};			   
				   
    curve.resize(9);
    curve.at(0 ) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}


void WaveletParams::getDefaultOpacityCurveW(std::vector<double> &curve) {
 double v[16]={ 0.00, 0.35, 0.35, 0.00,
                   0.35, 0.75, 0.35, 0.35,
                   0.60, 0.75, 0.35, 0.35,
                   1.00, 0.35, 0.00, 0.00 };
    curve.resize(17);
    curve.at(0) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}

void WaveletParams::getDefaultOpacityCurveWL(std::vector<double> &curve) {
	    double v[8]= { 0.0, 0.50,0.35,0.35,
                   1.00, 0.50,0.35,0.35};			   

 /*double v[12]={ 0.00, 0.53, 0.35, 0.00,
                   0.42, 0.53, 0.35, 0.35,
                   1.00, 0.15, 0.00, 0.00 };
				   */
    curve.resize(9);
    curve.at(0) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}


void WaveletParams::getDefaultCCWCurve(std::vector<double> &curve) {
    double v[12]= { 0.0, 0.25, 0.35, 0.35,
                   0.50, 0.75, 0.35, 0.35, 0.90, 0.0, 0.35, 0.35};	
			   
    curve.resize(13);
    curve.at(0 ) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
	
}

void WaveletParams::getCurves(WavCurve &cCurve, WavOpacityCurveRG &opacityCurveLUTRG, WavOpacityCurveBY &opacityCurveLUTBY, WavOpacityCurveW &opacityCurveLUTW, WavOpacityCurveWL &opacityCurveLUTWL) const {
		cCurve.Set(this->ccwcurve);
        opacityCurveLUTRG.Set(this->opacityCurveRG);
		opacityCurveLUTBY.Set(this->opacityCurveBY);
		opacityCurveLUTW.Set(this->opacityCurveW);		
		opacityCurveLUTWL.Set(this->opacityCurveWL);		
		
}

void WaveletParams::setDefaults() {
	getDefaultCCWCurve(ccwcurve);
    getDefaultOpacityCurveRG(opacityCurveRG);
    getDefaultOpacityCurveBY(opacityCurveBY);	
    getDefaultOpacityCurveW(opacityCurveW);	
    getDefaultOpacityCurveWL(opacityCurveWL);	
    enabled = false;  
    median = false;  
    medianlev = false;  
    linkedg = true;  
    cbenab = false;  
    lipst = false;  
    Medgreinf = "less"; //"none";
    avoid = false;
    tmr = false;
    strength = 100;
    balance = 0;
    iter = 0;
    wavclCurve.clear ();
    wavclCurve.push_back(DCT_Linear);

	Lmethod      	 = "4_";
	CHmethod      	 = "without";
	CHSLmethod      	 = "SL";
	EDmethod      	 = "CU";
	BAmethod      	 = "none";
	TMmethod      	 = "none";
	HSmethod      	 = "with";
	CLmethod      	 = "all";
	Backmethod      	 = "grey";
	Dirmethod      	 = "all";
	Tilesmethod      	 = "full";
	daubcoeffmethod      	 = "4_";
	rescon      = 0;
	resconH      = 0;
	reschro      = 0;
	tmrs      = 0;
	gamma      = 1;
	sky      	 = 0.;
	sup      	 = 0;
	thres      	 = 7;
	chroma       = 5;
	chro      	 = 0;
	contrast      	 = 0;
	edgrad 		 =15;
	edgval		 = 0;
	edgthresh    = 10;
	thr      	 = 35;
	thrH      	 = 65;
    skinprotect = 0.;
	hueskin.setValues(-5, 25, 170, 120); 
	hueskin2.setValues(-260, -250, -130, -140); 
    threshold=5;
    threshold2=4;
    edgedetect=80;
    edgedetectthr=20;
    edgedetectthr2=0;
	hllev.setValues(50, 75, 100, 98);     
	bllev.setValues(0, 2, 50, 25);     
	pastlev.setValues(0, 2, 30, 20);     
	satlev.setValues(30, 45, 130, 100);     
//	edgcont.setValues(bl, tl, br, tr);     
	edgcont.setValues(0, 10, 75, 40);     
    level0noise.setValues(0, 0);
    level1noise.setValues(0, 0);
    level2noise.setValues(0, 0);
    hhcurve.clear ();
    hhcurve.push_back(FCT_Linear);
    Chcurve.clear ();
    Chcurve.push_back(FCT_Linear);
	
    for(int i = 0; i < 9; i ++)
    {
        c[i] = 0;
    }
    for(int i = 0; i < 9; i ++)
    {
        ch[i] = 0;
    }

}


DirPyrDenoiseParams::DirPyrDenoiseParams () {
    setDefaults ();
}

void DirPyrDenoiseParams::getDefaultNoisCurve(std::vector<double> &curve) {
    double v[8]={ 0.05, 0.15, 0.35, 0.35,
                  0.55, 0.04, 0.35, 0.35};
    curve.resize(9);
    curve.at(0 ) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}

void DirPyrDenoiseParams::getDefaultCCCurve(std::vector<double> &curve) {
  //  double v[8]= { 0.15, 0.00,0.35,0.35,
  //                 0.60, 0.05,0.35,0.35};
    double v[8]= { 0.05, 0.50,0.35,0.35,
                   0.35, 0.05,0.35,0.35};

    curve.resize(9);
    curve.at(0 ) = double(FCT_MinMaxCPoints);
    for (size_t i=1; i<curve.size(); ++i)
        curve.at(i) = v[i-1];
}


void DirPyrDenoiseParams::setDefaults() {

    getDefaultNoisCurve(lcurve);
    getDefaultCCCurve(cccurve);

    enabled = false;
    enhance = false;
    median = false;
    autochroma = false;
    luma = 0;
    passes = 1;
    dmethod = "Lab";
    Lmethod = "SLI";//"CUR";// SLIDER method with value 0 is set as default, while the default Lcurve is populated via getDefaultNoisCurve and can be switched to by the user
    Cmethod = "MAN";
    C2method = "AUTO";
    smethod = "shal";
    medmethod = "soft";
    methodmed = "none";
    rgbmethod = "soft";
    Ldetail = 0;
    chroma = 15;
    redchro = 0;
    bluechro = 0;
    gamma = 1.7;
}

void DirPyrDenoiseParams::getCurves(NoiseCurve &lCurve, NoiseCurve &cCurve) const {
    lCurve.Set(this->lcurve);
    cCurve.Set(this->cccurve);
}

void ToneCurveParams::setDefaults() {
    autoexp       = false;
    clip          = 0.02;
    expcomp       = 0;
    brightness    = 0;
    contrast      = 0;
    saturation    = 0;
    black         = 0;
    hlcompr       = 0;
    hlcomprthresh = 33;
    shcompr       = 50;
    curve.clear ();
    curve.push_back(DCT_Linear);
    curve2.clear ();
    curve2.push_back(DCT_Linear);
    curveMode     = ToneCurveParams::TC_MODE_STD;
    curveMode2    = ToneCurveParams::TC_MODE_STD;
    hrenabled = false;
    method  = "Blend";
}

void LensProfParams::setDefaults() {
    lcpFile="";
    useDist=useVign=true;
    useCA=false;
}

void CoarseTransformParams::setDefaults() {
    rotate = 0;
    hflip = false;
    vflip = false;
}

void RAWParams::setDefaults() {
    bayersensor.method = RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::amaze];
    bayersensor.ccSteps = 0;
    bayersensor.dcb_iterations = 2;
    bayersensor.dcb_enhance = true;
    //bayersensor.all_enhance = false;
    bayersensor.lmmse_iterations = 2;
    bayersensor.black0 = 0.0;
    bayersensor.black1 = 0.0;
    bayersensor.black2 = 0.0;
    bayersensor.black3 = 0.0;
    bayersensor.twogreen = true;
    bayersensor.linenoise = 0;
    bayersensor.greenthresh = 0;

    xtranssensor.method = RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::threePass];
    xtranssensor.ccSteps = 0;
    xtranssensor.blackred = 0.0;
    xtranssensor.blackgreen = 0.0;
    xtranssensor.blackblue = 0.0;

    expos=1.0;
    preser=0.0;
    df_autoselect = false;
    ff_AutoSelect = false;
    ff_BlurRadius = 32;
    ff_BlurType = RAWParams::ff_BlurTypestring[RAWParams::area_ff];
    ff_AutoClipControl = false;
    ff_clipControl = 0;
    cared = 0;
    cablue = 0;
    ca_autocorrect = false;
    hotPixelFilter = false;
    deadPixelFilter = false;
    hotdeadpix_thresh = 100;
}

void ColorManagementParams::setDefaults() {
    input   = "(cameraICC)";
    blendCMSMatrix = false;
    toneCurve = true;
    applyLookTable = true;
    applyBaselineExposureOffset = true;
    applyHueSatMap = true;
    dcpIlluminant = 0;
    working = "ProPhoto";
    output  = "RT_sRGB";
    gamma  = "default";
    gampos =2.22;
    slpos=4.5;
    freegamma = false;
}

ProcParams::ProcParams () { 

    setDefaults (); 
}

void ProcParams::init () {

    WBParams::init();
}

void ProcParams::cleanup () {

    WBParams::cleanup();
}

ProcParams* ProcParams::create () {

    return new ProcParams();
}

void ProcParams::destroy (ProcParams* pp) {

    delete pp;
}

void ProcParams::setDefaults () {

    toneCurve.setDefaults();

    labCurve.brightness      = 0;
    labCurve.contrast        = 0;
    labCurve.chromaticity    = 0;
    labCurve.avoidcolorshift = false;
    labCurve.lcredsk = true;
    labCurve.rstprotection   = 0;
    labCurve.lcurve.clear ();
    labCurve.lcurve.push_back(DCT_Linear);
    labCurve.acurve.clear ();
    labCurve.acurve.push_back(DCT_Linear);
    labCurve.bcurve.clear ();
    labCurve.bcurve.push_back(DCT_Linear);
    labCurve.cccurve.clear ();
    labCurve.cccurve.push_back(DCT_Linear);
    labCurve.chcurve.clear ();
    labCurve.chcurve.push_back(FCT_Linear);
    labCurve.lhcurve.clear ();
    labCurve.lhcurve.push_back(FCT_Linear);
    labCurve.hhcurve.clear ();
    labCurve.hhcurve.push_back(FCT_Linear);

    labCurve.lccurve.clear ();
    labCurve.lccurve.push_back(DCT_Linear);
    labCurve.clcurve.clear ();
    labCurve.clcurve.push_back(DCT_Linear);

    rgbCurves.lumamode          = false;
    rgbCurves.rcurve.clear ();
    rgbCurves.rcurve.push_back(DCT_Linear);
    rgbCurves.gcurve.clear ();
    rgbCurves.gcurve.push_back(DCT_Linear);
    rgbCurves.bcurve.clear ();
    rgbCurves.bcurve.push_back(DCT_Linear);

    colorToning.setDefaults();

    sharpenEdge.enabled         = false;
    sharpenEdge.passes          = 2;
    sharpenEdge.amount          = 50.0;
    sharpenEdge.threechannels   = false;

    sharpenMicro.enabled        = false;
    sharpenMicro.amount         = 20.0;
    sharpenMicro.uniformity     = 50.0;
    sharpenMicro.matrix         = false;

    sharpening.enabled          = false;
    sharpening.radius           = 0.5;
    sharpening.amount           = 200;
    sharpening.threshold.setValues(20, 80, 2000, 1200);
    sharpening.edgesonly        = false;
    sharpening.edges_radius     = 1.9;
    sharpening.edges_tolerance  = 1800;
    sharpening.halocontrol      = false;
    sharpening.halocontrol_amount = 85;
    sharpening.method           = "usm";
    sharpening.deconvradius     = 0.75;
    sharpening.deconviter       = 30;
    sharpening.deconvdamping    = 20;
    sharpening.deconvamount     = 75;

    prsharpening.enabled          = false;
    prsharpening.radius           = 0.5;
    prsharpening.amount           = 200;
    prsharpening.threshold.setValues(20, 80, 2000, 1200);
    prsharpening.edgesonly        = false;
    prsharpening.edges_radius     = 1.9;
    prsharpening.edges_tolerance  = 1800;
    prsharpening.halocontrol      = false;
    prsharpening.halocontrol_amount = 85;
    prsharpening.method           = "usm";
    prsharpening.deconvradius     = 0.5;
    prsharpening.deconviter       = 100;
    prsharpening.deconvdamping    = 0;
    prsharpening.deconvamount     = 100;

    vibrance.enabled            = false;
    vibrance.pastels            = 0;
    vibrance.saturated          = 0;
    vibrance.psthreshold.setValues(0, 75);
    vibrance.protectskins       = false;
    vibrance.avoidcolorshift    = true;
    vibrance.pastsattog     	= true;
    vibrance.skintonescurve.clear ();
    vibrance.skintonescurve.push_back(DCT_Linear);

    wb.method       = "Camera";
    wb.temperature  = 6504;
    wb.green        = 1.0;
    wb.equal        = 1.0;
    colorappearance.enabled       = false;
    colorappearance.degree        = 90;
    colorappearance.autodegree    = true;
    colorappearance.surround      = "Average";
    colorappearance.adaplum       = 16;
    colorappearance.badpixsl       = 0;
    colorappearance.adapscen      = 2000.0;
    colorappearance.autoadapscen    = true;
    colorappearance.algo          = "No";
    colorappearance.wbmodel       = "RawT";
    colorappearance.jlight        = 0.0;
    colorappearance.qbright       = 0.0;
    colorappearance.chroma        = 0.0;
    colorappearance.schroma       = 0.0;
    colorappearance.mchroma       = 0.0;
    colorappearance.rstprotection = 0.0;
    colorappearance.contrast      = 0.0;
    colorappearance.qcontrast     = 0.0;
    colorappearance.colorh        = 0.0;
    colorappearance.surrsource    = false;
    colorappearance.gamut         = true;
//    colorappearance.badpix        = false;
    colorappearance.datacie       = false;
    colorappearance.tonecie       = false;
 //   colorappearance.sharpcie      = false;
    colorappearance.curve.clear ();
    colorappearance.curve.push_back(DCT_Linear);
    colorappearance.curve2.clear ();
    colorappearance.curve2.push_back(DCT_Linear);
    colorappearance.curveMode     =ColorAppearanceParams::TC_MODE_LIGHT;
    colorappearance.curveMode2    = ColorAppearanceParams::TC_MODE_LIGHT;
    colorappearance.curve3.clear ();
    colorappearance.curve3.push_back(DCT_Linear);
    colorappearance.curveMode3    = ColorAppearanceParams::TC_MODE_CHROMA;

    impulseDenoise.enabled      = false;
    impulseDenoise.thresh       = 50;

    defringe.enabled            = false;
    defringe.radius             = 2.0;
    defringe.threshold          = 13;
    defringe.huecurve.resize (25);
    defringe.huecurve.at(0)     = FCT_MinMaxCPoints;
    defringe.huecurve.at(1)     = 0.166666667;
    defringe.huecurve.at(2)     = 0.;
    defringe.huecurve.at(3)     = 0.35;
    defringe.huecurve.at(4)     = 0.35;
    defringe.huecurve.at(5)     = 0.347;
    defringe.huecurve.at(6)     = 0.;
    defringe.huecurve.at(7)     = 0.35;
    defringe.huecurve.at(8)     = 0.35;
    defringe.huecurve.at(9)     = 0.513667426;
    defringe.huecurve.at(10)    = 0;
    defringe.huecurve.at(11)    = 0.35;
    defringe.huecurve.at(12)    = 0.35;
    defringe.huecurve.at(13)    = 0.668944571;
    defringe.huecurve.at(14)    = 0.;
    defringe.huecurve.at(15)    = 0.35;
    defringe.huecurve.at(16)    = 0.35;
    defringe.huecurve.at(17)    = 0.8287775246;
    defringe.huecurve.at(18)    = 0.97835991;
    defringe.huecurve.at(19)    = 0.35;
    defringe.huecurve.at(20)    = 0.35;
    defringe.huecurve.at(21)    = 0.9908883827;
    defringe.huecurve.at(22)    = 0.;
    defringe.huecurve.at(23)    = 0.35;
    defringe.huecurve.at(24)    = 0.35;

    dirpyrDenoise.setDefaults();

    epd.enabled = false;
    epd.strength = 0.5;
    epd.gamma = 1.0;
    epd.edgeStopping = 1.4;
    epd.scale = 0.3;
    epd.reweightingIterates = 0;

    sh.enabled       = false;
    sh.hq            = false;
    sh.highlights    = 0;
    sh.htonalwidth   = 80;
    sh.shadows       = 0;
    sh.stonalwidth   = 80;
    sh.localcontrast = 0;
    sh.radius        = 40;
    
    crop.enabled    = false;
    crop.x          = -1;
    crop.y          = -1;
    crop.w          = 15000;
    crop.h          = 15000;
    crop.fixratio   = false;
    crop.ratio      = "3:2";
    crop.orientation= "As Image";
    crop.guide      = "Rule of thirds";
    
    coarse.setDefaults();
    
    commonTrans.autofill = true;

    rotate.degree       = 0;

    distortion.amount     = 0;
    
    perspective.horizontal = 0;
    perspective.vertical   = 0;

    gradient.enabled = false;
    gradient.degree = 0;
    gradient.feather = 25;
    gradient.strength = 0.60;
    gradient.centerX = 0;
    gradient.centerY = 0;

    pcvignette.enabled = false;
    pcvignette.strength = 0.60;
    pcvignette.feather = 50;
    pcvignette.roundness = 50;

    cacorrection.red  = 0;
    cacorrection.blue = 0;
    

    vignetting.amount = 0;
    vignetting.radius = 50;
    vignetting.strength = 1;
    vignetting.centerX = 0;
    vignetting.centerY = 0;

    lensProf.setDefaults();

    chmixer.red[0] = 100;
    chmixer.red[1] = 0;
    chmixer.red[2] = 0;
    chmixer.green[0] = 0;
    chmixer.green[1] = 100;
    chmixer.green[2] = 0;
    chmixer.blue[0] = 0;
    chmixer.blue[1] = 0;
    chmixer.blue[2] = 100;

    blackwhite.autoc  = false;	
    blackwhite.enabledcc  = true;
    blackwhite.enabled  = false;
    blackwhite.mixerRed  = 33;
    blackwhite.mixerGreen  = 33;
    blackwhite.mixerBlue  = 33;
    blackwhite.mixerOrange  = 33;
    blackwhite.mixerYellow  = 33;
    blackwhite.mixerCyan  = 33;
    blackwhite.mixerMagenta  = 33;
    blackwhite.mixerPurple  = 33;
    blackwhite.gammaRed  = 0;
    blackwhite.gammaGreen  = 0;
    blackwhite.gammaBlue  = 0;
    blackwhite.luminanceCurve.clear ();
    blackwhite.luminanceCurve.push_back (FCT_Linear);
    blackwhite.method = "Desaturation";
    blackwhite.filter = "None";
    blackwhite.setting = "NormalContrast";
    blackwhite.beforeCurve.clear ();
    blackwhite.beforeCurve.push_back(DCT_Linear);
    blackwhite.beforeCurveMode     = BlackWhiteParams::TC_MODE_STD_BW;
    blackwhite.afterCurve.clear ();
    blackwhite.afterCurve.push_back(DCT_Linear);
    blackwhite.afterCurveMode     = BlackWhiteParams::TC_MODE_STD_BW;
    blackwhite.algo          = "SP";

    resize.enabled = false;
    resize.scale  = 1.0;
    resize.appliesTo = "Cropped area";
    resize.method = "Lanczos";
    resize.dataspec = 3;
    resize.width = 900;
    resize.height = 900;

    icm.setDefaults();
  
    dirpyrequalizer.enabled = false;
    dirpyrequalizer.gamutlab = false;
    for(int i = 0; i < 6; i ++)
    {
        dirpyrequalizer.mult[i] = 1.0;
    }
    dirpyrequalizer.threshold = 0.2;
    dirpyrequalizer.skinprotect = 0.;
    dirpyrequalizer.hueskin.setValues(-5, 25, 170, 120);        //default (b_l 0, t_l 30, b_r 170, t_r 120);
 //   dirpyrequalizer.algo = "FI";

    hsvequalizer.hcurve.clear ();
    hsvequalizer.hcurve.push_back (FCT_Linear);
    hsvequalizer.scurve.clear ();
    hsvequalizer.scurve.push_back (FCT_Linear);
    hsvequalizer.vcurve.clear ();
    hsvequalizer.vcurve.push_back (FCT_Linear);

    filmSimulation.setDefaults();

    raw.setDefaults();

    exif.clear ();
    iptc.clear ();

    rank = 0;
    colorlabel = 0;
    inTrash = false;
    
    ppVersion = PPVERSION;
}

static Glib::ustring expandRelativePath(Glib::ustring procparams_fname, Glib::ustring prefix, Glib::ustring embedded_fname) {
    if (embedded_fname == "" || !Glib::path_is_absolute(procparams_fname)) {
        return embedded_fname;
    }
    if (prefix != "") {
        if (embedded_fname.length() < prefix.length() || embedded_fname.substr(0, prefix.length()) != prefix) {
            return embedded_fname;
        }
        embedded_fname = embedded_fname.substr(prefix.length());
    }
    if (Glib::path_is_absolute(embedded_fname)) {
        return prefix + embedded_fname;
    }
    Glib::ustring absPath = prefix + Glib::path_get_dirname(procparams_fname) + G_DIR_SEPARATOR_S + embedded_fname;
    return absPath;
}

static Glib::ustring relativePathIfInside(Glib::ustring procparams_fname, bool fnameAbsolute, Glib::ustring embedded_fname) {
    if (fnameAbsolute || embedded_fname == "" || !Glib::path_is_absolute(procparams_fname)) {
        return embedded_fname;
    }
    Glib::ustring prefix = "";
    if (embedded_fname.length() > 5 && embedded_fname.substr(0, 5) == "file:") {
        embedded_fname = embedded_fname.substr(5);
        prefix = "file:";
    }
    if (!Glib::path_is_absolute(embedded_fname)) {
        return prefix + embedded_fname;
    }
    Glib::ustring dir1 = Glib::path_get_dirname(procparams_fname) + G_DIR_SEPARATOR_S;
    Glib::ustring dir2 = Glib::path_get_dirname(embedded_fname) + G_DIR_SEPARATOR_S;
    if (dir2.substr(0, dir1.length()) != dir1) {
        // it's in a different directory, ie not inside
        return prefix + embedded_fname;
    }
    return prefix + embedded_fname.substr(dir1.length());
}

int ProcParams::save (Glib::ustring fname, Glib::ustring fname2, bool fnameAbsolute, ParamsEdited* pedited) {

    if (!fname.length() && !fname2.length())
        return 0;

    SafeKeyFile keyFile;

    keyFile.set_string  ("Version", "AppVersion", APPVERSION);
    keyFile.set_integer ("Version", "Version",    PPVERSION);

    if (!pedited || pedited->general.rank)       keyFile.set_integer ("General", "Rank",       rank);
    if (!pedited || pedited->general.colorlabel) keyFile.set_integer ("General", "ColorLabel", colorlabel);
    if (!pedited || pedited->general.intrash)    keyFile.set_boolean ("General", "InTrash",    inTrash);

    // save tone curve
    if (!pedited || pedited->toneCurve.autoexp)    keyFile.set_boolean ("Exposure", "Auto",           toneCurve.autoexp);
    if (!pedited || pedited->toneCurve.clip)       keyFile.set_double  ("Exposure", "Clip",           toneCurve.clip);
    if (!pedited || pedited->toneCurve.expcomp)    keyFile.set_double  ("Exposure", "Compensation",   toneCurve.expcomp);
    if (!pedited || pedited->toneCurve.brightness) keyFile.set_integer ("Exposure", "Brightness",     toneCurve.brightness);
    if (!pedited || pedited->toneCurve.contrast)   keyFile.set_integer ("Exposure", "Contrast",       toneCurve.contrast);
    if (!pedited || pedited->toneCurve.saturation) keyFile.set_integer ("Exposure", "Saturation",     toneCurve.saturation);
    if (!pedited || pedited->toneCurve.black)      keyFile.set_integer ("Exposure", "Black",          toneCurve.black);
    if (!pedited || pedited->toneCurve.hlcompr)    keyFile.set_integer ("Exposure", "HighlightCompr", toneCurve.hlcompr);
    if (!pedited || pedited->toneCurve.hlcomprthresh) keyFile.set_integer ("Exposure", "HighlightComprThreshold", toneCurve.hlcomprthresh);
    if (!pedited || pedited->toneCurve.shcompr)       keyFile.set_integer ("Exposure", "ShadowCompr",             toneCurve.shcompr);
    // save highlight recovery settings
    if (!pedited || pedited->toneCurve.hrenabled)     keyFile.set_boolean ("HLRecovery", "Enabled",  toneCurve.hrenabled);
    if (!pedited || pedited->toneCurve.method)      keyFile.set_string  ("HLRecovery", "Method",   toneCurve.method);
    if (!pedited || pedited->toneCurve.curveMode)  {
        Glib::ustring method;
        switch (toneCurve.curveMode) {
        case (ToneCurveParams::TC_MODE_STD):
            method = "Standard";
            break;
        case (ToneCurveParams::TC_MODE_FILMLIKE):
            method = "FilmLike";
            break;
        case (ToneCurveParams::TC_MODE_SATANDVALBLENDING):
            method = "SatAndValueBlending";
            break;
        case (ToneCurveParams::TC_MODE_WEIGHTEDSTD):
            method = "WeightedStd";
            break;
        case (ToneCurveParams::TC_MODE_LUMINANCE):
            method = "Luminance";
            break;
        }
        keyFile.set_string  ("Exposure", "CurveMode", method);
    }
    if (!pedited || pedited->toneCurve.curveMode2)  {
        Glib::ustring method;
        switch (toneCurve.curveMode2) {
        case (ToneCurveParams::TC_MODE_STD):
            method = "Standard";
            break;
        case (ToneCurveParams::TC_MODE_FILMLIKE):
            method = "FilmLike";
            break;
        case (ToneCurveParams::TC_MODE_SATANDVALBLENDING):
            method = "SatAndValueBlending";
            break;
        case (ToneCurveParams::TC_MODE_WEIGHTEDSTD):
            method = "WeightedStd";
            break;
        case (ToneCurveParams::TC_MODE_LUMINANCE):
            method = "Luminance";
            break;
        }
        keyFile.set_string  ("Exposure", "CurveMode2", method);
    }
    if (!pedited || pedited->toneCurve.curve) {
        Glib::ArrayHandle<double> tcurve = toneCurve.curve;
        keyFile.set_double_list("Exposure", "Curve", tcurve);
    }
    if (!pedited || pedited->toneCurve.curve2) {
        Glib::ArrayHandle<double> tcurve = toneCurve.curve2;
        keyFile.set_double_list("Exposure", "Curve2", tcurve);
    }

    // save channel mixer
    if (!pedited || pedited->chmixer.red[0] || pedited->chmixer.red[1] || pedited->chmixer.red[2]) {
        Glib::ArrayHandle<int> rmix (chmixer.red, 3, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Channel Mixer", "Red",   rmix);
    }
    if (!pedited || pedited->chmixer.green[0] || pedited->chmixer.green[1] || pedited->chmixer.green[2]) {
        Glib::ArrayHandle<int> gmix (chmixer.green, 3, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Channel Mixer", "Green", gmix);
    }
    if (!pedited || pedited->chmixer.blue[0] || pedited->chmixer.blue[1] || pedited->chmixer.blue[2]) {
        Glib::ArrayHandle<int> bmix (chmixer.blue, 3, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Channel Mixer", "Blue",  bmix);
    }

    //save Black & White
    if (!pedited || pedited->blackwhite.enabled)      keyFile.set_boolean ("Black & White", "Enabled",             blackwhite.enabled);
    if (!pedited || pedited->blackwhite.method)       keyFile.set_string  ("Black & White", "Method",              blackwhite.method );
    if (!pedited || pedited->blackwhite.autoc)        keyFile.set_boolean ("Black & White", "Auto",                blackwhite.autoc);
    if (!pedited || pedited->blackwhite.enabledcc)    keyFile.set_boolean ("Black & White", "ComplementaryColors", blackwhite.enabledcc);
    if (!pedited || pedited->blackwhite.setting)      keyFile.set_string  ("Black & White", "Setting",             blackwhite.setting );
    if (!pedited || pedited->blackwhite.filter)       keyFile.set_string  ("Black & White", "Filter",              blackwhite.filter );
    if (!pedited || pedited->blackwhite.mixerRed)     keyFile.set_integer ("Black & White", "MixerRed",            blackwhite.mixerRed);
    if (!pedited || pedited->blackwhite.mixerOrange)  keyFile.set_integer ("Black & White", "MixerOrange",         blackwhite.mixerOrange);
    if (!pedited || pedited->blackwhite.mixerYellow)  keyFile.set_integer ("Black & White", "MixerYellow",         blackwhite.mixerYellow);
    if (!pedited || pedited->blackwhite.mixerGreen)   keyFile.set_integer ("Black & White", "MixerGreen",          blackwhite.mixerGreen);
    if (!pedited || pedited->blackwhite.mixerCyan)    keyFile.set_integer ("Black & White", "MixerCyan",           blackwhite.mixerCyan);
    if (!pedited || pedited->blackwhite.mixerBlue)    keyFile.set_integer ("Black & White", "MixerBlue",           blackwhite.mixerBlue);
    if (!pedited || pedited->blackwhite.mixerMagenta) keyFile.set_integer ("Black & White", "MixerMagenta",        blackwhite.mixerMagenta);
    if (!pedited || pedited->blackwhite.mixerPurple)  keyFile.set_integer ("Black & White", "MixerPurple",         blackwhite.mixerPurple);
    if (!pedited || pedited->blackwhite.gammaRed)     keyFile.set_integer ("Black & White", "GammaRed",            blackwhite.gammaRed);
    if (!pedited || pedited->blackwhite.gammaGreen)   keyFile.set_integer ("Black & White", "GammaGreen",          blackwhite.gammaGreen);
    if (!pedited || pedited->blackwhite.gammaBlue)    keyFile.set_integer ("Black & White", "GammaBlue",           blackwhite.gammaBlue);
    if (!pedited || pedited->blackwhite.algo)         keyFile.set_string  ("Black & White", "Algorithm",     	   blackwhite.algo);

    if (!pedited || pedited->blackwhite.luminanceCurve) {
        Glib::ArrayHandle<double> luminanceCurve = blackwhite.luminanceCurve;
        keyFile.set_double_list("Black & White", "LuminanceCurve", luminanceCurve);
    }
    if (!pedited || pedited->blackwhite.beforeCurveMode) {
        Glib::ustring mode;
        switch (blackwhite.beforeCurveMode) {
        case (BlackWhiteParams::TC_MODE_STD_BW):
            mode  = "Standard";
            break;
        case (BlackWhiteParams::TC_MODE_FILMLIKE_BW):
            mode  = "FilmLike";
            break;
        case (BlackWhiteParams::TC_MODE_SATANDVALBLENDING_BW):
            mode = "SatAndValueBlending";
            break;
        case (BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW):
            mode = "WeightedStd";
            break;
        }
        keyFile.set_string  ("Black & White", "BeforeCurveMode", mode);
    }

    if (!pedited || pedited->blackwhite.afterCurveMode) {
        Glib::ustring mode;
        switch (blackwhite.afterCurveMode) {
        case (BlackWhiteParams::TC_MODE_STD_BW):
            mode = "Standard";
            break;
        case (BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW):
            mode = "WeightedStd";
            break;
        default:
            break;
        }
        keyFile.set_string  ("Black & White", "AfterCurveMode", mode);
    }

    if (!pedited || pedited->blackwhite.beforeCurve) {
        Glib::ArrayHandle<double> tcurvebw = blackwhite.beforeCurve;
        keyFile.set_double_list("Black & White", "BeforeCurve", tcurvebw);
    }
    if (!pedited || pedited->blackwhite.afterCurve) {
        Glib::ArrayHandle<double> tcurvebw = blackwhite.afterCurve;
        keyFile.set_double_list("Black & White", "AfterCurve", tcurvebw);
    }

    // save luma curve
    if (!pedited || pedited->labCurve.brightness)      keyFile.set_integer ("Luminance Curve", "Brightness",                 labCurve.brightness);
    if (!pedited || pedited->labCurve.contrast)        keyFile.set_integer ("Luminance Curve", "Contrast",                   labCurve.contrast);
    if (!pedited || pedited->labCurve.chromaticity)    keyFile.set_integer ("Luminance Curve", "Chromaticity",               labCurve.chromaticity);
    if (!pedited || pedited->labCurve.avoidcolorshift) keyFile.set_boolean ("Luminance Curve", "AvoidColorShift",            labCurve.avoidcolorshift);
    if (!pedited || pedited->labCurve.rstprotection)   keyFile.set_double  ("Luminance Curve", "RedAndSkinTonesProtection",  labCurve.rstprotection);
    if (!pedited || pedited->labCurve.lcredsk)         keyFile.set_boolean ("Luminance Curve", "LCredsk",                    labCurve.lcredsk);
	
    if (!pedited || pedited->labCurve.lcurve)  {
        Glib::ArrayHandle<double> lcurve = labCurve.lcurve;
        keyFile.set_double_list("Luminance Curve", "LCurve", lcurve);
    }
    if (!pedited || pedited->labCurve.acurve)  {
        Glib::ArrayHandle<double> acurve = labCurve.acurve;
        keyFile.set_double_list("Luminance Curve", "aCurve", acurve);
    }
    if (!pedited || pedited->labCurve.bcurve)  {
        Glib::ArrayHandle<double> bcurve = labCurve.bcurve;
        keyFile.set_double_list("Luminance Curve", "bCurve", bcurve);
    }
    if (!pedited || pedited->labCurve.cccurve)  {
        Glib::ArrayHandle<double> cccurve = labCurve.cccurve;
        keyFile.set_double_list("Luminance Curve", "ccCurve", cccurve);
    }
    if (!pedited || pedited->labCurve.chcurve)  {
        Glib::ArrayHandle<double> chcurve = labCurve.chcurve;
        keyFile.set_double_list("Luminance Curve", "chCurve", chcurve);
    }
    if (!pedited || pedited->labCurve.lhcurve)  {
        Glib::ArrayHandle<double> lhcurve = labCurve.lhcurve;
        keyFile.set_double_list("Luminance Curve", "lhCurve", lhcurve);
    }
    if (!pedited || pedited->labCurve.hhcurve)  {
        Glib::ArrayHandle<double> hhcurve = labCurve.hhcurve;
        keyFile.set_double_list("Luminance Curve", "hhCurve", hhcurve);
    }

    if (!pedited || pedited->labCurve.lccurve)  {
        Glib::ArrayHandle<double> lccurve = labCurve.lccurve;
        keyFile.set_double_list("Luminance Curve", "LcCurve", lccurve);
    }
    if (!pedited || pedited->labCurve.clcurve)  {
        Glib::ArrayHandle<double> clcurve = labCurve.clcurve;
        keyFile.set_double_list("Luminance Curve", "ClCurve", clcurve);
    }

    // save sharpening
    if (!pedited || pedited->sharpening.enabled)            keyFile.set_boolean ("Sharpening", "Enabled",             sharpening.enabled);
    if (!pedited || pedited->sharpening.method)             keyFile.set_string  ("Sharpening", "Method",              sharpening.method);
    if (!pedited || pedited->sharpening.radius)             keyFile.set_double  ("Sharpening", "Radius",              sharpening.radius);
    if (!pedited || pedited->sharpening.amount)             keyFile.set_integer ("Sharpening", "Amount",              sharpening.amount);
    if (!pedited || pedited->sharpening.threshold) {
        Glib::ArrayHandle<int> thresh (sharpening.threshold.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Sharpening",   "Threshold", thresh);
    }
    if (!pedited || pedited->sharpening.edgesonly)          keyFile.set_boolean ("Sharpening", "OnlyEdges",           sharpening.edgesonly);
    if (!pedited || pedited->sharpening.edges_radius)       keyFile.set_double  ("Sharpening", "EdgedetectionRadius", sharpening.edges_radius);
    if (!pedited || pedited->sharpening.edges_tolerance)    keyFile.set_integer ("Sharpening", "EdgeTolerance",       sharpening.edges_tolerance);
    if (!pedited || pedited->sharpening.halocontrol)        keyFile.set_boolean ("Sharpening", "HalocontrolEnabled",  sharpening.halocontrol);
    if (!pedited || pedited->sharpening.halocontrol_amount) keyFile.set_integer ("Sharpening", "HalocontrolAmount",   sharpening.halocontrol_amount);
    if (!pedited || pedited->sharpening.deconvradius)       keyFile.set_double  ("Sharpening", "DeconvRadius",        sharpening.deconvradius);
    if (!pedited || pedited->sharpening.deconvamount)       keyFile.set_integer ("Sharpening", "DeconvAmount",        sharpening.deconvamount);
    if (!pedited || pedited->sharpening.deconvdamping)      keyFile.set_integer ("Sharpening", "DeconvDamping",       sharpening.deconvdamping);
    if (!pedited || pedited->sharpening.deconviter)         keyFile.set_integer ("Sharpening", "DeconvIterations",    sharpening.deconviter);

    // save vibrance
    if (!pedited || pedited->vibrance.enabled)          keyFile.set_boolean ("Vibrance", "Enabled",         vibrance.enabled);
    if (!pedited || pedited->vibrance.pastels)          keyFile.set_integer ("Vibrance", "Pastels",         vibrance.pastels);
    if (!pedited || pedited->vibrance.saturated)        keyFile.set_integer ("Vibrance", "Saturated",       vibrance.saturated);
    if (!pedited || pedited->vibrance.psthreshold) {
        Glib::ArrayHandle<int> thresh (vibrance.psthreshold.value, 2, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Vibrance", "PSThreshold", thresh);
    }
    if (!pedited || pedited->vibrance.protectskins)     keyFile.set_boolean ("Vibrance", "ProtectSkins",    vibrance.protectskins);
    if (!pedited || pedited->vibrance.avoidcolorshift)  keyFile.set_boolean ("Vibrance", "AvoidColorShift", vibrance.avoidcolorshift);
    if (!pedited || pedited->vibrance.pastsattog)       keyFile.set_boolean ("Vibrance", "PastSatTog",      vibrance.pastsattog);
    if (!pedited || pedited->vibrance.skintonescurve)  {
        Glib::ArrayHandle<double> skintonescurve = vibrance.skintonescurve;
        keyFile.set_double_list("Vibrance", "SkinTonesCurve", skintonescurve);
    }

    //save edge sharpening
    if (!pedited || pedited->sharpenEdge.enabled)       keyFile.set_boolean ("SharpenEdge", "Enabled",       sharpenEdge.enabled);
    if (!pedited || pedited->sharpenEdge.passes)        keyFile.set_integer ("SharpenEdge", "Passes",        sharpenEdge.passes);
    if (!pedited || pedited->sharpenEdge.amount)        keyFile.set_double  ("SharpenEdge", "Strength",      sharpenEdge.amount);
    if (!pedited || pedited->sharpenEdge.threechannels) keyFile.set_boolean ("SharpenEdge", "ThreeChannels", sharpenEdge.threechannels);

    //save micro-contrast sharpening
    if (!pedited || pedited->sharpenMicro.enabled)      keyFile.set_boolean ("SharpenMicro", "Enabled",    sharpenMicro.enabled);
    if (!pedited || pedited->sharpenMicro.matrix)       keyFile.set_boolean ("SharpenMicro", "Matrix",     sharpenMicro.matrix);
    if (!pedited || pedited->sharpenMicro.amount)       keyFile.set_double  ("SharpenMicro", "Strength",   sharpenMicro.amount);
    if (!pedited || pedited->sharpenMicro.uniformity)   keyFile.set_double  ("SharpenMicro", "Uniformity", sharpenMicro.uniformity);

/*
    // save colorBoost
    if (!pedited || pedited->colorBoost.amount)                   keyFile.set_integer ("Color Boost", "Amount",             colorBoost.amount);
    if (!pedited || pedited->colorBoost.avoidclip)                keyFile.set_boolean ("Color Boost", "AvoidColorClipping", colorBoost.avoidclip);
    if (!pedited || pedited->colorBoost.enable_saturationlimiter) keyFile.set_boolean ("Color Boost", "SaturationLimiter",  colorBoost.enable_saturationlimiter);
    if (!pedited || pedited->colorBoost.saturationlimit)          keyFile.set_double  ("Color Boost", "SaturationLimit",    colorBoost.saturationlimit);
*/

    // save wb
    if (!pedited || pedited->wb.method)      keyFile.set_string  ("White Balance", "Setting",     wb.method);
    if (!pedited || pedited->wb.temperature) keyFile.set_integer ("White Balance", "Temperature", wb.temperature);
    if (!pedited || pedited->wb.green)       keyFile.set_double  ("White Balance", "Green",       wb.green);
    if (!pedited || pedited->wb.equal)       keyFile.set_double  ("White Balance", "Equal",       wb.equal);

/*
    // save colorShift
    if (!pedited || pedited->colorShift.a)   keyFile.set_double ("Color Shift", "ChannelA", colorShift.a);
    if (!pedited || pedited->colorShift.b)   keyFile.set_double ("Color Shift", "ChannelB", colorShift.b);
*/
    // save colorappearance
    if (!pedited || pedited->colorappearance.enabled)       keyFile.set_boolean ("Color appearance", "Enabled",       colorappearance.enabled);
    if (!pedited || pedited->colorappearance.degree)        keyFile.set_integer ("Color appearance", "Degree",        colorappearance.degree);
    if (!pedited || pedited->colorappearance.autodegree)    keyFile.set_boolean ("Color appearance", "AutoDegree",    colorappearance.autodegree);
    if (!pedited || pedited->colorappearance.surround)      keyFile.set_string  ("Color appearance", "Surround",      colorappearance.surround);
 // if (!pedited || pedited->colorappearance.backgrd)       keyFile.set_integer ("Color appearance", "Background",    colorappearance.backgrd);
    if (!pedited || pedited->colorappearance.adaplum)       keyFile.set_double  ("Color appearance", "AdaptLum",      colorappearance.adaplum);
    if (!pedited || pedited->colorappearance.badpixsl)      keyFile.set_integer ("Color appearance", "Badpixsl",      colorappearance.badpixsl);
    if (!pedited || pedited->colorappearance.wbmodel)       keyFile.set_string  ("Color appearance", "Model",         colorappearance.wbmodel);
    if (!pedited || pedited->colorappearance.algo)          keyFile.set_string  ("Color appearance", "Algorithm",     colorappearance.algo);

    if (!pedited || pedited->colorappearance.jlight)        keyFile.set_double  ("Color appearance", "J-Light",       colorappearance.jlight);
    if (!pedited || pedited->colorappearance.qbright)       keyFile.set_double  ("Color appearance", "Q-Bright",      colorappearance.qbright);
    if (!pedited || pedited->colorappearance.chroma)        keyFile.set_double  ("Color appearance", "C-Chroma",      colorappearance.chroma);
    if (!pedited || pedited->colorappearance.schroma)       keyFile.set_double  ("Color appearance", "S-Chroma",      colorappearance.schroma);
    if (!pedited || pedited->colorappearance.mchroma)       keyFile.set_double  ("Color appearance", "M-Chroma",      colorappearance.mchroma);
    if (!pedited || pedited->colorappearance.contrast)      keyFile.set_double  ("Color appearance", "J-Contrast",    colorappearance.contrast);
    if (!pedited || pedited->colorappearance.qcontrast)     keyFile.set_double  ("Color appearance", "Q-Contrast",    colorappearance.qcontrast);
    if (!pedited || pedited->colorappearance.colorh)        keyFile.set_double  ("Color appearance", "H-Hue",         colorappearance.colorh);
    if (!pedited || pedited->colorappearance.rstprotection) keyFile.set_double  ("Color appearance", "RSTProtection", colorappearance.rstprotection);

    if (!pedited || pedited->colorappearance.adapscen)      keyFile.set_double  ("Color appearance", "AdaptScene",    colorappearance.adapscen);
    if (!pedited || pedited->colorappearance.autoadapscen)  keyFile.set_boolean ("Color appearance", "AutoAdapscen",  colorappearance.autoadapscen);
    if (!pedited || pedited->colorappearance.surrsource)    keyFile.set_boolean ("Color appearance", "SurrSource",    colorappearance.surrsource);
    if (!pedited || pedited->colorappearance.gamut)         keyFile.set_boolean ("Color appearance", "Gamut",         colorappearance.gamut);
//    if (!pedited || pedited->colorappearance.badpix)        keyFile.set_boolean ("Color appearance", "Badpix",        colorappearance.badpix);
    if (!pedited || pedited->colorappearance.datacie)       keyFile.set_boolean ("Color appearance", "Datacie",       colorappearance.datacie);
    if (!pedited || pedited->colorappearance.tonecie)       keyFile.set_boolean ("Color appearance", "Tonecie",       colorappearance.tonecie);
//    if (!pedited || pedited->colorappearance.sharpcie)      keyFile.set_boolean ("Color appearance", "Sharpcie",      colorappearance.sharpcie);
    if (!pedited || pedited->colorappearance.curveMode)  {
        Glib::ustring method;
        switch (colorappearance.curveMode) {
        case (ColorAppearanceParams::TC_MODE_LIGHT):
            method = "Lightness";
            break;
        case (ColorAppearanceParams::TC_MODE_BRIGHT):
            method = "Brightness";
            break;
        }
        keyFile.set_string  ("Color appearance", "CurveMode", method);
    }
    if (!pedited || pedited->colorappearance.curveMode2)  {
        Glib::ustring method;
        switch (colorappearance.curveMode2) {
        case (ColorAppearanceParams::TC_MODE_LIGHT):
            method = "Lightness";
            break;
        case (ColorAppearanceParams::TC_MODE_BRIGHT):
            method = "Brightness";
            break;
        }
        keyFile.set_string  ("Color appearance", "CurveMode2", method);
    }
    if (!pedited || pedited->colorappearance.curveMode3)  {
        Glib::ustring method;
        switch (colorappearance.curveMode3) {
        case (ColorAppearanceParams::TC_MODE_CHROMA):
            method = "Chroma";
            break;
        case (ColorAppearanceParams::TC_MODE_SATUR):
            method = "Saturation";
            break;
        case (ColorAppearanceParams::TC_MODE_COLORF):
            method = "Colorfullness";
            break;
			
        }
        keyFile.set_string  ("Color appearance", "CurveMode3", method);
    }
	
    if (!pedited || pedited->colorappearance.curve) {
        Glib::ArrayHandle<double> tcurve = colorappearance.curve;
        keyFile.set_double_list("Color appearance", "Curve", tcurve);
    }
    if (!pedited || pedited->colorappearance.curve2) {
        Glib::ArrayHandle<double> tcurve = colorappearance.curve2;
        keyFile.set_double_list("Color appearance", "Curve2", tcurve);
    }
    if (!pedited || pedited->colorappearance.curve3) {
        Glib::ArrayHandle<double> tcurve = colorappearance.curve3;
        keyFile.set_double_list("Color appearance", "Curve3", tcurve);
    }


	
    // save impulseDenoise
    if (!pedited || pedited->impulseDenoise.enabled) keyFile.set_boolean ("Impulse Denoising", "Enabled",   impulseDenoise.enabled);
    if (!pedited || pedited->impulseDenoise.thresh)  keyFile.set_integer ("Impulse Denoising", "Threshold", impulseDenoise.thresh);

    // save defringe
    if (!pedited || pedited->defringe.enabled)       keyFile.set_boolean ("Defringing", "Enabled",   defringe.enabled);
    if (!pedited || pedited->defringe.radius)        keyFile.set_double  ("Defringing", "Radius",    defringe.radius);
    if (!pedited || pedited->defringe.threshold)     keyFile.set_integer ("Defringing", "Threshold", defringe.threshold);
    if (!pedited || pedited->defringe.huecurve)  {
        Glib::ArrayHandle<double> huecurve = defringe.huecurve;
        keyFile.set_double_list("Defringing", "HueCurve", huecurve);
    }

    // save dirpyrDenoise
    if (!pedited || pedited->dirpyrDenoise.enabled) keyFile.set_boolean ("Directional Pyramid Denoising", "Enabled", dirpyrDenoise.enabled);
    if (!pedited || pedited->dirpyrDenoise.enhance) keyFile.set_boolean ("Directional Pyramid Denoising", "Enhance", dirpyrDenoise.enhance);
    if (!pedited || pedited->dirpyrDenoise.median) keyFile.set_boolean ("Directional Pyramid Denoising", "Median", dirpyrDenoise.median);
    if (!pedited || pedited->dirpyrDenoise.autochroma) keyFile.set_boolean ("Directional Pyramid Denoising", "Auto", dirpyrDenoise.autochroma);
 //   if (!pedited || pedited->dirpyrDenoise.perform) keyFile.set_boolean ("Directional Pyramid Denoising", "Perform", dirpyrDenoise.perform);
    if (!pedited || pedited->dirpyrDenoise.luma)    keyFile.set_double ("Directional Pyramid Denoising", "Luma",    dirpyrDenoise.luma);
    if (!pedited || pedited->dirpyrDenoise.Ldetail) keyFile.set_double ("Directional Pyramid Denoising", "Ldetail", dirpyrDenoise.Ldetail);
    if (!pedited || pedited->dirpyrDenoise.chroma)  keyFile.set_double ("Directional Pyramid Denoising", "Chroma",  dirpyrDenoise.chroma);
    if (!pedited || pedited->dirpyrDenoise.dmethod)  keyFile.set_string  ("Directional Pyramid Denoising", "Method",  dirpyrDenoise.dmethod);
    if (!pedited || pedited->dirpyrDenoise.Lmethod)  keyFile.set_string  ("Directional Pyramid Denoising", "LMethod",  dirpyrDenoise.Lmethod);
    // never save 'auto chroma preview mode' to pp3
	if (!pedited || pedited->dirpyrDenoise.Cmethod) {
		if(dirpyrDenoise.Cmethod=="PRE")
			dirpyrDenoise.Cmethod = "MAN";
		keyFile.set_string  ("Directional Pyramid Denoising", "CMethod",  dirpyrDenoise.Cmethod);
    }
    if (!pedited || pedited->dirpyrDenoise.C2method) {
    	if(dirpyrDenoise.C2method=="PREV")
			dirpyrDenoise.C2method = "MANU";
		keyFile.set_string  ("Directional Pyramid Denoising", "C2Method",  dirpyrDenoise.C2method);
    }
    if (!pedited || pedited->dirpyrDenoise.smethod)  keyFile.set_string  ("Directional Pyramid Denoising", "SMethod",  dirpyrDenoise.smethod);
    if (!pedited || pedited->dirpyrDenoise.medmethod)  keyFile.set_string  ("Directional Pyramid Denoising", "MedMethod",  dirpyrDenoise.medmethod);
    if (!pedited || pedited->dirpyrDenoise.rgbmethod)  keyFile.set_string  ("Directional Pyramid Denoising", "RGBMethod",  dirpyrDenoise.rgbmethod);
    if (!pedited || pedited->dirpyrDenoise.methodmed)  keyFile.set_string  ("Directional Pyramid Denoising", "MethodMed",  dirpyrDenoise.methodmed);
	if (!pedited || pedited->dirpyrDenoise.redchro) keyFile.set_double ("Directional Pyramid Denoising", "Redchro",  dirpyrDenoise.redchro);
    if (!pedited || pedited->dirpyrDenoise.bluechro)keyFile.set_double ("Directional Pyramid Denoising", "Bluechro",  dirpyrDenoise.bluechro);
    if (!pedited || pedited->dirpyrDenoise.gamma)   keyFile.set_double  ("Directional Pyramid Denoising", "Gamma",   dirpyrDenoise.gamma);
    if (!pedited || pedited->dirpyrDenoise.passes)   keyFile.set_integer  ("Directional Pyramid Denoising", "Passes",   dirpyrDenoise.passes);
    if (!pedited || pedited->dirpyrDenoise.lcurve)  {
        Glib::ArrayHandle<double> lcurve = dirpyrDenoise.lcurve;
        keyFile.set_double_list("Directional Pyramid Denoising", "LCurve", lcurve);
    }
    if (!pedited || pedited->dirpyrDenoise.cccurve)  {
        Glib::ArrayHandle<double> cccurve = dirpyrDenoise.cccurve;
        keyFile.set_double_list("Directional Pyramid Denoising", "CCCurve", cccurve);
    }

    //Save epd.
    if (!pedited || pedited->epd.enabled)             keyFile.set_boolean ("EPD", "Enabled", epd.enabled);
    if (!pedited || pedited->epd.strength)            keyFile.set_double  ("EPD", "Strength", epd.strength);
    if (!pedited || pedited->epd.gamma)          	  keyFile.set_double  ("EPD", "Gamma", epd.gamma);
    if (!pedited || pedited->epd.edgeStopping)        keyFile.set_double  ("EPD", "EdgeStopping", epd.edgeStopping);
    if (!pedited || pedited->epd.scale)               keyFile.set_double  ("EPD", "Scale", epd.scale);
    if (!pedited || pedited->epd.reweightingIterates) keyFile.set_integer ("EPD", "ReweightingIterates", epd.reweightingIterates);

/*
    // save lumaDenoise
    if (!pedited || pedited->lumaDenoise.enabled)       keyFile.set_boolean ("Luminance Denoising", "Enabled",       lumaDenoise.enabled);
    if (!pedited || pedited->lumaDenoise.radius)        keyFile.set_double  ("Luminance Denoising", "Radius",        lumaDenoise.radius);
    if (!pedited || pedited->lumaDenoise.edgetolerance) keyFile.set_integer ("Luminance Denoising", "EdgeTolerance", lumaDenoise.edgetolerance);
*/

/*
    // save colorDenoise
    //if (!pedited || pedited->colorDenoise.enabled)      keyFile.set_boolean ("Chrominance Denoising", "Enabled", colorDenoise.enabled);
    if (!pedited || pedited->colorDenoise.amount)       keyFile.set_integer ("Chrominance Denoising", "Amount",  colorDenoise.amount);
*/

    // save sh
    if (!pedited || pedited->sh.enabled)       keyFile.set_boolean ("Shadows & Highlights", "Enabled",             sh.enabled);
    if (!pedited || pedited->sh.hq)            keyFile.set_boolean ("Shadows & Highlights", "HighQuality",         sh.hq);
    if (!pedited || pedited->sh.highlights)    keyFile.set_integer ("Shadows & Highlights", "Highlights",          sh.highlights);
    if (!pedited || pedited->sh.htonalwidth)   keyFile.set_integer ("Shadows & Highlights", "HighlightTonalWidth", sh.htonalwidth);
    if (!pedited || pedited->sh.shadows)       keyFile.set_integer ("Shadows & Highlights", "Shadows",             sh.shadows);
    if (!pedited || pedited->sh.stonalwidth)   keyFile.set_integer ("Shadows & Highlights", "ShadowTonalWidth",    sh.stonalwidth);
    if (!pedited || pedited->sh.localcontrast) keyFile.set_integer ("Shadows & Highlights", "LocalContrast",       sh.localcontrast);
    if (!pedited || pedited->sh.radius)        keyFile.set_integer ("Shadows & Highlights", "Radius",              sh.radius);

    // save crop
    if (!pedited || pedited->crop.enabled)     keyFile.set_boolean ("Crop", "Enabled",     crop.enabled);
    if (!pedited || pedited->crop.x)           keyFile.set_integer ("Crop", "X",           crop.x);
    if (!pedited || pedited->crop.y)           keyFile.set_integer ("Crop", "Y",           crop.y);
    if (!pedited || pedited->crop.w)           keyFile.set_integer ("Crop", "W",           crop.w);
    if (!pedited || pedited->crop.h)           keyFile.set_integer ("Crop", "H",           crop.h);
    if (!pedited || pedited->crop.fixratio)    keyFile.set_boolean ("Crop", "FixedRatio",  crop.fixratio);
    if (!pedited || pedited->crop.ratio)       keyFile.set_string  ("Crop", "Ratio",       crop.ratio);
    if (!pedited || pedited->crop.orientation) keyFile.set_string  ("Crop", "Orientation", crop.orientation);
    if (!pedited || pedited->crop.guide)       keyFile.set_string  ("Crop", "Guide",       crop.guide);
    
    // save coarse
    if (!pedited || pedited->coarse.rotate)    keyFile.set_integer ("Coarse Transformation", "Rotate",         coarse.rotate);
    if (!pedited || pedited->coarse.hflip)     keyFile.set_boolean ("Coarse Transformation", "HorizontalFlip", coarse.hflip);
    if (!pedited || pedited->coarse.vflip)     keyFile.set_boolean ("Coarse Transformation", "VerticalFlip",   coarse.vflip);
    
    // save commonTrans
    if (!pedited || pedited->commonTrans.autofill)   keyFile.set_boolean ("Common Properties for Transformations", "AutoFill", commonTrans.autofill);

    // save rotate
    if (!pedited || pedited->rotate.degree)          keyFile.set_double  ("Rotation", "Degree", rotate.degree);

    // save distortion
    if (!pedited || pedited->distortion.amount)      keyFile.set_double  ("Distortion", "Amount", distortion.amount);

    // lens profile
    if (!pedited || pedited->lensProf.lcpFile)       keyFile.set_string  ("LensProfile", "LCPFile", relativePathIfInside(fname, fnameAbsolute, lensProf.lcpFile));
    if (!pedited || pedited->lensProf.useDist)       keyFile.set_boolean  ("LensProfile", "UseDistortion", lensProf.useDist);
    if (!pedited || pedited->lensProf.useVign)       keyFile.set_boolean  ("LensProfile", "UseVignette", lensProf.useVign);
    if (!pedited || pedited->lensProf.useCA)         keyFile.set_boolean  ("LensProfile", "UseCA", lensProf.useCA);

    // save perspective correction
    if (!pedited || pedited->perspective.horizontal) keyFile.set_double  ("Perspective", "Horizontal", perspective.horizontal);
    if (!pedited || pedited->perspective.vertical)   keyFile.set_double  ("Perspective", "Vertical",   perspective.vertical);

    // save gradient
    if (!pedited || pedited->gradient.enabled)       keyFile.set_boolean ("Gradient", "Enabled", gradient.enabled);
    if (!pedited || pedited->gradient.degree)        keyFile.set_double  ("Gradient", "Degree", gradient.degree);
    if (!pedited || pedited->gradient.feather)       keyFile.set_integer ("Gradient", "Feather", gradient.feather);
    if (!pedited || pedited->gradient.strength)      keyFile.set_double  ("Gradient", "Strength", gradient.strength);
    if (!pedited || pedited->gradient.centerX)       keyFile.set_integer ("Gradient", "CenterX", gradient.centerX);
    if (!pedited || pedited->gradient.centerY)       keyFile.set_integer ("Gradient", "CenterY", gradient.centerY);

    // save post-crop vignette
    if (!pedited || pedited->pcvignette.enabled)     keyFile.set_boolean ("PCVignette", "Enabled", pcvignette.enabled);
    if (!pedited || pedited->pcvignette.strength)    keyFile.set_double  ("PCVignette", "Strength", pcvignette.strength);
    if (!pedited || pedited->pcvignette.feather)     keyFile.set_integer ("PCVignette", "Feather", pcvignette.feather);
    if (!pedited || pedited->pcvignette.roundness)   keyFile.set_integer ("PCVignette", "Roundness", pcvignette.roundness);

    // save C/A correction
    if (!pedited || pedited->cacorrection.red)       keyFile.set_double  ("CACorrection", "Red",  cacorrection.red);
    if (!pedited || pedited->cacorrection.blue)      keyFile.set_double  ("CACorrection", "Blue", cacorrection.blue);

    // save vignetting correction
    if (!pedited || pedited->vignetting.amount)      keyFile.set_integer ("Vignetting Correction", "Amount", vignetting.amount);
    if (!pedited || pedited->vignetting.radius)      keyFile.set_integer ("Vignetting Correction", "Radius", vignetting.radius);
    if (!pedited || pedited->vignetting.strength)    keyFile.set_integer ("Vignetting Correction", "Strength", vignetting.strength);
    if (!pedited || pedited->vignetting.centerX)     keyFile.set_integer ("Vignetting Correction", "CenterX", vignetting.centerX);
    if (!pedited || pedited->vignetting.centerY)     keyFile.set_integer ("Vignetting Correction", "CenterY", vignetting.centerY);


    if (!pedited || pedited->resize.enabled)         keyFile.set_boolean ("Resize", "Enabled",resize.enabled);
    if (!pedited || pedited->resize.scale)           keyFile.set_double  ("Resize", "Scale",  resize.scale);
    if (!pedited || pedited->resize.appliesTo)       keyFile.set_string  ("Resize", "AppliesTo", resize.appliesTo);
    if (!pedited || pedited->resize.method)          keyFile.set_string  ("Resize", "Method", resize.method);
    if (!pedited || pedited->resize.dataspec)        keyFile.set_integer ("Resize", "DataSpecified",  resize.dataspec);
    if (!pedited || pedited->resize.width)           keyFile.set_integer ("Resize", "Width",  resize.width);
    if (!pedited || pedited->resize.height)          keyFile.set_integer ("Resize", "Height", resize.height);

    if (!pedited || pedited->prsharpening.enabled)            keyFile.set_boolean ("PostResizeSharpening", "Enabled",             prsharpening.enabled);
    if (!pedited || pedited->prsharpening.method)             keyFile.set_string  ("PostResizeSharpening", "Method",              prsharpening.method);
    if (!pedited || pedited->prsharpening.radius)             keyFile.set_double  ("PostResizeSharpening", "Radius",              prsharpening.radius);
    if (!pedited || pedited->prsharpening.amount)             keyFile.set_integer ("PostResizeSharpening", "Amount",              prsharpening.amount);
    if (!pedited || pedited->prsharpening.threshold) {
        Glib::ArrayHandle<int> thresh (prsharpening.threshold.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("PostResizeSharpening",   "Threshold", thresh);
    }
    if (!pedited || pedited->prsharpening.edgesonly)          keyFile.set_boolean ("PostResizeSharpening", "OnlyEdges",           prsharpening.edgesonly);
    if (!pedited || pedited->prsharpening.edges_radius)       keyFile.set_double  ("PostResizeSharpening", "EdgedetectionRadius", prsharpening.edges_radius);
    if (!pedited || pedited->prsharpening.edges_tolerance)    keyFile.set_integer ("PostResizeSharpening", "EdgeTolerance",       prsharpening.edges_tolerance);
    if (!pedited || pedited->prsharpening.halocontrol)        keyFile.set_boolean ("PostResizeSharpening", "HalocontrolEnabled",  prsharpening.halocontrol);
    if (!pedited || pedited->prsharpening.halocontrol_amount) keyFile.set_integer ("PostResizeSharpening", "HalocontrolAmount",   prsharpening.halocontrol_amount);
    if (!pedited || pedited->prsharpening.deconvradius)       keyFile.set_double  ("PostResizeSharpening", "DeconvRadius",        prsharpening.deconvradius);
    if (!pedited || pedited->prsharpening.deconvamount)       keyFile.set_integer ("PostResizeSharpening", "DeconvAmount",        prsharpening.deconvamount);
    if (!pedited || pedited->prsharpening.deconvdamping)      keyFile.set_integer ("PostResizeSharpening", "DeconvDamping",       prsharpening.deconvdamping);
    if (!pedited || pedited->prsharpening.deconviter)         keyFile.set_integer ("PostResizeSharpening", "DeconvIterations",    prsharpening.deconviter);


    // save color management settings
    if (!pedited || pedited->icm.input)              keyFile.set_string  ("Color Management", "InputProfile",   relativePathIfInside(fname, fnameAbsolute, icm.input));
    if (!pedited || pedited->icm.toneCurve)          keyFile.set_boolean ("Color Management", "ToneCurve",   icm.toneCurve);
    if (!pedited || pedited->icm.applyLookTable)     keyFile.set_boolean ("Color Management", "ApplyLookTable",   icm.applyLookTable);
    if (!pedited || pedited->icm.applyBaselineExposureOffset)     keyFile.set_boolean ("Color Management", "ApplyBaselineExposureOffset",   icm.applyBaselineExposureOffset);
    if (!pedited || pedited->icm.applyHueSatMap)     keyFile.set_boolean ("Color Management", "ApplyHueSatMap",   icm.applyHueSatMap);
    if (!pedited || pedited->icm.blendCMSMatrix)     keyFile.set_boolean ("Color Management", "BlendCMSMatrix",   icm.blendCMSMatrix);
    if (!pedited || pedited->icm.dcpIlluminant)      keyFile.set_integer ("Color Management", "DCPIlluminant",   icm.dcpIlluminant);
    if (!pedited || pedited->icm.working)            keyFile.set_string  ("Color Management", "WorkingProfile", icm.working);
    if (!pedited || pedited->icm.output)             keyFile.set_string  ("Color Management", "OutputProfile",  icm.output);
    if (!pedited || pedited->icm.gamma)              keyFile.set_string  ("Color Management", "Gammafree",  icm.gamma);
    if (!pedited || pedited->icm.freegamma)          keyFile.set_boolean ("Color Management", "Freegamma",  icm.freegamma);
    if (!pedited || pedited->icm.gampos)             keyFile.set_double  ("Color Management", "GammaValue",  icm.gampos);
    if (!pedited || pedited->icm.slpos)              keyFile.set_double  ("Color Management", "GammaSlope",  icm.slpos);



    // save wavelet parameters
    if (!pedited || pedited->wavelet.enabled)    keyFile.set_boolean ("Wavelet", "Enabled", wavelet.enabled);
	if (!pedited || pedited->wavelet.strength)    keyFile.set_integer ("Wavelet", "Strength", wavelet.strength);
	if (!pedited || pedited->wavelet.balance)    keyFile.set_integer ("Wavelet", "Balance", wavelet.balance);
	if (!pedited || pedited->wavelet.iter)    keyFile.set_integer ("Wavelet", "Iter", wavelet.iter);
	if (!pedited || pedited->wavelet.thres)  keyFile.set_integer  ("Wavelet", "MaxLev",  wavelet.thres);
    if (!pedited || pedited->wavelet.Tilesmethod)  keyFile.set_string  ("Wavelet", "TilesMethod",  wavelet.Tilesmethod);
    if (!pedited || pedited->wavelet.daubcoeffmethod)  keyFile.set_string  ("Wavelet", "DaubMethod",  wavelet.daubcoeffmethod);
    if (!pedited || pedited->wavelet.CLmethod)  keyFile.set_string  ("Wavelet", "ChoiceLevMethod",  wavelet.CLmethod);
    if (!pedited || pedited->wavelet.Backmethod)  keyFile.set_string  ("Wavelet", "BackMethod",  wavelet.Backmethod);
    if (!pedited || pedited->wavelet.Lmethod)  keyFile.set_string  ("Wavelet", "LevMethod",  wavelet.Lmethod);
    if (!pedited || pedited->wavelet.Dirmethod)  keyFile.set_string  ("Wavelet", "DirMethod",  wavelet.Dirmethod);
	if (!pedited || pedited->wavelet.greenhigh)    keyFile.set_integer ("Wavelet", "CBgreenhigh", wavelet.greenhigh);
	if (!pedited || pedited->wavelet.greenmed)    keyFile.set_integer ("Wavelet", "CBgreenmed", wavelet.greenmed);
	if (!pedited || pedited->wavelet.greenlow)    keyFile.set_integer ("Wavelet", "CBgreenlow", wavelet.greenlow);
	if (!pedited || pedited->wavelet.bluehigh)    keyFile.set_integer ("Wavelet", "CBbluehigh", wavelet.bluehigh);
	if (!pedited || pedited->wavelet.bluemed)    keyFile.set_integer ("Wavelet", "CBbluemed", wavelet.bluemed);
	if (!pedited || pedited->wavelet.bluelow)    keyFile.set_integer ("Wavelet", "CBbluelow", wavelet.bluelow);
	
    for(int i = 0; i < 9; i++)
    {
        std::stringstream ss;
        ss << "Contrast" << (i+1);
	if (!pedited || pedited->wavelet.c[i])         keyFile.set_integer("Wavelet", ss.str(), wavelet.c[i]);
    }
    for(int i = 0; i < 9; i++)
    {
        std::stringstream ss;
        ss << "Chroma" << (i+1);
	if (!pedited || pedited->wavelet.ch[i])         keyFile.set_integer("Wavelet", ss.str(), wavelet.ch[i]);
    }
	
    if (!pedited || pedited->wavelet.sup)  keyFile.set_integer  ("Wavelet", "ContExtra",  wavelet.sup);	
    if (!pedited || pedited->wavelet.HSmethod)  keyFile.set_string  ("Wavelet", "HSMethod",  wavelet.HSmethod);
    if (!pedited || pedited->wavelet.hllev) {
        Glib::ArrayHandle<int> thresh (wavelet.hllev.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Wavelet",   "HLRange", thresh);
    }
    if (!pedited || pedited->wavelet.bllev) {
        Glib::ArrayHandle<int> thresh (wavelet.bllev.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Wavelet",   "SHRange", thresh);
    }
    if (!pedited || pedited->wavelet.edgcont) {
        Glib::ArrayHandle<int> thresh (wavelet.edgcont.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Wavelet",   "Edgcont", thresh);
    }
    if (!pedited || pedited->wavelet.level0noise) {
        Glib::ArrayHandle<double> thresh (wavelet.level0noise.value, 2, Glib::OWNERSHIP_NONE);
        keyFile.set_double_list("Wavelet",   "Level0noise", thresh);
    }
    if (!pedited || pedited->wavelet.level1noise) {
        Glib::ArrayHandle<double> thresh (wavelet.level1noise.value, 2, Glib::OWNERSHIP_NONE);
        keyFile.set_double_list("Wavelet",   "Level1noise", thresh);
    }
    if (!pedited || pedited->wavelet.level2noise) {
        Glib::ArrayHandle<double> thresh (wavelet.level2noise.value, 2, Glib::OWNERSHIP_NONE);
        keyFile.set_double_list("Wavelet",   "Level2noise", thresh);
    }
	
    if (!pedited || pedited->wavelet.threshold)  keyFile.set_integer  ("Wavelet", "ThresholdHighlight",  wavelet.threshold);
    if (!pedited || pedited->wavelet.threshold2)  keyFile.set_integer  ("Wavelet", "ThresholdShadow",  wavelet.threshold2);
    if (!pedited || pedited->wavelet.edgedetect)  keyFile.set_integer  ("Wavelet", "Edgedetect",  wavelet.edgedetect);
    if (!pedited || pedited->wavelet.edgedetectthr)  keyFile.set_integer  ("Wavelet", "Edgedetectthr",  wavelet.edgedetectthr);
    if (!pedited || pedited->wavelet.edgedetectthr2)  keyFile.set_integer  ("Wavelet", "EdgedetectthrHi",  wavelet.edgedetectthr2);
    if (!pedited || pedited->wavelet.chroma)  keyFile.set_integer  ("Wavelet", "ThresholdChroma",  wavelet.chroma);
    if (!pedited || pedited->wavelet.CHmethod)  keyFile.set_string  ("Wavelet", "CHromaMethod",  wavelet.CHmethod);
    if (!pedited || pedited->wavelet.Medgreinf)  keyFile.set_string  ("Wavelet", "Medgreinf",  wavelet.Medgreinf);
    if (!pedited || pedited->wavelet.CHSLmethod)  keyFile.set_string  ("Wavelet", "CHSLromaMethod",  wavelet.CHSLmethod);
    if (!pedited || pedited->wavelet.EDmethod)  keyFile.set_string  ("Wavelet", "EDMethod",  wavelet.EDmethod);
    if (!pedited || pedited->wavelet.BAmethod)  keyFile.set_string  ("Wavelet", "BAMethod",  wavelet.BAmethod);
    if (!pedited || pedited->wavelet.TMmethod)  keyFile.set_string  ("Wavelet", "TMMethod",  wavelet.TMmethod);
    if (!pedited || pedited->wavelet.chro)  keyFile.set_integer  ("Wavelet", "ChromaLink",  wavelet.chro);
    if (!pedited || pedited->wavelet.ccwcurve)  {
        Glib::ArrayHandle<double> ccwcurve = wavelet.ccwcurve;
        keyFile.set_double_list("Wavelet", "ContrastCurve", ccwcurve);
    }
    if (!pedited || pedited->wavelet.pastlev) {
        Glib::ArrayHandle<int> thresh (wavelet.pastlev.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Wavelet",   "Pastlev", thresh);
    }
    if (!pedited || pedited->wavelet.satlev) {
        Glib::ArrayHandle<int> thresh (wavelet.satlev.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Wavelet",   "Satlev", thresh);
    }
	
    if (!pedited || pedited->wavelet.opacityCurveRG) {
        Glib::ArrayHandle<double> curve = wavelet.opacityCurveRG;
        keyFile.set_double_list("Wavelet", "OpacityCurveRG", curve);
    }
    if (!pedited || pedited->wavelet.opacityCurveBY) {
        Glib::ArrayHandle<double> curve = wavelet.opacityCurveBY;
        keyFile.set_double_list("Wavelet", "OpacityCurveBY", curve);
    }
    if (!pedited || pedited->wavelet.opacityCurveW) {
        Glib::ArrayHandle<double> curve = wavelet.opacityCurveW;
        keyFile.set_double_list("Wavelet", "OpacityCurveW", curve);
    }
    if (!pedited || pedited->wavelet.opacityCurveWL) {
        Glib::ArrayHandle<double> curve = wavelet.opacityCurveWL;
        keyFile.set_double_list("Wavelet", "OpacityCurveWL", curve);
    }
	
    if (!pedited || pedited->wavelet.hhcurve) {
        Glib::ArrayHandle<double> curve = wavelet.hhcurve;
        keyFile.set_double_list("Wavelet", "HHcurve", curve);
    }
    if (!pedited || pedited->wavelet.Chcurve) {
        Glib::ArrayHandle<double> curve = wavelet.Chcurve;
        keyFile.set_double_list("Wavelet", "CHcurve", curve);
    }
    if (!pedited || pedited->wavelet.wavclCurve)  {
        Glib::ArrayHandle<double> wavclCurve = wavelet.wavclCurve;
        keyFile.set_double_list("Wavelet", "WavclCurve", wavclCurve);
    }
	
	
    if (!pedited || pedited->wavelet.median)    keyFile.set_boolean ("Wavelet", "Median", wavelet.median);
    if (!pedited || pedited->wavelet.medianlev)    keyFile.set_boolean ("Wavelet", "Medianlev", wavelet.medianlev);
    if (!pedited || pedited->wavelet.linkedg)    keyFile.set_boolean ("Wavelet", "Linkedg", wavelet.linkedg);
    if (!pedited || pedited->wavelet.cbenab)    keyFile.set_boolean ("Wavelet", "CBenab", wavelet.cbenab);
    if (!pedited || pedited->wavelet.lipst)    keyFile.set_boolean ("Wavelet", "Lipst", wavelet.lipst);
 //   if (!pedited || pedited->wavelet.edgreinf)    keyFile.set_boolean ("Wavelet", "Edgreinf", wavelet.edgreinf);
    if (!pedited || pedited->wavelet.skinprotect) keyFile.set_double ("Wavelet", "Skinprotect", wavelet.skinprotect);
     if (!pedited || pedited->wavelet.hueskin) {
        Glib::ArrayHandle<int> thresh (wavelet.hueskin.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Wavelet",   "Hueskin", thresh);
    }
    
	if (!pedited || pedited->wavelet.edgrad)  keyFile.set_integer  ("Wavelet", "Edgrad",  wavelet.edgrad);
    if (!pedited || pedited->wavelet.edgval)  keyFile.set_integer  ("Wavelet", "Edgval",  wavelet.edgval);
    if (!pedited || pedited->wavelet.edgthresh)  keyFile.set_integer  ("Wavelet", "ThrEdg",  wavelet.edgthresh);
 //   if (!pedited || pedited->wavelet.strength)  keyFile.set_integer  ("Wavelet", "Strength",  wavelet.strength);
  //  if (!pedited || pedited->wavelet.balance)  keyFile.set_integer  ("Wavelet", "Balance",  wavelet.balance);

	if (!pedited || pedited->wavelet.avoid)    keyFile.set_boolean ("Wavelet", "AvoidColorShift", wavelet.avoid);
	if (!pedited || pedited->wavelet.tmr)    keyFile.set_boolean ("Wavelet", "TMr", wavelet.tmr);
    if (!pedited || pedited->wavelet.rescon)  keyFile.set_integer  ("Wavelet", "ResidualcontShadow",  wavelet.rescon);
    if (!pedited || pedited->wavelet.resconH)  keyFile.set_integer  ("Wavelet", "ResidualcontHighlight",  wavelet.resconH);
    if (!pedited || pedited->wavelet.thr)  keyFile.set_integer  ("Wavelet", "ThresholdResidShadow",  wavelet.thr);
    if (!pedited || pedited->wavelet.thrH)  keyFile.set_integer  ("Wavelet", "ThresholdResidHighLight",  wavelet.thrH);
    if (!pedited || pedited->wavelet.reschro)  keyFile.set_integer  ("Wavelet", "Residualchroma",  wavelet.reschro);
    if (!pedited || pedited->wavelet.tmrs)  keyFile.set_double  ("Wavelet", "ResidualTM",  wavelet.tmrs);
    if (!pedited || pedited->wavelet.gamma)  keyFile.set_double  ("Wavelet", "Residualgamma",  wavelet.gamma);
    if (!pedited || pedited->wavelet.sky)  keyFile.set_double  ("Wavelet", "HueRangeResidual",  wavelet.sky);
    if (!pedited || pedited->wavelet.hueskin2) {
        Glib::ArrayHandle<int> thresh (wavelet.hueskin2.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Wavelet",   "HueRange", thresh);
    }
    if (!pedited || pedited->wavelet.contrast)  keyFile.set_integer  ("Wavelet", "Contrast",  wavelet.contrast);
	
    
    // save directional pyramid wavelet parameters
    if (!pedited || pedited->dirpyrequalizer.enabled) keyFile.set_boolean ("Directional Pyramid Equalizer", "Enabled", dirpyrequalizer.enabled);
    if (!pedited || pedited->dirpyrequalizer.gamutlab) keyFile.set_boolean ("Directional Pyramid Equalizer", "Gamutlab", dirpyrequalizer.gamutlab);
    for(int i = 0; i < 6; i++)
    {
        std::stringstream ss;
        ss << "Mult" << i;
        if (!pedited || pedited->dirpyrequalizer.mult[i]) keyFile.set_double("Directional Pyramid Equalizer", ss.str(), dirpyrequalizer.mult[i]);
    }
    if (!pedited || pedited->dirpyrequalizer.threshold) keyFile.set_double ("Directional Pyramid Equalizer", "Threshold", dirpyrequalizer.threshold);
    if (!pedited || pedited->dirpyrequalizer.skinprotect) keyFile.set_double ("Directional Pyramid Equalizer", "Skinprotect", dirpyrequalizer.skinprotect);
  //  if (!pedited || pedited->dirpyrequalizer.algo) keyFile.set_string ("Directional Pyramid Equalizer", "Algorithm", dirpyrequalizer.algo);
    if (!pedited || pedited->dirpyrequalizer.hueskin) {
        Glib::ArrayHandle<int> thresh (dirpyrequalizer.hueskin.value, 4, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("Directional Pyramid Equalizer",   "Hueskin", thresh);
    }

    // save hsv wavelet parameters
    if (!pedited || pedited->hsvequalizer.hcurve) {
        Glib::ArrayHandle<double> hcurve = hsvequalizer.hcurve;
        keyFile.set_double_list("HSV Equalizer", "HCurve", hcurve);
    }
    if (!pedited || pedited->hsvequalizer.scurve) {
        Glib::ArrayHandle<double> scurve = hsvequalizer.scurve;
        keyFile.set_double_list("HSV Equalizer", "SCurve", scurve);
    }
    if (!pedited || pedited->hsvequalizer.vcurve) {
        Glib::ArrayHandle<double> vcurve = hsvequalizer.vcurve;
        keyFile.set_double_list("HSV Equalizer", "VCurve", vcurve);
    }

    //save film simulation parameters
    if ( !pedited || pedited->filmSimulation.enabled )      keyFile.set_boolean( "Film Simulation", "Enabled", filmSimulation.enabled );
    if ( !pedited || pedited->filmSimulation.clutFilename ) keyFile.set_string ( "Film Simulation", "ClutFilename", filmSimulation.clutFilename );
    if ( !pedited || pedited->filmSimulation.strength )     keyFile.set_integer( "Film Simulation", "Strength", filmSimulation.strength );


    if (!pedited || pedited->rgbCurves.lumamode)     keyFile.set_boolean ("RGB Curves", "LumaMode",  rgbCurves.lumamode);

    if (!pedited || pedited->rgbCurves.rcurve) {
        Glib::ArrayHandle<double> RGBrcurve = rgbCurves.rcurve;
        keyFile.set_double_list("RGB Curves", "rCurve", RGBrcurve);
    }
    if (!pedited || pedited->rgbCurves.gcurve) {
        Glib::ArrayHandle<double> RGBgcurve = rgbCurves.gcurve;
        keyFile.set_double_list("RGB Curves", "gCurve", RGBgcurve);
    }
    if (!pedited || pedited->rgbCurves.bcurve) {
        Glib::ArrayHandle<double> RGBbcurve = rgbCurves.bcurve;
        keyFile.set_double_list("RGB Curves", "bCurve", RGBbcurve);
    }

    // save Color Toning
    if (!pedited || pedited->colorToning.enabled)    keyFile.set_boolean ("ColorToning", "Enabled", colorToning.enabled);
    if (!pedited || pedited->colorToning.method)     keyFile.set_string  ("ColorToning", "Method", colorToning.method);
    if (!pedited || pedited->colorToning.lumamode)   keyFile.set_boolean ("ColorToning", "Lumamode", colorToning.lumamode);
    if (!pedited || pedited->colorToning.twocolor)   keyFile.set_string  ("ColorToning", "Twocolor", colorToning.twocolor);
    if (!pedited || pedited->colorToning.redlow)     keyFile.set_double  ("ColorToning", "Redlow", colorToning.redlow);
    if (!pedited || pedited->colorToning.greenlow)   keyFile.set_double  ("ColorToning", "Greenlow", colorToning.greenlow);
    if (!pedited || pedited->colorToning.bluelow)    keyFile.set_double  ("ColorToning", "Bluelow", colorToning.bluelow);
    if (!pedited || pedited->colorToning.satlow)     keyFile.set_double  ("ColorToning", "Satlow", colorToning.satlow);
    if (!pedited || pedited->colorToning.balance)    keyFile.set_integer ("ColorToning", "Balance", colorToning.balance);
    if (!pedited || pedited->colorToning.sathigh)    keyFile.set_double  ("ColorToning", "Sathigh", colorToning.sathigh);
    if (!pedited || pedited->colorToning.redmed)     keyFile.set_double  ("ColorToning", "Redmed", colorToning.redmed);
    if (!pedited || pedited->colorToning.greenmed)   keyFile.set_double  ("ColorToning", "Greenmed", colorToning.greenmed);
    if (!pedited || pedited->colorToning.bluemed)    keyFile.set_double  ("ColorToning", "Bluemed", colorToning.bluemed);
    if (!pedited || pedited->colorToning.redhigh)    keyFile.set_double  ("ColorToning", "Redhigh", colorToning.redhigh);
    if (!pedited || pedited->colorToning.greenhigh)  keyFile.set_double  ("ColorToning", "Greenhigh", colorToning.greenhigh);
    if (!pedited || pedited->colorToning.bluehigh)   keyFile.set_double  ("ColorToning", "Bluehigh", colorToning.bluehigh);
    if (!pedited || pedited->colorToning.autosat)    keyFile.set_boolean ("ColorToning", "Autosat", colorToning.autosat);

    if (!pedited || pedited->colorToning.opacityCurve) {
        Glib::ArrayHandle<double> curve = colorToning.opacityCurve;
        keyFile.set_double_list("ColorToning", "OpacityCurve", curve);
    }
    if (!pedited || pedited->colorToning.colorCurve) {
        Glib::ArrayHandle<double> curve = colorToning.colorCurve;
        keyFile.set_double_list("ColorToning", "ColorCurve", curve);
    }
    if (!pedited || pedited->colorToning.satprotectionthreshold)  keyFile.set_integer ("ColorToning", "SatProtectionThreshold", colorToning.satProtectionThreshold );
    if (!pedited || pedited->colorToning.saturatedopacity)        keyFile.set_integer ("ColorToning", "SaturatedOpacity", colorToning.saturatedOpacity );
    if (!pedited || pedited->colorToning.strength)                keyFile.set_integer ("ColorToning", "Strength", colorToning.strength );

    if (!pedited || pedited->colorToning.hlColSat) {
        Glib::ArrayHandle<int> thresh (colorToning.hlColSat.value, 2, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("ColorToning", "HighlightsColorSaturation", thresh);
    }
    if (!pedited || pedited->colorToning.shadowsColSat) {
        Glib::ArrayHandle<int> thresh (colorToning.shadowsColSat.value, 2, Glib::OWNERSHIP_NONE);
        keyFile.set_integer_list("ColorToning", "ShadowsColorSaturation", thresh);
    }
    if (!pedited || pedited->colorToning.clcurve)  {
        Glib::ArrayHandle<double> clcurve = colorToning.clcurve;
        keyFile.set_double_list("ColorToning", "ClCurve", clcurve);
    }
    if (!pedited || pedited->colorToning.cl2curve)  {
        Glib::ArrayHandle<double> cl2curve = colorToning.cl2curve;
        keyFile.set_double_list("ColorToning", "Cl2Curve", cl2curve);
    }

    // save raw parameters
    if (!pedited || pedited->raw.darkFrame)            keyFile.set_string  ("RAW", "DarkFrame", relativePathIfInside(fname, fnameAbsolute, raw.dark_frame) );
    if (!pedited || pedited->raw.dfAuto)               keyFile.set_boolean ("RAW", "DarkFrameAuto", raw.df_autoselect );
    if (!pedited || pedited->raw.ff_file)              keyFile.set_string  ("RAW", "FlatFieldFile", relativePathIfInside(fname, fnameAbsolute, raw.ff_file) );
    if (!pedited || pedited->raw.ff_AutoSelect)        keyFile.set_boolean ("RAW", "FlatFieldAutoSelect", raw.ff_AutoSelect );
    if (!pedited || pedited->raw.ff_BlurRadius)        keyFile.set_integer ("RAW", "FlatFieldBlurRadius", raw.ff_BlurRadius );
    if (!pedited || pedited->raw.ff_BlurType)          keyFile.set_string  ("RAW", "FlatFieldBlurType", raw.ff_BlurType );
    if (!pedited || pedited->raw.ff_AutoClipControl)   keyFile.set_boolean ("RAW", "FlatFieldAutoClipControl", raw.ff_AutoClipControl );
    if (!pedited || pedited->raw.ff_clipControl)       keyFile.set_boolean ("RAW", "FlatFieldClipControl", raw.ff_clipControl );
    if (!pedited || pedited->raw.caCorrection)         keyFile.set_boolean ("RAW", "CA", raw.ca_autocorrect );
    if (!pedited || pedited->raw.caRed)                keyFile.set_double  ("RAW", "CARed", raw.cared );
    if (!pedited || pedited->raw.caBlue)               keyFile.set_double  ("RAW", "CABlue", raw.cablue );
    if (!pedited || pedited->raw.hotPixelFilter)       keyFile.set_boolean ("RAW", "HotPixelFilter", raw.hotPixelFilter );
    if (!pedited || pedited->raw.deadPixelFilter)      keyFile.set_boolean ("RAW", "DeadPixelFilter", raw.deadPixelFilter );
    if (!pedited || pedited->raw.hotDeadPixelThresh)   keyFile.set_integer ("RAW", "HotDeadPixelThresh", raw.hotdeadpix_thresh );

    if (!pedited || pedited->raw.bayersensor.method)          keyFile.set_string  ("RAW Bayer", "Method", raw.bayersensor.method );
    if (!pedited || pedited->raw.bayersensor.ccSteps)         keyFile.set_integer ("RAW Bayer", "CcSteps", raw.bayersensor.ccSteps);
    if (!pedited || pedited->raw.bayersensor.exBlack0)        keyFile.set_double  ("RAW Bayer", "PreBlack0", raw.bayersensor.black0 );
    if (!pedited || pedited->raw.bayersensor.exBlack1)        keyFile.set_double  ("RAW Bayer", "PreBlack1", raw.bayersensor.black1 );
    if (!pedited || pedited->raw.bayersensor.exBlack2)        keyFile.set_double  ("RAW Bayer", "PreBlack2", raw.bayersensor.black2 );
    if (!pedited || pedited->raw.bayersensor.exBlack3)        keyFile.set_double  ("RAW Bayer", "PreBlack3", raw.bayersensor.black3 );
    if (!pedited || pedited->raw.bayersensor.exTwoGreen)      keyFile.set_boolean ("RAW Bayer", "PreTwoGreen", raw.bayersensor.twogreen );
    if (!pedited || pedited->raw.bayersensor.linenoise)       keyFile.set_integer ("RAW Bayer", "LineDenoise", raw.bayersensor.linenoise);
    if (!pedited || pedited->raw.bayersensor.greenEq)         keyFile.set_integer ("RAW Bayer", "GreenEqThreshold", raw.bayersensor.greenthresh);
    if (!pedited || pedited->raw.bayersensor.dcbIterations)   keyFile.set_integer ("RAW Bayer", "DCBIterations", raw.bayersensor.dcb_iterations );
    if (!pedited || pedited->raw.bayersensor.dcbEnhance)      keyFile.set_boolean ("RAW Bayer", "DCBEnhance", raw.bayersensor.dcb_enhance );
    if (!pedited || pedited->raw.bayersensor.lmmseIterations) keyFile.set_integer ("RAW Bayer", "LMMSEIterations", raw.bayersensor.lmmse_iterations );
    //if (!pedited || pedited->raw.bayersensor.allEnhance)    keyFile.set_boolean ("RAW Bayer", "ALLEnhance", raw.bayersensor.all_enhance );

    if (!pedited || pedited->raw.xtranssensor.method)         keyFile.set_string  ("RAW X-Trans", "Method", raw.xtranssensor.method );
    if (!pedited || pedited->raw.xtranssensor.ccSteps)        keyFile.set_integer ("RAW X-Trans", "CcSteps", raw.xtranssensor.ccSteps);
    if (!pedited || pedited->raw.xtranssensor.exBlackRed)     keyFile.set_double  ("RAW X-Trans", "PreBlackRed", raw.xtranssensor.blackred );
    if (!pedited || pedited->raw.xtranssensor.exBlackGreen)   keyFile.set_double  ("RAW X-Trans", "PreBlackGreen", raw.xtranssensor.blackgreen );
    if (!pedited || pedited->raw.xtranssensor.exBlackBlue)    keyFile.set_double  ("RAW X-Trans", "PreBlackBlue", raw.xtranssensor.blackblue );


    // save raw exposition
    if (!pedited || pedited->raw.exPos)              keyFile.set_double  ("RAW", "PreExposure", raw.expos );
    if (!pedited || pedited->raw.exPreser)           keyFile.set_double  ("RAW", "PrePreserv", raw.preser );

    // save exif change list
    if (!pedited || pedited->exif) {
        for (ExifPairs::const_iterator i=exif.begin(); i!=exif.end(); i++)
            keyFile.set_string ("Exif", i->first, i->second);
    }

    // save iptc change list
    if (!pedited || pedited->iptc) {
        for (IPTCPairs::const_iterator i=iptc.begin(); i!=iptc.end(); i++) {
            Glib::ArrayHandle<Glib::ustring> values = i->second;
            keyFile.set_string_list ("IPTC", i->first, values);
        }
    }
    
    Glib::ustring sPParams = keyFile.to_data();

    int error1, error2;
    error1 = write (fname , sPParams);
    if (fname2.length()) {

        error2 = write (fname2, sPParams);
        // If at least one file has been saved, it's a success
        return error1 & error2;
    }
    else
        return error1;
}

int ProcParams::write (Glib::ustring &fname, Glib::ustring &content) const {

    int error = 0;
    if (fname.length()) {
        FILE *f;
        f = safe_g_fopen (fname, "wt");

        if (f==NULL)
            error = 1;
        else {
            fprintf (f, "%s", content.c_str());
            fclose (f);
        }
    }
    return error;
}

int ProcParams::load (Glib::ustring fname, ParamsEdited* pedited) {
    setlocale(LC_NUMERIC, "C"); // to set decimal point to "."
    if (fname.empty())
        return 1;

    SafeKeyFile keyFile;
    try {
        //setDefaults ();
        if (pedited)
            pedited->set(false);

        FILE* f = safe_g_fopen (fname, "rt");
        if (!f)
            return 1;
        char* buffer = new char[1024];
        std::ostringstream ostr;
        while (fgets (buffer, 1024, f))
            ostr << buffer << "\n";
        delete [] buffer;
        if (!keyFile.load_from_data (ostr.str())) 
            return 1;
        fclose (f);

    // load tonecurve:

ppVersion = PPVERSION;
appVersion = APPVERSION;
if (keyFile.has_group ("Version")) {    
    if (keyFile.has_key ("Version", "AppVersion")) appVersion = keyFile.get_string  ("Version", "AppVersion");
    if (keyFile.has_key ("Version", "Version"))    ppVersion  = keyFile.get_integer ("Version", "Version");
}
//printf("ProcParams::load called ppVersion=%i\n",ppVersion);

if (keyFile.has_group ("General")) {
    if (keyFile.has_key ("General", "Rank"))        { rank       = keyFile.get_integer ("General", "Rank"); if (pedited) pedited->general.rank = true; }
    if (keyFile.has_key ("General", "ColorLabel"))  { colorlabel = keyFile.get_integer ("General", "ColorLabel"); if (pedited) pedited->general.colorlabel = true; }
    if (keyFile.has_key ("General", "InTrash"))     { inTrash    = keyFile.get_boolean ("General", "InTrash"); if (pedited) pedited->general.intrash = true; }
}

if (keyFile.has_group ("Exposure")) {    
    if (ppVersion<PPVERSION_AEXP)
        toneCurve.autoexp = false; // prevent execution of autoexp when opening file created with earlier verions of autoexp algorithm
    else
        if (keyFile.has_key ("Exposure", "Auto"))           { toneCurve.autoexp       = keyFile.get_boolean ("Exposure", "Auto"); if (pedited) pedited->toneCurve.autoexp = true; }

    if (keyFile.has_key ("Exposure", "Clip"))           { toneCurve.clip          = keyFile.get_double  ("Exposure", "Clip"); if (pedited) pedited->toneCurve.clip = true; }
    if (keyFile.has_key ("Exposure", "Compensation"))   { toneCurve.expcomp       = keyFile.get_double  ("Exposure", "Compensation"); if (pedited) pedited->toneCurve.expcomp = true; }
    if (keyFile.has_key ("Exposure", "Brightness"))     { toneCurve.brightness    = keyFile.get_integer ("Exposure", "Brightness"); if (pedited) pedited->toneCurve.brightness = true; }
    if (keyFile.has_key ("Exposure", "Contrast"))       { toneCurve.contrast      = keyFile.get_integer ("Exposure", "Contrast"); if (pedited) pedited->toneCurve.contrast = true; }
    if (keyFile.has_key ("Exposure", "Saturation"))     { toneCurve.saturation    = keyFile.get_integer ("Exposure", "Saturation"); if (pedited) pedited->toneCurve.saturation = true; }
    if (keyFile.has_key ("Exposure", "Black"))          { toneCurve.black         = keyFile.get_integer ("Exposure", "Black"); if (pedited) pedited->toneCurve.black = true; }
    if (keyFile.has_key ("Exposure", "HighlightCompr")) { toneCurve.hlcompr       = keyFile.get_integer ("Exposure", "HighlightCompr"); if (pedited) pedited->toneCurve.hlcompr = true; }
    if (keyFile.has_key ("Exposure", "HighlightComprThreshold")) { toneCurve.hlcomprthresh = keyFile.get_integer ("Exposure", "HighlightComprThreshold"); if (pedited) pedited->toneCurve.hlcomprthresh = true; }
    if (keyFile.has_key ("Exposure", "ShadowCompr"))    { toneCurve.shcompr       = keyFile.get_integer ("Exposure", "ShadowCompr"); if (pedited) pedited->toneCurve.shcompr = true; }
    // load highlight recovery settings
    if (toneCurve.shcompr > 100) toneCurve.shcompr = 100; // older pp3 files can have values above 100.
    if (keyFile.has_key ("Exposure", "CurveMode"))      {
        Glib::ustring sMode = keyFile.get_string ("Exposure", "CurveMode");
        if      (sMode == "Standard")            toneCurve.curveMode = ToneCurveParams::TC_MODE_STD;
        else if (sMode == "FilmLike")            toneCurve.curveMode = ToneCurveParams::TC_MODE_FILMLIKE;
        else if (sMode == "SatAndValueBlending") toneCurve.curveMode = ToneCurveParams::TC_MODE_SATANDVALBLENDING;
        else if (sMode == "WeightedStd")         toneCurve.curveMode = ToneCurveParams::TC_MODE_WEIGHTEDSTD;
        else if (sMode == "Luminance")           toneCurve.curveMode = ToneCurveParams::TC_MODE_LUMINANCE;
        if (pedited) pedited->toneCurve.curveMode = true; 
    }
    if (keyFile.has_key ("Exposure", "CurveMode2"))      {
        Glib::ustring sMode = keyFile.get_string ("Exposure", "CurveMode2");
        if      (sMode == "Standard")            toneCurve.curveMode2 = ToneCurveParams::TC_MODE_STD;
        else if (sMode == "FilmLike")            toneCurve.curveMode2 = ToneCurveParams::TC_MODE_FILMLIKE;
        else if (sMode == "SatAndValueBlending") toneCurve.curveMode2 = ToneCurveParams::TC_MODE_SATANDVALBLENDING;
        else if (sMode == "WeightedStd")         toneCurve.curveMode2 = ToneCurveParams::TC_MODE_WEIGHTEDSTD;
        else if (sMode == "Luminance")           toneCurve.curveMode2 = ToneCurveParams::TC_MODE_LUMINANCE;
        if (pedited) pedited->toneCurve.curveMode2 = true;
    }
    if (ppVersion>200) {
    if (keyFile.has_key ("Exposure", "Curve"))          { toneCurve.curve         = keyFile.get_double_list ("Exposure", "Curve"); if (pedited) pedited->toneCurve.curve = true; }
    if (keyFile.has_key ("Exposure", "Curve2"))         { toneCurve.curve2        = keyFile.get_double_list ("Exposure", "Curve2"); if (pedited) pedited->toneCurve.curve2 = true; }
    }
}
if (keyFile.has_group ("HLRecovery")) {
    if (keyFile.has_key ("HLRecovery", "Enabled"))  { toneCurve.hrenabled  = keyFile.get_boolean ("HLRecovery", "Enabled"); if (pedited) pedited->toneCurve.hrenabled = true; }
    if (keyFile.has_key ("HLRecovery", "Method"))   { toneCurve.method   = keyFile.get_string  ("HLRecovery", "Method"); if (pedited) pedited->toneCurve.method = true; }
}

    // load channel mixer curve
if (keyFile.has_group ("Channel Mixer")) {
    if (keyFile.has_key ("Channel Mixer", "Red") && keyFile.has_key ("Channel Mixer", "Green") && keyFile.has_key ("Channel Mixer", "Blue")) {
        if (pedited) {
            pedited->chmixer.red[0]   = pedited->chmixer.red[1]   = pedited->chmixer.red[2] = true;
            pedited->chmixer.green[0] = pedited->chmixer.green[1] = pedited->chmixer.green[2] = true;
            pedited->chmixer.blue[0]  = pedited->chmixer.blue[1]  = pedited->chmixer.blue[2] = true;
        }

        Glib::ArrayHandle<int> rmix = keyFile.get_integer_list ("Channel Mixer", "Red");
        Glib::ArrayHandle<int> gmix = keyFile.get_integer_list ("Channel Mixer", "Green");
        Glib::ArrayHandle<int> bmix = keyFile.get_integer_list ("Channel Mixer", "Blue");
        memcpy (chmixer.red, rmix.data(), 3*sizeof(int));
        memcpy (chmixer.green, gmix.data(), 3*sizeof(int));
        memcpy (chmixer.blue, bmix.data(), 3*sizeof(int));
    }
}

    // load black & white
if (keyFile.has_group ("Black & White")) {
    if (keyFile.has_key ("Black & White", "Enabled"))             { blackwhite.enabled      = keyFile.get_boolean ("Black & White", "Enabled"); if (pedited) pedited->blackwhite.enabled = true; }
    if (keyFile.has_key ("Black & White", "Method"))              { blackwhite.method       = keyFile.get_string  ("Black & White", "Method"); if (pedited) pedited->blackwhite.method = true; }

    if (keyFile.has_key ("Black & White", "Auto"))                { blackwhite.autoc        = keyFile.get_boolean ("Black & White", "Auto"); if (pedited) pedited->blackwhite.autoc = true; }
    if (keyFile.has_key ("Black & White", "ComplementaryColors")) { blackwhite.enabledcc    = keyFile.get_boolean ("Black & White", "ComplementaryColors"); if (pedited) pedited->blackwhite.enabledcc = true; }
    if (keyFile.has_key ("Black & White", "MixerRed"))            { blackwhite.mixerRed     = keyFile.get_integer ("Black & White", "MixerRed"); if (pedited) pedited->blackwhite.mixerRed = true; }
    if (keyFile.has_key ("Black & White", "MixerOrange"))         { blackwhite.mixerOrange  = keyFile.get_integer ("Black & White", "MixerOrange"); if (pedited) pedited->blackwhite.mixerOrange = true; }
    if (keyFile.has_key ("Black & White", "MixerYellow"))         { blackwhite.mixerYellow  = keyFile.get_integer ("Black & White", "MixerYellow"); if (pedited) pedited->blackwhite.mixerYellow = true; }
    if (keyFile.has_key ("Black & White", "MixerGreen"))          { blackwhite.mixerGreen   = keyFile.get_integer ("Black & White", "MixerGreen"); if (pedited) pedited->blackwhite.mixerGreen = true; }
    if (keyFile.has_key ("Black & White", "MixerCyan"))           { blackwhite.mixerCyan    = keyFile.get_integer ("Black & White", "MixerCyan"); if (pedited) pedited->blackwhite.mixerCyan = true; }
    if (keyFile.has_key ("Black & White", "MixerBlue"))           { blackwhite.mixerBlue    = keyFile.get_integer ("Black & White", "MixerBlue"); if (pedited) pedited->blackwhite.mixerBlue = true; }
    if (keyFile.has_key ("Black & White", "MixerMagenta"))        { blackwhite.mixerMagenta = keyFile.get_integer ("Black & White", "MixerMagenta"); if (pedited) pedited->blackwhite.mixerMagenta = true; }
    if (keyFile.has_key ("Black & White", "MixerPurple"))         { blackwhite.mixerPurple  = keyFile.get_integer ("Black & White", "MixerPurple"); if (pedited) pedited->blackwhite.mixerPurple = true; }
    if (keyFile.has_key ("Black & White", "GammaRed"))            { blackwhite.gammaRed     = keyFile.get_integer ("Black & White", "GammaRed"); if (pedited) pedited->blackwhite.gammaRed = true; }
    if (keyFile.has_key ("Black & White", "GammaGreen"))          { blackwhite.gammaGreen   = keyFile.get_integer ("Black & White", "GammaGreen"); if (pedited) pedited->blackwhite.gammaGreen = true; }
    if (keyFile.has_key ("Black & White", "GammaBlue"))           { blackwhite.gammaBlue    = keyFile.get_integer ("Black & White", "GammaBlue"); if (pedited) pedited->blackwhite.gammaBlue = true; }
    if (keyFile.has_key ("Black & White", "Filter"))              { blackwhite.filter       = keyFile.get_string  ("Black & White", "Filter"); if (pedited) pedited->blackwhite.filter = true; }
    if (keyFile.has_key ("Black & White", "Setting"))             { blackwhite.setting      = keyFile.get_string  ("Black & White", "Setting"); if (pedited) pedited->blackwhite.setting = true; }
    if (keyFile.has_key ("Black & White", "LuminanceCurve"))      { blackwhite.luminanceCurve = keyFile.get_double_list ("Black & White", "LuminanceCurve"); if (pedited) pedited->blackwhite.luminanceCurve = true; }
    if (keyFile.has_key ("Black & White", "BeforeCurve"))         { blackwhite.beforeCurve    = keyFile.get_double_list ("Black & White", "BeforeCurve"); if (pedited) pedited->blackwhite.beforeCurve = true; }
    if (keyFile.has_key ("Black & White", "Algorithm"))    		  { blackwhite.algo          = keyFile.get_string  ("Black & White", "Algorithm"); if (pedited) pedited->blackwhite.algo = true; }
   if (keyFile.has_key ("Black & White", "BeforeCurveMode")) {
        Glib::ustring sMode = keyFile.get_string ("Black & White", "BeforeCurveMode");
        if      (sMode == "Standard")            blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_STD_BW;
        else if (sMode == "FilmLike")            blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_FILMLIKE_BW;
        else if (sMode == "SatAndValueBlending") blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_SATANDVALBLENDING_BW;
        else if (sMode == "WeightedStd")         blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW;
        if (pedited) pedited->blackwhite.beforeCurveMode = true;
    }
    if (keyFile.has_key ("Black & White", "AfterCurve"))          { blackwhite.afterCurve     = keyFile.get_double_list ("Black & White", "AfterCurve"); if (pedited) pedited->blackwhite.afterCurve = true; }
    if (keyFile.has_key ("Black & White", "AfterCurveMode"))      {
        Glib::ustring sMode2 = keyFile.get_string ("Black & White", "AfterCurveMode");
        if      (sMode2 == "Standard")            blackwhite.afterCurveMode = BlackWhiteParams::TC_MODE_STD_BW;
        else if (sMode2 == "WeightedStd")         blackwhite.afterCurveMode = BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW;
        if (pedited) pedited->blackwhite.afterCurveMode = true;
    }
}

    // load luma curve
if (keyFile.has_group ("Luminance Curve")) {
    if (keyFile.has_key ("Luminance Curve", "Brightness"))     { labCurve.brightness   = keyFile.get_integer ("Luminance Curve", "Brightness"); if (pedited) pedited->labCurve.brightness = true; }
    if (keyFile.has_key ("Luminance Curve", "Contrast"))       { labCurve.contrast     = keyFile.get_integer ("Luminance Curve", "Contrast"); if (pedited) pedited->labCurve.contrast = true; }

    if (ppVersion < 303) {
        // transform Saturation into Chromaticity
        // if Saturation == 0, should we set BWToning on?
        if (keyFile.has_key ("Luminance Curve", "Saturation"))                { labCurve.chromaticity    = keyFile.get_integer ("Luminance Curve", "Saturation"); if (pedited) pedited->labCurve.chromaticity = true; }
        // transform AvoidColorClipping into AvoidColorShift
        if (keyFile.has_key ("Luminance Curve", "AvoidColorClipping"))        { labCurve.avoidcolorshift = keyFile.get_boolean ("Luminance Curve", "AvoidColorClipping"); if (pedited) pedited->labCurve.avoidcolorshift = true; }
    }
    else {
        if (keyFile.has_key ("Luminance Curve", "Chromaticity"))              { labCurve.chromaticity    = keyFile.get_integer ("Luminance Curve", "Chromaticity"); if (pedited) pedited->labCurve.chromaticity = true; }
        if (keyFile.has_key ("Luminance Curve", "AvoidColorShift"))           { labCurve.avoidcolorshift = keyFile.get_boolean ("Luminance Curve", "AvoidColorShift"); if (pedited) pedited->labCurve.avoidcolorshift = true; }
        if (keyFile.has_key ("Luminance Curve", "RedAndSkinTonesProtection")) { labCurve.rstprotection   = keyFile.get_double  ("Luminance Curve", "RedAndSkinTonesProtection"); if (pedited) pedited->labCurve.rstprotection = true; }
    }
    if (keyFile.has_key ("Luminance Curve", "LCredsk"))         { labCurve.lcredsk            = keyFile.get_boolean     ("Luminance Curve", "LCredsk"); if (pedited) pedited->labCurve.lcredsk = true; }
    if (ppVersion < 314)
        // Backward compatibility: If BWtoning is true, Chromaticity has to be set to -100, which will produce the same effect
        // and will enable the b&w toning mode ('a' & 'b' curves)
        if (keyFile.has_key ("Luminance Curve", "BWtoning")) {
            if ( keyFile.get_boolean     ("Luminance Curve", "BWtoning")) {
                labCurve.chromaticity = -100;
                if (pedited) pedited->labCurve.chromaticity = true;
            }
        }
    if (keyFile.has_key ("Luminance Curve", "LCurve"))          { labCurve.lcurve             = keyFile.get_double_list ("Luminance Curve", "LCurve"); if (pedited) pedited->labCurve.lcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "aCurve"))          { labCurve.acurve             = keyFile.get_double_list ("Luminance Curve", "aCurve"); if (pedited) pedited->labCurve.acurve = true; }
    if (keyFile.has_key ("Luminance Curve", "bCurve"))          { labCurve.bcurve             = keyFile.get_double_list ("Luminance Curve", "bCurve"); if (pedited) pedited->labCurve.bcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "ccCurve"))         { labCurve.cccurve            = keyFile.get_double_list ("Luminance Curve", "ccCurve"); if (pedited) pedited->labCurve.cccurve = true; }
    if (keyFile.has_key ("Luminance Curve", "chCurve"))         { labCurve.chcurve            = keyFile.get_double_list ("Luminance Curve", "chCurve"); if (pedited) pedited->labCurve.chcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "lhCurve"))         { labCurve.lhcurve            = keyFile.get_double_list ("Luminance Curve", "lhCurve"); if (pedited) pedited->labCurve.lhcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "hhCurve"))         { labCurve.hhcurve            = keyFile.get_double_list ("Luminance Curve", "hhCurve"); if (pedited) pedited->labCurve.hhcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "LcCurve"))         { labCurve.lccurve            = keyFile.get_double_list ("Luminance Curve", "LcCurve"); if (pedited) pedited->labCurve.lccurve = true; }
    if (keyFile.has_key ("Luminance Curve", "ClCurve"))         { labCurve.clcurve            = keyFile.get_double_list ("Luminance Curve", "ClCurve"); if (pedited) pedited->labCurve.clcurve = true; }

    }

    // load sharpening
if (keyFile.has_group ("Sharpening")) {
    if (keyFile.has_key ("Sharpening", "Enabled"))              { sharpening.enabled          = keyFile.get_boolean ("Sharpening", "Enabled"); if (pedited) pedited->sharpening.enabled = true; }
    if (keyFile.has_key ("Sharpening", "Radius"))               { sharpening.radius           = keyFile.get_double  ("Sharpening", "Radius"); if (pedited) pedited->sharpening.radius = true; }
    if (keyFile.has_key ("Sharpening", "Amount"))               { sharpening.amount           = keyFile.get_integer ("Sharpening", "Amount"); if (pedited) pedited->sharpening.amount = true; }
    if (keyFile.has_key ("Sharpening", "Threshold"))            {
        if (ppVersion < 302) {
            int thresh = min(keyFile.get_integer ("Sharpening", "Threshold"), 2000);
            sharpening.threshold.setValues(thresh, thresh, 2000, 2000); // TODO: 2000 is the maximum value and is taken of rtgui/sharpening.cc ; should be changed by the tool modularization
        }
        else {
            Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Sharpening", "Threshold");
            sharpening.threshold.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 2000), min(thresh.data()[3], 2000));
        }
        if (pedited) pedited->sharpening.threshold = true;
    }
    if (keyFile.has_key ("Sharpening", "OnlyEdges"))            { sharpening.edgesonly        = keyFile.get_boolean ("Sharpening", "OnlyEdges"); if (pedited) pedited->sharpening.edgesonly = true; }
    if (keyFile.has_key ("Sharpening", "EdgedetectionRadius"))  { sharpening.edges_radius     = keyFile.get_double  ("Sharpening", "EdgedetectionRadius"); if (pedited) pedited->sharpening.edges_radius = true; }
    if (keyFile.has_key ("Sharpening", "EdgeTolerance"))        { sharpening.edges_tolerance  = keyFile.get_integer ("Sharpening", "EdgeTolerance"); if (pedited) pedited->sharpening.edges_tolerance = true; }
    if (keyFile.has_key ("Sharpening", "HalocontrolEnabled"))   { sharpening.halocontrol      = keyFile.get_boolean ("Sharpening", "HalocontrolEnabled"); if (pedited) pedited->sharpening.halocontrol = true; }
    if (keyFile.has_key ("Sharpening", "HalocontrolAmount"))    { sharpening.halocontrol_amount = keyFile.get_integer ("Sharpening", "HalocontrolAmount"); if (pedited) pedited->sharpening.halocontrol_amount = true; }
    if (keyFile.has_key ("Sharpening", "Method"))               { sharpening.method           = keyFile.get_string  ("Sharpening", "Method"); if (pedited) pedited->sharpening.method = true; }
    if (keyFile.has_key ("Sharpening", "DeconvRadius"))         { sharpening.deconvradius     = keyFile.get_double  ("Sharpening", "DeconvRadius"); if (pedited) pedited->sharpening.deconvradius = true; }
    if (keyFile.has_key ("Sharpening", "DeconvAmount"))         { sharpening.deconvamount     = keyFile.get_integer ("Sharpening", "DeconvAmount"); if (pedited) pedited->sharpening.deconvamount = true; }
    if (keyFile.has_key ("Sharpening", "DeconvDamping"))        { sharpening.deconvdamping    = keyFile.get_integer ("Sharpening", "DeconvDamping"); if (pedited) pedited->sharpening.deconvdamping = true; }
    if (keyFile.has_key ("Sharpening", "DeconvIterations"))     { sharpening.deconviter       = keyFile.get_integer ("Sharpening", "DeconvIterations"); if (pedited) pedited->sharpening.deconviter = true; }
}

    // load edge sharpening
if (keyFile.has_group ("SharpenEdge")) {
    if (keyFile.has_key ("SharpenEdge", "Enabled"))             { sharpenEdge.enabled         = keyFile.get_boolean ("SharpenEdge", "Enabled"); if (pedited) pedited->sharpenEdge.enabled = true; }
    if (keyFile.has_key ("SharpenEdge", "Passes"))              { sharpenEdge.passes          = keyFile.get_integer  ("SharpenEdge", "Passes"); if (pedited) pedited->sharpenEdge.passes = true; }
    if (keyFile.has_key ("SharpenEdge", "Strength"))            { sharpenEdge.amount          = keyFile.get_double  ("SharpenEdge", "Strength"); if (pedited) pedited->sharpenEdge.amount = true; }
    if (keyFile.has_key ("SharpenEdge", "ThreeChannels"))       { sharpenEdge.threechannels   = keyFile.get_boolean ("SharpenEdge", "ThreeChannels"); if (pedited) pedited->sharpenEdge.threechannels = true; }
}

    // load micro-contrast sharpening
if (keyFile.has_group ("SharpenMicro")) {
    if (keyFile.has_key ("SharpenMicro", "Enabled"))            { sharpenMicro.enabled        = keyFile.get_boolean ("SharpenMicro", "Enabled"); if (pedited) pedited->sharpenMicro.enabled = true; }
    if (keyFile.has_key ("SharpenMicro", "Matrix"))             { sharpenMicro.matrix         = keyFile.get_boolean ("SharpenMicro", "Matrix"); if (pedited) pedited->sharpenMicro.matrix = true; }
    if (keyFile.has_key ("SharpenMicro", "Strength"))           { sharpenMicro.amount         = keyFile.get_double  ("SharpenMicro", "Strength"); if (pedited) pedited->sharpenMicro.amount = true; }
    if (keyFile.has_key ("SharpenMicro", "Uniformity"))         { sharpenMicro.uniformity     = keyFile.get_double  ("SharpenMicro", "Uniformity"); if (pedited) pedited->sharpenMicro.uniformity = true; }
}

    // load vibrance
if (keyFile.has_group ("Vibrance")) {
    if (keyFile.has_key ("Vibrance", "Enabled"))                { vibrance.enabled            = keyFile.get_boolean ("Vibrance", "Enabled"); if (pedited) pedited->vibrance.enabled = true; }
    if (keyFile.has_key ("Vibrance", "Pastels"))                { vibrance.pastels            = keyFile.get_integer ("Vibrance", "Pastels"); if (pedited) pedited->vibrance.pastels = true; }
    if (keyFile.has_key ("Vibrance", "Saturated"))              { vibrance.saturated          = keyFile.get_integer ("Vibrance", "Saturated"); if (pedited) pedited->vibrance.saturated = true; }
    if (keyFile.has_key ("Vibrance", "PSThreshold"))            {
        if (ppVersion < 302) {
            int thresh = keyFile.get_integer ("Vibrance", "PSThreshold");
            vibrance.psthreshold.setValues(thresh, thresh);
        }
        else {
            Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Vibrance", "PSThreshold");
            vibrance.psthreshold.setValues(thresh.data()[0], thresh.data()[1]);
        }
        if (pedited) pedited->vibrance.psthreshold = true;
    }
    if (keyFile.has_key ("Vibrance", "ProtectSkins"))           { vibrance.protectskins       = keyFile.get_boolean ("Vibrance", "ProtectSkins"); if (pedited) pedited->vibrance.protectskins = true; }
    if (keyFile.has_key ("Vibrance", "AvoidColorShift"))        { vibrance.avoidcolorshift    = keyFile.get_boolean ("Vibrance", "AvoidColorShift"); if (pedited) pedited->vibrance.avoidcolorshift = true; }
    if (keyFile.has_key ("Vibrance", "PastSatTog"))             { vibrance.pastsattog         = keyFile.get_boolean ("Vibrance", "PastSatTog"); if (pedited) pedited->vibrance.pastsattog = true; }
    if (keyFile.has_key ("Vibrance", "SkinTonesCurve"))        	{ vibrance.skintonescurve     = keyFile.get_double_list ("Vibrance", "SkinTonesCurve"); if (pedited) pedited->vibrance.skintonescurve = true; }
}

    // load colorBoost
/*if (keyFile.has_group ("Color Boost")) {
    if (keyFile.has_key ("Color Boost", "Amount"))              { colorBoost.amount           = keyFile.get_integer ("Color Boost", "Amount"); if (pedited) pedited->colorBoost.amount = true; }
    else {
        int a=0, b=0;
        if (keyFile.has_key ("Color Boost", "ChannelA"))        { a                           = keyFile.get_integer ("Color Boost", "ChannelA"); }
        if (keyFile.has_key ("Color Boost", "ChannelB"))        { b                           = keyFile.get_integer ("Color Boost", "ChannelB"); }
        colorBoost.amount = (a+b) / 2;
        if (pedited) pedited->colorBoost.amount = true;
    }   
    if (keyFile.has_key ("Color Boost", "AvoidColorClipping"))  { colorBoost.avoidclip               = keyFile.get_boolean ("Color Boost", "AvoidColorClipping"); if (pedited) pedited->colorBoost.avoidclip = true; }
    if (keyFile.has_key ("Color Boost", "SaturationLimiter"))   { colorBoost.enable_saturationlimiter= keyFile.get_boolean ("Color Boost", "SaturationLimiter"); if (pedited) pedited->colorBoost.enable_saturationlimiter = true; }
    if (keyFile.has_key ("Color Boost", "SaturationLimit"))     { colorBoost.saturationlimit         = keyFile.get_double  ("Color Boost", "SaturationLimit"); if (pedited) pedited->colorBoost.saturationlimit = true; }
}*/

    // load wb
if (keyFile.has_group ("White Balance")) {
    if (keyFile.has_key ("White Balance", "Setting"))     { wb.method         = keyFile.get_string ("White Balance", "Setting"); if (pedited) pedited->wb.method = true; }
    if (keyFile.has_key ("White Balance", "Temperature")) { wb.temperature    = keyFile.get_integer ("White Balance", "Temperature"); if (pedited) pedited->wb.temperature = true; }
    if (keyFile.has_key ("White Balance", "Green"))       { wb.green          = keyFile.get_double ("White Balance", "Green"); if (pedited) pedited->wb.green = true; }
    if (keyFile.has_key ("White Balance", "Equal"))       { wb.equal          = keyFile.get_double ("White Balance", "Equal"); if (pedited) pedited->wb.equal = true; }
}

    // load colorShift
/*if (keyFile.has_group ("Color Shift")) {
    if (keyFile.has_key ("Color Shift", "ChannelA")) { colorShift.a = keyFile.get_double ("Color Shift", "ChannelA"); if (pedited) pedited->colorShift.a = true; }
    if (keyFile.has_key ("Color Shift", "ChannelB")) { colorShift.b = keyFile.get_double ("Color Shift", "ChannelB"); if (pedited) pedited->colorShift.b = true; }
}*/

    // load defringe
if (keyFile.has_group ("Defringing")) {
    if (keyFile.has_key ("Defringing", "Enabled"))        { defringe.enabled   = keyFile.get_boolean ("Defringing", "Enabled"); if (pedited) pedited->defringe.enabled = true; }
    if (keyFile.has_key ("Defringing", "Radius"))         { defringe.radius    = keyFile.get_double  ("Defringing", "Radius"); if (pedited) pedited->defringe.radius = true; }
    if (keyFile.has_key ("Defringing", "Threshold"))      { defringe.threshold = (float)keyFile.get_integer ("Defringing", "Threshold"); if (pedited) pedited->defringe.threshold = true; }
    if (ppVersion < 310) {
    	 defringe.threshold =  sqrt(defringe.threshold * 33.f/5.f);
    }
    if (keyFile.has_key ("Defringing", "HueCurve"))       { defringe.huecurve  = keyFile.get_double_list ("Defringing", "HueCurve"); if (pedited) pedited->defringe.huecurve = true; }
}
    // load colorappearance
if (keyFile.has_group ("Color appearance")) {
    if (keyFile.has_key ("Color appearance", "Enabled"))       {colorappearance.enabled       = keyFile.get_boolean ("Color appearance", "Enabled"); if (pedited) pedited->colorappearance.enabled = true; }
    if (keyFile.has_key ("Color appearance", "Degree"))        {colorappearance.degree        = keyFile.get_integer ("Color appearance", "Degree"); if (pedited) pedited->colorappearance.degree = true; }
    if (keyFile.has_key ("Color appearance", "AutoDegree"))    {colorappearance.autodegree    = keyFile.get_boolean ("Color appearance", "AutoDegree"); if (pedited) pedited->colorappearance.autodegree = true; }
    if (keyFile.has_key ("Color appearance", "Surround"))      {colorappearance.surround      = keyFile.get_string  ("Color appearance", "Surround"); if (pedited) pedited->colorappearance.surround = true; }
//  if (keyFile.has_key ("Color appearance", "Background"))    {colorappearance.backgrd       = keyFile.get_integer ("Color appearance", "Background"); if (pedited) pedited->colorappearance.backgrd = true; }
    if (keyFile.has_key ("Color appearance", "AdaptLum"))      {colorappearance.adaplum       = keyFile.get_double  ("Color appearance", "AdaptLum"); if (pedited) pedited->colorappearance.adaplum = true; }
    if (keyFile.has_key ("Color appearance", "Badpixsl"))      {colorappearance.badpixsl      = keyFile.get_integer ("Color appearance", "Badpixsl"); if (pedited) pedited->colorappearance.badpixsl = true; }
    if (keyFile.has_key ("Color appearance", "Model"))         {colorappearance.wbmodel       = keyFile.get_string  ("Color appearance", "Model"); if (pedited) pedited->colorappearance.wbmodel = true; }
    if (keyFile.has_key ("Color appearance", "Algorithm"))     {colorappearance.algo          = keyFile.get_string  ("Color appearance", "Algorithm"); if (pedited) pedited->colorappearance.algo = true; }
    if (keyFile.has_key ("Color appearance", "J-Light"))       {colorappearance.jlight        = keyFile.get_double  ("Color appearance", "J-Light"); if (pedited) pedited->colorappearance.jlight = true; }
    if (keyFile.has_key ("Color appearance", "Q-Bright"))      {colorappearance.qbright       = keyFile.get_double  ("Color appearance", "Q-Bright"); if (pedited) pedited->colorappearance.qbright = true; }
    if (keyFile.has_key ("Color appearance", "C-Chroma"))      {colorappearance.chroma        = keyFile.get_double  ("Color appearance", "C-Chroma"); if (pedited) pedited->colorappearance.chroma = true; }
    if (keyFile.has_key ("Color appearance", "S-Chroma"))      {colorappearance.schroma       = keyFile.get_double  ("Color appearance", "S-Chroma"); if (pedited) pedited->colorappearance.schroma = true; }
    if (keyFile.has_key ("Color appearance", "M-Chroma"))      {colorappearance.mchroma       = keyFile.get_double  ("Color appearance", "M-Chroma"); if (pedited) pedited->colorappearance.mchroma = true; }
    if (keyFile.has_key ("Color appearance", "RSTProtection")) {colorappearance.rstprotection = keyFile.get_double  ("Color appearance", "RSTProtection"); if (pedited) pedited->colorappearance.rstprotection = true; }
    if (keyFile.has_key ("Color appearance", "J-Contrast"))    {colorappearance.contrast      = keyFile.get_double  ("Color appearance", "J-Contrast"); if (pedited) pedited->colorappearance.contrast = true; }
    if (keyFile.has_key ("Color appearance", "Q-Contrast"))    {colorappearance.qcontrast     = keyFile.get_double  ("Color appearance", "Q-Contrast"); if (pedited) pedited->colorappearance.qcontrast = true; }
    if (keyFile.has_key ("Color appearance", "H-Hue"))         {colorappearance.colorh        = keyFile.get_double  ("Color appearance", "H-Hue"); if (pedited) pedited->colorappearance.colorh = true; }
    if (keyFile.has_key ("Color appearance", "AdaptScene"))    {colorappearance.adapscen      = keyFile.get_double  ("Color appearance", "AdaptScene"); if (pedited) pedited->colorappearance.adapscen = true; }
    if (keyFile.has_key ("Color appearance", "AutoAdapscen"))  {colorappearance.autoadapscen  = keyFile.get_boolean ("Color appearance", "AutoAdapscen"); if (pedited) pedited->colorappearance.autoadapscen = true; }
    if (keyFile.has_key ("Color appearance", "SurrSource"))    {colorappearance.surrsource    = keyFile.get_boolean ("Color appearance", "SurrSource"); if (pedited) pedited->colorappearance.surrsource = true; }
    if (keyFile.has_key ("Color appearance", "Gamut"))         {colorappearance.gamut         = keyFile.get_boolean ("Color appearance", "Gamut"); if (pedited) pedited->colorappearance.gamut = true; }
//    if (keyFile.has_key ("Color appearance", "Badpix"))        {colorappearance.badpix        = keyFile.get_boolean ("Color appearance", "Badpix"); if (pedited) pedited->colorappearance.badpix = true; }
    if (keyFile.has_key ("Color appearance", "Datacie"))       {colorappearance.datacie       = keyFile.get_boolean ("Color appearance", "Datacie"); if (pedited) pedited->colorappearance.datacie = true; }
    if (keyFile.has_key ("Color appearance", "Tonecie"))       {colorappearance.tonecie       = keyFile.get_boolean ("Color appearance", "Tonecie"); if (pedited) pedited->colorappearance.tonecie = true; }
//    if (keyFile.has_key ("Color appearance", "Sharpcie"))      {colorappearance.sharpcie      = keyFile.get_boolean ("Color appearance", "Sharpcie"); if (pedited) pedited->colorappearance.sharpcie = true; }
    if (keyFile.has_key ("Color appearance", "CurveMode"))      {
        Glib::ustring sMode = keyFile.get_string ("Color appearance", "CurveMode");
        if      (sMode == "Lightness")            colorappearance.curveMode = ColorAppearanceParams::TC_MODE_LIGHT;
        else if (sMode == "Brightness")           colorappearance.curveMode = ColorAppearanceParams::TC_MODE_BRIGHT;
        if (pedited) pedited->colorappearance.curveMode = true; 
    }
    if (keyFile.has_key ("Color appearance", "CurveMode2"))      {
        Glib::ustring sMode = keyFile.get_string ("Color appearance", "CurveMode2");
        if      (sMode == "Lightness")         	colorappearance.curveMode2 = ColorAppearanceParams::TC_MODE_LIGHT;
        else if (sMode == "Brightness")         colorappearance.curveMode2 = ColorAppearanceParams::TC_MODE_BRIGHT;
        if (pedited) pedited->colorappearance.curveMode2 = true;
    }
    if (keyFile.has_key ("Color appearance", "CurveMode3"))      {
        Glib::ustring sMode = keyFile.get_string ("Color appearance", "CurveMode3");
        if      (sMode == "Chroma")         	colorappearance.curveMode3 = ColorAppearanceParams::TC_MODE_CHROMA;
        else if (sMode == "Saturation")         colorappearance.curveMode3 = ColorAppearanceParams::TC_MODE_SATUR;
        else if (sMode == "Colorfullness")      colorappearance.curveMode3 = ColorAppearanceParams::TC_MODE_COLORF;

        if (pedited) pedited->colorappearance.curveMode3 = true;
    }
	
    if (ppVersion>200) {
    if (keyFile.has_key ("Color appearance", "Curve"))          { colorappearance.curve         = keyFile.get_double_list ("Color appearance", "Curve"); if (pedited) pedited->colorappearance.curve = true; }
    if (keyFile.has_key ("Color appearance", "Curve2"))         { colorappearance.curve2        = keyFile.get_double_list ("Color appearance", "Curve2"); if (pedited) pedited->colorappearance.curve2 = true; }
    if (keyFile.has_key ("Color appearance", "Curve3"))         { colorappearance.curve3        = keyFile.get_double_list ("Color appearance", "Curve3"); if (pedited) pedited->colorappearance.curve3 = true; }
	}
	
	}

    // load impulseDenoise
if (keyFile.has_group ("Impulse Denoising")) {
    if (keyFile.has_key ("Impulse Denoising", "Enabled"))   { impulseDenoise.enabled = keyFile.get_boolean ("Impulse Denoising", "Enabled"); if (pedited) pedited->impulseDenoise.enabled = true; }
    if (keyFile.has_key ("Impulse Denoising", "Threshold")) { impulseDenoise.thresh  = keyFile.get_integer ("Impulse Denoising", "Threshold"); if (pedited) pedited->impulseDenoise.thresh = true; }
}

    // load dirpyrDenoise
if (keyFile.has_group ("Directional Pyramid Denoising")) {//TODO: No longer an accurate description for FT denoise
    if (keyFile.has_key ("Directional Pyramid Denoising", "Enabled"))    { dirpyrDenoise.enabled = keyFile.get_boolean ("Directional Pyramid Denoising", "Enabled"); if (pedited) pedited->dirpyrDenoise.enabled = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Enhance"))    { dirpyrDenoise.enhance = keyFile.get_boolean ("Directional Pyramid Denoising", "Enhance"); if (pedited) pedited->dirpyrDenoise.enhance = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Median"))    { dirpyrDenoise.median = keyFile.get_boolean ("Directional Pyramid Denoising", "Median"); if (pedited) pedited->dirpyrDenoise.median = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Auto"))    { dirpyrDenoise.autochroma = keyFile.get_boolean ("Directional Pyramid Denoising", "Auto"); if (pedited) pedited->dirpyrDenoise.autochroma = true; }
 //   if (keyFile.has_key ("Directional Pyramid Denoising", "Perform"))    { dirpyrDenoise.perform = keyFile.get_boolean ("Directional Pyramid Denoising", "Perform"); if (pedited) pedited->dirpyrDenoise.perform = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Luma"))       { dirpyrDenoise.luma    = keyFile.get_double ("Directional Pyramid Denoising", "Luma"); if (pedited) pedited->dirpyrDenoise.luma = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Ldetail"))    { dirpyrDenoise.Ldetail = keyFile.get_double ("Directional Pyramid Denoising", "Ldetail"); if (pedited) pedited->dirpyrDenoise.Ldetail = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Chroma"))     { dirpyrDenoise.chroma  = keyFile.get_double ("Directional Pyramid Denoising", "Chroma"); if (pedited) pedited->dirpyrDenoise.chroma = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Method"))     {dirpyrDenoise.dmethod  = keyFile.get_string  ("Directional Pyramid Denoising", "Method"); if (pedited) pedited->dirpyrDenoise.dmethod = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "LMethod"))     {dirpyrDenoise.Lmethod  = keyFile.get_string  ("Directional Pyramid Denoising", "LMethod"); if (pedited) pedited->dirpyrDenoise.Lmethod = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "CMethod"))     {dirpyrDenoise.Cmethod  = keyFile.get_string  ("Directional Pyramid Denoising", "CMethod"); if (pedited) pedited->dirpyrDenoise.Cmethod = true; }
    // never load 'auto chroma preview mode' from pp3
	if(dirpyrDenoise.Cmethod=="PRE")
		dirpyrDenoise.Cmethod = "MAN";
	if (keyFile.has_key ("Directional Pyramid Denoising", "C2Method"))     {dirpyrDenoise.C2method  = keyFile.get_string  ("Directional Pyramid Denoising", "C2Method"); if (pedited) pedited->dirpyrDenoise.C2method = true; }
	if(dirpyrDenoise.C2method=="PREV")
		dirpyrDenoise.C2method = "MANU";
    if (keyFile.has_key ("Directional Pyramid Denoising", "SMethod"))     {dirpyrDenoise.smethod  = keyFile.get_string  ("Directional Pyramid Denoising", "SMethod"); if (pedited) pedited->dirpyrDenoise.smethod = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "MedMethod"))     {dirpyrDenoise.medmethod  = keyFile.get_string  ("Directional Pyramid Denoising", "MedMethod"); if (pedited) pedited->dirpyrDenoise.medmethod = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "MethodMed"))     {dirpyrDenoise.methodmed  = keyFile.get_string  ("Directional Pyramid Denoising", "MethodMed"); if (pedited) pedited->dirpyrDenoise.methodmed = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "RGBMethod"))     {dirpyrDenoise.rgbmethod  = keyFile.get_string  ("Directional Pyramid Denoising", "RGBMethod"); if (pedited) pedited->dirpyrDenoise.rgbmethod = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "LCurve"))          {dirpyrDenoise.lcurve             = keyFile.get_double_list ("Directional Pyramid Denoising", "LCurve"); if (pedited) pedited->dirpyrDenoise.lcurve = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "CCCurve"))          {dirpyrDenoise.cccurve             = keyFile.get_double_list ("Directional Pyramid Denoising", "CCCurve"); if (pedited) pedited->dirpyrDenoise.cccurve = true; }

    if (keyFile.has_key ("Directional Pyramid Denoising", "Redchro"))    { dirpyrDenoise.redchro  = keyFile.get_double ("Directional Pyramid Denoising", "Redchro"); if (pedited) pedited->dirpyrDenoise.redchro = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Bluechro"))   { dirpyrDenoise.bluechro  = keyFile.get_double ("Directional Pyramid Denoising", "Bluechro"); if (pedited) pedited->dirpyrDenoise.bluechro = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Gamma"))      { dirpyrDenoise.gamma   = keyFile.get_double ("Directional Pyramid Denoising", "Gamma"); if (pedited) pedited->dirpyrDenoise.gamma = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Passes"))      { dirpyrDenoise.passes   = keyFile.get_integer ("Directional Pyramid Denoising", "Passes"); if (pedited) pedited->dirpyrDenoise.passes = true; }
}

    //Load EPD.
if (keyFile.has_group ("EPD")) {
    if(keyFile.has_key("EPD", "Enabled"))             { epd.enabled = keyFile.get_boolean ("EPD", "Enabled"); if (pedited) pedited->epd.enabled = true; }
    if(keyFile.has_key("EPD", "Strength"))            { epd.strength = keyFile.get_double ("EPD", "Strength"); if (pedited) pedited->epd.strength = true; }
    if(keyFile.has_key("EPD", "Gamma"))         	   { epd.gamma = keyFile.get_double ("EPD", "Gamma"); if (pedited) pedited->epd.gamma = true; }
    if(keyFile.has_key("EPD", "EdgeStopping"))        { epd.edgeStopping = keyFile.get_double ("EPD", "EdgeStopping"); if (pedited) pedited->epd.edgeStopping = true; }
    if(keyFile.has_key("EPD", "Scale"))               { epd.scale = keyFile.get_double ("EPD", "Scale"); if (pedited) pedited->epd.scale = true; }
    if(keyFile.has_key("EPD", "ReweightingIterates")) { epd.reweightingIterates = keyFile.get_integer ("EPD", "ReweightingIterates"); if (pedited) pedited->epd.reweightingIterates = true; }
}

    // load lumaDenoise
/*if (keyFile.has_group ("Luminance Denoising")) {
    if (keyFile.has_key ("Luminance Denoising", "Enabled"))        { lumaDenoise.enabled       = keyFile.get_boolean ("Luminance Denoising", "Enabled"); if (pedited) pedited->lumaDenoise.enabled = true; }
    if (keyFile.has_key ("Luminance Denoising", "Radius"))         { lumaDenoise.radius        = keyFile.get_double  ("Luminance Denoising", "Radius"); if (pedited) pedited->lumaDenoise.radius = true; }
    if (keyFile.has_key ("Luminance Denoising", "EdgeTolerance"))  { lumaDenoise.edgetolerance = keyFile.get_integer ("Luminance Denoising", "EdgeTolerance"); if (pedited) pedited->lumaDenoise.edgetolerance = true; }
}*/

    // load colorDenoise
/*if (keyFile.has_group ("Chrominance Denoising")) {
    if (keyFile.has_key ("Chrominance Denoising", "Enabled"))      { colorDenoise.enabled       = keyFile.get_boolean 	("Chrominance Denoising", "Enabled"); if (pedited) pedited->colorDenoise.enabled = true; }
    // WARNING: radius doesn't exist anymore; is there any compatibility issue that require to keep the following line?
    if (keyFile.has_key ("Chrominance Denoising", "Radius"))       { colorDenoise.amount        = 10*keyFile.get_double ("Chrominance Denoising", "Radius"); }
    if (keyFile.has_key ("Chrominance Denoising", "Amount"))       { colorDenoise.amount        = keyFile.get_integer  	("Chrominance Denoising", "Amount"); if (pedited) pedited->colorDenoise.amount = true; }
}*/

    // load sh
if (keyFile.has_group ("Shadows & Highlights")) {
    if (keyFile.has_key ("Shadows & Highlights", "Enabled"))               { sh.enabled       = keyFile.get_boolean ("Shadows & Highlights", "Enabled"); if (pedited) pedited->sh.enabled = true; }
    if (keyFile.has_key ("Shadows & Highlights", "HighQuality"))           { sh.hq            = keyFile.get_boolean ("Shadows & Highlights", "HighQuality"); if (pedited) pedited->sh.hq = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Highlights"))            { sh.highlights    = keyFile.get_integer ("Shadows & Highlights", "Highlights"); if (pedited) pedited->sh.highlights = true; }
    if (keyFile.has_key ("Shadows & Highlights", "HighlightTonalWidth"))   { sh.htonalwidth   = keyFile.get_integer ("Shadows & Highlights", "HighlightTonalWidth"); if (pedited) pedited->sh.htonalwidth = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Shadows"))               { sh.shadows       = keyFile.get_integer ("Shadows & Highlights", "Shadows"); if (pedited) pedited->sh.shadows = true; }
    if (keyFile.has_key ("Shadows & Highlights", "ShadowTonalWidth"))      { sh.stonalwidth   = keyFile.get_integer ("Shadows & Highlights", "ShadowTonalWidth"); if (pedited) pedited->sh.stonalwidth = true; }
    if (keyFile.has_key ("Shadows & Highlights", "LocalContrast"))         { sh.localcontrast = keyFile.get_integer ("Shadows & Highlights", "LocalContrast"); if (pedited) pedited->sh.localcontrast = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Radius"))                { sh.radius        = keyFile.get_integer ("Shadows & Highlights", "Radius"); if (pedited) pedited->sh.radius = true; }
}
    
    // load crop
if (keyFile.has_group ("Crop")) {
    if (keyFile.has_key ("Crop", "Enabled"))    { crop.enabled    = keyFile.get_boolean ("Crop", "Enabled"); if (pedited) pedited->crop.enabled = true; }
    if (keyFile.has_key ("Crop", "X"))          { crop.x          = keyFile.get_integer ("Crop", "X"); if (pedited) pedited->crop.x = true; }
    if (keyFile.has_key ("Crop", "Y"))          { crop.y          = keyFile.get_integer ("Crop", "Y"); if (pedited) pedited->crop.y = true; }
    if (keyFile.has_key ("Crop", "W"))          { crop.w          = keyFile.get_integer ("Crop", "W"); if (pedited) pedited->crop.w = true; }
    if (keyFile.has_key ("Crop", "H"))          { crop.h          = keyFile.get_integer ("Crop", "H"); if (pedited) pedited->crop.h = true; }
    if (keyFile.has_key ("Crop", "FixedRatio")) { crop.fixratio   = keyFile.get_boolean ("Crop", "FixedRatio"); if (pedited) pedited->crop.fixratio = true; }
    if (keyFile.has_key ("Crop", "Ratio")) {
        crop.ratio      = keyFile.get_string  ("Crop", "Ratio");
        if (pedited) pedited->crop.ratio = true;
        //backwards compatibility for crop.ratio
        if (crop.ratio=="DIN")    crop.ratio = "1.414 - DIN EN ISO 216";
        if (crop.ratio=="8.5:11") crop.ratio = "8.5:11 - US Letter";
        if (crop.ratio=="11:17")  crop.ratio = "11:17 - Tabloid";
    }
    if (keyFile.has_key ("Crop", "Orientation"))  { crop.orientation= keyFile.get_string  ("Crop", "Orientation"); if (pedited) pedited->crop.orientation = true; }
    if (keyFile.has_key ("Crop", "Guide"))        { crop.guide      = keyFile.get_string  ("Crop", "Guide"); if (pedited) pedited->crop.guide = true; }
}

    // load coarse
if (keyFile.has_group ("Coarse Transformation")) {
    if (keyFile.has_key ("Coarse Transformation", "Rotate"))          { coarse.rotate = keyFile.get_integer ("Coarse Transformation", "Rotate"); if (pedited) pedited->coarse.rotate = true; }
    if (keyFile.has_key ("Coarse Transformation", "HorizontalFlip"))  { coarse.hflip  = keyFile.get_boolean ("Coarse Transformation", "HorizontalFlip"); if (pedited) pedited->coarse.hflip = true; }
    if (keyFile.has_key ("Coarse Transformation", "VerticalFlip"))    { coarse.vflip  = keyFile.get_boolean ("Coarse Transformation", "VerticalFlip"); if (pedited) pedited->coarse.vflip = true; }
}

    // load rotate
if (keyFile.has_group ("Rotation")) {
    if (keyFile.has_key ("Rotation", "Degree"))   { rotate.degree = keyFile.get_double ("Rotation", "Degree"); if (pedited) pedited->rotate.degree = true; }
}
    // load commonTrans
if (keyFile.has_group ("Common Properties for Transformations")) {
    if (keyFile.has_key ("Common Properties for Transformations", "AutoFill"))   { commonTrans.autofill = keyFile.get_boolean ("Common Properties for Transformations", "AutoFill"); if (pedited) pedited->commonTrans.autofill = true; }
}

    // load distortion
if (keyFile.has_group ("Distortion")) {
    if (keyFile.has_key ("Distortion", "Amount"))     { distortion.amount     = keyFile.get_double  ("Distortion", "Amount"); if (pedited) pedited->distortion.amount = true; }
}

    // lens profile
if (keyFile.has_group ("LensProfile")) {
    if (keyFile.has_key ("LensProfile", "LCPFile")) { lensProf.lcpFile = expandRelativePath(fname, "", keyFile.get_string ("LensProfile", "LCPFile")); if (pedited) pedited->lensProf.lcpFile = true; }
    if (keyFile.has_key ("LensProfile", "UseDistortion")) { lensProf.useDist = keyFile.get_boolean ("LensProfile", "UseDistortion"); if (pedited) pedited->lensProf.useDist = true; }
    if (keyFile.has_key ("LensProfile", "UseVignette")) { lensProf.useVign = keyFile.get_boolean ("LensProfile", "UseVignette"); if (pedited) pedited->lensProf.useVign = true; }
    if (keyFile.has_key ("LensProfile", "UseCA")) { lensProf.useCA = keyFile.get_boolean ("LensProfile", "UseCA"); if (pedited) pedited->lensProf.useCA = true; }
}
    
    // load perspective correction
if (keyFile.has_group ("Perspective")) {
    if (keyFile.has_key ("Perspective", "Horizontal"))  { perspective.horizontal = keyFile.get_double ("Perspective", "Horizontal"); if (pedited) pedited->perspective.horizontal = true; }
    if (keyFile.has_key ("Perspective", "Vertical"))    { perspective.vertical   = keyFile.get_double ("Perspective", "Vertical"); if (pedited) pedited->perspective.vertical = true; }
}

    // load gradient
if (keyFile.has_group ("Gradient")) {
    if (keyFile.has_key ("Gradient", "Enabled"))  { gradient.enabled  = keyFile.get_boolean ("Gradient", "Enabled"); if (pedited) pedited->gradient.enabled = true; }
    if (keyFile.has_key ("Gradient", "Degree"))   { gradient.degree   = keyFile.get_double  ("Gradient", "Degree");  if (pedited) pedited->gradient.degree = true; }
    if (keyFile.has_key ("Gradient", "Feather"))  { gradient.feather  = keyFile.get_integer ("Gradient", "Feather"); if (pedited) pedited->gradient.feather = true; }
    if (keyFile.has_key ("Gradient", "Strength")) { gradient.strength = keyFile.get_double  ("Gradient", "Strength");if (pedited) pedited->gradient.strength = true; }
    if (keyFile.has_key ("Gradient", "CenterX"))  { gradient.centerX  = keyFile.get_integer ("Gradient", "CenterX"); if (pedited) pedited->gradient.centerX = true; }
    if (keyFile.has_key ("Gradient", "CenterY"))  { gradient.centerY  = keyFile.get_integer ("Gradient", "CenterY"); if (pedited) pedited->gradient.centerY = true; }
}

if (keyFile.has_group ("PCVignette")) {
    if (keyFile.has_key ("PCVignette", "Enabled"))  { pcvignette.enabled  = keyFile.get_boolean ("PCVignette", "Enabled"); if (pedited) pedited->pcvignette.enabled = true; }
    if (keyFile.has_key ("PCVignette", "Strength")) { pcvignette.strength = keyFile.get_double  ("PCVignette", "Strength");if (pedited) pedited->pcvignette.strength = true; }
    if (keyFile.has_key ("PCVignette", "Feather"))  { pcvignette.feather  = keyFile.get_integer ("PCVignette", "Feather"); if (pedited) pedited->pcvignette.feather = true; }
    if (keyFile.has_key ("PCVignette", "Roundness"))  { pcvignette.roundness  = keyFile.get_integer ("PCVignette", "Roundness"); if (pedited) pedited->pcvignette.roundness = true; }
}

// load c/a correction
if (keyFile.has_group ("CACorrection")) {
    if (keyFile.has_key ("CACorrection", "Red"))  { cacorrection.red  = keyFile.get_double ("CACorrection", "Red"); if (pedited) pedited->cacorrection.red = true; }
    if (keyFile.has_key ("CACorrection", "Blue")) { cacorrection.blue = keyFile.get_double ("CACorrection", "Blue"); if (pedited) pedited->cacorrection.blue = true; }
}

    // load vignetting correction
if (keyFile.has_group ("Vignetting Correction")) {
    if (keyFile.has_key ("Vignetting Correction", "Amount"))   { vignetting.amount = keyFile.get_integer ("Vignetting Correction", "Amount"); if (pedited) pedited->vignetting.amount = true; }
    if (keyFile.has_key ("Vignetting Correction", "Radius"))   { vignetting.radius = keyFile.get_integer ("Vignetting Correction", "Radius"); if (pedited) pedited->vignetting.radius = true; }
    if (keyFile.has_key ("Vignetting Correction", "Strength")) { vignetting.strength = keyFile.get_integer ("Vignetting Correction", "Strength"); if (pedited) pedited->vignetting.strength = true; }
    if (keyFile.has_key ("Vignetting Correction", "CenterX"))  { vignetting.centerX = keyFile.get_integer ("Vignetting Correction", "CenterX"); if (pedited) pedited->vignetting.centerX = true; }
    if (keyFile.has_key ("Vignetting Correction", "CenterY"))  { vignetting.centerY = keyFile.get_integer ("Vignetting Correction", "CenterY"); if (pedited) pedited->vignetting.centerY = true; }
}

    // load resize settings
if (keyFile.has_group ("Resize")) {
    if (keyFile.has_key ("Resize", "Enabled"))       { resize.enabled   = keyFile.get_boolean ("Resize", "Enabled"); if (pedited) pedited->resize.enabled = true; }
    if (keyFile.has_key ("Resize", "Scale"))         { resize.scale     = keyFile.get_double ("Resize", "Scale"); if (pedited) pedited->resize.scale = true; }
    if (keyFile.has_key ("Resize", "AppliesTo"))     { resize.appliesTo = keyFile.get_string ("Resize", "AppliesTo"); if (pedited) pedited->resize.appliesTo = true; }
    if (keyFile.has_key ("Resize", "Method"))        { resize.method    = keyFile.get_string ("Resize", "Method"); if (pedited) pedited->resize.method = true; }
    if (keyFile.has_key ("Resize", "DataSpecified")) { resize.dataspec  = keyFile.get_integer ("Resize", "DataSpecified"); if (pedited) pedited->resize.dataspec = true; }
    if (keyFile.has_key ("Resize", "Width"))         { resize.width     = keyFile.get_integer ("Resize", "Width"); if (pedited) pedited->resize.width = true; }
    if (keyFile.has_key ("Resize", "Height"))        { resize.height    = keyFile.get_integer ("Resize", "Height"); if (pedited) pedited->resize.height = true; }
}

    // load post resize sharpening
if (keyFile.has_group ("PostResizeSharpening")) {
    if (keyFile.has_key ("PostResizeSharpening", "Enabled"))              { prsharpening.enabled          = keyFile.get_boolean ("PostResizeSharpening", "Enabled"); if (pedited) pedited->prsharpening.enabled = true; }
    if (keyFile.has_key ("PostResizeSharpening", "Radius"))               { prsharpening.radius           = keyFile.get_double  ("PostResizeSharpening", "Radius"); if (pedited) pedited->prsharpening.radius = true; }
    if (keyFile.has_key ("PostResizeSharpening", "Amount"))               { prsharpening.amount           = keyFile.get_integer ("PostResizeSharpening", "Amount"); if (pedited) pedited->prsharpening.amount = true; }
    if (keyFile.has_key ("PostResizeSharpening", "Threshold"))            {
        if (ppVersion < 302) {
            int thresh = min(keyFile.get_integer ("PostResizeSharpening", "Threshold"), 2000);
            prsharpening.threshold.setValues(thresh, thresh, 2000, 2000); // TODO: 2000 is the maximum value and is taken of rtgui/sharpening.cc ; should be changed by the tool modularization
        }
        else {
            Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("PostResizeSharpening", "Threshold");
            prsharpening.threshold.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 2000), min(thresh.data()[3], 2000));
        }
        if (pedited) pedited->prsharpening.threshold = true;
    }
    if (keyFile.has_key ("PostResizeSharpening", "OnlyEdges"))            { prsharpening.edgesonly        = keyFile.get_boolean ("PostResizeSharpening", "OnlyEdges"); if (pedited) pedited->prsharpening.edgesonly = true; }
    if (keyFile.has_key ("PostResizeSharpening", "EdgedetectionRadius"))  { prsharpening.edges_radius     = keyFile.get_double  ("PostResizeSharpening", "EdgedetectionRadius"); if (pedited) pedited->prsharpening.edges_radius = true; }
    if (keyFile.has_key ("PostResizeSharpening", "EdgeTolerance"))        { prsharpening.edges_tolerance  = keyFile.get_integer ("PostResizeSharpening", "EdgeTolerance"); if (pedited) pedited->prsharpening.edges_tolerance = true; }
    if (keyFile.has_key ("PostResizeSharpening", "HalocontrolEnabled"))   { prsharpening.halocontrol      = keyFile.get_boolean ("PostResizeSharpening", "HalocontrolEnabled"); if (pedited) pedited->prsharpening.halocontrol = true; }
    if (keyFile.has_key ("PostResizeSharpening", "HalocontrolAmount"))    { prsharpening.halocontrol_amount = keyFile.get_integer ("PostResizeSharpening", "HalocontrolAmount"); if (pedited) pedited->prsharpening.halocontrol_amount = true; }
    if (keyFile.has_key ("PostResizeSharpening", "Method"))               { prsharpening.method           = keyFile.get_string  ("PostResizeSharpening", "Method"); if (pedited) pedited->prsharpening.method = true; }
    if (keyFile.has_key ("PostResizeSharpening", "DeconvRadius"))         { prsharpening.deconvradius     = keyFile.get_double  ("PostResizeSharpening", "DeconvRadius"); if (pedited) pedited->prsharpening.deconvradius = true; }
    if (keyFile.has_key ("PostResizeSharpening", "DeconvAmount"))         { prsharpening.deconvamount     = keyFile.get_integer ("PostResizeSharpening", "DeconvAmount"); if (pedited) pedited->prsharpening.deconvamount = true; }
    if (keyFile.has_key ("PostResizeSharpening", "DeconvDamping"))        { prsharpening.deconvdamping    = keyFile.get_integer ("PostResizeSharpening", "DeconvDamping"); if (pedited) pedited->prsharpening.deconvdamping = true; }
    if (keyFile.has_key ("PostResizeSharpening", "DeconvIterations"))     { prsharpening.deconviter       = keyFile.get_integer ("PostResizeSharpening", "DeconvIterations"); if (pedited) pedited->prsharpening.deconviter = true; }
}

    // load color management settings
if (keyFile.has_group ("Color Management")) {
    if (keyFile.has_key ("Color Management", "InputProfile"))   { icm.input          = expandRelativePath(fname, "file:", keyFile.get_string ("Color Management", "InputProfile")); if (pedited) pedited->icm.input = true; }
    if (keyFile.has_key ("Color Management", "ToneCurve"))      { icm.toneCurve      = keyFile.get_boolean ("Color Management", "ToneCurve"); if (pedited) pedited->icm.toneCurve = true; }
    if (keyFile.has_key ("Color Management", "ApplyLookTable"))      { icm.applyLookTable      = keyFile.get_boolean ("Color Management", "ApplyLookTable"); if (pedited) pedited->icm.applyLookTable = true; }
    if (keyFile.has_key ("Color Management", "ApplyBaselineExposureOffset"))      { icm.applyBaselineExposureOffset      = keyFile.get_boolean ("Color Management", "ApplyBaselineExposureOffset"); if (pedited) pedited->icm.applyBaselineExposureOffset = true; }
    if (keyFile.has_key ("Color Management", "ApplyHueSatMap"))      { icm.applyHueSatMap      = keyFile.get_boolean ("Color Management", "ApplyHueSatMap"); if (pedited) pedited->icm.applyHueSatMap = true; }
    if (keyFile.has_key ("Color Management", "BlendCMSMatrix")) { icm.blendCMSMatrix = keyFile.get_boolean ("Color Management", "BlendCMSMatrix"); if (pedited) pedited->icm.blendCMSMatrix = true; }
    if (keyFile.has_key ("Color Management", "DCPIlluminant"))  { icm.dcpIlluminant  = keyFile.get_integer ("Color Management", "DCPIlluminant"); if (pedited) pedited->icm.dcpIlluminant = true; }
    if (keyFile.has_key ("Color Management", "WorkingProfile")) { icm.working        = keyFile.get_string ("Color Management", "WorkingProfile"); if (pedited) pedited->icm.working = true; }
    if (keyFile.has_key ("Color Management", "OutputProfile"))  { icm.output         = keyFile.get_string ("Color Management", "OutputProfile"); if (pedited) pedited->icm.output = true; }
    if (keyFile.has_key ("Color Management", "Gammafree"))      { icm.gamma          = keyFile.get_string ("Color Management", "Gammafree"); if (pedited) pedited->icm.gamfree = true; }
    if (keyFile.has_key ("Color Management", "Freegamma"))      { icm.freegamma      = keyFile.get_boolean ("Color Management", "Freegamma"); if (pedited) pedited->icm.freegamma = true; }
    if (keyFile.has_key ("Color Management", "GammaValue"))     { icm.gampos         = keyFile.get_double ("Color Management", "GammaValue"); if (pedited) pedited->icm.gampos = true; }
    if (keyFile.has_key ("Color Management", "GammaSlope"))     { icm.slpos          = keyFile.get_double ("Color Management", "GammaSlope"); if (pedited) pedited->icm.slpos = true; }

}
    // load wavelet wavelet parameters
if (keyFile.has_group ("Wavelet")) {
    if (keyFile.has_key ("Wavelet", "Enabled"))   { wavelet.enabled = keyFile.get_boolean ("Wavelet", "Enabled"); if (pedited) pedited->wavelet.enabled = true; }
	if (keyFile.has_key ("Wavelet", "Strength"))   { wavelet.strength = keyFile.get_integer ("Wavelet", "Strength"); if (pedited) pedited->wavelet.strength = true; }
	if (keyFile.has_key ("Wavelet", "Balance"))   { wavelet.balance = keyFile.get_integer ("Wavelet", "Balance"); if (pedited) pedited->wavelet.balance = true; }
	if (keyFile.has_key ("Wavelet", "Iter"))   { wavelet.iter = keyFile.get_integer ("Wavelet", "Iter"); if (pedited) pedited->wavelet.iter = true; }
    if (keyFile.has_key ("Wavelet", "Median")) {wavelet.median = keyFile.get_boolean ("Wavelet", "Median");if (pedited) pedited->wavelet.median = true;}
    if (keyFile.has_key ("Wavelet", "Medianlev")) {wavelet.medianlev = keyFile.get_boolean ("Wavelet", "Medianlev");if (pedited) pedited->wavelet.medianlev = true;}
    if (keyFile.has_key ("Wavelet", "Linkedg")) {wavelet.linkedg = keyFile.get_boolean ("Wavelet", "Linkedg");if (pedited) pedited->wavelet.linkedg = true;}
    if (keyFile.has_key ("Wavelet", "CBenab")) {wavelet.cbenab = keyFile.get_boolean ("Wavelet", "CBenab");if (pedited) pedited->wavelet.cbenab = true;}
	if (keyFile.has_key ("Wavelet", "CBgreenhigh"))   { wavelet.greenhigh = keyFile.get_integer ("Wavelet", "CBgreenhigh"); if (pedited) pedited->wavelet.greenhigh = true; }
	if (keyFile.has_key ("Wavelet", "CBgreenmed"))   { wavelet.greenmed = keyFile.get_integer ("Wavelet", "CBgreenmed"); if (pedited) pedited->wavelet.greenmed = true; }
	if (keyFile.has_key ("Wavelet", "CBgreenlow"))   { wavelet.greenlow = keyFile.get_integer ("Wavelet", "CBgreenlow"); if (pedited) pedited->wavelet.greenlow = true; }
	if (keyFile.has_key ("Wavelet", "CBbluehigh"))   { wavelet.bluehigh = keyFile.get_integer ("Wavelet", "CBbluehigh"); if (pedited) pedited->wavelet.bluehigh = true; }
	if (keyFile.has_key ("Wavelet", "CBbluemed"))   { wavelet.bluemed = keyFile.get_integer ("Wavelet", "CBbluemed"); if (pedited) pedited->wavelet.bluemed = true; }
	if (keyFile.has_key ("Wavelet", "CBbluelow"))   { wavelet.bluelow = keyFile.get_integer ("Wavelet", "CBbluelow"); if (pedited) pedited->wavelet.bluelow = true; }
 //   if (keyFile.has_key ("Wavelet", "Edgreinf")) {wavelet.edgreinf = keyFile.get_boolean ("Wavelet", "Edgreinf");if (pedited) pedited->wavelet.edgreinf = true;}
    if (keyFile.has_key ("Wavelet", "Lipst")) {wavelet.lipst = keyFile.get_boolean ("Wavelet", "Lipst");if (pedited) pedited->wavelet.lipst = true;}
    if (keyFile.has_key ("Wavelet", "AvoidColorShift")) {wavelet.avoid = keyFile.get_boolean ("Wavelet", "AvoidColorShift");if (pedited) pedited->wavelet.avoid = true;}
    if (keyFile.has_key ("Wavelet", "TMr")) {wavelet.tmr = keyFile.get_boolean ("Wavelet", "TMr");if (pedited) pedited->wavelet.tmr = true;}
    if (keyFile.has_key ("Wavelet", "LevMethod"))     {wavelet.Lmethod  = keyFile.get_string  ("Wavelet", "LevMethod"); if (pedited) pedited->wavelet.Lmethod = true; }
    if (keyFile.has_key ("Wavelet", "ChoiceLevMethod"))     {wavelet.CLmethod  = keyFile.get_string  ("Wavelet", "ChoiceLevMethod"); if (pedited) pedited->wavelet.CLmethod = true; }
    if (keyFile.has_key ("Wavelet", "BackMethod"))     {wavelet.Backmethod  = keyFile.get_string  ("Wavelet", "BackMethod"); if (pedited) pedited->wavelet.Backmethod = true; }
    if (keyFile.has_key ("Wavelet", "TilesMethod"))     {wavelet.Tilesmethod  = keyFile.get_string  ("Wavelet", "TilesMethod"); if (pedited) pedited->wavelet.Tilesmethod = true; }
    if (keyFile.has_key ("Wavelet", "DaubMethod"))     {wavelet.daubcoeffmethod  = keyFile.get_string  ("Wavelet", "DaubMethod"); if (pedited) pedited->wavelet.daubcoeffmethod = true; }
    if (keyFile.has_key ("Wavelet", "CHromaMethod"))     {wavelet.CHmethod  = keyFile.get_string  ("Wavelet", "CHromaMethod"); if (pedited) pedited->wavelet.CHmethod = true; }
    if (keyFile.has_key ("Wavelet", "Medgreinf"))     {wavelet.Medgreinf  = keyFile.get_string  ("Wavelet", "Medgreinf"); if (pedited) pedited->wavelet.Medgreinf = true; }
    if (keyFile.has_key ("Wavelet", "CHSLromaMethod"))     {wavelet.CHSLmethod  = keyFile.get_string  ("Wavelet", "CHSLromaMethod"); if (pedited) pedited->wavelet.CHSLmethod = true; }
    if (keyFile.has_key ("Wavelet", "EDMethod"))     {wavelet.EDmethod  = keyFile.get_string  ("Wavelet", "EDMethod"); if (pedited) pedited->wavelet.EDmethod = true; }
    if (keyFile.has_key ("Wavelet", "BAMethod"))     {wavelet.BAmethod  = keyFile.get_string  ("Wavelet", "BAMethod"); if (pedited) pedited->wavelet.BAmethod = true; }
    if (keyFile.has_key ("Wavelet", "TMMethod"))     {wavelet.TMmethod  = keyFile.get_string  ("Wavelet", "TMMethod"); if (pedited) pedited->wavelet.TMmethod = true; }
    if (keyFile.has_key ("Wavelet", "HSMethod"))     {wavelet.HSmethod  = keyFile.get_string  ("Wavelet", "HSMethod"); if (pedited) pedited->wavelet.HSmethod = true; }
    if (keyFile.has_key ("Wavelet", "DirMethod"))     {wavelet.Dirmethod  = keyFile.get_string  ("Wavelet", "DirMethod"); if (pedited) pedited->wavelet.Dirmethod = true; }
    if (keyFile.has_key ("Wavelet", "ResidualcontShadow"))     {wavelet.rescon  = keyFile.get_integer  ("Wavelet", "ResidualcontShadow"); if (pedited) pedited->wavelet.rescon = true; }
	if (keyFile.has_key ("Wavelet", "ResidualcontHighlight"))     {wavelet.resconH  = keyFile.get_integer  ("Wavelet", "ResidualcontHighlight"); if (pedited) pedited->wavelet.resconH = true; }   
	if (keyFile.has_key ("Wavelet", "Residualchroma"))     {wavelet.reschro  = keyFile.get_integer  ("Wavelet", "Residualchroma"); if (pedited) pedited->wavelet.reschro = true; }
	if (keyFile.has_key ("Wavelet", "ResidualTM"))     {wavelet.tmrs  = keyFile.get_double  ("Wavelet", "ResidualTM"); if (pedited) pedited->wavelet.tmrs = true; }
	if (keyFile.has_key ("Wavelet", "Residualgamma"))     {wavelet.gamma  = keyFile.get_double  ("Wavelet", "Residualgamma"); if (pedited) pedited->wavelet.gamma = true; }
    if (keyFile.has_key ("Wavelet", "ContExtra"))     {wavelet.sup  = keyFile.get_integer  ("Wavelet", "ContExtra"); if (pedited) pedited->wavelet.sup = true; }
    if (keyFile.has_key ("Wavelet", "HueRangeResidual"))     {wavelet.sky  = keyFile.get_double  ("Wavelet", "HueRangeResidual"); if (pedited) pedited->wavelet.sky = true; }
    if (keyFile.has_key ("Wavelet", "MaxLev"))     {wavelet.thres  = keyFile.get_integer  ("Wavelet", "MaxLev"); if (pedited) pedited->wavelet.thres = true; }
    if (keyFile.has_key ("Wavelet", "ThresholdHighLight"))     {wavelet.threshold  = keyFile.get_integer  ("Wavelet", "ThresholdHighLight"); if (pedited) pedited->wavelet.threshold = true; }
    if (keyFile.has_key ("Wavelet", "ThresholdShadow"))     {wavelet.threshold2  = keyFile.get_integer  ("Wavelet", "ThresholdShadow"); if (pedited) pedited->wavelet.threshold2 = true; }
    if (keyFile.has_key ("Wavelet", "Edgedetect"))     {wavelet.edgedetect  = keyFile.get_integer  ("Wavelet", "Edgedetect"); if (pedited) pedited->wavelet.edgedetect = true; }
    if (keyFile.has_key ("Wavelet", "Edgedetectthr"))     {wavelet.edgedetectthr  = keyFile.get_integer  ("Wavelet", "Edgedetectthr"); if (pedited) pedited->wavelet.edgedetectthr = true; }
    if (keyFile.has_key ("Wavelet", "EdgedetectthrHi"))     {wavelet.edgedetectthr2  = keyFile.get_integer  ("Wavelet", "EdgedetectthrHi"); if (pedited) pedited->wavelet.edgedetectthr2 = true; }
    if (keyFile.has_key ("Wavelet", "ThresholdChroma"))     {wavelet.chroma  = keyFile.get_integer  ("Wavelet", "ThresholdChroma"); if (pedited) pedited->wavelet.chroma = true; }
    if (keyFile.has_key ("Wavelet", "ChromaLink"))     {wavelet.chro  = keyFile.get_integer  ("Wavelet", "ChromaLink"); if (pedited) pedited->wavelet.chro = true; }
    if (keyFile.has_key ("Wavelet", "Contrast"))     {wavelet.contrast  = keyFile.get_integer  ("Wavelet", "Contrast"); if (pedited) pedited->wavelet.contrast = true; }
    if (keyFile.has_key ("Wavelet", "Edgrad"))     {wavelet.edgrad  = keyFile.get_integer  ("Wavelet", "Edgrad"); if (pedited) pedited->wavelet.edgrad = true; }
    if (keyFile.has_key ("Wavelet", "Edgval"))     {wavelet.edgval  = keyFile.get_integer  ("Wavelet", "Edgval"); if (pedited) pedited->wavelet.edgval = true; }
    if (keyFile.has_key ("Wavelet", "ThrEdg"))     {wavelet.edgthresh  = keyFile.get_integer  ("Wavelet", "ThrEdg"); if (pedited) pedited->wavelet.edgthresh = true; }
    if (keyFile.has_key ("Wavelet", "ThresholdResidShadow"))     {wavelet.thr  = keyFile.get_integer  ("Wavelet", "ThresholdResidShadow"); if (pedited) pedited->wavelet.thr = true; }
	if (keyFile.has_key ("Wavelet", "ThresholdResidHighLight"))     {wavelet.thrH  = keyFile.get_integer  ("Wavelet", "ThresholdResidHighLight"); if (pedited) pedited->wavelet.thrH = true; }
    if (keyFile.has_key ("Wavelet", "ContrastCurve")) {wavelet.ccwcurve = keyFile.get_double_list ("Wavelet", "ContrastCurve"); if (pedited) pedited->wavelet.ccwcurve = true; }
    if (keyFile.has_key ("Wavelet", "OpacityCurveRG"))    { wavelet.opacityCurveRG = keyFile.get_double_list ("Wavelet", "OpacityCurveRG"); if (pedited) pedited->wavelet.opacityCurveRG = true; }
    if (keyFile.has_key ("Wavelet", "OpacityCurveBY"))    { wavelet.opacityCurveBY = keyFile.get_double_list ("Wavelet", "OpacityCurveBY"); if (pedited) pedited->wavelet.opacityCurveBY = true; }
    if (keyFile.has_key ("Wavelet", "OpacityCurveW"))    { wavelet.opacityCurveW = keyFile.get_double_list ("Wavelet", "OpacityCurveW"); if (pedited) pedited->wavelet.opacityCurveW = true; }
    if (keyFile.has_key ("Wavelet", "OpacityCurveWL"))    { wavelet.opacityCurveWL = keyFile.get_double_list ("Wavelet", "OpacityCurveWL"); if (pedited) pedited->wavelet.opacityCurveWL = true; }
    if (keyFile.has_key ("Wavelet", "HHcurve"))    { wavelet.hhcurve = keyFile.get_double_list ("Wavelet", "HHcurve"); if (pedited) pedited->wavelet.hhcurve = true; }
    if (keyFile.has_key ("Wavelet", "CHcurve"))    { wavelet.Chcurve = keyFile.get_double_list ("Wavelet", "CHcurve"); if (pedited) pedited->wavelet.Chcurve = true; }
    if (keyFile.has_key ("Wavelet", "WavclCurve"))    { wavelet.wavclCurve = keyFile.get_double_list ("Wavelet", "WavclCurve"); if (pedited) pedited->wavelet.wavclCurve = true; }
	if (keyFile.has_key ("Wavelet", "Hueskin"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Wavelet", "Hueskin");
        wavelet.hueskin.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->wavelet.hueskin = true;
    }
	if (keyFile.has_key ("Wavelet", "HueRange"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Wavelet", "HueRange");
        wavelet.hueskin2.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->wavelet.hueskin2 = true;
    }
	
    if (keyFile.has_key ("Wavelet", "HLRange"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Wavelet", "HLRange");
        wavelet.hllev.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->wavelet.hllev = true;
    }
    if (keyFile.has_key ("Wavelet", "SHRange"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Wavelet", "SHRange");
        wavelet.bllev.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->wavelet.bllev = true;
    }
    if (keyFile.has_key ("Wavelet", "Edgcont"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Wavelet", "Edgcont");
        wavelet.edgcont.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->wavelet.edgcont = true;
    }
    if (keyFile.has_key ("Wavelet", "Level0noise"))   {
        Glib::ArrayHandle<double> thresh = keyFile.get_double_list ("Wavelet", "Level0noise");
        wavelet.level0noise.setValues(thresh.data()[0], thresh.data()[1]);
        if (pedited) pedited->wavelet.level0noise = true;
    }
    if (keyFile.has_key ("Wavelet", "Level1noise"))   {
        Glib::ArrayHandle<double> thresh = keyFile.get_double_list ("Wavelet", "Level1noise");
        wavelet.level1noise.setValues(thresh.data()[0], thresh.data()[1]);
        if (pedited) pedited->wavelet.level1noise = true;
    }
    if (keyFile.has_key ("Wavelet", "Level2noise"))   {
        Glib::ArrayHandle<double> thresh = keyFile.get_double_list ("Wavelet", "Level2noise");
        wavelet.level2noise.setValues(thresh.data()[0], thresh.data()[1]);
        if (pedited) pedited->wavelet.level2noise = true;
    }
	
    if (keyFile.has_key ("Wavelet", "Pastlev"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Wavelet", "Pastlev");
        wavelet.pastlev.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->wavelet.pastlev = true;
    }
    if (keyFile.has_key ("Wavelet", "Satlev"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Wavelet", "Satlev");
        wavelet.satlev.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->wavelet.satlev = true;
    }
	
	
    if(keyFile.has_key ("Wavelet", "Skinprotect")) { wavelet.skinprotect = keyFile.get_double ("Wavelet", "Skinprotect"); if (pedited) pedited->wavelet.skinprotect = true; }
    for(int i = 0; i < 9; i ++)
    {
        std::stringstream ss;
        ss << "Contrast" << (i+1);
        if(keyFile.has_key ("Wavelet", ss.str())) {wavelet.c[i] = keyFile.get_integer ("Wavelet", ss.str()); if (pedited) pedited->wavelet.c[i]   = true;} 
    }
    for(int i = 0; i < 9; i ++)
    {
        std::stringstream ss;
        ss << "Chroma" << (i+1);
        if(keyFile.has_key ("Wavelet", ss.str())) {wavelet.ch[i] = keyFile.get_integer ("Wavelet", ss.str()); if (pedited) pedited->wavelet.ch[i]   = true;} 
    }
	
}

    // load directional pyramid equalizer parameters
if (keyFile.has_group ("Directional Pyramid Equalizer")) {
    if (keyFile.has_key ("Directional Pyramid Equalizer", "Enabled"))   { dirpyrequalizer.enabled = keyFile.get_boolean ("Directional Pyramid Equalizer", "Enabled"); if (pedited) pedited->dirpyrequalizer.enabled = true; }
    if (keyFile.has_key ("Directional Pyramid Equalizer", "Gamutlab"))  { dirpyrequalizer.gamutlab = keyFile.get_boolean ("Directional Pyramid Equalizer", "Gamutlab"); if (pedited) pedited->dirpyrequalizer.gamutlab = true; }
 //   if (keyFile.has_key ("Directional Pyramid Equalizer", "Algorithm")) { dirpyrequalizer.algo = keyFile.get_string ("Directional Pyramid Equalizer", "Algorithm"); if (pedited) pedited->dirpyrequalizer.algo = true; }
    if (keyFile.has_key ("Directional Pyramid Equalizer", "Hueskin"))   {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Directional Pyramid Equalizer", "Hueskin");
        dirpyrequalizer.hueskin.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 300), min(thresh.data()[3], 300));
        if (pedited) pedited->dirpyrequalizer.hueskin = true;
    }

    if (ppVersion < 316) {
        for(int i = 0; i < 5; i ++) {
            std::stringstream ss;
            ss << "Mult" << i;
            if(keyFile.has_key ("Directional Pyramid Equalizer", ss.str())) {
                if(i==4) { dirpyrequalizer.threshold = keyFile.get_double ("Directional Pyramid Equalizer", ss.str()); if (pedited) pedited->dirpyrequalizer.threshold = true; }
                else     { dirpyrequalizer.mult[i]   = keyFile.get_double ("Directional Pyramid Equalizer", ss.str()); if (pedited) pedited->dirpyrequalizer.mult[i]   = true; }
            }
        }
        dirpyrequalizer.mult[4] = 1.0;
    }
    else {
        // 5 level wavelet + dedicated threshold parameter
        for(int i = 0; i < 6; i ++) {
            std::stringstream ss;
            ss << "Mult" << i;
            if(keyFile.has_key ("Directional Pyramid Equalizer", ss.str())) { dirpyrequalizer.mult[i]  = keyFile.get_double ("Directional Pyramid Equalizer", ss.str()); if (pedited) pedited->dirpyrequalizer.mult[i]   = true; }
        }
        if(keyFile.has_key ("Directional Pyramid Equalizer", "Threshold"))   { dirpyrequalizer.threshold = keyFile.get_double ("Directional Pyramid Equalizer", "Threshold"); if (pedited) pedited->dirpyrequalizer.threshold = true; }
        if(keyFile.has_key ("Directional Pyramid Equalizer", "Skinprotect")) { dirpyrequalizer.skinprotect = keyFile.get_double ("Directional Pyramid Equalizer", "Skinprotect"); if (pedited) pedited->dirpyrequalizer.skinprotect = true; }
    }
}

    // load CLUT parameters
if ( keyFile.has_group( "Film Simulation" ) )
{
    if ( keyFile.has_key( "Film Simulation", "Enabled" ) )      { filmSimulation.enabled = keyFile.get_boolean( "Film Simulation", "Enabled" ); if ( pedited ) pedited->filmSimulation.enabled = true; }
    if ( keyFile.has_key( "Film Simulation", "ClutFilename" ) ) { filmSimulation.clutFilename = keyFile.get_string( "Film Simulation", "ClutFilename" ); if ( pedited ) pedited->filmSimulation.clutFilename = true; }
    if ( keyFile.has_key( "Film Simulation", "Strength" ) )     {
        if (ppVersion < 321)
            filmSimulation.strength = int(keyFile.get_double( "Film Simulation", "Strength" )*100 + 0.1);
        else
            filmSimulation.strength = keyFile.get_integer( "Film Simulation", "Strength" );
        if ( pedited ) pedited->filmSimulation.strength = true;
    }
}

    // load HSV wavelet parameters
if (keyFile.has_group ("HSV Equalizer")) {
    if (ppVersion>=300) {
        if (keyFile.has_key ("HSV Equalizer", "HCurve")) { hsvequalizer.hcurve = keyFile.get_double_list ("HSV Equalizer", "HCurve"); if (pedited) pedited->hsvequalizer.hcurve = true; }
        if (keyFile.has_key ("HSV Equalizer", "SCurve")) { hsvequalizer.scurve = keyFile.get_double_list ("HSV Equalizer", "SCurve"); if (pedited) pedited->hsvequalizer.scurve = true; }
        if (keyFile.has_key ("HSV Equalizer", "VCurve")) { hsvequalizer.vcurve = keyFile.get_double_list ("HSV Equalizer", "VCurve"); if (pedited) pedited->hsvequalizer.vcurve = true; }
    }
}

    // load RGB curves
if (keyFile.has_group ("RGB Curves")) {
    if (keyFile.has_key ("RGB Curves", "LumaMode"))  { rgbCurves.lumamode = keyFile.get_boolean ("RGB Curves", "LumaMode"); if (pedited) pedited->rgbCurves.lumamode = true; }
    if (keyFile.has_key ("RGB Curves", "rCurve")) { rgbCurves.rcurve = keyFile.get_double_list ("RGB Curves", "rCurve"); if (pedited) pedited->rgbCurves.rcurve = true; }
    if (keyFile.has_key ("RGB Curves", "gCurve")) { rgbCurves.gcurve = keyFile.get_double_list ("RGB Curves", "gCurve"); if (pedited) pedited->rgbCurves.gcurve = true; }
    if (keyFile.has_key ("RGB Curves", "bCurve")) { rgbCurves.bcurve  = keyFile.get_double_list ("RGB Curves", "bCurve"); if (pedited) pedited->rgbCurves.bcurve = true; }
}

    // load Color Toning
if (keyFile.has_group ("ColorToning")) {
    if (keyFile.has_key ("ColorToning", "Enabled"))       { colorToning.enabled = keyFile.get_boolean ("ColorToning", "Enabled"); if (pedited) pedited->colorToning.enabled = true; }
    if (keyFile.has_key ("ColorToning", "Method"))        { colorToning.method = keyFile.get_string ("ColorToning", "Method"); if (pedited) pedited->colorToning.method = true; }
    if (keyFile.has_key ("ColorToning", "Lumamode"))      { colorToning.lumamode = keyFile.get_boolean ("ColorToning", "Lumamode"); if (pedited) pedited->colorToning.lumamode = true; }
    if (keyFile.has_key ("ColorToning", "Twocolor"))      { colorToning.twocolor = keyFile.get_string ("ColorToning", "Twocolor"); if (pedited) pedited->colorToning.twocolor = true; }
    if (keyFile.has_key ("ColorToning", "OpacityCurve"))  { colorToning.opacityCurve = keyFile.get_double_list ("ColorToning", "OpacityCurve"); if (pedited) pedited->colorToning.opacityCurve = true; }
    if (keyFile.has_key ("ColorToning", "ColorCurve"))    { colorToning.colorCurve = keyFile.get_double_list ("ColorToning", "ColorCurve"); if (pedited) pedited->colorToning.colorCurve = true; }
    if (keyFile.has_key ("ColorToning", "Autosat"))       { colorToning.autosat = keyFile.get_boolean ("ColorToning", "Autosat"); if (pedited) pedited->colorToning.autosat = true; }
    if (keyFile.has_key ("ColorToning", "SatProtectionThreshold"))    { colorToning.satProtectionThreshold = keyFile.get_integer ("ColorToning", "SatProtectionThreshold"); if (pedited) pedited->colorToning.satprotectionthreshold = true; }
    if (keyFile.has_key ("ColorToning", "SaturatedOpacity"))          { colorToning.saturatedOpacity = keyFile.get_integer ("ColorToning", "SaturatedOpacity"); if (pedited) pedited->colorToning.saturatedopacity = true; }
    if (keyFile.has_key ("ColorToning", "Strength"))                  { colorToning.strength = keyFile.get_integer ("ColorToning", "Strength"); if (pedited) pedited->colorToning.strength = true; }
    if (keyFile.has_key ("ColorToning", "HighlightsColorSaturation")) {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("ColorToning", "HighlightsColorSaturation");
        colorToning.hlColSat.setValues(thresh.data()[0], thresh.data()[1]);
        if (pedited) pedited->colorToning.hlColSat = true;
    }
    if (keyFile.has_key ("ColorToning", "ShadowsColorSaturation")) {
        Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("ColorToning", "ShadowsColorSaturation");
        colorToning.shadowsColSat.setValues(thresh.data()[0], thresh.data()[1]);
        if (pedited) pedited->colorToning.shadowsColSat = true;
    }
    if (keyFile.has_key ("ColorToning", "ClCurve"))       { colorToning.clcurve = keyFile.get_double_list ("ColorToning", "ClCurve"); if (pedited) pedited->colorToning.clcurve = true; }
    if (keyFile.has_key ("ColorToning", "Cl2Curve"))      { colorToning.cl2curve = keyFile.get_double_list ("ColorToning", "Cl2Curve"); if (pedited) pedited->colorToning.cl2curve = true; }
    if (keyFile.has_key ("ColorToning", "Redlow"))        { colorToning.redlow = keyFile.get_double ("ColorToning", "Redlow"); if (pedited) pedited->colorToning.redlow = true; }
    if (keyFile.has_key ("ColorToning", "Greenlow"))      { colorToning.greenlow = keyFile.get_double ("ColorToning", "Greenlow"); if (pedited) pedited->colorToning.greenlow = true; }
    if (keyFile.has_key ("ColorToning", "Bluelow"))       { colorToning.bluelow = keyFile.get_double ("ColorToning", "Bluelow"); if (pedited) pedited->colorToning.bluelow = true; }
    if (keyFile.has_key ("ColorToning", "Satlow"))        { colorToning.satlow = keyFile.get_double ("ColorToning", "Satlow"); if (pedited) pedited->colorToning.satlow = true; }
    if (keyFile.has_key ("ColorToning", "Balance"))       { colorToning.balance = keyFile.get_integer ("ColorToning", "Balance"); if (pedited) pedited->colorToning.balance = true; }
    if (keyFile.has_key ("ColorToning", "Sathigh"))       { colorToning.sathigh = keyFile.get_double ("ColorToning", "Sathigh"); if (pedited) pedited->colorToning.sathigh = true; }
    if (keyFile.has_key ("ColorToning", "Redmed"))        { colorToning.redmed = keyFile.get_double ("ColorToning", "Redmed"); if (pedited) pedited->colorToning.redmed = true; }
    if (keyFile.has_key ("ColorToning", "Greenmed"))      { colorToning.greenmed = keyFile.get_double ("ColorToning", "Greenmed"); if (pedited) pedited->colorToning.greenmed = true; }
    if (keyFile.has_key ("ColorToning", "Bluemed"))       { colorToning.bluemed = keyFile.get_double ("ColorToning", "Bluemed"); if (pedited) pedited->colorToning.bluemed = true; }
    if (keyFile.has_key ("ColorToning", "Redhigh"))       { colorToning.redhigh = keyFile.get_double ("ColorToning", "Redhigh"); if (pedited) pedited->colorToning.redhigh = true; }
    if (keyFile.has_key ("ColorToning", "Greenhigh"))     { colorToning.greenhigh = keyFile.get_double ("ColorToning", "Greenhigh"); if (pedited) pedited->colorToning.greenhigh = true; }
    if (keyFile.has_key ("ColorToning", "Bluehigh"))      { colorToning.bluehigh = keyFile.get_double ("ColorToning", "Bluehigh"); if (pedited) pedited->colorToning.bluehigh = true; }
}

    // load raw settings
if (keyFile.has_group ("RAW")) {
    if (keyFile.has_key ("RAW", "DarkFrame"))                { raw.dark_frame = expandRelativePath(fname, "", keyFile.get_string  ("RAW", "DarkFrame" )); if (pedited) pedited->raw.darkFrame = true; }
    if (keyFile.has_key ("RAW", "DarkFrameAuto"))            { raw.df_autoselect = keyFile.get_boolean ("RAW", "DarkFrameAuto" ); if (pedited) pedited->raw.dfAuto = true; }
    if (keyFile.has_key ("RAW", "FlatFieldFile"))            { raw.ff_file = expandRelativePath(fname, "", keyFile.get_string  ("RAW", "FlatFieldFile" )); if (pedited) pedited->raw.ff_file = true; }
    if (keyFile.has_key ("RAW", "FlatFieldAutoSelect"))      { raw.ff_AutoSelect = keyFile.get_boolean  ("RAW", "FlatFieldAutoSelect" );  if (pedited) pedited->raw.ff_AutoSelect = true; }
    if (keyFile.has_key ("RAW", "FlatFieldBlurRadius"))      { raw.ff_BlurRadius = keyFile.get_integer  ("RAW", "FlatFieldBlurRadius" ); if (pedited) pedited->raw.ff_BlurRadius = true; }
    if (keyFile.has_key ("RAW", "FlatFieldBlurType"))        { raw.ff_BlurType = keyFile.get_string  ("RAW", "FlatFieldBlurType" ); if (pedited) pedited->raw.ff_BlurType = true; }
    if (keyFile.has_key ("RAW", "FlatFieldAutoClipControl")) { raw.ff_AutoClipControl = keyFile.get_boolean  ("RAW", "FlatFieldAutoClipControl" );  if (pedited) pedited->raw.ff_AutoClipControl = true; }
    if (keyFile.has_key ("RAW", "FlatFieldClipControl"))     { raw.ff_clipControl = keyFile.get_boolean  ("RAW", "FlatFieldClipControl" );  if (pedited) pedited->raw.ff_clipControl = true; }
    if (keyFile.has_key ("RAW", "CA"))                       { raw.ca_autocorrect = keyFile.get_boolean ("RAW", "CA" ); if (pedited) pedited->raw.caCorrection = true; }
    if (keyFile.has_key ("RAW", "CARed"))                    { raw.cared = keyFile.get_double ("RAW", "CARed" ); if (pedited) pedited->raw.caRed = true; }
    if (keyFile.has_key ("RAW", "CABlue"))                   { raw.cablue = keyFile.get_double ("RAW", "CABlue" ); if (pedited) pedited->raw.caBlue = true; }
    // for compatibility to elder pp3 versions
    if (keyFile.has_key ("RAW", "HotDeadPixels"))            { raw.deadPixelFilter = raw.hotPixelFilter = keyFile.get_boolean ("RAW", "HotDeadPixels" ); if (pedited) pedited->raw.hotPixelFilter = pedited->raw.deadPixelFilter = true; }
    if (keyFile.has_key ("RAW", "HotPixelFilter"))           { raw.hotPixelFilter = keyFile.get_boolean ("RAW", "HotPixelFilter" ); if (pedited) pedited->raw.hotPixelFilter = true; }
    if (keyFile.has_key ("RAW", "DeadPixelFilter"))          { raw.deadPixelFilter = keyFile.get_boolean ("RAW", "DeadPixelFilter" ); if (pedited) pedited->raw.deadPixelFilter = true; }
    
    if (keyFile.has_key ("RAW", "HotDeadPixelThresh"))       { raw.hotdeadpix_thresh = keyFile.get_integer ("RAW", "HotDeadPixelThresh" ); if (pedited) pedited->raw.hotDeadPixelThresh = true; }
    if (keyFile.has_key ("RAW", "PreExposure"))              { raw.expos =keyFile.get_double("RAW", "PreExposure"); if (pedited) pedited->raw.exPos = true; }
    if (keyFile.has_key ("RAW", "PrePreserv"))               { raw.preser =keyFile.get_double("RAW", "PrePreserv"); if (pedited) pedited->raw.exPreser = true; }

    if (ppVersion < 320) {
        if (keyFile.has_key ("RAW", "Method"))           { raw.bayersensor.method = keyFile.get_string ("RAW", "Method"); if (pedited) pedited->raw.bayersensor.method = true; }
        if (keyFile.has_key ("RAW", "CcSteps"))          { raw.bayersensor.ccSteps  = keyFile.get_integer ("RAW", "CcSteps"); if (pedited) pedited->raw.bayersensor.ccSteps = true; }
        if (keyFile.has_key ("RAW", "LineDenoise"))      { raw.bayersensor.linenoise = keyFile.get_integer ("RAW", "LineDenoise" ); if (pedited) pedited->raw.bayersensor.linenoise = true; }
        if (keyFile.has_key ("RAW", "GreenEqThreshold")) { raw.bayersensor.greenthresh= keyFile.get_integer ("RAW", "GreenEqThreshold"); if (pedited) pedited->raw.bayersensor.greenEq = true; }
        if (keyFile.has_key ("RAW", "DCBIterations"))    { raw.bayersensor.dcb_iterations = keyFile.get_integer("RAW", "DCBIterations"); if (pedited) pedited->raw.bayersensor.dcbIterations = true; }
        if (keyFile.has_key ("RAW", "DCBEnhance"))       { raw.bayersensor.dcb_enhance = keyFile.get_boolean("RAW", "DCBEnhance"); if (pedited) pedited->raw.bayersensor.dcbEnhance = true; }
        if (keyFile.has_key ("RAW", "LMMSEIterations"))  { raw.bayersensor.lmmse_iterations = keyFile.get_integer("RAW", "LMMSEIterations"); if (pedited) pedited->raw.bayersensor.lmmseIterations = true; }
        if (keyFile.has_key ("RAW", "PreBlackzero"))     { raw.bayersensor.black0 = keyFile.get_double("RAW", "PreBlackzero"); if (pedited) pedited->raw.bayersensor.exBlack0 = true; }
        if (keyFile.has_key ("RAW", "PreBlackone"))      { raw.bayersensor.black1 = keyFile.get_double("RAW", "PreBlackone"); if (pedited) pedited->raw.bayersensor.exBlack1 = true; }
        if (keyFile.has_key ("RAW", "PreBlacktwo"))      { raw.bayersensor.black2 = keyFile.get_double("RAW", "PreBlacktwo"); if (pedited) pedited->raw.bayersensor.exBlack2 = true; }
        if (keyFile.has_key ("RAW", "PreBlackthree"))    { raw.bayersensor.black3 = keyFile.get_double("RAW", "PreBlackthree"); if (pedited) pedited->raw.bayersensor.exBlack3 = true; }
        if (keyFile.has_key ("RAW", "PreTwoGreen"))      { raw.bayersensor.twogreen = keyFile.get_boolean("RAW", "PreTwoGreen"); if (pedited) pedited->raw.bayersensor.exTwoGreen = true; }
        //if (keyFile.has_key ("RAW", "ALLEnhance"))     { raw.bayersensor.all_enhance = keyFile.get_boolean("RAW", "ALLEnhance"); if (pedited) pedited->raw.bayersensor.allEnhance = true; }
    }
}

// load Bayer sensors' raw settings
if (keyFile.has_group ("RAW Bayer")) {
    if (keyFile.has_key ("RAW Bayer", "Method"))           { raw.bayersensor.method = keyFile.get_string ("RAW Bayer", "Method"); if (pedited) pedited->raw.bayersensor.method = true; }
    if (keyFile.has_key ("RAW Bayer", "CcSteps"))          { raw.bayersensor.ccSteps  = keyFile.get_integer ("RAW Bayer", "CcSteps"); if (pedited) pedited->raw.bayersensor.ccSteps = true; }
    if (keyFile.has_key ("RAW Bayer", "PreBlack0"))        { raw.bayersensor.black0 = keyFile.get_double("RAW Bayer", "PreBlack0"); if (pedited) pedited->raw.bayersensor.exBlack0 = true; }
    if (keyFile.has_key ("RAW Bayer", "PreBlack1"))        { raw.bayersensor.black1 = keyFile.get_double("RAW Bayer", "PreBlack1"); if (pedited) pedited->raw.bayersensor.exBlack1 = true; }
    if (keyFile.has_key ("RAW Bayer", "PreBlack2"))        { raw.bayersensor.black2 = keyFile.get_double("RAW Bayer", "PreBlack2"); if (pedited) pedited->raw.bayersensor.exBlack2 = true; }
    if (keyFile.has_key ("RAW Bayer", "PreBlack3"))        { raw.bayersensor.black3 = keyFile.get_double("RAW Bayer", "PreBlack3"); if (pedited) pedited->raw.bayersensor.exBlack3 = true; }
    if (keyFile.has_key ("RAW Bayer", "PreTwoGreen"))      { raw.bayersensor.twogreen = keyFile.get_boolean("RAW Bayer", "PreTwoGreen"); if (pedited) pedited->raw.bayersensor.exTwoGreen = true; }
    if (keyFile.has_key ("RAW Bayer", "LineDenoise"))      { raw.bayersensor.linenoise = keyFile.get_integer ("RAW Bayer", "LineDenoise" ); if (pedited) pedited->raw.bayersensor.linenoise = true; }
    if (keyFile.has_key ("RAW Bayer", "GreenEqThreshold")) { raw.bayersensor.greenthresh= keyFile.get_integer ("RAW Bayer", "GreenEqThreshold"); if (pedited) pedited->raw.bayersensor.greenEq = true; }
    if (keyFile.has_key ("RAW Bayer", "DCBIterations"))    { raw.bayersensor.dcb_iterations = keyFile.get_integer("RAW Bayer", "DCBIterations"); if (pedited) pedited->raw.bayersensor.dcbIterations = true; }
    if (keyFile.has_key ("RAW Bayer", "DCBEnhance"))       { raw.bayersensor.dcb_enhance = keyFile.get_boolean("RAW Bayer", "DCBEnhance"); if (pedited) pedited->raw.bayersensor.dcbEnhance = true; }
    if (keyFile.has_key ("RAW Bayer", "LMMSEIterations"))  { raw.bayersensor.lmmse_iterations = keyFile.get_integer("RAW Bayer", "LMMSEIterations"); if (pedited) pedited->raw.bayersensor.lmmseIterations = true; }
    //if (keyFile.has_key ("RAW Bayer", "ALLEnhance"))     { raw.bayersensor.all_enhance = keyFile.get_boolean("RAW Bayer", "ALLEnhance"); if (pedited) pedited->raw.bayersensor.allEnhance = true; }
}

// load X-Trans sensors' raw settings
if (keyFile.has_group ("RAW X-Trans")) {
    if (keyFile.has_key ("RAW X-Trans", "Method"))           { raw.xtranssensor.method = keyFile.get_string ("RAW X-Trans", "Method"); if (pedited) pedited->raw.xtranssensor.method = true; }
    if (keyFile.has_key ("RAW X-Trans", "CcSteps"))          { raw.xtranssensor.ccSteps  = keyFile.get_integer ("RAW X-Trans", "CcSteps"); if (pedited) pedited->raw.xtranssensor.ccSteps = true; }
    if (keyFile.has_key ("RAW X-Trans", "PreBlackRed"))      { raw.xtranssensor.blackred = keyFile.get_double("RAW X-Trans", "PreBlackRed"); if (pedited) pedited->raw.xtranssensor.exBlackRed = true; }
    if (keyFile.has_key ("RAW X-Trans", "PreBlackGreen"))    { raw.xtranssensor.blackgreen = keyFile.get_double("RAW X-Trans", "PreBlackGreen"); if (pedited) pedited->raw.xtranssensor.exBlackGreen = true; }
    if (keyFile.has_key ("RAW X-Trans", "PreBlackBlue"))     { raw.xtranssensor.blackblue = keyFile.get_double("RAW X-Trans", "PreBlackBlue"); if (pedited) pedited->raw.xtranssensor.exBlackBlue = true; }
}

    // load exif change settings
if (keyFile.has_group ("Exif")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("Exif");
    for (int i=0; i<(int)keys.size(); i++) {
        Glib::ustring tmpStr = keyFile.get_string ("Exif", keys[i]);
        exif[keys[i]] = keyFile.get_string ("Exif", keys[i]);
        if (pedited) pedited->exif = true;
    }
}

    /*
     * Load iptc change settings
     *
     * Existing values are preserved, and the stored values
     * are added to the list. To reset a field, the user has to
     * save the profile with the field leaved empty, but still
     * terminated by a semi-column ";"
     *
     * Please note that the old Keywords and SupplementalCategories
     * tag content is fully replaced by the new one,
     * i.e. they don't merge
     */
if (keyFile.has_group ("IPTC")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("IPTC");
    IPTCPairs::iterator element;
    for (unsigned int i=0; i<keys.size(); i++) {
        // does this key already exist?
        element = iptc.find(keys[i]);
        if (element != iptc.end()) {
            // it already exist so we cleanup the values
            element->second.clear();
        }

        // TODO: look out if merging Keywords and SupplementalCategories from the procparams chain would be interesting
        std::vector<Glib::ustring> currIptc = keyFile.get_string_list ("IPTC", keys[i]);
        for (
            std::vector<Glib::ustring>::iterator currLoadedTagValue=currIptc.begin();
            currLoadedTagValue!=currIptc.end();
            currLoadedTagValue++)
        {
            iptc[keys[i]].push_back(currLoadedTagValue->data());
        }
        if (pedited) pedited->iptc = true;
    }
}


        return 0;
    }
    catch (const Glib::Error& e) {
        printf ("-->%s\n", e.what().c_str());
        return 1;
    }
    catch (...) {
        printf ("-->unknown exception!\n");
        return 1;
    }
    return 0;
}

const Glib::ustring ColorManagementParams::NoICMString = Glib::ustring("No ICM: sRGB output");

bool operator==(const WaveletParams & a, const WaveletParams & b) {
    if(a.enabled != b.enabled)
        return false;

    for(int i = 0; i < 9; i++) {
        if(a.c[i] != b.c[i])
            return false;
    }
    for(int i = 0; i < 9; i++) {
        if(a.ch[i] != b.ch[i])
            return false;
    }
	
    return true;
}



bool operator==(const DirPyrEqualizerParams & a, const DirPyrEqualizerParams & b) {
	if(a.enabled != b.enabled)
		return false;
		
	for(int i = 0; i < 6; i++) {
		if(a.mult[i] != b.mult[i])
			return false;
	}
	if (a.threshold != b.threshold)
		return false;

	return true;
}

/*bool operator==(const ExifPairs& a, const ExifPairs& b) {

    return a.field == b.field && a.value == b.value;
}

bool operator==(const IPTCPairs& a, const IPTCPairs& b) {

    return a.field == b.field && a.values == b.values;
}*/
bool ProcParams::operator== (const ProcParams& other) {

	return
		toneCurve.curve == other.toneCurve.curve
		&& toneCurve.curve2 == other.toneCurve.curve2
		&& toneCurve.brightness == other.toneCurve.brightness
		&& toneCurve.black == other.toneCurve.black
		&& toneCurve.contrast == other.toneCurve.contrast
		&& toneCurve.saturation == other.toneCurve.saturation
		&& toneCurve.shcompr == other.toneCurve.shcompr
		&& toneCurve.hlcompr == other.toneCurve.hlcompr
		&& toneCurve.hlcomprthresh == other.toneCurve.hlcomprthresh
		&& toneCurve.autoexp == other.toneCurve.autoexp
		&& toneCurve.clip == other.toneCurve.clip
		&& toneCurve.expcomp == other.toneCurve.expcomp
		&& toneCurve.curveMode == other.toneCurve.curveMode
		&& toneCurve.curveMode2 == other.toneCurve.curveMode2
		&& toneCurve.hrenabled == other.toneCurve.hrenabled
		&& toneCurve.method == other.toneCurve.method
		&& labCurve.lcurve == other.labCurve.lcurve
		&& labCurve.acurve == other.labCurve.acurve
		&& labCurve.bcurve == other.labCurve.bcurve
		&& labCurve.cccurve == other.labCurve.cccurve
		&& labCurve.chcurve == other.labCurve.chcurve
		&& labCurve.lhcurve == other.labCurve.lhcurve
		&& labCurve.hhcurve == other.labCurve.hhcurve
		&& labCurve.lccurve == other.labCurve.lccurve
		&& labCurve.clcurve == other.labCurve.clcurve
		&& labCurve.brightness == other.labCurve.brightness
		&& labCurve.contrast == other.labCurve.contrast
		&& labCurve.chromaticity == other.labCurve.chromaticity
		&& labCurve.avoidcolorshift == other.labCurve.avoidcolorshift
		&& labCurve.rstprotection == other.labCurve.rstprotection
		&& labCurve.lcredsk == other.labCurve.lcredsk
		&& sharpenEdge.enabled == other.sharpenEdge.enabled
		&& sharpenEdge.passes == other.sharpenEdge.passes
		&& sharpenEdge.amount == other.sharpenEdge.amount
		&& sharpenEdge.threechannels == other.sharpenEdge.threechannels
		&& sharpenMicro.enabled == other.sharpenMicro.enabled
		&& sharpenMicro.matrix == other.sharpenMicro.matrix
		&& sharpenMicro.amount == other.sharpenMicro.amount
		&& sharpenMicro.uniformity == other.sharpenMicro.uniformity
		&& sharpening.enabled == other.sharpening.enabled
		&& sharpening.radius == other.sharpening.radius
		&& sharpening.amount == other.sharpening.amount
		&& sharpening.threshold == other.sharpening.threshold
		&& sharpening.edgesonly == other.sharpening.edgesonly
		&& sharpening.edges_radius == other.sharpening.edges_radius
		&& sharpening.edges_tolerance == other.sharpening.edges_tolerance
		&& sharpening.halocontrol == other.sharpening.halocontrol
		&& sharpening.halocontrol_amount== other.sharpening.halocontrol_amount
		&& sharpening.method == other.sharpening.method
		&& sharpening.deconvamount == other.sharpening.deconvamount
		&& sharpening.deconvradius == other.sharpening.deconvradius
		&& sharpening.deconviter == other.sharpening.deconviter
		&& sharpening.deconvdamping == other.sharpening.deconvdamping
		&& prsharpening.enabled == other.prsharpening.enabled
		&& prsharpening.radius == other.prsharpening.radius
		&& prsharpening.amount == other.prsharpening.amount
		&& prsharpening.threshold == other.prsharpening.threshold
		&& prsharpening.edgesonly == other.prsharpening.edgesonly
		&& prsharpening.edges_radius == other.prsharpening.edges_radius
		&& prsharpening.edges_tolerance == other.prsharpening.edges_tolerance
		&& prsharpening.halocontrol == other.prsharpening.halocontrol
		&& prsharpening.halocontrol_amount== other.prsharpening.halocontrol_amount
		&& prsharpening.method == other.prsharpening.method
		&& prsharpening.deconvamount == other.prsharpening.deconvamount
		&& prsharpening.deconvradius == other.prsharpening.deconvradius
		&& prsharpening.deconviter == other.prsharpening.deconviter
		&& prsharpening.deconvdamping == other.prsharpening.deconvdamping
		&& vibrance.enabled == other.vibrance.enabled
		&& vibrance.pastels == other.vibrance.pastels
		&& vibrance.saturated == other.vibrance.saturated
		&& vibrance.psthreshold == other.vibrance.psthreshold
		&& vibrance.protectskins == other.vibrance.protectskins
		&& vibrance.avoidcolorshift == other.vibrance.avoidcolorshift
		&& vibrance.pastsattog == other.vibrance.pastsattog
		&& vibrance.skintonescurve == other.vibrance.skintonescurve
		//&& colorBoost.amount == other.colorBoost.amount
		//&& colorBoost.avoidclip == other.colorBoost.avoidclip
		//&& colorBoost.enable_saturationlimiter == other.colorBoost.enable_saturationlimiter
		//&& colorBoost.saturationlimit == other.colorBoost.saturationlimit
		&& wb.method == other.wb.method
		&& wb.green == other.wb.green
		&& wb.temperature == other.wb.temperature
		&& wb.equal == other.wb.equal
		//&& colorShift.a == other.colorShift.a
		//&& colorShift.b == other.colorShift.b
		&& colorappearance.enabled == other.colorappearance.enabled
		&& colorappearance.degree == other.colorappearance.degree
		&& colorappearance.autodegree == other.colorappearance.autodegree
		&& colorappearance.surround == other.colorappearance.surround
		&& colorappearance.adapscen == other.colorappearance.adapscen
		&& colorappearance.autoadapscen == other.colorappearance.autoadapscen
		&& colorappearance.adaplum == other.colorappearance.adaplum
		&& colorappearance.badpixsl == other.colorappearance.badpixsl
		&& colorappearance.wbmodel == other.colorappearance.wbmodel
		&& colorappearance.algo == other.colorappearance.algo
		&& colorappearance.curveMode == other.colorappearance.curveMode
		&& colorappearance.curveMode2 == other.colorappearance.curveMode2
		&& colorappearance.curveMode3 == other.colorappearance.curveMode3
		&& colorappearance.jlight == other.colorappearance.jlight
		&& colorappearance.qbright == other.colorappearance.qbright
		&& colorappearance.chroma == other.colorappearance.chroma
		&& colorappearance.schroma == other.colorappearance.schroma
		&& colorappearance.mchroma == other.colorappearance.mchroma
		&& colorappearance.rstprotection == other.colorappearance.rstprotection
		&& colorappearance.contrast == other.colorappearance.contrast
		&& colorappearance.qcontrast == other.colorappearance.qcontrast
		&& colorappearance.colorh == other.colorappearance.colorh
		&& impulseDenoise.enabled == other.impulseDenoise.enabled
		&& impulseDenoise.thresh == other.impulseDenoise.thresh
		&& dirpyrDenoise.enabled == other.dirpyrDenoise.enabled
		&& dirpyrDenoise.enhance == other.dirpyrDenoise.enhance
		&& dirpyrDenoise.median == other.dirpyrDenoise.median
		&& dirpyrDenoise.autochroma == other.dirpyrDenoise.autochroma
//		&& dirpyrDenoise.perform == other.dirpyrDenoise.perform
		&& dirpyrDenoise.luma == other.dirpyrDenoise.luma
		&& dirpyrDenoise.lcurve == other.dirpyrDenoise.lcurve
		&& dirpyrDenoise.cccurve == other.dirpyrDenoise.cccurve
		&& dirpyrDenoise.Ldetail == other.dirpyrDenoise.Ldetail
		&& dirpyrDenoise.chroma == other.dirpyrDenoise.chroma
		&& dirpyrDenoise.dmethod == other.dirpyrDenoise.dmethod
		&& dirpyrDenoise.Lmethod == other.dirpyrDenoise.Lmethod
		&& dirpyrDenoise.Cmethod == other.dirpyrDenoise.Cmethod
		&& dirpyrDenoise.C2method == other.dirpyrDenoise.C2method
		&& dirpyrDenoise.smethod == other.dirpyrDenoise.smethod
		&& dirpyrDenoise.medmethod == other.dirpyrDenoise.medmethod
		&& dirpyrDenoise.methodmed == other.dirpyrDenoise.methodmed
		&& dirpyrDenoise.rgbmethod == other.dirpyrDenoise.rgbmethod
		&& dirpyrDenoise.redchro == other.dirpyrDenoise.redchro
		&& dirpyrDenoise.bluechro == other.dirpyrDenoise.bluechro
		&& dirpyrDenoise.gamma == other.dirpyrDenoise.gamma
		&& dirpyrDenoise.passes == other.dirpyrDenoise.passes
		&& epd.enabled == other.epd.enabled
		&& epd.strength == other.epd.strength
		&& epd.gamma == other.epd.gamma
		&& epd.edgeStopping == other.epd.edgeStopping
		&& epd.scale == other.epd.scale
		&& epd.reweightingIterates == other.epd.reweightingIterates
		&& defringe.enabled == other.defringe.enabled
		&& defringe.radius == other.defringe.radius
		&& defringe.threshold == other.defringe.threshold
		&& defringe.huecurve == other.defringe.huecurve
		
		//&& lumaDenoise.enabled == other.lumaDenoise.enabled
		//&& lumaDenoise.radius == other.lumaDenoise.radius
		//&& lumaDenoise.edgetolerance == other.lumaDenoise.edgetolerance
		//&& colorDenoise.enabled == other.colorDenoise.enabled
		//&& colorDenoise.edgetolerance == other.colorDenoise.edgetolerance
		//&& colorDenoise.edgesensitive == other.colorDenoise.edgesensitive
		&& sh.enabled == other.sh.enabled
		&& sh.hq == other.sh.hq
		&& sh.highlights == other.sh.highlights
		&& sh.htonalwidth == other.sh.htonalwidth
		&& sh.shadows == other.sh.shadows
		&& sh.stonalwidth == other.sh.stonalwidth
		&& sh.localcontrast == other.sh.localcontrast
		&& sh.radius == other.sh.radius
		&& crop.enabled == other.crop.enabled
		&& crop.x == other.crop.x
		&& crop.y == other.crop.y
		&& crop.w == other.crop.w
		&& crop.h == other.crop.h
		&& crop.fixratio == other.crop.fixratio
		&& crop.ratio == other.crop.ratio
		&& crop.orientation == other.crop.orientation
		&& crop.guide == other.crop.guide
		&& coarse.rotate == other.coarse.rotate
		&& coarse.hflip == other.coarse.hflip
		&& coarse.vflip == other.coarse.vflip
		&& rotate.degree == other.rotate.degree
		&& commonTrans.autofill == other.commonTrans.autofill
		&& distortion.amount == other.distortion.amount
		&& lensProf.lcpFile == other.lensProf.lcpFile
		&& lensProf.useDist == other.lensProf.useDist
		&& lensProf.useVign == other.lensProf.useVign
		&& lensProf.useCA == other.lensProf.useCA
		&& perspective.horizontal == other.perspective.horizontal
		&& perspective.vertical == other.perspective.vertical
		&& gradient.enabled == other.gradient.enabled
		&& gradient.degree == other.gradient.degree
		&& gradient.feather == other.gradient.feather
		&& gradient.strength == other.gradient.strength
		&& gradient.centerX == other.gradient.centerX
		&& gradient.centerY == other.gradient.centerY
		&& pcvignette.enabled == other.pcvignette.enabled
		&& pcvignette.strength == other.pcvignette.strength
		&& pcvignette.feather == other.pcvignette.feather
		&& pcvignette.roundness == other.pcvignette.roundness
		&& cacorrection.red == other.cacorrection.red
		&& cacorrection.blue == other.cacorrection.blue
		&& vignetting.amount == other.vignetting.amount
		&& vignetting.radius == other.vignetting.radius
		&& vignetting.strength == other.vignetting.strength
		&& vignetting.centerX == other.vignetting.centerX
		&& vignetting.centerY == other.vignetting.centerY
		&& !memcmp (&chmixer.red, &other.chmixer.red, 3*sizeof(int))
		&& !memcmp (&chmixer.green, &other.chmixer.green, 3*sizeof(int))
		&& !memcmp (&chmixer.blue, &other.chmixer.blue, 3*sizeof(int))
		&& blackwhite.mixerRed == other.blackwhite.mixerRed
		&& blackwhite.mixerOrange == other.blackwhite.mixerOrange
		&& blackwhite.mixerYellow == other.blackwhite.mixerYellow
		&& blackwhite.mixerGreen == other.blackwhite.mixerGreen
		&& blackwhite.mixerCyan == other.blackwhite.mixerCyan
		&& blackwhite.mixerBlue == other.blackwhite.mixerBlue
		&& blackwhite.mixerMagenta == other.blackwhite.mixerMagenta
		&& blackwhite.mixerPurple == other.blackwhite.mixerPurple
		&& blackwhite.gammaRed == other.blackwhite.gammaRed
		&& blackwhite.gammaGreen == other.blackwhite.gammaGreen
		&& blackwhite.gammaBlue == other.blackwhite.gammaBlue
		&& blackwhite.filter == other.blackwhite.filter
		&& blackwhite.setting == other.blackwhite.setting
		&& blackwhite.method == other.blackwhite.method
		&& blackwhite.luminanceCurve == other.blackwhite.luminanceCurve
		&& blackwhite.beforeCurve == other.blackwhite.beforeCurve
		&& blackwhite.afterCurve == other.blackwhite.afterCurve
		&& blackwhite.beforeCurveMode == other.blackwhite.beforeCurveMode
		&& blackwhite.afterCurveMode == other.blackwhite.afterCurveMode
		&& blackwhite.autoc == other.blackwhite.autoc
		&& blackwhite.algo == other.blackwhite.algo
		&& resize.scale == other.resize.scale
		&& resize.appliesTo == other.resize.appliesTo
		&& resize.method == other.resize.method
		&& resize.dataspec == other.resize.dataspec
		&& resize.width == other.resize.width
		&& resize.height == other.resize.height
		&& raw.bayersensor.method == other.raw.bayersensor.method
		&& raw.bayersensor.ccSteps == other.raw.bayersensor.ccSteps
		&& raw.bayersensor.black0==other.raw.bayersensor.black0
		&& raw.bayersensor.black1==other.raw.bayersensor.black1
		&& raw.bayersensor.black2==other.raw.bayersensor.black2
		&& raw.bayersensor.black3==other.raw.bayersensor.black3
		&& raw.bayersensor.twogreen==other.raw.bayersensor.twogreen
		&& raw.bayersensor.greenthresh == other.raw.bayersensor.greenthresh
		&& raw.bayersensor.linenoise == other.raw.bayersensor.linenoise
		&& raw.bayersensor.dcb_enhance == other.raw.bayersensor.dcb_enhance
		&& raw.bayersensor.dcb_iterations == other.raw.bayersensor.dcb_iterations
		&& raw.xtranssensor.method == other.raw.xtranssensor.method
		&& raw.xtranssensor.ccSteps == other.raw.xtranssensor.ccSteps
		&& raw.xtranssensor.blackred==other.raw.xtranssensor.blackred
		&& raw.xtranssensor.blackgreen==other.raw.xtranssensor.blackgreen
		&& raw.xtranssensor.blackblue==other.raw.xtranssensor.blackblue
		&& raw.dark_frame == other.raw.dark_frame
		&& raw.df_autoselect == other.raw.df_autoselect
		&& raw.ff_file == other.raw.ff_file
		&& raw.ff_AutoSelect == other.raw.ff_AutoSelect
		&& raw.ff_BlurRadius == other.raw.ff_BlurRadius
		&& raw.ff_BlurType == other.raw.ff_BlurType
		&& raw.ff_AutoClipControl == other.raw.ff_AutoClipControl
		&& raw.ff_clipControl == other.raw.ff_clipControl
		&& raw.expos==other.raw.expos
		&& raw.preser==other.raw.preser
		&& raw.ca_autocorrect == other.raw.ca_autocorrect
		&& raw.cared == other.raw.cared
		&& raw.cablue == other.raw.cablue
		&& raw.hotPixelFilter == other.raw.hotPixelFilter
		&& raw.deadPixelFilter == other.raw.deadPixelFilter
		&& raw.hotdeadpix_thresh == other.raw.hotdeadpix_thresh
		&& icm.input == other.icm.input
		&& icm.toneCurve == other.icm.toneCurve
		&& icm.applyLookTable == other.icm.applyLookTable
		&& icm.applyBaselineExposureOffset == other.icm.applyBaselineExposureOffset
		&& icm.applyHueSatMap == other.icm.applyHueSatMap
		&& icm.blendCMSMatrix == other.icm.blendCMSMatrix
		&& icm.dcpIlluminant == other.icm.dcpIlluminant
		&& icm.working == other.icm.working
		&& icm.output == other.icm.output
		&& icm.gamma == other.icm.gamma
		&& icm.freegamma == other.icm.freegamma
		&& icm.gampos == other.icm.gampos
		&& icm.slpos == other.icm.slpos
		&& wavelet == other.wavelet
		&& wavelet.Lmethod == other.wavelet.Lmethod
		&& wavelet.CLmethod == other.wavelet.CLmethod
		&& wavelet.Backmethod == other.wavelet.Backmethod
		&& wavelet.Tilesmethod == other.wavelet.Tilesmethod
		&& wavelet.daubcoeffmethod == other.wavelet.daubcoeffmethod
		&& wavelet.CHmethod == other.wavelet.CHmethod
		&& wavelet.CHSLmethod == other.wavelet.CHSLmethod
		&& wavelet.EDmethod == other.wavelet.EDmethod
		&& wavelet.BAmethod == other.wavelet.BAmethod
		&& wavelet.TMmethod == other.wavelet.TMmethod
		&& wavelet.HSmethod == other.wavelet.HSmethod
		&& wavelet.Dirmethod == other.wavelet.Dirmethod
		&& wavelet.rescon == other.wavelet.rescon
		&& wavelet.resconH == other.wavelet.resconH
		&& wavelet.reschro == other.wavelet.reschro
		&& wavelet.tmrs == other.wavelet.tmrs
		&& wavelet.gamma == other.wavelet.gamma
		&& wavelet.sup == other.wavelet.sup
		&& wavelet.sky == other.wavelet.sky
		&& wavelet.thres == other.wavelet.thres
		&& wavelet.threshold == other.wavelet.threshold
		&& wavelet.chroma == other.wavelet.chroma
		&& wavelet.chro == other.wavelet.chro
		&& wavelet.tmr == other.wavelet.tmr
		&& wavelet.contrast == other.wavelet.contrast
		&& wavelet.median == other.wavelet.median
		&& wavelet.medianlev == other.wavelet.medianlev
		&& wavelet.linkedg == other.wavelet.linkedg
		&& wavelet.cbenab == other.wavelet.cbenab
		&& wavelet.lipst == other.wavelet.lipst
		&& wavelet.Medgreinf == other.wavelet.Medgreinf
		&& wavelet.edgrad == other.wavelet.edgrad
		&& wavelet.edgval == other.wavelet.edgval
		&& wavelet.edgthresh == other.wavelet.edgthresh
		&& wavelet.thr == other.wavelet.thr
		&& wavelet.thrH == other.wavelet.thrH
		&& wavelet.threshold == other.wavelet.threshold
		&& wavelet.threshold2 == other.wavelet.threshold2
		&& wavelet.edgedetect == other.wavelet.edgedetect
		&& wavelet.edgedetectthr == other.wavelet.edgedetectthr
		&& wavelet.edgedetectthr2 == other.wavelet.edgedetectthr2
		&& wavelet.hueskin == other.wavelet.hueskin
		&& wavelet.hueskin2 == other.wavelet.hueskin2
		&& wavelet.hllev == other.wavelet.hllev
		&& wavelet.bllev == other.wavelet.bllev
		&& wavelet.edgcont == other.wavelet.edgcont
		&& wavelet.level0noise == other.wavelet.level0noise
		&& wavelet.level1noise == other.wavelet.level1noise
		&& wavelet.level2noise == other.wavelet.level2noise
		&& wavelet.pastlev == other.wavelet.pastlev
		&& wavelet.satlev == other.wavelet.satlev
		&& wavelet.opacityCurveRG == other.wavelet.opacityCurveRG
		&& wavelet.opacityCurveBY == other.wavelet.opacityCurveBY
		&& wavelet.opacityCurveW == other.wavelet.opacityCurveW
		&& wavelet.opacityCurveWL == other.wavelet.opacityCurveWL
		&& wavelet.hhcurve == other.wavelet.hhcurve
		&& wavelet.Chcurve == other.wavelet.Chcurve
		&& wavelet.ccwcurve == other.wavelet.ccwcurve
		&& wavelet.wavclCurve == other.wavelet.wavclCurve
		&& wavelet.skinprotect == other.wavelet.skinprotect
		&& wavelet.strength == other.wavelet.strength
		&& wavelet.balance == other.wavelet.balance
		&& wavelet.greenhigh == other.wavelet.greenhigh
		&& wavelet.greenmed == other.wavelet.greenmed
		&& wavelet.greenlow == other.wavelet.greenlow
		&& wavelet.bluehigh == other.wavelet.bluehigh
		&& wavelet.bluemed == other.wavelet.bluemed
		&& wavelet.bluelow == other.wavelet.bluelow
		&& wavelet.iter == other.wavelet.iter
		&& dirpyrequalizer == other.dirpyrequalizer
	//	&& dirpyrequalizer.algo == other.dirpyrequalizer.algo
		&& dirpyrequalizer.hueskin == other.dirpyrequalizer.hueskin
		&& dirpyrequalizer.threshold == other.dirpyrequalizer.threshold
		&& dirpyrequalizer.skinprotect == other.dirpyrequalizer.skinprotect
		&& hsvequalizer.hcurve == other.hsvequalizer.hcurve
		&& hsvequalizer.scurve == other.hsvequalizer.scurve
		&& hsvequalizer.vcurve == other.hsvequalizer.vcurve
		&& filmSimulation.enabled == other.filmSimulation.enabled
		&& filmSimulation.clutFilename == other.filmSimulation.clutFilename
		&& filmSimulation.strength == other.filmSimulation.strength
		&& rgbCurves.rcurve == other.rgbCurves.rcurve
		&& rgbCurves.gcurve == other.rgbCurves.gcurve
		&& rgbCurves.bcurve == other.rgbCurves.bcurve
		&& colorToning.enabled == other.colorToning.enabled
		&& colorToning.twocolor == other.colorToning.twocolor
		&& colorToning.method == other.colorToning.method
		&& colorToning.colorCurve == other.colorToning.colorCurve
		&& colorToning.opacityCurve == other.colorToning.opacityCurve
		&& colorToning.autosat == other.colorToning.autosat
		&& colorToning.satProtectionThreshold == other.colorToning.satProtectionThreshold
		&& colorToning.saturatedOpacity == other.colorToning.saturatedOpacity
		&& colorToning.strength == other.colorToning.strength
		&& colorToning.hlColSat == other.colorToning.hlColSat
		&& colorToning.shadowsColSat == other.colorToning.shadowsColSat
		&& colorToning.balance == other.colorToning.balance
		&& colorToning.clcurve == other.colorToning.clcurve
		&& colorToning.cl2curve == other.colorToning.cl2curve
		&& colorToning.redlow == other.colorToning.redlow
		&& colorToning.greenlow == other.colorToning.greenlow
		&& colorToning.bluelow == other.colorToning.bluelow
		&& colorToning.satlow == other.colorToning.satlow
		&& colorToning.sathigh == other.colorToning.sathigh
		&& colorToning.redmed == other.colorToning.redmed
		&& colorToning.greenmed == other.colorToning.greenmed
		&& colorToning.bluemed == other.colorToning.bluemed
		&& colorToning.redhigh == other.colorToning.redhigh
		&& colorToning.greenhigh == other.colorToning.greenhigh
		&& colorToning.bluehigh == other.colorToning.bluehigh
		&& exif==other.exif
		&& iptc==other.iptc;
}

bool ProcParams::operator!= (const ProcParams& other) {

    return !(*this==other);
}

PartialProfile::PartialProfile(bool createInstance, bool paramsEditedValue) {
    if (createInstance) {
        pparams = new ProcParams();
        pedited = new ParamsEdited(paramsEditedValue);
    }
    else {
        pparams = NULL;
        pedited = NULL;
    }
}

PartialProfile::PartialProfile(ProcParams* pp, ParamsEdited* pe, bool fullCopy) {
    if (fullCopy && pp) {
        pparams = new ProcParams(*pp);
    }
    else
        pparams = pp;

    if (fullCopy && pe) {
        pedited = new ParamsEdited(*pe);
    }
    else
        pedited = pe;
}

PartialProfile::PartialProfile(const ProcParams* pp, const ParamsEdited* pe) {
    if (pp) {
        pparams = new ProcParams(*pp);
    }
    else
        pparams = NULL;

    if (pe) {
        pedited = new ParamsEdited(*pe);
    }
    else
        pedited = NULL;
}

int PartialProfile::load (Glib::ustring fName) {
    if (!pparams) pparams = new ProcParams();
    if (!pedited) pedited = new ParamsEdited();
    if (fName == DEFPROFILE_INTERNAL)
        return 0;
    else
        return pparams->load(fName, pedited);
}

void PartialProfile::deleteInstance () {
    if (pparams) { delete pparams; pparams = NULL; }
    if (pedited) { delete pedited; pedited = NULL; }
}

/*
 * Set the all values of the General section to false
 * in order to preserve them in applyTo
 */
void PartialProfile::clearGeneral () {
    if (pedited) {
        pedited->general.colorlabel = false;
        pedited->general.intrash = false;
        pedited->general.rank = false;
    }
}

const void PartialProfile::applyTo(ProcParams *destParams) const {
    if (destParams && pparams && pedited) {
        pedited->combine(*destParams, *pparams, true);
    }
}

void PartialProfile::set(bool v) {
    if (pedited) pedited->set(v);
}

}
}

