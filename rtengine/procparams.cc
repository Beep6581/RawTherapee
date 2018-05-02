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

#include <map>

#include <locale.h>

#include <glib/gstdio.h>

#include "curves.h"
#include "procparams.h"

#include "../rtgui/multilangmgr.h"
#include "../rtgui/options.h"
#include "../rtgui/paramsedited.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/version.h"

using namespace std;
extern Options options;

namespace
{

Glib::ustring expandRelativePath(const Glib::ustring &procparams_fname, const Glib::ustring &prefix, Glib::ustring embedded_fname)
{
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

Glib::ustring relativePathIfInside(const Glib::ustring &procparams_fname, bool fnameAbsolute, Glib::ustring embedded_fname)
{
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

void avoidEmptyCurve(std::vector<double> &curve)
{
    if (curve.empty()) {
        curve.push_back(FCT_Linear);
    }
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    int& value
)
{
    value = keyfile.get_integer(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    double& value
)
{
    value = keyfile.get_double(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool& value
)
{
    value = keyfile.get_boolean(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    Glib::ustring& value
)
{
    value = keyfile.get_string(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    std::vector<double>& value
)
{
    value = keyfile.get_double_list(group_name, key);
    avoidEmptyCurve(value);
}

template<typename T>
bool assignFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool has_params_edited,
    T& value,
    bool& params_edited_value
)
{
    if (keyfile.has_key(group_name, key)) {
        getFromKeyfile(keyfile, group_name, key, value);

        if (has_params_edited) {
            params_edited_value = true;
        }

        return true;
    }

    return false;
}

template<typename T, typename = typename std::enable_if<std::is_enum<T>::value>::type>
bool assignFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool has_params_edited,
    const std::map<std::string, T>& mapping,
    T& value,
    bool& params_edited_value
)
{
    if (keyfile.has_key(group_name, key)) {
        Glib::ustring v;
        getFromKeyfile(keyfile, group_name, key, v);

        const typename std::map<std::string, T>::const_iterator m = mapping.find(v);

        if (m != mapping.end()) {
            value = m->second;
        } else {
            return false;
        }

        if (has_params_edited) {
            params_edited_value = true;
        }

        return true;
    }

    return false;
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    int value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_integer(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    double value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_double(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_boolean(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const Glib::ustring& value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_string(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::vector<int>& value,
    Glib::KeyFile& keyfile
)
{
    const Glib::ArrayHandle<int> list = value;
    keyfile.set_integer_list(group_name, key, list);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::vector<double>& value,
    Glib::KeyFile& keyfile
)
{
    const Glib::ArrayHandle<double> list = value;
    keyfile.set_double_list(group_name, key, list);
}

template<typename T>
bool saveToKeyfile(
    bool save,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const T& value,
    Glib::KeyFile& keyfile
)
{
    if (save) {
        putToKeyfile(group_name, key, value, keyfile);
        return true;
    }

    return false;
}

template<typename T, typename = typename std::enable_if<std::is_enum<T>::value>::type>
bool saveToKeyfile(
    bool save,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::map<T, const char*>& mapping,
    const T& value,
    Glib::KeyFile& keyfile
)
{
    if (save) {
        const typename std::map<T, const char*>::const_iterator m = mapping.find(value);

        if (m != mapping.end()) {
            keyfile.set_string(group_name, key, m->second);
            return true;
        }
    }

    return false;
}

}

namespace rtengine
{

namespace procparams
{

ToneCurveParams::ToneCurveParams() :
    autoexp(false),
    clip(0.02),
    hrenabled(false),
    method("Blend"),
    expcomp(0),
    curve{
    DCT_Linear
},
curve2{
    DCT_Linear
},
curveMode(ToneCurveParams::TcMode::STD),
          curveMode2(ToneCurveParams::TcMode::STD),
          brightness(0),
          black(0),
          contrast(0),
          saturation(0),
          shcompr(50),
          hlcompr(0),
          hlcomprthresh(33),
    histmatching(false),
    clampOOG(true)
{
}

bool ToneCurveParams::operator ==(const ToneCurveParams& other) const
{
    return
        autoexp == other.autoexp
        && clip == other.clip
        && hrenabled == other.hrenabled
        && method == other.method
        && expcomp == other.expcomp
        && curve == other.curve
        && curve2 == other.curve2
        && curveMode == other.curveMode
        && curveMode2 == other.curveMode2
        && brightness == other.brightness
        && black == other.black
        && contrast == other.contrast
        && saturation == other.saturation
        && shcompr == other.shcompr
        && hlcompr == other.hlcompr
        && hlcomprthresh == other.hlcomprthresh
        && histmatching == other.histmatching
        && clampOOG == other.clampOOG;
}

bool ToneCurveParams::operator !=(const ToneCurveParams& other) const
{
    return !(*this == other);
}

RetinexParams::RetinexParams() :
    enabled(false),
    cdcurve{
    DCT_Linear
},
cdHcurve{
    DCT_Linear
},
lhcurve{
    DCT_Linear
},
transmissionCurve{
    FCT_MinMaxCPoints,
    0.00,
    0.50,
    0.35,
    0.35,
    0.60,
    0.75,
    0.35,
    0.35,
    1.00,
    0.50,
    0.35,
    0.35
},
gaintransmissionCurve{
    FCT_MinMaxCPoints,
    0.00,
    0.1,
    0.35,
    0.00,
    0.25,
    0.25,
    0.35,
    0.35,
    0.70,
    0.25,
    0.35,
    0.35,
    1.00,
    0.1,
    0.00,
    0.00
},
mapcurve{
    DCT_Linear
},
str(20),
scal(3),
iter(1),
grad(1),
grads(1),
gam(1.30),
slope(3.),
neigh(80),
offs(0),
highlights(0),
htonalwidth(80),
shadows(0),
stonalwidth(80),
radius(40),
retinexMethod("high"),
retinexcolorspace("Lab"),
gammaretinex("none"),
mapMethod("none"),
viewMethod("none"),
vart(200),
limd(8),
highl(4),
skal(3),
medianmap(false)
{
}

bool RetinexParams::operator ==(const RetinexParams& other) const
{
    return
        enabled == other.enabled
        && cdcurve == other.cdcurve
        && cdHcurve == other.cdHcurve
        && lhcurve == other.lhcurve
        && transmissionCurve == other.transmissionCurve
        && gaintransmissionCurve == other.gaintransmissionCurve
        && mapcurve == other.mapcurve
        && str == other.str
        && scal == other.scal
        && iter == other.iter
        && grad == other.grad
        && grads == other.grads
        && gam == other.gam
        && slope == other.slope
        && neigh == other.neigh
        && offs == other.offs
        && highlights == other.highlights
        && htonalwidth == other.htonalwidth
        && shadows == other.shadows
        && stonalwidth == other.stonalwidth
        && radius == other.radius
        && retinexMethod == other.retinexMethod
        && retinexcolorspace == other.retinexcolorspace
        && gammaretinex == other.gammaretinex
        && mapMethod == other.mapMethod
        && viewMethod == other.viewMethod
        && vart == other.vart
        && limd == other.limd
        && highl == other.highl
        && skal == other.skal
        && medianmap == other.medianmap;
}

bool RetinexParams::operator !=(const RetinexParams& other) const
{
    return !(*this == other);
}

void RetinexParams::getCurves(RetinextransmissionCurve &transmissionCurveLUT, RetinexgaintransmissionCurve &gaintransmissionCurveLUT) const
{
    transmissionCurveLUT.Set(this->transmissionCurve);
    gaintransmissionCurveLUT.Set(this->gaintransmissionCurve);

}

LCurveParams::LCurveParams() :
    enabled(false),
    lcurve{
    DCT_Linear
},
acurve{
    DCT_Linear
},
bcurve{
    DCT_Linear
},
cccurve{
    DCT_Linear
},
chcurve{
    FCT_Linear
},
lhcurve{
    FCT_Linear
},
hhcurve{
    FCT_Linear
},
lccurve{
    DCT_Linear
},
clcurve{
    DCT_Linear
},
brightness(0),
contrast(0),
chromaticity(0),
avoidcolorshift(false),
rstprotection(0),
lcredsk(true)
{
}

bool LCurveParams::operator ==(const LCurveParams& other) const
{
    return
        enabled == other.enabled
        && lcurve == other.lcurve
        && acurve == other.acurve
        && bcurve == other.bcurve
        && cccurve == other.cccurve
        && chcurve == other.chcurve
        && lhcurve == other.lhcurve
        && hhcurve == other.hhcurve
        && lccurve == other.lccurve
        && clcurve == other.clcurve
        && brightness == other.brightness
        && contrast == other.contrast
        && chromaticity == other.chromaticity
        && avoidcolorshift == other.avoidcolorshift
        && rstprotection == other.rstprotection
        && lcredsk == other.lcredsk;
}

bool LCurveParams::operator !=(const LCurveParams& other) const
{
    return !(*this == other);
}

RGBCurvesParams::RGBCurvesParams() :
    enabled(false),
    lumamode(false),
    rcurve{
    DCT_Linear
},
gcurve{
    DCT_Linear
},
bcurve{
    DCT_Linear
}
{
}

bool RGBCurvesParams::operator ==(const RGBCurvesParams& other) const
{
    return
        enabled == other.enabled
        && lumamode == other.lumamode
        && rcurve == other.rcurve
        && gcurve == other.gcurve
        && bcurve == other.bcurve;
}

bool RGBCurvesParams::operator !=(const RGBCurvesParams& other) const
{
    return !(*this == other);
}


LocalContrastParams::LocalContrastParams():
    enabled(false),
    radius(80),
    amount(0.2),
    darkness(1.0),
    lightness(1.0)
{
}


bool LocalContrastParams::operator==(const LocalContrastParams &other) const
{
    return
        enabled == other.enabled
        && radius == other.radius
        && amount == other.amount
        && darkness == other.darkness
        && lightness == other.lightness;
}


bool LocalContrastParams::operator!=(const LocalContrastParams &other) const
{
    return !(*this == other);
}


const double ColorToningParams::LABGRID_CORR_MAX = 12000.f;

ColorToningParams::ColorToningParams() :
    enabled(false),
    autosat(true),
    opacityCurve{
    FCT_MinMaxCPoints,
    0.00,
    0.3,
    0.35,
    0.00,
    0.25,
    0.8,
    0.35,
    0.35,
    0.70,
    0.8,
    0.35,
    0.35,
    1.00,
    0.3,
    0.00,
    0.00
},
colorCurve{
    FCT_MinMaxCPoints,
    0.050,
    0.62,
    0.25,
    0.25,
    0.585,
    0.11,
    0.25,
    0.25
},
satProtectionThreshold(30),
                       saturatedOpacity(80),
                       strength(50),
                       balance(0),
                       hlColSat(60, 80, false),
                       shadowsColSat(80, 208, false),
clcurve{
    DCT_NURBS,
    0.00,
    0.00,
    0.35,
    0.65,
    1.00,
    1.00
},
cl2curve{
    DCT_NURBS,
    0.00,
    0.00,
    0.35,
    0.65,
    1.00,
    1.00
},
method("Lab"),
twocolor("Std"),
redlow(0.0),
greenlow(0.0),
bluelow(0.0),
redmed(0.0),
greenmed(0.0),
bluemed(0.0),
redhigh(0.0),
greenhigh(0.0),
bluehigh(0.0),
satlow(0.0),
sathigh(0.0),
lumamode(true),
labgridALow(0.0),
labgridBLow(0.0),
labgridAHigh(0.0),
labgridBHigh(0.0)
{
}

bool ColorToningParams::operator ==(const ColorToningParams& other) const
{
    return
        enabled == other.enabled
        && autosat == other.autosat
        && opacityCurve == other.opacityCurve
        && colorCurve == other.colorCurve
        && satProtectionThreshold == other.satProtectionThreshold
        && saturatedOpacity == other.saturatedOpacity
        && strength == other.strength
        && balance == other.balance
        && hlColSat == other.hlColSat
        && shadowsColSat == other.shadowsColSat
        && clcurve == other.clcurve
        && cl2curve == other.cl2curve
        && method == other.method
        && twocolor == other.twocolor
        && redlow == other.redlow
        && greenlow == other.greenlow
        && bluelow == other.bluelow
        && redmed == other.redmed
        && greenmed == other.greenmed
        && bluemed == other.bluemed
        && redhigh == other.redhigh
        && greenhigh == other.greenhigh
        && bluehigh == other.bluehigh
        && satlow == other.satlow
        && sathigh == other.sathigh
        && lumamode == other.lumamode
        && labgridALow == other.labgridALow
        && labgridBLow == other.labgridBLow
        && labgridAHigh == other.labgridAHigh
        && labgridBHigh == other.labgridBHigh;
}

bool ColorToningParams::operator !=(const ColorToningParams& other) const
{
    return !(*this == other);
}

void ColorToningParams::mixerToCurve(std::vector<double>& colorCurve, std::vector<double>& opacityCurve) const
{
    // check if non null first
    if (!redlow && !greenlow && !bluelow && !redmed && !greenmed && !bluemed && !redhigh && !greenhigh && !bluehigh) {
        colorCurve.resize(1);
        colorCurve.at(0) = FCT_Linear;
        opacityCurve.resize(1);
        opacityCurve.at(0) = FCT_Linear;
        return;
    }

    float low[3]; // RGB color for shadows
    float med[3]; // RGB color for mid-tones
    float high[3]; // RGB color for highlights
    float lowSat = 0.f;
    float medSat = 0.f;
    float highSat = 0.f;
    float minTmp, maxTmp;

// Fill the shadow mixer values of the Color TOning tool
    low[0] = float (redlow) / 100.f;  // [-1. ; +1.]
    low[1] = float (greenlow) / 100.f; // [-1. ; +1.]
    low[2] = float (bluelow) / 100.f;  // [-1. ; +1.]
    minTmp = min<float> (low[0], low[1], low[2]);
    maxTmp = max<float> (low[0], low[1], low[2]);

    if (maxTmp - minTmp > 0.005f) {
        float v[3];
        lowSat = (maxTmp - minTmp) / 2.f;

        if (low[0] == minTmp) {
            v[0] = 0.f;
        } else if (low[1] == minTmp) {
            v[1] = 0.f;
        } else if (low[2] == minTmp) {
            v[2] = 0.f;
        }

        if (low[0] == maxTmp) {
            v[0] = 1.f;
        } else if (low[1] == maxTmp) {
            v[1] = 1.f;
        } else if (low[2] == maxTmp) {
            v[2] = 1.f;
        }

        if (low[0] != minTmp && low[0] != maxTmp) {
            v[0] = (low[0] - minTmp) / (maxTmp - minTmp);
        } else if (low[1] != minTmp && low[1] != maxTmp) {
            v[1] = (low[1] - minTmp) / (maxTmp - minTmp);
        } else if (low[2] != minTmp && low[2] != maxTmp) {
            v[2] = (low[2] - minTmp) / (maxTmp - minTmp);
        }

        low[0] = v[0];
        low[1] = v[1];
        low[2] = v[2];
    } else {
        low[0] = low[1] = low[2] = 1.f;
    }

// Fill the mid-tones mixer values of the Color TOning tool
    med[0] = float (redmed) / 100.f;  // [-1. ; +1.]
    med[1] = float (greenmed) / 100.f; // [-1. ; +1.]
    med[2] = float (bluemed) / 100.f;  // [-1. ; +1.]
    minTmp = min<float> (med[0], med[1], med[2]);
    maxTmp = max<float> (med[0], med[1], med[2]);

    if (maxTmp - minTmp > 0.005f) {
        float v[3];
        medSat = (maxTmp - minTmp) / 2.f;

        if (med[0] == minTmp) {
            v[0] = 0.f;
        } else if (med[1] == minTmp) {
            v[1] = 0.f;
        } else if (med[2] == minTmp) {
            v[2] = 0.f;
        }

        if (med[0] == maxTmp) {
            v[0] = 1.f;
        } else if (med[1] == maxTmp) {
            v[1] = 1.f;
        } else if (med[2] == maxTmp) {
            v[2] = 1.f;
        }

        if (med[0] != minTmp && med[0] != maxTmp) {
            v[0] = (med[0] - minTmp) / (maxTmp - minTmp);
        } else if (med[1] != minTmp && med[1] != maxTmp) {
            v[1] = (med[1] - minTmp) / (maxTmp - minTmp);
        } else if (med[2] != minTmp && med[2] != maxTmp) {
            v[2] = (med[2] - minTmp) / (maxTmp - minTmp);
        }

        med[0] = v[0];
        med[1] = v[1];
        med[2] = v[2];
    } else {
        med[0] = med[1] = med[2] = 1.f;
    }

    // Fill the highlight mixer values of the Color TOning tool
    high[0] = float (redhigh) / 100.f;   // [-1. ; +1.]
    high[1] = float (greenhigh) / 100.f; // [-1. ; +1.]
    high[2] = float (bluehigh) / 100.f;  // [-1. ; +1.]
    minTmp = min<float> (high[0], high[1], high[2]);
    maxTmp = max<float> (high[0], high[1], high[2]);

    if (maxTmp - minTmp > 0.005f) {
        float v[3];
        highSat = (maxTmp - minTmp) / 2.f;

        if (high[0] == minTmp) {
            v[0] = 0.f;
        } else if (high[1] == minTmp) {
            v[1] = 0.f;
        } else if (high[2] == minTmp) {
            v[2] = 0.f;
        }

        if (high[0] == maxTmp) {
            v[0] = 1.f;
        } else if (high[1] == maxTmp) {
            v[1] = 1.f;
        } else if (high[2] == maxTmp) {
            v[2] = 1.f;
        }

        if (high[0] != minTmp && high[0] != maxTmp) {
            v[0] = (high[0] - minTmp) / (maxTmp - minTmp);
        } else if (high[1] != minTmp && high[1] != maxTmp) {
            v[1] = (high[1] - minTmp) / (maxTmp - minTmp);
        } else if (high[2] != minTmp && high[2] != maxTmp) {
            v[2] = (high[2] - minTmp) / (maxTmp - minTmp);
        }

        high[0] = v[0];
        high[1] = v[1];
        high[2] = v[2];
    } else {
        high[0] = high[1] = high[2] = 1.f;
    }

    const double xPosLow  = 0.1;
    const double xPosMed  = 0.4;
    const double xPosHigh = 0.7;

    colorCurve.resize(medSat != 0.f ? 13 : 9);
    colorCurve.at(0) = FCT_MinMaxCPoints;
    opacityCurve.resize(13);
    opacityCurve.at(0) = FCT_MinMaxCPoints;

    float h, s, l;
    int idx = 1;

    if (lowSat == 0.f) {
        if (medSat != 0.f) {
            Color::rgb2hsl(med[0], med[1], med[2], h, s, l);
        } else { // highSat can't be null if the 2 other ones are!
            Color::rgb2hsl(high[0], high[1], high[2], h, s, l);
        }
    } else {
        Color::rgb2hsl(low[0], low[1], low[2], h, s, l);
    }

    colorCurve.at(idx++) = xPosLow;
    colorCurve.at(idx++) = h;
    colorCurve.at(idx++) = 0.35;
    colorCurve.at(idx++) = 0.35;

    if (medSat != 0.f) {
        Color::rgb2hsl(med[0], med[1], med[2], h, s, l);
        colorCurve.at(idx++) = xPosMed;
        colorCurve.at(idx++) = h;
        colorCurve.at(idx++) = 0.35;
        colorCurve.at(idx++) = 0.35;
    }

    if (highSat == 0.f) {
        if (medSat != 0.f) {
            Color::rgb2hsl(med[0], med[1], med[2], h, s, l);
        } else { // lowSat can't be null if the 2 other ones are!
            Color::rgb2hsl(low[0], low[1], low[2], h, s, l);
        }
    } else {
        Color::rgb2hsl(high[0], high[1], high[2], h, s, l);
    }

    colorCurve.at(idx++) = xPosHigh;
    colorCurve.at(idx++) = h;
    colorCurve.at(idx++) = 0.35;
    colorCurve.at(idx)   = 0.35;

    opacityCurve.at(1)  = xPosLow;
    opacityCurve.at(2)  = double (lowSat);
    opacityCurve.at(3)  = 0.35;
    opacityCurve.at(4)  = 0.35;
    opacityCurve.at(5)  = xPosMed;
    opacityCurve.at(6)  = double (medSat);
    opacityCurve.at(7)  = 0.35;
    opacityCurve.at(8)  = 0.35;
    opacityCurve.at(9)  = xPosHigh;
    opacityCurve.at(10) = double (highSat);
    opacityCurve.at(11) = 0.35;
    opacityCurve.at(12) = 0.35;
}

void ColorToningParams::slidersToCurve(std::vector<double>& colorCurve, std::vector<double>& opacityCurve) const
{
    if (hlColSat.getBottom() == 0 && shadowsColSat.getBottom() == 0) { // if both opacity are null, set both curves to Linear
        colorCurve.resize(1);
        colorCurve.at(0) = FCT_Linear;
        opacityCurve.resize(1);
        opacityCurve.at(0) = FCT_Linear;
        return;
    }

    colorCurve.resize(9);
    colorCurve.at(0) = FCT_MinMaxCPoints;
    colorCurve.at(1) = 0.26 + 0.12 * double (balance) / 100.;
    colorCurve.at(2) = double (shadowsColSat.getTop()) / 360.;
    colorCurve.at(3) = 0.35;
    colorCurve.at(4) = 0.35;
    colorCurve.at(5) = 0.64 + 0.12 * double (balance) / 100.;
    colorCurve.at(6) = double (hlColSat.getTop()) / 360.;
    colorCurve.at(7) = 0.35;
    colorCurve.at(8) = 0.35;

    opacityCurve.resize(9);
    opacityCurve.at(0) = FCT_MinMaxCPoints;
    opacityCurve.at(1) = colorCurve.at(1);
    opacityCurve.at(2) = double (shadowsColSat.getBottom()) / 100.;
    opacityCurve.at(3) = 0.35;
    opacityCurve.at(4) = 0.35;
    opacityCurve.at(5) = colorCurve.at(5);
    opacityCurve.at(6) = double (hlColSat.getBottom()) / 100.;
    opacityCurve.at(7) = 0.35;
    opacityCurve.at(8) = 0.35;
}

void ColorToningParams::getCurves(ColorGradientCurve& colorCurveLUT, OpacityCurve& opacityCurveLUT, const double xyz_rgb[3][3], bool& opautili) const
{
    float satur = 0.8f;
    float lumin = 0.5f; //middle of luminance for optimization of gamut - no real importance...as we work in XYZ and gamut control

    // Transform slider values to control points
    std::vector<double> cCurve, oCurve;

    if (method == "RGBSliders" || method == "Splitlr") {
        slidersToCurve(cCurve, oCurve);
    } else if (method == "Splitco") {
        mixerToCurve(cCurve, oCurve);
    } else {
        cCurve = this->colorCurve;
        oCurve = this->opacityCurve;
    }

    if (method == "Lab") {
        if (twocolor == "Separ") {
            satur = 0.9f;
        }

        if (twocolor == "All" || twocolor == "Two") {
            satur = 0.9f;
        }

        colorCurveLUT.SetXYZ(cCurve, xyz_rgb, satur, lumin);
        opacityCurveLUT.Set(oCurve, opautili);
    } else if (method == "Splitlr" || method == "Splitco") {
        colorCurveLUT.SetXYZ(cCurve, xyz_rgb, satur, lumin);
        opacityCurveLUT.Set(oCurve, opautili);
    } else if (method.substr(0, 3) == "RGB") {
        colorCurveLUT.SetRGB(cCurve);
        opacityCurveLUT.Set(oCurve, opautili);
    }
}

SharpeningParams::SharpeningParams() :
    enabled(false),
    radius(0.5),
    amount(200),
    threshold(20, 80, 2000, 1200, false),
    edgesonly(false),
    edges_radius(1.9),
    edges_tolerance(1800),
    halocontrol(false),
    halocontrol_amount(85),
    method("usm"),
    deconvamount(75),
    deconvradius(0.75),
    deconviter(30),
    deconvdamping(20)
{
}

bool SharpeningParams::operator ==(const SharpeningParams& other) const
{
    return
        enabled == other.enabled
        && radius == other.radius
        && amount == other.amount
        && threshold == other.threshold
        && edgesonly == other.edgesonly
        && edges_radius == other.edges_radius
        && edges_tolerance == other.edges_tolerance
        && halocontrol == other.halocontrol
        && halocontrol_amount == other.halocontrol_amount
        && method == other.method
        && deconvamount == other.deconvamount
        && deconvradius == other.deconvradius
        && deconviter == other.deconviter
        && deconvdamping == other.deconvdamping;
}

bool SharpeningParams::operator !=(const SharpeningParams& other) const
{
    return !(*this == other);
}

SharpenEdgeParams::SharpenEdgeParams() :
    enabled(false),
    passes(2),
    amount(50.0),
    threechannels(false)
{
}

bool SharpenEdgeParams::operator ==(const SharpenEdgeParams& other) const
{
    return
        enabled == other.enabled
        && passes == other.passes
        && amount == other.amount
        && threechannels == other.threechannels;
}

bool SharpenEdgeParams::operator !=(const SharpenEdgeParams& other) const
{
    return !(*this == other);
}

SharpenMicroParams::SharpenMicroParams() :
    enabled(false),
    matrix(false),
    amount(20.0),
    uniformity(50.0)
{
}

bool SharpenMicroParams::operator ==(const SharpenMicroParams& other) const
{
    return
        enabled == other.enabled
        && matrix == other.matrix
        && amount == other.amount
        && uniformity == other.uniformity;
}

bool SharpenMicroParams::operator !=(const SharpenMicroParams& other) const
{
    return !(*this == other);
}

VibranceParams::VibranceParams() :
    enabled(false),
    pastels(0),
    saturated(0),
    psthreshold(0, 75, false),
    protectskins(false),
    avoidcolorshift(true),
    pastsattog(true),
    skintonescurve{
    DCT_Linear
}
{
}

bool VibranceParams::operator ==(const VibranceParams& other) const
{
    return
        enabled == other.enabled
        && pastels == other.pastels
        && saturated == other.saturated
        && psthreshold == other.psthreshold
        && protectskins == other.protectskins
        && avoidcolorshift == other.avoidcolorshift
        && pastsattog == other.pastsattog
        && skintonescurve == other.skintonescurve;
}

bool VibranceParams::operator !=(const VibranceParams& other) const
{
    return !(*this == other);
}

WBParams::WBParams() :
    enabled(true),
    method("Camera"),
    temperature(6504),
    green(1.0),
    equal(1.0),
    tempBias(0.0)
{
}

bool WBParams::operator ==(const WBParams& other) const
{
    return
        enabled == other.enabled
        && method == other.method
        && temperature == other.temperature
        && green == other.green
        && equal == other.equal
        && tempBias == other.tempBias;
}

bool WBParams::operator !=(const WBParams& other) const
{
    return !(*this == other);
}

const std::vector<WBEntry>& WBParams::getWbEntries()
{
    static const std::vector<WBEntry> wb_entries = {
        {"Camera",               WBEntry::Type::CAMERA,      M("TP_WBALANCE_CAMERA"),         0, 1.f,   1.f,   0.f},
        {"Auto",                 WBEntry::Type::AUTO,        M("TP_WBALANCE_AUTO"),           0, 1.f,   1.f,   0.f},
        {"Daylight",             WBEntry::Type::DAYLIGHT,    M("TP_WBALANCE_DAYLIGHT"),    5300, 1.f,   1.f,   0.f},
        {"Cloudy",               WBEntry::Type::CLOUDY,      M("TP_WBALANCE_CLOUDY"),      6200, 1.f,   1.f,   0.f},
        {"Shade",                WBEntry::Type::SHADE,       M("TP_WBALANCE_SHADE"),       7600, 1.f,   1.f,   0.f},
        {"Water 1",              WBEntry::Type::WATER,       M("TP_WBALANCE_WATER1"),     35000, 0.3f,  1.1f,  0.f},
        {"Water 2",              WBEntry::Type::WATER,       M("TP_WBALANCE_WATER2"),     48000, 0.63f, 1.38f, 0.f},
        {"Tungsten",             WBEntry::Type::TUNGSTEN,    M("TP_WBALANCE_TUNGSTEN"),    2856, 1.f,   1.f,   0.f},
        {"Fluo F1",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO1"),       6430, 1.f,   1.f,   0.f},
        {"Fluo F2",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO2"),       4230, 1.f,   1.f,   0.f},
        {"Fluo F3",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO3"),       3450, 1.f,   1.f,   0.f},
        {"Fluo F4",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO4"),       2940, 1.f,   1.f,   0.f},
        {"Fluo F5",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO5"),       6350, 1.f,   1.f,   0.f},
        {"Fluo F6",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO6"),       4150, 1.f,   1.f,   0.f},
        {"Fluo F7",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO7"),       6500, 1.f,   1.f,   0.f},
        {"Fluo F8",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO8"),       5020, 1.f,   1.f,   0.f},
        {"Fluo F9",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO9"),       4330, 1.f,   1.f,   0.f},
        {"Fluo F10",             WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO10"),      5300, 1.f,   1.f,   0.f},
        {"Fluo F11",             WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO11"),      4000, 1.f,   1.f,   0.f},
        {"Fluo F12",             WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO12"),      3000, 1.f,   1.f,   0.f},
        {"HMI Lamp",             WBEntry::Type::LAMP,        M("TP_WBALANCE_HMI"),         4800, 1.f,   1.f,   0.f},
        {"GTI Lamp",             WBEntry::Type::LAMP,        M("TP_WBALANCE_GTI"),         5000, 1.f,   1.f,   0.f},
        {"JudgeIII Lamp",        WBEntry::Type::LAMP,        M("TP_WBALANCE_JUDGEIII"),    5100, 1.f,   1.f,   0.f},
        {"Solux Lamp 3500K",     WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX35"),     3480, 1.f,   1.f,   0.f},
        {"Solux Lamp 4100K",     WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX41"),     3930, 1.f,   1.f,   0.f},
        {"Solux Lamp 4700K",     WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX47"),     4700, 1.f,   1.f,   0.f},
        {"NG Solux Lamp 4700K",  WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX47_NG"),  4480, 1.f,   1.f,   0.f},
        {"LED LSI Lumelex 2040", WBEntry::Type::LED,         M("TP_WBALANCE_LED_LSI"),     2970, 1.f,   1.f,   0.f},
        {"LED CRS SP12 WWMR16",  WBEntry::Type::LED,         M("TP_WBALANCE_LED_CRS"),     3050, 1.f,   1.f,   0.f},
        {"Flash 5500K",          WBEntry::Type::FLASH,       M("TP_WBALANCE_FLASH55"),     5500, 1.f,   1.f,   0.f},
        {"Flash 6000K",          WBEntry::Type::FLASH,       M("TP_WBALANCE_FLASH60"),     6000, 1.f,   1.f,   0.f},
        {"Flash 6500K",          WBEntry::Type::FLASH,       M("TP_WBALANCE_FLASH65"),     6500, 1.f,   1.f,   0.f},
        // Should remain the last one
        {"Custom",               WBEntry::Type::CUSTOM,      M("TP_WBALANCE_CUSTOM"),        0, 1.f,   1.f,   0.f}
    };

    return wb_entries;
}

ColorAppearanceParams::ColorAppearanceParams() :
    enabled(false),
    degree(90),
    autodegree(true),
    degreeout(90),
    autodegreeout(true),
    curve{
    DCT_Linear
},
curve2{
    DCT_Linear
},
curve3{
    DCT_Linear
},
curveMode(TcMode::LIGHT),
curveMode2(TcMode::LIGHT),
curveMode3(CtcMode::CHROMA),
surround("Average"),
surrsrc("Average"),
adapscen(2000.0),
autoadapscen(true),
ybscen(18),
autoybscen(true),
adaplum(16),
badpixsl(0),
wbmodel("RawT"),
algo("No"),
contrast(0.0),
qcontrast(0.0),
jlight(0.0),
qbright(0.0),
chroma(0.0),
schroma(0.0),
mchroma(0.0),
colorh(0.0),
rstprotection(0.0),
surrsource(false),
gamut(true),
datacie(false),
tonecie(false),
tempout(5000),
ybout(18),
greenout(1.0),
tempsc(5000),
greensc(1.0)
{
}

bool ColorAppearanceParams::operator ==(const ColorAppearanceParams& other) const
{
    return
        enabled == other.enabled
        && degree == other.degree
        && autodegree == other.autodegree
        && degreeout == other.degreeout
        && autodegreeout == other.autodegreeout
        && curve == other.curve
        && curve2 == other.curve2
        && curve3 == other.curve3
        && curveMode == other.curveMode
        && curveMode2 == other.curveMode2
        && curveMode3 == other.curveMode3
        && surround == other.surround
        && surrsrc == other.surrsrc
        && adapscen == other.adapscen
        && autoadapscen == other.autoadapscen
        && ybscen == other.ybscen
        && autoybscen == other.autoybscen
        && adaplum == other.adaplum
        && badpixsl == other.badpixsl
        && wbmodel == other.wbmodel
        && algo == other.algo
        && contrast == other.contrast
        && qcontrast == other.qcontrast
        && jlight == other.jlight
        && qbright == other.qbright
        && chroma == other.chroma
        && schroma == other.schroma
        && mchroma == other.mchroma
        && colorh == other.colorh
        && rstprotection == other.rstprotection
        && surrsource == other.surrsource
        && gamut == other.gamut
        && datacie == other.datacie
        && tonecie == other.tonecie
        && tempout == other.tempout
        && ybout == other.ybout
        && greenout == other.greenout
        && tempsc == other.tempsc
        && greensc == other.greensc;
}

bool ColorAppearanceParams::operator !=(const ColorAppearanceParams& other) const
{
    return !(*this == other);
}

DefringeParams::DefringeParams() :
    enabled(false),
    radius(2.0),
    threshold(13),
    huecurve{
    FCT_MinMaxCPoints,
    0.166666667,
    0.,
    0.35,
    0.35,
    0.347,
    0.,
    0.35,
    0.35,
    0.513667426,
    0,
    0.35,
    0.35,
    0.668944571,
    0.,
    0.35,
    0.35,
    0.8287775246,
    0.97835991,
    0.35,
    0.35,
    0.9908883827,
    0.,
    0.35,
    0.35
}
{
}

bool DefringeParams::operator ==(const DefringeParams& other) const
{
    return
        enabled == other.enabled
        && radius == other.radius
        && threshold == other.threshold
        && huecurve == other.huecurve;
}

bool DefringeParams::operator !=(const DefringeParams& other) const
{
    return !(*this == other);
}

ImpulseDenoiseParams::ImpulseDenoiseParams() :
    enabled(false),
    thresh(50)
{
}

bool ImpulseDenoiseParams::operator ==(const ImpulseDenoiseParams& other) const
{
    return
        enabled == other.enabled
        && thresh == other.thresh;
}

bool ImpulseDenoiseParams::operator !=(const ImpulseDenoiseParams& other) const
{
    return !(*this == other);
}

DirPyrDenoiseParams::DirPyrDenoiseParams() :
    lcurve{
    FCT_MinMaxCPoints,
    0.05,
    0.15,
    0.35,
    0.35,
    0.55,
    0.04,
    0.35,
    0.35
},
cccurve{
    FCT_MinMaxCPoints,
    0.05,
    0.50,
    0.35,
    0.35,
    0.35,
    0.05,
    0.35,
    0.35
},
enabled(false),
enhance(false),
median(false),
perform(false),
luma(0),
Ldetail(0),
chroma(15),
redchro(0),
bluechro(0),
gamma(1.7),
dmethod("Lab"),
Lmethod("SLI"),
Cmethod("MAN"),
C2method("AUTO"),
smethod("shal"),
medmethod("soft"),
methodmed("none"),
rgbmethod("soft"),
passes(1)
{
}

bool DirPyrDenoiseParams::operator ==(const DirPyrDenoiseParams& other) const
{
    return
        lcurve == other.lcurve
        && cccurve == other.cccurve
        && enabled == other.enabled
        && enhance == other.enhance
        && median == other.median
        && perform == other.perform
        && luma == other.luma
        && Ldetail == other.Ldetail
        && chroma == other.chroma
        && redchro == other.redchro
        && bluechro == other.bluechro
        && gamma == other.gamma
        && dmethod == other.dmethod
        && Lmethod == other.Lmethod
        && Cmethod == other.Cmethod
        && C2method == other.C2method
        && smethod == other.smethod
        && medmethod == other.medmethod
        && methodmed == other.methodmed
        && rgbmethod == other.rgbmethod
        && passes == other.passes;
}

bool DirPyrDenoiseParams::operator !=(const DirPyrDenoiseParams& other) const
{
    return !(*this == other);
}

void DirPyrDenoiseParams::getCurves(NoiseCurve &lCurve, NoiseCurve &cCurve) const
{
    lCurve.Set(this->lcurve);
    cCurve.Set(this->cccurve);
}

EPDParams::EPDParams() :
    enabled(false),
    strength(0.5),
    gamma(1.0),
    edgeStopping(1.4),
    scale(1.0),
    reweightingIterates(0)
{
}

bool EPDParams::operator ==(const EPDParams& other) const
{
    return
        enabled == other.enabled
        && strength == other.strength
        && gamma == other.gamma
        && edgeStopping == other.edgeStopping
        && scale == other.scale
        && reweightingIterates == other.reweightingIterates;
}

bool EPDParams::operator !=(const EPDParams& other) const
{
    return !(*this == other);
}

FattalToneMappingParams::FattalToneMappingParams() :
    enabled(false),
    threshold(0),
    amount(30),
    anchor(50)
{
}

bool FattalToneMappingParams::operator ==(const FattalToneMappingParams& other) const
{
    return
        enabled == other.enabled
        && threshold == other.threshold
        && amount == other.amount
        && anchor == other.anchor;
}

bool FattalToneMappingParams::operator !=(const FattalToneMappingParams& other) const
{
    return !(*this == other);
}

SHParams::SHParams() :
    enabled(false),
    highlights(0),
    htonalwidth(70),
    shadows(0),
    stonalwidth(30),
    radius(40)
{
}

bool SHParams::operator ==(const SHParams& other) const
{
    return
        enabled == other.enabled
        && highlights == other.highlights
        && htonalwidth == other.htonalwidth
        && shadows == other.shadows
        && stonalwidth == other.stonalwidth
        && radius == other.radius;
}

bool SHParams::operator !=(const SHParams& other) const
{
    return !(*this == other);
}

CropParams::CropParams() :
    enabled(false),
    x(-1),
    y(-1),
    w(15000),
    h(15000),
    fixratio(true),
    ratio("As Image"),
    orientation("As Image"),
    guide("Frame")
{
}

bool CropParams::operator ==(const CropParams& other) const
{
    return
        enabled == other.enabled
        && x == other.x
        && y == other.y
        && w == other.w
        && h == other.h
        && fixratio == other.fixratio
        && ratio == other.ratio
        && orientation == other.orientation
        && guide == other.guide;
}

bool CropParams::operator !=(const CropParams& other) const
{
    return !(*this == other);
}

void CropParams::mapToResized(int resizedWidth, int resizedHeight, int scale, int& x1, int& x2, int& y1, int& y2) const
{
    x1 = 0, x2 = resizedWidth, y1 = 0, y2 = resizedHeight;

    if (enabled) {
        x1 = min(resizedWidth - 1, max(0, x / scale));
        y1 = min(resizedHeight - 1, max(0, y / scale));
        x2 = min(resizedWidth, max(0, (x + w) / scale));
        y2 = min(resizedHeight, max(0, (y + h) / scale));
    }
}

CoarseTransformParams::CoarseTransformParams() :
    rotate(0),
    hflip(false),
    vflip(false)
{
}

bool CoarseTransformParams::operator ==(const CoarseTransformParams& other) const
{
    return
        rotate == other.rotate
        && hflip == other.hflip
        && vflip == other.vflip;
}

bool CoarseTransformParams::operator !=(const CoarseTransformParams& other) const
{
    return !(*this == other);
}

CommonTransformParams::CommonTransformParams() :
    autofill(true)
{
}

bool CommonTransformParams::operator ==(const CommonTransformParams& other) const
{
    return autofill == other.autofill;
}

bool CommonTransformParams::operator !=(const CommonTransformParams& other) const
{
    return !(*this == other);
}

RotateParams::RotateParams() :
    degree(0.0)
{
}

bool RotateParams::operator ==(const RotateParams& other) const
{
    return degree == other.degree;
}

bool RotateParams::operator !=(const RotateParams& other) const
{
    return !(*this == other);
}

DistortionParams::DistortionParams() :
    amount(0.0)
{
}

bool DistortionParams::operator ==(const DistortionParams& other) const
{
    return amount == other.amount;
}

bool DistortionParams::operator !=(const DistortionParams& other) const
{
    return !(*this == other);
}

LensProfParams::LensProfParams() :
    lcMode(LcMode::NONE),
    useDist(true),
    useVign(true),
    useCA(false)
{
}

bool LensProfParams::operator ==(const LensProfParams& other) const
{
    return
        lcMode == other.lcMode
        && lcpFile == other.lcpFile
        && useCA == other.useCA
        && lfCameraMake == other.lfCameraMake
        && lfCameraModel == other.lfCameraModel
        && lfLens == other.lfLens;
}

bool LensProfParams::operator !=(const LensProfParams& other) const
{
    return !(*this == other);
}

bool LensProfParams::useLensfun() const
{
    return lcMode == LcMode::LENSFUNAUTOMATCH || lcMode == LcMode::LENSFUNMANUAL;
}

bool LensProfParams::lfAutoMatch() const
{
    return lcMode == LcMode::LENSFUNAUTOMATCH;
}

bool LensProfParams::useLcp() const
{
    return lcMode == LcMode::LCP && lcpFile.length() > 0;
}

bool LensProfParams::lfManual() const
{
    return lcMode == LcMode::LENSFUNMANUAL;
}

const std::vector<const char*>& LensProfParams::getMethodStrings() const
{
    static const std::vector<const char*> method_strings = {
        "none",
        "lfauto",
        "lfmanual",
        "lcp"
    };
    return method_strings;
}

Glib::ustring LensProfParams::getMethodString(LcMode mode) const
{
    return getMethodStrings()[toUnderlying(mode)];
}

LensProfParams::LcMode LensProfParams::getMethodNumber(const Glib::ustring& mode) const
{
    for (std::vector<const char*>::size_type i = 0; i < getMethodStrings().size(); ++i) {
        if (getMethodStrings()[i] == mode) {
            return static_cast<LcMode>(i);
        }
    }

    return LcMode::NONE;
}

PerspectiveParams::PerspectiveParams() :
    horizontal(0.0),
    vertical(0.0)
{
}

bool PerspectiveParams::operator ==(const PerspectiveParams& other) const
{
    return
        horizontal == other.horizontal
        && vertical == other.vertical;
}

bool PerspectiveParams::operator !=(const PerspectiveParams& other) const
{
    return !(*this == other);
}

GradientParams::GradientParams() :
    enabled(false),
    degree(0.0),
    feather(25),
    strength(0.60),
    centerX(0),
    centerY(0)
{
}

bool GradientParams::operator ==(const GradientParams& other) const
{
    return
        enabled == other.enabled
        && degree == other.degree
        && feather == other.feather
        && strength == other.strength
        && centerX == other.centerX
        && centerY == other.centerY;
}

bool GradientParams::operator !=(const GradientParams& other) const
{
    return !(*this == other);
}

PCVignetteParams::PCVignetteParams() :
    enabled(false),
    strength(0.60),
    feather(50),
    roundness(50)
{
}

bool PCVignetteParams::operator ==(const PCVignetteParams& other) const
{
    return
        enabled == other.enabled
        && strength == other.strength
        && feather == other.feather
        && roundness == other.roundness;
}

bool PCVignetteParams::operator !=(const PCVignetteParams& other) const
{
    return !(*this == other);
}

VignettingParams::VignettingParams() :
    amount(0),
    radius(50),
    strength(1),
    centerX(0),
    centerY(0)
{
}

bool VignettingParams::operator ==(const VignettingParams& other) const
{
    return
        amount == other.amount
        && radius == other.radius
        && strength == other.strength
        && centerX == other.centerX
        && centerY == other.centerY;
}

bool VignettingParams::operator !=(const VignettingParams& other) const
{
    return !(*this == other);
}

ChannelMixerParams::ChannelMixerParams() :
    enabled(false),
    red{
    100,
    0,
    0
},
green{
    0,
    100,
    0
},
blue{
    0,
    0,
    100
}
{
}

bool ChannelMixerParams::operator ==(const ChannelMixerParams& other) const
{
    if (enabled != other.enabled) {
        return false;
    }

    for (unsigned int i = 0; i < 3; ++i) {
        if (
            red[i] != other.red[i]
            || green[i] != other.green[i]
            || blue[i] != other.blue[i]
        ) {
            return false;
        }
    }

    return true;
}

bool ChannelMixerParams::operator !=(const ChannelMixerParams& other) const
{
    return !(*this == other);
}

BlackWhiteParams::BlackWhiteParams() :
    beforeCurve{
    DCT_Linear
},
beforeCurveMode(BlackWhiteParams::TcMode::STD_BW),
afterCurve{
    DCT_Linear
},
afterCurveMode(BlackWhiteParams::TcMode::STD_BW),
algo("SP"),
luminanceCurve{
    FCT_Linear
},
autoc(false),
enabledcc(true),
enabled(false),
filter("None"),
setting("NormalContrast"),
method("Desaturation"),
mixerRed(33),
mixerOrange(33),
mixerYellow(33),
mixerGreen(33),
mixerCyan(33),
mixerBlue(33),
mixerMagenta(33),
mixerPurple(33),
gammaRed(0),
gammaGreen(0),
gammaBlue(0)
{
}

bool BlackWhiteParams::operator ==(const BlackWhiteParams& other) const
{
    return
        beforeCurve == other.beforeCurve
        && beforeCurveMode == other.beforeCurveMode
        && afterCurve == other.afterCurve
        && afterCurveMode == other.afterCurveMode
        && algo == other.algo
        && luminanceCurve == other.luminanceCurve
        && autoc == other.autoc
        && enabledcc == other.enabledcc
        && enabled == other.enabled
        && filter == other.filter
        && setting == other.setting
        && method == other.method
        && mixerRed == other.mixerRed
        && mixerOrange == other.mixerOrange
        && mixerYellow == other.mixerYellow
        && mixerGreen == other.mixerGreen
        && mixerCyan == other.mixerCyan
        && mixerBlue == other.mixerBlue
        && mixerMagenta == other.mixerMagenta
        && mixerPurple == other.mixerPurple
        && gammaRed == other.gammaRed
        && gammaGreen == other.gammaGreen
        && gammaBlue == other.gammaBlue;
}

bool BlackWhiteParams::operator !=(const BlackWhiteParams& other) const
{
    return !(*this == other);
}

CACorrParams::CACorrParams() :
    red(0.0),
    blue(0.0)
{
}

bool CACorrParams::operator ==(const CACorrParams& other) const
{
    return
        red == other.red
        && blue == other.blue;
}

bool CACorrParams::operator !=(const CACorrParams& other) const
{
    return !(*this == other);
}

ResizeParams::ResizeParams() :
    enabled(false),
    scale(1.0),
    appliesTo("Cropped area"),
    method("Lanczos"),
    dataspec(3),
    width(900),
    height(900)
{
}

bool ResizeParams::operator ==(const ResizeParams& other) const
{
    return
        enabled == other.enabled
        && scale == other.scale
        && appliesTo == other.appliesTo
        && method == other.method
        && dataspec == other.dataspec
        && width == other.width
        && height == other.height;
}

bool ResizeParams::operator !=(const ResizeParams& other) const
{
    return !(*this == other);
}

const Glib::ustring ColorManagementParams::NoICMString = Glib::ustring("No ICM: sRGB output");

ColorManagementParams::ColorManagementParams() :
    input("(cameraICC)"),
    toneCurve(false),
    applyLookTable(false),
    applyBaselineExposureOffset(true),
    applyHueSatMap(true),
    dcpIlluminant(0),
    working("ProPhoto"),
    output("RT_sRGB"),
    outputIntent(RI_RELATIVE),
    outputBPC(true),
    gamma("Free"),
    gampos(2.4),
    slpos(12.92310),
	predx(0.6400),
	predy(0.3300),
	pgrex(0.3000),
	pgrey(0.6000),
	pblux(0.1500),
	pbluy(0.0600),
    gamm(2.4),
    slop(12.92),
	
    wprimaries("sRGB"),
    wprofile("none"),
    wtemp("DEF"),
    freegamma(false),
    wtrcin("none")
	
{
}

bool ColorManagementParams::operator ==(const ColorManagementParams& other) const
{
    return
        input == other.input
        && toneCurve == other.toneCurve
        && applyLookTable == other.applyLookTable
        && applyBaselineExposureOffset == other.applyBaselineExposureOffset
        && applyHueSatMap == other.applyHueSatMap
        && dcpIlluminant == other.dcpIlluminant
        && working == other.working
        && output == other.output
        && outputIntent == other.outputIntent
        && outputBPC == other.outputBPC
        && gamma == other.gamma
        && gampos == other.gampos
        && slpos == other.slpos
        && predx == other.predx
        && predy == other.predy
        && pgrex == other.pgrex
        && pgrey == other.pgrey
        && pblux == other.pblux
        && pbluy == other.pbluy
        && gamm == other.gamm
        && slop == other.slop
        && wprimaries == other.wprimaries
        && wprofile == other.wprofile
        && wtemp == other.wtemp
        && wtrcin == other.wtrcin
        && freegamma == other.freegamma;
}

bool ColorManagementParams::operator !=(const ColorManagementParams& other) const
{
    return !(*this == other);
}

WaveletParams::WaveletParams() :
    ccwcurve{
    static_cast<double>(FCT_MinMaxCPoints),
    0.0,
    0.25,
    0.35,
    0.35,
    0.50,
    0.75,
    0.35,
    0.35,
    0.90,
    0.0,
    0.35,
    0.35
},
opacityCurveRG{
    static_cast<double>(FCT_MinMaxCPoints),
    0.0,
    0.50,
    0.35,
    0.35,
    1.00,
    0.50,
    0.35,
    0.35
},
opacityCurveBY{
    static_cast<double>(FCT_MinMaxCPoints),
    0.0,
    0.50,
    0.35,
    0.35,
    1.00,
    0.50,
    0.35,
    0.35
},
opacityCurveW{
    static_cast<double>(FCT_MinMaxCPoints),
    0.00,
    0.35,
    0.35,
    0.00,
    0.35,
    0.75,
    0.35,
    0.35,
    0.60,
    0.75,
    0.35,
    0.35,
    1.00,
    0.35,
    0.00,
    0.00
},
opacityCurveWL{
    static_cast<double>(FCT_MinMaxCPoints),
    0.0,
    0.50,
    0.35,
    0.35,
    1.00,
    0.50,
    0.35,
    0.35
},
hhcurve{
    FCT_Linear
},
Chcurve{
    FCT_Linear
},
wavclCurve {
    DCT_Linear
},
enabled(false),
        median(false),
        medianlev(false),
        linkedg(true),
        cbenab(false),
        greenlow(0),
        bluelow(0),
        greenmed(0),
        bluemed(0),
        greenhigh(0),
        bluehigh(0),
        lipst(false),
        avoid(false),
        tmr(false),
        strength(100),
        balance(0),
        iter(0),
        expcontrast(false),
        expchroma(false),
        c{},
        ch{},
        expedge(false),
        expresid(false),
        expfinal(false),
        exptoning(false),
        expnoise(false),
        Lmethod(4),
        CLmethod("all"),
        Backmethod("grey"),
        Tilesmethod("full"),
        daubcoeffmethod("4_"),
        CHmethod("without"),
        Medgreinf("less"),
        CHSLmethod("SL"),
        EDmethod("CU"),
        NPmethod("none"),
        BAmethod("none"),
        TMmethod("cont"),
        Dirmethod("all"),
        HSmethod("with"),
        rescon(0),
        resconH(0),
        reschro(0),
        tmrs(0),
        gamma(1),
        sup(0),
        sky(0.0),
        thres(7),
        chroma(5),
        chro(0),
        threshold(5),
        threshold2(4),
        edgedetect(90),
        edgedetectthr(20),
        edgedetectthr2(0),
        edgesensi(60),
        edgeampli(10),
        contrast(0),
        edgrad(15),
        edgval(0),
        edgthresh(10),
        thr(35),
        thrH(65),
        skinprotect(0.0),
        hueskin(-5, 25, 170, 120, false),
        hueskin2(-260, -250, -130, -140, false),
        hllev(50, 75, 100, 98, false),
        bllev(0, 2, 50, 25, false),
        pastlev(0, 2, 30, 20, false),
        satlev(30, 45, 130, 100, false),
        edgcont(0, 10, 75, 40, false),
        level0noise(0, 0, false),
        level1noise(0, 0, false),
        level2noise(0, 0, false),
        level3noise(0, 0, false)
{
}

bool WaveletParams::operator ==(const WaveletParams& other) const
{
    return
        ccwcurve == other.ccwcurve
        && opacityCurveRG == other.opacityCurveRG
        && opacityCurveBY == other.opacityCurveBY
        && opacityCurveW == other.opacityCurveW
        && opacityCurveWL == other.opacityCurveWL
        && hhcurve == other.hhcurve
        && Chcurve == other.Chcurve
        && wavclCurve == other.wavclCurve
        && enabled == other.enabled
        && median == other.median
        && medianlev == other.medianlev
        && linkedg == other.linkedg
        && cbenab == other.cbenab
        && greenlow == other.greenlow
        && bluelow == other.bluelow
        && greenmed == other.greenmed
        && bluemed == other.bluemed
        && greenhigh == other.greenhigh
        && bluehigh == other.bluehigh
        && lipst == other.lipst
        && avoid == other.avoid
        && tmr == other.tmr
        && strength == other.strength
        && balance == other.balance
        && iter == other.iter
        && expcontrast == other.expcontrast
        && expchroma == other.expchroma
    && [this, &other]() -> bool {
        for (unsigned int i = 0; i < 9; ++i)
        {
            if (c[i] != other.c[i] || ch[i] != other.ch[i]) {
                return false;
            }
        }

        return true;
    }()
    && expedge == other.expedge
    && expresid == other.expresid
    && expfinal == other.expfinal
    && exptoning == other.exptoning
    && expnoise == other.expnoise
    && Lmethod == other.Lmethod
    && CLmethod == other.CLmethod
    && Backmethod == other.Backmethod
    && Tilesmethod == other.Tilesmethod
    && daubcoeffmethod == other.daubcoeffmethod
    && CHmethod == other.CHmethod
    && Medgreinf == other.Medgreinf
    && CHSLmethod == other.CHSLmethod
    && EDmethod == other.EDmethod
    && NPmethod == other.NPmethod
    && BAmethod == other.BAmethod
    && TMmethod == other.TMmethod
    && Dirmethod == other.Dirmethod
    && HSmethod == other.HSmethod
    && rescon == other.rescon
    && resconH == other.resconH
    && reschro == other.reschro
    && tmrs == other.tmrs
    && gamma == other.gamma
    && sup == other.sup
    && sky == other.sky
    && thres == other.thres
    && chroma == other.chroma
    && chro == other.chro
    && threshold == other.threshold
    && threshold2 == other.threshold2
    && edgedetect == other.edgedetect
    && edgedetectthr == other.edgedetectthr
    && edgedetectthr2 == other.edgedetectthr2
    && edgesensi == other.edgesensi
    && edgeampli == other.edgeampli
    && contrast == other.contrast
    && edgrad == other.edgrad
    && edgval == other.edgval
    && edgthresh == other.edgthresh
    && thr == other.thr
    && thrH == other.thrH
    && skinprotect == other.skinprotect
    && hueskin == other.hueskin
    && hueskin2 == other.hueskin2
    && hllev == other.hllev
    && bllev == other.bllev
    && pastlev == other.pastlev
    && satlev == other.satlev
    && edgcont == other.edgcont
    && level0noise == other.level0noise
    && level1noise == other.level1noise
    && level2noise == other.level2noise
    && level3noise == other.level3noise;
}

bool WaveletParams::operator !=(const WaveletParams& other) const
{
    return !(*this == other);
}

void WaveletParams::getCurves(
    WavCurve& cCurve,
    WavOpacityCurveRG& opacityCurveLUTRG,
    WavOpacityCurveBY& opacityCurveLUTBY,
    WavOpacityCurveW& opacityCurveLUTW,
    WavOpacityCurveWL& opacityCurveLUTWL
) const
{
    cCurve.Set(this->ccwcurve);
    opacityCurveLUTRG.Set(this->opacityCurveRG);
    opacityCurveLUTBY.Set(this->opacityCurveBY);
    opacityCurveLUTW.Set(this->opacityCurveW);
    opacityCurveLUTWL.Set(this->opacityCurveWL);

}

DirPyrEqualizerParams::DirPyrEqualizerParams() :
    enabled(false),
    gamutlab(false),
    mult{
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0
},
threshold(0.2),
skinprotect(0.0),
hueskin(-5, 25, 170, 120, false),
cbdlMethod("bef")
{
}

bool DirPyrEqualizerParams::operator ==(const DirPyrEqualizerParams& other) const
{
    return
        enabled == other.enabled
        && gamutlab == other.gamutlab
    && [this, &other]() -> bool {
        for (unsigned int i = 0; i < 6; ++i)
        {
            if (mult[i] != other.mult[i]) {
                return false;
            }
        }

        return true;
    }()
    && threshold == other.threshold
    && skinprotect == other.skinprotect
    && hueskin == other.hueskin
    && cbdlMethod == other.cbdlMethod;
}

bool DirPyrEqualizerParams::operator !=(const DirPyrEqualizerParams& other) const
{
    return !(*this == other);
}

HSVEqualizerParams::HSVEqualizerParams() :
    enabled(false),
    hcurve{
    FCT_Linear
},
scurve{
    FCT_Linear
},
vcurve{
    FCT_Linear
}
{
}

bool HSVEqualizerParams::operator ==(const HSVEqualizerParams& other) const
{
    return
        enabled == other.enabled
        && hcurve == other.hcurve
        && scurve == other.scurve
        && vcurve == other.vcurve;
}

bool HSVEqualizerParams::operator !=(const HSVEqualizerParams& other) const
{
    return !(*this == other);
}

FilmSimulationParams::FilmSimulationParams() :
    enabled(false),
    strength(100)
{
}

bool FilmSimulationParams::operator ==(const FilmSimulationParams& other) const
{
    return
        enabled == other.enabled
        && clutFilename == other.clutFilename
        && strength == other.strength;
}

bool FilmSimulationParams::operator !=(const FilmSimulationParams& other) const
{
    return !(*this == other);
}

RAWParams::BayerSensor::BayerSensor() :
    method(getMethodString(Method::AMAZE)),
    imageNum(0),
    ccSteps(0),
    black0(0.0),
    black1(0.0),
    black2(0.0),
    black3(0.0),
    twogreen(true),
    linenoise(0),
    linenoiseDirection(LineNoiseDirection::BOTH),
    greenthresh(0),
    dcb_iterations(2),
    lmmse_iterations(2),
    pixelShiftMotionCorrectionMethod(PSMotionCorrectionMethod::AUTO),
    pixelShiftEperIso(0.0),
    pixelShiftSigma(1.0),
    pixelShiftShowMotion(false),
    pixelShiftShowMotionMaskOnly(false),
    pixelShiftHoleFill(true),
    pixelShiftMedian(false),
    pixelShiftGreen(true),
    pixelShiftBlur(true),
    pixelShiftSmoothFactor(0.7),
    pixelShiftLmmse(false),
    pixelShiftEqualBright(false),
    pixelShiftEqualBrightChannel(false),
    pixelShiftNonGreenCross(true),
    dcb_enhance(true),
    pdafLinesFilter(false)
{
}

bool RAWParams::BayerSensor::operator ==(const BayerSensor& other) const
{
    return
        method == other.method
        && imageNum == other.imageNum
        && ccSteps == other.ccSteps
        && black0 == other.black0
        && black1 == other.black1
        && black2 == other.black2
        && black3 == other.black3
        && twogreen == other.twogreen
        && linenoise == other.linenoise
        && linenoiseDirection == other.linenoiseDirection
        && greenthresh == other.greenthresh
        && dcb_iterations == other.dcb_iterations
        && lmmse_iterations == other.lmmse_iterations
        && pixelShiftMotionCorrectionMethod == other.pixelShiftMotionCorrectionMethod
        && pixelShiftEperIso == other.pixelShiftEperIso
        && pixelShiftSigma == other.pixelShiftSigma
        && pixelShiftShowMotion == other.pixelShiftShowMotion
        && pixelShiftShowMotionMaskOnly == other.pixelShiftShowMotionMaskOnly
        && pixelShiftHoleFill == other.pixelShiftHoleFill
        && pixelShiftMedian == other.pixelShiftMedian
        && pixelShiftGreen == other.pixelShiftGreen
        && pixelShiftBlur == other.pixelShiftBlur
        && pixelShiftSmoothFactor == other.pixelShiftSmoothFactor
        && pixelShiftLmmse == other.pixelShiftLmmse
        && pixelShiftEqualBright == other.pixelShiftEqualBright
        && pixelShiftEqualBrightChannel == other.pixelShiftEqualBrightChannel
        && pixelShiftNonGreenCross == other.pixelShiftNonGreenCross
        && dcb_enhance == other.dcb_enhance
        && pdafLinesFilter == other.pdafLinesFilter;
}

bool RAWParams::BayerSensor::operator !=(const BayerSensor& other) const
{
    return !(*this == other);
}

void RAWParams::BayerSensor::setPixelShiftDefaults()
{
    pixelShiftMotionCorrectionMethod = RAWParams::BayerSensor::PSMotionCorrectionMethod::AUTO;
    pixelShiftEperIso = 0.0;
    pixelShiftSigma = 1.0;
    pixelShiftHoleFill = true;
    pixelShiftMedian = false;
    pixelShiftGreen = true;
    pixelShiftBlur = true;
    pixelShiftSmoothFactor = 0.7;
    pixelShiftLmmse = false;
    pixelShiftEqualBright = false;
    pixelShiftEqualBrightChannel = false;
    pixelShiftNonGreenCross = true;
}

const std::vector<const char*>& RAWParams::BayerSensor::getMethodStrings()
{
    static const std::vector<const char*> method_strings {
        "amaze",
        "igv",
        "lmmse",
        "eahd",
        "hphd",
        "vng4",
        "dcb",
        "ahd",
        "rcd",
        "fast",
        "mono",
        "none",
        "pixelshift"
    };
    return method_strings;
}

Glib::ustring RAWParams::BayerSensor::getMethodString(Method method)
{
    return getMethodStrings()[toUnderlying(method)];
}

RAWParams::XTransSensor::XTransSensor() :
    method(getMethodString(Method::THREE_PASS)),
    ccSteps(0),
    blackred(0.0),
    blackgreen(0.0),
    blackblue(0.0)
{
}

bool RAWParams::XTransSensor::operator ==(const XTransSensor& other) const
{
    return
        method == other.method
        && ccSteps == other.ccSteps
        && blackred == other.blackred
        && blackgreen == other.blackgreen
        && blackblue == other.blackblue;
}

bool RAWParams::XTransSensor::operator !=(const XTransSensor& other) const
{
    return !(*this == other);
}

const std::vector<const char*>& RAWParams::XTransSensor::getMethodStrings()
{
    static const std::vector<const char*> method_strings {
        "3-pass (best)",
        "1-pass (medium)",
        "fast",
        "mono",
        "none"
    };
    return method_strings;
}

Glib::ustring RAWParams::XTransSensor::getMethodString(Method method)
{
    return getMethodStrings()[toUnderlying(method)];
}

RAWParams::RAWParams() :
    df_autoselect(false),
    ff_AutoSelect(false),
    ff_BlurRadius(32),
    ff_BlurType(getFlatFieldBlurTypeString(FlatFieldBlurType::AREA)),
    ff_AutoClipControl(false),
    ff_clipControl(0),
    ca_autocorrect(false),
    cared(0.0),
    cablue(0.0),
    expos(1.0),
    preser(0.0),
    hotPixelFilter(false),
    deadPixelFilter(false),
    hotdeadpix_thresh(100)
{
}

bool RAWParams::operator ==(const RAWParams& other) const
{
    return
        bayersensor == other.bayersensor
        && xtranssensor == other.xtranssensor
        && dark_frame == other.dark_frame
        && df_autoselect == other.df_autoselect
        && ff_file == other.ff_file
        && ff_AutoSelect == other.ff_AutoSelect
        && ff_BlurRadius == other.ff_BlurRadius
        && ff_BlurType == other.ff_BlurType
        && ff_AutoClipControl == other.ff_AutoClipControl
        && ff_clipControl == other.ff_clipControl
        && ca_autocorrect == other.ca_autocorrect
        && cared == other.cared
        && cablue == other.cablue
        && expos == other.expos
        && preser == other.preser
        && hotPixelFilter == other.hotPixelFilter
        && deadPixelFilter == other.deadPixelFilter
        && hotdeadpix_thresh == other.hotdeadpix_thresh;
}

bool RAWParams::operator !=(const RAWParams& other) const
{
    return !(*this == other);
}

const std::vector<const char*>& RAWParams::getFlatFieldBlurTypeStrings()
{
    static const std::vector<const char*> blur_type_strings {
        "Area Flatfield",
        "Vertical Flatfield",
        "Horizontal Flatfield",
        "V+H Flatfield"
    };
    return blur_type_strings;
}

Glib::ustring RAWParams::getFlatFieldBlurTypeString(FlatFieldBlurType type)
{
    return getFlatFieldBlurTypeStrings()[toUnderlying(type)];
}


MetaDataParams::MetaDataParams():
    mode(MetaDataParams::TUNNEL)
{
}

bool MetaDataParams::operator==(const MetaDataParams &other) const
{
    return mode == other.mode;
}

bool MetaDataParams::operator!=(const MetaDataParams &other) const
{
    return !(*this == other);
}


ProcParams::ProcParams()
{
    setDefaults();
}

void ProcParams::setDefaults()
{
    toneCurve = ToneCurveParams();

    labCurve = LCurveParams();

    rgbCurves = RGBCurvesParams();

    localContrast = LocalContrastParams();

    colorToning = ColorToningParams();

    sharpenEdge = SharpenEdgeParams();

    sharpenMicro = SharpenMicroParams();

    sharpening = SharpeningParams();

    prsharpening = SharpeningParams();
    prsharpening.method = "rld";
    prsharpening.deconvamount = 100;
    prsharpening.deconvradius = 0.45;
    prsharpening.deconviter = 100;
    prsharpening.deconvdamping = 0;

    vibrance = VibranceParams();

    wb = WBParams();

    colorappearance = ColorAppearanceParams();

    defringe = DefringeParams();

    impulseDenoise = ImpulseDenoiseParams();

    dirpyrDenoise = DirPyrDenoiseParams();

    epd = EPDParams();

    fattal = FattalToneMappingParams();

    sh = SHParams();

    crop = CropParams();

    coarse = CoarseTransformParams();

    commonTrans = CommonTransformParams();

    rotate = RotateParams();

    distortion = DistortionParams();

    lensProf = LensProfParams();

    perspective = PerspectiveParams();

    gradient = GradientParams();

    pcvignette = PCVignetteParams();

    vignetting = VignettingParams();

    chmixer = ChannelMixerParams();

    blackwhite = BlackWhiteParams();

    cacorrection = CACorrParams();

    resize = ResizeParams();

    icm = ColorManagementParams();

    wavelet = WaveletParams();

    dirpyrequalizer = DirPyrEqualizerParams();

    hsvequalizer = HSVEqualizerParams();

    filmSimulation = FilmSimulationParams();

    raw = RAWParams();

    metadata = MetaDataParams();
    exif.clear();
    iptc.clear();

    rank = 0;
    colorlabel = 0;
    inTrash = false;

    ppVersion = PPVERSION;
}

int ProcParams::save(const Glib::ustring& fname, const Glib::ustring& fname2, bool fnameAbsolute, ParamsEdited* pedited)
{
    if (fname.empty() && fname2.empty()) {
        return 0;
    }

    Glib::ustring sPParams;

    try {
        Glib::KeyFile keyFile;

// Version
        keyFile.set_string("Version", "AppVersion", RTVERSION);
        keyFile.set_integer("Version", "Version", PPVERSION);

        saveToKeyfile(!pedited || pedited->general.rank, "General", "Rank", rank, keyFile);
        saveToKeyfile(!pedited || pedited->general.colorlabel, "General", "ColorLabel", colorlabel, keyFile);
        saveToKeyfile(!pedited || pedited->general.intrash, "General", "InTrash", inTrash, keyFile);

// Tone curve
        saveToKeyfile(!pedited || pedited->toneCurve.autoexp, "Exposure", "Auto", toneCurve.autoexp, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.clip, "Exposure", "Clip", toneCurve.clip, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.expcomp, "Exposure", "Compensation", toneCurve.expcomp, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.brightness, "Exposure", "Brightness", toneCurve.brightness, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.contrast, "Exposure", "Contrast", toneCurve.contrast, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.saturation, "Exposure", "Saturation", toneCurve.saturation, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.black, "Exposure", "Black", toneCurve.black, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.hlcompr, "Exposure", "HighlightCompr", toneCurve.hlcompr, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.hlcomprthresh, "Exposure", "HighlightComprThreshold", toneCurve.hlcomprthresh, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.shcompr, "Exposure", "ShadowCompr", toneCurve.shcompr, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.histmatching, "Exposure", "HistogramMatching", toneCurve.histmatching, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.clampOOG, "Exposure", "ClampOOG", toneCurve.clampOOG, keyFile);

// Highlight recovery
        saveToKeyfile(!pedited || pedited->toneCurve.hrenabled, "HLRecovery", "Enabled", toneCurve.hrenabled, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.method, "HLRecovery", "Method", toneCurve.method, keyFile);

        const std::map<ToneCurveParams::TcMode, const char*> tc_mapping = {
            {ToneCurveParams::TcMode::STD, "Standard"},
            {ToneCurveParams::TcMode::FILMLIKE, "FilmLike"},
            {ToneCurveParams::TcMode::SATANDVALBLENDING, "SatAndValueBlending"},
            {ToneCurveParams::TcMode::WEIGHTEDSTD, "WeightedStd"},
            {ToneCurveParams::TcMode::LUMINANCE, "Luminance"},
            {ToneCurveParams::TcMode::PERCEPTUAL, "Perceptual"}
        };

        saveToKeyfile(!pedited || pedited->toneCurve.curveMode, "Exposure", "CurveMode", tc_mapping, toneCurve.curveMode, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.curveMode2, "Exposure", "CurveMode2", tc_mapping, toneCurve.curveMode2, keyFile);

        saveToKeyfile(!pedited || pedited->toneCurve.curve, "Exposure", "Curve", toneCurve.curve, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.curve2, "Exposure", "Curve2", toneCurve.curve2, keyFile);

// Retinex
        saveToKeyfile(!pedited || pedited->retinex.enabled, "Retinex", "Enabled", retinex.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.str, "Retinex", "Str", retinex.str, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.scal, "Retinex", "Scal", retinex.scal, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.iter, "Retinex", "Iter", retinex.iter, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.grad, "Retinex", "Grad", retinex.grad, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.grads, "Retinex", "Grads", retinex.grads, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.gam, "Retinex", "Gam", retinex.gam, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.slope, "Retinex", "Slope", retinex.slope, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.medianmap, "Retinex", "Median", retinex.medianmap, keyFile);

        saveToKeyfile(!pedited || pedited->retinex.neigh, "Retinex", "Neigh", retinex.neigh, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.offs, "Retinex", "Offs", retinex.offs, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.vart, "Retinex", "Vart", retinex.vart, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.limd, "Retinex", "Limd", retinex.limd, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.highl, "Retinex", "highl", retinex.highl, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.skal, "Retinex", "skal", retinex.skal, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.retinexMethod, "Retinex", "RetinexMethod", retinex.retinexMethod, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.mapMethod, "Retinex", "mapMethod", retinex.mapMethod, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.viewMethod, "Retinex", "viewMethod", retinex.viewMethod, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.retinexcolorspace, "Retinex", "Retinexcolorspace", retinex.retinexcolorspace, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.gammaretinex, "Retinex", "Gammaretinex", retinex.gammaretinex, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.cdcurve, "Retinex", "CDCurve", retinex.cdcurve, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.mapcurve, "Retinex", "MAPCurve", retinex.mapcurve, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.cdHcurve, "Retinex", "CDHCurve", retinex.cdHcurve, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.lhcurve, "Retinex", "LHCurve", retinex.lhcurve, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.highlights, "Retinex", "Highlights", retinex.highlights, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.htonalwidth, "Retinex", "HighlightTonalWidth", retinex.htonalwidth, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.shadows, "Retinex", "Shadows", retinex.shadows, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.stonalwidth, "Retinex", "ShadowTonalWidth", retinex.stonalwidth, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.radius, "Retinex", "Radius", retinex.radius, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.transmissionCurve, "Retinex", "TransmissionCurve", retinex.transmissionCurve, keyFile);
        saveToKeyfile(!pedited || pedited->retinex.gaintransmissionCurve, "Retinex", "GainTransmissionCurve", retinex.gaintransmissionCurve, keyFile);

// Local contrast
        saveToKeyfile(!pedited || pedited->localContrast.enabled, "Local Contrast", "Enabled", localContrast.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->localContrast.radius, "Local Contrast", "Radius", localContrast.radius, keyFile);
        saveToKeyfile(!pedited || pedited->localContrast.amount, "Local Contrast", "Amount", localContrast.amount, keyFile);
        saveToKeyfile(!pedited || pedited->localContrast.darkness, "Local Contrast", "Darkness", localContrast.darkness, keyFile);
        saveToKeyfile(!pedited || pedited->localContrast.lightness, "Local Contrast", "Lightness", localContrast.lightness, keyFile);


// Channel mixer
        saveToKeyfile(!pedited || pedited->chmixer.enabled, "Channel Mixer", "Enabled", chmixer.enabled, keyFile);

        if (!pedited || pedited->chmixer.red[0] || pedited->chmixer.red[1] || pedited->chmixer.red[2]) {
            Glib::ArrayHandle<int> rmix(chmixer.red, 3, Glib::OWNERSHIP_NONE);
            keyFile.set_integer_list("Channel Mixer", "Red", rmix);
        }

        if (!pedited || pedited->chmixer.green[0] || pedited->chmixer.green[1] || pedited->chmixer.green[2]) {
            Glib::ArrayHandle<int> gmix(chmixer.green, 3, Glib::OWNERSHIP_NONE);
            keyFile.set_integer_list("Channel Mixer", "Green", gmix);
        }

        if (!pedited || pedited->chmixer.blue[0] || pedited->chmixer.blue[1] || pedited->chmixer.blue[2]) {
            Glib::ArrayHandle<int> bmix(chmixer.blue, 3, Glib::OWNERSHIP_NONE);
            keyFile.set_integer_list("Channel Mixer", "Blue", bmix);
        }

// Black & White
        saveToKeyfile(!pedited || pedited->blackwhite.enabled, "Black & White", "Enabled", blackwhite.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.method, "Black & White", "Method", blackwhite.method, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.autoc, "Black & White", "Auto", blackwhite.autoc, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.enabledcc, "Black & White", "ComplementaryColors", blackwhite.enabledcc, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.setting, "Black & White", "Setting", blackwhite.setting, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.filter, "Black & White", "Filter", blackwhite.filter, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerRed, "Black & White", "MixerRed", blackwhite.mixerRed, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerOrange, "Black & White", "MixerOrange", blackwhite.mixerOrange, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerYellow, "Black & White", "MixerYellow", blackwhite.mixerYellow, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerGreen, "Black & White", "MixerGreen", blackwhite.mixerGreen, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerCyan, "Black & White", "MixerCyan", blackwhite.mixerCyan, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerBlue, "Black & White", "MixerBlue", blackwhite.mixerBlue, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerMagenta, "Black & White", "MixerMagenta", blackwhite.mixerMagenta, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.mixerPurple, "Black & White", "MixerPurple", blackwhite.mixerPurple, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.gammaRed, "Black & White", "GammaRed", blackwhite.gammaRed, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.gammaGreen, "Black & White", "GammaGreen", blackwhite.gammaGreen, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.gammaBlue, "Black & White", "GammaBlue", blackwhite.gammaBlue, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.algo, "Black & White", "Algorithm", blackwhite.algo, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.luminanceCurve, "Black & White", "LuminanceCurve", blackwhite.luminanceCurve, keyFile);
        saveToKeyfile(
            !pedited || pedited->blackwhite.beforeCurveMode,
            "Black & White",
        "BeforeCurveMode", {
            {BlackWhiteParams::TcMode::STD_BW, "Standard"},
            {BlackWhiteParams::TcMode::FILMLIKE_BW, "FilmLike"},
            {BlackWhiteParams::TcMode::SATANDVALBLENDING_BW, "SatAndValueBlending"},
            {BlackWhiteParams::TcMode::WEIGHTEDSTD_BW, "WeightedStd"}

        },
        blackwhite.beforeCurveMode,
        keyFile
        );
        saveToKeyfile(
            !pedited || pedited->blackwhite.afterCurveMode,
            "Black & White",
        "AfterCurveMode", {
            {BlackWhiteParams::TcMode::STD_BW, "Standard"},
            {BlackWhiteParams::TcMode::WEIGHTEDSTD_BW, "WeightedStd"}

        },
        blackwhite.afterCurveMode,
        keyFile
        );
        saveToKeyfile(!pedited || pedited->blackwhite.beforeCurve, "Black & White", "BeforeCurve", blackwhite.beforeCurve, keyFile);
        saveToKeyfile(!pedited || pedited->blackwhite.afterCurve, "Black & White", "AfterCurve", blackwhite.afterCurve, keyFile);

// Luma curve
        saveToKeyfile(!pedited || pedited->labCurve.enabled, "Luminance Curve", "Enabled", labCurve.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.brightness, "Luminance Curve", "Brightness", labCurve.brightness, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.contrast, "Luminance Curve", "Contrast", labCurve.contrast, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.chromaticity, "Luminance Curve", "Chromaticity", labCurve.chromaticity, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.avoidcolorshift, "Luminance Curve", "AvoidColorShift", labCurve.avoidcolorshift, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.rstprotection, "Luminance Curve", "RedAndSkinTonesProtection", labCurve.rstprotection, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.lcredsk, "Luminance Curve", "LCredsk", labCurve.lcredsk, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.lcurve, "Luminance Curve", "LCurve", labCurve.lcurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.acurve, "Luminance Curve", "aCurve", labCurve.acurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.bcurve, "Luminance Curve", "bCurve", labCurve.bcurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.cccurve, "Luminance Curve", "ccCurve", labCurve.cccurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.chcurve, "Luminance Curve", "chCurve", labCurve.chcurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.lhcurve, "Luminance Curve", "lhCurve", labCurve.lhcurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.hhcurve, "Luminance Curve", "hhCurve", labCurve.hhcurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.lccurve, "Luminance Curve", "LcCurve", labCurve.lccurve, keyFile);
        saveToKeyfile(!pedited || pedited->labCurve.clcurve, "Luminance Curve", "ClCurve", labCurve.clcurve, keyFile);

// Sharpening
        saveToKeyfile(!pedited || pedited->sharpening.enabled, "Sharpening", "Enabled", sharpening.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.method, "Sharpening", "Method", sharpening.method, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.radius, "Sharpening", "Radius", sharpening.radius, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.amount, "Sharpening", "Amount", sharpening.amount, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.threshold, "Sharpening", "Threshold", sharpening.threshold.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.edgesonly, "Sharpening", "OnlyEdges", sharpening.edgesonly, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.edges_radius, "Sharpening", "EdgedetectionRadius", sharpening.edges_radius, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.edges_tolerance, "Sharpening", "EdgeTolerance", sharpening.edges_tolerance, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.halocontrol, "Sharpening", "HalocontrolEnabled", sharpening.halocontrol, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.halocontrol_amount, "Sharpening", "HalocontrolAmount", sharpening.halocontrol_amount, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.deconvradius, "Sharpening", "DeconvRadius", sharpening.deconvradius, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.deconvamount, "Sharpening", "DeconvAmount", sharpening.deconvamount, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.deconvdamping, "Sharpening", "DeconvDamping", sharpening.deconvdamping, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.deconviter, "Sharpening", "DeconvIterations", sharpening.deconviter, keyFile);

// Vibrance
        saveToKeyfile(!pedited || pedited->vibrance.enabled, "Vibrance", "Enabled", vibrance.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->vibrance.pastels, "Vibrance", "Pastels", vibrance.pastels, keyFile);
        saveToKeyfile(!pedited || pedited->vibrance.saturated, "Vibrance", "Saturated", vibrance.saturated, keyFile);
        saveToKeyfile(!pedited || pedited->vibrance.psthreshold, "Vibrance", "PSThreshold", vibrance.psthreshold.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->vibrance.protectskins, "Vibrance", "ProtectSkins", vibrance.protectskins, keyFile);
        saveToKeyfile(!pedited || pedited->vibrance.avoidcolorshift, "Vibrance", "AvoidColorShift", vibrance.avoidcolorshift, keyFile);
        saveToKeyfile(!pedited || pedited->vibrance.pastsattog, "Vibrance", "PastSatTog", vibrance.pastsattog, keyFile);
        saveToKeyfile(!pedited || pedited->vibrance.skintonescurve, "Vibrance", "SkinTonesCurve", vibrance.skintonescurve, keyFile);

// Edge sharpening
        saveToKeyfile(!pedited || pedited->sharpenEdge.enabled, "SharpenEdge", "Enabled", sharpenEdge.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->sharpenEdge.passes, "SharpenEdge", "Passes", sharpenEdge.passes, keyFile);
        saveToKeyfile(!pedited || pedited->sharpenEdge.amount, "SharpenEdge", "Strength", sharpenEdge.amount, keyFile);
        saveToKeyfile(!pedited || pedited->sharpenEdge.threechannels, "SharpenEdge", "ThreeChannels", sharpenEdge.threechannels, keyFile);

// Micro-contrast sharpening
        saveToKeyfile(!pedited || pedited->sharpenMicro.enabled, "SharpenMicro", "Enabled", sharpenMicro.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->sharpenMicro.matrix, "SharpenMicro", "Matrix", sharpenMicro.matrix, keyFile);
        saveToKeyfile(!pedited || pedited->sharpenMicro.amount, "SharpenMicro", "Strength", sharpenMicro.amount, keyFile);
        saveToKeyfile(!pedited || pedited->sharpenMicro.uniformity, "SharpenMicro", "Uniformity", sharpenMicro.uniformity, keyFile);

// WB
        saveToKeyfile(!pedited || pedited->wb.enabled, "White Balance", "Enabled", wb.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->wb.method, "White Balance", "Setting", wb.method, keyFile);
        saveToKeyfile(!pedited || pedited->wb.temperature, "White Balance", "Temperature", wb.temperature, keyFile);
        saveToKeyfile(!pedited || pedited->wb.green, "White Balance", "Green", wb.green, keyFile);
        saveToKeyfile(!pedited || pedited->wb.equal, "White Balance", "Equal", wb.equal, keyFile);
        saveToKeyfile(!pedited || pedited->wb.tempBias, "White Balance", "TemperatureBias", wb.tempBias, keyFile);

// Colorappearance
        saveToKeyfile(!pedited || pedited->colorappearance.enabled, "Color appearance", "Enabled", colorappearance.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.degree, "Color appearance", "Degree", colorappearance.degree, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.autodegree, "Color appearance", "AutoDegree", colorappearance.autodegree, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.degreeout, "Color appearance", "Degreeout",        colorappearance.degreeout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.autodegreeout, "Color appearance", "AutoDegreeout",    colorappearance.autodegreeout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.surround, "Color appearance", "Surround", colorappearance.surround, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.surrsrc, "Color appearance", "Surrsrc", colorappearance.surrsrc, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.adaplum, "Color appearance", "AdaptLum", colorappearance.adaplum, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.badpixsl, "Color appearance", "Badpixsl", colorappearance.badpixsl, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.wbmodel, "Color appearance", "Model", colorappearance.wbmodel, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.algo, "Color appearance", "Algorithm", colorappearance.algo, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.jlight, "Color appearance", "J-Light", colorappearance.jlight, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.qbright, "Color appearance", "Q-Bright", colorappearance.qbright, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.chroma, "Color appearance", "C-Chroma", colorappearance.chroma, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.schroma, "Color appearance", "S-Chroma", colorappearance.schroma, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.mchroma, "Color appearance", "M-Chroma", colorappearance.mchroma, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.contrast, "Color appearance", "J-Contrast", colorappearance.contrast, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.qcontrast, "Color appearance", "Q-Contrast", colorappearance.qcontrast, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.colorh, "Color appearance", "H-Hue", colorappearance.colorh, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.rstprotection, "Color appearance", "RSTProtection", colorappearance.rstprotection, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.adapscen, "Color appearance", "AdaptScene", colorappearance.adapscen, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.autoadapscen, "Color appearance", "AutoAdapscen", colorappearance.autoadapscen, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.ybscen, "Color appearance", "YbScene", colorappearance.ybscen, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.autoybscen, "Color appearance", "Autoybscen", colorappearance.autoybscen, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.surrsource, "Color appearance", "SurrSource", colorappearance.surrsource, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.gamut, "Color appearance", "Gamut", colorappearance.gamut, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.tempout, "Color appearance", "Tempout", colorappearance.tempout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.greenout, "Color appearance", "Greenout", colorappearance.greenout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.tempsc, "Color appearance", "Tempsc", colorappearance.tempsc, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.greensc, "Color appearance", "Greensc", colorappearance.greensc, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.ybout, "Color appearance", "Ybout", colorappearance.ybout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.datacie, "Color appearance", "Datacie", colorappearance.datacie, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.tonecie, "Color appearance", "Tonecie", colorappearance.tonecie, keyFile);

        const std::map<ColorAppearanceParams::TcMode, const char*> ca_mapping = {
            {ColorAppearanceParams::TcMode::LIGHT, "Lightness"},
            {ColorAppearanceParams::TcMode::BRIGHT, "Brightness"}
        };

        saveToKeyfile(!pedited || pedited->colorappearance.curveMode, "Color appearance", "CurveMode", ca_mapping, colorappearance.curveMode, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.curveMode2, "Color appearance", "CurveMode2", ca_mapping, colorappearance.curveMode2, keyFile);
        saveToKeyfile(
            !pedited || pedited->colorappearance.curveMode3,
            "Color appearance",
        "CurveMode3", {
            {ColorAppearanceParams::CtcMode::CHROMA, "Chroma"},
            {ColorAppearanceParams::CtcMode::SATUR, "Saturation"},
            {ColorAppearanceParams::CtcMode::COLORF, "Colorfullness"}

        },
        colorappearance.curveMode3,
        keyFile
        );
        saveToKeyfile(!pedited || pedited->colorappearance.curve, "Color appearance", "Curve", colorappearance.curve, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.curve2, "Color appearance", "Curve2", colorappearance.curve2, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.curve3, "Color appearance", "Curve3", colorappearance.curve3, keyFile);

// Impulse denoise
        saveToKeyfile(!pedited || pedited->impulseDenoise.enabled, "Impulse Denoising", "Enabled", impulseDenoise.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->impulseDenoise.thresh, "Impulse Denoising", "Threshold", impulseDenoise.thresh, keyFile);

// Defringe
        saveToKeyfile(!pedited || pedited->defringe.enabled, "Defringing", "Enabled", defringe.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->defringe.radius, "Defringing", "Radius", defringe.radius, keyFile);
        saveToKeyfile(!pedited || pedited->defringe.threshold, "Defringing", "Threshold", defringe.threshold, keyFile);
        saveToKeyfile(!pedited || pedited->defringe.huecurve, "Defringing", "HueCurve", defringe.huecurve, keyFile);

// Directional pyramid denoising
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.enabled, "Directional Pyramid Denoising", "Enabled", dirpyrDenoise.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.enhance, "Directional Pyramid Denoising", "Enhance", dirpyrDenoise.enhance, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.median, "Directional Pyramid Denoising", "Median", dirpyrDenoise.median, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.luma, "Directional Pyramid Denoising", "Luma", dirpyrDenoise.luma, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.Ldetail, "Directional Pyramid Denoising", "Ldetail", dirpyrDenoise.Ldetail, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.chroma, "Directional Pyramid Denoising", "Chroma", dirpyrDenoise.chroma, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.dmethod, "Directional Pyramid Denoising", "Method", dirpyrDenoise.dmethod, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.Lmethod, "Directional Pyramid Denoising", "LMethod", dirpyrDenoise.Lmethod, keyFile);

        if (dirpyrDenoise.Cmethod == "PRE") {
            dirpyrDenoise.Cmethod = "MAN"; // Never save 'auto chroma preview mode' to pp3
        }

        saveToKeyfile(!pedited || pedited->dirpyrDenoise.Cmethod, "Directional Pyramid Denoising", "CMethod", dirpyrDenoise.Cmethod, keyFile);

        if (dirpyrDenoise.C2method == "PREV") {
            dirpyrDenoise.C2method = "MANU";
        }

        saveToKeyfile(!pedited || pedited->dirpyrDenoise.C2method, "Directional Pyramid Denoising", "C2Method", dirpyrDenoise.C2method, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.smethod, "Directional Pyramid Denoising", "SMethod", dirpyrDenoise.smethod, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.medmethod, "Directional Pyramid Denoising", "MedMethod", dirpyrDenoise.medmethod, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.rgbmethod, "Directional Pyramid Denoising", "RGBMethod", dirpyrDenoise.rgbmethod, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.methodmed, "Directional Pyramid Denoising", "MethodMed", dirpyrDenoise.methodmed, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.redchro, "Directional Pyramid Denoising", "Redchro", dirpyrDenoise.redchro, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.bluechro, "Directional Pyramid Denoising", "Bluechro", dirpyrDenoise.bluechro, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.gamma, "Directional Pyramid Denoising", "Gamma", dirpyrDenoise.gamma, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.passes, "Directional Pyramid Denoising", "Passes", dirpyrDenoise.passes, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.lcurve, "Directional Pyramid Denoising", "LCurve", dirpyrDenoise.lcurve, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrDenoise.cccurve, "Directional Pyramid Denoising", "CCCurve", dirpyrDenoise.cccurve, keyFile);

// EPD
        saveToKeyfile(!pedited || pedited->epd.enabled, "EPD", "Enabled", epd.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->epd.strength, "EPD", "Strength", epd.strength, keyFile);
        saveToKeyfile(!pedited || pedited->epd.gamma, "EPD", "Gamma", epd.gamma, keyFile);
        saveToKeyfile(!pedited || pedited->epd.edgeStopping, "EPD", "EdgeStopping", epd.edgeStopping, keyFile);
        saveToKeyfile(!pedited || pedited->epd.scale, "EPD", "Scale", epd.scale, keyFile);
        saveToKeyfile(!pedited || pedited->epd.reweightingIterates, "EPD", "ReweightingIterates", epd.reweightingIterates, keyFile);

// Fattal
        saveToKeyfile(!pedited || pedited->fattal.enabled, "FattalToneMapping", "Enabled", fattal.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->fattal.threshold, "FattalToneMapping", "Threshold", fattal.threshold, keyFile);
        saveToKeyfile(!pedited || pedited->fattal.amount, "FattalToneMapping", "Amount", fattal.amount, keyFile);
        saveToKeyfile(!pedited || pedited->fattal.anchor, "FattalToneMapping", "Anchor", fattal.anchor, keyFile);

// Shadows & highlights
        saveToKeyfile(!pedited || pedited->sh.enabled, "Shadows & Highlights", "Enabled", sh.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->sh.highlights, "Shadows & Highlights", "Highlights", sh.highlights, keyFile);
        saveToKeyfile(!pedited || pedited->sh.htonalwidth, "Shadows & Highlights", "HighlightTonalWidth", sh.htonalwidth, keyFile);
        saveToKeyfile(!pedited || pedited->sh.shadows, "Shadows & Highlights", "Shadows", sh.shadows, keyFile);
        saveToKeyfile(!pedited || pedited->sh.stonalwidth, "Shadows & Highlights", "ShadowTonalWidth", sh.stonalwidth, keyFile);
        saveToKeyfile(!pedited || pedited->sh.radius, "Shadows & Highlights", "Radius", sh.radius, keyFile);

// Crop
        saveToKeyfile(!pedited || pedited->crop.enabled, "Crop", "Enabled", crop.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->crop.x, "Crop", "X", crop.x, keyFile);
        saveToKeyfile(!pedited || pedited->crop.y, "Crop", "Y", crop.y, keyFile);
        saveToKeyfile(!pedited || pedited->crop.w, "Crop", "W", crop.w, keyFile);
        saveToKeyfile(!pedited || pedited->crop.h, "Crop", "H", crop.h, keyFile);
        saveToKeyfile(!pedited || pedited->crop.fixratio, "Crop", "FixedRatio", crop.fixratio, keyFile);
        saveToKeyfile(!pedited || pedited->crop.ratio, "Crop", "Ratio", crop.ratio, keyFile);
        saveToKeyfile(!pedited || pedited->crop.orientation, "Crop", "Orientation", crop.orientation, keyFile);
        saveToKeyfile(!pedited || pedited->crop.guide, "Crop", "Guide", crop.guide, keyFile);

// Coarse transformation
        saveToKeyfile(!pedited || pedited->coarse.rotate, "Coarse Transformation", "Rotate", coarse.rotate, keyFile);
        saveToKeyfile(!pedited || pedited->coarse.hflip, "Coarse Transformation", "HorizontalFlip", coarse.hflip, keyFile);
        saveToKeyfile(!pedited || pedited->coarse.vflip, "Coarse Transformation", "VerticalFlip", coarse.vflip, keyFile);

// Common properties for transformations
        saveToKeyfile(!pedited || pedited->commonTrans.autofill, "Common Properties for Transformations", "AutoFill", commonTrans.autofill, keyFile);

// Rotation
        saveToKeyfile(!pedited || pedited->rotate.degree, "Rotation", "Degree", rotate.degree, keyFile);

// Distortion
        saveToKeyfile(!pedited || pedited->distortion.amount, "Distortion", "Amount", distortion.amount, keyFile);

// Lens profile
        saveToKeyfile(!pedited || pedited->lensProf.lcMode, "LensProfile", "LcMode", lensProf.getMethodString(lensProf.lcMode), keyFile);
        saveToKeyfile(!pedited || pedited->lensProf.lcpFile, "LensProfile", "LCPFile", relativePathIfInside(fname, fnameAbsolute, lensProf.lcpFile), keyFile);
        saveToKeyfile(!pedited || pedited->lensProf.useDist, "LensProfile", "UseDistortion", lensProf.useDist, keyFile);
        saveToKeyfile(!pedited || pedited->lensProf.useVign, "LensProfile", "UseVignette", lensProf.useVign, keyFile);
        saveToKeyfile(!pedited || pedited->lensProf.useCA, "LensProfile", "UseCA", lensProf.useCA, keyFile);
        saveToKeyfile(!pedited || pedited->lensProf.lfCameraMake, "LensProfile", "LFCameraMake", lensProf.lfCameraMake, keyFile);
        saveToKeyfile(!pedited || pedited->lensProf.lfCameraModel, "LensProfile", "LFCameraModel", lensProf.lfCameraModel, keyFile);
        saveToKeyfile(!pedited || pedited->lensProf.lfLens, "LensProfile", "LFLens", lensProf.lfLens, keyFile);

// Perspective correction
        saveToKeyfile(!pedited || pedited->perspective.horizontal, "Perspective", "Horizontal", perspective.horizontal, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.vertical, "Perspective", "Vertical", perspective.vertical, keyFile);

// Gradient
        saveToKeyfile(!pedited || pedited->gradient.enabled, "Gradient", "Enabled", gradient.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.degree, "Gradient", "Degree", gradient.degree, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.feather, "Gradient", "Feather", gradient.feather, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.strength, "Gradient", "Strength", gradient.strength, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.centerX, "Gradient", "CenterX", gradient.centerX, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.centerY, "Gradient", "CenterY", gradient.centerY, keyFile);

// Post-crop vignette
        saveToKeyfile(!pedited || pedited->pcvignette.enabled, "PCVignette", "Enabled", pcvignette.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->pcvignette.strength, "PCVignette", "Strength", pcvignette.strength, keyFile);
        saveToKeyfile(!pedited || pedited->pcvignette.feather, "PCVignette", "Feather", pcvignette.feather, keyFile);
        saveToKeyfile(!pedited || pedited->pcvignette.roundness, "PCVignette", "Roundness", pcvignette.roundness, keyFile);

// C/A correction
        saveToKeyfile(!pedited || pedited->cacorrection.red, "CACorrection", "Red", cacorrection.red, keyFile);
        saveToKeyfile(!pedited || pedited->cacorrection.blue, "CACorrection", "Blue", cacorrection.blue, keyFile);

// Vignetting correction
        saveToKeyfile(!pedited || pedited->vignetting.amount, "Vignetting Correction", "Amount", vignetting.amount, keyFile);
        saveToKeyfile(!pedited || pedited->vignetting.radius, "Vignetting Correction", "Radius", vignetting.radius, keyFile);
        saveToKeyfile(!pedited || pedited->vignetting.strength, "Vignetting Correction", "Strength", vignetting.strength, keyFile);
        saveToKeyfile(!pedited || pedited->vignetting.centerX, "Vignetting Correction", "CenterX", vignetting.centerX, keyFile);
        saveToKeyfile(!pedited || pedited->vignetting.centerY, "Vignetting Correction", "CenterY", vignetting.centerY, keyFile);

// Resize
        saveToKeyfile(!pedited || pedited->resize.enabled, "Resize", "Enabled", resize.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->resize.scale, "Resize", "Scale", resize.scale, keyFile);
        saveToKeyfile(!pedited || pedited->resize.appliesTo, "Resize", "AppliesTo", resize.appliesTo, keyFile);
        saveToKeyfile(!pedited || pedited->resize.method, "Resize", "Method", resize.method, keyFile);
        saveToKeyfile(!pedited || pedited->resize.dataspec, "Resize", "DataSpecified", resize.dataspec, keyFile);
        saveToKeyfile(!pedited || pedited->resize.width, "Resize", "Width", resize.width, keyFile);
        saveToKeyfile(!pedited || pedited->resize.height, "Resize", "Height", resize.height, keyFile);

// Post resize sharpening
        saveToKeyfile(!pedited || pedited->prsharpening.enabled, "PostResizeSharpening", "Enabled", prsharpening.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.method, "PostResizeSharpening", "Method", prsharpening.method, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.radius, "PostResizeSharpening", "Radius", prsharpening.radius, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.amount, "PostResizeSharpening", "Amount", prsharpening.amount, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.threshold, "PostResizeSharpening", "Threshold", prsharpening.threshold.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.edgesonly, "PostResizeSharpening", "OnlyEdges", prsharpening.edgesonly, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.edges_radius, "PostResizeSharpening", "EdgedetectionRadius", prsharpening.edges_radius, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.edges_tolerance, "PostResizeSharpening", "EdgeTolerance", prsharpening.edges_tolerance, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.halocontrol, "PostResizeSharpening", "HalocontrolEnabled", prsharpening.halocontrol, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.halocontrol_amount, "PostResizeSharpening", "HalocontrolAmount", prsharpening.halocontrol_amount, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.deconvradius, "PostResizeSharpening", "DeconvRadius", prsharpening.deconvradius, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.deconvamount, "PostResizeSharpening", "DeconvAmount", prsharpening.deconvamount, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.deconvdamping, "PostResizeSharpening", "DeconvDamping", prsharpening.deconvdamping, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.deconviter, "PostResizeSharpening", "DeconvIterations", prsharpening.deconviter, keyFile);

// Color management
        saveToKeyfile(!pedited || pedited->icm.input, "Color Management", "InputProfile", relativePathIfInside(fname, fnameAbsolute, icm.input), keyFile);
        saveToKeyfile(!pedited || pedited->icm.toneCurve, "Color Management", "ToneCurve", icm.toneCurve, keyFile);
        saveToKeyfile(!pedited || pedited->icm.applyLookTable, "Color Management", "ApplyLookTable", icm.applyLookTable, keyFile);
        saveToKeyfile(!pedited || pedited->icm.applyBaselineExposureOffset, "Color Management", "ApplyBaselineExposureOffset", icm.applyBaselineExposureOffset, keyFile);
        saveToKeyfile(!pedited || pedited->icm.applyHueSatMap, "Color Management", "ApplyHueSatMap", icm.applyHueSatMap, keyFile);
        saveToKeyfile(!pedited || pedited->icm.dcpIlluminant, "Color Management", "DCPIlluminant", icm.dcpIlluminant, keyFile);
        saveToKeyfile(!pedited || pedited->icm.working, "Color Management", "WorkingProfile", icm.working, keyFile);
        saveToKeyfile(!pedited || pedited->icm.output, "Color Management", "OutputProfile", icm.output, keyFile);
        saveToKeyfile(
            !pedited || pedited->icm.outputIntent,
            "Color Management",
        "OutputProfileIntent", {
            {RI_PERCEPTUAL, "Perceptual"},
            {RI_RELATIVE, "Relative"},
            {RI_SATURATION, "Saturation"},
            {RI_ABSOLUTE, "Absolute"}

        },
        icm.outputIntent,
        keyFile
        );
        saveToKeyfile(!pedited || pedited->icm.outputBPC, "Color Management", "OutputBPC", icm.outputBPC, keyFile);
        saveToKeyfile(!pedited || pedited->icm.gamma, "Color Management", "GammaCustom", icm.gamma, keyFile);
        saveToKeyfile(!pedited || pedited->icm.freegamma, "Color Management", "Freegamma", icm.freegamma, keyFile);
        saveToKeyfile(!pedited || pedited->icm.gampos, "Color Management", "GammaValue", icm.gampos, keyFile);
        saveToKeyfile(!pedited || pedited->icm.slpos, "Color Management", "GammaSlope", icm.slpos, keyFile);
        saveToKeyfile(!pedited || pedited->icm.predx, "Color Management", "GammaPredx", icm.predx, keyFile);
        saveToKeyfile(!pedited || pedited->icm.predy, "Color Management", "GammaPredy", icm.predy, keyFile);
        saveToKeyfile(!pedited || pedited->icm.pgrex, "Color Management", "GammaPgrex", icm.pgrex, keyFile);
        saveToKeyfile(!pedited || pedited->icm.pgrey, "Color Management", "GammaPgrey", icm.pgrey, keyFile);
        saveToKeyfile(!pedited || pedited->icm.pblux, "Color Management", "GammaPblux", icm.pblux, keyFile);
        saveToKeyfile(!pedited || pedited->icm.pbluy, "Color Management", "GammaPbluy", icm.pbluy, keyFile);
        saveToKeyfile(!pedited || pedited->icm.gamm, "Color Management", "GammaValueIn", icm.gamm, keyFile);
        saveToKeyfile(!pedited || pedited->icm.slop, "Color Management", "GammaSlopeIn", icm.slop, keyFile);

        saveToKeyfile(!pedited || pedited->icm.wprimaries, "Color Management", "GammaPrimaries", icm.wprimaries, keyFile);
        saveToKeyfile(!pedited || pedited->icm.wtemp, "Color Management", "GammaTemp", icm.wtemp, keyFile);
        saveToKeyfile(!pedited || pedited->icm.wprofile, "Color Management", "GammaProfile", icm.wprofile, keyFile);
        saveToKeyfile(!pedited || pedited->icm.wtrcin, "Color Management", "GammaTRCIN", icm.wtrcin, keyFile);

// Wavelet
        saveToKeyfile(!pedited || pedited->wavelet.enabled, "Wavelet", "Enabled", wavelet.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.strength, "Wavelet", "Strength", wavelet.strength, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.balance, "Wavelet", "Balance", wavelet.balance, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.iter, "Wavelet", "Iter", wavelet.iter, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.thres, "Wavelet", "MaxLev", wavelet.thres, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.Tilesmethod, "Wavelet", "TilesMethod", wavelet.Tilesmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.daubcoeffmethod, "Wavelet", "DaubMethod", wavelet.daubcoeffmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.CLmethod, "Wavelet", "ChoiceLevMethod", wavelet.CLmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.Backmethod, "Wavelet", "BackMethod", wavelet.Backmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.Lmethod, "Wavelet", "LevMethod", wavelet.Lmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.Dirmethod, "Wavelet", "DirMethod", wavelet.Dirmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.greenhigh, "Wavelet", "CBgreenhigh", wavelet.greenhigh, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.greenmed, "Wavelet", "CBgreenmed", wavelet.greenmed, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.greenlow, "Wavelet", "CBgreenlow", wavelet.greenlow, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.bluehigh, "Wavelet", "CBbluehigh", wavelet.bluehigh, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.bluemed, "Wavelet", "CBbluemed", wavelet.bluemed, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.bluelow, "Wavelet", "CBbluelow", wavelet.bluelow, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expcontrast, "Wavelet", "Expcontrast", wavelet.expcontrast, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expchroma, "Wavelet", "Expchroma", wavelet.expchroma, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expedge, "Wavelet", "Expedge", wavelet.expedge, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expresid, "Wavelet", "Expresid", wavelet.expresid, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expfinal, "Wavelet", "Expfinal", wavelet.expfinal, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.exptoning, "Wavelet", "Exptoning", wavelet.exptoning, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expnoise, "Wavelet", "Expnoise", wavelet.expnoise, keyFile);

        for (int i = 0; i < 9; i++) {
            std::stringstream ss;
            ss << "Contrast" << (i + 1);

            saveToKeyfile(!pedited || pedited->wavelet.c[i], "Wavelet", ss.str(), wavelet.c[i], keyFile);
        }

        for (int i = 0; i < 9; i++) {
            std::stringstream ss;
            ss << "Chroma" << (i + 1);

            saveToKeyfile(!pedited || pedited->wavelet.ch[i], "Wavelet", ss.str(), wavelet.ch[i], keyFile);
        }

        saveToKeyfile(!pedited || pedited->wavelet.sup, "Wavelet", "ContExtra", wavelet.sup, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.HSmethod, "Wavelet", "HSMethod", wavelet.HSmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.hllev, "Wavelet", "HLRange", wavelet.hllev.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.bllev, "Wavelet", "SHRange", wavelet.bllev.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgcont, "Wavelet", "Edgcont", wavelet.edgcont.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.level0noise, "Wavelet", "Level0noise", wavelet.level0noise.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.level1noise, "Wavelet", "Level1noise", wavelet.level1noise.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.level2noise, "Wavelet", "Level2noise", wavelet.level2noise.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.level3noise, "Wavelet", "Level3noise", wavelet.level3noise.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.threshold, "Wavelet", "ThresholdHighlight", wavelet.threshold, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.threshold2, "Wavelet", "ThresholdShadow", wavelet.threshold2, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgedetect, "Wavelet", "Edgedetect", wavelet.edgedetect, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgedetectthr, "Wavelet", "Edgedetectthr", wavelet.edgedetectthr, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgedetectthr2, "Wavelet", "EdgedetectthrHi", wavelet.edgedetectthr2, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgesensi, "Wavelet", "Edgesensi", wavelet.edgesensi, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgeampli, "Wavelet", "Edgeampli", wavelet.edgeampli, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.chroma, "Wavelet", "ThresholdChroma", wavelet.chroma, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.CHmethod, "Wavelet", "CHromaMethod", wavelet.CHmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.Medgreinf, "Wavelet", "Medgreinf", wavelet.Medgreinf, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.CHSLmethod, "Wavelet", "CHSLromaMethod", wavelet.CHSLmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.EDmethod, "Wavelet", "EDMethod", wavelet.EDmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.NPmethod, "Wavelet", "NPMethod", wavelet.NPmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.BAmethod, "Wavelet", "BAMethod", wavelet.BAmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.TMmethod, "Wavelet", "TMMethod", wavelet.TMmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.chro, "Wavelet", "ChromaLink", wavelet.chro, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.ccwcurve, "Wavelet", "ContrastCurve", wavelet.ccwcurve, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.pastlev, "Wavelet", "Pastlev", wavelet.pastlev.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.satlev, "Wavelet", "Satlev", wavelet.satlev.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.opacityCurveRG, "Wavelet", "OpacityCurveRG", wavelet.opacityCurveRG, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.opacityCurveBY, "Wavelet", "OpacityCurveBY", wavelet.opacityCurveBY, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.opacityCurveW, "Wavelet", "OpacityCurveW", wavelet.opacityCurveW, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.opacityCurveWL, "Wavelet", "OpacityCurveWL", wavelet.opacityCurveWL, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.hhcurve, "Wavelet", "HHcurve", wavelet.hhcurve, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.Chcurve, "Wavelet", "CHcurve", wavelet.Chcurve, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.wavclCurve, "Wavelet", "WavclCurve", wavelet.wavclCurve, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.median, "Wavelet", "Median", wavelet.median, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.medianlev, "Wavelet", "Medianlev", wavelet.medianlev, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.linkedg, "Wavelet", "Linkedg", wavelet.linkedg, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.cbenab, "Wavelet", "CBenab", wavelet.cbenab, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.lipst, "Wavelet", "Lipst", wavelet.lipst, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.skinprotect, "Wavelet", "Skinprotect", wavelet.skinprotect, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.hueskin, "Wavelet", "Hueskin", wavelet.hueskin.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgrad, "Wavelet", "Edgrad", wavelet.edgrad, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgval, "Wavelet", "Edgval", wavelet.edgval, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgthresh, "Wavelet", "ThrEdg", wavelet.edgthresh, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.avoid, "Wavelet", "AvoidColorShift", wavelet.avoid, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.tmr, "Wavelet", "TMr", wavelet.tmr, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.rescon, "Wavelet", "ResidualcontShadow", wavelet.rescon, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.resconH, "Wavelet", "ResidualcontHighlight", wavelet.resconH, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.thr, "Wavelet", "ThresholdResidShadow", wavelet.thr, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.thrH, "Wavelet", "ThresholdResidHighLight", wavelet.thrH, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.reschro, "Wavelet", "Residualchroma", wavelet.reschro, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.tmrs, "Wavelet", "ResidualTM", wavelet.tmrs, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.gamma, "Wavelet", "Residualgamma", wavelet.gamma, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.sky, "Wavelet", "HueRangeResidual", wavelet.sky, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.hueskin2, "Wavelet", "HueRange", wavelet.hueskin2.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.contrast, "Wavelet", "Contrast", wavelet.contrast, keyFile);

// Directional pyramid equalizer
        saveToKeyfile(!pedited || pedited->dirpyrequalizer.enabled, "Directional Pyramid Equalizer", "Enabled", dirpyrequalizer.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrequalizer.gamutlab, "Directional Pyramid Equalizer", "Gamutlab", dirpyrequalizer.gamutlab, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrequalizer.cbdlMethod, "Directional Pyramid Equalizer", "cbdlMethod", dirpyrequalizer.cbdlMethod, keyFile);

        for (int i = 0; i < 6; i++) {
            std::stringstream ss;
            ss << "Mult" << i;

            saveToKeyfile(!pedited || pedited->dirpyrequalizer.mult[i], "Directional Pyramid Equalizer", ss.str(), dirpyrequalizer.mult[i], keyFile);
        }

        saveToKeyfile(!pedited || pedited->dirpyrequalizer.threshold, "Directional Pyramid Equalizer", "Threshold", dirpyrequalizer.threshold, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrequalizer.skinprotect, "Directional Pyramid Equalizer", "Skinprotect", dirpyrequalizer.skinprotect, keyFile);
        saveToKeyfile(!pedited || pedited->dirpyrequalizer.hueskin, "Directional Pyramid Equalizer", "Hueskin", dirpyrequalizer.hueskin.toVector(), keyFile);

// HSV Equalizer
        saveToKeyfile(!pedited || pedited->hsvequalizer.enabled, "HSV Equalizer", "Enabled", hsvequalizer.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->hsvequalizer.hcurve, "HSV Equalizer", "HCurve", hsvequalizer.hcurve, keyFile);
        saveToKeyfile(!pedited || pedited->hsvequalizer.scurve, "HSV Equalizer", "SCurve", hsvequalizer.scurve, keyFile);
        saveToKeyfile(!pedited || pedited->hsvequalizer.vcurve, "HSV Equalizer", "VCurve", hsvequalizer.vcurve, keyFile);

// Film simulation
        saveToKeyfile(!pedited || pedited->filmSimulation.enabled, "Film Simulation", "Enabled", filmSimulation.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->filmSimulation.clutFilename, "Film Simulation", "ClutFilename", filmSimulation.clutFilename, keyFile);
        saveToKeyfile(!pedited || pedited->filmSimulation.strength, "Film Simulation", "Strength", filmSimulation.strength, keyFile);

        saveToKeyfile(!pedited || pedited->rgbCurves.enabled, "RGB Curves", "Enabled", rgbCurves.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->rgbCurves.lumamode, "RGB Curves", "LumaMode", rgbCurves.lumamode, keyFile);
        saveToKeyfile(!pedited || pedited->rgbCurves.rcurve, "RGB Curves", "rCurve", rgbCurves.rcurve, keyFile);
        saveToKeyfile(!pedited || pedited->rgbCurves.gcurve, "RGB Curves", "gCurve", rgbCurves.gcurve, keyFile);
        saveToKeyfile(!pedited || pedited->rgbCurves.bcurve, "RGB Curves", "bCurve", rgbCurves.bcurve, keyFile);

// Color toning
        saveToKeyfile(!pedited || pedited->colorToning.enabled, "ColorToning", "Enabled", colorToning.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.method, "ColorToning", "Method", colorToning.method, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.lumamode, "ColorToning", "Lumamode", colorToning.lumamode, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.twocolor, "ColorToning", "Twocolor", colorToning.twocolor, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.redlow, "ColorToning", "Redlow", colorToning.redlow, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.greenlow, "ColorToning", "Greenlow", colorToning.greenlow, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.bluelow, "ColorToning", "Bluelow", colorToning.bluelow, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.satlow, "ColorToning", "Satlow", colorToning.satlow, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.balance, "ColorToning", "Balance", colorToning.balance, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.sathigh, "ColorToning", "Sathigh", colorToning.sathigh, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.redmed, "ColorToning", "Redmed", colorToning.redmed, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.greenmed, "ColorToning", "Greenmed", colorToning.greenmed, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.bluemed, "ColorToning", "Bluemed", colorToning.bluemed, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.redhigh, "ColorToning", "Redhigh", colorToning.redhigh, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.greenhigh, "ColorToning", "Greenhigh", colorToning.greenhigh, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.bluehigh, "ColorToning", "Bluehigh", colorToning.bluehigh, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.autosat, "ColorToning", "Autosat", colorToning.autosat, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.opacityCurve, "ColorToning", "OpacityCurve", colorToning.opacityCurve, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.colorCurve, "ColorToning", "ColorCurve", colorToning.colorCurve, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.satprotectionthreshold, "ColorToning", "SatProtectionThreshold", colorToning.satProtectionThreshold, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.saturatedopacity, "ColorToning", "SaturatedOpacity", colorToning.saturatedOpacity, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.strength, "ColorToning", "Strength", colorToning.strength, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.hlColSat, "ColorToning", "HighlightsColorSaturation", colorToning.hlColSat.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.shadowsColSat, "ColorToning", "ShadowsColorSaturation", colorToning.shadowsColSat.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.clcurve, "ColorToning", "ClCurve", colorToning.clcurve, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.cl2curve, "ColorToning", "Cl2Curve", colorToning.cl2curve, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.labgridALow, "ColorToning", "LabGridALow", colorToning.labgridALow, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.labgridBLow, "ColorToning", "LabGridBLow", colorToning.labgridBLow, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.labgridAHigh, "ColorToning", "LabGridAHigh", colorToning.labgridAHigh, keyFile);
        saveToKeyfile(!pedited || pedited->colorToning.labgridBHigh, "ColorToning", "LabGridBHigh", colorToning.labgridBHigh, keyFile);

// Raw
        saveToKeyfile(!pedited || pedited->raw.darkFrame, "RAW", "DarkFrame", relativePathIfInside(fname, fnameAbsolute, raw.dark_frame), keyFile);
        saveToKeyfile(!pedited || pedited->raw.df_autoselect, "RAW", "DarkFrameAuto", raw.df_autoselect, keyFile);
        saveToKeyfile(!pedited || pedited->raw.ff_file, "RAW", "FlatFieldFile", relativePathIfInside(fname, fnameAbsolute, raw.ff_file), keyFile);
        saveToKeyfile(!pedited || pedited->raw.ff_AutoSelect, "RAW", "FlatFieldAutoSelect", raw.ff_AutoSelect, keyFile);
        saveToKeyfile(!pedited || pedited->raw.ff_BlurRadius, "RAW", "FlatFieldBlurRadius", raw.ff_BlurRadius, keyFile);
        saveToKeyfile(!pedited || pedited->raw.ff_BlurType, "RAW", "FlatFieldBlurType", raw.ff_BlurType, keyFile);
        saveToKeyfile(!pedited || pedited->raw.ff_AutoClipControl, "RAW", "FlatFieldAutoClipControl", raw.ff_AutoClipControl, keyFile);
        saveToKeyfile(!pedited || pedited->raw.ff_clipControl, "RAW", "FlatFieldClipControl", raw.ff_clipControl, keyFile);
        saveToKeyfile(!pedited || pedited->raw.ca_autocorrect, "RAW", "CA", raw.ca_autocorrect, keyFile);
        saveToKeyfile(!pedited || pedited->raw.cared, "RAW", "CARed", raw.cared, keyFile);
        saveToKeyfile(!pedited || pedited->raw.cablue, "RAW", "CABlue", raw.cablue, keyFile);
        saveToKeyfile(!pedited || pedited->raw.hotPixelFilter, "RAW", "HotPixelFilter", raw.hotPixelFilter, keyFile);
        saveToKeyfile(!pedited || pedited->raw.deadPixelFilter, "RAW", "DeadPixelFilter", raw.deadPixelFilter, keyFile);
        saveToKeyfile(!pedited || pedited->raw.hotdeadpix_thresh, "RAW", "HotDeadPixelThresh", raw.hotdeadpix_thresh, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.method, "RAW Bayer", "Method", raw.bayersensor.method, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.imageNum, "RAW Bayer", "ImageNum", raw.bayersensor.imageNum + 1, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.ccSteps, "RAW Bayer", "CcSteps", raw.bayersensor.ccSteps, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.exBlack0, "RAW Bayer", "PreBlack0", raw.bayersensor.black0, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.exBlack1, "RAW Bayer", "PreBlack1", raw.bayersensor.black1, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.exBlack2, "RAW Bayer", "PreBlack2", raw.bayersensor.black2, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.exBlack3, "RAW Bayer", "PreBlack3", raw.bayersensor.black3, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.exTwoGreen, "RAW Bayer", "PreTwoGreen", raw.bayersensor.twogreen, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.linenoise, "RAW Bayer", "LineDenoise", raw.bayersensor.linenoise, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.linenoise, "RAW Bayer", "LineDenoiseDirection", toUnderlying(raw.bayersensor.linenoiseDirection), keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.greenEq, "RAW Bayer", "GreenEqThreshold", raw.bayersensor.greenthresh, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.dcbIterations, "RAW Bayer", "DCBIterations", raw.bayersensor.dcb_iterations, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.dcbEnhance, "RAW Bayer", "DCBEnhance", raw.bayersensor.dcb_enhance, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.lmmseIterations, "RAW Bayer", "LMMSEIterations", raw.bayersensor.lmmse_iterations, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod, "RAW Bayer", "PixelShiftMotionCorrectionMethod", toUnderlying(raw.bayersensor.pixelShiftMotionCorrectionMethod), keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftEperIso, "RAW Bayer", "PixelShiftEperIso", raw.bayersensor.pixelShiftEperIso, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftSigma, "RAW Bayer", "PixelShiftSigma", raw.bayersensor.pixelShiftSigma, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftShowMotion, "RAW Bayer", "PixelShiftShowMotion", raw.bayersensor.pixelShiftShowMotion, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftShowMotionMaskOnly, "RAW Bayer", "PixelShiftShowMotionMaskOnly", raw.bayersensor.pixelShiftShowMotionMaskOnly, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftHoleFill, "RAW Bayer", "pixelShiftHoleFill", raw.bayersensor.pixelShiftHoleFill, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftMedian, "RAW Bayer", "pixelShiftMedian", raw.bayersensor.pixelShiftMedian, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftGreen, "RAW Bayer", "pixelShiftGreen", raw.bayersensor.pixelShiftGreen, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftBlur, "RAW Bayer", "pixelShiftBlur", raw.bayersensor.pixelShiftBlur, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftSmooth, "RAW Bayer", "pixelShiftSmoothFactor", raw.bayersensor.pixelShiftSmoothFactor, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftLmmse, "RAW Bayer", "pixelShiftLmmse", raw.bayersensor.pixelShiftLmmse, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftEqualBright, "RAW Bayer", "pixelShiftEqualBright", raw.bayersensor.pixelShiftEqualBright, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftEqualBrightChannel, "RAW Bayer", "pixelShiftEqualBrightChannel", raw.bayersensor.pixelShiftEqualBrightChannel, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftNonGreenCross, "RAW Bayer", "pixelShiftNonGreenCross", raw.bayersensor.pixelShiftNonGreenCross, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pdafLinesFilter, "RAW Bayer", "PDAFLinesFilter", raw.bayersensor.pdafLinesFilter, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.method, "RAW X-Trans", "Method", raw.xtranssensor.method, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.ccSteps, "RAW X-Trans", "CcSteps", raw.xtranssensor.ccSteps, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.exBlackRed, "RAW X-Trans", "PreBlackRed", raw.xtranssensor.blackred, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.exBlackGreen, "RAW X-Trans", "PreBlackGreen", raw.xtranssensor.blackgreen, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.exBlackBlue, "RAW X-Trans", "PreBlackBlue", raw.xtranssensor.blackblue, keyFile);

// Raw exposition
        saveToKeyfile(!pedited || pedited->raw.exPos, "RAW", "PreExposure", raw.expos, keyFile);
        saveToKeyfile(!pedited || pedited->raw.exPreser, "RAW", "PrePreserv", raw.preser, keyFile);

// MetaData
        saveToKeyfile(!pedited || pedited->metadata.mode, "MetaData", "Mode", metadata.mode, keyFile);

// EXIF change list
        if (!pedited || pedited->exif) {
            for (ExifPairs::const_iterator i = exif.begin(); i != exif.end(); ++i) {
                keyFile.set_string("Exif", i->first, i->second);
            }
        }

// IPTC change list
        if (!pedited || pedited->iptc) {
            for (IPTCPairs::const_iterator i = iptc.begin(); i != iptc.end(); ++i) {
                Glib::ArrayHandle<Glib::ustring> values = i->second;
                keyFile.set_string_list("IPTC", i->first, values);
            }
        }

        sPParams = keyFile.to_data();

    } catch (Glib::KeyFileError&) {}

    if (sPParams.empty()) {
        return 1;
    }

    int error1, error2;
    error1 = write(fname, sPParams);

    if (!fname2.empty()) {

        error2 = write(fname2, sPParams);
        // If at least one file has been saved, it's a success
        return error1 & error2;
    } else {
        return error1;
    }
}

int ProcParams::load(const Glib::ustring& fname, ParamsEdited* pedited)
{
    setlocale(LC_NUMERIC, "C");  // to set decimal point to "."

    if (fname.empty()) {
        return 1;
    }

    Glib::KeyFile keyFile;

    try {
        if (pedited) {
            pedited->set(false);
        }

        if (!Glib::file_test(fname, Glib::FILE_TEST_EXISTS) ||
                !keyFile.load_from_file(fname)) {
            return 1;
        }

        ppVersion = PPVERSION;
        appVersion = RTVERSION;

        if (keyFile.has_group("Version")) {
            if (keyFile.has_key("Version", "AppVersion")) {
                appVersion = keyFile.get_string("Version", "AppVersion");
            }

            if (keyFile.has_key("Version", "Version")) {
                ppVersion = keyFile.get_integer("Version", "Version");
            }
        }

        if (keyFile.has_group("General")) {
            assignFromKeyfile(keyFile, "General", "Rank", pedited, rank, pedited->general.rank);
            assignFromKeyfile(keyFile, "General", "ColorLabel", pedited, colorlabel, pedited->general.colorlabel);
            assignFromKeyfile(keyFile, "General", "InTrash", pedited, inTrash, pedited->general.intrash);
        }

        if (keyFile.has_group("Exposure")) {
            if (ppVersion < PPVERSION_AEXP) {
                toneCurve.autoexp = false; // prevent execution of autoexp when opening file created with earlier versions of autoexp algorithm
            } else {
                assignFromKeyfile(keyFile, "Exposure", "Auto", pedited, toneCurve.autoexp, pedited->toneCurve.autoexp);
            }

            assignFromKeyfile(keyFile, "Exposure", "Clip", pedited, toneCurve.clip, pedited->toneCurve.clip);
            assignFromKeyfile(keyFile, "Exposure", "Compensation", pedited, toneCurve.expcomp, pedited->toneCurve.expcomp);
            assignFromKeyfile(keyFile, "Exposure", "Brightness", pedited, toneCurve.brightness, pedited->toneCurve.brightness);
            assignFromKeyfile(keyFile, "Exposure", "Contrast", pedited, toneCurve.contrast, pedited->toneCurve.contrast);
            assignFromKeyfile(keyFile, "Exposure", "Saturation", pedited, toneCurve.saturation, pedited->toneCurve.saturation);
            assignFromKeyfile(keyFile, "Exposure", "Black", pedited, toneCurve.black, pedited->toneCurve.black);
            assignFromKeyfile(keyFile, "Exposure", "HighlightCompr", pedited, toneCurve.hlcompr, pedited->toneCurve.hlcompr);
            assignFromKeyfile(keyFile, "Exposure", "HighlightComprThreshold", pedited, toneCurve.hlcomprthresh, pedited->toneCurve.hlcomprthresh);
            assignFromKeyfile(keyFile, "Exposure", "ShadowCompr", pedited, toneCurve.shcompr, pedited->toneCurve.shcompr);

            if (toneCurve.shcompr > 100) {
                toneCurve.shcompr = 100; // older pp3 files can have values above 100.
            }

            const std::map<std::string, ToneCurveParams::TcMode> tc_mapping = {
                {"Standard", ToneCurveParams::TcMode::STD},
                {"FilmLike", ToneCurveParams::TcMode::FILMLIKE},
                {"SatAndValueBlending", ToneCurveParams::TcMode::SATANDVALBLENDING},
                {"WeightedStd", ToneCurveParams::TcMode::WEIGHTEDSTD},
                {"Luminance", ToneCurveParams::TcMode::LUMINANCE},
                {"Perceptual", ToneCurveParams::TcMode::PERCEPTUAL}
            };

            assignFromKeyfile(keyFile, "Exposure", "CurveMode", pedited, tc_mapping, toneCurve.curveMode, pedited->toneCurve.curveMode);
            assignFromKeyfile(keyFile, "Exposure", "CurveMode2", pedited, tc_mapping, toneCurve.curveMode2, pedited->toneCurve.curveMode2);

            if (ppVersion > 200) {
                assignFromKeyfile(keyFile, "Exposure", "Curve", pedited, toneCurve.curve, pedited->toneCurve.curve);
                assignFromKeyfile(keyFile, "Exposure", "Curve2", pedited, toneCurve.curve2, pedited->toneCurve.curve2);
            }

            assignFromKeyfile(keyFile, "Exposure", "HistogramMatching", pedited, toneCurve.histmatching, pedited->toneCurve.histmatching);
            assignFromKeyfile(keyFile, "Exposure", "ClampOOG", pedited, toneCurve.clampOOG, pedited->toneCurve.clampOOG);
        }

        if (keyFile.has_group("HLRecovery")) {
            assignFromKeyfile(keyFile, "HLRecovery", "Enabled", pedited, toneCurve.hrenabled, pedited->toneCurve.hrenabled);
            assignFromKeyfile(keyFile, "HLRecovery", "Method", pedited, toneCurve.method, pedited->toneCurve.method);
        }

        if (keyFile.has_group("Channel Mixer")) {
            if (ppVersion >= 329) {
                assignFromKeyfile(keyFile, "Channel Mixer", "Enabled", pedited, chmixer.enabled, pedited->chmixer.enabled);
            } else {
                chmixer.enabled = true;

                if (pedited) {
                    pedited->chmixer.enabled = true;
                }
            }

            if (keyFile.has_key("Channel Mixer", "Red") && keyFile.has_key("Channel Mixer", "Green") && keyFile.has_key("Channel Mixer", "Blue")) {
                const std::vector<int> rmix = keyFile.get_integer_list("Channel Mixer", "Red");
                const std::vector<int> gmix = keyFile.get_integer_list("Channel Mixer", "Green");
                const std::vector<int> bmix = keyFile.get_integer_list("Channel Mixer", "Blue");

                if (rmix.size() == 3 && gmix.size() == 3 && bmix.size() == 3) {
                    memcpy(chmixer.red,   rmix.data(), 3 * sizeof(int));
                    memcpy(chmixer.green, gmix.data(), 3 * sizeof(int));
                    memcpy(chmixer.blue,  bmix.data(), 3 * sizeof(int));
                }

                if (pedited) {
                    pedited->chmixer.red[0] =   pedited->chmixer.red[1] =   pedited->chmixer.red[2] = true;
                    pedited->chmixer.green[0] = pedited->chmixer.green[1] = pedited->chmixer.green[2] = true;
                    pedited->chmixer.blue[0] =  pedited->chmixer.blue[1] =  pedited->chmixer.blue[2] = true;
                }
            }
        }

        if (keyFile.has_group("Black & White")) {
            assignFromKeyfile(keyFile, "Black & White", "Enabled", pedited, blackwhite.enabled, pedited->blackwhite.enabled);
            assignFromKeyfile(keyFile, "Black & White", "Method", pedited, blackwhite.method, pedited->blackwhite.method);
            assignFromKeyfile(keyFile, "Black & White", "Auto", pedited, blackwhite.autoc, pedited->blackwhite.autoc);
            assignFromKeyfile(keyFile, "Black & White", "ComplementaryColors", pedited, blackwhite.enabledcc, pedited->blackwhite.enabledcc);
            assignFromKeyfile(keyFile, "Black & White", "MixerRed", pedited, blackwhite.mixerRed, pedited->blackwhite.mixerRed);
            assignFromKeyfile(keyFile, "Black & White", "MixerOrange", pedited, blackwhite.mixerOrange, pedited->blackwhite.mixerOrange);
            assignFromKeyfile(keyFile, "Black & White", "MixerYellow", pedited, blackwhite.mixerYellow, pedited->blackwhite.mixerYellow);
            assignFromKeyfile(keyFile, "Black & White", "MixerGreen", pedited, blackwhite.mixerGreen, pedited->blackwhite.mixerGreen);
            assignFromKeyfile(keyFile, "Black & White", "MixerCyan", pedited, blackwhite.mixerCyan, pedited->blackwhite.mixerCyan);
            assignFromKeyfile(keyFile, "Black & White", "MixerBlue", pedited, blackwhite.mixerBlue, pedited->blackwhite.mixerBlue);
            assignFromKeyfile(keyFile, "Black & White", "MixerMagenta", pedited, blackwhite.mixerMagenta, pedited->blackwhite.mixerMagenta);
            assignFromKeyfile(keyFile, "Black & White", "MixerPurple", pedited, blackwhite.mixerPurple, pedited->blackwhite.mixerPurple);
            assignFromKeyfile(keyFile, "Black & White", "GammaRed", pedited, blackwhite.gammaRed, pedited->blackwhite.gammaRed);
            assignFromKeyfile(keyFile, "Black & White", "GammaGreen", pedited, blackwhite.gammaGreen, pedited->blackwhite.gammaGreen);
            assignFromKeyfile(keyFile, "Black & White", "GammaBlue", pedited, blackwhite.gammaBlue, pedited->blackwhite.gammaBlue);
            assignFromKeyfile(keyFile, "Black & White", "Filter", pedited, blackwhite.filter, pedited->blackwhite.filter);
            assignFromKeyfile(keyFile, "Black & White", "Setting", pedited, blackwhite.setting, pedited->blackwhite.setting);
            assignFromKeyfile(keyFile, "Black & White", "LuminanceCurve", pedited, blackwhite.luminanceCurve, pedited->blackwhite.luminanceCurve);

            assignFromKeyfile(keyFile, "Black & White", "BeforeCurve", pedited, blackwhite.beforeCurve, pedited->blackwhite.beforeCurve);

            assignFromKeyfile(keyFile, "Black & White", "Algorithm", pedited, blackwhite.algo, pedited->blackwhite.algo);
            assignFromKeyfile(
                keyFile,
                "Black & White",
                "BeforeCurveMode",
            pedited, {
                {"Standard", BlackWhiteParams::TcMode::STD_BW},
                {"FilmLike", BlackWhiteParams::TcMode::FILMLIKE_BW},
                {"SatAndValueBlending", BlackWhiteParams::TcMode::SATANDVALBLENDING_BW},
                {"WeightedStd", BlackWhiteParams::TcMode::WEIGHTEDSTD_BW}
            },
            blackwhite.beforeCurveMode,
            pedited->blackwhite.beforeCurveMode
            );

            assignFromKeyfile(keyFile, "Black & White", "AfterCurve", pedited, blackwhite.afterCurve, pedited->blackwhite.afterCurve);
            assignFromKeyfile(
                keyFile,
                "Black & White",
                "AfterCurveMode",
            pedited, {
                {"Standard", BlackWhiteParams::TcMode::STD_BW},
                {"WeightedStd", BlackWhiteParams::TcMode::WEIGHTEDSTD_BW}
            },
            blackwhite.afterCurveMode,
            pedited->blackwhite.afterCurveMode
            );
        }

        if (keyFile.has_group("Retinex")) {
            assignFromKeyfile(keyFile, "Retinex", "Median", pedited, retinex.medianmap, pedited->retinex.medianmap);
            assignFromKeyfile(keyFile, "Retinex", "RetinexMethod", pedited, retinex.retinexMethod, pedited->retinex.retinexMethod);
            assignFromKeyfile(keyFile, "Retinex", "mapMethod", pedited, retinex.mapMethod, pedited->retinex.mapMethod);
            assignFromKeyfile(keyFile, "Retinex", "viewMethod", pedited, retinex.viewMethod, pedited->retinex.viewMethod);

            assignFromKeyfile(keyFile, "Retinex", "Retinexcolorspace", pedited, retinex.retinexcolorspace, pedited->retinex.retinexcolorspace);
            assignFromKeyfile(keyFile, "Retinex", "Gammaretinex", pedited, retinex.gammaretinex, pedited->retinex.gammaretinex);
            assignFromKeyfile(keyFile, "Retinex", "Enabled", pedited, retinex.enabled, pedited->retinex.enabled);
            assignFromKeyfile(keyFile, "Retinex", "Neigh", pedited, retinex.neigh, pedited->retinex.neigh);
            assignFromKeyfile(keyFile, "Retinex", "Str", pedited, retinex.str, pedited->retinex.str);
            assignFromKeyfile(keyFile, "Retinex", "Scal", pedited, retinex.scal, pedited->retinex.scal);
            assignFromKeyfile(keyFile, "Retinex", "Iter", pedited, retinex.iter, pedited->retinex.iter);
            assignFromKeyfile(keyFile, "Retinex", "Grad", pedited, retinex.grad, pedited->retinex.grad);
            assignFromKeyfile(keyFile, "Retinex", "Grads", pedited, retinex.grads, pedited->retinex.grads);
            assignFromKeyfile(keyFile, "Retinex", "Gam", pedited, retinex.gam, pedited->retinex.gam);
            assignFromKeyfile(keyFile, "Retinex", "Slope", pedited, retinex.slope, pedited->retinex.slope);
            assignFromKeyfile(keyFile, "Retinex", "Offs", pedited, retinex.offs, pedited->retinex.offs);
            assignFromKeyfile(keyFile, "Retinex", "Vart", pedited, retinex.vart, pedited->retinex.vart);
            assignFromKeyfile(keyFile, "Retinex", "Limd", pedited, retinex.limd, pedited->retinex.limd);
            assignFromKeyfile(keyFile, "Retinex", "highl", pedited, retinex.highl, pedited->retinex.highl);
            assignFromKeyfile(keyFile, "Retinex", "skal", pedited, retinex.skal, pedited->retinex.skal);
            assignFromKeyfile(keyFile, "Retinex", "CDCurve", pedited, retinex.cdcurve, pedited->retinex.cdcurve);

            assignFromKeyfile(keyFile, "Retinex", "MAPCurve", pedited, retinex.mapcurve, pedited->retinex.mapcurve);

            assignFromKeyfile(keyFile, "Retinex", "CDHCurve", pedited, retinex.cdHcurve, pedited->retinex.cdHcurve);

            assignFromKeyfile(keyFile, "Retinex", "LHCurve", pedited, retinex.lhcurve, pedited->retinex.lhcurve);

            assignFromKeyfile(keyFile, "Retinex", "Highlights", pedited, retinex.highlights, pedited->retinex.highlights);
            assignFromKeyfile(keyFile, "Retinex", "HighlightTonalWidth", pedited, retinex.htonalwidth, pedited->retinex.htonalwidth);
            assignFromKeyfile(keyFile, "Retinex", "Shadows", pedited, retinex.shadows, pedited->retinex.shadows);
            assignFromKeyfile(keyFile, "Retinex", "ShadowTonalWidth", pedited, retinex.stonalwidth, pedited->retinex.stonalwidth);

            assignFromKeyfile(keyFile, "Retinex", "Radius", pedited, retinex.radius, pedited->retinex.radius);

            assignFromKeyfile(keyFile, "Retinex", "TransmissionCurve", pedited, retinex.transmissionCurve, pedited->retinex.transmissionCurve);

            assignFromKeyfile(keyFile, "Retinex", "GainTransmissionCurve", pedited, retinex.gaintransmissionCurve, pedited->retinex.gaintransmissionCurve);
        }

        if (keyFile.has_group("Local Contrast")) {
            assignFromKeyfile(keyFile, "Local Contrast", "Enabled", pedited, localContrast.enabled, pedited->localContrast.enabled);
            assignFromKeyfile(keyFile, "Local Contrast", "Radius", pedited, localContrast.radius, pedited->localContrast.radius);
            assignFromKeyfile(keyFile, "Local Contrast", "Amount", pedited, localContrast.amount, pedited->localContrast.amount);
            assignFromKeyfile(keyFile, "Local Contrast", "Darkness", pedited, localContrast.darkness, pedited->localContrast.darkness);
            assignFromKeyfile(keyFile, "Local Contrast", "Lightness", pedited, localContrast.lightness, pedited->localContrast.lightness);
        }

        if (keyFile.has_group("Luminance Curve")) {
            if (ppVersion >= 329) {
                assignFromKeyfile(keyFile, "Luminance Curve", "Enabled", pedited, labCurve.enabled, pedited->labCurve.enabled);
            } else {
                labCurve.enabled = true;

                if (pedited) {
                    pedited->labCurve.enabled = true;
                }
            }

            assignFromKeyfile(keyFile, "Luminance Curve", "Brightness", pedited, labCurve.brightness, pedited->labCurve.brightness);
            assignFromKeyfile(keyFile, "Luminance Curve", "Contrast", pedited, labCurve.contrast, pedited->labCurve.contrast);

            if (ppVersion < 303) {
                // transform Saturation into Chromaticity
                // if Saturation == 0, should we set BWToning on?
                assignFromKeyfile(keyFile, "Luminance Curve", "Saturation", pedited, labCurve.chromaticity, pedited->labCurve.chromaticity);
                // transform AvoidColorClipping into AvoidColorShift
                assignFromKeyfile(keyFile, "Luminance Curve", "AvoidColorClipping", pedited, labCurve.avoidcolorshift, pedited->labCurve.avoidcolorshift);
            } else {
                if (keyFile.has_key("Luminance Curve", "Chromaticity")) {
                    labCurve.chromaticity = keyFile.get_integer("Luminance Curve", "Chromaticity");

                    if (ppVersion >= 303 && ppVersion < 314 && labCurve.chromaticity == -100) {
                        blackwhite.enabled = true;
                    }

                    if (pedited) {
                        pedited->labCurve.chromaticity = true;
                    }
                }

                assignFromKeyfile(keyFile, "Luminance Curve", "AvoidColorShift", pedited, labCurve.avoidcolorshift, pedited->labCurve.avoidcolorshift);
                assignFromKeyfile(keyFile, "Luminance Curve", "RedAndSkinTonesProtection", pedited, labCurve.rstprotection, pedited->labCurve.rstprotection);
            }

            assignFromKeyfile(keyFile, "Luminance Curve", "LCredsk", pedited, labCurve.lcredsk, pedited->labCurve.lcredsk);

            if (ppVersion < 314) {
                // Backward compatibility: If BWtoning is true, Chromaticity has to be set to -100, which will produce the same effect
                // and will enable the b&w toning mode ('a' & 'b' curves)
                if (keyFile.has_key("Luminance Curve", "BWtoning")) {
                    if (keyFile.get_boolean("Luminance Curve", "BWtoning")) {
                        labCurve.chromaticity = -100;

                        if (pedited) {
                            pedited->labCurve.chromaticity = true;
                        }
                    }
                }
            }

            assignFromKeyfile(keyFile, "Luminance Curve", "LCurve", pedited, labCurve.lcurve, pedited->labCurve.lcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "aCurve", pedited, labCurve.acurve, pedited->labCurve.acurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "bCurve", pedited, labCurve.bcurve, pedited->labCurve.bcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "ccCurve", pedited, labCurve.cccurve, pedited->labCurve.cccurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "chCurve", pedited, labCurve.chcurve, pedited->labCurve.chcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "lhCurve", pedited, labCurve.lhcurve, pedited->labCurve.lhcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "hhCurve", pedited, labCurve.hhcurve, pedited->labCurve.hhcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "LcCurve", pedited, labCurve.lccurve, pedited->labCurve.lccurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "ClCurve", pedited, labCurve.clcurve, pedited->labCurve.clcurve);
        }

        if (keyFile.has_group("Sharpening")) {
            assignFromKeyfile(keyFile, "Sharpening", "Enabled", pedited, sharpening.enabled, pedited->sharpening.enabled);
            assignFromKeyfile(keyFile, "Sharpening", "Radius", pedited, sharpening.radius, pedited->sharpening.radius);
            assignFromKeyfile(keyFile, "Sharpening", "Amount", pedited, sharpening.amount, pedited->sharpening.amount);

            if (keyFile.has_key("Sharpening", "Threshold")) {
                if (ppVersion < 302) {
                    int thresh = min(keyFile.get_integer("Sharpening", "Threshold"), 2000);
                    sharpening.threshold.setValues(thresh, thresh, 2000, 2000);  // TODO: 2000 is the maximum value and is taken of rtgui/sharpening.cc ; should be changed by the tool modularization
                } else {
                    const std::vector<int> thresh = keyFile.get_integer_list("Sharpening", "Threshold");

                    if (thresh.size() >= 4) {
                        sharpening.threshold.setValues(thresh[0], thresh[1], min(thresh[2], 2000), min(thresh[3], 2000));
                    }
                }

                if (pedited) {
                    pedited->sharpening.threshold = true;
                }
            }

            assignFromKeyfile(keyFile, "Sharpening", "OnlyEdges", pedited, sharpening.edgesonly, pedited->sharpening.edgesonly);
            assignFromKeyfile(keyFile, "Sharpening", "EdgedetectionRadius", pedited, sharpening.edges_radius, pedited->sharpening.edges_radius);
            assignFromKeyfile(keyFile, "Sharpening", "EdgeTolerance", pedited, sharpening.edges_tolerance, pedited->sharpening.edges_tolerance);
            assignFromKeyfile(keyFile, "Sharpening", "HalocontrolEnabled", pedited, sharpening.halocontrol, pedited->sharpening.halocontrol);
            assignFromKeyfile(keyFile, "Sharpening", "HalocontrolAmount", pedited, sharpening.halocontrol_amount, pedited->sharpening.halocontrol_amount);
            assignFromKeyfile(keyFile, "Sharpening", "Method", pedited, sharpening.method, pedited->sharpening.method);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvRadius", pedited, sharpening.deconvradius, pedited->sharpening.deconvradius);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvAmount", pedited, sharpening.deconvamount, pedited->sharpening.deconvamount);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvDamping", pedited, sharpening.deconvdamping, pedited->sharpening.deconvdamping);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvIterations", pedited, sharpening.deconviter, pedited->sharpening.deconviter);
        }

        if (keyFile.has_group("SharpenEdge")) {
            assignFromKeyfile(keyFile, "SharpenEdge", "Enabled", pedited, sharpenEdge.enabled, pedited->sharpenEdge.enabled);
            assignFromKeyfile(keyFile, "SharpenEdge", "Passes", pedited, sharpenEdge.passes, pedited->sharpenEdge.passes);
            assignFromKeyfile(keyFile, "SharpenEdge", "Strength", pedited, sharpenEdge.amount, pedited->sharpenEdge.amount);
            assignFromKeyfile(keyFile, "SharpenEdge", "ThreeChannels", pedited, sharpenEdge.threechannels, pedited->sharpenEdge.threechannels);
        }

        if (keyFile.has_group("SharpenMicro")) {
            assignFromKeyfile(keyFile, "SharpenMicro", "Enabled", pedited, sharpenMicro.enabled, pedited->sharpenMicro.enabled);
            assignFromKeyfile(keyFile, "SharpenMicro", "Matrix", pedited, sharpenMicro.matrix, pedited->sharpenMicro.matrix);
            assignFromKeyfile(keyFile, "SharpenMicro", "Strength", pedited, sharpenMicro.amount, pedited->sharpenMicro.amount);
            assignFromKeyfile(keyFile, "SharpenMicro", "Uniformity", pedited, sharpenMicro.uniformity, pedited->sharpenMicro.uniformity);
        }

        if (keyFile.has_group("Vibrance")) {
            assignFromKeyfile(keyFile, "Vibrance", "Enabled", pedited, vibrance.enabled, pedited->vibrance.enabled);
            assignFromKeyfile(keyFile, "Vibrance", "Pastels", pedited, vibrance.pastels, pedited->vibrance.pastels);
            assignFromKeyfile(keyFile, "Vibrance", "Saturated", pedited, vibrance.saturated, pedited->vibrance.saturated);

            if (keyFile.has_key("Vibrance", "PSThreshold")) {
                if (ppVersion < 302) {
                    int thresh = keyFile.get_integer("Vibrance", "PSThreshold");
                    vibrance.psthreshold.setValues(thresh, thresh);
                } else {
                    const std::vector<int> thresh = keyFile.get_integer_list("Vibrance", "PSThreshold");

                    if (thresh.size() >= 2) {
                        vibrance.psthreshold.setValues(thresh[0], thresh[1]);
                    }
                }

                if (pedited) {
                    pedited->vibrance.psthreshold = true;
                }
            }

            assignFromKeyfile(keyFile, "Vibrance", "ProtectSkins", pedited, vibrance.protectskins, pedited->vibrance.protectskins);
            assignFromKeyfile(keyFile, "Vibrance", "AvoidColorShift", pedited, vibrance.avoidcolorshift, pedited->vibrance.avoidcolorshift);
            assignFromKeyfile(keyFile, "Vibrance", "PastSatTog", pedited, vibrance.pastsattog, pedited->vibrance.pastsattog);
            assignFromKeyfile(keyFile, "Vibrance", "SkinTonesCurve", pedited, vibrance.skintonescurve, pedited->vibrance.skintonescurve);
        }

        if (keyFile.has_group("White Balance")) {
            assignFromKeyfile(keyFile, "White Balance", "Enabled", pedited, wb.enabled, pedited->wb.enabled);
            assignFromKeyfile(keyFile, "White Balance", "Setting", pedited, wb.method, pedited->wb.method);
            assignFromKeyfile(keyFile, "White Balance", "Temperature", pedited, wb.temperature, pedited->wb.temperature);
            assignFromKeyfile(keyFile, "White Balance", "Green", pedited, wb.green, pedited->wb.green);
            assignFromKeyfile(keyFile, "White Balance", "Equal", pedited, wb.equal, pedited->wb.equal);
            assignFromKeyfile(keyFile, "White Balance", "TemperatureBias", pedited, wb.tempBias, pedited->wb.tempBias);
        }

        if (keyFile.has_group("Defringing")) {
            assignFromKeyfile(keyFile, "Defringing", "Enabled", pedited, defringe.enabled, pedited->defringe.enabled);
            assignFromKeyfile(keyFile, "Defringing", "Radius", pedited, defringe.radius, pedited->defringe.radius);

            if (keyFile.has_key("Defringing", "Threshold")) {
                defringe.threshold = (float)keyFile.get_integer("Defringing", "Threshold");

                if (pedited) {
                    pedited->defringe.threshold = true;
                }
            }

            if (ppVersion < 310) {
                defringe.threshold = sqrt(defringe.threshold * 33.f / 5.f);
            }

            assignFromKeyfile(keyFile, "Defringing", "HueCurve", pedited, defringe.huecurve, pedited->defringe.huecurve);
        }

        if (keyFile.has_group("Color appearance")) {
            assignFromKeyfile(keyFile, "Color appearance", "Enabled", pedited, colorappearance.enabled, pedited->colorappearance.enabled);
            assignFromKeyfile(keyFile, "Color appearance", "Degree", pedited, colorappearance.degree, pedited->colorappearance.degree);
            assignFromKeyfile(keyFile, "Color appearance", "AutoDegree", pedited, colorappearance.autodegree, pedited->colorappearance.autodegree);
            assignFromKeyfile(keyFile, "Color appearance", "Degreeout", pedited, colorappearance.degreeout, pedited->colorappearance.degreeout);

            assignFromKeyfile(keyFile, "Color appearance", "AutoDegreeout", pedited, colorappearance.autodegreeout, pedited->colorappearance.autodegreeout);

            assignFromKeyfile(keyFile, "Color appearance", "Surround", pedited, colorappearance.surround, pedited->colorappearance.surround);
            assignFromKeyfile(keyFile, "Color appearance", "Surrsrc", pedited, colorappearance.surrsrc, pedited->colorappearance.surrsrc);
            assignFromKeyfile(keyFile, "Color appearance", "AdaptLum", pedited, colorappearance.adaplum, pedited->colorappearance.adaplum);
            assignFromKeyfile(keyFile, "Color appearance", "Badpixsl", pedited, colorappearance.badpixsl, pedited->colorappearance.badpixsl);
            assignFromKeyfile(keyFile, "Color appearance", "Model", pedited, colorappearance.wbmodel, pedited->colorappearance.wbmodel);
            assignFromKeyfile(keyFile, "Color appearance", "Algorithm", pedited, colorappearance.algo, pedited->colorappearance.algo);
            assignFromKeyfile(keyFile, "Color appearance", "J-Light", pedited, colorappearance.jlight, pedited->colorappearance.jlight);
            assignFromKeyfile(keyFile, "Color appearance", "Q-Bright", pedited, colorappearance.qbright, pedited->colorappearance.qbright);
            assignFromKeyfile(keyFile, "Color appearance", "C-Chroma", pedited, colorappearance.chroma, pedited->colorappearance.chroma);
            assignFromKeyfile(keyFile, "Color appearance", "S-Chroma", pedited, colorappearance.schroma, pedited->colorappearance.schroma);
            assignFromKeyfile(keyFile, "Color appearance", "M-Chroma", pedited, colorappearance.mchroma, pedited->colorappearance.mchroma);
            assignFromKeyfile(keyFile, "Color appearance", "RSTProtection", pedited, colorappearance.rstprotection, pedited->colorappearance.rstprotection);
            assignFromKeyfile(keyFile, "Color appearance", "J-Contrast", pedited, colorappearance.contrast, pedited->colorappearance.contrast);
            assignFromKeyfile(keyFile, "Color appearance", "Q-Contrast", pedited, colorappearance.qcontrast, pedited->colorappearance.qcontrast);
            assignFromKeyfile(keyFile, "Color appearance", "H-Hue", pedited, colorappearance.colorh, pedited->colorappearance.colorh);
            assignFromKeyfile(keyFile, "Color appearance", "AdaptScene", pedited, colorappearance.adapscen, pedited->colorappearance.adapscen);
            assignFromKeyfile(keyFile, "Color appearance", "AutoAdapscen", pedited, colorappearance.autoadapscen, pedited->colorappearance.autoadapscen);
            assignFromKeyfile(keyFile, "Color appearance", "YbScene", pedited, colorappearance.ybscen, pedited->colorappearance.ybscen);
            assignFromKeyfile(keyFile, "Color appearance", "Autoybscen", pedited, colorappearance.autoybscen, pedited->colorappearance.autoybscen);
            assignFromKeyfile(keyFile, "Color appearance", "SurrSource", pedited, colorappearance.surrsource, pedited->colorappearance.surrsource);
            assignFromKeyfile(keyFile, "Color appearance", "Gamut", pedited, colorappearance.gamut, pedited->colorappearance.gamut);
            assignFromKeyfile(keyFile, "Color appearance", "Tempout", pedited, colorappearance.tempout, pedited->colorappearance.tempout);
            assignFromKeyfile(keyFile, "Color appearance", "Greenout", pedited, colorappearance.greenout, pedited->colorappearance.greenout);
            assignFromKeyfile(keyFile, "Color appearance", "Tempsc", pedited, colorappearance.tempsc, pedited->colorappearance.tempsc);
            assignFromKeyfile(keyFile, "Color appearance", "Greensc", pedited, colorappearance.greensc, pedited->colorappearance.greensc);
            assignFromKeyfile(keyFile, "Color appearance", "Ybout", pedited, colorappearance.ybout, pedited->colorappearance.ybout);
            assignFromKeyfile(keyFile, "Color appearance", "Datacie", pedited, colorappearance.datacie, pedited->colorappearance.datacie);
            assignFromKeyfile(keyFile, "Color appearance", "Tonecie", pedited, colorappearance.tonecie, pedited->colorappearance.tonecie);

            const std::map<std::string, ColorAppearanceParams::TcMode> tc_mapping = {
                {"Lightness", ColorAppearanceParams::TcMode::LIGHT},
                {"Brightness", ColorAppearanceParams::TcMode::BRIGHT}
            };
            assignFromKeyfile(keyFile, "Color appearance", "CurveMode", pedited, tc_mapping, colorappearance.curveMode, pedited->colorappearance.curveMode);
            assignFromKeyfile(keyFile, "Color appearance", "CurveMode2", pedited, tc_mapping, colorappearance.curveMode2, pedited->colorappearance.curveMode2);

            assignFromKeyfile(
                keyFile,
                "Color appearance",
                "CurveMode3",
            pedited, {
                {"Chroma", ColorAppearanceParams::CtcMode::CHROMA},
                {"Saturation", ColorAppearanceParams::CtcMode::SATUR},
                {"Colorfullness", ColorAppearanceParams::CtcMode::COLORF}
            },
            colorappearance.curveMode3,
            pedited->colorappearance.curveMode3
            );

            if (ppVersion > 200) {
                assignFromKeyfile(keyFile, "Color appearance", "Curve", pedited, colorappearance.curve, pedited->colorappearance.curve);
                assignFromKeyfile(keyFile, "Color appearance", "Curve2", pedited, colorappearance.curve2, pedited->colorappearance.curve2);
                assignFromKeyfile(keyFile, "Color appearance", "Curve3", pedited, colorappearance.curve3, pedited->colorappearance.curve3);
            }

        }

        if (keyFile.has_group("Impulse Denoising")) {
            assignFromKeyfile(keyFile, "Impulse Denoising", "Enabled", pedited, impulseDenoise.enabled, pedited->impulseDenoise.enabled);
            assignFromKeyfile(keyFile, "Impulse Denoising", "Threshold", pedited, impulseDenoise.thresh, pedited->impulseDenoise.thresh);
        }

        if (keyFile.has_group("Directional Pyramid Denoising")) { //TODO: No longer an accurate description for FT denoise
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Enabled", pedited, dirpyrDenoise.enabled, pedited->dirpyrDenoise.enabled);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Enhance", pedited, dirpyrDenoise.enhance, pedited->dirpyrDenoise.enhance);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Median", pedited, dirpyrDenoise.median, pedited->dirpyrDenoise.median);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Luma", pedited, dirpyrDenoise.luma, pedited->dirpyrDenoise.luma);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Ldetail", pedited, dirpyrDenoise.Ldetail, pedited->dirpyrDenoise.Ldetail);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Chroma", pedited, dirpyrDenoise.chroma, pedited->dirpyrDenoise.chroma);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Method", pedited, dirpyrDenoise.dmethod, pedited->dirpyrDenoise.dmethod);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "LMethod", pedited, dirpyrDenoise.Lmethod, pedited->dirpyrDenoise.Lmethod);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "CMethod", pedited, dirpyrDenoise.Cmethod, pedited->dirpyrDenoise.Cmethod);

            if (dirpyrDenoise.Cmethod == "PRE") {
                dirpyrDenoise.Cmethod = "MAN"; // Never load 'auto chroma preview mode' from pp3
            }

            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "C2Method", pedited, dirpyrDenoise.C2method, pedited->dirpyrDenoise.C2method);

            if (dirpyrDenoise.C2method == "PREV") {
                dirpyrDenoise.C2method = "MANU";
            }

            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "SMethod", pedited, dirpyrDenoise.smethod, pedited->dirpyrDenoise.smethod);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "MedMethod", pedited, dirpyrDenoise.medmethod, pedited->dirpyrDenoise.medmethod);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "MethodMed", pedited, dirpyrDenoise.methodmed, pedited->dirpyrDenoise.methodmed);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "RGBMethod", pedited, dirpyrDenoise.rgbmethod, pedited->dirpyrDenoise.rgbmethod);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "LCurve", pedited, dirpyrDenoise.lcurve, pedited->dirpyrDenoise.lcurve);

            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "CCCurve", pedited, dirpyrDenoise.cccurve, pedited->dirpyrDenoise.cccurve);

            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Redchro", pedited, dirpyrDenoise.redchro, pedited->dirpyrDenoise.redchro);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Bluechro", pedited, dirpyrDenoise.bluechro, pedited->dirpyrDenoise.bluechro);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Gamma", pedited, dirpyrDenoise.gamma, pedited->dirpyrDenoise.gamma);
            assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Passes", pedited, dirpyrDenoise.passes, pedited->dirpyrDenoise.passes);
        }

        if (keyFile.has_group("EPD")) {
            assignFromKeyfile(keyFile, "EPD", "Enabled", pedited, epd.enabled, pedited->epd.enabled);
            assignFromKeyfile(keyFile, "EPD", "Strength", pedited, epd.strength, pedited->epd.strength);
            assignFromKeyfile(keyFile, "EPD", "Gamma", pedited, epd.gamma, pedited->epd.gamma);
            assignFromKeyfile(keyFile, "EPD", "EdgeStopping", pedited, epd.edgeStopping, pedited->epd.edgeStopping);
            assignFromKeyfile(keyFile, "EPD", "Scale", pedited, epd.scale, pedited->epd.scale);
            assignFromKeyfile(keyFile, "EPD", "ReweightingIterates", pedited, epd.reweightingIterates, pedited->epd.reweightingIterates);
        }

        if (keyFile.has_group("FattalToneMapping")) {
            assignFromKeyfile(keyFile, "FattalToneMapping", "Enabled", pedited, fattal.enabled, pedited->fattal.enabled);
            assignFromKeyfile(keyFile, "FattalToneMapping", "Threshold", pedited, fattal.threshold, pedited->fattal.threshold);
            assignFromKeyfile(keyFile, "FattalToneMapping", "Amount", pedited, fattal.amount, pedited->fattal.amount);
            assignFromKeyfile(keyFile, "FattalToneMapping", "Anchor", pedited, fattal.anchor, pedited->fattal.anchor);
        }

        if (keyFile.has_group ("Shadows & Highlights") && ppVersion >= 333) {
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Enabled", pedited, sh.enabled, pedited->sh.enabled);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Highlights", pedited, sh.highlights, pedited->sh.highlights);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "HighlightTonalWidth", pedited, sh.htonalwidth, pedited->sh.htonalwidth);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Shadows", pedited, sh.shadows, pedited->sh.shadows);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "ShadowTonalWidth", pedited, sh.stonalwidth, pedited->sh.stonalwidth);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Radius", pedited, sh.radius, pedited->sh.radius);

            if (keyFile.has_key("Shadows & Highlights", "LocalContrast") && ppVersion < 329) {
                int lc = keyFile.get_integer("Shadows & Highlights", "LocalContrast");
                localContrast.amount = float(lc) / 30.;

                if (pedited) {
                    pedited->localContrast.amount = true;
                }

                localContrast.enabled = sh.enabled;

                if (pedited) {
                    pedited->localContrast.enabled = true;
                }

                localContrast.radius = sh.radius;

                if (pedited) {
                    pedited->localContrast.radius = true;
                }
            }
        }

        if (keyFile.has_group("Crop")) {
            assignFromKeyfile(keyFile, "Crop", "Enabled", pedited, crop.enabled, pedited->crop.enabled);
            assignFromKeyfile(keyFile, "Crop", "X", pedited, crop.x, pedited->crop.x);
            assignFromKeyfile(keyFile, "Crop", "Y", pedited, crop.y, pedited->crop.y);

            if (keyFile.has_key("Crop", "W")) {
                crop.w = std::max(keyFile.get_integer("Crop", "W"), 1);

                if (pedited) {
                    pedited->crop.w = true;
                }
            }

            if (keyFile.has_key("Crop", "H")) {
                crop.h = std::max(keyFile.get_integer("Crop", "H"), 1);

                if (pedited) {
                    pedited->crop.h = true;
                }
            }

            assignFromKeyfile(keyFile, "Crop", "FixedRatio", pedited, crop.fixratio, pedited->crop.fixratio);

            if (assignFromKeyfile(keyFile, "Crop", "Ratio", pedited, crop.ratio, pedited->crop.ratio)) {
                //backwards compatibility for crop.ratio
                if (crop.ratio == "DIN") {
                    crop.ratio = "1.414 - DIN EN ISO 216";
                }

                if (crop.ratio == "8.5:11") {
                    crop.ratio = "8.5:11 - US Letter";
                }

                if (crop.ratio == "11:17") {
                    crop.ratio = "11:17 - Tabloid";
                }
            }

            assignFromKeyfile(keyFile, "Crop", "Orientation", pedited, crop.orientation, pedited->crop.orientation);
            assignFromKeyfile(keyFile, "Crop", "Guide", pedited, crop.guide, pedited->crop.guide);
        }

        if (keyFile.has_group("Coarse Transformation")) {
            assignFromKeyfile(keyFile, "Coarse Transformation", "Rotate", pedited, coarse.rotate, pedited->coarse.rotate);
            assignFromKeyfile(keyFile, "Coarse Transformation", "HorizontalFlip", pedited, coarse.hflip, pedited->coarse.hflip);
            assignFromKeyfile(keyFile, "Coarse Transformation", "VerticalFlip", pedited, coarse.vflip, pedited->coarse.vflip);
        }

        if (keyFile.has_group("Rotation")) {
            assignFromKeyfile(keyFile, "Rotation", "Degree", pedited, rotate.degree, pedited->rotate.degree);
        }

        if (keyFile.has_group("Common Properties for Transformations")) {
            assignFromKeyfile(keyFile, "Common Properties for Transformations", "AutoFill", pedited, commonTrans.autofill, pedited->commonTrans.autofill);
        }

        if (keyFile.has_group("Distortion")) {
            assignFromKeyfile(keyFile, "Distortion", "Amount", pedited, distortion.amount, pedited->distortion.amount);
        }

        if (keyFile.has_group("LensProfile")) {
            if (keyFile.has_key("LensProfile", "LcMode")) {
                lensProf.lcMode = lensProf.getMethodNumber(keyFile.get_string("LensProfile", "LcMode"));

                if (pedited) {
                    pedited->lensProf.lcMode = true;
                }
            }

            if (keyFile.has_key("LensProfile", "LCPFile")) {
                lensProf.lcpFile = expandRelativePath(fname, "", keyFile.get_string("LensProfile", "LCPFile"));

                if (pedited) {
                    pedited->lensProf.lcpFile = true;
                }

                if (ppVersion < 327 && !lensProf.lcpFile.empty()) {
                    lensProf.lcMode = LensProfParams::LcMode::LCP;
                }
            }

            assignFromKeyfile(keyFile, "LensProfile", "UseDistortion", pedited, lensProf.useDist, pedited->lensProf.useDist);
            assignFromKeyfile(keyFile, "LensProfile", "UseVignette", pedited, lensProf.useVign, pedited->lensProf.useVign);
            assignFromKeyfile(keyFile, "LensProfile", "UseCA", pedited, lensProf.useCA, pedited->lensProf.useCA);

            if (keyFile.has_key("LensProfile", "LFCameraMake")) {
                lensProf.lfCameraMake = keyFile.get_string("LensProfile", "LFCameraMake");

                if (pedited) {
                    pedited->lensProf.lfCameraMake = true;
                }
            }

            if (keyFile.has_key("LensProfile", "LFCameraModel")) {
                lensProf.lfCameraModel = keyFile.get_string("LensProfile", "LFCameraModel");

                if (pedited) {
                    pedited->lensProf.lfCameraModel = true;
                }
            }

            if (keyFile.has_key("LensProfile", "LFLens")) {
                lensProf.lfLens = keyFile.get_string("LensProfile", "LFLens");

                if (pedited) {
                    pedited->lensProf.lfLens = true;
                }
            }
        }

        if (keyFile.has_group("Perspective")) {
            assignFromKeyfile(keyFile, "Perspective", "Horizontal", pedited, perspective.horizontal, pedited->perspective.horizontal);
            assignFromKeyfile(keyFile, "Perspective", "Vertical", pedited, perspective.vertical, pedited->perspective.vertical);
        }

        if (keyFile.has_group("Gradient")) {
            assignFromKeyfile(keyFile, "Gradient", "Enabled", pedited, gradient.enabled, pedited->gradient.enabled);
            assignFromKeyfile(keyFile, "Gradient", "Degree", pedited, gradient.degree, pedited->gradient.degree);
            assignFromKeyfile(keyFile, "Gradient", "Feather", pedited, gradient.feather, pedited->gradient.feather);
            assignFromKeyfile(keyFile, "Gradient", "Strength", pedited, gradient.strength, pedited->gradient.strength);
            assignFromKeyfile(keyFile, "Gradient", "CenterX", pedited, gradient.centerX, pedited->gradient.centerX);
            assignFromKeyfile(keyFile, "Gradient", "CenterY", pedited, gradient.centerY, pedited->gradient.centerY);
        }

        if (keyFile.has_group("PCVignette")) {
            assignFromKeyfile(keyFile, "PCVignette", "Enabled", pedited, pcvignette.enabled, pedited->pcvignette.enabled);
            assignFromKeyfile(keyFile, "PCVignette", "Strength", pedited, pcvignette.strength, pedited->pcvignette.strength);
            assignFromKeyfile(keyFile, "PCVignette", "Feather", pedited, pcvignette.feather, pedited->pcvignette.feather);
            assignFromKeyfile(keyFile, "PCVignette", "Roundness", pedited, pcvignette.roundness, pedited->pcvignette.roundness);
        }

        if (keyFile.has_group("CACorrection")) {
            assignFromKeyfile(keyFile, "CACorrection", "Red", pedited, cacorrection.red, pedited->cacorrection.red);
            assignFromKeyfile(keyFile, "CACorrection", "Blue", pedited, cacorrection.blue, pedited->cacorrection.blue);
        }

        if (keyFile.has_group("Vignetting Correction")) {
            assignFromKeyfile(keyFile, "Vignetting Correction", "Amount", pedited, vignetting.amount, pedited->vignetting.amount);
            assignFromKeyfile(keyFile, "Vignetting Correction", "Radius", pedited, vignetting.radius, pedited->vignetting.radius);
            assignFromKeyfile(keyFile, "Vignetting Correction", "Strength", pedited, vignetting.strength, pedited->vignetting.strength);
            assignFromKeyfile(keyFile, "Vignetting Correction", "CenterX", pedited, vignetting.centerX, pedited->vignetting.centerX);
            assignFromKeyfile(keyFile, "Vignetting Correction", "CenterY", pedited, vignetting.centerY, pedited->vignetting.centerY);
        }

        if (keyFile.has_group("Resize")) {
            assignFromKeyfile(keyFile, "Resize", "Enabled", pedited, resize.enabled, pedited->resize.enabled);
            assignFromKeyfile(keyFile, "Resize", "Scale", pedited, resize.scale, pedited->resize.scale);
            assignFromKeyfile(keyFile, "Resize", "AppliesTo", pedited, resize.appliesTo, pedited->resize.appliesTo);
            assignFromKeyfile(keyFile, "Resize", "Method", pedited, resize.method, pedited->resize.method);
            assignFromKeyfile(keyFile, "Resize", "DataSpecified", pedited, resize.dataspec, pedited->resize.dataspec);
            assignFromKeyfile(keyFile, "Resize", "Width", pedited, resize.width, pedited->resize.width);
            assignFromKeyfile(keyFile, "Resize", "Height", pedited, resize.height, pedited->resize.height);
        }

        if (keyFile.has_group("PostResizeSharpening")) {
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Enabled", pedited, prsharpening.enabled, pedited->prsharpening.enabled);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Radius", pedited, prsharpening.radius, pedited->prsharpening.radius);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Amount", pedited, prsharpening.amount, pedited->prsharpening.amount);

            if (keyFile.has_key("PostResizeSharpening", "Threshold")) {
                if (ppVersion < 302) {
                    int thresh = min(keyFile.get_integer("PostResizeSharpening", "Threshold"), 2000);
                    prsharpening.threshold.setValues(thresh, thresh, 2000, 2000);  // TODO: 2000 is the maximum value and is taken of rtgui/sharpening.cc ; should be changed by the tool modularization
                } else {
                    const std::vector<int> thresh = keyFile.get_integer_list("PostResizeSharpening", "Threshold");

                    if (thresh.size() >= 4) {
                        prsharpening.threshold.setValues(thresh[0], thresh[1], min(thresh[2], 2000), min(thresh[3], 2000));
                    }
                }

                if (pedited) {
                    pedited->prsharpening.threshold = true;
                }
            }

            assignFromKeyfile(keyFile, "PostResizeSharpening", "OnlyEdges", pedited, prsharpening.edgesonly, pedited->prsharpening.edgesonly);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "EdgedetectionRadius", pedited, prsharpening.edges_radius, pedited->prsharpening.edges_radius);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "EdgeTolerance", pedited, prsharpening.edges_tolerance, pedited->prsharpening.edges_tolerance);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "HalocontrolEnabled", pedited, prsharpening.halocontrol, pedited->prsharpening.halocontrol);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "HalocontrolAmount", pedited, prsharpening.halocontrol_amount, pedited->prsharpening.halocontrol_amount);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Method", pedited, prsharpening.method, pedited->prsharpening.method);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvRadius", pedited, prsharpening.deconvradius, pedited->prsharpening.deconvradius);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvAmount", pedited, prsharpening.deconvamount, pedited->prsharpening.deconvamount);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvDamping", pedited, prsharpening.deconvdamping, pedited->prsharpening.deconvdamping);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvIterations", pedited, prsharpening.deconviter, pedited->prsharpening.deconviter);
        }

        if (keyFile.has_group("Color Management")) {
            if (keyFile.has_key("Color Management", "InputProfile")) {
                icm.input = expandRelativePath(fname, "file:", keyFile.get_string("Color Management", "InputProfile"));

                if (pedited) {
                    pedited->icm.input = true;
                }
            }

            assignFromKeyfile(keyFile, "Color Management", "ToneCurve", pedited, icm.toneCurve, pedited->icm.toneCurve);
            assignFromKeyfile(keyFile, "Color Management", "ApplyLookTable", pedited, icm.applyLookTable, pedited->icm.applyLookTable);
            assignFromKeyfile(keyFile, "Color Management", "ApplyBaselineExposureOffset", pedited, icm.applyBaselineExposureOffset, pedited->icm.applyBaselineExposureOffset);
            assignFromKeyfile(keyFile, "Color Management", "ApplyHueSatMap", pedited, icm.applyHueSatMap, pedited->icm.applyHueSatMap);
            assignFromKeyfile(keyFile, "Color Management", "DCPIlluminant", pedited, icm.dcpIlluminant, pedited->icm.dcpIlluminant);
            assignFromKeyfile(keyFile, "Color Management", "WorkingProfile", pedited, icm.working, pedited->icm.working);
            assignFromKeyfile(keyFile, "Color Management", "OutputProfile", pedited, icm.output, pedited->icm.output);

            if (keyFile.has_key("Color Management", "OutputProfileIntent")) {
                Glib::ustring intent = keyFile.get_string("Color Management", "OutputProfileIntent");

                if (intent == "Perceptual") {
                    icm.outputIntent = RI_PERCEPTUAL;
                } else if (intent == "Relative") {
                    icm.outputIntent = RI_RELATIVE;
                } else if (intent == "Saturation") {
                    icm.outputIntent = RI_SATURATION;
                } else if (intent == "Absolute") {
                    icm.outputIntent = RI_ABSOLUTE;
                }

                if (pedited) {
                    pedited->icm.outputIntent = true;
                }
            }

            assignFromKeyfile(keyFile, "Color Management", "OutputBPC", pedited, icm.outputBPC, pedited->icm.outputBPC);
            assignFromKeyfile(keyFile, "Color Management", "GammaCustom", pedited, icm.gamma, pedited->icm.gamma);
            assignFromKeyfile(keyFile, "Color Management", "Freegamma", pedited, icm.freegamma, pedited->icm.freegamma);
            assignFromKeyfile(keyFile, "Color Management", "GammaValue", pedited, icm.gampos, pedited->icm.gampos);
            assignFromKeyfile(keyFile, "Color Management", "GammaSlope", pedited, icm.slpos, pedited->icm.slpos);
            assignFromKeyfile(keyFile, "Color Management", "GammaPredx", pedited, icm.predx, pedited->icm.predx);
            assignFromKeyfile(keyFile, "Color Management", "GammaPredy", pedited, icm.predy, pedited->icm.predy);
            assignFromKeyfile(keyFile, "Color Management", "GammaPgrex", pedited, icm.pgrex, pedited->icm.pgrex);
            assignFromKeyfile(keyFile, "Color Management", "GammaPgrey", pedited, icm.pgrey, pedited->icm.pgrey);
            assignFromKeyfile(keyFile, "Color Management", "GammaPblux", pedited, icm.pblux, pedited->icm.pblux);
            assignFromKeyfile(keyFile, "Color Management", "GammaPbluy", pedited, icm.pbluy, pedited->icm.pbluy);

            assignFromKeyfile(keyFile, "Color Management", "GammaValueIn", pedited, icm.gamm, pedited->icm.gamm);
            assignFromKeyfile(keyFile, "Color Management", "GammaSlopeIn", pedited, icm.slop, pedited->icm.slop);
			
            assignFromKeyfile(keyFile, "Color Management", "GammaPrimaries", pedited, icm.wprimaries, pedited->icm.wprimaries);
            assignFromKeyfile(keyFile, "Color Management", "GammaProfile", pedited, icm.wprofile, pedited->icm.wprofile);
            assignFromKeyfile(keyFile, "Color Management", "GammaTemp", pedited, icm.wtemp, pedited->icm.wtemp);
            assignFromKeyfile(keyFile, "Color Management", "GammaTRCIN", pedited, icm.wtrcin, pedited->icm.wtrcin);
        }

        if (keyFile.has_group("Wavelet")) {
            assignFromKeyfile(keyFile, "Wavelet", "Enabled", pedited, wavelet.enabled, pedited->wavelet.enabled);
            assignFromKeyfile(keyFile, "Wavelet", "Strength", pedited, wavelet.strength, pedited->wavelet.strength);
            assignFromKeyfile(keyFile, "Wavelet", "Balance", pedited, wavelet.balance, pedited->wavelet.balance);
            assignFromKeyfile(keyFile, "Wavelet", "Iter", pedited, wavelet.iter, pedited->wavelet.iter);
            assignFromKeyfile(keyFile, "Wavelet", "Median", pedited, wavelet.median, pedited->wavelet.median);
            assignFromKeyfile(keyFile, "Wavelet", "Medianlev", pedited, wavelet.medianlev, pedited->wavelet.medianlev);
            assignFromKeyfile(keyFile, "Wavelet", "Linkedg", pedited, wavelet.linkedg, pedited->wavelet.linkedg);
            assignFromKeyfile(keyFile, "Wavelet", "CBenab", pedited, wavelet.cbenab, pedited->wavelet.cbenab);
            assignFromKeyfile(keyFile, "Wavelet", "CBgreenhigh", pedited, wavelet.greenhigh, pedited->wavelet.greenhigh);
            assignFromKeyfile(keyFile, "Wavelet", "CBgreenmed", pedited, wavelet.greenmed, pedited->wavelet.greenmed);
            assignFromKeyfile(keyFile, "Wavelet", "CBgreenlow", pedited, wavelet.greenlow, pedited->wavelet.greenlow);
            assignFromKeyfile(keyFile, "Wavelet", "CBbluehigh", pedited, wavelet.bluehigh, pedited->wavelet.bluehigh);
            assignFromKeyfile(keyFile, "Wavelet", "CBbluemed", pedited, wavelet.bluemed, pedited->wavelet.bluemed);
            assignFromKeyfile(keyFile, "Wavelet", "CBbluelow", pedited, wavelet.bluelow, pedited->wavelet.bluelow);
            assignFromKeyfile(keyFile, "Wavelet", "Lipst", pedited, wavelet.lipst, pedited->wavelet.lipst);
            assignFromKeyfile(keyFile, "Wavelet", "AvoidColorShift", pedited, wavelet.avoid, pedited->wavelet.avoid);
            assignFromKeyfile(keyFile, "Wavelet", "TMr", pedited, wavelet.tmr, pedited->wavelet.tmr);

            if (ppVersion < 331) { // wavelet.Lmethod was a string before version 331
                Glib::ustring temp;
                assignFromKeyfile(keyFile, "Wavelet", "LevMethod", pedited, temp, pedited->wavelet.Lmethod);

                if (!temp.empty()) {
                    wavelet.Lmethod = std::stoi(temp);
                }
            } else {
                assignFromKeyfile(keyFile, "Wavelet", "LevMethod", pedited, wavelet.Lmethod, pedited->wavelet.Lmethod);
            }

            assignFromKeyfile(keyFile, "Wavelet", "ChoiceLevMethod", pedited, wavelet.CLmethod, pedited->wavelet.CLmethod);
            assignFromKeyfile(keyFile, "Wavelet", "BackMethod", pedited, wavelet.Backmethod, pedited->wavelet.Backmethod);
            assignFromKeyfile(keyFile, "Wavelet", "TilesMethod", pedited, wavelet.Tilesmethod, pedited->wavelet.Tilesmethod);
            assignFromKeyfile(keyFile, "Wavelet", "DaubMethod", pedited, wavelet.daubcoeffmethod, pedited->wavelet.daubcoeffmethod);
            assignFromKeyfile(keyFile, "Wavelet", "CHromaMethod", pedited, wavelet.CHmethod, pedited->wavelet.CHmethod);
            assignFromKeyfile(keyFile, "Wavelet", "Medgreinf", pedited, wavelet.Medgreinf, pedited->wavelet.Medgreinf);
            assignFromKeyfile(keyFile, "Wavelet", "CHSLromaMethod", pedited, wavelet.CHSLmethod, pedited->wavelet.CHSLmethod);
            assignFromKeyfile(keyFile, "Wavelet", "EDMethod", pedited, wavelet.EDmethod, pedited->wavelet.EDmethod);
            assignFromKeyfile(keyFile, "Wavelet", "NPMethod", pedited, wavelet.NPmethod, pedited->wavelet.NPmethod);
            assignFromKeyfile(keyFile, "Wavelet", "BAMethod", pedited, wavelet.BAmethod, pedited->wavelet.BAmethod);
            assignFromKeyfile(keyFile, "Wavelet", "TMMethod", pedited, wavelet.TMmethod, pedited->wavelet.TMmethod);
            assignFromKeyfile(keyFile, "Wavelet", "HSMethod", pedited, wavelet.HSmethod, pedited->wavelet.HSmethod);
            assignFromKeyfile(keyFile, "Wavelet", "DirMethod", pedited, wavelet.Dirmethod, pedited->wavelet.Dirmethod);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualcontShadow", pedited, wavelet.rescon, pedited->wavelet.rescon);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualcontHighlight", pedited, wavelet.resconH, pedited->wavelet.resconH);
            assignFromKeyfile(keyFile, "Wavelet", "Residualchroma", pedited, wavelet.reschro, pedited->wavelet.reschro);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualTM", pedited, wavelet.tmrs, pedited->wavelet.tmrs);
            assignFromKeyfile(keyFile, "Wavelet", "Residualgamma", pedited, wavelet.gamma, pedited->wavelet.gamma);
            assignFromKeyfile(keyFile, "Wavelet", "ContExtra", pedited, wavelet.sup, pedited->wavelet.sup);
            assignFromKeyfile(keyFile, "Wavelet", "HueRangeResidual", pedited, wavelet.sky, pedited->wavelet.sky);
            assignFromKeyfile(keyFile, "Wavelet", "MaxLev", pedited, wavelet.thres, pedited->wavelet.thres);
            assignFromKeyfile(keyFile, "Wavelet", "ThresholdHighlight", pedited, wavelet.threshold, pedited->wavelet.threshold);
            assignFromKeyfile(keyFile, "Wavelet", "ThresholdShadow", pedited, wavelet.threshold2, pedited->wavelet.threshold2);
            assignFromKeyfile(keyFile, "Wavelet", "Edgedetect", pedited, wavelet.edgedetect, pedited->wavelet.edgedetect);
            assignFromKeyfile(keyFile, "Wavelet", "Edgedetectthr", pedited, wavelet.edgedetectthr, pedited->wavelet.edgedetectthr);
            assignFromKeyfile(keyFile, "Wavelet", "EdgedetectthrHi", pedited, wavelet.edgedetectthr2, pedited->wavelet.edgedetectthr2);
            assignFromKeyfile(keyFile, "Wavelet", "Edgesensi", pedited, wavelet.edgesensi, pedited->wavelet.edgesensi);
            assignFromKeyfile(keyFile, "Wavelet", "Edgeampli", pedited, wavelet.edgeampli, pedited->wavelet.edgeampli);
            assignFromKeyfile(keyFile, "Wavelet", "ThresholdChroma", pedited, wavelet.chroma, pedited->wavelet.chroma);
            assignFromKeyfile(keyFile, "Wavelet", "ChromaLink", pedited, wavelet.chro, pedited->wavelet.chro);
            assignFromKeyfile(keyFile, "Wavelet", "Contrast", pedited, wavelet.contrast, pedited->wavelet.contrast);
            assignFromKeyfile(keyFile, "Wavelet", "Edgrad", pedited, wavelet.edgrad, pedited->wavelet.edgrad);
            assignFromKeyfile(keyFile, "Wavelet", "Edgval", pedited, wavelet.edgval, pedited->wavelet.edgval);
            assignFromKeyfile(keyFile, "Wavelet", "ThrEdg", pedited, wavelet.edgthresh, pedited->wavelet.edgthresh);
            assignFromKeyfile(keyFile, "Wavelet", "ThresholdResidShadow", pedited, wavelet.thr, pedited->wavelet.thr);
            assignFromKeyfile(keyFile, "Wavelet", "ThresholdResidHighLight", pedited, wavelet.thrH, pedited->wavelet.thrH);
            assignFromKeyfile(keyFile, "Wavelet", "ContrastCurve", pedited, wavelet.ccwcurve, pedited->wavelet.ccwcurve);
            assignFromKeyfile(keyFile, "Wavelet", "OpacityCurveRG", pedited, wavelet.opacityCurveRG, pedited->wavelet.opacityCurveRG);
            assignFromKeyfile(keyFile, "Wavelet", "OpacityCurveBY", pedited, wavelet.opacityCurveBY, pedited->wavelet.opacityCurveBY);
            assignFromKeyfile(keyFile, "Wavelet", "OpacityCurveW", pedited, wavelet.opacityCurveW, pedited->wavelet.opacityCurveW);
            assignFromKeyfile(keyFile, "Wavelet", "OpacityCurveWL", pedited, wavelet.opacityCurveWL, pedited->wavelet.opacityCurveWL);
            assignFromKeyfile(keyFile, "Wavelet", "HHcurve", pedited, wavelet.hhcurve, pedited->wavelet.hhcurve);
            assignFromKeyfile(keyFile, "Wavelet", "CHcurve", pedited, wavelet.Chcurve, pedited->wavelet.Chcurve);
            assignFromKeyfile(keyFile, "Wavelet", "WavclCurve", pedited, wavelet.wavclCurve, pedited->wavelet.wavclCurve);

            if (keyFile.has_key("Wavelet", "Hueskin")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Wavelet", "Hueskin");

                if (thresh.size() >= 4) {
                    wavelet.hueskin.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->wavelet.hueskin = true;
                }
            }

            if (keyFile.has_key("Wavelet", "HueRange")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Wavelet", "HueRange");

                if (thresh.size() >= 4) {
                    wavelet.hueskin2.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->wavelet.hueskin2 = true;
                }
            }

            if (keyFile.has_key("Wavelet", "HLRange")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Wavelet", "HLRange");

                if (thresh.size() >= 4) {
                    wavelet.hllev.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->wavelet.hllev = true;
                }
            }

            if (keyFile.has_key("Wavelet", "SHRange")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Wavelet", "SHRange");

                if (thresh.size() >= 4) {
                    wavelet.bllev.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->wavelet.bllev = true;
                }
            }

            if (keyFile.has_key("Wavelet", "Edgcont")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Wavelet", "Edgcont");

                if (thresh.size() >= 4) {
                    wavelet.edgcont.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->wavelet.edgcont = true;
                }
            }

            if (keyFile.has_key("Wavelet", "Level0noise")) {
                const std::vector<double> thresh = keyFile.get_double_list("Wavelet", "Level0noise");

                if (thresh.size() >= 2) {
                    wavelet.level0noise.setValues(thresh[0], thresh[1]);
                }

                if (pedited) {
                    pedited->wavelet.level0noise = true;
                }
            }

            if (keyFile.has_key("Wavelet", "Level1noise")) {
                const std::vector<double> thresh = keyFile.get_double_list("Wavelet", "Level1noise");

                if (thresh.size() >= 2) {
                    wavelet.level1noise.setValues(thresh[0], thresh[1]);
                }

                if (pedited) {
                    pedited->wavelet.level1noise = true;
                }
            }

            if (keyFile.has_key("Wavelet", "Level2noise")) {
                const std::vector<double> thresh = keyFile.get_double_list("Wavelet", "Level2noise");

                if (thresh.size() >= 2) {
                    wavelet.level2noise.setValues(thresh[0], thresh[1]);
                }

                if (pedited) {
                    pedited->wavelet.level2noise = true;
                }
            }

            if (keyFile.has_key("Wavelet", "Level3noise")) {
                const std::vector<double> thresh = keyFile.get_double_list("Wavelet", "Level3noise");

                if (thresh.size() >= 2) {
                    wavelet.level3noise.setValues(thresh[0], thresh[1]);
                }

                if (pedited) {
                    pedited->wavelet.level3noise = true;
                }
            }

            if (keyFile.has_key("Wavelet", "Pastlev")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Wavelet", "Pastlev");

                if (thresh.size() >= 4) {
                    wavelet.pastlev.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->wavelet.pastlev = true;
                }
            }

            if (keyFile.has_key("Wavelet", "Satlev")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Wavelet", "Satlev");

                if (thresh.size() >= 4) {
                    wavelet.satlev.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->wavelet.satlev = true;
                }
            }

            assignFromKeyfile(keyFile, "Wavelet", "Skinprotect", pedited, wavelet.skinprotect, pedited->wavelet.skinprotect);
            assignFromKeyfile(keyFile, "Wavelet", "Expcontrast", pedited, wavelet.expcontrast, pedited->wavelet.expcontrast);
            assignFromKeyfile(keyFile, "Wavelet", "Expchroma", pedited, wavelet.expchroma, pedited->wavelet.expchroma);

            for (int i = 0; i < 9; ++i) {
                std::stringstream ss;
                ss << "Contrast" << (i + 1);

                if (keyFile.has_key("Wavelet", ss.str())) {
                    wavelet.c[i] = keyFile.get_integer("Wavelet", ss.str());

                    if (pedited) {
                        pedited->wavelet.c[i] = true;
                    }
                }
            }

            for (int i = 0; i < 9; ++i) {
                std::stringstream ss;
                ss << "Chroma" << (i + 1);

                if (keyFile.has_key("Wavelet", ss.str())) {
                    wavelet.ch[i] = keyFile.get_integer("Wavelet", ss.str());

                    if (pedited) {
                        pedited->wavelet.ch[i] = true;
                    }
                }
            }

            assignFromKeyfile(keyFile, "Wavelet", "Expedge", pedited, wavelet.expedge, pedited->wavelet.expedge);
            assignFromKeyfile(keyFile, "Wavelet", "Expresid", pedited, wavelet.expresid, pedited->wavelet.expresid);
            assignFromKeyfile(keyFile, "Wavelet", "Expfinal", pedited, wavelet.expfinal, pedited->wavelet.expfinal);
            assignFromKeyfile(keyFile, "Wavelet", "Exptoning", pedited, wavelet.exptoning, pedited->wavelet.exptoning);
            assignFromKeyfile(keyFile, "Wavelet", "Expnoise", pedited, wavelet.expnoise, pedited->wavelet.expnoise);
        }

        if (keyFile.has_group("Directional Pyramid Equalizer")) {
            assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Enabled", pedited, dirpyrequalizer.enabled, pedited->dirpyrequalizer.enabled);
            assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Gamutlab", pedited, dirpyrequalizer.gamutlab, pedited->dirpyrequalizer.gamutlab);
            assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "cbdlMethod", pedited, dirpyrequalizer.cbdlMethod, pedited->dirpyrequalizer.cbdlMethod);

            if (keyFile.has_key("Directional Pyramid Equalizer", "Hueskin")) {
                const std::vector<int> thresh = keyFile.get_integer_list("Directional Pyramid Equalizer", "Hueskin");

                if (thresh.size() >= 4) {
                    dirpyrequalizer.hueskin.setValues(thresh[0], thresh[1], min(thresh[2], 300), min(thresh[3], 300));
                }

                if (pedited) {
                    pedited->dirpyrequalizer.hueskin = true;
                }
            }

            if (ppVersion < 316) {
                for (int i = 0; i < 5; i ++) {
                    std::stringstream ss;
                    ss << "Mult" << i;

                    if (keyFile.has_key("Directional Pyramid Equalizer", ss.str())) {
                        if (i == 4) {
                            dirpyrequalizer.threshold = keyFile.get_double("Directional Pyramid Equalizer", ss.str());

                            if (pedited) {
                                pedited->dirpyrequalizer.threshold = true;
                            }
                        } else {
                            dirpyrequalizer.mult[i] = keyFile.get_double("Directional Pyramid Equalizer", ss.str());

                            if (pedited) {
                                pedited->dirpyrequalizer.mult[i] = true;
                            }
                        }
                    }
                }

                dirpyrequalizer.mult[4] = 1.0;
            } else {
                // 5 level wavelet + dedicated threshold parameter
                for (int i = 0; i < 6; i ++) {
                    std::stringstream ss;
                    ss << "Mult" << i;

                    if (keyFile.has_key("Directional Pyramid Equalizer", ss.str())) {
                        dirpyrequalizer.mult[i] = keyFile.get_double("Directional Pyramid Equalizer", ss.str());

                        if (pedited) {
                            pedited->dirpyrequalizer.mult[i] = true;
                        }
                    }
                }

                assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Threshold", pedited, dirpyrequalizer.threshold, pedited->dirpyrequalizer.threshold);
                assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Skinprotect", pedited, dirpyrequalizer.skinprotect, pedited->dirpyrequalizer.skinprotect);
            }
        }

        if (keyFile.has_group("Film Simulation")) {
            assignFromKeyfile(keyFile, "Film Simulation", "Enabled", pedited, filmSimulation.enabled, pedited->filmSimulation.enabled);
            assignFromKeyfile(keyFile, "Film Simulation", "ClutFilename", pedited, filmSimulation.clutFilename, pedited->filmSimulation.clutFilename);

            if (keyFile.has_key("Film Simulation", "Strength")) {
                if (ppVersion < 321) {
                    filmSimulation.strength = keyFile.get_double("Film Simulation", "Strength") * 100 + 0.1;
                } else {
                    filmSimulation.strength = keyFile.get_integer("Film Simulation", "Strength");
                }

                if (pedited) {
                    pedited->filmSimulation.strength = true;
                }
            }
        }

        if (keyFile.has_group("HSV Equalizer")) {
            if (ppVersion >= 329) {
                assignFromKeyfile(keyFile, "HSV Equalizer", "Enabled", pedited, hsvequalizer.enabled, pedited->hsvequalizer.enabled);
            } else {
                hsvequalizer.enabled = true;

                if (pedited) {
                    pedited->hsvequalizer.enabled = true;
                }
            }

            if (ppVersion >= 300) {
                assignFromKeyfile(keyFile, "HSV Equalizer", "HCurve", pedited, hsvequalizer.hcurve, pedited->hsvequalizer.hcurve);
                assignFromKeyfile(keyFile, "HSV Equalizer", "SCurve", pedited, hsvequalizer.scurve, pedited->hsvequalizer.scurve);
                assignFromKeyfile(keyFile, "HSV Equalizer", "VCurve", pedited, hsvequalizer.vcurve, pedited->hsvequalizer.vcurve);
            }
        }

        if (keyFile.has_group("RGB Curves")) {
            if (ppVersion >= 329) {
                assignFromKeyfile(keyFile, "RGB Curves", "Enabled", pedited, rgbCurves.enabled, pedited->rgbCurves.enabled);
            } else {
                rgbCurves.enabled = true;

                if (pedited) {
                    pedited->rgbCurves.enabled = true;
                }
            }

            assignFromKeyfile(keyFile, "RGB Curves", "LumaMode", pedited, rgbCurves.lumamode, pedited->rgbCurves.lumamode);
            assignFromKeyfile(keyFile, "RGB Curves", "rCurve", pedited, rgbCurves.rcurve, pedited->rgbCurves.rcurve);
            assignFromKeyfile(keyFile, "RGB Curves", "gCurve", pedited, rgbCurves.gcurve, pedited->rgbCurves.gcurve);
            assignFromKeyfile(keyFile, "RGB Curves", "bCurve", pedited, rgbCurves.bcurve, pedited->rgbCurves.bcurve);
        }

        if (keyFile.has_group("ColorToning")) {
            assignFromKeyfile(keyFile, "ColorToning", "Enabled", pedited, colorToning.enabled, pedited->colorToning.enabled);
            assignFromKeyfile(keyFile, "ColorToning", "Method", pedited, colorToning.method, pedited->colorToning.method);
            assignFromKeyfile(keyFile, "ColorToning", "Lumamode", pedited, colorToning.lumamode, pedited->colorToning.lumamode);
            assignFromKeyfile(keyFile, "ColorToning", "Twocolor", pedited, colorToning.twocolor, pedited->colorToning.twocolor);
            assignFromKeyfile(keyFile, "ColorToning", "OpacityCurve", pedited, colorToning.opacityCurve, pedited->colorToning.opacityCurve);
            assignFromKeyfile(keyFile, "ColorToning", "ColorCurve", pedited, colorToning.colorCurve, pedited->colorToning.colorCurve);
            assignFromKeyfile(keyFile, "ColorToning", "Autosat", pedited, colorToning.autosat, pedited->colorToning.autosat);
            assignFromKeyfile(keyFile, "ColorToning", "SatProtectionThreshold", pedited, colorToning.satProtectionThreshold, pedited->colorToning.satprotectionthreshold);
            assignFromKeyfile(keyFile, "ColorToning", "SaturatedOpacity", pedited, colorToning.saturatedOpacity, pedited->colorToning.saturatedopacity);
            assignFromKeyfile(keyFile, "ColorToning", "Strength", pedited, colorToning.strength, pedited->colorToning.strength);

            if (keyFile.has_key("ColorToning", "HighlightsColorSaturation")) {
                const std::vector<int> thresh = keyFile.get_integer_list("ColorToning", "HighlightsColorSaturation");

                if (thresh.size() >= 2) {
                    colorToning.hlColSat.setValues(thresh[0], thresh[1]);
                }

                if (pedited) {
                    pedited->colorToning.hlColSat = true;
                }
            }

            if (keyFile.has_key("ColorToning", "ShadowsColorSaturation")) {
                const std::vector<int> thresh = keyFile.get_integer_list("ColorToning", "ShadowsColorSaturation");

                if (thresh.size() >= 2) {
                    colorToning.shadowsColSat.setValues(thresh[0], thresh[1]);
                }

                if (pedited) {
                    pedited->colorToning.shadowsColSat = true;
                }
            }

            assignFromKeyfile(keyFile, "ColorToning", "ClCurve", pedited, colorToning.clcurve, pedited->colorToning.clcurve);
            assignFromKeyfile(keyFile, "ColorToning", "Cl2Curve", pedited, colorToning.cl2curve, pedited->colorToning.cl2curve);
            assignFromKeyfile(keyFile, "ColorToning", "Redlow", pedited, colorToning.redlow, pedited->colorToning.redlow);
            assignFromKeyfile(keyFile, "ColorToning", "Greenlow", pedited, colorToning.greenlow, pedited->colorToning.greenlow);
            assignFromKeyfile(keyFile, "ColorToning", "Bluelow", pedited, colorToning.bluelow, pedited->colorToning.bluelow);
            assignFromKeyfile(keyFile, "ColorToning", "Satlow", pedited, colorToning.satlow, pedited->colorToning.satlow);
            assignFromKeyfile(keyFile, "ColorToning", "Balance", pedited, colorToning.balance, pedited->colorToning.balance);
            assignFromKeyfile(keyFile, "ColorToning", "Sathigh", pedited, colorToning.sathigh, pedited->colorToning.sathigh);
            assignFromKeyfile(keyFile, "ColorToning", "Redmed", pedited, colorToning.redmed, pedited->colorToning.redmed);
            assignFromKeyfile(keyFile, "ColorToning", "Greenmed", pedited, colorToning.greenmed, pedited->colorToning.greenmed);
            assignFromKeyfile(keyFile, "ColorToning", "Bluemed", pedited, colorToning.bluemed, pedited->colorToning.bluemed);
            assignFromKeyfile(keyFile, "ColorToning", "Redhigh", pedited, colorToning.redhigh, pedited->colorToning.redhigh);
            assignFromKeyfile(keyFile, "ColorToning", "Greenhigh", pedited, colorToning.greenhigh, pedited->colorToning.greenhigh);
            assignFromKeyfile(keyFile, "ColorToning", "Bluehigh", pedited, colorToning.bluehigh, pedited->colorToning.bluehigh);

            assignFromKeyfile(keyFile, "ColorToning", "LabGridALow", pedited, colorToning.labgridALow, pedited->colorToning.labgridALow);
            assignFromKeyfile(keyFile, "ColorToning", "LabGridBLow", pedited, colorToning.labgridBLow, pedited->colorToning.labgridBLow);
            assignFromKeyfile(keyFile, "ColorToning", "LabGridAHigh", pedited, colorToning.labgridAHigh, pedited->colorToning.labgridAHigh);
            assignFromKeyfile(keyFile, "ColorToning", "LabGridBHigh", pedited, colorToning.labgridBHigh, pedited->colorToning.labgridBHigh);
        }

        if (keyFile.has_group("RAW")) {
            if (keyFile.has_key("RAW", "DarkFrame")) {
                raw.dark_frame = expandRelativePath(fname, "", keyFile.get_string("RAW", "DarkFrame"));

                if (pedited) {
                    pedited->raw.darkFrame = true;
                }
            }

            assignFromKeyfile(keyFile, "RAW", "DarkFrameAuto", pedited, raw.df_autoselect, pedited->raw.df_autoselect);

            if (keyFile.has_key("RAW", "FlatFieldFile")) {
                raw.ff_file = expandRelativePath(fname, "", keyFile.get_string("RAW", "FlatFieldFile"));

                if (pedited) {
                    pedited->raw.ff_file = true;
                }
            }

            assignFromKeyfile(keyFile, "RAW", "FlatFieldAutoSelect", pedited, raw.ff_AutoSelect, pedited->raw.ff_AutoSelect);
            assignFromKeyfile(keyFile, "RAW", "FlatFieldBlurRadius", pedited, raw.ff_BlurRadius, pedited->raw.ff_BlurRadius);
            assignFromKeyfile(keyFile, "RAW", "FlatFieldBlurType", pedited, raw.ff_BlurType, pedited->raw.ff_BlurType);
            assignFromKeyfile(keyFile, "RAW", "FlatFieldAutoClipControl", pedited, raw.ff_AutoClipControl, pedited->raw.ff_AutoClipControl);

            if (ppVersion < 328) {
                // With ppversion < 328 this value was stored as a boolean, which is nonsense.
                // To avoid annoying warnings we skip reading and assume 0.
                raw.ff_clipControl = 0;
            } else {
                assignFromKeyfile(keyFile, "RAW", "FlatFieldClipControl", pedited, raw.ff_clipControl, pedited->raw.ff_clipControl);
            }

            assignFromKeyfile(keyFile, "RAW", "CA", pedited, raw.ca_autocorrect, pedited->raw.ca_autocorrect);
            assignFromKeyfile(keyFile, "RAW", "CARed", pedited, raw.cared, pedited->raw.cared);
            assignFromKeyfile(keyFile, "RAW", "CABlue", pedited, raw.cablue, pedited->raw.cablue);
            // For compatibility to elder pp3 versions
            assignFromKeyfile(keyFile, "RAW", "HotDeadPixels", pedited, raw.hotPixelFilter, pedited->raw.hotPixelFilter);
            raw.deadPixelFilter = raw.hotPixelFilter;

            if (pedited) {
                pedited->raw.deadPixelFilter = pedited->raw.hotPixelFilter;
            }

            assignFromKeyfile(keyFile, "RAW", "HotPixelFilter", pedited, raw.hotPixelFilter, pedited->raw.hotPixelFilter);
            assignFromKeyfile(keyFile, "RAW", "DeadPixelFilter", pedited, raw.deadPixelFilter, pedited->raw.deadPixelFilter);
            assignFromKeyfile(keyFile, "RAW", "HotDeadPixelThresh", pedited, raw.hotdeadpix_thresh, pedited->raw.hotdeadpix_thresh);
            assignFromKeyfile(keyFile, "RAW", "PreExposure", pedited, raw.expos, pedited->raw.exPos);
            assignFromKeyfile(keyFile, "RAW", "PrePreserv", pedited, raw.preser, pedited->raw.exPreser);

            if (ppVersion < 320) {
                assignFromKeyfile(keyFile, "RAW", "Method", pedited, raw.bayersensor.method, pedited->raw.bayersensor.method);
                assignFromKeyfile(keyFile, "RAW", "CcSteps", pedited, raw.bayersensor.ccSteps, pedited->raw.bayersensor.ccSteps);
                assignFromKeyfile(keyFile, "RAW", "LineDenoise", pedited, raw.bayersensor.linenoise, pedited->raw.bayersensor.linenoise);
                assignFromKeyfile(keyFile, "RAW", "GreenEqThreshold", pedited, raw.bayersensor.greenthresh, pedited->raw.bayersensor.greenEq);
                assignFromKeyfile(keyFile, "RAW", "DCBIterations", pedited, raw.bayersensor.dcb_iterations, pedited->raw.bayersensor.dcbIterations);
                assignFromKeyfile(keyFile, "RAW", "DCBEnhance", pedited, raw.bayersensor.dcb_enhance, pedited->raw.bayersensor.dcbEnhance);
                assignFromKeyfile(keyFile, "RAW", "LMMSEIterations", pedited, raw.bayersensor.lmmse_iterations, pedited->raw.bayersensor.lmmseIterations);
                assignFromKeyfile(keyFile, "RAW", "PreBlackzero", pedited, raw.bayersensor.black0, pedited->raw.bayersensor.exBlack0);
                assignFromKeyfile(keyFile, "RAW", "PreBlackone", pedited, raw.bayersensor.black1, pedited->raw.bayersensor.exBlack1);
                assignFromKeyfile(keyFile, "RAW", "PreBlacktwo", pedited, raw.bayersensor.black2, pedited->raw.bayersensor.exBlack2);
                assignFromKeyfile(keyFile, "RAW", "PreBlackthree", pedited, raw.bayersensor.black3, pedited->raw.bayersensor.exBlack3);
                assignFromKeyfile(keyFile, "RAW", "PreTwoGreen", pedited, raw.bayersensor.twogreen, pedited->raw.bayersensor.exTwoGreen);
            }
        }

        if (keyFile.has_group("RAW Bayer")) {
            assignFromKeyfile(keyFile, "RAW Bayer", "Method", pedited, raw.bayersensor.method, pedited->raw.bayersensor.method);

            if (keyFile.has_key("RAW Bayer", "ImageNum")) {
                raw.bayersensor.imageNum = keyFile.get_integer("RAW Bayer", "ImageNum") - 1;

                if (pedited) {
                    pedited->raw.bayersensor.imageNum = true;
                }
            }

            assignFromKeyfile(keyFile, "RAW Bayer", "CcSteps", pedited, raw.bayersensor.ccSteps, pedited->raw.bayersensor.ccSteps);
            assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack0", pedited, raw.bayersensor.black0, pedited->raw.bayersensor.exBlack0);
            assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack1", pedited, raw.bayersensor.black1, pedited->raw.bayersensor.exBlack1);
            assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack2", pedited, raw.bayersensor.black2, pedited->raw.bayersensor.exBlack2);
            assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack3", pedited, raw.bayersensor.black3, pedited->raw.bayersensor.exBlack3);
            assignFromKeyfile(keyFile, "RAW Bayer", "PreTwoGreen", pedited, raw.bayersensor.twogreen, pedited->raw.bayersensor.exTwoGreen);
            assignFromKeyfile(keyFile, "RAW Bayer", "LineDenoise", pedited, raw.bayersensor.linenoise, pedited->raw.bayersensor.linenoise);

            if (keyFile.has_key("RAW Bayer", "LineDenoiseDirection")) {
                raw.bayersensor.linenoiseDirection = RAWParams::BayerSensor::LineNoiseDirection(keyFile.get_integer("RAW Bayer", "LineDenoiseDirection"));

                if (pedited) {
                    pedited->raw.bayersensor.linenoiseDirection = true;
                }
            }

            assignFromKeyfile(keyFile, "RAW Bayer", "GreenEqThreshold", pedited, raw.bayersensor.greenthresh, pedited->raw.bayersensor.greenEq);
            assignFromKeyfile(keyFile, "RAW Bayer", "DCBIterations", pedited, raw.bayersensor.dcb_iterations, pedited->raw.bayersensor.dcbIterations);
            assignFromKeyfile(keyFile, "RAW Bayer", "DCBEnhance", pedited, raw.bayersensor.dcb_enhance, pedited->raw.bayersensor.dcbEnhance);
            assignFromKeyfile(keyFile, "RAW Bayer", "LMMSEIterations", pedited, raw.bayersensor.lmmse_iterations, pedited->raw.bayersensor.lmmseIterations);

            if (keyFile.has_key("RAW Bayer", "PixelShiftMotionCorrectionMethod")) {
                raw.bayersensor.pixelShiftMotionCorrectionMethod = (RAWParams::BayerSensor::PSMotionCorrectionMethod)keyFile.get_integer("RAW Bayer", "PixelShiftMotionCorrectionMethod");

                if (pedited) {
                    pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod = true;
                }
            }

            assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftEperIso", pedited, raw.bayersensor.pixelShiftEperIso, pedited->raw.bayersensor.pixelShiftEperIso);
            if (ppVersion < 332) {
                raw.bayersensor.pixelShiftEperIso += 1.0;
            }
            assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftSigma", pedited, raw.bayersensor.pixelShiftSigma, pedited->raw.bayersensor.pixelShiftSigma);
            assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftShowMotion", pedited, raw.bayersensor.pixelShiftShowMotion, pedited->raw.bayersensor.pixelShiftShowMotion);
            assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftShowMotionMaskOnly", pedited, raw.bayersensor.pixelShiftShowMotionMaskOnly, pedited->raw.bayersensor.pixelShiftShowMotionMaskOnly);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftHoleFill", pedited, raw.bayersensor.pixelShiftHoleFill, pedited->raw.bayersensor.pixelShiftHoleFill);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftMedian", pedited, raw.bayersensor.pixelShiftMedian, pedited->raw.bayersensor.pixelShiftMedian);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftGreen", pedited, raw.bayersensor.pixelShiftGreen, pedited->raw.bayersensor.pixelShiftGreen);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftBlur", pedited, raw.bayersensor.pixelShiftBlur, pedited->raw.bayersensor.pixelShiftBlur);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftSmoothFactor", pedited, raw.bayersensor.pixelShiftSmoothFactor, pedited->raw.bayersensor.pixelShiftSmooth);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftLmmse", pedited, raw.bayersensor.pixelShiftLmmse, pedited->raw.bayersensor.pixelShiftLmmse);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftEqualBright", pedited, raw.bayersensor.pixelShiftEqualBright, pedited->raw.bayersensor.pixelShiftEqualBright);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftEqualBrightChannel", pedited, raw.bayersensor.pixelShiftEqualBrightChannel, pedited->raw.bayersensor.pixelShiftEqualBrightChannel);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftNonGreenCross", pedited, raw.bayersensor.pixelShiftNonGreenCross, pedited->raw.bayersensor.pixelShiftNonGreenCross);
            assignFromKeyfile(keyFile, "RAW Bayer", "PDAFLinesFilter", pedited, raw.bayersensor.pdafLinesFilter, pedited->raw.bayersensor.pdafLinesFilter);
        }

        if (keyFile.has_group("RAW X-Trans")) {
            assignFromKeyfile(keyFile, "RAW X-Trans", "Method", pedited, raw.xtranssensor.method, pedited->raw.xtranssensor.method);
            assignFromKeyfile(keyFile, "RAW X-Trans", "CcSteps", pedited, raw.xtranssensor.ccSteps, pedited->raw.xtranssensor.ccSteps);
            assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackRed", pedited, raw.xtranssensor.blackred, pedited->raw.xtranssensor.exBlackRed);
            assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackGreen", pedited, raw.xtranssensor.blackgreen, pedited->raw.xtranssensor.exBlackGreen);
            assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackBlue", pedited, raw.xtranssensor.blackblue, pedited->raw.xtranssensor.exBlackBlue);
        }

        if (keyFile.has_group("MetaData")) {
            int mode = int(MetaDataParams::TUNNEL);
            assignFromKeyfile(keyFile, "MetaData", "Mode", pedited, mode, pedited->metadata.mode);

            if (mode >= int(MetaDataParams::TUNNEL) && mode <= int(MetaDataParams::STRIP)) {
                metadata.mode = static_cast<MetaDataParams::Mode>(mode);
            }
        }

        if (keyFile.has_group("Exif")) {
            std::vector<Glib::ustring> keys = keyFile.get_keys("Exif");

            for (const auto& key : keyFile.get_keys("Exif")) {
                exif[key] = keyFile.get_string("Exif", key);

                if (pedited) {
                    pedited->exif = true;
                }
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
        if (keyFile.has_group("IPTC")) {
            for (const auto& key : keyFile.get_keys("IPTC")) {
                // does this key already exist?
                const IPTCPairs::iterator element = iptc.find(key);

                if (element != iptc.end()) {
                    // it already exist so we cleanup the values
                    element->second.clear();
                }

                // TODO: look out if merging Keywords and SupplementalCategories from the procparams chain would be interesting
                for (const auto& currLoadedTagValue : keyFile.get_string_list("IPTC", key)) {
                    iptc[key].push_back(currLoadedTagValue);
                }

                if (pedited) {
                    pedited->iptc = true;
                }
            }
        }

        return 0;
    } catch (const Glib::Error& e) {
        printf("-->%s\n", e.what().c_str());
        setDefaults();
        return 1;
    } catch (...) {
        printf("-->unknown exception!\n");
        setDefaults();
        return 1;
    }

    return 0;
}

ProcParams* ProcParams::create()
{
    return new ProcParams();
}

void ProcParams::destroy(ProcParams* pp)
{
    delete pp;
}

bool ProcParams::operator ==(const ProcParams& other) const
{
    return
        toneCurve == other.toneCurve
        && retinex == other.retinex
        && localContrast == other.localContrast
        && labCurve == other.labCurve
        && sharpenEdge == other.sharpenEdge
        && sharpenMicro == other.sharpenMicro
        && sharpening == other.sharpening
        && prsharpening == other.prsharpening
        && vibrance == other.vibrance
        && wb == other.wb
        && colorappearance == other.colorappearance
        && impulseDenoise == other.impulseDenoise
        && dirpyrDenoise == other.dirpyrDenoise
        && epd == other.epd
        && fattal == other.fattal
        && defringe == other.defringe
        && sh == other.sh
        && crop == other.crop
        && coarse == other.coarse
        && rotate == other.rotate
        && commonTrans == other.commonTrans
        && distortion == other.distortion
        && lensProf == other.lensProf
        && perspective == other.perspective
        && gradient == other.gradient
        && pcvignette == other.pcvignette
        && cacorrection == other.cacorrection
        && vignetting == other.vignetting
        && chmixer == other.chmixer
        && blackwhite == other.blackwhite
        && resize == other.resize
        && raw == other.raw
        && icm == other.icm
        && wavelet == other.wavelet
        && dirpyrequalizer == other.dirpyrequalizer
        && hsvequalizer == other.hsvequalizer
        && filmSimulation == other.filmSimulation
        && rgbCurves == other.rgbCurves
        && colorToning == other.colorToning
        && metadata == other.metadata
        && exif == other.exif
        && iptc == other.iptc;
}

bool ProcParams::operator !=(const ProcParams& other) const
{
    return !(*this == other);
}

void ProcParams::init()
{
}

void ProcParams::cleanup()
{
}

int ProcParams::write(const Glib::ustring& fname, const Glib::ustring& content) const
{
    int error = 0;

    if (fname.length()) {
        FILE *f;
        f = g_fopen(fname.c_str(), "wt");

        if (f == nullptr) {
            error = 1;
        } else {
            fprintf(f, "%s", content.c_str());
            fclose(f);
        }
    }

    return error;
}

PartialProfile::PartialProfile(bool createInstance, bool paramsEditedValue)
{
    if (createInstance) {
        pparams = new ProcParams();
        pedited = new ParamsEdited(paramsEditedValue);
    } else {
        pparams = nullptr;
        pedited = nullptr;
    }
}

PartialProfile::PartialProfile(ProcParams* pp, ParamsEdited* pe, bool fullCopy)
{
    if (fullCopy && pp) {
        pparams = new ProcParams(*pp);
    } else {
        pparams = pp;
    }

    if (fullCopy && pe) {
        pedited = new ParamsEdited(*pe);
    } else {
        pedited = pe;
    }
}

PartialProfile::PartialProfile(const ProcParams* pp, const ParamsEdited* pe)
{
    if (pp) {
        pparams = new ProcParams(*pp);
    } else {
        pparams = nullptr;
    }

    if (pe) {
        pedited = new ParamsEdited(*pe);
    } else {
        pedited = nullptr;
    }
}

void PartialProfile::deleteInstance()
{
    if (pparams) {
        delete pparams;
        pparams = nullptr;
    }

    if (pedited) {
        delete pedited;
        pedited = nullptr;
    }
}

void PartialProfile::clearGeneral()
{
    if (pedited) {
        pedited->general.colorlabel = false;
        pedited->general.intrash = false;
        pedited->general.rank = false;
    }
}

int PartialProfile::load(const Glib::ustring& fName)
{
    if (!pparams) {
        pparams = new ProcParams();
    }

    if (!pedited) {
        pedited = new ParamsEdited();
    }

    if (fName == DEFPROFILE_INTERNAL) {
        return 0;
    } else if (fName == DEFPROFILE_DYNAMIC) {
        return -1; // should not happen here
    } else {
        return pparams->load(fName, pedited);
    }
}

/*
 * Set the all values of the General section to false
 * in order to preserve them in applyTo
 */
void PartialProfile::set(bool v)
{
    if (pedited) {
        pedited->set(v);
    }
}

void PartialProfile::applyTo(ProcParams* destParams) const
{
    if (destParams && pparams && pedited) {
        pedited->combine(*destParams, *pparams, true);
    }
}

AutoPartialProfile::AutoPartialProfile() :
    PartialProfile(true)
{
}

AutoPartialProfile::~AutoPartialProfile()
{
    deleteInstance();
}

}

}
