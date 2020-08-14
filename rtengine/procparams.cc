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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <map>

#include <locale.h>

#include <glib/gstdio.h>
#include <glibmm/fileutils.h>
#include <glibmm/miscutils.h>
#include <glibmm/keyfile.h>

#include "color.h"
#include "curves.h"
#include "procparams.h"
#include "utils.h"

#include "../rtgui/multilangmgr.h"
#include "../rtgui/options.h"
#include "../rtgui/paramsedited.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/version.h"

using namespace std;

namespace
{

Glib::ustring expandRelativePath(const Glib::ustring &procparams_fname, const Glib::ustring &prefix, Glib::ustring embedded_fname)
{
    if (embedded_fname.empty() || !Glib::path_is_absolute(procparams_fname)) {
        return embedded_fname;
    }

    if (!prefix.empty()) {
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
    if (fnameAbsolute || embedded_fname.empty() || !Glib::path_is_absolute(procparams_fname)) {
        return embedded_fname;
    }

    Glib::ustring prefix;

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
    std::vector<int>& value
)
{
    value = keyfile.get_integer_list(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    std::vector<double>& value
)
{
    value = keyfile.get_double_list(group_name, key);
    rtengine::sanitizeCurve(value);
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
    curveMode(ToneCurveMode::STD),
    curveMode2(ToneCurveMode::STD),
    brightness(0),
    black(0),
    contrast(0),
    saturation(0),
    shcompr(50),
    hlcompr(0),
    hlcomprthresh(0),
    histmatching(false),
    fromHistMatching(false),
    clampOOG(true)
{
}

bool ToneCurveParams::isPanningRelatedChange(const ToneCurveParams& other) const
{
    return !
        (autoexp == other.autoexp
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
        && clampOOG == other.clampOOG);
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
        && fromHistMatching == other.fromHistMatching
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
    complexmethod("normal"),
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
        && complexmethod == other.complexmethod
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
const double ColorToningParams::LABGRID_CORR_SCALE = 3.f;

ColorToningParams::LabCorrectionRegion::LabCorrectionRegion():
    a(0),
    b(0),
    saturation(0),
    slope(1),
    offset(0),
    power(1),
    hueMask{
        FCT_MinMaxCPoints,
        0.166666667,
        1.,
        0.35,
        0.35,
        0.8287775246,
        1.,
        0.35,
        0.35
    },
    chromaticityMask{
        FCT_MinMaxCPoints,
        0.,
        1.,
        0.35,
        0.35,
        1.,
        1.,
        0.35,
        0.35
    },
    lightnessMask{
        FCT_MinMaxCPoints,
        0.,
        1.,
        0.35,
        0.35,
        1.,
        1.,
        0.35,
        0.35
    },
    maskBlur(0),
    channel(ColorToningParams::LabCorrectionRegion::CHAN_ALL)
{
}

bool ColorToningParams::LabCorrectionRegion::operator==(const LabCorrectionRegion &other) const
{
    return
        a == other.a
        && b == other.b
        && saturation == other.saturation
        && slope == other.slope
        && offset == other.offset
        && power == other.power
        && hueMask == other.hueMask
        && chromaticityMask == other.chromaticityMask
        && lightnessMask == other.lightnessMask
        && maskBlur == other.maskBlur
        && channel == other.channel;
}

bool ColorToningParams::LabCorrectionRegion::operator!=(const LabCorrectionRegion &other) const
{
    return !(*this == other);
}

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
    method("LabRegions"),
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
    labgridBHigh(0.0),
    labregions{LabCorrectionRegion()},
    labregionsShowMask(-1)
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
        && labgridBHigh == other.labgridBHigh
        && labregions == other.labregions
        && labregionsShowMask == other.labregionsShowMask;
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
    contrast(20.0),
    autoContrast(false),
    blurradius(0.2),
    gamma(1.0),
    radius(0.5),
    amount(200),
    threshold(20, 80, 2000, 1200, false),
    edgesonly(false),
    edges_radius(1.9),
    edges_tolerance(1800),
    halocontrol(false),
    halocontrol_amount(85),
    method("usm"),
    deconvamount(100),
    deconvradius(0.75),
    deconviter(30),
    deconvdamping(0)
{
}

bool SharpeningParams::operator ==(const SharpeningParams& other) const
{
    return
        enabled == other.enabled
        && contrast == other.contrast
        && blurradius == other.blurradius
        && gamma == other.gamma
        && radius == other.radius
        && amount == other.amount
        && threshold == other.threshold
        && autoContrast == other.autoContrast
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

CaptureSharpeningParams::CaptureSharpeningParams() :
    enabled(false),
    autoContrast(true),
    autoRadius(true),
    contrast(10.0),
    deconvradius(0.75),
    deconvradiusOffset(0.0),
    deconviter(20),
    deconvitercheck(true)
{
}

bool CaptureSharpeningParams::operator ==(const CaptureSharpeningParams& other) const
{
    return
        enabled == other.enabled
        && contrast == other.contrast
        && autoContrast == other.autoContrast
        && autoRadius == other.autoRadius
        && deconvradius == other.deconvradius
        && deconvitercheck == other.deconvitercheck
        && deconvradiusOffset == other.deconvradiusOffset
        && deconviter == other.deconviter;
}

bool CaptureSharpeningParams::operator !=(const CaptureSharpeningParams& other) const
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
    contrast(20.0),
    uniformity(5)
{
}

bool SharpenMicroParams::operator ==(const SharpenMicroParams& other) const
{
    return
        enabled == other.enabled
        && matrix == other.matrix
        && amount == other.amount
        && contrast == other.contrast
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

bool WBParams::isPanningRelatedChange(const WBParams& other) const
{
    return
        !(
            enabled == other.enabled
            && (
                (
                    method == "Camera" 
                    && other.method == "Camera"
                )
            || (
                method == other.method
                && temperature == other.temperature
                && green == other.green
                && equal == other.equal
                && tempBias == other.tempBias
            )
        )
    );
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
        {"autitcgreen",       	 WBEntry::Type::AUTO,        M("TP_WBALANCE_AUTOITCGREEN"),   0, 1.f,    1.f,    0.f},
        {"autold",               WBEntry::Type::AUTO,        M("TP_WBALANCE_AUTOOLD"),        0, 1.f,   1.f,   0.f},
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
    illum("i50"),
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
    autotempout(true),
    ybout(18),
    greenout(1.0),
    tempsc(5003),
    greensc(1.0),
    presetcat02(false)
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
        && illum == other.illum
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
        && autotempout == other.autotempout
        && ybout == other.ybout
        && greenout == other.greenout
        && tempsc == other.tempsc
        && greensc == other.greensc
        && presetcat02 == other.presetcat02;
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
    threshold(30),
    amount(20),
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
    radius(40),
    lab(false)
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
        && radius == other.radius
        && lab == other.lab;
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
    method("log"),
    autofill(true)
{
}

bool CommonTransformParams::operator ==(const CommonTransformParams& other) const
{
    return method == other.method && autofill == other.autofill;
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
        && lfLens == other.lfLens
        && useDist == other.useDist
        && useVign == other.useVign;
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
    method("simple"),
    horizontal(0.0),
    vertical(0.0),
    camera_crop_factor(0.0),
    camera_focal_length(0.0),
    camera_pitch(0.0),
    camera_roll(0.0),
    camera_shift_horiz(0.0),
    camera_shift_vert(0.0),
    camera_yaw(0.0),
    projection_pitch(0.0),
    projection_rotate(0.0),
    projection_shift_horiz(0.0),
    projection_shift_vert(0.0),
    projection_yaw(0.0)
{
}

bool PerspectiveParams::operator ==(const PerspectiveParams& other) const
{
    return
        method == other.method
        && horizontal == other.horizontal
        && vertical == other.vertical
        && camera_focal_length == other.camera_focal_length
        && camera_crop_factor == other.camera_crop_factor
        && camera_pitch == other.camera_pitch
        && camera_roll == other.camera_roll
        && camera_shift_horiz == other.camera_shift_horiz
        && camera_shift_vert == other.camera_shift_vert
        && camera_yaw == other.camera_yaw
        && projection_shift_horiz == other.projection_shift_horiz
        && projection_shift_vert == other.projection_shift_vert
        && projection_rotate == other.projection_rotate
        && projection_pitch == other.projection_pitch
        && projection_yaw == other.projection_yaw;
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
        1000,
        0,
        0
    },
    green{
        0,
        1000,
        0
    },
    blue{
        0,
        0,
        1000
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
    setting("RGB-Rel"),
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
    height(900),
    allowUpscaling(false)
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
        && height == other.height
        && allowUpscaling == other.allowUpscaling;
}

bool ResizeParams::operator !=(const ResizeParams& other) const
{
    return !(*this == other);
}

const Glib::ustring ColorManagementParams::NoICMString = Glib::ustring("No ICM: sRGB output");

ColorManagementParams::ColorManagementParams() :
    inputProfile("(cameraICC)"),
    toneCurve(false),
    applyLookTable(false),
    applyBaselineExposureOffset(true),
    applyHueSatMap(true),
    dcpIlluminant(0),
    workingProfile("ProPhoto"),
    workingTRC("none"),
    workingTRCGamma(2.4),
    workingTRCSlope(12.92310),
    outputProfile(options.rtSettings.srgb),
    outputIntent(RI_RELATIVE),
    outputBPC(true)
{
}

bool ColorManagementParams::operator ==(const ColorManagementParams& other) const
{
    return
        inputProfile == other.inputProfile
        && toneCurve == other.toneCurve
        && applyLookTable == other.applyLookTable
        && applyBaselineExposureOffset == other.applyBaselineExposureOffset
        && applyHueSatMap == other.applyHueSatMap
        && dcpIlluminant == other.dcpIlluminant
        && workingProfile == other.workingProfile
        && workingTRC == other.workingTRC
        && workingTRCGamma == other.workingTRCGamma
        && workingTRCSlope == other.workingTRCSlope
        && outputProfile == other.outputProfile
        && outputIntent == other.outputIntent
        && outputBPC == other.outputBPC;
}

bool ColorManagementParams::operator !=(const ColorManagementParams& other) const
{
    return !(*this == other);
}

const double WaveletParams::LABGRID_CORR_MAX = 12800.f;
const double WaveletParams::LABGRID_CORR_SCALE = 3.276f;
const double WaveletParams::LABGRIDL_DIRECT_SCALE = 41950.;

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
    blcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0, 
        0.0, 
        0.0, 
        0.35, 
        0.5, 
        0., 
        0.35,
        0.35,
        1.0,
        0.0,
        0.35,
        0.35
/*      
        0.0,
        0.35, 
        0.35, 
        1.0, 
        0.0, 
        0.35, 
        0.35
*/
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
    opacityCurveSH{
        static_cast<double>(FCT_MinMaxCPoints),
        0.,
        1.,
        0.35,
        0.35,
        0.15,
        0.9,
        0.35,
        0.35,
        0.4,
        0.8,
        0.35,
        0.35,
        0.4,
        0.5,
        0.35,
        0.35,
        0.5,
        0.5,
        0.35,
        0.35,
        0.5,
        0.2,
        0.35,
        0.35,
        0.8,
        0.1,
        0.35,
        0.35,
        1.0,
        0.,
        0.35,
        0.35
    },
/*
    opacityCurveSH{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.,
        0.35,
        0.35,
        0.4,
        0.5,
        0.35,
        0.35,
        0.5,
        0.5,
        0.35,
        0.35,
        1.,
        0.,
        0.35,
        0.35
    },
*/
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
    linkedg(false),
    cbenab(false),
    greenlow(0),
    bluelow(0),
    greenmed(0),
    bluemed(0),
    greenhigh(0),
    bluehigh(0),
    ballum(7.),
    balchrom(0.),
    chromfi(0.),
    chromco(0.),
    mergeL(20.),
    mergeC(20.),
    softrad(0.),
    softradend(0.),
    lipst(false),
    avoid(false),
    showmask(false),
    oldsh(true),
    tmr(false),
    strength(100),
    balance(0),
    sigmafin(1.0),
    sigmaton(1.0),
    sigmacol(1.0),
    sigmadir(1.0),
    rangeab(20.0),
    protab(0.0),
    iter(0),
    expcontrast(false),
    expchroma(false),
    c{},
    ch{},
    expedge(false),
    expbl(false),
    expresid(false),
    expfinal(false),
    exptoning(false),
    expnoise(false),
    expclari(false),
    labgridALow(0.0),
    labgridBLow(0.0),
    labgridAHigh(0.0),
    labgridBHigh(0.0),
    Lmethod(4),
    CLmethod("all"),
    Backmethod("grey"),
    Tilesmethod("full"),
    complexmethod("normal"),
    daubcoeffmethod("4_"),
    CHmethod("without"),
    Medgreinf("less"),
    ushamethod("clari"),
    CHSLmethod("SL"),
    EDmethod("CU"),
    NPmethod("none"),
    BAmethod("none"),
    TMmethod("cont"),
    Dirmethod("all"),
    HSmethod("with"),
    sigma(1.0),
    offset(1.0),
    lowthr(40.0),
    rescon(0),
    resconH(0),
    reschro(0),
    resblur(0),
    resblurc(0),
    tmrs(0),
    edgs(1.4),
    scale(1.),
    gamma(1),
    sup(0),
    sky(0.0),
    thres(7),
    chroma(5),
    chro(0),
    threshold(4),
    threshold2(5),
    edgedetect(90),
    edgedetectthr(20),
    edgedetectthr2(0),
    edgesensi(60),
    edgeampli(10),
    contrast(0),
    edgrad(15),
    edgeffect(1.0),
    edgval(0),
    edgthresh(10),
    thr(30),
    thrH(70),
    radius(40),
    skinprotect(0.0),
    chrwav(0.),
    bluwav(1.0),
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
        && blcurve == other.blcurve
        && opacityCurveRG == other.opacityCurveRG
        && opacityCurveSH == other.opacityCurveSH
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
        && ballum == other.ballum
        && balchrom == other.balchrom
        && chromfi == other.chromfi
        && chromco == other.chromco
        && mergeL == other.mergeL
        && mergeC == other.mergeC
        && softrad == other.softrad
        && softradend == other.softradend
        && lipst == other.lipst
        && avoid == other.avoid
        && showmask == other.showmask
        && oldsh == other.oldsh
        && tmr == other.tmr
        && strength == other.strength
        && balance == other.balance
        && sigmafin == other.sigmafin
        && sigmaton == other.sigmaton
        && sigmacol == other.sigmacol
        && sigmadir == other.sigmadir
        && rangeab == other.rangeab
        && protab == other.protab
        && iter == other.iter
        && labgridALow == other.labgridALow
        && labgridBLow == other.labgridBLow
        && labgridAHigh == other.labgridAHigh
        && labgridBHigh == other.labgridBHigh
        && expcontrast == other.expcontrast
        && expchroma == other.expchroma
        && [this, &other]() -> bool
            {
                for (unsigned int i = 0; i < 9; ++i) {
                    if (c[i] != other.c[i] || ch[i] != other.ch[i]) {
                        return false;
                    }
                }
                return true;
            }()
        && expedge == other.expedge
        && expbl == other.expbl
        && expresid == other.expresid
        && expfinal == other.expfinal
        && expclari == other.expclari
        && exptoning == other.exptoning
        && expnoise == other.expnoise
        && Lmethod == other.Lmethod
        && CLmethod == other.CLmethod
        && Backmethod == other.Backmethod
        && Tilesmethod == other.Tilesmethod
        && complexmethod == other.complexmethod
        && daubcoeffmethod == other.daubcoeffmethod
        && CHmethod == other.CHmethod
        && Medgreinf == other.Medgreinf
        && ushamethod == other.ushamethod
        && CHSLmethod == other.CHSLmethod
        && EDmethod == other.EDmethod
        && NPmethod == other.NPmethod
        && BAmethod == other.BAmethod
        && TMmethod == other.TMmethod
        && Dirmethod == other.Dirmethod
        && HSmethod == other.HSmethod
        && sigma == other.sigma
        && offset == other.offset
        && lowthr == other.lowthr
        && rescon == other.rescon
        && resconH == other.resconH
        && reschro == other.reschro
        && resblur == other.resblur
        && resblurc == other.resblurc
        && tmrs == other.tmrs
        && edgs == other.edgs
        && scale == other.scale
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
        && edgeffect == other.edgeffect
        && edgval == other.edgval
        && edgthresh == other.edgthresh
        && thr == other.thr
        && thrH == other.thrH
        && radius == other.radius
        && skinprotect == other.skinprotect
        && chrwav == other.chrwav
        && bluwav == other.bluwav
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
    Wavblcurve& tCurve,
    WavOpacityCurveRG& opacityCurveLUTRG,
    WavOpacityCurveSH& opacityCurveLUTSH,
    WavOpacityCurveBY& opacityCurveLUTBY,
    WavOpacityCurveW& opacityCurveLUTW,
    WavOpacityCurveWL& opacityCurveLUTWL
) const
{
    cCurve.Set(this->ccwcurve);
    tCurve.Set(this->blcurve);
    opacityCurveLUTRG.Set(this->opacityCurveRG);
    opacityCurveLUTSH.Set(this->opacityCurveSH);
    opacityCurveLUTBY.Set(this->opacityCurveBY);
    opacityCurveLUTW.Set(this->opacityCurveW);
    opacityCurveLUTWL.Set(this->opacityCurveWL);

}

LocallabParams::LocallabSpot::LocallabSpot() :
    // Control spot settings
    name(""),
    isvisible(true),
    prevMethod("hide"),
    shape("ELI"),
    spotMethod("norm"),
    wavMethod("D4"),
    sensiexclu(12),
    structexclu(0),
    struc(4.0),
    shapeMethod("IND"),
    loc{150, 150, 150, 150},
    centerX(0),
    centerY(0),
    circrad(18),
    qualityMethod("enh"),
    complexMethod("mod"),
    transit(60.),
    feather(25.),
    thresh(2.0),
    iter(2.0),
    balan(1.0),
    balanh(1.0),
    colorde(5.0),
    colorscope(30.0),
    transitweak(1.0),
    transitgrad(0.0),
    activ(true),
    avoid(false),
    blwh(false),
    recurs(false),
    laplac(true),
    deltae(true),
    shortc(false),
    savrest(false),
    scopemask(60),
    lumask(10),
    // Color & Light
    visicolor(false),
    expcolor(false),
    complexcolor(2),
    curvactiv(false),
    lightness(0),
    contrast(0),
    chroma(0),
    labgridALow(0.0),
    labgridBLow(0.0),
    labgridAHigh(0.0),
    labgridBHigh(0.0),
    labgridALowmerg(0.0),
    labgridBLowmerg(0.0),
    labgridAHighmerg(-3500.0),
    labgridBHighmerg(-4600.0),
    strengthgrid(30),
    sensi(15),
    structcol(0),
    strcol(0.),
    strcolab(0.),
    strcolh(0.),
    angcol(0.),
    blurcolde(5),
    blurcol(0.2),
    contcol(0.),
    blendmaskcol(0),
    radmaskcol(0.0),
    chromaskcol(0.0),
    gammaskcol(1.0),
    slomaskcol(0.0),
    shadmaskcol(0),
    strumaskcol(0.),
    lapmaskcol(0.0),
    qualitycurveMethod("none"),
    gridMethod("one"),
    merMethod("mone"),
    toneMethod("fou"),
    mergecolMethod("one"),
    llcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    lccurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    cccurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    clcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    rgbcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    LHcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.50,
        0.35,
        0.35,
        0.166,
        0.50,
        0.35,
        0.35,
        0.333,
        0.50,
        0.35,
        0.35,
        0.50,
        0.50,
        0.35,
        0.35,
        0.666,
        0.50,
        0.35,
        0.35,
        0.833,
        0.50,
        0.35,
        0.35
    },
    HHcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.50,
        0.35,
        0.35,
        0.166,
        0.50,
        0.35,
        0.35,
        0.333,
        0.50,
        0.35,
        0.35,
        0.50,
        0.50,
        0.35,
        0.35,
        0.666,
        0.50,
        0.35,
        0.35,
        0.833,
        0.50,
        0.35,
        0.35
    },
    CHcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.50,
        0.35,
        0.35,
        0.166,
        0.50,
        0.35,
        0.35,
        0.333,
        0.50,
        0.35,
        0.35,
        0.50,
        0.50,
        0.35,
        0.35,
        0.666,
        0.50,
        0.35,
        0.35,
        0.833,
        0.50,
        0.35,
        0.35
    },
    invers(false),
    special(false),
    toolcol(true),
    enaColorMask(false),
    fftColorMask(true),
    CCmaskcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.00,
        1.0,
        0.35,
        0.35
    },
    LLmaskcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.00,
        1.0,
        0.35,
        0.35
    },
    HHmaskcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.00,
        1.0,
        0.35,
        0.35
    },
    HHhmaskcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.5,
        0.35,
        0.35,
        0.50,
        0.5,
        0.35,
        0.35,
        1.00,
        0.5,
        0.35,
        0.35
    },
    softradiuscol(0.0),
    opacol(60.0),
    mercol(18.0),
    merlucol(32.0),
    conthrcol(0.0),
    Lmaskcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    LLmaskcolcurvewav{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.5,
        0.35,
        0.35,
        1.,
        0.5,
        0.35,
        0.35
    },
    csthresholdcol(0, 0, 6, 5, false),
    // Exposure
    visiexpose(false),
    expexpose(false),
    complexexpose(0),
    expcomp(0.0),
    hlcompr(20),
    hlcomprthresh(0),
    black(0),
    shadex(0),
    shcompr(50),
    expchroma(5),
    sensiex(60),
    structexp(0),
    blurexpde(5),
    strexp(0.),
    angexp(0.),
    excurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    inversex(false),
    enaExpMask(false),
    enaExpMaskaft(false),
    CCmaskexpcurve{
        static_cast<double>(FCT_MinMaxCPoints),0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmaskexpcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmaskexpcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    blendmaskexp(0),
    radmaskexp(0.0),
    chromaskexp(0.0),
    gammaskexp(1.0),
    slomaskexp(0.0),
    lapmaskexp(0.0),
    strmaskexp(0.0),
    angmaskexp(0.0),
    softradiusexp(0.0),
    Lmaskexpcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    expMethod("std"),
    exnoiseMethod("none"),
    laplacexp(0.0),
    balanexp(1.0),
    linear(0.05),
    gamm(0.4),
    fatamount(1.0),
    fatdetail(40.0),
    fatanchor(1.0),
    fatlevel(1.),
    // Shadow highlight
    visishadhigh(false),
    expshadhigh(false),
    complexshadhigh(0),
    shMethod("tone"),
    multsh{0, 0, 0, 0, 0},
    highlights(0),
    h_tonalwidth(70),
    shadows(0),
    s_tonalwidth(30),
    sh_radius(40),
    sensihs(15),
    enaSHMask(false),
    CCmaskSHcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmaskSHcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmaskSHcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    blendmaskSH(0),
    radmaskSH(0.0),
    blurSHde(5),
    strSH(0.),
    angSH(0.),
    inverssh(false),
    chromaskSH(0.0),
    gammaskSH(1.0),
    slomaskSH(0.0),
    lapmaskSH(0.0),
    detailSH(0),
    LmaskSHcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    fatamountSH(1.0),
    fatanchorSH(50.0),
    gamSH(2.4),
    sloSH(12.92),
    // Vibrance
    visivibrance(false),
    expvibrance(false),
    complexvibrance(0),
    saturated(0),
    pastels(0),
    warm(0),
    psthreshold({0, 75, false}),
    protectskins(false),
    avoidcolorshift(true),
    pastsattog(true),
    sensiv(15),
    skintonescurve{
        static_cast<double>(DCT_Linear)
    },
    CCmaskvibcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmaskvibcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmaskvibcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    enavibMask(false),
    blendmaskvib(0),
    radmaskvib(0.0),
    chromaskvib(0.0),
    gammaskvib(1.0),
    slomaskvib(0.0),
    lapmaskvib(0.0),
    strvib(0.0),
    strvibab(0.0),
    strvibh(0.0),
    angvib(0.0),
    Lmaskvibcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    // Soft Light
    visisoft(false),
    expsoft(false),
    complexsoft(0),
    streng(0),
    sensisf(30),
    laplace(25.),
    softMethod("soft"),
    // Blur & Noise
    visiblur(false),
    expblur(false),
    complexblur(0),
    radius(1.5),
    strength(0),
    sensibn(40),
    itera(1),
    guidbl(0),
    strbl(50),
    isogr(400),
    strengr(0),
    scalegr(100),
    epsbl(0),
    blMethod("blur"),
    chroMethod("lum"),
    blurMethod("norm"),
    medMethod("33"),
    activlum(true),
    noiselumf(0.),
    noiselumf0(0.),
    noiselumf2(0.),
    noiselumc(0.),
    noiselumdetail(0.),
    noiselequal(7),
    noisechrof(0.),
    noisechroc(0.),
    noisechrodetail(0.),
    adjblur(0),
    bilateral(0),
    sensiden(60),
    detailthr(0),
    locwavcurveden{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.0,
        0.0,
        0.35,
        0.5,
        0.,
        0.35,
        0.35,
        1.0,
        0.0,
        0.35,
        0.35
    },
    showmaskblMethodtyp("blur"),
    CCmaskblcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmaskblcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmaskblcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    enablMask(false),
    fftwbl(false),
    toolbl(false),
    blendmaskbl(0),
    radmaskbl(0.0),
    chromaskbl(0.0),
    gammaskbl(1.0),
    slomaskbl(0.0),
    lapmaskbl(0.0),
    shadmaskbl(0),
    shadmaskblsha(0),
    strumaskbl(0.),
    Lmaskblcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    LLmaskblcurvewav{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.5,
        0.35,
        0.35,
        1.,
        0.5,
        0.35,
        0.35
    },
    csthresholdblur(0, 0, 6, 5, false),
    // Tone Mapping
    visitonemap(false),
    exptonemap(false),
    complextonemap(0),
    stren(0.5),
    gamma(1.0),
    estop(1.4),
    scaltm(1.0),
    rewei(0),
    satur(0.),
    sensitm(60),
    softradiustm(0.0),
    amount(95.),
    equiltm(true),
    CCmasktmcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmasktmcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmasktmcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    enatmMask(false),
    enatmMaskaft(false),
    blendmasktm(0),
    radmasktm(0.0),
    chromasktm(0.0),
    gammasktm(1.0),
    slomasktm(0.0),
    lapmasktm(0.0),
    Lmasktmcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    // Retinex
    visireti(false),
    expreti(false),
    complexreti(0),
    retinexMethod("high"),
    str(0.),
    chrrt(0.0),
    neigh(50.0),
    vart(150.0),
    offs(0.0),
    dehaz(0),
    depth(25),
    sensih(60),
    localTgaincurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.12,
        0.35,
        0.35,
        0.70,
        0.50,
        0.35,
        0.35,
        1.00,
        0.12,
        0.35,
        0.35
    },
    localTtranscurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.50,
        0.35,
        0.35,
        0.5,
        0.5,
        0.35,
        0.35,
        1.00,
        0.50,
        0.35,
        0.35
    },
    inversret(false),
    equilret(true),
    loglin(false),
    lumonly(false),
    softradiusret(40.0),
    CCmaskreticurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmaskreticurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmaskreticurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    enaretiMask(false),
    enaretiMasktmap(true),
    blendmaskreti(0),
    radmaskreti(0.0),
    chromaskreti(0.0),
    gammaskreti(1.0),
    slomaskreti(0.0),
    lapmaskreti(0.0),
    scalereti(2.0),
    darkness(2.0),
    lightnessreti(1.0),
    limd(8.0),
    cliptm(1.0),
    fftwreti(false),
    Lmaskreticurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    // Sharpening
    visisharp(false),
    expsharp(false),
    complexsharp(0),
    sharcontrast(20),
    sharradius(0.75),
    sharamount(100),
    shardamping(0),
    shariter(30),
    sharblur(0.2),
    sensisha(40),
    inverssha(false),
    // Local Contrast
    visicontrast(false),
    expcontrast(false),
    complexcontrast(0),
    lcradius(80),
    lcamount(0.0),
    lcdarkness(1.0),
    lclightness(1.0),
    sigmalc(1.0),
    levelwav(4),
    residcont(0.0),
    residsha(0.0),
    residshathr(30.0),
    residhi(0.0),
    residhithr(70.0),
    residblur(0.0),
    levelblur(0.0),
    sigmabl(1.0),
    residchro(0.0),
    residcomp(0.0),
    sigma(1.0),
    offset(1.0),
    sigmadr(1.0),
    threswav(1.4),
    chromalev(1.0),
    chromablu(0.0),
    sigmadc(1.0),
    deltad(0.0),
    fatres(0.0),
    clarilres(0.0),
    claricres(0.0),
    clarisoft(1.0),
    sigmalc2(1.0),
    strwav(0.0),
    angwav(0.0),
    strengthw(0.0),
    sigmaed(1.0),
    radiusw(15.0),
    detailw(10.0),
    gradw(90.0),
    tloww(20.0),
    thigw(0.0),
    edgw(60.0),
    basew(10.0),
    sensilc(60),
    fftwlc(false),
    blurlc(true),
    wavblur(false),
    wavedg(false),
    waveshow(false),
    wavcont(false),
    wavcomp(false),
    wavgradl(false),
    wavcompre(false),
    origlc(false),
    localcontMethod("loc"),
    localedgMethod("thr"),
    localneiMethod("low"),
    locwavcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.5,
        0.35,
        0.35,
        1.,
        0.5,
        0.35,
        0.35
    },
    csthreshold(0, 0, 6, 6, false),
    loclevwavcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.0,
        0.0,
        0.35,
        0.5,
        0.,
        0.35,
        0.35,
        1.0,
        0.0,
        0.35,
        0.35
    },
    locconwavcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.5,
        0.35,
        0.35,
        1.,
        0.5,
        0.35,
        0.35
    },
    loccompwavcurve{
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
    loccomprewavcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.75,
        0.35,
        0.35,
        1.,
        0.75,
        0.35,
        0.35
    },
    locedgwavcurve{
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
    CCmasklccurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmasklccurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmasklccurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    enalcMask(false),
    blendmasklc(0),
    radmasklc(0.0),
    chromasklc(0.0),
    Lmasklccurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    // Contrast by detail levels
    visicbdl(false),
    expcbdl(false),
    complexcbdl(0),
    mult{1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
    chromacbdl(0.),
    threshold(0.2),
    sensicb(60),
    clarityml(0.1),
    contresid(0),
    softradiuscb(0.0),
    enacbMask(false),
    CCmaskcbcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    LLmaskcbcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    HHmaskcbcurve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.0,
        1.0,
        0.35,
        0.35
    },
    blendmaskcb(0),
    radmaskcb(0.0),
    chromaskcb(0.0),
    gammaskcb(1.0),
    slomaskcb(0.0),
    lapmaskcb(0.0),
    Lmaskcbcurve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    // Log encoding
    visilog(false),
    explog(false),
    autocompute(false),
    sourceGray(10.),
    targetGray(18.),
    Autogray(true),
    fullimage(true),
    blackEv(-5.0),
    whiteEv(10.0),
    detail(0.6),
    sensilog(60),
    baselog(2.),
    strlog(0.0),
    anglog(0.0),
    // mask
    visimask(false),
    complexmask(0),
    expmask(false),
    sensimask(60),
    blendmask(-10.),
    blendmaskab(-10.),
    softradiusmask(1.0),
    enamask(false),
    fftmask(true),
    blurmask(0.2),
    contmask(0.),
    CCmask_curve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.00,
        1.0,
        0.35,
        0.35
    },
    LLmask_curve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.00,
        1.0,
        0.35,
        0.35
    },
    HHmask_curve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        1.0,
        0.35,
        0.35,
        0.50,
        1.0,
        0.35,
        0.35,
        1.00,
        1.0,
        0.35,
        0.35
    },
    strumaskmask(0.),
    toolmask(true),
    radmask(0.0),
    lapmask(0.0),
    chromask(0.0),
    gammask(1.0),
    slopmask(0.0),
    shadmask(0.0),
    str_mask(0),
    ang_mask(0),
    HHhmask_curve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.5,
        0.35,
        0.35,
        0.50,
        0.5,
        0.35,
        0.35,
        1.00,
        0.5,
        0.35,
        0.35
    },
    Lmask_curve{
        static_cast<double>(DCT_NURBS),
        0.0,
        0.0,
        1.0,
        1.0
    },
    LLmask_curvewav{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.5,
        0.35,
        0.35,
        1.,
        0.5,
        0.35,
        0.35
    },
    csthresholdmask(0, 0, 6, 5, false)

{
}

bool LocallabParams::LocallabSpot::operator ==(const LocallabSpot& other) const
{
    return
        // Control spot settings
        name == other.name
        && isvisible == other.isvisible
        && prevMethod == other.prevMethod
        && shape == other.shape
        && spotMethod == other.spotMethod
        && wavMethod == other.wavMethod
        && sensiexclu == other.sensiexclu
        && structexclu == other.structexclu
        && struc == other.struc
        && shapeMethod == other.shapeMethod
        && loc == other.loc
        && centerX == other.centerX
        && centerY == other.centerY
        && circrad == other.circrad
        && qualityMethod == other.qualityMethod
        && complexMethod == other.complexMethod
        && transit == other.transit
        && feather == other.feather
        && thresh == other.thresh
        && iter == other.iter
        && balan == other.balan
        && balanh == other.balanh
        && colorde == other.colorde
        && colorscope == other.colorscope
        && transitweak == other.transitweak
        && transitgrad == other.transitgrad
        && activ == other.activ
        && avoid == other.avoid
        && blwh == other.blwh
        && recurs == other.recurs
        && laplac == other.laplac
        && deltae == other.deltae
        && shortc == other.shortc
        && savrest == other.savrest
        && scopemask == other.scopemask
        && lumask == other.lumask
        // Color & Light
        && visicolor == other.visicolor
        && expcolor == other.expcolor
        && complexcolor == other.complexcolor
        && curvactiv == other.curvactiv
        && lightness == other.lightness
        && contrast == other.contrast
        && chroma == other.chroma
        && labgridALow == other.labgridALow
        && labgridBLow == other.labgridBLow
        && labgridAHigh == other.labgridAHigh
        && labgridBHigh == other.labgridBHigh
        && labgridALowmerg == other.labgridALowmerg
        && labgridBLowmerg == other.labgridBLowmerg
        && labgridAHighmerg == other.labgridAHighmerg
        && labgridBHighmerg == other.labgridBHighmerg
        && strengthgrid == other.strengthgrid
        && sensi == other.sensi
        && structcol == other.structcol
        && strcol == other.strcol
        && strcolab == other.strcolab
        && strcolh == other.strcolh
        && angcol == other.angcol
        && blurcolde == other.blurcolde
        && blurcol == other.blurcol
        && contcol == other.contcol
        && blendmaskcol == other.blendmaskcol
        && radmaskcol == other.radmaskcol
        && chromaskcol == other.chromaskcol
        && gammaskcol == other.gammaskcol
        && slomaskcol == other.slomaskcol
        && shadmaskcol == other.shadmaskcol
        && strumaskcol == other.strumaskcol
        && lapmaskcol == other.lapmaskcol
        && qualitycurveMethod == other.qualitycurveMethod
        && gridMethod == other.gridMethod
        && merMethod == other.merMethod
        && toneMethod == other.toneMethod
        && mergecolMethod == other.mergecolMethod
        && llcurve == other.llcurve
        && lccurve == other.lccurve
        && cccurve == other.cccurve
        && clcurve == other.clcurve
        && rgbcurve == other.rgbcurve
        && LHcurve == other.LHcurve
        && HHcurve == other.HHcurve
        && CHcurve == other.CHcurve
        && invers == other.invers
        && special == other.special
        && toolcol == other.toolcol
        && enaColorMask == other.enaColorMask
        && fftColorMask == other.fftColorMask
        && CCmaskcurve == other.CCmaskcurve
        && LLmaskcurve == other.LLmaskcurve
        && HHmaskcurve == other.HHmaskcurve
        && HHhmaskcurve == other.HHhmaskcurve
        && softradiuscol == other.softradiuscol
        && opacol == other.opacol
        && mercol == other.mercol
        && merlucol == other.merlucol
        && conthrcol == other.conthrcol
        && Lmaskcurve == other.Lmaskcurve
        && LLmaskcolcurvewav == other.LLmaskcolcurvewav
        && csthresholdcol == other.csthresholdcol
        // Exposure
        && visiexpose == other.visiexpose
        && expexpose == other.expexpose
        && complexexpose == other.complexexpose
        && expcomp == other.expcomp
        && hlcompr == other.hlcompr
        && hlcomprthresh == other.hlcomprthresh
        && black == other.black
        && shadex == other.shadex
        && shcompr == other.shcompr
        && expchroma == other.expchroma
        && sensiex == other.sensiex
        && structexp == other.structexp
        && blurexpde == other.blurexpde
        && strexp == other.strexp
        && angexp == other.angexp
        && excurve == other.excurve
        && inversex == other.inversex
        && enaExpMask == other.enaExpMask
        && enaExpMaskaft == other.enaExpMaskaft
        && CCmaskexpcurve == other.CCmaskexpcurve
        && LLmaskexpcurve == other.LLmaskexpcurve
        && HHmaskexpcurve == other.HHmaskexpcurve
        && blendmaskexp == other.blendmaskexp
        && radmaskexp == other.radmaskexp
        && chromaskexp == other.chromaskexp
        && gammaskexp == other.gammaskexp
        && slomaskexp == other.slomaskexp
        && lapmaskexp == other.lapmaskexp
        && strmaskexp == other.strmaskexp
        && angmaskexp == other.angmaskexp
        && softradiusexp == other.softradiusexp
        && Lmaskexpcurve == other.Lmaskexpcurve
        && expMethod == other.expMethod
        && exnoiseMethod == other.exnoiseMethod
        && laplacexp == other.laplacexp
        && balanexp == other.balanexp
        && linear == other.linear
        && gamm == other.gamm
        && fatamount == other.fatamount
        && fatdetail == other.fatdetail
        && fatanchor == other.fatanchor
        && fatlevel == other.fatlevel
        // Shadow highlight
        && visishadhigh == other.visishadhigh
        && expshadhigh == other.expshadhigh
        && complexshadhigh == other.complexshadhigh
        && shMethod == other.shMethod
        && [this, &other]() -> bool
            {
                for (int i = 0; i < 5; ++i) {
                    if (multsh[i] != other.multsh[i]) {
                        return false;
                    }
                }
                return true;
            }()
        && highlights == other.highlights
        && h_tonalwidth == other.h_tonalwidth
        && shadows == other.shadows
        && s_tonalwidth == other.s_tonalwidth
        && sh_radius == other.sh_radius
        && sensihs == other.sensihs
        && enaSHMask == other.enaSHMask
        && CCmaskSHcurve == other.CCmaskSHcurve
        && LLmaskSHcurve == other.LLmaskSHcurve
        && HHmaskSHcurve == other.HHmaskSHcurve
        && blendmaskSH == other.blendmaskSH
        && radmaskSH == other.radmaskSH
        && blurSHde == other.blurSHde
        && strSH == other.strSH
        && angSH == other.angSH
        && inverssh == other.inverssh
        && chromaskSH == other.chromaskSH
        && gammaskSH == other.gammaskSH
        && slomaskSH == other.slomaskSH
        && lapmaskSH == other.lapmaskSH
        && detailSH == other.detailSH
        && LmaskSHcurve == other.LmaskSHcurve
        && fatamountSH == other.fatamountSH
        && fatanchorSH == other.fatanchorSH
        && gamSH == other.gamSH
        && sloSH == other.sloSH
        // Vibrance
        && visivibrance == other.visivibrance
        && expvibrance == other.expvibrance
        && complexvibrance == other.complexvibrance
        && saturated == other.saturated
        && pastels == other.pastels
        && warm == other.warm
        && psthreshold == other.psthreshold
        && protectskins == other.protectskins
        && avoidcolorshift == other.avoidcolorshift
        && pastsattog == other.pastsattog
        && sensiv == other.sensiv
        && skintonescurve == other.skintonescurve
        && CCmaskvibcurve == other.CCmaskvibcurve
        && LLmaskvibcurve == other.LLmaskvibcurve
        && HHmaskvibcurve == other.HHmaskvibcurve
        && enavibMask == other.enavibMask
        && blendmaskvib == other.blendmaskvib
        && radmaskvib == other.radmaskvib
        && chromaskvib == other.chromaskvib
        && gammaskvib == other.gammaskvib
        && slomaskvib == other.slomaskvib
        && lapmaskvib == other.lapmaskvib
        && strvib == other.strvib
        && strvibab == other.strvibab
        && strvibh == other.strvibh
        && angvib == other.angvib
        && Lmaskvibcurve == other.Lmaskvibcurve
        // Soft Light
        && visisoft == other.visisoft
        && expsoft == other.expsoft
        && complexsoft == other.complexsoft
        && streng == other.streng
        && sensisf == other.sensisf
        && laplace == other.laplace
        && softMethod == other.softMethod
        // Blur & Noise
        && visiblur == other.visiblur
        && expblur == other.expblur
        && complexblur == other.complexblur
        && radius == other.radius
        && strength == other.strength
        && sensibn == other.sensibn
        && itera == other.itera
        && guidbl == other.guidbl
        && strbl == other.strbl
        && isogr == other.isogr
        && strengr == other.strengr
        && scalegr == other.scalegr
        && epsbl == other.epsbl
        && blMethod == other.blMethod
        && chroMethod == other.chroMethod
        && blurMethod == other.blurMethod
        && medMethod == other.medMethod
        && activlum == other.activlum
        && noiselumf == other.noiselumf
        && noiselumf0 == other.noiselumf0
        && noiselumf2 == other.noiselumf2
        && noiselumc == other.noiselumc
        && noiselumdetail == other.noiselumdetail
        && noiselequal == other.noiselequal
        && noisechrof == other.noisechrof
        && noisechroc == other.noisechroc
        && noisechrodetail == other.noisechrodetail
        && adjblur == other.adjblur
        && bilateral == other.bilateral
        && sensiden == other.sensiden
        && detailthr == other.detailthr
        && locwavcurveden == other.locwavcurveden
        && showmaskblMethodtyp == other.showmaskblMethodtyp
        && CCmaskblcurve == other.CCmaskblcurve
        && LLmaskblcurve == other.LLmaskblcurve
        && HHmaskblcurve == other.HHmaskblcurve
        && enablMask == other.enablMask
        && fftwbl == other.fftwbl
        && toolbl == other.toolbl
        && blendmaskbl == other.blendmaskbl
        && radmaskbl == other.radmaskbl
        && chromaskbl == other.chromaskbl
        && gammaskbl == other.gammaskbl
        && slomaskbl == other.slomaskbl
        && lapmaskbl == other.lapmaskbl
        && shadmaskbl == other.shadmaskbl
        && shadmaskblsha == other.shadmaskblsha
        && strumaskbl == other.strumaskbl
        && Lmaskblcurve == other.Lmaskblcurve
        && LLmaskblcurvewav == other.LLmaskblcurvewav
        && csthresholdblur == other.csthresholdblur
        // Tone Mapping
        && visitonemap == other.visitonemap
        && exptonemap == other.exptonemap
        && complextonemap == other.complextonemap
        && stren == other.stren
        && gamma == other.gamma
        && estop == other.estop
        && scaltm == other.scaltm
        && rewei == other.rewei
        && satur == other.satur
        && sensitm == other.sensitm
        && softradiustm == other.softradiustm
        && amount == other.amount
        && equiltm == other.equiltm
        && CCmasktmcurve == other.CCmasktmcurve
        && LLmasktmcurve == other.LLmasktmcurve
        && HHmasktmcurve == other.HHmasktmcurve
        && enatmMask == other.enatmMask
        && enatmMaskaft == other.enatmMaskaft
        && blendmasktm == other.blendmasktm
        && radmasktm == other.radmasktm
        && chromasktm == other.chromasktm
        && gammasktm == other.gammasktm
        && slomasktm == other.slomasktm
        && lapmasktm == other.lapmasktm
        && Lmasktmcurve == other.Lmasktmcurve
        // Retinex
        && visireti == other.visireti
        && expreti == other.expreti
        && complexreti == other.complexreti
        && retinexMethod == other.retinexMethod
        && str == other.str
        && chrrt == other.chrrt
        && neigh == other.neigh
        && vart == other.vart
        && offs == other.offs
        && dehaz == other.dehaz
        && depth == other.depth
        && sensih == other.sensih
        && localTgaincurve == other.localTgaincurve
        && localTtranscurve == other.localTtranscurve
        && inversret == other.inversret
        && equilret == other.equilret
        && loglin == other.loglin
        && lumonly == other.lumonly
        && softradiusret == other.softradiusret
        && CCmaskreticurve == other.CCmaskreticurve
        && LLmaskreticurve == other.LLmaskreticurve
        && HHmaskreticurve == other.HHmaskreticurve
        && enaretiMask == other.enaretiMask
        && enaretiMasktmap == other.enaretiMasktmap
        && blendmaskreti == other.blendmaskreti
        && radmaskreti == other.radmaskreti
        && chromaskreti == other.chromaskreti
        && gammaskreti == other.gammaskreti
        && slomaskreti == other.slomaskreti
        && lapmaskreti == other.lapmaskreti
        && scalereti == other.scalereti
        && darkness == other.darkness
        && lightnessreti == other.lightnessreti
        && limd == other.limd
        && cliptm == other.cliptm
        && fftwreti == other.fftwreti
        && Lmaskreticurve == other.Lmaskreticurve
        // Sharpening
        && visisharp == other.visisharp
        && expsharp == other.expsharp
        && complexsharp == other.complexsharp
        && sharcontrast == other.sharcontrast
        && sharradius == other.sharradius
        && sharamount == other.sharamount
        && shardamping == other.shardamping
        && shariter == other.shariter
        && sharblur == other.sharblur
        && sensisha == other.sensisha
        && inverssha == other.inverssha
        // Local contrast
        && visicontrast == other.visicontrast
        && expcontrast == other.expcontrast
        && complexcontrast == other.complexcontrast
        && lcradius == other.lcradius
        && lcamount == other.lcamount
        && lcdarkness == other.lcdarkness
        && lclightness == other.lclightness
        && sigmalc == other.sigmalc
        && levelwav == other.levelwav
        && residcont == other.residcont
        && residsha == other.residsha
        && residshathr == other.residshathr
        && residhi == other.residhi
        && residhithr == other.residhithr
        && residblur == other.residblur
        && levelblur == other.levelblur
        && sigmabl == other.sigmabl
        && residchro == other.residchro
        && residcomp == other.residcomp
        && sigma == other.sigma
        && offset == other.offset
        && sigmadr == other.sigmadr
        && threswav == other.threswav
        && chromalev == other.chromalev
        && chromablu == other.chromablu
        && sigmadc == other.sigmadc
        && deltad == other.deltad
        && fatres == other.fatres
        && clarilres == other.clarilres
        && claricres == other.claricres
        && clarisoft == other.clarisoft
        && sigmalc2 == other.sigmalc2
        && strwav == other.strwav
        && angwav == other.angwav
        && strengthw == other.strengthw
        && sigmaed == other.sigmaed
        && radiusw == other.radiusw
        && detailw == other.detailw
        && gradw == other.gradw
        && tloww == other.tloww
        && thigw == other.thigw
        && edgw == other.edgw
        && basew == other.basew
        && sensilc == other.sensilc
        && fftwlc == other.fftwlc
        && blurlc == other.blurlc
        && wavblur == other.wavblur
        && wavedg == other.wavedg
        && waveshow == other.waveshow
        && wavcont == other.wavcont
        && wavcomp == other.wavcomp
        && wavgradl == other.wavgradl
        && wavcompre == other.wavcompre
        && origlc == other.origlc
        && localcontMethod == other.localcontMethod
        && localedgMethod == other.localedgMethod
        && localneiMethod == other.localneiMethod
        && locwavcurve == other.locwavcurve
        && csthreshold == other.csthreshold
        && loclevwavcurve == other.loclevwavcurve
        && locconwavcurve == other.locconwavcurve
        && loccompwavcurve == other.loccompwavcurve
        && loccomprewavcurve == other.loccomprewavcurve
        && locedgwavcurve == other.locedgwavcurve
        && CCmasklccurve == other.CCmasklccurve
        && LLmasklccurve == other.LLmasklccurve
        && HHmasklccurve == other.HHmasklccurve
        && enalcMask == other.enalcMask
        && blendmasklc == other.blendmasklc
        && radmasklc == other.radmasklc
        && chromasklc == other.chromasklc
        && Lmasklccurve == other.Lmasklccurve
        // Contrast by detail levels
        && visicbdl == other.visicbdl
        && expcbdl == other.expcbdl
        && complexcbdl == other.complexcbdl
        && [this, &other]() -> bool
            {
                for (int i = 0; i < 6; ++i) {
                    if (mult[i] != other.mult[i]) {
                        return false;
                    }
                }
                return true;
            }()
        && chromacbdl == other.chromacbdl
        && threshold == other.threshold
        && sensicb == other.sensicb
        && clarityml == other.clarityml
        && contresid == other.contresid
        && softradiuscb == other.softradiuscb
        && enacbMask == other.enacbMask
        && CCmaskcbcurve == other.CCmaskcbcurve
        && LLmaskcbcurve == other.LLmaskcbcurve
        && HHmaskcbcurve == other.HHmaskcbcurve
        && blendmaskcb == other.blendmaskcb
        && radmaskcb == other.radmaskcb
        && chromaskcb == other.chromaskcb
        && gammaskcb == other.gammaskcb
        && slomaskcb == other.slomaskcb
        && lapmaskcb == other.lapmaskcb
        && Lmaskcbcurve == other.Lmaskcbcurve
        // Log encoding
        && visilog == other.visilog
        && explog == other.explog
        && autocompute == other.autocompute
        && sourceGray == other.sourceGray
        && targetGray == other.targetGray
        && Autogray == other.Autogray
        && fullimage == other.fullimage
        && blackEv == other.blackEv
        && whiteEv == other.whiteEv
        && detail == other.detail
        && sensilog == other.sensilog
        && baselog == other.baselog
        && strlog == other.strlog
        && anglog == other.anglog
        // mask
        && visimask == other.visimask
        && complexmask == other.complexmask
        && expmask == other.expmask
        && sensimask == other.sensimask
        && blendmask == other.blendmask
        && blendmaskab == other.blendmaskab
        && softradiusmask == other.softradiusmask
        && enamask == other.enamask
        && fftmask == other.fftmask
        && blurmask == other.blurmask
        && contmask == other.contmask
        && CCmask_curve == other.CCmask_curve
        && LLmask_curve == other.LLmask_curve
        && HHmask_curve == other.HHmask_curve
        && strumaskmask == other.strumaskmask
        && toolmask == other.toolmask
        && radmask == other.radmask
        && lapmask == other.lapmask
        && chromask == other.chromask
        && gammask == other.gammask
        && slopmask == other.slopmask
        && shadmask == other.shadmask
        && str_mask == other.str_mask
        && ang_mask == other.ang_mask
        && HHhmask_curve == other.HHhmask_curve
        && Lmask_curve == other.Lmask_curve
        && LLmask_curvewav == other.LLmask_curvewav
        && csthresholdmask == other.csthresholdmask;

}

bool LocallabParams::LocallabSpot::operator !=(const LocallabSpot& other) const
{
    return !(*this == other);
}

const double LocallabParams::LABGRIDL_CORR_MAX = 12800.;
const double LocallabParams::LABGRIDL_CORR_SCALE = 3.276;
const double LocallabParams::LABGRIDL_DIRECT_SCALE = 41950.;

LocallabParams::LocallabParams() :
    enabled(false),
    selspot(0),
    spots()
{
}

bool LocallabParams::operator ==(const LocallabParams& other) const
{
    return
        enabled == other.enabled
        && selspot == other.selspot
        && spots == other.spots;
}

bool LocallabParams::operator !=(const LocallabParams& other) const
{
    return !(*this == other);
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
        && [this, &other]() -> bool
            {
                for (unsigned int i = 0; i < 6; ++i) {
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


SoftLightParams::SoftLightParams() :
    enabled(false),
    strength(30)
{
}

bool SoftLightParams::operator ==(const SoftLightParams& other) const
{
    return
        enabled == other.enabled
        && strength == other.strength;
}

bool SoftLightParams::operator !=(const SoftLightParams& other) const
{
    return !(*this == other);
}


DehazeParams::DehazeParams() :
    enabled(false),
    strength(50),
    showDepthMap(false),
    depth(25),
    luminance(false)
{
}

bool DehazeParams::operator ==(const DehazeParams& other) const
{
    return
        enabled == other.enabled
        && strength == other.strength
        && showDepthMap == other.showDepthMap
        && depth == other.depth
        && luminance == other.luminance;
}

bool DehazeParams::operator !=(const DehazeParams& other) const
{
    return !(*this == other);
}


RAWParams::BayerSensor::BayerSensor() :
    method(getMethodString(Method::AMAZE)),
    border(4),
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
    dualDemosaicAutoContrast(true),
    dualDemosaicContrast(20),
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
    pixelShiftEqualBright(false),
    pixelShiftEqualBrightChannel(false),
    pixelShiftNonGreenCross(true),
    pixelShiftDemosaicMethod(getPSDemosaicMethodString(PSDemosaicMethod::AMAZE)),
    dcb_enhance(true),
    pdafLinesFilter(false)
{
}

bool RAWParams::BayerSensor::operator ==(const BayerSensor& other) const
{
    return
        method == other.method
        && border == other.border
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
        && dualDemosaicAutoContrast == other.dualDemosaicAutoContrast
        && dualDemosaicContrast == other.dualDemosaicContrast
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
        && pixelShiftEqualBright == other.pixelShiftEqualBright
        && pixelShiftEqualBrightChannel == other.pixelShiftEqualBrightChannel
        && pixelShiftNonGreenCross == other.pixelShiftNonGreenCross
        && pixelShiftDemosaicMethod == other.pixelShiftDemosaicMethod
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
    pixelShiftEqualBright = false;
    pixelShiftEqualBrightChannel = false;
    pixelShiftNonGreenCross = true;
    pixelShiftDemosaicMethod = getPSDemosaicMethodString(PSDemosaicMethod::AMAZE);
}

const std::vector<const char*>& RAWParams::BayerSensor::getMethodStrings()
{
    static const std::vector<const char*> method_strings {
        "amaze",
        "amazevng4",
        "rcd",
        "rcdvng4",
        "dcb",
        "dcbvng4",
        "lmmse",
        "igv",
        "ahd",
        "eahd",
        "hphd",
        "vng4",
        "fast",
        "mono",
        "pixelshift",
        "none"
    };
    return method_strings;
}

Glib::ustring RAWParams::BayerSensor::getMethodString(Method method)
{
    return getMethodStrings()[toUnderlying(method)];
}

const std::vector<const char*>& RAWParams::BayerSensor::getPSDemosaicMethodStrings()
{
    static const std::vector<const char*> method_strings {
        "amaze",
        "amazevng4",
        "rcdvng4",
        "lmmse"
    };
    return method_strings;
}

Glib::ustring RAWParams::BayerSensor::getPSDemosaicMethodString(PSDemosaicMethod method)
{
    return getPSDemosaicMethodStrings()[toUnderlying(method)];
}

RAWParams::XTransSensor::XTransSensor() :
    method(getMethodString(Method::THREE_PASS)),
    dualDemosaicAutoContrast(true),
    dualDemosaicContrast(20),
    border(7),
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
        && dualDemosaicAutoContrast == other.dualDemosaicAutoContrast
        && dualDemosaicContrast == other.dualDemosaicContrast
        && border == other.border
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
        "4-pass",
        "3-pass (best)",
        "2-pass",
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


RAWParams::PreprocessWB::PreprocessWB() :
    mode(Mode::AUTO)
{
}

bool RAWParams::PreprocessWB::operator ==(const PreprocessWB& other) const
{
    return mode == other.mode;
}

bool RAWParams::PreprocessWB::operator !=(const PreprocessWB& other) const
{
    return !(*this == other);
}


RAWParams::RAWParams() :
    df_autoselect(false),
    ff_AutoSelect(false),
    ff_BlurRadius(32),
    ff_BlurType(getFlatFieldBlurTypeString(FlatFieldBlurType::AREA)),
    ff_AutoClipControl(false),
    ff_clipControl(0),
    ca_autocorrect(false),
    ca_avoidcolourshift(true),
    caautoiterations(2),
    cared(0.0),
    cablue(0.0),
    expos(1.0),
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
        && ca_avoidcolourshift == other.ca_avoidcolourshift
        && caautoiterations == other.caautoiterations
        && cared == other.cared
        && cablue == other.cablue
        && expos == other.expos
        && preprocessWB == other.preprocessWB
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

FilmNegativeParams::FilmNegativeParams() :
    enabled(false),
    redRatio(1.36),
    greenExp(1.5),
    blueRatio(0.86),
    redBase(0),
    greenBase(0),
    blueBase(0)
{
}

bool FilmNegativeParams::operator ==(const FilmNegativeParams& other) const
{
    return
        enabled == other.enabled
        && redRatio == other.redRatio
        && greenExp == other.greenExp
        && blueRatio == other.blueRatio
        && redBase == other.redBase
        && greenBase == other.greenBase
        && blueBase == other.blueBase;
}

bool FilmNegativeParams::operator !=(const FilmNegativeParams& other) const
{
    return !(*this == other);
}

ProcParams::ProcParams()
{
    setDefaults();
}

void ProcParams::setDefaults()
{
    toneCurve = {};

    labCurve = {};

    rgbCurves = {};

    localContrast = {};

    colorToning = {};

    sharpenEdge = {};

    sharpenMicro = {};

    sharpening = {};

    prsharpening = {};
    prsharpening.contrast = 15.0;
    prsharpening.method = "rld";
    prsharpening.deconvamount = 100;
    prsharpening.deconvradius = 0.45;
    prsharpening.deconviter = 100;
    prsharpening.deconvdamping = 0;

    pdsharpening = {};

    vibrance = {};

    wb = {};

    colorappearance = {};

    defringe = {};

    impulseDenoise = {};

    dirpyrDenoise = {};

    epd = {};

    fattal = {};

    sh = {};

    crop = {};

    coarse = {};

    commonTrans = {};

    rotate = {};

    distortion = {};

    lensProf = {};

    perspective = {};

    gradient = {};

    pcvignette = {};

    vignetting = {};

    locallab = {};

    chmixer = {};

    blackwhite = {};

    cacorrection = {};

    resize = {};

    icm = {};

    wavelet = {};

    dirpyrequalizer = {};

    hsvequalizer = {};

    filmSimulation = {};

    softlight = {};

    dehaze = {};

    raw = {};

    metadata = {};
    exif.clear();
    iptc.clear();

    // -1 means that there's no pp3 data with rank yet. In this case, the
    // embedded Rating metadata should take precedence. -1 should never be
    // written to pp3 on disk.
    rank = -1;
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
        saveToKeyfile(!pedited || pedited->toneCurve.fromHistMatching, "Exposure", "CurveFromHistogramMatching", toneCurve.fromHistMatching, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.clampOOG, "Exposure", "ClampOOG", toneCurve.clampOOG, keyFile);

// Highlight recovery
        saveToKeyfile(!pedited || pedited->toneCurve.hrenabled, "HLRecovery", "Enabled", toneCurve.hrenabled, keyFile);
        saveToKeyfile(!pedited || pedited->toneCurve.method, "HLRecovery", "Method", toneCurve.method, keyFile);

        const std::map<ToneCurveMode, const char*> tc_mapping = {
            {ToneCurveMode::STD, "Standard"},
            {ToneCurveMode::FILMLIKE, "FilmLike"},
            {ToneCurveMode::SATANDVALBLENDING, "SatAndValueBlending"},
            {ToneCurveMode::WEIGHTEDSTD, "WeightedStd"},
            {ToneCurveMode::LUMINANCE, "Luminance"},
            {ToneCurveMode::PERCEPTUAL, "Perceptual"}
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
        saveToKeyfile(!pedited || pedited->retinex.complexmethod, "Retinex", "complexMethod", retinex.complexmethod, keyFile);
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
            "BeforeCurveMode",
            {
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
            "AfterCurveMode",
            {
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
        saveToKeyfile(!pedited || pedited->sharpening.contrast, "Sharpening", "Contrast", sharpening.contrast, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.method, "Sharpening", "Method", sharpening.method, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.radius, "Sharpening", "Radius", sharpening.radius, keyFile);
        saveToKeyfile(!pedited || pedited->sharpening.blurradius, "Sharpening", "BlurRadius", sharpening.blurradius, keyFile);
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
        saveToKeyfile(!pedited || pedited->sharpenMicro.contrast, "SharpenMicro", "Contrast", sharpenMicro.contrast, keyFile);
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
        saveToKeyfile(!pedited || pedited->colorappearance.illum, "Color appearance", "Illum", colorappearance.illum, keyFile);
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
        saveToKeyfile(!pedited || pedited->colorappearance.autotempout, "Color appearance", "Autotempout", colorappearance.autotempout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.greenout, "Color appearance", "Greenout", colorappearance.greenout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.tempsc, "Color appearance", "Tempsc", colorappearance.tempsc, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.greensc, "Color appearance", "Greensc", colorappearance.greensc, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.ybout, "Color appearance", "Ybout", colorappearance.ybout, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.datacie, "Color appearance", "Datacie", colorappearance.datacie, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.tonecie, "Color appearance", "Tonecie", colorappearance.tonecie, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.presetcat02, "Color appearance", "Presetcat02", colorappearance.presetcat02, keyFile);

        const std::map<ColorAppearanceParams::TcMode, const char*> ca_mapping = {
            {ColorAppearanceParams::TcMode::LIGHT, "Lightness"},
            {ColorAppearanceParams::TcMode::BRIGHT, "Brightness"}
        };

        saveToKeyfile(!pedited || pedited->colorappearance.curveMode, "Color appearance", "CurveMode", ca_mapping, colorappearance.curveMode, keyFile);
        saveToKeyfile(!pedited || pedited->colorappearance.curveMode2, "Color appearance", "CurveMode2", ca_mapping, colorappearance.curveMode2, keyFile);
        saveToKeyfile(
            !pedited || pedited->colorappearance.curveMode3,
            "Color appearance",
            "CurveMode3",
            {
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

// Dehaze
        saveToKeyfile(!pedited || pedited->dehaze.enabled, "Dehaze", "Enabled", dehaze.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->dehaze.strength, "Dehaze", "Strength", dehaze.strength, keyFile);        
        saveToKeyfile(!pedited || pedited->dehaze.showDepthMap, "Dehaze", "ShowDepthMap", dehaze.showDepthMap, keyFile);        
        saveToKeyfile(!pedited || pedited->dehaze.depth, "Dehaze", "Depth", dehaze.depth, keyFile);        
        saveToKeyfile(!pedited || pedited->dehaze.depth, "Dehaze", "Luminance", dehaze.luminance, keyFile);

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
        saveToKeyfile(!pedited || pedited->sh.lab, "Shadows & Highlights", "Lab", sh.lab, keyFile);

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
        saveToKeyfile(!pedited || pedited->commonTrans.method, "Common Properties for Transformations", "Method", commonTrans.method, keyFile);
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
        saveToKeyfile(!pedited || pedited->perspective.method, "Perspective", "Method", perspective.method, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.horizontal, "Perspective", "Horizontal", perspective.horizontal, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.vertical, "Perspective", "Vertical", perspective.vertical, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.camera_crop_factor, "Perspective", "CameraCropFactor", perspective.camera_crop_factor, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.camera_focal_length, "Perspective", "CameraFocalLength", perspective.camera_focal_length, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.camera_pitch, "Perspective", "CameraPitch", perspective.camera_pitch, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.camera_roll, "Perspective", "CameraRoll", perspective.camera_roll, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.camera_shift_horiz, "Perspective", "CameraShiftHorizontal", perspective.camera_shift_horiz, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.camera_shift_vert, "Perspective", "CameraShiftVertical", perspective.camera_shift_vert, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.camera_yaw, "Perspective", "CameraYaw", perspective.camera_yaw, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.projection_shift_horiz, "Perspective", "ProjectionShiftHorizontal", perspective.projection_shift_horiz, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.projection_pitch, "Perspective", "ProjectionPitch", perspective.projection_pitch, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.projection_rotate, "Perspective", "ProjectionRotate", perspective.projection_rotate, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.projection_shift_horiz, "Perspective", "ProjectionShiftHorizontal", perspective.projection_shift_horiz, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.projection_shift_vert, "Perspective", "ProjectionShiftVertical", perspective.projection_shift_vert, keyFile);
        saveToKeyfile(!pedited || pedited->perspective.projection_yaw, "Perspective", "ProjectionYaw", perspective.projection_yaw, keyFile);

// Gradient
        saveToKeyfile(!pedited || pedited->gradient.enabled, "Gradient", "Enabled", gradient.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.degree, "Gradient", "Degree", gradient.degree, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.feather, "Gradient", "Feather", gradient.feather, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.strength, "Gradient", "Strength", gradient.strength, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.centerX, "Gradient", "CenterX", gradient.centerX, keyFile);
        saveToKeyfile(!pedited || pedited->gradient.centerY, "Gradient", "CenterY", gradient.centerY, keyFile);

// Locallab
        saveToKeyfile(!pedited || pedited->locallab.enabled, "Locallab", "Enabled", locallab.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->locallab.selspot, "Locallab", "Selspot", locallab.selspot, keyFile);

        for (size_t i = 0; i < locallab.spots.size(); ++i) {
            if (!pedited || i < pedited->locallab.spots.size()) {
                const LocallabParams::LocallabSpot& spot = locallab.spots.at(i);
                const LocallabParamsEdited::LocallabSpotEdited* const spot_edited =
                    pedited
                        ? &pedited->locallab.spots.at(i)
                        : nullptr;
                const std::string index_str = std::to_string(i);
                // Control spot settings
                saveToKeyfile(!pedited || spot_edited->name, "Locallab", "Name_" + index_str, spot.name, keyFile);
                saveToKeyfile(!pedited || spot_edited->isvisible, "Locallab", "Isvisible_" + index_str, spot.isvisible, keyFile);
                saveToKeyfile(!pedited || spot_edited->prevMethod, "Locallab", "PrevMethod_" + index_str, spot.prevMethod, keyFile);
                saveToKeyfile(!pedited || spot_edited->shape, "Locallab", "Shape_" + index_str, spot.shape, keyFile);
                saveToKeyfile(!pedited || spot_edited->spotMethod, "Locallab", "SpotMethod_" + index_str, spot.spotMethod, keyFile);
                saveToKeyfile(!pedited || spot_edited->wavMethod, "Locallab", "WavMethod_" + index_str, spot.wavMethod, keyFile);
                saveToKeyfile(!pedited || spot_edited->sensiexclu, "Locallab", "SensiExclu_" + index_str, spot.sensiexclu, keyFile);
                saveToKeyfile(!pedited || spot_edited->structexclu, "Locallab", "StructExclu_" + index_str, spot.structexclu, keyFile);
                saveToKeyfile(!pedited || spot_edited->struc, "Locallab", "Struc_" + index_str, spot.struc, keyFile);
                saveToKeyfile(!pedited || spot_edited->shapeMethod, "Locallab", "ShapeMethod_" + index_str, spot.shapeMethod, keyFile);
                saveToKeyfile(!pedited || spot_edited->loc, "Locallab", "Loc_" + index_str, spot.loc, keyFile);
                saveToKeyfile(!pedited || spot_edited->centerX, "Locallab", "CenterX_" + index_str, spot.centerX, keyFile);
                saveToKeyfile(!pedited || spot_edited->centerY, "Locallab", "CenterY_" + index_str, spot.centerY, keyFile);
                saveToKeyfile(!pedited || spot_edited->circrad, "Locallab", "Circrad_" + index_str, spot.circrad, keyFile);
                saveToKeyfile(!pedited || spot_edited->qualityMethod, "Locallab", "QualityMethod_" + index_str, spot.qualityMethod, keyFile);
                saveToKeyfile(!pedited || spot_edited->complexMethod, "Locallab", "ComplexMethod_" + index_str, spot.complexMethod, keyFile);
                saveToKeyfile(!pedited || spot_edited->transit, "Locallab", "Transit_" + index_str, spot.transit, keyFile);
                saveToKeyfile(!pedited || spot_edited->feather, "Locallab", "Feather_" + index_str, spot.feather, keyFile);
                saveToKeyfile(!pedited || spot_edited->thresh, "Locallab", "Thresh_" + index_str, spot.thresh, keyFile);
                saveToKeyfile(!pedited || spot_edited->iter, "Locallab", "Iter_" + index_str, spot.iter, keyFile);
                saveToKeyfile(!pedited || spot_edited->balan, "Locallab", "Balan_" + index_str, spot.balan, keyFile);
                saveToKeyfile(!pedited || spot_edited->balanh, "Locallab", "Balanh_" + index_str, spot.balanh, keyFile);
                saveToKeyfile(!pedited || spot_edited->colorde, "Locallab", "Colorde_" + index_str, spot.colorde, keyFile);
                saveToKeyfile(!pedited || spot_edited->colorscope, "Locallab", "Colorscope_" + index_str, spot.colorscope, keyFile);
                saveToKeyfile(!pedited || spot_edited->transitweak, "Locallab", "Transitweak_" + index_str, spot.transitweak, keyFile);
                saveToKeyfile(!pedited || spot_edited->transitgrad, "Locallab", "Transitgrad_" + index_str, spot.transitgrad, keyFile);
                saveToKeyfile(!pedited || spot_edited->activ, "Locallab", "Activ_" + index_str, spot.activ, keyFile);
                saveToKeyfile(!pedited || spot_edited->avoid, "Locallab", "Avoid_" + index_str, spot.avoid, keyFile);
                saveToKeyfile(!pedited || spot_edited->blwh, "Locallab", "Blwh_" + index_str, spot.blwh, keyFile);
                saveToKeyfile(!pedited || spot_edited->recurs, "Locallab", "Recurs_" + index_str, spot.recurs, keyFile);
                saveToKeyfile(!pedited || spot_edited->laplac, "Locallab", "Laplac_" + index_str, spot.laplac, keyFile);
                saveToKeyfile(!pedited || spot_edited->deltae, "Locallab", "Deltae_" + index_str, spot.deltae, keyFile);
                saveToKeyfile(!pedited || spot_edited->shortc, "Locallab", "Shortc_" + index_str, spot.shortc, keyFile);
                saveToKeyfile(!pedited || spot_edited->savrest, "Locallab", "Savrest_" + index_str, spot.savrest, keyFile);
                saveToKeyfile(!pedited || spot_edited->scopemask, "Locallab", "Scopemask_" + index_str, spot.scopemask, keyFile);
                saveToKeyfile(!pedited || spot_edited->lumask, "Locallab", "Lumask_" + index_str, spot.lumask, keyFile);
                // Color & Light
                if ((!pedited || spot_edited->visicolor) && spot.visicolor) {
                    saveToKeyfile(!pedited || spot_edited->expcolor, "Locallab", "Expcolor_" + index_str, spot.expcolor, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexcolor, "Locallab", "Complexcolor_" + index_str, spot.complexcolor, keyFile);
                    saveToKeyfile(!pedited || spot_edited->curvactiv, "Locallab", "Curvactiv_" + index_str, spot.curvactiv, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lightness, "Locallab", "Lightness_" + index_str, spot.lightness, keyFile);
                    saveToKeyfile(!pedited || spot_edited->contrast, "Locallab", "Contrast_" + index_str, spot.contrast, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chroma, "Locallab", "Chroma_" + index_str, spot.chroma, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridALow, "Locallab", "labgridALow_" + index_str, spot.labgridALow, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridBLow, "Locallab", "labgridBLow_" + index_str, spot.labgridBLow, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridAHigh, "Locallab", "labgridAHigh_" + index_str, spot.labgridAHigh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridBHigh, "Locallab", "labgridBHigh_" + index_str, spot.labgridBHigh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridALowmerg, "Locallab", "labgridALowmerg_" + index_str, spot.labgridALowmerg, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridBLowmerg, "Locallab", "labgridBLowmerg_" + index_str, spot.labgridBLowmerg, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridAHighmerg, "Locallab", "labgridAHighmerg_" + index_str, spot.labgridAHighmerg, keyFile);
                    saveToKeyfile(!pedited || spot_edited->labgridBHighmerg, "Locallab", "labgridBHighmerg_" + index_str, spot.labgridBHighmerg, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strengthgrid, "Locallab", "Strengthgrid_" + index_str, spot.strengthgrid, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensi, "Locallab", "Sensi_" + index_str, spot.sensi, keyFile);
                    saveToKeyfile(!pedited || spot_edited->structcol, "Locallab", "Structcol_" + index_str, spot.structcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strcol, "Locallab", "Strcol_" + index_str, spot.strcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strcolab, "Locallab", "Strcolab_" + index_str, spot.strcolab, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strcolh, "Locallab", "Strcolh_" + index_str, spot.strcolh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->angcol, "Locallab", "Angcol_" + index_str, spot.angcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blurcolde, "Locallab", "Blurcolde_" + index_str, spot.blurcolde, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blurcol, "Locallab", "Blurcol_" + index_str, spot.blurcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->contcol, "Locallab", "Contcol_" + index_str, spot.contcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskcol, "Locallab", "Blendmaskcol_" + index_str, spot.blendmaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmaskcol, "Locallab", "Radmaskcol_" + index_str, spot.radmaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromaskcol, "Locallab", "Chromaskcol_" + index_str, spot.chromaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammaskcol, "Locallab", "Gammaskcol_" + index_str, spot.gammaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomaskcol, "Locallab", "Slomaskcol_" + index_str, spot.slomaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shadmaskcol, "Locallab", "shadmaskcol_" + index_str, spot.shadmaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strumaskcol, "Locallab", "strumaskcol_" + index_str, spot.strumaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmaskcol, "Locallab", "Lapmaskcol_" + index_str, spot.lapmaskcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->qualitycurveMethod, "Locallab", "QualityCurveMethod_" + index_str, spot.qualitycurveMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gridMethod, "Locallab", "gridMethod_" + index_str, spot.gridMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->merMethod, "Locallab", "Merg_Method_" + index_str, spot.merMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->toneMethod, "Locallab", "ToneMethod_" + index_str, spot.toneMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->mergecolMethod, "Locallab", "mergecolMethod_" + index_str, spot.mergecolMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->llcurve, "Locallab", "LLCurve_" + index_str, spot.llcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lccurve, "Locallab", "LCCurve_" + index_str, spot.lccurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->cccurve, "Locallab", "CCCurve_" + index_str, spot.cccurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->clcurve, "Locallab", "CLCurve_" + index_str, spot.clcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->rgbcurve, "Locallab", "RGBCurve_" + index_str, spot.rgbcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LHcurve, "Locallab", "LHCurve_" + index_str, spot.LHcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHcurve, "Locallab", "HHCurve_" + index_str, spot.HHcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CHcurve, "Locallab", "CHCurve_" + index_str, spot.CHcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->invers, "Locallab", "Invers_" + index_str, spot.invers, keyFile);
                    saveToKeyfile(!pedited || spot_edited->special, "Locallab", "Special_" + index_str, spot.special, keyFile);
                    saveToKeyfile(!pedited || spot_edited->toolcol, "Locallab", "Toolcol_" + index_str, spot.toolcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enaColorMask, "Locallab", "EnaColorMask_" + index_str, spot.enaColorMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fftColorMask, "Locallab", "FftColorMask_" + index_str, spot.fftColorMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmaskcurve, "Locallab", "CCmaskCurve_" + index_str, spot.CCmaskcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskcurve, "Locallab", "LLmaskCurve_" + index_str, spot.LLmaskcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmaskcurve, "Locallab", "HHmaskCurve_" + index_str, spot.HHmaskcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHhmaskcurve, "Locallab", "HHhmaskCurve_" + index_str, spot.HHhmaskcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->softradiuscol, "Locallab", "Softradiuscol_" + index_str, spot.softradiuscol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->opacol, "Locallab", "Opacol_" + index_str, spot.opacol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->mercol, "Locallab", "Mercol_" + index_str, spot.mercol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->merlucol, "Locallab", "Merlucol_" + index_str, spot.merlucol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->conthrcol, "Locallab", "Conthrcol_" + index_str, spot.conthrcol, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmaskcurve, "Locallab", "LmaskCurve_" + index_str, spot.Lmaskcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskcolcurvewav, "Locallab", "LLmaskcolCurvewav_" + index_str, spot.LLmaskcolcurvewav, keyFile);
                    saveToKeyfile(!pedited || spot_edited->csthresholdcol, "Locallab", "CSThresholdcol_" + index_str, spot.csthresholdcol.toVector(), keyFile);
                }
                // Exposure
                if ((!pedited || spot_edited->visiexpose) && spot.visiexpose) {
                    saveToKeyfile(!pedited || spot_edited->expexpose, "Locallab", "Expexpose_" + index_str, spot.expexpose, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexexpose, "Locallab", "Complexexpose_" + index_str, spot.complexexpose, keyFile);
                    saveToKeyfile(!pedited || spot_edited->expcomp, "Locallab", "Expcomp_" + index_str, spot.expcomp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->hlcompr, "Locallab", "Hlcompr_" + index_str, spot.hlcompr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->hlcomprthresh, "Locallab", "Hlcomprthresh_" + index_str, spot.hlcomprthresh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->black, "Locallab", "Black_" + index_str, spot.black, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shadex, "Locallab", "Shadex_" + index_str, spot.shadex, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shcompr, "Locallab", "Shcompr_" + index_str, spot.shcompr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->expchroma, "Locallab", "Expchroma_" + index_str, spot.expchroma, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensiex, "Locallab", "Sensiex_" + index_str, spot.sensiex, keyFile);
                    saveToKeyfile(!pedited || spot_edited->structexp, "Locallab", "Structexp_" + index_str, spot.structexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blurexpde, "Locallab", "Blurexpde_" + index_str, spot.blurexpde, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strexp, "Locallab", "Strexp_" + index_str, spot.strexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->angexp, "Locallab", "Angexp_" + index_str, spot.angexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->excurve, "Locallab", "ExCurve_" + index_str, spot.excurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->inversex, "Locallab", "Inversex_" + index_str, spot.inversex, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enaExpMask, "Locallab", "EnaExpMask_" + index_str, spot.enaExpMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enaExpMaskaft, "Locallab", "EnaExpMaskaft_" + index_str, spot.enaExpMaskaft, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmaskexpcurve, "Locallab", "CCmaskexpCurve_" + index_str, spot.CCmaskexpcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskexpcurve, "Locallab", "LLmaskexpCurve_" + index_str, spot.LLmaskexpcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmaskexpcurve, "Locallab", "HHmaskexpCurve_" + index_str, spot.HHmaskexpcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskexp, "Locallab", "Blendmaskexp_" + index_str, spot.blendmaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmaskexp, "Locallab", "Radmaskexp_" + index_str, spot.radmaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromaskexp, "Locallab", "Chromaskexp_" + index_str, spot.chromaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammaskexp, "Locallab", "Gammaskexp_" + index_str, spot.gammaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomaskexp, "Locallab", "Slomaskexp_" + index_str, spot.slomaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmaskexp, "Locallab", "Lapmaskexp_" + index_str, spot.lapmaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strmaskexp, "Locallab", "Strmaskexp_" + index_str, spot.strmaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->angmaskexp, "Locallab", "Angmaskexp_" + index_str, spot.angmaskexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->softradiusexp, "Locallab", "Softradiusexp_" + index_str, spot.softradiusexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmaskexpcurve, "Locallab", "LmaskexpCurve_" + index_str, spot.Lmaskexpcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->expMethod, "Locallab", "ExpMethod_" + index_str, spot.expMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->exnoiseMethod, "Locallab", "ExnoiseMethod_" + index_str, spot.exnoiseMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->laplacexp, "Locallab", "Laplacexp_" + index_str, spot.laplacexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->balanexp, "Locallab", "Balanexp_" + index_str, spot.balanexp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->linear, "Locallab", "Linearexp_" + index_str, spot.linear, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gamm, "Locallab", "Gamm_" + index_str, spot.gamm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fatamount, "Locallab", "Fatamount_" + index_str, spot.fatamount, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fatdetail, "Locallab", "Fatdetail_" + index_str, spot.fatdetail, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fatanchor, "Locallab", "Fatanchor_" + index_str, spot.fatanchor, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fatlevel, "Locallab", "Fatlevel_" + index_str, spot.fatlevel, keyFile);
                }
                // Shadow highlight
                if ((!pedited || spot_edited->visishadhigh) && spot.visishadhigh) {
                    saveToKeyfile(!pedited || spot_edited->expshadhigh, "Locallab", "Expshadhigh_" + index_str, spot.expshadhigh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexshadhigh, "Locallab", "Complexshadhigh_" + index_str, spot.complexshadhigh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shMethod, "Locallab", "ShMethod_" + index_str, spot.shMethod, keyFile);

                    for (int j = 0; j < 5; j++) {
                        saveToKeyfile(!pedited || spot_edited->multsh[j], "Locallab", "Multsh" + std::to_string(j) + "_" + index_str, spot.multsh[j], keyFile);
                    }

                    saveToKeyfile(!pedited || spot_edited->highlights, "Locallab", "highlights_" + index_str, spot.highlights, keyFile);
                    saveToKeyfile(!pedited || spot_edited->h_tonalwidth, "Locallab", "h_tonalwidth_" + index_str, spot.h_tonalwidth, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shadows, "Locallab", "shadows_" + index_str, spot.shadows, keyFile);
                    saveToKeyfile(!pedited || spot_edited->s_tonalwidth, "Locallab", "s_tonalwidth_" + index_str, spot.s_tonalwidth, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sh_radius, "Locallab", "sh_radius_" + index_str, spot.sh_radius, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensihs, "Locallab", "sensihs_" + index_str, spot.sensihs, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enaSHMask, "Locallab", "EnaSHMask_" + index_str, spot.enaSHMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmaskSHcurve, "Locallab", "CCmaskSHCurve_" + index_str, spot.CCmaskSHcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskSHcurve, "Locallab", "LLmaskSHCurve_" + index_str, spot.LLmaskSHcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmaskSHcurve, "Locallab", "HHmaskSHCurve_" + index_str, spot.HHmaskSHcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskSH, "Locallab", "BlendmaskSH_" + index_str, spot.blendmaskSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmaskSH, "Locallab", "RadmaskSH_" + index_str, spot.radmaskSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blurSHde, "Locallab", "BlurSHde_" + index_str, spot.blurSHde, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strSH, "Locallab", "StrSH_" + index_str, spot.strSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->angSH, "Locallab", "AngSH_" + index_str, spot.angSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->inverssh, "Locallab", "Inverssh_" + index_str, spot.inverssh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromaskSH, "Locallab", "ChromaskSH_" + index_str, spot.chromaskSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammaskSH, "Locallab", "GammaskSH_" + index_str, spot.gammaskSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomaskSH, "Locallab", "SlomaskSH_" + index_str, spot.slomaskSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->detailSH, "Locallab", "DetailSH_" + index_str, spot.detailSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LmaskSHcurve, "Locallab", "LmaskSHCurve_" + index_str, spot.LmaskSHcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fatamountSH, "Locallab", "FatamountSH_" + index_str, spot.fatamountSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fatanchorSH, "Locallab", "FatanchorSH_" + index_str, spot.fatanchorSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gamSH, "Locallab", "GamSH_" + index_str, spot.gamSH, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sloSH, "Locallab", "SloSH_" + index_str, spot.sloSH, keyFile);
                }
                // Vibrance
                if ((!pedited || spot_edited->visivibrance) && spot.visivibrance) {
                    saveToKeyfile(!pedited || spot_edited->expvibrance, "Locallab", "Expvibrance_" + index_str, spot.expvibrance, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexvibrance, "Locallab", "Complexvibrance_" + index_str, spot.complexvibrance, keyFile);
                    saveToKeyfile(!pedited || spot_edited->saturated, "Locallab", "Saturated_" + index_str, spot.saturated, keyFile);
                    saveToKeyfile(!pedited || spot_edited->pastels, "Locallab", "Pastels_" + index_str, spot.pastels, keyFile);
                    saveToKeyfile(!pedited || spot_edited->warm, "Locallab", "Warm_" + index_str, spot.warm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->psthreshold, "Locallab", "PSThreshold_" + index_str, spot.psthreshold.toVector(), keyFile);
                    saveToKeyfile(!pedited || spot_edited->protectskins, "Locallab", "ProtectSkins_" + index_str, spot.protectskins, keyFile);
                    saveToKeyfile(!pedited || spot_edited->avoidcolorshift, "Locallab", "AvoidColorShift_" + index_str, spot.avoidcolorshift, keyFile);
                    saveToKeyfile(!pedited || spot_edited->pastsattog, "Locallab", "PastSatTog_" + index_str, spot.pastsattog, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensiv, "Locallab", "Sensiv_" + index_str, spot.sensiv, keyFile);
                    saveToKeyfile(!pedited || spot_edited->skintonescurve, "Locallab", "SkinTonesCurve_" + index_str, spot.skintonescurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmaskvibcurve, "Locallab", "CCmaskvibCurve_" + index_str, spot.CCmaskvibcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskvibcurve, "Locallab", "LLmaskvibCurve_" + index_str, spot.LLmaskvibcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmaskvibcurve, "Locallab", "HHmaskvibCurve_" + index_str, spot.HHmaskvibcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enavibMask, "Locallab", "EnavibMask_" + index_str, spot.enavibMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskvib, "Locallab", "Blendmaskvib_" + index_str, spot.blendmaskvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmaskvib, "Locallab", "Radmaskvib_" + index_str, spot.radmaskvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromaskvib, "Locallab", "Chromaskvib_" + index_str, spot.chromaskvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammaskvib, "Locallab", "Gammaskvib_" + index_str, spot.gammaskvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomaskvib, "Locallab", "Slomaskvib_" + index_str, spot.slomaskvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmaskvib, "Locallab", "Lapmaskvib_" + index_str, spot.lapmaskvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strvib, "Locallab", "Strvib_" + index_str, spot.strvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strvibab, "Locallab", "Strvibab_" + index_str, spot.strvibab, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strvibh, "Locallab", "Strvibh_" + index_str, spot.strvibh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->angvib, "Locallab", "Angvib_" + index_str, spot.angvib, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmaskvibcurve, "Locallab", "LmaskvibCurve_" + index_str, spot.Lmaskvibcurve, keyFile);
                }
                // Soft Light
                if ((!pedited || spot_edited->visisoft) && spot.visisoft) {
                    saveToKeyfile(!pedited || spot_edited->expsoft, "Locallab", "Expsoft_" + index_str, spot.expsoft, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexsoft, "Locallab", "Complexsoft_" + index_str, spot.complexsoft, keyFile);
                    saveToKeyfile(!pedited || spot_edited->streng, "Locallab", "Streng_" + index_str, spot.streng, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensisf, "Locallab", "Sensisf_" + index_str, spot.sensisf, keyFile);
                    saveToKeyfile(!pedited || spot_edited->laplace, "Locallab", "Laplace_" + index_str, spot.laplace, keyFile);
                    saveToKeyfile(!pedited || spot_edited->softMethod, "Locallab", "SoftMethod_" + index_str, spot.softMethod, keyFile);
                }
                // Blur & Noise
                if ((!pedited || spot_edited->visiblur) && spot.visiblur) {
                    saveToKeyfile(!pedited || spot_edited->expblur, "Locallab", "Expblur_" + index_str, spot.expblur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexblur, "Locallab", "Complexblur_" + index_str, spot.complexblur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radius, "Locallab", "Radius_" + index_str, spot.radius, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strength, "Locallab", "Strength_" + index_str, spot.strength, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensibn, "Locallab", "Sensibn_" + index_str, spot.sensibn, keyFile);
                    saveToKeyfile(!pedited || spot_edited->itera, "Locallab", "Iteramed_" + index_str, spot.itera, keyFile);
                    saveToKeyfile(!pedited || spot_edited->guidbl, "Locallab", "Guidbl_" + index_str, spot.guidbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strbl, "Locallab", "Strbl_" + index_str, spot.strbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->isogr, "Locallab", "Isogr_" + index_str, spot.isogr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strengr, "Locallab", "Strengr_" + index_str, spot.strengr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->scalegr, "Locallab", "Scalegr_" + index_str, spot.scalegr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->epsbl, "Locallab", "Epsbl_" + index_str, spot.epsbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blMethod, "Locallab", "BlMethod_" + index_str, spot.blMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chroMethod, "Locallab", "ChroMethod_" + index_str, spot.chroMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blurMethod, "Locallab", "BlurMethod_" + index_str, spot.blurMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->medMethod, "Locallab", "MedMethod_" + index_str, spot.medMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->activlum, "Locallab", "activlum_" + index_str, spot.activlum, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noiselumf, "Locallab", "noiselumf_" + index_str, spot.noiselumf, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noiselumf0, "Locallab", "noiselumf0_" + index_str, spot.noiselumf0, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noiselumf2, "Locallab", "noiselumf2_" + index_str, spot.noiselumf2, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noiselumc, "Locallab", "noiselumc_" + index_str, spot.noiselumc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noiselumdetail, "Locallab", "noiselumdetail_" + index_str, spot.noiselumdetail, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noiselequal, "Locallab", "noiselequal_" + index_str, spot.noiselequal, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noisechrof, "Locallab", "noisechrof_" + index_str, spot.noisechrof, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noisechroc, "Locallab", "noisechroc_" + index_str, spot.noisechroc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->noisechrodetail, "Locallab", "noisechrodetail_" + index_str, spot.noisechrodetail, keyFile);
                    saveToKeyfile(!pedited || spot_edited->adjblur, "Locallab", "Adjblur_" + index_str, spot.adjblur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->bilateral, "Locallab", "Bilateral_" + index_str, spot.bilateral, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensiden, "Locallab", "Sensiden_" + index_str, spot.sensiden, keyFile);
                    saveToKeyfile(!pedited || spot_edited->detailthr, "Locallab", "Detailthr_" + index_str, spot.detailthr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->locwavcurveden, "Locallab", "LocwavCurveden_" + index_str, spot.locwavcurveden, keyFile);
                    saveToKeyfile(!pedited || spot_edited->showmaskblMethodtyp, "Locallab", "Showmasktyp_" + index_str, spot.showmaskblMethodtyp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmaskblcurve, "Locallab", "CCmaskblCurve_" + index_str, spot.CCmaskblcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskblcurve, "Locallab", "LLmaskblCurve_" + index_str, spot.LLmaskblcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmaskblcurve, "Locallab", "HHmaskblCurve_" + index_str, spot.HHmaskblcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enablMask, "Locallab", "EnablMask_" + index_str, spot.enablMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fftwbl, "Locallab", "Fftwbl_" + index_str, spot.fftwbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->toolbl, "Locallab", "Toolbl_" + index_str, spot.toolbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskbl, "Locallab", "Blendmaskbl_" + index_str, spot.blendmaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmaskbl, "Locallab", "Radmaskbl_" + index_str, spot.radmaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromaskbl, "Locallab", "Chromaskbl_" + index_str, spot.chromaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammaskbl, "Locallab", "Gammaskbl_" + index_str, spot.gammaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomaskbl, "Locallab", "Slomaskbl_" + index_str, spot.slomaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmaskbl, "Locallab", "Lapmaskbl_" + index_str, spot.lapmaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shadmaskbl, "Locallab", "shadmaskbl_" + index_str, spot.shadmaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shadmaskblsha, "Locallab", "shadmaskblsha_" + index_str, spot.shadmaskblsha, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strumaskbl, "Locallab", "strumaskbl_" + index_str, spot.strumaskbl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmaskblcurve, "Locallab", "LmaskblCurve_" + index_str, spot.Lmaskblcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskblcurvewav, "Locallab", "LLmaskblCurvewav_" + index_str, spot.LLmaskblcurvewav, keyFile);
                    saveToKeyfile(!pedited || spot_edited->csthresholdblur, "Locallab", "CSThresholdblur_" + index_str, spot.csthresholdblur.toVector(), keyFile);
                }
                // Tone Mapping
                if ((!pedited || spot_edited->visitonemap) && spot.visitonemap) {
                    saveToKeyfile(!pedited || spot_edited->exptonemap, "Locallab", "Exptonemap_" + index_str, spot.exptonemap, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complextonemap, "Locallab", "Complextonemap_" + index_str, spot.complextonemap, keyFile);
                    saveToKeyfile(!pedited || spot_edited->stren, "Locallab", "Stren_" + index_str, spot.stren, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gamma, "Locallab", "Gamma_" + index_str, spot.gamma, keyFile);
                    saveToKeyfile(!pedited || spot_edited->estop, "Locallab", "Estop_" + index_str, spot.estop, keyFile);
                    saveToKeyfile(!pedited || spot_edited->scaltm, "Locallab", "Scaltm_" + index_str, spot.scaltm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->rewei, "Locallab", "Rewei_" + index_str, spot.rewei, keyFile);
                    saveToKeyfile(!pedited || spot_edited->satur, "Locallab", "Satur_" + index_str, spot.satur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensitm, "Locallab", "Sensitm_" + index_str, spot.sensitm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->softradiustm, "Locallab", "Softradiustm_" + index_str, spot.softradiustm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->amount, "Locallab", "Amount_" + index_str, spot.amount, keyFile);
                    saveToKeyfile(!pedited || spot_edited->equiltm, "Locallab", "Equiltm_" + index_str, spot.equiltm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmasktmcurve, "Locallab", "CCmasktmCurve_" + index_str, spot.CCmasktmcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmasktmcurve, "Locallab", "LLmasktmCurve_" + index_str, spot.LLmasktmcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmasktmcurve, "Locallab", "HHmasktmCurve_" + index_str, spot.HHmasktmcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enatmMask, "Locallab", "EnatmMask_" + index_str, spot.enatmMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enatmMaskaft, "Locallab", "EnatmMaskaft_" + index_str, spot.enatmMaskaft, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmasktm, "Locallab", "Blendmasktm_" + index_str, spot.blendmasktm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmasktm, "Locallab", "Radmasktm_" + index_str, spot.radmasktm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromasktm, "Locallab", "Chromasktm_" + index_str, spot.chromasktm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammasktm, "Locallab", "Gammasktm_" + index_str, spot.gammasktm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomasktm, "Locallab", "Slomasktm_" + index_str, spot.slomasktm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmasktm, "Locallab", "Lapmasktm_" + index_str, spot.lapmasktm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmasktmcurve, "Locallab", "LmasktmCurve_" + index_str, spot.Lmasktmcurve, keyFile);
                }
                // Retinex
                if ((!pedited || spot_edited->visireti) && spot.visireti) {
                    saveToKeyfile(!pedited || spot_edited->expreti, "Locallab", "Expreti_" + index_str, spot.expreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexreti, "Locallab", "Complexreti_" + index_str, spot.complexreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->retinexMethod, "Locallab", "retinexMethod_" + index_str, spot.retinexMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->str, "Locallab", "Str_" + index_str, spot.str, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chrrt, "Locallab", "Chrrt_" + index_str, spot.chrrt, keyFile);
                    saveToKeyfile(!pedited || spot_edited->neigh, "Locallab", "Neigh_" + index_str, spot.neigh, keyFile);
                    saveToKeyfile(!pedited || spot_edited->vart, "Locallab", "Vart_" + index_str, spot.vart, keyFile);
                    saveToKeyfile(!pedited || spot_edited->offs, "Locallab", "Offs_" + index_str, spot.offs, keyFile);
                    saveToKeyfile(!pedited || spot_edited->dehaz, "Locallab", "Dehaz_" + index_str, spot.dehaz, keyFile);
                    saveToKeyfile(!pedited || spot_edited->depth, "Locallab", "Depth_" + index_str, spot.depth, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensih, "Locallab", "Sensih_" + index_str, spot.sensih, keyFile);
                    saveToKeyfile(!pedited || spot_edited->localTgaincurve, "Locallab", "TgainCurve_" + index_str, spot.localTgaincurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->localTtranscurve, "Locallab", "TtransCurve_" + index_str, spot.localTtranscurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->inversret, "Locallab", "Inversret_" + index_str, spot.inversret, keyFile);
                    saveToKeyfile(!pedited || spot_edited->equilret, "Locallab", "Equilret_" + index_str, spot.equilret, keyFile);
                    saveToKeyfile(!pedited || spot_edited->loglin, "Locallab", "Loglin_" + index_str, spot.loglin, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lumonly, "Locallab", "Lumonly_" + index_str, spot.lumonly, keyFile);
                    saveToKeyfile(!pedited || spot_edited->softradiusret, "Locallab", "Softradiusret_" + index_str, spot.softradiusret, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmaskreticurve, "Locallab", "CCmaskretiCurve_" + index_str, spot.CCmaskreticurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskreticurve, "Locallab", "LLmaskretiCurve_" + index_str, spot.LLmaskreticurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmaskreticurve, "Locallab", "HHmaskretiCurve_" + index_str, spot.HHmaskreticurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enaretiMask, "Locallab", "EnaretiMask_" + index_str, spot.enaretiMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enaretiMasktmap, "Locallab", "EnaretiMasktmap_" + index_str, spot.enaretiMasktmap, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskreti, "Locallab", "Blendmaskreti_" + index_str, spot.blendmaskreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmaskreti, "Locallab", "Radmaskreti_" + index_str, spot.radmaskreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromaskreti, "Locallab", "Chromaskreti_" + index_str, spot.chromaskreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammaskreti, "Locallab", "Gammaskreti_" + index_str, spot.gammaskreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomaskreti, "Locallab", "Slomaskreti_" + index_str, spot.slomaskreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmaskreti, "Locallab", "Lapmaskreti_" + index_str, spot.lapmaskreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->scalereti, "Locallab", "Scalereti_" + index_str, spot.scalereti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->darkness, "Locallab", "Darkness_" + index_str, spot.darkness, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lightnessreti, "Locallab", "Lightnessreti_" + index_str, spot.lightnessreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->limd, "Locallab", "Limd_" + index_str, spot.limd, keyFile);
                    saveToKeyfile(!pedited || spot_edited->cliptm, "Locallab", "Cliptm_" + index_str, spot.cliptm, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fftwreti, "Locallab", "Fftwreti_" + index_str, spot.fftwreti, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmaskreticurve, "Locallab", "LmaskretiCurve_" + index_str, spot.Lmaskreticurve, keyFile);
                }
                // Sharpening
                if ((!pedited || spot_edited->visisharp) && spot.visisharp) {
                    saveToKeyfile(!pedited || spot_edited->expsharp, "Locallab", "Expsharp_" + index_str, spot.expsharp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexsharp, "Locallab", "Complexsharp_" + index_str, spot.complexsharp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sharcontrast, "Locallab", "Sharcontrast_" + index_str, spot.sharcontrast, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sharradius, "Locallab", "Sharradius_" + index_str, spot.sharradius, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sharamount, "Locallab", "Sharamount_" + index_str, spot.sharamount, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shardamping, "Locallab", "Shardamping_" + index_str, spot.shardamping, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shariter, "Locallab", "Shariter_" + index_str, spot.shariter, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sharblur, "Locallab", "Sharblur_" + index_str, spot.sharblur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensisha, "Locallab", "Sensisha_" + index_str, spot.sensisha, keyFile);
                    saveToKeyfile(!pedited || spot_edited->inverssha, "Locallab", "Inverssha_" + index_str, spot.inverssha, keyFile);
                }
                // Local Contrast
                if ((!pedited || spot_edited->visicontrast) && spot.visicontrast) {
                    saveToKeyfile(!pedited || spot_edited->expcontrast, "Locallab", "Expcontrast_" + index_str, spot.expcontrast, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexcontrast, "Locallab", "Complexcontrast_" + index_str, spot.complexcontrast, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lcradius, "Locallab", "Lcradius_" + index_str, spot.lcradius, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lcamount, "Locallab", "Lcamount_" + index_str, spot.lcamount, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lcdarkness, "Locallab", "Lcdarkness_" + index_str, spot.lcdarkness, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lclightness, "Locallab", "Lclightness_" + index_str, spot.lclightness, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sigmalc, "Locallab", "Sigmalc_" + index_str, spot.sigmalc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->levelwav, "Locallab", "Levelwav_" + index_str, spot.levelwav, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residcont, "Locallab", "Residcont_" + index_str, spot.residcont, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residsha, "Locallab", "Residsha_" + index_str, spot.residsha, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residshathr, "Locallab", "Residshathr_" + index_str, spot.residshathr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residhi, "Locallab", "Residhi_" + index_str, spot.residhi, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residhithr, "Locallab", "Residhithr_" + index_str, spot.residhithr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residblur, "Locallab", "Residblur_" + index_str, spot.residblur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->levelblur, "Locallab", "Levelblur_" + index_str, spot.levelblur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sigmabl, "Locallab", "Sigmabl_" + index_str, spot.sigmabl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residchro, "Locallab", "Residchro_" + index_str, spot.residchro, keyFile);
                    saveToKeyfile(!pedited || spot_edited->residcomp, "Locallab", "Residcomp_" + index_str, spot.residcomp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sigma, "Locallab", "Sigma_" + index_str, spot.sigma, keyFile);
                    saveToKeyfile(!pedited || spot_edited->offset, "Locallab", "Offset_" + index_str, spot.offset, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sigmadr, "Locallab", "Sigmadr_" + index_str, spot.sigmadr, keyFile);
                    saveToKeyfile(!pedited || spot_edited->threswav, "Locallab", "Threswav_" + index_str, spot.threswav, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromalev, "Locallab", "Chromalev_" + index_str, spot.chromalev, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromablu, "Locallab", "Chromablu_" + index_str, spot.chromablu, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sigmadc, "Locallab", "sigmadc_" + index_str, spot.sigmadc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->deltad, "Locallab", "deltad_" + index_str, spot.deltad, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fatres, "Locallab", "Fatres_" + index_str, spot.fatres, keyFile);
                    saveToKeyfile(!pedited || spot_edited->clarilres, "Locallab", "ClariLres_" + index_str, spot.clarilres, keyFile);
                    saveToKeyfile(!pedited || spot_edited->claricres, "Locallab", "ClariCres_" + index_str, spot.claricres, keyFile);
                    saveToKeyfile(!pedited || spot_edited->clarisoft, "Locallab", "Clarisoft_" + index_str, spot.clarisoft, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sigmalc2, "Locallab", "Sigmalc2_" + index_str, spot.sigmalc2, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strwav, "Locallab", "Strwav_" + index_str, spot.strwav, keyFile);
                    saveToKeyfile(!pedited || spot_edited->angwav, "Locallab", "Angwav_" + index_str, spot.angwav, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strengthw, "Locallab", "Strengthw_" + index_str, spot.strengthw, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sigmaed, "Locallab", "Sigmaed_" + index_str, spot.sigmaed, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radiusw, "Locallab", "Radiusw_" + index_str, spot.radiusw, keyFile);
                    saveToKeyfile(!pedited || spot_edited->detailw, "Locallab", "Detailw_" + index_str, spot.detailw, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gradw, "Locallab", "Gradw_" + index_str, spot.gradw, keyFile);
                    saveToKeyfile(!pedited || spot_edited->tloww, "Locallab", "Tloww_" + index_str, spot.tloww, keyFile);
                    saveToKeyfile(!pedited || spot_edited->thigw, "Locallab", "Thigw_" + index_str, spot.thigw, keyFile);
                    saveToKeyfile(!pedited || spot_edited->edgw, "Locallab", "Edgw_" + index_str, spot.edgw, keyFile);
                    saveToKeyfile(!pedited || spot_edited->basew, "Locallab", "Basew_" + index_str, spot.basew, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensilc, "Locallab", "Sensilc_" + index_str, spot.sensilc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fftwlc, "Locallab", "Fftwlc_" + index_str, spot.fftwlc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blurlc, "Locallab", "Blurlc_" + index_str, spot.blurlc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->wavblur, "Locallab", "Wavblur_" + index_str, spot.wavblur, keyFile);
                    saveToKeyfile(!pedited || spot_edited->wavedg, "Locallab", "Wavedg_" + index_str, spot.wavedg, keyFile);
                    saveToKeyfile(!pedited || spot_edited->waveshow, "Locallab", "Waveshow_" + index_str, spot.waveshow, keyFile);
                    saveToKeyfile(!pedited || spot_edited->wavcont, "Locallab", "Wavcont_" + index_str, spot.wavcont, keyFile);
                    saveToKeyfile(!pedited || spot_edited->wavcomp, "Locallab", "Wavcomp_" + index_str, spot.wavcomp, keyFile);
                    saveToKeyfile(!pedited || spot_edited->wavgradl, "Locallab", "Wavgradl_" + index_str, spot.wavgradl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->wavcompre, "Locallab", "Wavcompre_" + index_str, spot.wavcompre, keyFile);
                    saveToKeyfile(!pedited || spot_edited->origlc, "Locallab", "Origlc_" + index_str, spot.origlc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->localcontMethod, "Locallab", "localcontMethod_" + index_str, spot.localcontMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->localedgMethod, "Locallab", "localedgMethod_" + index_str, spot.localedgMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->localneiMethod, "Locallab", "localneiMethod_" + index_str, spot.localneiMethod, keyFile);
                    saveToKeyfile(!pedited || spot_edited->locwavcurve, "Locallab", "LocwavCurve_" + index_str, spot.locwavcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->csthreshold, "Locallab", "CSThreshold_" + index_str, spot.csthreshold.toVector(), keyFile);
                    saveToKeyfile(!pedited || spot_edited->loclevwavcurve, "Locallab", "LoclevwavCurve_" + index_str, spot.loclevwavcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->locconwavcurve, "Locallab", "LocconwavCurve_" + index_str, spot.locconwavcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->loccompwavcurve, "Locallab", "LoccompwavCurve_" + index_str, spot.loccompwavcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->loccomprewavcurve, "Locallab", "LoccomprewavCurve_" + index_str, spot.loccomprewavcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->locedgwavcurve, "Locallab", "LocedgwavCurve_" + index_str, spot.locedgwavcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmasklccurve, "Locallab", "CCmasklcCurve_" + index_str, spot.CCmasklccurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmasklccurve, "Locallab", "LLmasklcCurve_" + index_str, spot.LLmasklccurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmasklccurve, "Locallab", "HHmasklcCurve_" + index_str, spot.HHmasklccurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enalcMask, "Locallab", "EnalcMask_" + index_str, spot.enalcMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmasklc, "Locallab", "Blendmasklc_" + index_str, spot.blendmasklc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmasklc, "Locallab", "Radmasklc_" + index_str, spot.radmasklc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromasklc, "Locallab", "Chromasklc_" + index_str, spot.chromasklc, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmasklccurve, "Locallab", "LmasklcCurve_" + index_str, spot.Lmasklccurve, keyFile);
                }
                // Contrast by detail levels
                if ((!pedited || spot_edited->visicbdl) && spot.visicbdl) {
                    saveToKeyfile(!pedited || spot_edited->expcbdl, "Locallab", "Expcbdl_" + index_str, spot.expcbdl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexcbdl, "Locallab", "Complexcbdl_" + index_str, spot.complexcbdl, keyFile);

                    for (int j = 0; j < 6; j++) {
                        saveToKeyfile(!pedited || spot_edited->mult[j], "Locallab", "Mult" + std::to_string(j) + "_" + index_str, spot.mult[j], keyFile);
                    }

                    saveToKeyfile(!pedited || spot_edited->chromacbdl, "Locallab", "Chromacbdl_" + index_str, spot.chromacbdl, keyFile);
                    saveToKeyfile(!pedited || spot_edited->threshold, "Locallab", "Threshold_" + index_str, spot.threshold, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensicb, "Locallab", "Sensicb_" + index_str, spot.sensicb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->clarityml, "Locallab", "Clarityml_" + index_str, spot.clarityml, keyFile);
                    saveToKeyfile(!pedited || spot_edited->contresid, "Locallab", "Contresid_" + index_str, spot.contresid, keyFile);
                    saveToKeyfile(!pedited || spot_edited->softradiuscb, "Locallab", "Softradiuscb_" + index_str, spot.softradiuscb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enacbMask, "Locallab", "EnacbMask_" + index_str, spot.enacbMask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmaskcbcurve, "Locallab", "CCmaskcbCurve_" + index_str, spot.CCmaskcbcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmaskcbcurve, "Locallab", "LLmaskcbCurve_" + index_str, spot.LLmaskcbcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmaskcbcurve, "Locallab", "HHmaskcbCurve_" + index_str, spot.HHmaskcbcurve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskcb, "Locallab", "Blendmaskcb_" + index_str, spot.blendmaskcb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmaskcb, "Locallab", "Radmaskcb_" + index_str, spot.radmaskcb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromaskcb, "Locallab", "Chromaskcb_" + index_str, spot.chromaskcb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammaskcb, "Locallab", "Gammaskcb_" + index_str, spot.gammaskcb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slomaskcb, "Locallab", "Slomaskcb_" + index_str, spot.slomaskcb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmaskcb, "Locallab", "Lapmaskcb_" + index_str, spot.lapmaskcb, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmaskcbcurve, "Locallab", "LmaskcbCurve_" + index_str, spot.Lmaskcbcurve, keyFile);
                }
                // Log encoding
                if ((!pedited || spot_edited->visilog) && spot.visilog) {
                    saveToKeyfile(!pedited || spot_edited->explog, "Locallab", "Explog_" + index_str, spot.explog, keyFile);
                    saveToKeyfile(!pedited || spot_edited->autocompute, "Locallab", "Autocompute_" + index_str, spot.autocompute, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sourceGray, "Locallab", "SourceGray_" + index_str, spot.sourceGray, keyFile);
                    saveToKeyfile(!pedited || spot_edited->targetGray, "Locallab", "TargetGray_" + index_str, spot.targetGray, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Autogray, "Locallab", "Autogray_" + index_str, spot.Autogray, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fullimage, "Locallab", "Fullimage_" + index_str, spot.fullimage, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blackEv, "Locallab", "BlackEv_" + index_str, spot.blackEv, keyFile);
                    saveToKeyfile(!pedited || spot_edited->whiteEv, "Locallab", "WhiteEv_" + index_str, spot.whiteEv, keyFile);
                    saveToKeyfile(!pedited || spot_edited->detail, "Locallab", "Detail_" + index_str, spot.detail, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensilog, "Locallab", "Sensilog_" + index_str, spot.sensilog, keyFile);
                    saveToKeyfile(!pedited || spot_edited->baselog, "Locallab", "Baselog_" + index_str, spot.baselog, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strlog, "Locallab", "Strlog_" + index_str, spot.strlog, keyFile);
                    saveToKeyfile(!pedited || spot_edited->anglog, "Locallab", "Anglog_" + index_str, spot.anglog, keyFile);
                }
                //mask
                if ((!pedited || spot_edited->visimask) && spot.visimask) {
                    saveToKeyfile(!pedited || spot_edited->expmask, "Locallab", "Expmask_" + index_str, spot.expmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->complexmask, "Locallab", "Complexmask_" + index_str, spot.complexmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->sensimask, "Locallab", "Sensimask_" + index_str, spot.sensimask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmask, "Locallab", "Blendmaskmask_" + index_str, spot.blendmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blendmaskab, "Locallab", "Blendmaskmaskab_" + index_str, spot.blendmaskab, keyFile);
                    saveToKeyfile(!pedited || spot_edited->softradiusmask, "Locallab", "Softradiusmask_" + index_str, spot.softradiusmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->enamask, "Locallab", "Enamask_" + index_str, spot.enamask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->fftmask, "Locallab", "Fftmask_" + index_str, spot.fftmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->blurmask, "Locallab", "Blurmask_" + index_str, spot.blurmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->contmask, "Locallab", "Contmask_" + index_str, spot.contmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->CCmask_curve, "Locallab", "CCmask_Curve_" + index_str, spot.CCmask_curve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmask_curve, "Locallab", "LLmask_Curve_" + index_str, spot.LLmask_curve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHmask_curve, "Locallab", "HHmask_Curve_" + index_str, spot.HHmask_curve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->strumaskmask, "Locallab", "Strumaskmask_" + index_str, spot.strumaskmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->toolmask, "Locallab", "Toolmask_" + index_str, spot.toolmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->radmask, "Locallab", "Radmask_" + index_str, spot.radmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->lapmask, "Locallab", "Lapmask_" + index_str, spot.lapmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->chromask, "Locallab", "Chromask_" + index_str, spot.chromask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->gammask, "Locallab", "Gammask_" + index_str, spot.gammask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->slopmask, "Locallab", "Slopmask_" + index_str, spot.slopmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->shadmask, "Locallab", "Shadmask_" + index_str, spot.shadmask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->str_mask, "Locallab", "Str_mask_" + index_str, spot.str_mask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->ang_mask, "Locallab", "Ang_mask_" + index_str, spot.ang_mask, keyFile);
                    saveToKeyfile(!pedited || spot_edited->HHhmask_curve, "Locallab", "HHhmask_Curve_" + index_str, spot.HHhmask_curve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->Lmask_curve, "Locallab", "Lmask_Curve_" + index_str, spot.Lmask_curve, keyFile);
                    saveToKeyfile(!pedited || spot_edited->LLmask_curvewav, "Locallab", "LLmask_Curvewav_" + index_str, spot.LLmask_curvewav, keyFile);
                    saveToKeyfile(!pedited || spot_edited->csthresholdmask, "Locallab", "CSThresholdmask_" + index_str, spot.csthresholdmask.toVector(), keyFile);
                }
            }
        }

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
        saveToKeyfile(!pedited || pedited->resize.allowUpscaling, "Resize", "AllowUpscaling", resize.allowUpscaling, keyFile);

// Post demosaic sharpening
        saveToKeyfile(!pedited || pedited->pdsharpening.enabled, "PostDemosaicSharpening", "Enabled", pdsharpening.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->pdsharpening.contrast, "PostDemosaicSharpening", "Contrast", pdsharpening.contrast, keyFile);
        saveToKeyfile(!pedited || pedited->pdsharpening.autoContrast, "PostDemosaicSharpening", "AutoContrast", pdsharpening.autoContrast, keyFile);
        saveToKeyfile(!pedited || pedited->pdsharpening.autoRadius, "PostDemosaicSharpening", "AutoRadius", pdsharpening.autoRadius, keyFile);
        saveToKeyfile(!pedited || pedited->pdsharpening.deconvradius, "PostDemosaicSharpening", "DeconvRadius", pdsharpening.deconvradius, keyFile);
        saveToKeyfile(!pedited || pedited->pdsharpening.deconvradiusOffset, "PostDemosaicSharpening", "DeconvRadiusOffset", pdsharpening.deconvradiusOffset, keyFile);
        saveToKeyfile(!pedited || pedited->pdsharpening.deconvitercheck, "PostDemosaicSharpening", "DeconvIterCheck", pdsharpening.deconvitercheck, keyFile);
        saveToKeyfile(!pedited || pedited->pdsharpening.deconviter, "PostDemosaicSharpening", "DeconvIterations", pdsharpening.deconviter, keyFile);

// Post resize sharpening
        saveToKeyfile(!pedited || pedited->prsharpening.enabled, "PostResizeSharpening", "Enabled", prsharpening.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->prsharpening.contrast, "PostResizeSharpening", "Contrast", prsharpening.contrast, keyFile);
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
        saveToKeyfile(!pedited || pedited->icm.inputProfile, "Color Management", "InputProfile", relativePathIfInside(fname, fnameAbsolute, icm.inputProfile), keyFile);
        saveToKeyfile(!pedited || pedited->icm.toneCurve, "Color Management", "ToneCurve", icm.toneCurve, keyFile);
        saveToKeyfile(!pedited || pedited->icm.applyLookTable, "Color Management", "ApplyLookTable", icm.applyLookTable, keyFile);
        saveToKeyfile(!pedited || pedited->icm.applyBaselineExposureOffset, "Color Management", "ApplyBaselineExposureOffset", icm.applyBaselineExposureOffset, keyFile);
        saveToKeyfile(!pedited || pedited->icm.applyHueSatMap, "Color Management", "ApplyHueSatMap", icm.applyHueSatMap, keyFile);
        saveToKeyfile(!pedited || pedited->icm.dcpIlluminant, "Color Management", "DCPIlluminant", icm.dcpIlluminant, keyFile);
        saveToKeyfile(!pedited || pedited->icm.workingProfile, "Color Management", "WorkingProfile", icm.workingProfile, keyFile);
        saveToKeyfile(!pedited || pedited->icm.workingTRC, "Color Management", "WorkingTRC", icm.workingTRC, keyFile);
        saveToKeyfile(!pedited || pedited->icm.workingTRCGamma, "Color Management", "WorkingTRCGamma", icm.workingTRCGamma, keyFile);
        saveToKeyfile(!pedited || pedited->icm.workingTRCSlope, "Color Management", "WorkingTRCSlope", icm.workingTRCSlope, keyFile);
        saveToKeyfile(!pedited || pedited->icm.outputProfile, "Color Management", "OutputProfile", icm.outputProfile, keyFile);
        saveToKeyfile(
            !pedited || pedited->icm.outputIntent,
            "Color Management",
            "OutputProfileIntent",
            {
                {RI_PERCEPTUAL, "Perceptual"},
                {RI_RELATIVE, "Relative"},
                {RI_SATURATION, "Saturation"},
                {RI_ABSOLUTE, "Absolute"}

            },
            icm.outputIntent,
            keyFile
        );
        saveToKeyfile(!pedited || pedited->icm.outputBPC, "Color Management", "OutputBPC", icm.outputBPC, keyFile);

// Wavelet
        saveToKeyfile(!pedited || pedited->wavelet.enabled, "Wavelet", "Enabled", wavelet.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.strength, "Wavelet", "Strength", wavelet.strength, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.balance, "Wavelet", "Balance", wavelet.balance, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.sigmafin, "Wavelet", "Sigmafin", wavelet.sigmafin, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.sigmaton, "Wavelet", "Sigmaton", wavelet.sigmaton, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.sigmacol, "Wavelet", "Sigmacol", wavelet.sigmacol, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.sigmadir, "Wavelet", "Sigmadir", wavelet.sigmadir, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.rangeab, "Wavelet", "Rangeab", wavelet.rangeab, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.protab, "Wavelet", "Protab", wavelet.protab, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.iter, "Wavelet", "Iter", wavelet.iter, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.thres, "Wavelet", "MaxLev", wavelet.thres, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.Tilesmethod, "Wavelet", "TilesMethod", wavelet.Tilesmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.complexmethod, "Wavelet", "complexMethod", wavelet.complexmethod, keyFile);
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
        saveToKeyfile(!pedited || pedited->wavelet.ballum, "Wavelet", "Ballum", wavelet.ballum, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.balchrom, "Wavelet", "Balchrom", wavelet.balchrom, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.chromfi, "Wavelet", "Chromfine", wavelet.chromfi, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.chromco, "Wavelet", "Chromcoarse", wavelet.chromco, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.mergeL, "Wavelet", "MergeL", wavelet.mergeL, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.mergeC, "Wavelet", "MergeC", wavelet.mergeC, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.softrad, "Wavelet", "Softrad", wavelet.softrad, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.softradend, "Wavelet", "Softradend", wavelet.softradend, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expcontrast, "Wavelet", "Expcontrast", wavelet.expcontrast, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expchroma, "Wavelet", "Expchroma", wavelet.expchroma, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expedge, "Wavelet", "Expedge", wavelet.expedge, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expbl, "Wavelet", "expbl", wavelet.expbl, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expresid, "Wavelet", "Expresid", wavelet.expresid, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expfinal, "Wavelet", "Expfinal", wavelet.expfinal, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.exptoning, "Wavelet", "Exptoning", wavelet.exptoning, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expnoise, "Wavelet", "Expnoise", wavelet.expnoise, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.expclari, "Wavelet", "Expclari", wavelet.expclari, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.labgridALow, "Wavelet", "LabGridALow", wavelet.labgridALow, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.labgridBLow, "Wavelet", "LabGridBLow", wavelet.labgridBLow, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.labgridAHigh, "Wavelet", "LabGridAHigh", wavelet.labgridAHigh, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.labgridBHigh, "Wavelet", "LabGridBHigh", wavelet.labgridBHigh, keyFile);

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
        saveToKeyfile(!pedited || pedited->wavelet.ushamethod, "Wavelet", "Ushamethod", wavelet.ushamethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.CHSLmethod, "Wavelet", "CHSLromaMethod", wavelet.CHSLmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.EDmethod, "Wavelet", "EDMethod", wavelet.EDmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.NPmethod, "Wavelet", "NPMethod", wavelet.NPmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.BAmethod, "Wavelet", "BAMethod", wavelet.BAmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.TMmethod, "Wavelet", "TMMethod", wavelet.TMmethod, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.chro, "Wavelet", "ChromaLink", wavelet.chro, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.ccwcurve, "Wavelet", "ContrastCurve", wavelet.ccwcurve, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.blcurve, "Wavelet", "blcurve", wavelet.blcurve, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.pastlev, "Wavelet", "Pastlev", wavelet.pastlev.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.satlev, "Wavelet", "Satlev", wavelet.satlev.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.opacityCurveRG, "Wavelet", "OpacityCurveRG", wavelet.opacityCurveRG, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.opacityCurveSH, "Wavelet", "Levalshc", wavelet.opacityCurveSH, keyFile);
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
        saveToKeyfile(!pedited || pedited->wavelet.chrwav, "Wavelet", "chrwav", wavelet.chrwav, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.bluwav, "Wavelet", "bluwav", wavelet.bluwav, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.hueskin, "Wavelet", "Hueskin", wavelet.hueskin.toVector(), keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgrad, "Wavelet", "Edgrad", wavelet.edgrad, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgeffect, "Wavelet", "Edgeffect", wavelet.edgeffect, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgval, "Wavelet", "Edgval", wavelet.edgval, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgthresh, "Wavelet", "ThrEdg", wavelet.edgthresh, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.avoid, "Wavelet", "AvoidColorShift", wavelet.avoid, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.showmask, "Wavelet", "Showmask", wavelet.showmask, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.oldsh, "Wavelet", "Oldsh", wavelet.oldsh, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.tmr, "Wavelet", "TMr", wavelet.tmr, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.sigma, "Wavelet", "Sigma", wavelet.sigma, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.offset, "Wavelet", "Offset", wavelet.offset, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.lowthr, "Wavelet", "Lowthr", wavelet.lowthr, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.rescon, "Wavelet", "ResidualcontShadow", wavelet.rescon, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.resconH, "Wavelet", "ResidualcontHighlight", wavelet.resconH, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.thr, "Wavelet", "ThresholdResidShadow", wavelet.thr, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.thrH, "Wavelet", "ThresholdResidHighLight", wavelet.thrH, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.radius, "Wavelet", "Residualradius", wavelet.radius, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.reschro, "Wavelet", "Residualchroma", wavelet.reschro, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.resblur, "Wavelet", "Residualblur", wavelet.resblur, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.resblurc, "Wavelet", "Residualblurc", wavelet.resblurc, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.tmrs, "Wavelet", "ResidualTM", wavelet.tmrs, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.edgs, "Wavelet", "ResidualEDGS", wavelet.edgs, keyFile);
        saveToKeyfile(!pedited || pedited->wavelet.scale, "Wavelet", "ResidualSCALE", wavelet.scale, keyFile);
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

// Soft Light
        saveToKeyfile(!pedited || pedited->softlight.enabled, "SoftLight", "Enabled", softlight.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->softlight.strength, "SoftLight", "Strength", softlight.strength, keyFile);

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
        if (!pedited || pedited->colorToning.labregions) {
            for (size_t j = 0; j < colorToning.labregions.size(); ++j) {
                std::string n = std::to_string(j+1);
                auto &l = colorToning.labregions[j];
                putToKeyfile("ColorToning", Glib::ustring("LabRegionA_") + n, l.a, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionB_") + n, l.b, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionSaturation_") + n, l.saturation, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionSlope_") + n, l.slope, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionOffset_") + n, l.offset, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionPower_") + n, l.power, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionHueMask_") + n, l.hueMask, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionChromaticityMask_") + n, l.chromaticityMask, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionLightnessMask_") + n, l.lightnessMask, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionMaskBlur_") + n, l.maskBlur, keyFile);
                putToKeyfile("ColorToning", Glib::ustring("LabRegionChannel_") + n, l.channel, keyFile);
            }
        }
        saveToKeyfile(!pedited || pedited->colorToning.labregionsShowMask, "ColorToning", "LabRegionsShowMask", colorToning.labregionsShowMask, keyFile);

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
        saveToKeyfile(!pedited || pedited->raw.ca_avoidcolourshift, "RAW", "CAAvoidColourshift", raw.ca_avoidcolourshift, keyFile);
        saveToKeyfile(!pedited || pedited->raw.caautoiterations, "RAW", "CAAutoIterations", raw.caautoiterations, keyFile);
        saveToKeyfile(!pedited || pedited->raw.cared, "RAW", "CARed", raw.cared, keyFile);
        saveToKeyfile(!pedited || pedited->raw.cablue, "RAW", "CABlue", raw.cablue, keyFile);
        saveToKeyfile(!pedited || pedited->raw.hotPixelFilter, "RAW", "HotPixelFilter", raw.hotPixelFilter, keyFile);
        saveToKeyfile(!pedited || pedited->raw.deadPixelFilter, "RAW", "DeadPixelFilter", raw.deadPixelFilter, keyFile);
        saveToKeyfile(!pedited || pedited->raw.hotdeadpix_thresh, "RAW", "HotDeadPixelThresh", raw.hotdeadpix_thresh, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.method, "RAW Bayer", "Method", raw.bayersensor.method, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.border, "RAW Bayer", "Border", raw.bayersensor.border, keyFile);
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
        saveToKeyfile(!pedited || pedited->raw.bayersensor.dualDemosaicAutoContrast, "RAW Bayer", "DualDemosaicAutoContrast", raw.bayersensor.dualDemosaicAutoContrast, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.dualDemosaicContrast, "RAW Bayer", "DualDemosaicContrast", raw.bayersensor.dualDemosaicContrast, keyFile);
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
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftEqualBright, "RAW Bayer", "pixelShiftEqualBright", raw.bayersensor.pixelShiftEqualBright, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftEqualBrightChannel, "RAW Bayer", "pixelShiftEqualBrightChannel", raw.bayersensor.pixelShiftEqualBrightChannel, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftNonGreenCross, "RAW Bayer", "pixelShiftNonGreenCross", raw.bayersensor.pixelShiftNonGreenCross, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pixelShiftDemosaicMethod, "RAW Bayer", "pixelShiftDemosaicMethod", raw.bayersensor.pixelShiftDemosaicMethod, keyFile);
        saveToKeyfile(!pedited || pedited->raw.bayersensor.pdafLinesFilter, "RAW Bayer", "PDAFLinesFilter", raw.bayersensor.pdafLinesFilter, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.method, "RAW X-Trans", "Method", raw.xtranssensor.method, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.dualDemosaicAutoContrast, "RAW X-Trans", "DualDemosaicAutoContrast", raw.xtranssensor.dualDemosaicAutoContrast, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.dualDemosaicContrast, "RAW X-Trans", "DualDemosaicContrast", raw.xtranssensor.dualDemosaicContrast, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.border, "RAW X-Trans", "Border", raw.xtranssensor.border, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.ccSteps, "RAW X-Trans", "CcSteps", raw.xtranssensor.ccSteps, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.exBlackRed, "RAW X-Trans", "PreBlackRed", raw.xtranssensor.blackred, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.exBlackGreen, "RAW X-Trans", "PreBlackGreen", raw.xtranssensor.blackgreen, keyFile);
        saveToKeyfile(!pedited || pedited->raw.xtranssensor.exBlackBlue, "RAW X-Trans", "PreBlackBlue", raw.xtranssensor.blackblue, keyFile);

// Raw exposition
        saveToKeyfile(!pedited || pedited->raw.exPos, "RAW", "PreExposure", raw.expos, keyFile);

// MetaData
        saveToKeyfile(!pedited || pedited->metadata.mode, "MetaData", "Mode", metadata.mode, keyFile);

// Film negative
        saveToKeyfile(!pedited || pedited->filmNegative.enabled, "Film Negative", "Enabled", filmNegative.enabled, keyFile);
        saveToKeyfile(!pedited || pedited->filmNegative.redRatio, "Film Negative", "RedRatio", filmNegative.redRatio, keyFile);
        saveToKeyfile(!pedited || pedited->filmNegative.greenExp, "Film Negative", "GreenExponent", filmNegative.greenExp, keyFile);
        saveToKeyfile(!pedited || pedited->filmNegative.blueRatio, "Film Negative", "BlueRatio", filmNegative.blueRatio, keyFile);
        saveToKeyfile(!pedited || pedited->filmNegative.baseValues, "Film Negative", "RedBase", filmNegative.redBase, keyFile);
        saveToKeyfile(!pedited || pedited->filmNegative.baseValues, "Film Negative", "GreenBase", filmNegative.greenBase, keyFile);
        saveToKeyfile(!pedited || pedited->filmNegative.baseValues, "Film Negative", "BlueBase", filmNegative.blueBase, keyFile);

// Preprocess WB
        saveToKeyfile(!pedited || pedited->raw.preprocessWB.mode, "RAW Preprocess WB", "Mode", toUnderlying(raw.preprocessWB.mode), keyFile);

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

            const std::map<std::string, ToneCurveMode> tc_mapping = {
                {"Standard", ToneCurveMode::STD},
                {"FilmLike", ToneCurveMode::FILMLIKE},
                {"SatAndValueBlending", ToneCurveMode::SATANDVALBLENDING},
                {"WeightedStd", ToneCurveMode::WEIGHTEDSTD},
                {"Luminance", ToneCurveMode::LUMINANCE},
                {"Perceptual", ToneCurveMode::PERCEPTUAL}
            };

            assignFromKeyfile(keyFile, "Exposure", "CurveMode", pedited, tc_mapping, toneCurve.curveMode, pedited->toneCurve.curveMode);
            assignFromKeyfile(keyFile, "Exposure", "CurveMode2", pedited, tc_mapping, toneCurve.curveMode2, pedited->toneCurve.curveMode2);

            if (ppVersion > 200) {
                assignFromKeyfile(keyFile, "Exposure", "Curve", pedited, toneCurve.curve, pedited->toneCurve.curve);
                assignFromKeyfile(keyFile, "Exposure", "Curve2", pedited, toneCurve.curve2, pedited->toneCurve.curve2);
            }

            assignFromKeyfile(keyFile, "Exposure", "HistogramMatching", pedited, toneCurve.histmatching, pedited->toneCurve.histmatching);
            if (ppVersion < 340) {
                toneCurve.fromHistMatching = false;
                if (pedited) {
                    pedited->toneCurve.fromHistMatching = true;
                }
            } else {
                assignFromKeyfile(keyFile, "Exposure", "CurveFromHistogramMatching", pedited, toneCurve.fromHistMatching, pedited->toneCurve.fromHistMatching);
            }
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
                if (ppVersion < 338) {
                    for (int i = 0; i < 3; ++i) {
                        chmixer.red[i] *= 10;
                        chmixer.green[i] *= 10;
                        chmixer.blue[i] *= 10;
                    }
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
                pedited,
                {
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
                pedited,
                {
                    {"Standard", BlackWhiteParams::TcMode::STD_BW},
                    {"WeightedStd", BlackWhiteParams::TcMode::WEIGHTEDSTD_BW}
                },
                blackwhite.afterCurveMode,
                pedited->blackwhite.afterCurveMode
            );
        }

        if (keyFile.has_group("Retinex")) {
            assignFromKeyfile(keyFile, "Retinex", "Median", pedited, retinex.medianmap, pedited->retinex.medianmap);
            assignFromKeyfile(keyFile, "Retinex", "complexMethod", pedited, retinex.complexmethod, pedited->retinex.complexmethod);
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

            if (ppVersion >= 334) {
                assignFromKeyfile(keyFile, "Sharpening", "Contrast", pedited, sharpening.contrast, pedited->sharpening.contrast);
            } else {
                sharpening.contrast = 0;

                if (pedited) {
                    pedited->sharpening.contrast = true;
                }
            }

            assignFromKeyfile(keyFile, "Sharpening", "Radius", pedited, sharpening.radius, pedited->sharpening.radius);
            assignFromKeyfile(keyFile, "Sharpening", "BlurRadius", pedited, sharpening.blurradius, pedited->sharpening.blurradius);
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

            if (ppVersion >= 334) {
                assignFromKeyfile(keyFile, "SharpenMicro", "Contrast", pedited, sharpenMicro.contrast, pedited->sharpenMicro.contrast);
            } else {
                sharpenMicro.contrast = 0;

                if (pedited) {
                    pedited->sharpenMicro.contrast = true;
                }
            }
            if (ppVersion >= 346) {
                assignFromKeyfile(keyFile, "SharpenMicro", "Uniformity", pedited, sharpenMicro.uniformity, pedited->sharpenMicro.uniformity);
            } else {
                double temp = 50.0;
                assignFromKeyfile(keyFile, "SharpenMicro", "Uniformity", pedited, temp, pedited->sharpenMicro.uniformity);
                sharpenMicro.uniformity = temp / 10;
            }
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
            if (wb.method == "Auto") {
                wb.method = "autold";
            }
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
            assignFromKeyfile(keyFile, "Color appearance", "Illum", pedited, colorappearance.illum, pedited->colorappearance.illum);
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
            assignFromKeyfile(keyFile, "Color appearance", "Autotempout", pedited, colorappearance.autotempout, pedited->colorappearance.autotempout);
            assignFromKeyfile(keyFile, "Color appearance", "Greenout", pedited, colorappearance.greenout, pedited->colorappearance.greenout);
            assignFromKeyfile(keyFile, "Color appearance", "Tempsc", pedited, colorappearance.tempsc, pedited->colorappearance.tempsc);
            assignFromKeyfile(keyFile, "Color appearance", "Greensc", pedited, colorappearance.greensc, pedited->colorappearance.greensc);
            assignFromKeyfile(keyFile, "Color appearance", "Ybout", pedited, colorappearance.ybout, pedited->colorappearance.ybout);
            assignFromKeyfile(keyFile, "Color appearance", "Datacie", pedited, colorappearance.datacie, pedited->colorappearance.datacie);
            assignFromKeyfile(keyFile, "Color appearance", "Tonecie", pedited, colorappearance.tonecie, pedited->colorappearance.tonecie);
            assignFromKeyfile(keyFile, "Color appearance", "Presetcat02", pedited, colorappearance.presetcat02, pedited->colorappearance.presetcat02);

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
                pedited,
                {
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

        if (keyFile.has_group("Shadows & Highlights") && ppVersion >= 333) {
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Enabled", pedited, sh.enabled, pedited->sh.enabled);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Highlights", pedited, sh.highlights, pedited->sh.highlights);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "HighlightTonalWidth", pedited, sh.htonalwidth, pedited->sh.htonalwidth);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Shadows", pedited, sh.shadows, pedited->sh.shadows);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "ShadowTonalWidth", pedited, sh.stonalwidth, pedited->sh.stonalwidth);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Radius", pedited, sh.radius, pedited->sh.radius);
            if (ppVersion >= 344) {
                assignFromKeyfile(keyFile, "Shadows & Highlights", "Lab", pedited, sh.lab, pedited->sh.lab);
            } else {
                sh.lab = true;
            }

            if (keyFile.has_key("Shadows & Highlights", "LocalContrast") && ppVersion < 329) {
                int lc = keyFile.get_integer("Shadows & Highlights", "LocalContrast");
                localContrast.amount = float(lc) / 30.f;

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
            if (keyFile.has_key("Common Properties for Transformations", "Method")) {
                assignFromKeyfile(keyFile, "Common Properties for Transformations", "Method", pedited, commonTrans.method, pedited->commonTrans.method);
            } else {
                commonTrans.method = "lin";
            }
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
            assignFromKeyfile(keyFile, "Perspective", "Method", pedited, perspective.method, pedited->perspective.method);
            assignFromKeyfile(keyFile, "Perspective", "Horizontal", pedited, perspective.horizontal, pedited->perspective.horizontal);
            assignFromKeyfile(keyFile, "Perspective", "Vertical", pedited, perspective.vertical, pedited->perspective.vertical);
            assignFromKeyfile(keyFile, "Perspective", "CameraShiftHorizontal", pedited, perspective.camera_shift_horiz, pedited->perspective.camera_shift_horiz);
            assignFromKeyfile(keyFile, "Perspective", "CameraShiftVertical", pedited, perspective.camera_shift_vert, pedited->perspective.camera_shift_vert);
            assignFromKeyfile(keyFile, "Perspective", "CameraPitch", pedited, perspective.camera_pitch, pedited->perspective.camera_pitch);
            assignFromKeyfile(keyFile, "Perspective", "CameraRoll", pedited, perspective.camera_roll, pedited->perspective.camera_roll);
            assignFromKeyfile(keyFile, "Perspective", "CameraCropFactor", pedited, perspective.camera_crop_factor, pedited->perspective.camera_crop_factor);
            assignFromKeyfile(keyFile, "Perspective", "CameraFocalLength", pedited, perspective.camera_focal_length, pedited->perspective.camera_focal_length);
            assignFromKeyfile(keyFile, "Perspective", "CameraYaw", pedited, perspective.camera_yaw, pedited->perspective.camera_yaw);
            assignFromKeyfile(keyFile, "Perspective", "ProjectionPitch", pedited, perspective.projection_pitch, pedited->perspective.projection_pitch);
            assignFromKeyfile(keyFile, "Perspective", "ProjectionRotate", pedited, perspective.projection_rotate, pedited->perspective.projection_rotate);
            assignFromKeyfile(keyFile, "Perspective", "ProjectionShiftHorizontal", pedited, perspective.projection_shift_horiz, pedited->perspective.projection_shift_horiz);
            assignFromKeyfile(keyFile, "Perspective", "ProjectionShiftVertical", pedited, perspective.projection_shift_vert, pedited->perspective.projection_shift_vert);
            assignFromKeyfile(keyFile, "Perspective", "ProjectionYaw", pedited, perspective.projection_yaw, pedited->perspective.projection_yaw);
        }

        if (keyFile.has_group("Gradient")) {
            assignFromKeyfile(keyFile, "Gradient", "Enabled", pedited, gradient.enabled, pedited->gradient.enabled);
            assignFromKeyfile(keyFile, "Gradient", "Degree", pedited, gradient.degree, pedited->gradient.degree);
            assignFromKeyfile(keyFile, "Gradient", "Feather", pedited, gradient.feather, pedited->gradient.feather);
            assignFromKeyfile(keyFile, "Gradient", "Strength", pedited, gradient.strength, pedited->gradient.strength);
            assignFromKeyfile(keyFile, "Gradient", "CenterX", pedited, gradient.centerX, pedited->gradient.centerX);
            assignFromKeyfile(keyFile, "Gradient", "CenterY", pedited, gradient.centerY, pedited->gradient.centerY);
        }

        if (keyFile.has_group("Locallab")) {
            assignFromKeyfile(keyFile, "Locallab", "Enabled", pedited, locallab.enabled, pedited->locallab.enabled);
            assignFromKeyfile(keyFile, "Locallab", "Selspot", pedited, locallab.selspot, pedited->locallab.selspot);

            Glib::ustring ppName;
            bool peName;
            int i = 0;

            while (assignFromKeyfile(keyFile, "Locallab", "Name_" + std::to_string(i), pedited, ppName, peName)) {
                const std::string index_str = std::to_string(i);

                // Create new LocallabSpot and LocallabParamsEdited
                LocallabParams::LocallabSpot spot;
                spot.name = ppName;
                LocallabParamsEdited::LocallabSpotEdited spotEdited(false);
                spotEdited.name = peName;

                // Control spot settings
                assignFromKeyfile(keyFile, "Locallab", "Isvisible_" + index_str, pedited, spot.isvisible, spotEdited.isvisible);
                assignFromKeyfile(keyFile, "Locallab", "PrevMethod_" + index_str, pedited, spot.prevMethod, spotEdited.prevMethod);
                assignFromKeyfile(keyFile, "Locallab", "Shape_" + index_str, pedited, spot.shape, spotEdited.shape);
                assignFromKeyfile(keyFile, "Locallab", "SpotMethod_" + index_str, pedited, spot.spotMethod, spotEdited.spotMethod);
                assignFromKeyfile(keyFile, "Locallab", "wavMethod_" + index_str, pedited, spot.wavMethod, spotEdited.wavMethod);
                assignFromKeyfile(keyFile, "Locallab", "SensiExclu_" + index_str, pedited, spot.sensiexclu, spotEdited.sensiexclu);
                assignFromKeyfile(keyFile, "Locallab", "StructExclu_" + index_str, pedited, spot.structexclu, spotEdited.structexclu);
                assignFromKeyfile(keyFile, "Locallab", "Struc_" + index_str, pedited, spot.struc, spotEdited.struc);
                assignFromKeyfile(keyFile, "Locallab", "ShapeMethod_" + index_str, pedited, spot.shapeMethod, spotEdited.shapeMethod);
                assignFromKeyfile(keyFile, "Locallab", "Loc_" + index_str, pedited, spot.loc, spotEdited.loc);
                assignFromKeyfile(keyFile, "Locallab", "CenterX_" + index_str, pedited, spot.centerX, spotEdited.centerX);
                assignFromKeyfile(keyFile, "Locallab", "CenterY_" + index_str, pedited, spot.centerY, spotEdited.centerY);
                assignFromKeyfile(keyFile, "Locallab", "Circrad_" + index_str, pedited, spot.circrad, spotEdited.circrad);
                assignFromKeyfile(keyFile, "Locallab", "QualityMethod_" + index_str, pedited, spot.qualityMethod, spotEdited.qualityMethod);
                assignFromKeyfile(keyFile, "Locallab", "ComplexMethod_" + index_str, pedited, spot.complexMethod, spotEdited.complexMethod);
                assignFromKeyfile(keyFile, "Locallab", "Transit_" + index_str, pedited, spot.transit, spotEdited.transit);
                assignFromKeyfile(keyFile, "Locallab", "Feather_" + index_str, pedited, spot.feather, spotEdited.feather);
                assignFromKeyfile(keyFile, "Locallab", "Thresh_" + index_str, pedited, spot.thresh, spotEdited.thresh);
                assignFromKeyfile(keyFile, "Locallab", "Iter_" + index_str, pedited, spot.iter, spotEdited.iter);
                assignFromKeyfile(keyFile, "Locallab", "Balan_" + index_str, pedited, spot.balan, spotEdited.balan);
                assignFromKeyfile(keyFile, "Locallab", "Balanh_" + index_str, pedited, spot.balanh, spotEdited.balanh);
                assignFromKeyfile(keyFile, "Locallab", "Colorde_" + index_str, pedited, spot.colorde, spotEdited.colorde);
                assignFromKeyfile(keyFile, "Locallab", "Colorscope_" + index_str, pedited, spot.colorscope, spotEdited.colorscope);
                assignFromKeyfile(keyFile, "Locallab", "Transitweak_" + index_str, pedited, spot.transitweak, spotEdited.transitweak);
                assignFromKeyfile(keyFile, "Locallab", "Transitgrad_" + index_str, pedited, spot.transitgrad, spotEdited.transitgrad);
                assignFromKeyfile(keyFile, "Locallab", "Activ_" + index_str, pedited, spot.activ, spotEdited.activ);
                assignFromKeyfile(keyFile, "Locallab", "Avoid_" + index_str, pedited, spot.avoid, spotEdited.avoid);
                assignFromKeyfile(keyFile, "Locallab", "Blwh_" + index_str, pedited, spot.blwh, spotEdited.blwh);
                assignFromKeyfile(keyFile, "Locallab", "Recurs_" + index_str, pedited, spot.recurs, spotEdited.recurs);
                assignFromKeyfile(keyFile, "Locallab", "Laplac_" + index_str, pedited, spot.laplac, spotEdited.laplac);
                assignFromKeyfile(keyFile, "Locallab", "Deltae_" + index_str, pedited, spot.deltae, spotEdited.deltae);
                assignFromKeyfile(keyFile, "Locallab", "Shortc_" + index_str, pedited, spot.shortc, spotEdited.shortc);
                assignFromKeyfile(keyFile, "Locallab", "Savrest_" + index_str, pedited, spot.savrest, spotEdited.savrest);
                assignFromKeyfile(keyFile, "Locallab", "Scopemask_" + index_str, pedited, spot.scopemask, spotEdited.scopemask);
                assignFromKeyfile(keyFile, "Locallab", "Lumask_" + index_str, pedited, spot.lumask, spotEdited.lumask);
                // Color & Light
                spot.visicolor = assignFromKeyfile(keyFile, "Locallab", "Expcolor_" + index_str, pedited, spot.expcolor, spotEdited.expcolor);

                if (spot.visicolor) {
                    spotEdited.visicolor = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexcolor_" + index_str, pedited, spot.complexcolor, spotEdited.complexcolor);
                assignFromKeyfile(keyFile, "Locallab", "Curvactiv_" + index_str, pedited, spot.curvactiv, spotEdited.curvactiv);
                assignFromKeyfile(keyFile, "Locallab", "Lightness_" + index_str, pedited, spot.lightness, spotEdited.lightness);
                assignFromKeyfile(keyFile, "Locallab", "Contrast_" + index_str, pedited, spot.contrast, spotEdited.contrast);
                assignFromKeyfile(keyFile, "Locallab", "Chroma_" + index_str, pedited, spot.chroma, spotEdited.chroma);
                assignFromKeyfile(keyFile, "Locallab", "labgridALow_" + index_str, pedited, spot.labgridALow, spotEdited.labgridALow);
                assignFromKeyfile(keyFile, "Locallab", "labgridBLow_" + index_str, pedited, spot.labgridBLow, spotEdited.labgridBLow);
                assignFromKeyfile(keyFile, "Locallab", "labgridAHigh_" + index_str, pedited, spot.labgridAHigh, spotEdited.labgridAHigh);
                assignFromKeyfile(keyFile, "Locallab", "labgridBHigh_" + index_str, pedited, spot.labgridBHigh, spotEdited.labgridBHigh);
                assignFromKeyfile(keyFile, "Locallab", "labgridALowmerg_" + index_str, pedited, spot.labgridALowmerg, spotEdited.labgridALowmerg);
                assignFromKeyfile(keyFile, "Locallab", "labgridBLowmerg_" + index_str, pedited, spot.labgridBLowmerg, spotEdited.labgridBLowmerg);
                assignFromKeyfile(keyFile, "Locallab", "labgridAHighmerg_" + index_str, pedited, spot.labgridAHighmerg, spotEdited.labgridAHighmerg);
                assignFromKeyfile(keyFile, "Locallab", "labgridBHighmerg_" + index_str, pedited, spot.labgridBHighmerg, spotEdited.labgridBHighmerg);
                assignFromKeyfile(keyFile, "Locallab", "Strengthgrid_" + index_str, pedited, spot.strengthgrid, spotEdited.strengthgrid);
                assignFromKeyfile(keyFile, "Locallab", "Sensi_" + index_str, pedited, spot.sensi, spotEdited.sensi);
                assignFromKeyfile(keyFile, "Locallab", "Structcol_" + index_str, pedited, spot.structcol, spotEdited.structcol);
                assignFromKeyfile(keyFile, "Locallab", "Strcol_" + index_str, pedited, spot.strcol, spotEdited.strcol);
                assignFromKeyfile(keyFile, "Locallab", "Strcolab_" + index_str, pedited, spot.strcolab, spotEdited.strcolab);
                assignFromKeyfile(keyFile, "Locallab", "Strcolh_" + index_str, pedited, spot.strcolh, spotEdited.strcolh);
                assignFromKeyfile(keyFile, "Locallab", "Angcol_" + index_str, pedited, spot.angcol, spotEdited.angcol);
                assignFromKeyfile(keyFile, "Locallab", "Blurcolde_" + index_str, pedited, spot.blurcolde, spotEdited.blurcolde);
                assignFromKeyfile(keyFile, "Locallab", "Blurcol_" + index_str, pedited, spot.blurcol, spotEdited.blurcol);
                assignFromKeyfile(keyFile, "Locallab", "Contcol_" + index_str, pedited, spot.contcol, spotEdited.contcol);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskcol_" + index_str, pedited, spot.blendmaskcol, spotEdited.blendmaskcol);
                assignFromKeyfile(keyFile, "Locallab", "Radmaskcol_" + index_str, pedited, spot.radmaskcol, spotEdited.radmaskcol);
                assignFromKeyfile(keyFile, "Locallab", "Chromaskcol_" + index_str, pedited, spot.chromaskcol, spotEdited.chromaskcol);
                assignFromKeyfile(keyFile, "Locallab", "Gammaskcol_" + index_str, pedited, spot.gammaskcol, spotEdited.gammaskcol);
                assignFromKeyfile(keyFile, "Locallab", "Slomaskcol_" + index_str, pedited, spot.slomaskcol, spotEdited.slomaskcol);
                assignFromKeyfile(keyFile, "Locallab", "shadmaskcol_" + index_str, pedited, spot.shadmaskcol, spotEdited.shadmaskcol);
                assignFromKeyfile(keyFile, "Locallab", "strumaskcol_" + index_str, pedited, spot.strumaskcol, spotEdited.strumaskcol);
                assignFromKeyfile(keyFile, "Locallab", "Lapmaskcol_" + index_str, pedited, spot.lapmaskcol, spotEdited.lapmaskcol);
                assignFromKeyfile(keyFile, "Locallab", "QualityCurveMethod_" + index_str, pedited, spot.qualitycurveMethod, spotEdited.qualitycurveMethod);
                assignFromKeyfile(keyFile, "Locallab", "gridMethod_" + index_str, pedited, spot.gridMethod, spotEdited.gridMethod);
                assignFromKeyfile(keyFile, "Locallab", "Merg_Method_" + index_str, pedited, spot.merMethod, spotEdited.merMethod);
                assignFromKeyfile(keyFile, "Locallab", "ToneMethod_" + index_str, pedited, spot.toneMethod, spotEdited.toneMethod);
                assignFromKeyfile(keyFile, "Locallab", "mergecolMethod_" + index_str, pedited, spot.mergecolMethod, spotEdited.mergecolMethod);
                assignFromKeyfile(keyFile, "Locallab", "LLCurve_" + index_str, pedited, spot.llcurve, spotEdited.llcurve);
                assignFromKeyfile(keyFile, "Locallab", "LCCurve_" + index_str, pedited, spot.lccurve, spotEdited.lccurve);
                assignFromKeyfile(keyFile, "Locallab", "CCCurve_" + index_str, pedited, spot.cccurve, spotEdited.cccurve);
                assignFromKeyfile(keyFile, "Locallab", "CLCurve_" + index_str, pedited, spot.clcurve, spotEdited.clcurve);
                assignFromKeyfile(keyFile, "Locallab", "RGBCurve_" + index_str, pedited, spot.rgbcurve, spotEdited.rgbcurve);
                assignFromKeyfile(keyFile, "Locallab", "LHCurve_" + index_str, pedited, spot.LHcurve, spotEdited.LHcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHCurve_" + index_str, pedited, spot.HHcurve, spotEdited.HHcurve);
                assignFromKeyfile(keyFile, "Locallab", "CHCurve_" + index_str, pedited, spot.CHcurve, spotEdited.CHcurve);
                assignFromKeyfile(keyFile, "Locallab", "Invers_" + index_str, pedited, spot.invers, spotEdited.invers);
                assignFromKeyfile(keyFile, "Locallab", "Special_" + index_str, pedited, spot.special, spotEdited.special);
                assignFromKeyfile(keyFile, "Locallab", "Toolcol_" + index_str, pedited, spot.toolcol, spotEdited.toolcol);
                assignFromKeyfile(keyFile, "Locallab", "EnaColorMask_" + index_str, pedited, spot.enaColorMask, spotEdited.enaColorMask);
                assignFromKeyfile(keyFile, "Locallab", "FftColorMask_" + index_str, pedited, spot.fftColorMask, spotEdited.fftColorMask);
                assignFromKeyfile(keyFile, "Locallab", "CCmaskCurve_" + index_str, pedited, spot.CCmaskcurve, spotEdited.CCmaskcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskCurve_" + index_str, pedited, spot.LLmaskcurve, spotEdited.LLmaskcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmaskCurve_" + index_str, pedited, spot.HHmaskcurve, spotEdited.HHmaskcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHhmaskCurve_" + index_str, pedited, spot.HHhmaskcurve, spotEdited.HHhmaskcurve);
                assignFromKeyfile(keyFile, "Locallab", "Softradiuscol_" + index_str, pedited, spot.softradiuscol, spotEdited.softradiuscol);
                assignFromKeyfile(keyFile, "Locallab", "Opacol_" + index_str, pedited, spot.opacol, spotEdited.opacol);
                assignFromKeyfile(keyFile, "Locallab", "Mercol_" + index_str, pedited, spot.mercol, spotEdited.mercol);
                assignFromKeyfile(keyFile, "Locallab", "Merlucol_" + index_str, pedited, spot.merlucol, spotEdited.merlucol);
                assignFromKeyfile(keyFile, "Locallab", "Conthrcol_" + index_str, pedited, spot.conthrcol, spotEdited.conthrcol);
                assignFromKeyfile(keyFile, "Locallab", "LmaskCurve_" + index_str, pedited, spot.Lmaskcurve, spotEdited.Lmaskcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskcolCurvewav_" + index_str, pedited, spot.LLmaskcolcurvewav, spotEdited.LLmaskcolcurvewav);

                if (keyFile.has_key("Locallab", "CSThresholdcol_" + index_str)) {
                    const std::vector<int> thresh = keyFile.get_integer_list("Locallab", "CSThresholdcol_" + index_str);

                    if (thresh.size() >= 4) {
                        spot.csthresholdcol.setValues(thresh[0], thresh[1], min(thresh[2], 10), min(thresh[3], 10));
                    }

                    spotEdited.csthresholdcol = true;
                }

                // Exposure
                spot.visiexpose = assignFromKeyfile(keyFile, "Locallab", "Expexpose_" + index_str, pedited, spot.expexpose, spotEdited.expexpose);

                if (spot.visiexpose) {
                    spotEdited.visiexpose = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexexpose_" + index_str, pedited, spot.complexexpose, spotEdited.complexexpose);
                assignFromKeyfile(keyFile, "Locallab", "Expcomp_" + index_str, pedited, spot.expcomp, spotEdited.expcomp);
                assignFromKeyfile(keyFile, "Locallab", "Hlcompr_" + index_str, pedited, spot.hlcompr, spotEdited.hlcompr);
                assignFromKeyfile(keyFile, "Locallab", "Hlcomprthresh_" + index_str, pedited, spot.hlcomprthresh, spotEdited.hlcomprthresh);
                assignFromKeyfile(keyFile, "Locallab", "Black_" + index_str, pedited, spot.black, spotEdited.black);
                assignFromKeyfile(keyFile, "Locallab", "Shadex_" + index_str, pedited, spot.shadex, spotEdited.shadex);
                assignFromKeyfile(keyFile, "Locallab", "Shcompr_" + index_str, pedited, spot.shcompr, spotEdited.shcompr);
                assignFromKeyfile(keyFile, "Locallab", "Expchroma_" + index_str, pedited, spot.expchroma, spotEdited.expchroma);
                assignFromKeyfile(keyFile, "Locallab", "Sensiex_" + index_str, pedited, spot.sensiex, spotEdited.sensiex);
                assignFromKeyfile(keyFile, "Locallab", "Structexp_" + index_str, pedited, spot.structexp, spotEdited.structexp);
                assignFromKeyfile(keyFile, "Locallab", "Blurexpde_" + index_str, pedited, spot.blurexpde, spotEdited.blurexpde);
                assignFromKeyfile(keyFile, "Locallab", "Strexp_" + index_str, pedited, spot.strexp, spotEdited.strexp);
                assignFromKeyfile(keyFile, "Locallab", "Angexp_" + index_str, pedited, spot.angexp, spotEdited.angexp);
                assignFromKeyfile(keyFile, "Locallab", "ExCurve_" + index_str, pedited, spot.excurve, spotEdited.excurve);
                assignFromKeyfile(keyFile, "Locallab", "Inversex_" + index_str, pedited, spot.inversex, spotEdited.inversex);
                assignFromKeyfile(keyFile, "Locallab", "EnaExpMask_" + index_str, pedited, spot.enaExpMask, spotEdited.enaExpMask);
                assignFromKeyfile(keyFile, "Locallab", "EnaExpMaskaft_" + index_str, pedited, spot.enaExpMaskaft, spotEdited.enaExpMaskaft);
                assignFromKeyfile(keyFile, "Locallab", "CCmaskexpCurve_" + index_str, pedited, spot.CCmaskexpcurve, spotEdited.CCmaskexpcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskexpCurve_" + index_str, pedited, spot.LLmaskexpcurve, spotEdited.LLmaskexpcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmaskexpCurve_" + index_str, pedited, spot.HHmaskexpcurve, spotEdited.HHmaskexpcurve);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskexp_" + index_str, pedited, spot.blendmaskexp, spotEdited.blendmaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Radmaskexp_" + index_str, pedited, spot.radmaskexp, spotEdited.radmaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Chromaskexp_" + index_str, pedited, spot.chromaskexp, spotEdited.chromaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Gammaskexp_" + index_str, pedited, spot.gammaskexp, spotEdited.gammaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Slomaskexp_" + index_str, pedited, spot.slomaskexp, spotEdited.slomaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Lapmaskexp_" + index_str, pedited, spot.lapmaskexp, spotEdited.lapmaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Strmaskexp_" + index_str, pedited, spot.strmaskexp, spotEdited.strmaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Angmaskexp_" + index_str, pedited, spot.angmaskexp, spotEdited.angmaskexp);
                assignFromKeyfile(keyFile, "Locallab", "Softradiusexp_" + index_str, pedited, spot.softradiusexp, spotEdited.softradiusexp);
                assignFromKeyfile(keyFile, "Locallab", "LmaskexpCurve_" + index_str, pedited, spot.Lmaskexpcurve, spotEdited.Lmaskexpcurve);
                assignFromKeyfile(keyFile, "Locallab", "ExpMethod_" + index_str, pedited, spot.expMethod, spotEdited.expMethod);
                assignFromKeyfile(keyFile, "Locallab", "ExnoiseMethod_" + index_str, pedited, spot.exnoiseMethod, spotEdited.exnoiseMethod);
                assignFromKeyfile(keyFile, "Locallab", "Laplacexp_" + index_str, pedited, spot.laplacexp, spotEdited.laplacexp);
                assignFromKeyfile(keyFile, "Locallab", "Balanexp_" + index_str, pedited, spot.balanexp, spotEdited.balanexp);
                assignFromKeyfile(keyFile, "Locallab", "Linearexp_" + index_str, pedited, spot.linear, spotEdited.linear);
                assignFromKeyfile(keyFile, "Locallab", "Gamm_" + index_str, pedited, spot.gamm, spotEdited.gamm);
                assignFromKeyfile(keyFile, "Locallab", "Fatamount_" + index_str, pedited, spot.fatamount, spotEdited.fatamount);
                assignFromKeyfile(keyFile, "Locallab", "Fatdetail_" + index_str, pedited, spot.fatdetail, spotEdited.fatdetail);
                assignFromKeyfile(keyFile, "Locallab", "Fatanchor_" + index_str, pedited, spot.fatanchor, spotEdited.fatanchor);
                assignFromKeyfile(keyFile, "Locallab", "Fatlevel_" + index_str, pedited, spot.fatlevel, spotEdited.fatlevel);
                // Shadow highlight
                spot.visishadhigh = assignFromKeyfile(keyFile, "Locallab", "Expshadhigh_" + index_str, pedited, spot.expshadhigh, spotEdited.expshadhigh);

                if (spot.visishadhigh) {
                    spotEdited.visishadhigh = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexshadhigh_" + index_str, pedited, spot.complexshadhigh, spotEdited.complexshadhigh);
                assignFromKeyfile(keyFile, "Locallab", "ShMethod_" + index_str, pedited, spot.shMethod, spotEdited.shMethod);

                for (int j = 0; j < 5; j ++) {
                    assignFromKeyfile(keyFile, "Locallab", "Multsh" + std::to_string(j) + "_" + index_str, pedited, spot.multsh[j], spotEdited.multsh[j]);
                }

                assignFromKeyfile(keyFile, "Locallab", "Expshadhigh_" + index_str, pedited, spot.expshadhigh, spotEdited.expshadhigh);
                assignFromKeyfile(keyFile, "Locallab", "highlights_" + index_str, pedited, spot.highlights, spotEdited.highlights);
                assignFromKeyfile(keyFile, "Locallab", "h_tonalwidth_" + index_str, pedited, spot.h_tonalwidth, spotEdited.h_tonalwidth);
                assignFromKeyfile(keyFile, "Locallab", "shadows_" + index_str, pedited, spot.shadows, spotEdited.shadows);
                assignFromKeyfile(keyFile, "Locallab", "s_tonalwidth_" + index_str, pedited, spot.s_tonalwidth, spotEdited.s_tonalwidth);
                assignFromKeyfile(keyFile, "Locallab", "sh_radius_" + index_str, pedited, spot.sh_radius, spotEdited.sh_radius);
                assignFromKeyfile(keyFile, "Locallab", "sensihs_" + index_str, pedited, spot.sensihs, spotEdited.sensihs);
                assignFromKeyfile(keyFile, "Locallab", "EnaSHMask_" + index_str, pedited, spot.enaSHMask, spotEdited.enaSHMask);
                assignFromKeyfile(keyFile, "Locallab", "CCmaskSHCurve_" + index_str, pedited, spot.CCmaskSHcurve, spotEdited.CCmaskSHcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskSHCurve_" + index_str, pedited, spot.LLmaskSHcurve, spotEdited.LLmaskSHcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmaskSHCurve_" + index_str, pedited, spot.HHmaskSHcurve, spotEdited.HHmaskSHcurve);
                assignFromKeyfile(keyFile, "Locallab", "BlendmaskSH_" + index_str, pedited, spot.blendmaskSH, spotEdited.blendmaskSH);
                assignFromKeyfile(keyFile, "Locallab", "RadmaskSH_" + index_str, pedited, spot.radmaskSH, spotEdited.radmaskSH);
                assignFromKeyfile(keyFile, "Locallab", "BlurSHde_" + index_str, pedited, spot.blurSHde, spotEdited.blurSHde);
                assignFromKeyfile(keyFile, "Locallab", "StrSH_" + index_str, pedited, spot.strSH, spotEdited.strSH);
                assignFromKeyfile(keyFile, "Locallab", "AngSH_" + index_str, pedited, spot.angSH, spotEdited.angSH);
                assignFromKeyfile(keyFile, "Locallab", "Inverssh_" + index_str, pedited, spot.inverssh, spotEdited.inverssh);
                assignFromKeyfile(keyFile, "Locallab", "ChromaskSH_" + index_str, pedited, spot.chromaskSH, spotEdited.chromaskSH);
                assignFromKeyfile(keyFile, "Locallab", "GammaskSH_" + index_str, pedited, spot.gammaskSH, spotEdited.gammaskSH);
                assignFromKeyfile(keyFile, "Locallab", "SlomaskSH_" + index_str, pedited, spot.slomaskSH, spotEdited.slomaskSH);
                assignFromKeyfile(keyFile, "Locallab", "LapmaskSH_" + index_str, pedited, spot.lapmaskSH, spotEdited.lapmaskSH);
                assignFromKeyfile(keyFile, "Locallab", "DetailSH_" + index_str, pedited, spot.detailSH, spotEdited.detailSH);
                assignFromKeyfile(keyFile, "Locallab", "LmaskSHCurve_" + index_str, pedited, spot.LmaskSHcurve, spotEdited.LmaskSHcurve);
                assignFromKeyfile(keyFile, "Locallab", "FatamountSH_" + index_str, pedited, spot.fatamountSH, spotEdited.fatamountSH);
                assignFromKeyfile(keyFile, "Locallab", "FatanchorSH_" + index_str, pedited, spot.fatanchorSH, spotEdited.fatanchorSH);
                assignFromKeyfile(keyFile, "Locallab", "GamSH_" + index_str, pedited, spot.gamSH, spotEdited.gamSH);
                assignFromKeyfile(keyFile, "Locallab", "SloSH_" + index_str, pedited, spot.sloSH, spotEdited.sloSH);
                // Vibrance
                spot.visivibrance = assignFromKeyfile(keyFile, "Locallab", "Expvibrance_" + index_str, pedited, spot.expvibrance, spotEdited.expvibrance);

                if (spot.visivibrance) {
                    spotEdited.visivibrance = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexvibrance_" + index_str, pedited, spot.complexvibrance, spotEdited.complexvibrance);
                assignFromKeyfile(keyFile, "Locallab", "Saturated_" + index_str, pedited, spot.saturated, spotEdited.saturated);
                assignFromKeyfile(keyFile, "Locallab", "Pastels_" + index_str, pedited, spot.pastels, spotEdited.pastels);
                assignFromKeyfile(keyFile, "Locallab", "Warm_" + index_str, pedited, spot.warm, spotEdited.warm);

                if (keyFile.has_key("Locallab", "PSThreshold_" + index_str)) {
                    const std::vector<int> thresh = keyFile.get_integer_list("Locallab", "PSThreshold_" + index_str);

                    if (thresh.size() >= 2) {
                        spot.psthreshold.setValues(thresh[0], thresh[1]);
                    }

                    spotEdited.psthreshold = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "ProtectSkins_" + index_str, pedited, spot.protectskins, spotEdited.protectskins);
                assignFromKeyfile(keyFile, "Locallab", "AvoidColorShift_" + index_str, pedited, spot.avoidcolorshift, spotEdited.avoidcolorshift);
                assignFromKeyfile(keyFile, "Locallab", "PastSatTog_" + index_str, pedited, spot.pastsattog, spotEdited.pastsattog);
                assignFromKeyfile(keyFile, "Locallab", "Sensiv_" + index_str, pedited, spot.sensiv, spotEdited.sensiv);
                assignFromKeyfile(keyFile, "Locallab", "SkinTonesCurve_" + index_str, pedited, spot.skintonescurve, spotEdited.skintonescurve);
                assignFromKeyfile(keyFile, "Locallab", "CCmaskvibCurve_" + index_str, pedited, spot.CCmaskvibcurve, spotEdited.CCmaskvibcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskvibCurve_" + index_str, pedited, spot.LLmaskvibcurve, spotEdited.LLmaskvibcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmaskvibCurve_" + index_str, pedited, spot.HHmaskvibcurve, spotEdited.HHmaskvibcurve);
                assignFromKeyfile(keyFile, "Locallab", "EnavibMask_" + index_str, pedited, spot.enavibMask, spotEdited.enavibMask);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskvib_" + index_str, pedited, spot.blendmaskvib, spotEdited.blendmaskvib);
                assignFromKeyfile(keyFile, "Locallab", "Radmaskvib_" + index_str, pedited, spot.radmaskvib, spotEdited.radmaskvib);
                assignFromKeyfile(keyFile, "Locallab", "Chromaskvib_" + index_str, pedited, spot.chromaskvib, spotEdited.chromaskvib);
                assignFromKeyfile(keyFile, "Locallab", "Gammaskvib_" + index_str, pedited, spot.gammaskvib, spotEdited.gammaskvib);
                assignFromKeyfile(keyFile, "Locallab", "Slomaskvib_" + index_str, pedited, spot.slomaskvib, spotEdited.slomaskvib);
                assignFromKeyfile(keyFile, "Locallab", "Lapmaskvib_" + index_str, pedited, spot.lapmaskvib, spotEdited.lapmaskvib);
                assignFromKeyfile(keyFile, "Locallab", "Strvib_" + index_str, pedited, spot.strvib, spotEdited.strvib);
                assignFromKeyfile(keyFile, "Locallab", "Strvibab_" + index_str, pedited, spot.strvibab, spotEdited.strvibab);
                assignFromKeyfile(keyFile, "Locallab", "Strvibh_" + index_str, pedited, spot.strvibh, spotEdited.strvibh);
                assignFromKeyfile(keyFile, "Locallab", "Angvib_" + index_str, pedited, spot.angvib, spotEdited.angvib);
                assignFromKeyfile(keyFile, "Locallab", "LmaskvibCurve_" + index_str, pedited, spot.Lmaskvibcurve, spotEdited.Lmaskvibcurve);
                // Soft Light
                spot.visisoft = assignFromKeyfile(keyFile, "Locallab", "Expsoft_" + index_str, pedited, spot.expsoft, spotEdited.expsoft);

                if (spot.visisoft) {
                    spotEdited.visisoft = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexsoft_" + index_str, pedited, spot.complexsoft, spotEdited.complexsoft);
                assignFromKeyfile(keyFile, "Locallab", "Streng_" + index_str, pedited, spot.streng, spotEdited.streng);
                assignFromKeyfile(keyFile, "Locallab", "Sensisf_" + index_str, pedited, spot.sensisf, spotEdited.sensisf);
                assignFromKeyfile(keyFile, "Locallab", "Laplace_" + index_str, pedited, spot.laplace, spotEdited.laplace);
                assignFromKeyfile(keyFile, "Locallab", "SoftMethod_" + index_str, pedited, spot.softMethod, spotEdited.softMethod);
                // Blur & Noise
                spot.visiblur = assignFromKeyfile(keyFile, "Locallab", "Expblur_" + index_str, pedited, spot.expblur, spotEdited.expblur);

                if (spot.visiblur) {
                    spotEdited.visiblur = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexblur_" + index_str, pedited, spot.complexblur, spotEdited.complexblur);
                assignFromKeyfile(keyFile, "Locallab", "Radius_" + index_str, pedited, spot.radius, spotEdited.radius);
                assignFromKeyfile(keyFile, "Locallab", "Strength_" + index_str, pedited, spot.strength, spotEdited.strength);
                assignFromKeyfile(keyFile, "Locallab", "Sensibn_" + index_str, pedited, spot.sensibn, spotEdited.sensibn);
                assignFromKeyfile(keyFile, "Locallab", "Iteramed_" + index_str, pedited, spot.itera, spotEdited.itera);
                assignFromKeyfile(keyFile, "Locallab", "Guidbl_" + index_str, pedited, spot.guidbl, spotEdited.guidbl);
                assignFromKeyfile(keyFile, "Locallab", "Strbl_" + index_str, pedited, spot.strbl, spotEdited.strbl);
                assignFromKeyfile(keyFile, "Locallab", "Isogr_" + index_str, pedited, spot.isogr, spotEdited.isogr);
                assignFromKeyfile(keyFile, "Locallab", "Strengr_" + index_str, pedited, spot.strengr, spotEdited.strengr);
                assignFromKeyfile(keyFile, "Locallab", "Scalegr_" + index_str, pedited, spot.scalegr, spotEdited.scalegr);
                assignFromKeyfile(keyFile, "Locallab", "Epsbl_" + index_str, pedited, spot.epsbl, spotEdited.epsbl);
                assignFromKeyfile(keyFile, "Locallab", "BlMethod_" + index_str, pedited, spot.blMethod, spotEdited.blMethod);
                assignFromKeyfile(keyFile, "Locallab", "ChroMethod_" + index_str, pedited, spot.chroMethod, spotEdited.chroMethod);
                assignFromKeyfile(keyFile, "Locallab", "BlurMethod_" + index_str, pedited, spot.blurMethod, spotEdited.blurMethod);
                assignFromKeyfile(keyFile, "Locallab", "MedMethod_" + index_str, pedited, spot.medMethod, spotEdited.medMethod);
                assignFromKeyfile(keyFile, "Locallab", "activlum_" + index_str, pedited, spot.activlum, spotEdited.activlum);
                assignFromKeyfile(keyFile, "Locallab", "noiselumf_" + index_str, pedited, spot.noiselumf, spotEdited.noiselumf);
                assignFromKeyfile(keyFile, "Locallab", "noiselumf0_" + index_str, pedited, spot.noiselumf0, spotEdited.noiselumf0);
                assignFromKeyfile(keyFile, "Locallab", "noiselumf2_" + index_str, pedited, spot.noiselumf2, spotEdited.noiselumf2);
                assignFromKeyfile(keyFile, "Locallab", "noiselumc_" + index_str, pedited, spot.noiselumc, spotEdited.noiselumc);
                assignFromKeyfile(keyFile, "Locallab", "noiselumdetail_" + index_str, pedited, spot.noiselumdetail, spotEdited.noiselumdetail);
                assignFromKeyfile(keyFile, "Locallab", "noiselequal_" + index_str, pedited, spot.noiselequal, spotEdited.noiselequal);
                assignFromKeyfile(keyFile, "Locallab", "noisechrof_" + index_str, pedited, spot.noisechrof, spotEdited.noisechrof);
                assignFromKeyfile(keyFile, "Locallab", "noisechroc_" + index_str, pedited, spot.noisechroc, spotEdited.noisechroc);
                assignFromKeyfile(keyFile, "Locallab", "noisechrodetail_" + index_str, pedited, spot.noisechrodetail, spotEdited.noisechrodetail);
                assignFromKeyfile(keyFile, "Locallab", "Adjblur_" + index_str, pedited, spot.adjblur, spotEdited.adjblur);
                assignFromKeyfile(keyFile, "Locallab", "Bilateral_" + index_str, pedited, spot.bilateral, spotEdited.bilateral);
                assignFromKeyfile(keyFile, "Locallab", "Sensiden_" + index_str, pedited, spot.sensiden, spotEdited.sensiden);
                assignFromKeyfile(keyFile, "Locallab", "Detailthr_" + index_str, pedited, spot.detailthr, spotEdited.detailthr);
                assignFromKeyfile(keyFile, "Locallab", "LocwavCurveden_" + index_str, pedited, spot.locwavcurveden, spotEdited.locwavcurveden);
                assignFromKeyfile(keyFile, "Locallab", "Showmasktyp_" + index_str, pedited, spot.showmaskblMethodtyp, spotEdited.showmaskblMethodtyp);
                assignFromKeyfile(keyFile, "Locallab", "CCmaskblCurve_" + index_str, pedited, spot.CCmaskblcurve, spotEdited.CCmaskblcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskblCurve_" + index_str, pedited, spot.LLmaskblcurve, spotEdited.LLmaskblcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmaskblCurve_" + index_str, pedited, spot.HHmaskblcurve, spotEdited.HHmaskblcurve);
                assignFromKeyfile(keyFile, "Locallab", "EnablMask_" + index_str, pedited, spot.enablMask, spotEdited.enablMask);
                assignFromKeyfile(keyFile, "Locallab", "Fftwbl_" + index_str, pedited, spot.fftwbl, spotEdited.fftwbl);
                assignFromKeyfile(keyFile, "Locallab", "Toolbl_" + index_str, pedited, spot.toolbl, spotEdited.toolbl);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskbl_" + index_str, pedited, spot.blendmaskbl, spotEdited.blendmaskbl);
                assignFromKeyfile(keyFile, "Locallab", "Radmaskbl_" + index_str, pedited, spot.radmaskbl, spotEdited.radmaskbl);
                assignFromKeyfile(keyFile, "Locallab", "Chromaskbl_" + index_str, pedited, spot.chromaskbl, spotEdited.chromaskbl);
                assignFromKeyfile(keyFile, "Locallab", "Gammaskbl_" + index_str, pedited, spot.gammaskbl, spotEdited.gammaskbl);
                assignFromKeyfile(keyFile, "Locallab", "Slomaskbl_" + index_str, pedited, spot.slomaskbl, spotEdited.slomaskbl);
                assignFromKeyfile(keyFile, "Locallab", "Lapmaskbl_" + index_str, pedited, spot.lapmaskbl, spotEdited.lapmaskbl);
                assignFromKeyfile(keyFile, "Locallab", "shadmaskbl_" + index_str, pedited, spot.shadmaskbl, spotEdited.shadmaskbl);
                assignFromKeyfile(keyFile, "Locallab", "shadmaskblsha_" + index_str, pedited, spot.shadmaskblsha, spotEdited.shadmaskblsha);
                assignFromKeyfile(keyFile, "Locallab", "strumaskbl_" + index_str, pedited, spot.strumaskbl, spotEdited.strumaskbl);
                assignFromKeyfile(keyFile, "Locallab", "LmaskblCurve_" + index_str, pedited, spot.Lmaskblcurve, spotEdited.Lmaskblcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskblCurvewav_" + index_str, pedited, spot.LLmaskblcurvewav, spotEdited.LLmaskblcurvewav);

                if (keyFile.has_key("Locallab", "CSThresholdblur_" + index_str)) {
                    const std::vector<int> thresh = keyFile.get_integer_list("Locallab", "CSThresholdblur_" + index_str);

                    if (thresh.size() >= 4) {
                        spot.csthresholdblur.setValues(thresh[0], thresh[1], min(thresh[2], 10), min(thresh[3], 10));
                    }

                    spotEdited.csthresholdblur = true;
                }
                // Tone Mapping
                spot.visitonemap = assignFromKeyfile(keyFile, "Locallab", "Exptonemap_" + index_str, pedited, spot.exptonemap, spotEdited.exptonemap);

                if (spot.visitonemap) {
                    spotEdited.visitonemap = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complextonemap_" + index_str, pedited, spot.complextonemap, spotEdited.complextonemap);
                assignFromKeyfile(keyFile, "Locallab", "Stren_" + index_str, pedited, spot.stren, spotEdited.stren);
                assignFromKeyfile(keyFile, "Locallab", "Gamma_" + index_str, pedited, spot.gamma, spotEdited.gamma);
                assignFromKeyfile(keyFile, "Locallab", "Estop_" + index_str, pedited, spot.estop, spotEdited.estop);
                assignFromKeyfile(keyFile, "Locallab", "Scaltm_" + index_str, pedited, spot.scaltm, spotEdited.scaltm);
                assignFromKeyfile(keyFile, "Locallab", "Rewei_" + index_str, pedited, spot.rewei, spotEdited.rewei);
                assignFromKeyfile(keyFile, "Locallab", "Satur_" + index_str, pedited, spot.satur, spotEdited.satur);
                assignFromKeyfile(keyFile, "Locallab", "Sensitm_" + index_str, pedited, spot.sensitm, spotEdited.sensitm);
                assignFromKeyfile(keyFile, "Locallab", "Softradiustm_" + index_str, pedited, spot.softradiustm, spotEdited.softradiustm);
                assignFromKeyfile(keyFile, "Locallab", "Amount_" + index_str, pedited, spot.amount, spotEdited.amount);
                assignFromKeyfile(keyFile, "Locallab", "Equiltm_" + index_str, pedited, spot.equiltm, spotEdited.equiltm);
                assignFromKeyfile(keyFile, "Locallab", "CCmasktmCurve_" + index_str, pedited, spot.CCmasktmcurve, spotEdited.CCmasktmcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmasktmCurve_" + index_str, pedited, spot.LLmasktmcurve, spotEdited.LLmasktmcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmasktmCurve_" + index_str, pedited, spot.HHmasktmcurve, spotEdited.HHmasktmcurve);
                assignFromKeyfile(keyFile, "Locallab", "EnatmMask_" + index_str, pedited, spot.enatmMask, spotEdited.enatmMask);
                assignFromKeyfile(keyFile, "Locallab", "EnatmMaskaft_" + index_str, pedited, spot.enatmMaskaft, spotEdited.enatmMaskaft);
                assignFromKeyfile(keyFile, "Locallab", "Blendmasktm_" + index_str, pedited, spot.blendmasktm, spotEdited.blendmasktm);
                assignFromKeyfile(keyFile, "Locallab", "Radmasktm_" + index_str, pedited, spot.radmasktm, spotEdited.radmasktm);
                assignFromKeyfile(keyFile, "Locallab", "Chromasktm_" + index_str, pedited, spot.chromasktm, spotEdited.chromasktm);
                assignFromKeyfile(keyFile, "Locallab", "Gammasktm_" + index_str, pedited, spot.gammasktm, spotEdited.gammasktm);
                assignFromKeyfile(keyFile, "Locallab", "Slomasktm_" + index_str, pedited, spot.slomasktm, spotEdited.slomasktm);
                assignFromKeyfile(keyFile, "Locallab", "Lapmasktm_" + index_str, pedited, spot.lapmasktm, spotEdited.lapmasktm);
                assignFromKeyfile(keyFile, "Locallab", "LmasktmCurve_" + index_str, pedited, spot.Lmasktmcurve, spotEdited.Lmasktmcurve);
                // Retinex
                spot.visireti = assignFromKeyfile(keyFile, "Locallab", "Expreti_" + index_str, pedited, spot.expreti, spotEdited.expreti);

                if (spot.visireti) {
                    spotEdited.visireti = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexreti_" + index_str, pedited, spot.complexreti, spotEdited.complexreti);
                assignFromKeyfile(keyFile, "Locallab", "retinexMethod_" + index_str, pedited, spot.retinexMethod, spotEdited.retinexMethod);
                assignFromKeyfile(keyFile, "Locallab", "Str_" + index_str, pedited, spot.str, spotEdited.str);
                assignFromKeyfile(keyFile, "Locallab", "Chrrt_" + index_str, pedited, spot.chrrt, spotEdited.chrrt);
                assignFromKeyfile(keyFile, "Locallab", "Neigh_" + index_str, pedited, spot.neigh, spotEdited.neigh);
                assignFromKeyfile(keyFile, "Locallab", "Vart_" + index_str, pedited, spot.vart, spotEdited.vart);
                assignFromKeyfile(keyFile, "Locallab", "Offs_" + index_str, pedited, spot.offs, spotEdited.offs);
                assignFromKeyfile(keyFile, "Locallab", "Dehaz_" + index_str, pedited, spot.dehaz, spotEdited.dehaz);
                assignFromKeyfile(keyFile, "Locallab", "Depth_" + index_str, pedited, spot.depth, spotEdited.depth);
                assignFromKeyfile(keyFile, "Locallab", "Sensih_" + index_str, pedited, spot.sensih, spotEdited.sensih);
                assignFromKeyfile(keyFile, "Locallab", "TgainCurve_" + index_str, pedited, spot.localTgaincurve, spotEdited.localTgaincurve);
                assignFromKeyfile(keyFile, "Locallab", "TtransCurve_" + index_str, pedited, spot.localTtranscurve, spotEdited.localTtranscurve);
                assignFromKeyfile(keyFile, "Locallab", "Inversret_" + index_str, pedited, spot.inversret, spotEdited.inversret);
                assignFromKeyfile(keyFile, "Locallab", "Equilret_" + index_str, pedited, spot.equilret, spotEdited.equilret);
                assignFromKeyfile(keyFile, "Locallab", "Loglin_" + index_str, pedited, spot.loglin, spotEdited.loglin);
                assignFromKeyfile(keyFile, "Locallab", "Lumonly_" + index_str, pedited, spot.lumonly, spotEdited.lumonly);
                assignFromKeyfile(keyFile, "Locallab", "Softradiusret_" + index_str, pedited, spot.softradiusret, spotEdited.softradiusret);
                assignFromKeyfile(keyFile, "Locallab", "CCmaskretiCurve_" + index_str, pedited, spot.CCmaskreticurve, spotEdited.CCmaskreticurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskretiCurve_" + index_str, pedited, spot.LLmaskreticurve, spotEdited.LLmaskreticurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmaskretiCurve_" + index_str, pedited, spot.HHmaskreticurve, spotEdited.HHmaskreticurve);
                assignFromKeyfile(keyFile, "Locallab", "EnaretiMask_" + index_str, pedited, spot.enaretiMask, spotEdited.enaretiMask);
                assignFromKeyfile(keyFile, "Locallab", "EnaretiMasktmap_" + index_str, pedited, spot.enaretiMasktmap, spotEdited.enaretiMasktmap);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskreti_" + index_str, pedited, spot.blendmaskreti, spotEdited.blendmaskreti);
                assignFromKeyfile(keyFile, "Locallab", "Radmaskreti_" + index_str, pedited, spot.radmaskreti, spotEdited.radmaskreti);
                assignFromKeyfile(keyFile, "Locallab", "Chromaskreti_" + index_str, pedited, spot.chromaskreti, spotEdited.chromaskreti);
                assignFromKeyfile(keyFile, "Locallab", "Gammaskreti_" + index_str, pedited, spot.gammaskreti, spotEdited.gammaskreti);
                assignFromKeyfile(keyFile, "Locallab", "Slomaskreti_" + index_str, pedited, spot.slomaskreti, spotEdited.slomaskreti);
                assignFromKeyfile(keyFile, "Locallab", "Lapmaskreti_" + index_str, pedited, spot.lapmaskreti, spotEdited.lapmaskreti);
                assignFromKeyfile(keyFile, "Locallab", "Scalereti_" + index_str, pedited, spot.scalereti, spotEdited.scalereti);
                assignFromKeyfile(keyFile, "Locallab", "Darkness_" + index_str, pedited, spot.darkness, spotEdited.darkness);
                assignFromKeyfile(keyFile, "Locallab", "Lightnessreti_" + index_str, pedited, spot.lightnessreti, spotEdited.lightnessreti);
                assignFromKeyfile(keyFile, "Locallab", "Limd_" + index_str, pedited, spot.limd, spotEdited.limd);
                assignFromKeyfile(keyFile, "Locallab", "Cliptm_" + index_str, pedited, spot.cliptm, spotEdited.cliptm);
                assignFromKeyfile(keyFile, "Locallab", "Fftwreti_" + index_str, pedited, spot.fftwreti, spotEdited.fftwreti);
                assignFromKeyfile(keyFile, "Locallab", "LmaskretiCurve_" + index_str, pedited, spot.Lmaskreticurve, spotEdited.Lmaskreticurve);
                // Sharpening
                spot.visisharp = assignFromKeyfile(keyFile, "Locallab", "Expsharp_" + index_str, pedited, spot.expsharp, spotEdited.expsharp);

                if (spot.visisharp) {
                    spotEdited.visisharp = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexsharp_" + index_str, pedited, spot.complexsharp, spotEdited.complexsharp);
                assignFromKeyfile(keyFile, "Locallab", "Sharcontrast_" + index_str, pedited, spot.sharcontrast, spotEdited.sharcontrast);
                assignFromKeyfile(keyFile, "Locallab", "Sharradius_" + index_str, pedited, spot.sharradius, spotEdited.sharradius);
                assignFromKeyfile(keyFile, "Locallab", "Sharamount_" + index_str, pedited, spot.sharamount, spotEdited.sharamount);
                assignFromKeyfile(keyFile, "Locallab", "Shardamping_" + index_str, pedited, spot.shardamping, spotEdited.shardamping);
                assignFromKeyfile(keyFile, "Locallab", "Shariter_" + index_str, pedited, spot.shariter, spotEdited.shariter);
                assignFromKeyfile(keyFile, "Locallab", "Sharblur_" + index_str, pedited, spot.sharblur, spotEdited.sharblur);
                assignFromKeyfile(keyFile, "Locallab", "Sensisha_" + index_str, pedited, spot.sensisha, spotEdited.sensisha);
                assignFromKeyfile(keyFile, "Locallab", "Inverssha_" + index_str, pedited, spot.inverssha, spotEdited.inverssha);
                // Local Contrast
                spot.visicontrast = assignFromKeyfile(keyFile, "Locallab", "Expcontrast_" + index_str, pedited, spot.expcontrast, spotEdited.expcontrast);

                if (spot.visicontrast) {
                    spotEdited.visicontrast = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexcontrast_" + index_str, pedited, spot.complexcontrast, spotEdited.complexcontrast);
                assignFromKeyfile(keyFile, "Locallab", "Lcradius_" + index_str, pedited, spot.lcradius, spotEdited.lcradius);
                assignFromKeyfile(keyFile, "Locallab", "Lcamount_" + index_str, pedited, spot.lcamount, spotEdited.lcamount);
                assignFromKeyfile(keyFile, "Locallab", "Lcdarkness_" + index_str, pedited, spot.lcdarkness, spotEdited.lcdarkness);
                assignFromKeyfile(keyFile, "Locallab", "Lclightness_" + index_str, pedited, spot.lclightness, spotEdited.lclightness);
                assignFromKeyfile(keyFile, "Locallab", "Sigmalc_" + index_str, pedited, spot.sigmalc, spotEdited.sigmalc);
                assignFromKeyfile(keyFile, "Locallab", "Levelwav_" + index_str, pedited, spot.levelwav, spotEdited.levelwav);
                assignFromKeyfile(keyFile, "Locallab", "Residcont_" + index_str, pedited, spot.residcont, spotEdited.residcont);
                assignFromKeyfile(keyFile, "Locallab", "Residsha_" + index_str, pedited, spot.residsha, spotEdited.residsha);
                assignFromKeyfile(keyFile, "Locallab", "Residshathr_" + index_str, pedited, spot.residshathr, spotEdited.residshathr);
                assignFromKeyfile(keyFile, "Locallab", "Residhi_" + index_str, pedited, spot.residhi, spotEdited.residhi);
                assignFromKeyfile(keyFile, "Locallab", "Residhithr_" + index_str, pedited, spot.residhithr, spotEdited.residhithr);
                assignFromKeyfile(keyFile, "Locallab", "Residblur_" + index_str, pedited, spot.residblur, spotEdited.residblur);
                assignFromKeyfile(keyFile, "Locallab", "Levelblur_" + index_str, pedited, spot.levelblur, spotEdited.levelblur);
                assignFromKeyfile(keyFile, "Locallab", "Sigmabl_" + index_str, pedited, spot.sigmabl, spotEdited.sigmabl);
                assignFromKeyfile(keyFile, "Locallab", "Residchro_" + index_str, pedited, spot.residchro, spotEdited.residchro);
                assignFromKeyfile(keyFile, "Locallab", "Residcomp_" + index_str, pedited, spot.residcomp, spotEdited.residcomp);
                assignFromKeyfile(keyFile, "Locallab", "Sigma_" + index_str, pedited, spot.sigma, spotEdited.sigma);
                assignFromKeyfile(keyFile, "Locallab", "Offset_" + index_str, pedited, spot.offset, spotEdited.offset);
                assignFromKeyfile(keyFile, "Locallab", "Sigmadr_" + index_str, pedited, spot.sigmadr, spotEdited.sigmadr);
                assignFromKeyfile(keyFile, "Locallab", "Threswav_" + index_str, pedited, spot.threswav, spotEdited.threswav);
                assignFromKeyfile(keyFile, "Locallab", "Chromalev_" + index_str, pedited, spot.chromalev, spotEdited.chromalev);
                assignFromKeyfile(keyFile, "Locallab", "Chromablu_" + index_str, pedited, spot.chromablu, spotEdited.chromablu);
                assignFromKeyfile(keyFile, "Locallab", "sigmadc_" + index_str, pedited, spot.sigmadc, spotEdited.sigmadc);
                assignFromKeyfile(keyFile, "Locallab", "deltad_" + index_str, pedited, spot.deltad, spotEdited.deltad);
                assignFromKeyfile(keyFile, "Locallab", "Fatres_" + index_str, pedited, spot.fatres, spotEdited.fatres);
                assignFromKeyfile(keyFile, "Locallab", "ClariLres_" + index_str, pedited, spot.clarilres, spotEdited.clarilres);
                assignFromKeyfile(keyFile, "Locallab", "ClariCres_" + index_str, pedited, spot.claricres, spotEdited.claricres);
                assignFromKeyfile(keyFile, "Locallab", "Clarisoft_" + index_str, pedited, spot.clarisoft, spotEdited.clarisoft);
                assignFromKeyfile(keyFile, "Locallab", "Sigmalc2_" + index_str, pedited, spot.sigmalc2, spotEdited.sigmalc2);
                assignFromKeyfile(keyFile, "Locallab", "Strwav_" + index_str, pedited, spot.strwav, spotEdited.strwav);
                assignFromKeyfile(keyFile, "Locallab", "Angwav_" + index_str, pedited, spot.angwav, spotEdited.angwav);
                assignFromKeyfile(keyFile, "Locallab", "Strengthw_" + index_str, pedited, spot.strengthw, spotEdited.strengthw);
                assignFromKeyfile(keyFile, "Locallab", "Sigmaed_" + index_str, pedited, spot.sigmaed, spotEdited.sigmaed);
                assignFromKeyfile(keyFile, "Locallab", "Radiusw_" + index_str, pedited, spot.radiusw, spotEdited.radiusw);
                assignFromKeyfile(keyFile, "Locallab", "Detailw_" + index_str, pedited, spot.detailw, spotEdited.detailw);
                assignFromKeyfile(keyFile, "Locallab", "Gradw_" + index_str, pedited, spot.gradw, spotEdited.gradw);
                assignFromKeyfile(keyFile, "Locallab", "Tloww_" + index_str, pedited, spot.tloww, spotEdited.tloww);
                assignFromKeyfile(keyFile, "Locallab", "Thigw_" + index_str, pedited, spot.thigw, spotEdited.thigw);
                assignFromKeyfile(keyFile, "Locallab", "Edgw_" + index_str, pedited, spot.edgw, spotEdited.edgw);
                assignFromKeyfile(keyFile, "Locallab", "Basew_" + index_str, pedited, spot.basew, spotEdited.basew);
                assignFromKeyfile(keyFile, "Locallab", "Sensilc_" + index_str, pedited, spot.sensilc, spotEdited.sensilc);
                assignFromKeyfile(keyFile, "Locallab", "Fftwlc_" + index_str, pedited, spot.fftwlc, spotEdited.fftwlc);
                assignFromKeyfile(keyFile, "Locallab", "Blurlc_" + index_str, pedited, spot.blurlc, spotEdited.blurlc);
                assignFromKeyfile(keyFile, "Locallab", "Wavblur_" + index_str, pedited, spot.wavblur, spotEdited.wavblur);
                assignFromKeyfile(keyFile, "Locallab", "Wavedg_" + index_str, pedited, spot.wavedg, spotEdited.wavedg);
                assignFromKeyfile(keyFile, "Locallab", "Waveshow_" + index_str, pedited, spot.waveshow, spotEdited.waveshow);
                assignFromKeyfile(keyFile, "Locallab", "Wavcont_" + index_str, pedited, spot.wavcont, spotEdited.wavcont);
                assignFromKeyfile(keyFile, "Locallab", "Wavcomp_" + index_str, pedited, spot.wavcomp, spotEdited.wavcomp);
                assignFromKeyfile(keyFile, "Locallab", "Wavgradl_" + index_str, pedited, spot.wavgradl, spotEdited.wavgradl);
                assignFromKeyfile(keyFile, "Locallab", "Wavcompre_" + index_str, pedited, spot.wavcompre, spotEdited.wavcompre);
                assignFromKeyfile(keyFile, "Locallab", "Origlc_" + index_str, pedited, spot.origlc, spotEdited.origlc);
                assignFromKeyfile(keyFile, "Locallab", "localcontMethod_" + index_str, pedited, spot.localcontMethod, spotEdited.localcontMethod);
                assignFromKeyfile(keyFile, "Locallab", "localedgMethod_" + index_str, pedited, spot.localedgMethod, spotEdited.localedgMethod);
                assignFromKeyfile(keyFile, "Locallab", "localneiMethod_" + index_str, pedited, spot.localneiMethod, spotEdited.localneiMethod);
                assignFromKeyfile(keyFile, "Locallab", "LocwavCurve_" + index_str, pedited, spot.locwavcurve, spotEdited.locwavcurve);
                assignFromKeyfile(keyFile, "Locallab", "LoclevwavCurve_" + index_str, pedited, spot.loclevwavcurve, spotEdited.loclevwavcurve);
                assignFromKeyfile(keyFile, "Locallab", "LocconwavCurve_" + index_str, pedited, spot.locconwavcurve, spotEdited.locconwavcurve);
                assignFromKeyfile(keyFile, "Locallab", "LoccompwavCurve_" + index_str, pedited, spot.loccompwavcurve, spotEdited.loccompwavcurve);
                assignFromKeyfile(keyFile, "Locallab", "LoccomprewavCurve_" + index_str, pedited, spot.loccomprewavcurve, spotEdited.loccomprewavcurve);
                assignFromKeyfile(keyFile, "Locallab", "LocedgwavCurve_" + index_str, pedited, spot.locedgwavcurve, spotEdited.locedgwavcurve);

                if (keyFile.has_key("Locallab", "CSThreshold_" + index_str)) {

                    const std::vector<int> thresh = keyFile.get_integer_list("Locallab", "CSThreshold_" + index_str);

                    if (thresh.size() >= 4) {
                        spot.csthreshold.setValues(thresh[0], thresh[1], min(thresh[2], 10), min(thresh[3], 10));
                    }

                    spotEdited.csthreshold = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "CCmasklcCurve_" + index_str, pedited, spot.CCmasklccurve, spotEdited.CCmasklccurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmasklcCurve_" + index_str, pedited, spot.LLmasklccurve, spotEdited.LLmasklccurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmasklcCurve_" + index_str, pedited, spot.HHmasklccurve, spotEdited.HHmasklccurve);
                assignFromKeyfile(keyFile, "Locallab", "EnalcMask_" + index_str, pedited, spot.enalcMask, spotEdited.enalcMask);
                assignFromKeyfile(keyFile, "Locallab", "Blendmasklc_" + index_str, pedited, spot.blendmasklc, spotEdited.blendmasklc);
                assignFromKeyfile(keyFile, "Locallab", "Radmasklc_" + index_str, pedited, spot.radmasklc, spotEdited.radmasklc);
                assignFromKeyfile(keyFile, "Locallab", "Chromasklc_" + index_str, pedited, spot.chromasklc, spotEdited.chromasklc);
                assignFromKeyfile(keyFile, "Locallab", "LmasklcCurve_" + index_str, pedited, spot.Lmasklccurve, spotEdited.Lmasklccurve);
                // Contrast by detail levels
                spot.visicbdl = assignFromKeyfile(keyFile, "Locallab", "Expcbdl_" + index_str, pedited, spot.expcbdl, spotEdited.expcbdl);

                if (spot.visicbdl) {
                    spotEdited.visicbdl = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Complexcbdl_" + index_str, pedited, spot.complexcbdl, spotEdited.complexcbdl);

                for (int j = 0; j < 6; j ++) {
                    assignFromKeyfile(keyFile, "Locallab", "Mult" + std::to_string(j) + "_" + index_str, pedited, spot.mult[j], spotEdited.mult[j]);
                }

                assignFromKeyfile(keyFile, "Locallab", "Chromacbdl_" + index_str, pedited, spot.chromacbdl, spotEdited.chromacbdl);
                assignFromKeyfile(keyFile, "Locallab", "Threshold_" + index_str, pedited, spot.threshold, spotEdited.threshold);
                assignFromKeyfile(keyFile, "Locallab", "Sensicb_" + index_str, pedited, spot.sensicb, spotEdited.sensicb);
                assignFromKeyfile(keyFile, "Locallab", "Clarityml_" + index_str, pedited, spot.clarityml, spotEdited.clarityml);
                assignFromKeyfile(keyFile, "Locallab", "Contresid_" + index_str, pedited, spot.contresid, spotEdited.contresid);
                assignFromKeyfile(keyFile, "Locallab", "Softradiuscb_" + index_str, pedited, spot.softradiuscb, spotEdited.softradiuscb);
                assignFromKeyfile(keyFile, "Locallab", "EnacbMask_" + index_str, pedited, spot.enacbMask, spotEdited.enacbMask);
                assignFromKeyfile(keyFile, "Locallab", "CCmaskcbCurve_" + index_str, pedited, spot.CCmaskcbcurve, spotEdited.CCmaskcbcurve);
                assignFromKeyfile(keyFile, "Locallab", "LLmaskcbCurve_" + index_str, pedited, spot.LLmaskcbcurve, spotEdited.LLmaskcbcurve);
                assignFromKeyfile(keyFile, "Locallab", "HHmaskcbCurve_" + index_str, pedited, spot.HHmaskcbcurve, spotEdited.HHmaskcbcurve);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskcb_" + index_str, pedited, spot.blendmaskcb, spotEdited.blendmaskcb);
                assignFromKeyfile(keyFile, "Locallab", "Radmaskcb_" + index_str, pedited, spot.radmaskcb, spotEdited.radmaskcb);
                assignFromKeyfile(keyFile, "Locallab", "Chromaskcb_" + index_str, pedited, spot.chromaskcb, spotEdited.chromaskcb);
                assignFromKeyfile(keyFile, "Locallab", "Gammaskcb_" + index_str, pedited, spot.gammaskcb, spotEdited.gammaskcb);
                assignFromKeyfile(keyFile, "Locallab", "Slomaskcb_" + index_str, pedited, spot.slomaskcb, spotEdited.slomaskcb);
                assignFromKeyfile(keyFile, "Locallab", "Lapmaskcb_" + index_str, pedited, spot.lapmaskcb, spotEdited.lapmaskcb);
                assignFromKeyfile(keyFile, "Locallab", "LmaskcbCurve_" + index_str, pedited, spot.Lmaskcbcurve, spotEdited.Lmaskcbcurve);
                // Log encoding
                spot.visilog = assignFromKeyfile(keyFile, "Locallab", "Explog_" + index_str, pedited, spot.explog, spotEdited.explog);

                if (spot.visilog) {
                    spotEdited.visilog = true;
                }

                assignFromKeyfile(keyFile, "Locallab", "Autocompute_" + index_str, pedited, spot.autocompute, spotEdited.autocompute);
                assignFromKeyfile(keyFile, "Locallab", "SourceGray_" + index_str, pedited, spot.sourceGray, spotEdited.sourceGray);
                assignFromKeyfile(keyFile, "Locallab", "TargetGray_" + index_str, pedited, spot.targetGray, spotEdited.targetGray);
                assignFromKeyfile(keyFile, "Locallab", "AutoGray_" + index_str, pedited, spot.Autogray, spotEdited.Autogray);
                assignFromKeyfile(keyFile, "Locallab", "Fullimage_" + index_str, pedited, spot.fullimage, spotEdited.fullimage);
                assignFromKeyfile(keyFile, "Locallab", "BlackEv_" + index_str, pedited, spot.blackEv, spotEdited.blackEv);
                assignFromKeyfile(keyFile, "Locallab", "WhiteEv_" + index_str, pedited, spot.whiteEv, spotEdited.whiteEv);
                assignFromKeyfile(keyFile, "Locallab", "Detail_" + index_str, pedited, spot.detail, spotEdited.detail);
                assignFromKeyfile(keyFile, "Locallab", "Sensilog_" + index_str, pedited, spot.sensilog, spotEdited.sensilog);
                assignFromKeyfile(keyFile, "Locallab", "Baselog_" + index_str, pedited, spot.baselog, spotEdited.baselog);
                assignFromKeyfile(keyFile, "Locallab", "Strlog_" + index_str, pedited, spot.strlog, spotEdited.strlog);
                assignFromKeyfile(keyFile, "Locallab", "Anglog_" + index_str, pedited, spot.anglog, spotEdited.anglog);
                // mask
                spot.visimask = assignFromKeyfile(keyFile, "Locallab", "Expmask_" + index_str, pedited, spot.expmask, spotEdited.expmask);
                assignFromKeyfile(keyFile, "Locallab", "Complexmask_" + index_str, pedited, spot.complexmask, spotEdited.complexmask);
                assignFromKeyfile(keyFile, "Locallab", "Sensimask_" + index_str, pedited, spot.sensimask, spotEdited.sensimask);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskmask_" + index_str, pedited, spot.blendmask, spotEdited.blendmask);
                assignFromKeyfile(keyFile, "Locallab", "Blendmaskmaskab_" + index_str, pedited, spot.blendmaskab, spotEdited.blendmaskab);
                assignFromKeyfile(keyFile, "Locallab", "Softradiusmask_" + index_str, pedited, spot.softradiusmask, spotEdited.softradiusmask);
                assignFromKeyfile(keyFile, "Locallab", "Enamask_" + index_str, pedited, spot.enamask, spotEdited.enamask);
                assignFromKeyfile(keyFile, "Locallab", "Fftmask_" + index_str, pedited, spot.fftmask, spotEdited.fftmask);
                assignFromKeyfile(keyFile, "Locallab", "Blurmask_" + index_str, pedited, spot.blurmask, spotEdited.blurmask);
                assignFromKeyfile(keyFile, "Locallab", "Contmask_" + index_str, pedited, spot.contmask, spotEdited.contmask);
                assignFromKeyfile(keyFile, "Locallab", "CCmask_Curve_" + index_str, pedited, spot.CCmask_curve, spotEdited.CCmask_curve);
                assignFromKeyfile(keyFile, "Locallab", "LLmask_Curve_" + index_str, pedited, spot.LLmask_curve, spotEdited.LLmask_curve);
                assignFromKeyfile(keyFile, "Locallab", "HHmask_Curve_" + index_str, pedited, spot.HHmask_curve, spotEdited.HHmask_curve);
                assignFromKeyfile(keyFile, "Locallab", "Strumaskmask_" + index_str, pedited, spot.strumaskmask, spotEdited.strumaskmask);
                assignFromKeyfile(keyFile, "Locallab", "Toolmask_" + index_str, pedited, spot.toolmask, spotEdited.toolmask);
                assignFromKeyfile(keyFile, "Locallab", "Radmask_" + index_str, pedited, spot.radmask, spotEdited.radmask);
                assignFromKeyfile(keyFile, "Locallab", "Lapmask_" + index_str, pedited, spot.lapmask, spotEdited.lapmask);
                assignFromKeyfile(keyFile, "Locallab", "Chromask_" + index_str, pedited, spot.chromask, spotEdited.chromask);
                assignFromKeyfile(keyFile, "Locallab", "Gammask_" + index_str, pedited, spot.gammask, spotEdited.gammask);
                assignFromKeyfile(keyFile, "Locallab", "Slopmask_" + index_str, pedited, spot.slopmask, spotEdited.slopmask);
                assignFromKeyfile(keyFile, "Locallab", "Shadmask_" + index_str, pedited, spot.shadmask, spotEdited.shadmask);
                assignFromKeyfile(keyFile, "Locallab", "Str_mask_" + index_str, pedited, spot.str_mask, spotEdited.str_mask);
                assignFromKeyfile(keyFile, "Locallab", "Ang_mask_" + index_str, pedited, spot.ang_mask, spotEdited.ang_mask);
                assignFromKeyfile(keyFile, "Locallab", "HHhmask_Curve_" + index_str, pedited, spot.HHhmask_curve, spotEdited.HHhmask_curve);
                assignFromKeyfile(keyFile, "Locallab", "Lmask_Curve_" + index_str, pedited, spot.Lmask_curve, spotEdited.Lmask_curve);
                assignFromKeyfile(keyFile, "Locallab", "LLmask_Curvewav_" + index_str, pedited, spot.LLmask_curvewav, spotEdited.LLmask_curvewav);

                if (keyFile.has_key("Locallab", "CSThresholdmask_" + index_str)) {
                    const std::vector<int> thresh = keyFile.get_integer_list("Locallab", "CSThresholdmask_" + index_str);

                    if (thresh.size() >= 4) {
                        spot.csthresholdmask.setValues(thresh[0], thresh[1], min(thresh[2], 10), min(thresh[3], 10));
                    }

                    spotEdited.csthresholdmask = true;
                }

                if (spot.visimask) {
                    spotEdited.visimask = true;
                }

                // Append LocallabSpot and LocallabParamsEdited
                locallab.spots.push_back(spot);

                if (pedited) {
                    pedited->locallab.spots.push_back(spotEdited);
                }

                // Update increment
                ++i;
            }
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
            if (ppVersion >= 339) {
                assignFromKeyfile(keyFile, "Resize", "AllowUpscaling", pedited, resize.allowUpscaling, pedited->resize.allowUpscaling);
            } else {
                resize.allowUpscaling = false;
                if (pedited) {
                    pedited->resize.allowUpscaling = true;
                }
            }
        }

        if (keyFile.has_group("PostDemosaicSharpening")) {
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "Enabled", pedited, pdsharpening.enabled, pedited->pdsharpening.enabled);
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "Contrast", pedited, pdsharpening.contrast, pedited->pdsharpening.contrast);
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "AutoContrast", pedited, pdsharpening.autoContrast, pedited->pdsharpening.autoContrast);
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "AutoRadius", pedited, pdsharpening.autoRadius, pedited->pdsharpening.autoRadius);
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "DeconvRadius", pedited, pdsharpening.deconvradius, pedited->pdsharpening.deconvradius);
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "DeconvRadiusOffset", pedited, pdsharpening.deconvradiusOffset, pedited->pdsharpening.deconvradiusOffset);
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "DeconvIterCheck", pedited, pdsharpening.deconvitercheck, pedited->pdsharpening.deconvitercheck);
            assignFromKeyfile(keyFile, "PostDemosaicSharpening", "DeconvIterations", pedited, pdsharpening.deconviter, pedited->pdsharpening.deconviter);
        }

        if (keyFile.has_group("PostResizeSharpening")) {
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Enabled", pedited, prsharpening.enabled, pedited->prsharpening.enabled);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Contrast", pedited, prsharpening.contrast, pedited->prsharpening.contrast);
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
                icm.inputProfile = expandRelativePath(fname, "file:", keyFile.get_string("Color Management", "InputProfile"));

                if (pedited) {
                    pedited->icm.inputProfile = true;
                }
            }

            assignFromKeyfile(keyFile, "Color Management", "ToneCurve", pedited, icm.toneCurve, pedited->icm.toneCurve);
            assignFromKeyfile(keyFile, "Color Management", "ApplyLookTable", pedited, icm.applyLookTable, pedited->icm.applyLookTable);
            assignFromKeyfile(keyFile, "Color Management", "ApplyBaselineExposureOffset", pedited, icm.applyBaselineExposureOffset, pedited->icm.applyBaselineExposureOffset);
            assignFromKeyfile(keyFile, "Color Management", "ApplyHueSatMap", pedited, icm.applyHueSatMap, pedited->icm.applyHueSatMap);
            assignFromKeyfile(keyFile, "Color Management", "DCPIlluminant", pedited, icm.dcpIlluminant, pedited->icm.dcpIlluminant);
            assignFromKeyfile(keyFile, "Color Management", "WorkingProfile", pedited, icm.workingProfile, pedited->icm.workingProfile);
            assignFromKeyfile(keyFile, "Color Management", "WorkingTRC", pedited, icm.workingTRC, pedited->icm.workingTRC);
            assignFromKeyfile(keyFile, "Color Management", "WorkingTRCGamma", pedited, icm.workingTRCGamma, pedited->icm.workingTRCGamma);
            assignFromKeyfile(keyFile, "Color Management", "WorkingTRCSlope", pedited, icm.workingTRCSlope, pedited->icm.workingTRCSlope);

            assignFromKeyfile(keyFile, "Color Management", "OutputProfile", pedited, icm.outputProfile, pedited->icm.outputProfile);
            if (ppVersion < 341) {
                if (icm.outputProfile == "RT_Medium_gsRGB") {
                    icm.outputProfile = "RTv4_Medium";
                } else if (icm.outputProfile == "RT_Large_gBT709" || icm.outputProfile == "RT_Large_g10" || icm.outputProfile == "RT_Large_gsRGB") {
                    icm.outputProfile = "RTv4_Large";
                } else if (icm.outputProfile == "WideGamutRGB") {
                    icm.outputProfile = "RTv4_Wide";
                } else if (icm.outputProfile == "RT_sRGB_gBT709" || icm.outputProfile == "RT_sRGB_g10" || icm.outputProfile == "RT_sRGB") {
                    icm.outputProfile = "RTv4_sRGB";
                } else if (icm.outputProfile == "BetaRGB") { // Have we ever provided this profile ? Should we convert this filename ?
                    icm.outputProfile = "RTv4_Beta";
                } else if (icm.outputProfile == "BestRGB") { // Have we ever provided this profile ? Should we convert this filename ?
                    icm.outputProfile = "RTv4_Best";
                } else if (icm.outputProfile == "Rec2020") {
                    icm.outputProfile = "RTv4_Rec2020";
                } else if (icm.outputProfile == "Bruce") { // Have we ever provided this profile ? Should we convert this filename ?
                    icm.outputProfile = "RTv4_Bruce";
                } else if (icm.outputProfile == "ACES") {
                    icm.outputProfile = "RTv4_ACES-AP0";
                }
            }
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
        }

        if (keyFile.has_group("Wavelet")) {
            assignFromKeyfile(keyFile, "Wavelet", "Enabled", pedited, wavelet.enabled, pedited->wavelet.enabled);
            assignFromKeyfile(keyFile, "Wavelet", "Strength", pedited, wavelet.strength, pedited->wavelet.strength);
            assignFromKeyfile(keyFile, "Wavelet", "Balance", pedited, wavelet.balance, pedited->wavelet.balance);
            assignFromKeyfile(keyFile, "Wavelet", "Sigmafin", pedited, wavelet.sigmafin, pedited->wavelet.sigmafin);
            assignFromKeyfile(keyFile, "Wavelet", "Sigmaton", pedited, wavelet.sigmaton, pedited->wavelet.sigmaton);
            assignFromKeyfile(keyFile, "Wavelet", "Sigmacol", pedited, wavelet.sigmacol, pedited->wavelet.sigmacol);
            assignFromKeyfile(keyFile, "Wavelet", "Sigmadir", pedited, wavelet.sigmadir, pedited->wavelet.sigmadir);
            assignFromKeyfile(keyFile, "Wavelet", "Rangeab", pedited, wavelet.rangeab, pedited->wavelet.rangeab);
            assignFromKeyfile(keyFile, "Wavelet", "Protab", pedited, wavelet.protab, pedited->wavelet.protab);
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
            assignFromKeyfile(keyFile, "Wavelet", "Ballum", pedited, wavelet.ballum, pedited->wavelet.ballum);
            assignFromKeyfile(keyFile, "Wavelet", "Balchrom", pedited, wavelet.balchrom, pedited->wavelet.balchrom);
            assignFromKeyfile(keyFile, "Wavelet", "Chromfine", pedited, wavelet.chromfi, pedited->wavelet.chromfi);
            assignFromKeyfile(keyFile, "Wavelet", "Chromcoarse", pedited, wavelet.chromco, pedited->wavelet.chromco);
            assignFromKeyfile(keyFile, "Wavelet", "MergeL", pedited, wavelet.mergeL, pedited->wavelet.mergeL);
            assignFromKeyfile(keyFile, "Wavelet", "MergeC", pedited, wavelet.mergeC, pedited->wavelet.mergeC);
            assignFromKeyfile(keyFile, "Wavelet", "Softrad", pedited, wavelet.softrad, pedited->wavelet.softrad);
            assignFromKeyfile(keyFile, "Wavelet", "Softradend", pedited, wavelet.softradend, pedited->wavelet.softradend);
            assignFromKeyfile(keyFile, "Wavelet", "Lipst", pedited, wavelet.lipst, pedited->wavelet.lipst);
            assignFromKeyfile(keyFile, "Wavelet", "AvoidColorShift", pedited, wavelet.avoid, pedited->wavelet.avoid);
            assignFromKeyfile(keyFile, "Wavelet", "Showmask", pedited, wavelet.showmask, pedited->wavelet.showmask);
            assignFromKeyfile(keyFile, "Wavelet", "Oldsh", pedited, wavelet.oldsh, pedited->wavelet.oldsh);
            assignFromKeyfile(keyFile, "Wavelet", "TMr", pedited, wavelet.tmr, pedited->wavelet.tmr);
            assignFromKeyfile(keyFile, "Wavelet", "LabGridALow", pedited, wavelet.labgridALow, pedited->wavelet.labgridALow);
            assignFromKeyfile(keyFile, "Wavelet", "LabGridBLow", pedited, wavelet.labgridBLow, pedited->wavelet.labgridBLow);
            assignFromKeyfile(keyFile, "Wavelet", "LabGridAHigh", pedited, wavelet.labgridAHigh, pedited->wavelet.labgridAHigh);
            assignFromKeyfile(keyFile, "Wavelet", "LabGridBHigh", pedited, wavelet.labgridBHigh, pedited->wavelet.labgridBHigh);

            if (ppVersion < 331) { // wavelet.Lmethod was a string before version 331
                Glib::ustring temp;
                assignFromKeyfile(keyFile, "Wavelet", "LevMethod", pedited, temp, pedited->wavelet.Lmethod);

                try {
                    wavelet.Lmethod = std::stoi(temp);
                } catch (...) {
                }
            } else {
                assignFromKeyfile(keyFile, "Wavelet", "LevMethod", pedited, wavelet.Lmethod, pedited->wavelet.Lmethod);
            }

            assignFromKeyfile(keyFile, "Wavelet", "ChoiceLevMethod", pedited, wavelet.CLmethod, pedited->wavelet.CLmethod);
            assignFromKeyfile(keyFile, "Wavelet", "BackMethod", pedited, wavelet.Backmethod, pedited->wavelet.Backmethod);
            assignFromKeyfile(keyFile, "Wavelet", "TilesMethod", pedited, wavelet.Tilesmethod, pedited->wavelet.Tilesmethod);
            assignFromKeyfile(keyFile, "Wavelet", "complexMethod", pedited, wavelet.complexmethod, pedited->wavelet.complexmethod);
            assignFromKeyfile(keyFile, "Wavelet", "DaubMethod", pedited, wavelet.daubcoeffmethod, pedited->wavelet.daubcoeffmethod);
            assignFromKeyfile(keyFile, "Wavelet", "CHromaMethod", pedited, wavelet.CHmethod, pedited->wavelet.CHmethod);
            assignFromKeyfile(keyFile, "Wavelet", "Medgreinf", pedited, wavelet.Medgreinf, pedited->wavelet.Medgreinf);
            assignFromKeyfile(keyFile, "Wavelet", "Ushamethod", pedited, wavelet.ushamethod, pedited->wavelet.ushamethod);
            assignFromKeyfile(keyFile, "Wavelet", "CHSLromaMethod", pedited, wavelet.CHSLmethod, pedited->wavelet.CHSLmethod);
            assignFromKeyfile(keyFile, "Wavelet", "EDMethod", pedited, wavelet.EDmethod, pedited->wavelet.EDmethod);
            assignFromKeyfile(keyFile, "Wavelet", "NPMethod", pedited, wavelet.NPmethod, pedited->wavelet.NPmethod);
            assignFromKeyfile(keyFile, "Wavelet", "BAMethod", pedited, wavelet.BAmethod, pedited->wavelet.BAmethod);
            assignFromKeyfile(keyFile, "Wavelet", "TMMethod", pedited, wavelet.TMmethod, pedited->wavelet.TMmethod);
            assignFromKeyfile(keyFile, "Wavelet", "HSMethod", pedited, wavelet.HSmethod, pedited->wavelet.HSmethod);
            assignFromKeyfile(keyFile, "Wavelet", "DirMethod", pedited, wavelet.Dirmethod, pedited->wavelet.Dirmethod);
            assignFromKeyfile(keyFile, "Wavelet", "Sigma", pedited, wavelet.sigma, pedited->wavelet.sigma);
            assignFromKeyfile(keyFile, "Wavelet", "Offset", pedited, wavelet.offset, pedited->wavelet.offset);
            assignFromKeyfile(keyFile, "Wavelet", "Lowthr", pedited, wavelet.lowthr, pedited->wavelet.lowthr);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualcontShadow", pedited, wavelet.rescon, pedited->wavelet.rescon);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualcontHighlight", pedited, wavelet.resconH, pedited->wavelet.resconH);
            assignFromKeyfile(keyFile, "Wavelet", "Residualchroma", pedited, wavelet.reschro, pedited->wavelet.reschro);
            assignFromKeyfile(keyFile, "Wavelet", "Residualblur", pedited, wavelet.resblur, pedited->wavelet.resblur);
            assignFromKeyfile(keyFile, "Wavelet", "Residualblurc", pedited, wavelet.resblurc, pedited->wavelet.resblurc);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualTM", pedited, wavelet.tmrs, pedited->wavelet.tmrs);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualEDGS", pedited, wavelet.edgs, pedited->wavelet.edgs);
            assignFromKeyfile(keyFile, "Wavelet", "ResidualSCALE", pedited, wavelet.scale, pedited->wavelet.scale);
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
            assignFromKeyfile(keyFile, "Wavelet", "Edgeffect", pedited, wavelet.edgeffect, pedited->wavelet.edgeffect);
            assignFromKeyfile(keyFile, "Wavelet", "Edgval", pedited, wavelet.edgval, pedited->wavelet.edgval);
            assignFromKeyfile(keyFile, "Wavelet", "ThrEdg", pedited, wavelet.edgthresh, pedited->wavelet.edgthresh);
            assignFromKeyfile(keyFile, "Wavelet", "ThresholdResidShadow", pedited, wavelet.thr, pedited->wavelet.thr);
            assignFromKeyfile(keyFile, "Wavelet", "ThresholdResidHighLight", pedited, wavelet.thrH, pedited->wavelet.thrH);
            assignFromKeyfile(keyFile, "Wavelet", "Residualradius", pedited, wavelet.radius, pedited->wavelet.radius);
            assignFromKeyfile(keyFile, "Wavelet", "ContrastCurve", pedited, wavelet.ccwcurve, pedited->wavelet.ccwcurve);
            assignFromKeyfile(keyFile, "Wavelet", "blcurve", pedited, wavelet.blcurve, pedited->wavelet.blcurve);
            assignFromKeyfile(keyFile, "Wavelet", "OpacityCurveRG", pedited, wavelet.opacityCurveRG, pedited->wavelet.opacityCurveRG);
            assignFromKeyfile(keyFile, "Wavelet", "Levalshc", pedited, wavelet.opacityCurveSH, pedited->wavelet.opacityCurveSH);
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
            assignFromKeyfile(keyFile, "Wavelet", "chrwav", pedited, wavelet.chrwav, pedited->wavelet.chrwav);
            assignFromKeyfile(keyFile, "Wavelet", "bluwav", pedited, wavelet.bluwav, pedited->wavelet.bluwav);
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
            assignFromKeyfile(keyFile, "Wavelet", "expbl", pedited, wavelet.expbl, pedited->wavelet.expbl);
            assignFromKeyfile(keyFile, "Wavelet", "Expresid", pedited, wavelet.expresid, pedited->wavelet.expresid);
            assignFromKeyfile(keyFile, "Wavelet", "Expfinal", pedited, wavelet.expfinal, pedited->wavelet.expfinal);
            assignFromKeyfile(keyFile, "Wavelet", "Exptoning", pedited, wavelet.exptoning, pedited->wavelet.exptoning);
            assignFromKeyfile(keyFile, "Wavelet", "Expnoise", pedited, wavelet.expnoise, pedited->wavelet.expnoise);
            assignFromKeyfile(keyFile, "Wavelet", "Expclari", pedited, wavelet.expclari, pedited->wavelet.expclari);
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

        if (keyFile.has_group("SoftLight")) {
            assignFromKeyfile(keyFile, "SoftLight", "Enabled", pedited, softlight.enabled, pedited->softlight.enabled);
            assignFromKeyfile(keyFile, "SoftLight", "Strength", pedited, softlight.strength, pedited->softlight.strength);
        }

        if (keyFile.has_group("Dehaze")) {
            assignFromKeyfile(keyFile, "Dehaze", "Enabled", pedited, dehaze.enabled, pedited->dehaze.enabled);
            assignFromKeyfile(keyFile, "Dehaze", "Strength", pedited, dehaze.strength, pedited->dehaze.strength);
            assignFromKeyfile(keyFile, "Dehaze", "ShowDepthMap", pedited, dehaze.showDepthMap, pedited->dehaze.showDepthMap);
            assignFromKeyfile(keyFile, "Dehaze", "Depth", pedited, dehaze.depth, pedited->dehaze.depth);
            assignFromKeyfile(keyFile, "Dehaze", "Luminance", pedited, dehaze.luminance, pedited->dehaze.luminance);
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
            if (ppVersion < 337) {
                const double scale = ColorToningParams::LABGRID_CORR_SCALE;
                colorToning.labgridALow *= scale;
                colorToning.labgridAHigh *= scale;
                colorToning.labgridBLow *= scale;
                colorToning.labgridBHigh *= scale;
            }
            std::vector<ColorToningParams::LabCorrectionRegion> lg;
            bool found = false;
            bool done = false;
            for (int i = 1; !done; ++i) {
                ColorToningParams::LabCorrectionRegion cur;
                done = true;
                std::string n = std::to_string(i);
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionA_") + n, pedited, cur.a, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionB_") + n, pedited, cur.b, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionSaturation_") + n, pedited, cur.saturation, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionSlope_") + n, pedited, cur.slope, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionOffset_") + n, pedited, cur.offset, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionPower_") + n, pedited, cur.power, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionHueMask_") + n, pedited, cur.hueMask, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionChromaticityMask_") + n, pedited, cur.chromaticityMask, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionLightnessMask_") + n, pedited, cur.lightnessMask, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionMaskBlur_") + n, pedited, cur.maskBlur, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "ColorToning", Glib::ustring("LabRegionChannel_") + n, pedited, cur.channel, pedited->colorToning.labregions)) {
                    found = true;
                    done = false;
                }
                if (!done) {
                    lg.emplace_back(cur);
                }
            }
            if (found) {
                colorToning.labregions = std::move(lg);
            }
            assignFromKeyfile(keyFile, "ColorToning", "LabRegionsShowMask", pedited, colorToning.labregionsShowMask, pedited->colorToning.labregionsShowMask);
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
            if (ppVersion >= 342) {
                assignFromKeyfile(keyFile, "RAW", "CAAutoIterations", pedited, raw.caautoiterations, pedited->raw.caautoiterations);
            } else {
                raw.caautoiterations = 1;
            }

            if (ppVersion >= 343) {
                assignFromKeyfile(keyFile, "RAW", "CAAvoidColourshift", pedited, raw.ca_avoidcolourshift, pedited->raw.ca_avoidcolourshift);
            } else {
                raw.ca_avoidcolourshift = false;
            }
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
            assignFromKeyfile(keyFile, "RAW Bayer", "Border", pedited, raw.bayersensor.border, pedited->raw.bayersensor.border);

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
            assignFromKeyfile(keyFile, "RAW Bayer", "DualDemosaicAutoContrast", pedited, raw.bayersensor.dualDemosaicAutoContrast, pedited->raw.bayersensor.dualDemosaicAutoContrast);
            if (ppVersion < 345) {
                raw.bayersensor.dualDemosaicAutoContrast = false;
                if (pedited) {
                    pedited->raw.bayersensor.dualDemosaicAutoContrast = true;
                }
            }
            assignFromKeyfile(keyFile, "RAW Bayer", "DualDemosaicContrast", pedited, raw.bayersensor.dualDemosaicContrast, pedited->raw.bayersensor.dualDemosaicContrast);

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
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftEqualBright", pedited, raw.bayersensor.pixelShiftEqualBright, pedited->raw.bayersensor.pixelShiftEqualBright);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftEqualBrightChannel", pedited, raw.bayersensor.pixelShiftEqualBrightChannel, pedited->raw.bayersensor.pixelShiftEqualBrightChannel);
            assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftNonGreenCross", pedited, raw.bayersensor.pixelShiftNonGreenCross, pedited->raw.bayersensor.pixelShiftNonGreenCross);

            if (ppVersion < 336) {
                if (keyFile.has_key("RAW Bayer", "pixelShiftLmmse")) {
                    const bool useLmmse = keyFile.get_boolean("RAW Bayer", "pixelShiftLmmse");

                    if (useLmmse) {
                        raw.bayersensor.pixelShiftDemosaicMethod = raw.bayersensor.getPSDemosaicMethodString(RAWParams::BayerSensor::PSDemosaicMethod::LMMSE);
                    } else {
                        raw.bayersensor.pixelShiftDemosaicMethod = raw.bayersensor.getPSDemosaicMethodString(RAWParams::BayerSensor::PSDemosaicMethod::AMAZE);
                    }

                    if (pedited) {
                        pedited->raw.bayersensor.pixelShiftDemosaicMethod = true;
                    }
                }
            } else {
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftDemosaicMethod", pedited, raw.bayersensor.pixelShiftDemosaicMethod, pedited->raw.bayersensor.pixelShiftDemosaicMethod);
            }

            assignFromKeyfile(keyFile, "RAW Bayer", "PDAFLinesFilter", pedited, raw.bayersensor.pdafLinesFilter, pedited->raw.bayersensor.pdafLinesFilter);
        }

        if (keyFile.has_group("RAW X-Trans")) {
            assignFromKeyfile(keyFile, "RAW X-Trans", "Method", pedited, raw.xtranssensor.method, pedited->raw.xtranssensor.method);
            assignFromKeyfile(keyFile, "RAW X-Trans", "DualDemosaicAutoContrast", pedited, raw.xtranssensor.dualDemosaicAutoContrast, pedited->raw.xtranssensor.dualDemosaicAutoContrast);
            if (ppVersion < 345) {
                raw.xtranssensor.dualDemosaicAutoContrast = false;
                if (pedited) {
                    pedited->raw.xtranssensor.dualDemosaicAutoContrast = true;
                }
            }
            assignFromKeyfile(keyFile, "RAW X-Trans", "DualDemosaicContrast", pedited, raw.xtranssensor.dualDemosaicContrast, pedited->raw.xtranssensor.dualDemosaicContrast);
            assignFromKeyfile(keyFile, "RAW X-Trans", "Border", pedited, raw.xtranssensor.border, pedited->raw.xtranssensor.border);
            assignFromKeyfile(keyFile, "RAW X-Trans", "CcSteps", pedited, raw.xtranssensor.ccSteps, pedited->raw.xtranssensor.ccSteps);
            assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackRed", pedited, raw.xtranssensor.blackred, pedited->raw.xtranssensor.exBlackRed);
            assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackGreen", pedited, raw.xtranssensor.blackgreen, pedited->raw.xtranssensor.exBlackGreen);
            assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackBlue", pedited, raw.xtranssensor.blackblue, pedited->raw.xtranssensor.exBlackBlue);
        }

        if (keyFile.has_group("Film Negative")) {
            assignFromKeyfile(keyFile, "Film Negative", "Enabled", pedited, filmNegative.enabled, pedited->filmNegative.enabled);
            assignFromKeyfile(keyFile, "Film Negative", "RedRatio", pedited, filmNegative.redRatio, pedited->filmNegative.redRatio);
            assignFromKeyfile(keyFile, "Film Negative", "GreenExponent", pedited, filmNegative.greenExp, pedited->filmNegative.greenExp);
            assignFromKeyfile(keyFile, "Film Negative", "BlueRatio", pedited, filmNegative.blueRatio, pedited->filmNegative.blueRatio);
            if (ppVersion >= 347) {
                bool r, g, b;
                assignFromKeyfile(keyFile, "Film Negative", "RedBase", pedited, filmNegative.redBase, r);
                assignFromKeyfile(keyFile, "Film Negative", "GreenBase", pedited, filmNegative.greenBase, g);
                assignFromKeyfile(keyFile, "Film Negative", "BlueBase", pedited, filmNegative.blueBase, b);
                if (pedited) {
                    pedited->filmNegative.baseValues = r || g || b;
                }
            } else {
                // Backwards compatibility with film negative in RT 5.7: use special film base value -1,
                // to signal that the old channel scaling method should be used.
                filmNegative.redBase = -1.f;
                filmNegative.greenBase = -1.f;
                filmNegative.blueBase = -1.f;
                if (pedited) {
                    pedited->filmNegative.baseValues = true;
                }
            }
        }

        if (keyFile.has_group("RAW Preprocess WB")) {
            if (keyFile.has_key("RAW Preprocess WB", "Mode")) {
                raw.preprocessWB.mode = RAWParams::PreprocessWB::Mode(keyFile.get_integer("RAW Preprocess WB", "Mode"));

                if (pedited) {
                    pedited->raw.preprocessWB.mode = true;
                }
            }
        }

        if (keyFile.has_group("MetaData")) {
            int mode = int(MetaDataParams::TUNNEL);
            assignFromKeyfile(keyFile, "MetaData", "Mode", pedited, mode, pedited->metadata.mode);

            if (mode >= int(MetaDataParams::TUNNEL) && mode <= int(MetaDataParams::STRIP)) {
                metadata.mode = static_cast<MetaDataParams::Mode>(mode);
            }
        }

        if (keyFile.has_group("Exif")) {
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
        && locallab == other.locallab
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
        && softlight == other.softlight
        && rgbCurves == other.rgbCurves
        && colorToning == other.colorToning
        && metadata == other.metadata
        && exif == other.exif
        && iptc == other.iptc
        && dehaze == other.dehaze
        && filmNegative == other.filmNegative;
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

void PartialProfile::applyTo(ProcParams* destParams, bool fromLastSave) const
{
    if (destParams && pparams && pedited) {
        bool fromHistMatching = fromLastSave && destParams->toneCurve.histmatching && pparams->toneCurve.histmatching;
        pedited->combine(*destParams, *pparams, true);
        if (!fromLastSave) {
            destParams->toneCurve.fromHistMatching = fromHistMatching;
        }
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
