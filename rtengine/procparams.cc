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
#include <glib/gstdio.h>
#include "procparams.h"
#include <glibmm.h>
#include <sstream>
#include <string.h>

#include <safekeyfile.h>

namespace rtengine {

ProcParams defaultProcParams;

void ProcParams::setFloat  (const String& key, float value) {

	floatParams[key] = value;
}

float ProcParams::getFloat (const String& key) {

	if (floatParams.count (key))
		return floatParams[key];
	else
		return defaultProcParams.floatParams[key];
}

void ProcParams::setInteger (const String& key, int value) {

	intParams[key] = value;
}

int ProcParams::getInteger (const String& key) {

	if (intParams.count (key))
		return intParams[key];
	else
		return defaultProcParams.intParams[key];
}

void ProcParams::setBoolean (const String& key, bool value) {

	boolParams[key] = value;
}

bool ProcParams::getBoolean (const String& key) {

	if (boolParams.count (key))
		return boolParams[key];
	else
		return defaultProcParams.boolParams[key];
}

void ProcParams::setString  (const String& key, const String& value) {

	stringParams[key] = value;
}

String ProcParams::getString  (const String& key) {

	if (stringParams.count (key))
		return stringParams[key];
	else
		return defaultProcParams.stringParams[key];
}

void ProcParams::setFloatList (const String& key, const FloatList& value) {

	floatListParams[key] = value;
}

FloatList& ProcParams::getFloatList (const String& key) {

	if (floatListParams.count (key))
		return floatListParams[key];
	else
		return defaultProcParams.floatListParams[key];
}

void ProcParams::setIntegerList (const String& key, const IntList& value) {

	intListParams[key] = value;
}

IntList& ProcParams::getIntegerList (const String& key) {

	if (intListParams.count (key))
		return intListParams[key];
	else
		return defaultProcParams.intListParams[key];
}

void ProcParams::setStringList (const String& key, const StringList& value) {

	stringListParams[key] = value;
}

StringList& ProcParams::getStringList (const String& key) {

	if (stringListParams.count (key))
		return stringListParams[key];
	else
		return defaultProcParams.stringListParams[key];
}


ProcParams::ProcParams () { 

    setDefaults (); 
}       

ProcParams* ProcParams::create () {

    return new ProcParams();
}

void ProcParams::destroy (ProcParams* pp) {

    delete pp;
}

void ProcParams::setDefaults () {

}

int ProcParams::save (Glib::ustring fname) const {

	return 0;
}

int ProcParams::load (Glib::ustring fname) {

	return 0;
}

bool operator==(const ExifPair& a, const ExifPair& b) {

    return a.field == b.field && a.value == b.value;
}

bool operator==(const IPTCPair& a, const IPTCPair& b) {

    return a.field == b.field && a.values == b.values;
}
bool ProcParams::operator== (const ProcParams& other) {

	return false;
/*
    return 
           toneCurve.curve      == other.toneCurve.curve
        && toneCurve.brightness == other.toneCurve.brightness
        && toneCurve.black      == other.toneCurve.black
        && toneCurve.contrast   == other.toneCurve.contrast
        && toneCurve.shcompr    == other.toneCurve.shcompr
        && toneCurve.hlcompr    == other.toneCurve.hlcompr
        && toneCurve.autoexp    == other.toneCurve.autoexp
        && toneCurve.clip       == other.toneCurve.clip
        && toneCurve.expcomp    == other.toneCurve.expcomp
        && lumaCurve.curve      == other.lumaCurve.curve
        && lumaCurve.brightness == other.lumaCurve.brightness
        && lumaCurve.contrast   == other.lumaCurve.contrast
        && sharpening.enabled   == other.sharpening.enabled
        && sharpening.radius    == other.sharpening.radius
        && sharpening.amount    == other.sharpening.amount
        && sharpening.threshold     == other.sharpening.threshold
        && sharpening.edgesonly     == other.sharpening.edgesonly
        && sharpening.edges_radius  == other.sharpening.edges_radius
        && sharpening.edges_tolerance   == other.sharpening.edges_tolerance
        && sharpening.halocontrol   == other.sharpening.halocontrol
        && sharpening.halocontrol_amount== other.sharpening.halocontrol_amount
        && sharpening.method        == other.sharpening.method
        && sharpening.deconvamount  == other.sharpening.deconvamount
        && sharpening.deconvradius  == other.sharpening.deconvradius
        && sharpening.deconviter    == other.sharpening.deconviter
        && sharpening.deconvdamping == other.sharpening.deconvdamping
        && colorBoost.amount        == other.colorBoost.amount
        && colorBoost.avoidclip == other.colorBoost.avoidclip
        && colorBoost.enable_saturationlimiter == other.colorBoost.enable_saturationlimiter
        && colorBoost.saturationlimit == other.colorBoost.saturationlimit
        && wb.method        == other.wb.method
        && wb.green         == other.wb.green
        && wb.temperature   == other.wb.temperature
        && colorShift.a     == other.colorShift.a
        && colorShift.b     == other.colorShift.b
        && lumaDenoise.enabled      == other.lumaDenoise.enabled
        && lumaDenoise.radius       == other.lumaDenoise.radius
        && lumaDenoise.edgetolerance == other.lumaDenoise.edgetolerance
        && colorDenoise.enabled      == other.colorDenoise.enabled
        && colorDenoise.radius       == other.colorDenoise.radius
        && colorDenoise.edgetolerance == other.colorDenoise.edgetolerance
        && colorDenoise.edgesensitive == other.colorDenoise.edgesensitive
        && sh.enabled       == other.sh.enabled
        && sh.hq            == other.sh.hq
        && sh.highlights    == other.sh.highlights
        && sh.htonalwidth   == other.sh.htonalwidth
        && sh.shadows       == other.sh.shadows
        && sh.stonalwidth   == other.sh.stonalwidth
        && sh.localcontrast == other.sh.localcontrast
        && sh.radius        == other.sh.radius
        && crop.enabled == other.crop.enabled
        && crop.x       == other.crop.x
        && crop.y       == other.crop.y
        && crop.w       == other.crop.w
        && crop.h       == other.crop.h
        && crop.fixratio == other.crop.fixratio
        && crop.ratio   == other.crop.ratio
        && crop.orientation == other.crop.orientation
        && crop.guide   == other.crop.guide
        && coarse.rotate == other.coarse.rotate
        && coarse.hflip == other.coarse.hflip
        && coarse.vflip == other.coarse.vflip
        && rotate.degree == other.rotate.degree
        && commonTrans.autofill == other.commonTrans.autofill
        && distortion.uselensfun == other.distortion.uselensfun
        && distortion.amount == other.distortion.amount
        && perspective.horizontal == other.perspective.horizontal
        && perspective.vertical == other.perspective.vertical
        && cacorrection.red == other.cacorrection.red
        && cacorrection.blue == other.cacorrection.blue
        && vignetting.amount == other.vignetting.amount
        && vignetting.radius == other.vignetting.radius
        && !memcmp (&chmixer.red, &other.chmixer.red, 3*sizeof(int))
        && !memcmp (&chmixer.green, &other.chmixer.green, 3*sizeof(int))
        && !memcmp (&chmixer.blue, &other.chmixer.blue, 3*sizeof(int))
        && hlrecovery.enabled   == other.hlrecovery.enabled
        && hlrecovery.method    == other.hlrecovery.method
        && resize.scale     == other.resize.scale
        && resize.method    == other.resize.method
        && resize.dataspec  == other.resize.dataspec
        && resize.width     == other.resize.width
        && resize.height    == other.resize.height
        && icm.input        == other.icm.input
        && icm.gammaOnInput == other.icm.gammaOnInput
        && icm.working      == other.icm.working
        && icm.output       == other.icm.output
        && exif==other.exif
        && iptc==other.iptc;*/
}

bool ProcParams::operator!= (const ProcParams& other) {

    return !(*this==other);
}
}

