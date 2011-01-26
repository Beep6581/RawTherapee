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

	defaultProcParams.setBoolean ("CropEnabled", false);
	defaultProcParams.setInteger ("CropRectX", 0);
	defaultProcParams.setInteger ("CropRectY", 0);
	defaultProcParams.setInteger ("CropRectW", 100000);
	defaultProcParams.setInteger ("CropRectH", 100000);
	defaultProcParams.setBoolean ("CropFixRectRatio", true);
	defaultProcParams.setString  ("CropRectRatio", "3:2");
	defaultProcParams.setString  ("CropRectOrientation", "Landscape");
	defaultProcParams.setString  ("CropGuide", "None");

	defaultProcParams.setString  ("ColorManagementInputProfile", "");
	defaultProcParams.setBoolean ("ColorManagementGammaOnInput", false);
	defaultProcParams.setString  ("ColorManagementWorkingProfile", "sRGB");
	defaultProcParams.setString  ("ColorManagementOutputProfile", "sRGB");

    setDefaults (); 
}       

ProcParams* ProcParams::create () {

    return new ProcParams();
}

void ProcParams::destroy (ProcParams* pp) {

    delete pp;
}

void ProcParams::setDefaults () {

	floatParams = defaultProcParams.floatParams;
	intParams = defaultProcParams.intParams;
	boolParams = defaultProcParams.boolParams;
	stringParams = defaultProcParams.stringParams;
	floatListParams = defaultProcParams.floatListParams;
	intListParams = defaultProcParams.intListParams;
	stringListParams = defaultProcParams.stringListParams;
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

	return floatParams == other.floatParams
		&& intParams == other.intParams
		&& boolParams == other.boolParams
		&& stringParams == other.stringParams
		&& floatListParams == other.floatListParams
		&& intListParams == other.intListParams
		&& stringListParams == other.stringListParams;
}

bool ProcParams::operator!= (const ProcParams& other) {

    return !(*this==other);
}
}

