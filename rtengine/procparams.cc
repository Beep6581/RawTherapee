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
#include <fstream>
#include <string.h>
#include <exiv2/exiv2.hpp>

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

void ProcParams::setDefaults () {

	floatParams = defaultProcParams.floatParams;
	intParams = defaultProcParams.intParams;
	boolParams = defaultProcParams.boolParams;
	stringParams = defaultProcParams.stringParams;
	floatListParams = defaultProcParams.floatListParams;
	intListParams = defaultProcParams.intListParams;
	stringListParams = defaultProcParams.stringListParams;
}

int ProcParams::save (const String& fname) const {
	
	try {
		// create xmp data and register our namespace
		Exiv2::XmpData xmpData;

		// save float parameters. An "F" is appended to the end of the keys to indicate that these are floats.
		for (std::map<String,float>::const_iterator i=floatParams.begin(); i!=floatParams.end(); i++)
			xmpData[Glib::ustring::compose("Xmp.rt.%1F", i->first)] = Glib::ustring::format (i->second);

		// save integer parameters. An "I" is appended to the end of the keys to indicate that these are integers.
		for (std::map<String,int>::const_iterator i=intParams.begin(); i!=intParams.end(); i++)
			xmpData[Glib::ustring::compose("Xmp.rt.%1I", i->first)] = Glib::ustring::format (i->second);
			
		// save boolean parameters. A "B" is appended to the end of the keys to indicate that these are booleans.
		for (std::map<String,bool>::const_iterator i=boolParams.begin(); i!=boolParams.end(); i++)
			xmpData[Glib::ustring::compose("Xmp.rt.%1B", i->first)] = i->second ? "true" : "false";

		// save string parameters. An "S" is appended to the end of the keys to indicate that these are strings.
		for (std::map<String,String>::const_iterator i=stringParams.begin(); i!=stringParams.end(); i++)
			xmpData[Glib::ustring::compose("Xmp.rt.%1S", i->first)] = i->second;

		// save float list parameters. An "FL" is appended to the end of the keys to indicate that these are floatlists.
		for (std::map<String,FloatList>::const_iterator i=floatListParams.begin(); i!=floatListParams.end(); i++) {
			Exiv2::Value::AutoPtr arr = Exiv2::Value::create (Exiv2::xmpSeq);
			for (int j=0; j<i->second.size(); j++)
				arr->read (Glib::ustring::format (i->second[j]));
			xmpData.add (Exiv2::XmpKey (Glib::ustring::compose("Xmp.rt.%1FL", i->first)), arr.get());
		}

		// save int list parameters. An "IL" is appended to the end of the keys to indicate that these are intlists.
		for (std::map<String,IntList>::const_iterator i=intListParams.begin(); i!=intListParams.end(); i++) {
			Exiv2::Value::AutoPtr arr = Exiv2::Value::create (Exiv2::xmpSeq);
			for (int j=0; j<i->second.size(); j++)
				arr->read (Glib::ustring::format (i->second[j]));
			xmpData.add (Exiv2::XmpKey (Glib::ustring::compose("Xmp.rt.%1IL", i->first)), arr.get());
		}

		// save string list parameters. An "SL" is appended to the end of the keys to indicate that these are stringlists.
		for (std::map<String,StringList>::const_iterator i=stringListParams.begin(); i!=stringListParams.end(); i++) {
			Exiv2::Value::AutoPtr arr = Exiv2::Value::create (Exiv2::xmpSeq);
			for (int j=0; j<i->second.size(); j++)
				arr->read (i->second[j]);
			xmpData.add (Exiv2::XmpKey (Glib::ustring::compose("Xmp.rt.%1SL", i->first)), arr.get());
		}
		
		// sort and create xmp packet
		xmpData.sortByKey ();
		std::string xmpPacket;
		if (Exiv2::XmpParser::encode(xmpPacket, xmpData)) 
			return 1;

		// save to file
		std::ofstream f (fname.c_str ());
		if (f.is_open ()) {
			f << xmpPacket;
			f.close ();
		}
		else
			return 2;
	}
	catch (Exiv2::AnyError& e) {
		return 3;
	}
	return 0;
}

int ProcParams::load (const String& fname) {

	try {
		// open file and read content
		std::ifstream f (fname.c_str ());
		std::string xmpPacket;
		if (f.is_open ()) {
			xmpPacket = std::string (std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
			f.close ();
		}
		else
			return 2;
			
		// create xmp data and register our namespace
		Exiv2::XmpData xmpData;
		if (Exiv2::XmpParser::decode(xmpData, xmpPacket)) 
			return 1;
			
		// read all values
		for (Exiv2::XmpData::const_iterator i=xmpData.begin(); i!=xmpData.end(); i++) {
			std::string key = i->key ();
			if (key.size()>=9) {
				// find out type of the item
				char last = key[key.size()-1];
				std::string typeID;
				std::string propName;
				if (last=='L') {
					typeID = std::string (key, key.size()-2, 2);
					propName = std::string (key, 7, key.size()-9);
				}
				else {
					typeID = std::string (&last, 1);
					propName = std::string (key, 7, key.size()-8);
				}
				// store values
				if (typeID == "F")
					setFloat (propName, i->toFloat());
				else if (typeID == "I")
					setInteger (propName, i->toLong());
				else if (typeID == "S")
					setString (propName, i->toString());
				else if (typeID == "B")
					setBoolean (propName, i->toString()=="true");
				else if (typeID == "FL") {
					FloatList arr;
					arr.resize (i->count ());
					for (int j=0; j<arr.size(); j++)
						arr[j] = i->toFloat (j);
					setFloatList (propName, arr);
				}
				else if (typeID == "IL") {
					IntList arr;
					arr.resize (i->count ());
					for (int j=0; j<arr.size(); j++)
						arr[j] = i->toLong (j);
					setIntegerList (propName, arr);
				}
				else if (typeID == "SL") {
					StringList arr;
					arr.resize (i->count ());
					for (int j=0; j<arr.size(); j++)
						arr[j] = i->toString (j);
					setStringList (propName, arr);
				}
			}
		}
	}
	catch (Exiv2::AnyError& e) {
		return 3;
	}
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

