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
#include "procparams.h"
#include <string.h>

#ifdef QTBUILD
#include <QFile>
#include <QTextStream>
#define Compose2(a,b,c) 	(QString(a).arg(b).arg(c)).toStdString()
#define Number2String(a)	QString::number(a).toStdString()
#else
#include <fstream>
#define Compose2(a,b,c) 	Glib::ustring::compose(a, b, c)
#define Number2String(a)	Glib::ustring::format(a)
#endif

namespace rtengine {

ProcParams defaultProcParams;

void ProcParams::setFloat  (const String& group, const String& key, float value) {

	String key_ = group + '/' + key;
	floatParams[key_] = value;
}

float ProcParams::getFloat (const String& group, const String& key) {

	String key_ = group + '/' + key;
	if (floatParams.count (key_))
		return floatParams[key_];
	else
		return defaultProcParams.floatParams[key_];
}

void ProcParams::setInteger (const String& group, const String& key, int value) {

	String key_ = group + '/' + key;
	intParams[key_] = value;
}

int ProcParams::getInteger (const String& group, const String& key) {

	String key_ = group + '/' + key;
	if (intParams.count (key_))
		return intParams[key_];
	else
		return defaultProcParams.intParams[key_];
}

void ProcParams::setBoolean (const String& group, const String& key, bool value) {

	String key_ = group + '/' + key;
	boolParams[key_] = value;
}

bool ProcParams::getBoolean (const String& group, const String& key) {

	String key_ = group + '/' + key;
	if (boolParams.count (key_))
		return boolParams[key_];
	else
		return defaultProcParams.boolParams[key_];
}

void ProcParams::setString  (const String& group, const String& key, const String& value) {

	String key_ = group + '/' + key;
	stringParams[key_] = value;
}

String ProcParams::getString  (const String& group, const String& key) {

	String key_ = group + '/' + key;
	if (stringParams.count (key_))
		return stringParams[key_];
	else
		return defaultProcParams.stringParams[key_];
}

void ProcParams::setFloatVector (const String& group, const String& key, const FloatVector& value) {

	String key_ = group + '/' + key;
	floatVectorParams[key_] = value;
}

FloatVector& ProcParams::getFloatVector (const String& group, const String& key) {

	String key_ = group + '/' + key;
	if (floatVectorParams.count (key_))
		return floatVectorParams[key_];
	else
		return defaultProcParams.floatVectorParams[key_];
}

void ProcParams::setIntegerVector (const String& group, const String& key, const IntVector& value) {

	String key_ = group + '/' + key;
	intVectorParams[key_] = value;
}

IntVector& ProcParams::getIntegerVector (const String& group, const String& key) {

	String key_ = group + '/' + key;
	if (intVectorParams.count (key_))
		return intVectorParams[key_];
	else
		return defaultProcParams.intVectorParams[key_];
}

void ProcParams::setStringList (const String& group, const String& key, const StringList& value) {

	String key_ = group + '/' + key;
	stringListParams[key_] = value;
}

StringList& ProcParams::getStringList (const String& group, const String& key) {

	String key_ = group + '/' + key;
	if (stringListParams.count (key_))
		return stringListParams[key_];
	else
		return defaultProcParams.stringListParams[key_];
}


ProcParams::ProcParams () { 

	defaultProcParams.setBoolean ("Crop", "Enabled", false);
	defaultProcParams.setInteger ("Crop", "RectX", 0);
	defaultProcParams.setInteger ("Crop", "RectY", 0);
	defaultProcParams.setInteger ("Crop", "RectW", 100000);
	defaultProcParams.setInteger ("Crop", "RectH", 100000);
	defaultProcParams.setBoolean ("Crop", "FixRectRatio", true);
	defaultProcParams.setString  ("Crop", "RectRatio", "3:2");
	defaultProcParams.setString  ("Crop", "RectOrientation", "Landscape");
	defaultProcParams.setString  ("Crop", "Guide", "None");

	defaultProcParams.setString  ("ColorManagement", "InputProfile", "");
	defaultProcParams.setBoolean ("ColorManagement", "GammaOnInput", false);
	defaultProcParams.setString  ("ColorManagement", "WorkingProfile", "sRGB");
	defaultProcParams.setString  ("ColorManagement", "OutputProfile", "sRGB");

    setDefaults (); 
}       

void ProcParams::setDefaults () {

	floatParams = defaultProcParams.floatParams;
	intParams = defaultProcParams.intParams;
	boolParams = defaultProcParams.boolParams;
	stringParams = defaultProcParams.stringParams;
	floatVectorParams = defaultProcParams.floatVectorParams;
	intVectorParams = defaultProcParams.intVectorParams;
	stringListParams = defaultProcParams.stringListParams;
}

String ProcParams::getKey (const String& skey) const {

#ifdef QTBUILD
	int pos = skey.indexOf ('/');
	if (pos>=0)
		return skey.mid (pos+1);
	else
		return String ("");
#else
	int pos = skey.find ('/');
	if (pos!=String::npos)
		return skey.substr (pos+1);
	else
		return String ("");
#endif
}

String ProcParams::getGroup (const String& skey) const{

#ifdef QTBUILD
	int pos = skey.indexOf ('/');
	if (pos>=0)
		return skey.left (pos);
	else
		return String ("");
#else
	int pos = skey.find ('/');
	if (pos!=String::npos)
		return skey.substr (0, pos);
	else
		return String ("");
#endif
}

String ProcParams::removeQualifier (const String& skey) const {

	if (skey.size() > 3)
#ifdef QTBUILD
		return skey.mid (3);
#else
	return skey.substr (3);
#endif
	else
		return String ("");
}

void ProcParams::addToXmp (Exiv2::XmpData& xmpData) const {

	try {
		// save float paramereplaceters. An "F" is appended to the end of the keys to indicate that these are floats.
		for (std::map<String,float>::const_iterator i=floatParams.begin(); i!=floatParams.end(); i++)
			xmpData[Compose2("Xmp.rt.%1/rt:%2F", getGroup(i->first), getKey(i->first))] = Number2String(i->second);

		// save integer parameters. An "I" is appended to the end of the keys to indicate that these are integers.
		for (std::map<String,int>::const_iterator i=intParams.begin(); i!=intParams.end(); i++)
			xmpData[Compose2("Xmp.rt.%1/rt:%2I", getGroup(i->first), getKey(i->first))] = Number2String(i->second);
			
		// save boolean parameters. A "B" is appended to the end of the keys to indicate that these are booleans.
		for (std::map<String,bool>::const_iterator i=boolParams.begin(); i!=boolParams.end(); i++)
			xmpData[Compose2("Xmp.rt.%1/rt:%2B", getGroup(i->first), getKey(i->first))] = i->second ? "true" : "false";

		// save string parameters. An "S" is appended to the end of the keys to indicate that these are strings.
		for (std::map<String,String>::const_iterator i=stringParams.begin(); i!=stringParams.end(); i++)
			xmpData[Compose2("Xmp.rt.%1/rt:%2S", getGroup(i->first), getKey(i->first))] = String2StdString(i->second);

		// save float list parameters. An "FL" is appended to the end of the keys to indicate that these are floatlists.
		for (std::map<String,FloatVector>::const_iterator i=floatVectorParams.begin(); i!=floatVectorParams.end(); i++) {
			Exiv2::Value::AutoPtr arr = Exiv2::Value::create (Exiv2::xmpSeq);
			for (int j=0; j<i->second.size(); j++)
				arr->read (Number2String(i->second[j]));
			xmpData.add (Exiv2::XmpKey (Compose2("Xmp.rt.%1/rt:%2FL", getGroup(i->first), getKey(i->first))), arr.get());
		}

		// save int list parameters. An "IL" is appended to the end of the keys to indicate that these are intlists.
		for (std::map<String,IntVector>::const_iterator i=intVectorParams.begin(); i!=intVectorParams.end(); i++) {
			Exiv2::Value::AutoPtr arr = Exiv2::Value::create (Exiv2::xmpSeq);
			for (int j=0; j<i->second.size(); j++)
				arr->read (Number2String(i->second[j]));
			xmpData.add (Exiv2::XmpKey (Compose2("Xmp.rt.%1/rt:%2IL", getGroup(i->first), getKey(i->first))), arr.get());
		}

		// save string list parameters. An "SL" is appended to the end of the keys to indicate that these are stringlists.
		for (std::map<String,StringList>::const_iterator i=stringListParams.begin(); i!=stringListParams.end(); i++) {
			Exiv2::Value::AutoPtr arr = Exiv2::Value::create (Exiv2::xmpSeq);
			for (int j=0; j<i->second.size(); j++)
				arr->read (String2StdString(i->second[j]));
			xmpData.add (Exiv2::XmpKey (Compose2("Xmp.rt.%1/rt:%2SL", getGroup(i->first), getKey(i->first))), arr.get());
		}
	}
	catch (Exiv2::AnyError& e) {
	}
}


int ProcParams::save (const String& fname) const {
	
	try {
		// create xmp data
		Exiv2::XmpData xmpData;
		addToXmp (xmpData);
		
		// sort and create xmp packet
		xmpData.sortByKey ();
		std::string xmpPacket;
		if (Exiv2::XmpParser::encode(xmpPacket, xmpData)) 
			return 1;

		// save to file
#ifdef QTBUILD
		QFile f (fname);
		 if (f.open (QFile::WriteOnly | QFile::Truncate)) {
		     QTextStream out(&f);
		     out << StdString2String(xmpPacket);
		     f.close ();
		 }
		 else
			 return 2;
#else
		std::ofstream f (fname.c_str ());
		if (f.is_open ()) {
			f << xmpPacket;
			f.close ();
		}
		else
			return 2;
#endif
	}
	catch (Exiv2::AnyError& e) {
		return 3;
	}
	return 0;
}

int ProcParams::load (const String& fname) {

	try {
		// open file and read content
		std::string xmpPacket;
#ifdef QTBUILD
		QFile f (fname);
		 if (f.open (QFile::ReadOnly)) {
		     QTextStream in(&f);
		     xmpPacket = String2StdString(in.readAll());
		     f.close ();
		 }
		 else
			 return 2;
#else
		std::ifstream f (fname.c_str ());
		if (f.is_open ()) {
			xmpPacket = std::string (std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
			f.close ();
		}
		else
			return 2;
#endif
			
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
				String propName;
				if (last=='L') {
					typeID = std::string (key, key.size()-2, 2);
					propName = StdString2String(std::string (key, 7, key.size()-9));
				}
				else {
					typeID = std::string (&last, 1);
					propName = StdString2String(std::string (key, 7, key.size()-8));
				}
				// store values
				if (typeID == "F")
					setFloat (getGroup (propName), removeQualifier (getKey (propName)), i->toFloat());
				else if (typeID == "I")
					setInteger (getGroup (propName), removeQualifier (getKey (propName)), i->toLong());
				else if (typeID == "S")
					setString (getGroup (propName), removeQualifier (getKey (propName)), StdString2String(i->toString()));
				else if (typeID == "B")
					setBoolean (getGroup (propName), removeQualifier (getKey (propName)), i->toString()=="true");
				else if (typeID == "FL") {
					FloatVector arr;
					for (int j=0; j<arr.size(); j++)
						arr.push_back (i->toFloat (j));
					setFloatVector (getGroup (propName), removeQualifier (getKey (propName)), arr);
				}
				else if (typeID == "IL") {
					IntVector arr;
					for (int j=0; j<arr.size(); j++)
						arr.push_back (i->toLong (j));
					setIntegerVector (getGroup (propName), removeQualifier (getKey (propName)), arr);
				}
				else if (typeID == "SL") {
					StringList arr;
					for (int j=0; j<arr.size(); j++)
						arr.push_back (StdString2String(i->toString (j)));
					setStringList (getGroup (propName), removeQualifier (getKey (propName)), arr);
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
		&& floatVectorParams == other.floatVectorParams
		&& intVectorParams == other.intVectorParams
		&& stringListParams == other.stringListParams;
}

bool ProcParams::operator!= (const ProcParams& other) {

    return !(*this==other);
}
}

