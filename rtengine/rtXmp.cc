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

#include "rtXmp.h"

namespace rtengine {

const char *kXmpProcessing="Processing";
const char *kXmpSnapshotId="SnapshotId";
const char *kXmpSnapshot="SnapshotName";
const char *kXmpQueued="Queued";
const char *kXmpSaved="Saved";
const char *kXmpOutput="OutputFile";
const char *kXmpVersion="Version";
const char *kXmpEnabled="Enabled";
const char *kXmpMerged="XmpMerged";

bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, bool &var)
{
	Exiv2::XmpData::const_iterator iter = xmpData.findKey(Exiv2::XmpKey(key));
	if( iter == xmpData.end()){
		return false;
	}
	try{
		var = (iter->getValue()->toString().compare("True")==0);
	}
	catch( Exiv2::AnyError &e){
		return false;
	}
	return true;
}

bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, int &var)
{
	Exiv2::XmpData::const_iterator iter = xmpData.findKey(Exiv2::XmpKey(key));
	if( iter == xmpData.end()){
		return false;
	}
	try{
		var = int(iter->getValue()->toLong());
	}
	catch( Exiv2::AnyError &e){
		return false;
	}
	return true;
}

bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, double &var)
{
	Exiv2::XmpData::const_iterator iter = xmpData.findKey(Exiv2::XmpKey(key));
	if( iter == xmpData.end()){
		return false;
	}
	try{
		var = double(iter->getValue()->toFloat());
	}
	catch( Exiv2::AnyError &e){
		return false;
	}
	return true;
}

bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, float &var)
{
	Exiv2::XmpData::const_iterator iter = xmpData.findKey(Exiv2::XmpKey(key));
	if( iter == xmpData.end()){
		return false;
	}
	try{
		var = iter->getValue()->toFloat();
	}
	catch( Exiv2::AnyError &e){
		return false;
	}
	return true;
}

bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, std::string &var)
{
	Exiv2::XmpData::const_iterator iter = xmpData.findKey(Exiv2::XmpKey(key));
	if( iter == xmpData.end()){
		return false;
	}
	try{
		var = iter->getValue()->toString();
	}
	catch( Exiv2::AnyError &e){
		return false;
	}
	return true;
}

bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, Glib::ustring &var)
{
	Exiv2::XmpData::const_iterator iter = xmpData.findKey(Exiv2::XmpKey(key));
	if( iter == xmpData.end()){
		return false;
	}
	try{
		var = iter->getValue()->toString();
	}
	catch( Exiv2::AnyError &e){
		return false;
	}
	return true;
}

bool readIPTCFromXmp( const Exiv2::XmpData &xmpData, IPTCMeta &tag, std::vector< Glib::ustring> &v )
{
	Exiv2::XmpData::const_iterator iter = xmpData.findKey(Exiv2::XmpKey( tag.getXmpKey() ));
	if( iter == xmpData.end()){
		return false;
	}
	v.clear();
	switch( tag.arrType ){
	case Exiv2::xmpText:
		v.push_back( iter->getValue()->toString() );
		break;
	case Exiv2::xmpBag:
	case Exiv2::xmpSeq:
	{
		for( int i=0; i<iter->getValue()->count();i++ )
			v.push_back( iter->getValue()->toString(i) );
		break;
	}
	case Exiv2::langAlt:
	{
		Exiv2::LangAltValue x(iter->getValue()->toString());
		v.push_back( x.toString("x-default") );
		break;
	}
	default:
		return false;
	}
	return true;
}

}

