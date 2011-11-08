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
 #ifndef _RTXMP_
 #define _RTXMP_
 
#include <glibmm.h>
#include <vector>
#include <exiv2/exiv2.hpp>
#include <iptcmeta.h>
 
namespace rtengine {

extern const char *kXmpProcessing;
extern const char *kXmpSnapshotId;
extern const char *kXmpSnapshot;
extern const char *kXmpQueued;
extern const char *kXmpSaved;
extern const char *kXmpOutput;
extern const char *kXmpVersion;
extern const char *kXmpEnabled;
extern const char *kXmpMerged;

template< class T > std::string serializeVector( const std::vector<T> &v)
{
	std::ostringstream os;
	for( size_t i = 0; i < v.size(); i++)
		os << v[i]<<";";
	return os.str();
}

template< class T > std::string serializeArray( T* v,size_t n)
{
	std::ostringstream os;
	for( size_t i = 0; i < n; i++)
		os << v[i]<<";";
	return os.str();
}
 
/*
 *  This is a quick workaround: reading different data types should be integrated in XMP class
 */
template <class T> int readVarFromXmp(Exiv2::XmpData &xmpData, const std::string& key, std::vector<T> &v )
{
	if( xmpData.findKey(Exiv2::XmpKey(key)) == xmpData.end())
		return 1;
	v.clear();
	std::istringstream s( xmpData[key].value().toString() );
	while( !s.eof()){
		char separator;
		T val;
	    s >> val;
	    if( !s.eof() ){
	    	s >> separator;
	    	v.push_back(val);
	    }
	}
	return 0;
}

template <class T> int readVarFromXmp(Exiv2::XmpData &xmpData, const std::string& key, T *v,size_t n )
{
	if( xmpData.findKey(Exiv2::XmpKey(key)) == xmpData.end())
		return 1;
	std::istringstream s( xmpData[key].value().toString() );
	size_t i=0;
	while( !s.eof() && i<n){
		char separator;
		T val;
	    s >> val;
	    if( !s.eof() ){
	    	s >> separator;
	    	v[i++]=val;
	    }
	}
	return 0;
}

bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, bool &var);
bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, int &var);
bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, double &var);
bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, std::string &var);
bool readVarFromXmp( const Exiv2::XmpData &xmpData, const Glib::ustring &key, Glib::ustring &var);
bool readIPTCFromXmp( const Exiv2::XmpData &xmpData, IPTCMeta &tag, std::vector< Glib::ustring> &v );

} 
#endif
