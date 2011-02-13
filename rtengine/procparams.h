/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2011 Gabor Horvath <hgabor@rawtherapee.com>
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
#ifndef _PROCPARAMS_H_
#define _PROCPARAMS_H_

#include <vector>
#include "rtcommon.h"
#include <map>
#include <vector>
#include <exiv2/exiv2.hpp>

namespace rtengine {

/**
  * A class representing a key/value for the exif metadata information
  */
class ExifPair {

    public:
		String field;
		String value;
};

/**
  * The IPTC key/value pairs
  */
class IPTCPair {

    public:
		String field;
        StringList values;
};

/**
  * This class holds all the processing parameters applied on the images
  */

class ProcParams {

		// to be replaced by QHash later
		std::map<String,float>   floatParams;
		std::map<String,int> 	 intParams;
		std::map<String,bool> 	 boolParams;
		std::map<String,String> stringParams;
		std::map<String,FloatVector> floatVectorParams;
		std::map<String,IntVector> intVectorParams;
		std::map<String,StringList> stringListParams;

		String getKey   (const String& skey) const;
		String getGroup (const String& skey) const;
		String removeQualifier (const String& skey) const;

    public:

		void   setFloat   (const String& group, const String& key, float value);
		float  getFloat   (const String& group, const String& key);

		void   setInteger (const String& group, const String& key, int value);
		int    getInteger (const String& group, const String& key);

		void   setBoolean (const String& group, const String& key, bool value);
		bool   getBoolean (const String& group, const String& key);

		void   setString  (const String& group, const String& key, const String& value);
		String getString  (const String& group, const String& key);

		void   		 setFloatVector (const String& group, const String& key, const FloatVector& value);
		FloatVector& getFloatVector (const String& group, const String& key);

		void   		setIntegerVector (const String& group, const String& key, const IntVector& value);
		IntVector& 	getIntegerVector (const String& group, const String& key);

		void   		setStringList (const String& group, const String& key, const StringList& value);
		StringList& getStringList (const String& group, const String& key);

		// --------------8<------------------ to be removed when all filters are rewritten ------

        std::vector<ExifPair> exif;             ///< List of modifications appplied on the exif tags of the input image
        std::vector<IPTCPair> iptc;             ///< The IPTC tags and values to be saved to the output image
        int version;                            ///< Version of the file from which the parameters have been read
        
      /**
        * The constructor only sets the hand-wired defaults.
        */
        ProcParams          ();
      /**
        * Sets the hand-wired defaults parameters.
        */
        void    setDefaults ();
      /**
        * Saves the parameters to a file.
        * @param fname the name of the file
        * @return Error code (=0 if no error)
        */
        int     save        (const String& fname) const;
      /**
        * Loads the parameters from a file.
        * @param fname the name of the file
        * @return Error code (=0 if no error)
        */
        int     load        (const String& fname);

		void 	addToXmp 	(Exiv2::XmpData& xmpData) const;

        bool operator== (const ProcParams& other);
        bool operator!= (const ProcParams& other);
};

extern ProcParams defaultProcParams;
}
#endif
