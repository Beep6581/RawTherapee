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
#ifndef __IMAGEDATA_H__
#define __IMAGEDATA_H__

#include <time.h>
#include <string>
#include "rtengine.h"
#include <exiv2/exiv2.hpp>

namespace rtengine {

class ImageData : public ImageMetaData {

  protected:
	Exiv2::ExifData exifData;
	Exiv2::IptcData iptcData;
	Exiv2::XmpData  xmpData;

    struct tm time;
    int iso;
    float fNumber;
    float focalLen;
    int defRot;
    Exiv2::Rational expTime;
    
    std::string make, model;
    std::string lens;

    void extractInfo ();
    
  public:

    ImageData (const String& fname);

	const Exiv2::ExifData& 	getExifData () const { return exifData; }
	const Exiv2::IptcData& 	getIptcData () const { return iptcData; }
	const Exiv2::XmpData& 	getXmpData () const  { return xmpData; }

    struct tm   	getDateTime () const { return time;      }
    int         	getISO 		() const { return iso; }
    float       	getFNumber  () const { return fNumber;  }
    float       	getFocalLen () const { return focalLen;   }
    Exiv2::Rational getExposureTime () const { return expTime;   }
    int				getDefaultRotation () const { return defRot; }

    std::string getMake     () const { return make;      }
    std::string getModel    () const { return model;     }
    std::string getLens     () const { return lens;      }
    
};
};
#endif
