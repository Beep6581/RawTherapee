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
//
// A class representing a 16 bit rgb image with separate planes and 16 byte aligned data
//
#ifndef _IMAGE_
#define _IMAGE_

#include "rtcommon.h"
#include "finalimage16.h"
#include <lcms2.h>
#include <FreeImage.h>
#include <exiv2/exiv2.hpp>

namespace rtengine {

class Image : public FinalImage16 {

		Exiv2::ExifData exifData;
		Exiv2::IptcData iptcData;
		Exiv2::XmpData	xmpData;
	
		FIBITMAP* bitmap;
	
		FIBITMAP* convertTo24bpp ();
		void	  writeMetadata (const String& fname);	
				  Image (FIBITMAP* img);

	public:
	
		Image (int width, int height);
		~Image ();

        void getEmbeddedICCProfile (int& length, unsigned char*& pdata);
        void setEmbeddedICCProfile (int length, unsigned char* pdata);
        
        static Image* load (const String& fname, bool fast=false);
        int save (const String& fname);
        
        void setMetadata (const Exiv2::ExifData& ed, const Exiv2::IptcData& id, const Exiv2::XmpData& xd);

        // Implementing FinalImage16 interface
        int getWidth ();
        int getHeight ();
        int getScanLineSize ();
        unsigned char* getData ();
        int saveAsPNG  (const String& fname, PNGCompression compr = PNGZDefaultCompression, bool bps16=true);
        int saveAsJPEG (const String& fname, int quality = -1, JPEGSubSampling ss = JPEGSubSampling_420);
        int saveAsTIFF (const String& fname, TIFFCompression compr = TIFFLZWCompression, bool bps16=true);
};

};
#endif

