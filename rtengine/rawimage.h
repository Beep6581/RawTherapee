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
#ifndef _RAWIMAGE_H__
#define _RAWIMAGE_H__

#include "rtcommon.h"
#include "matrix33.h"
#include "colortemp.h"

namespace rtengine {

class RawImage {

	public:

		int width;					// width of the image
		int height;					// height of the image

		int fujiWidth;			

		unsigned filter;			// Bayer pattern

		ColorTemp camSpaceTemp;		// color temperature in camera color space
		ColorTemp rgbSpaceTemp;		// color temperature in srgb D50 color space ??? (Warning! These 2 member variables do not make sense. I'm not good in color theory, however this does what it should for sure.
  
		float defgain;				// this is the multiplier by which each pixel values shall be multiplied to get the image without HL recovery

		unsigned short* allocation;	// the array that holds pixel values
		unsigned short** data;      // pointer array pointing to the rows

		Matrix33 cam_srgb;			// transformation matrix from the camera color space to the srgb D50 color space
		Matrix33 srgb_cam;			// inverse of the transformation matrix
  
		int profileLength;			// length of the embedded profile
		char* profileData;			// the icc profile embedded to the raw file

	public:

		RawImage ();
		~RawImage ();

		int load (const String& fname);

		inline bool isRed (int row, int col) {
			return (filter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==0;
		}
		inline bool isGreen (int row, int col) {
			return (filter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==1;
		}
		inline bool isBlue (int row, int col) {
			return (filter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==2;
		}
};
}
#endif
