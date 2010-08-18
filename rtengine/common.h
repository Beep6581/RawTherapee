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
#ifndef _COMMON_
#define _COMMON_

#define ISRED(image,row,col) \
	((image->filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==0)
#define ISGREEN(image,row,col) \
	((image->filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==1)
#define ISBLUE(image,row,col) \
	((image->filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==2)

#define FISRED(filter,row,col) \
	((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==0 || !filter)
#define FISGREEN(filter,row,col) \
	((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==1 || !filter)
#define FISBLUE(filter,row,col) \
	((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==2 || !filter)

#define CMAXVAL 65535

#include <time.h>
#include <glibmm.h>

namespace rtengine {

struct badPix
{
	int x;
	int y;
	badPix(	int xc, int yc ):x(xc),y(yc){}
};

struct RawImage {

  Glib::ustring fname; // complete filename
  int width;  // with of the image as reported by dcraw
  int height; // height of the image as reported by dcraw

  unsigned filters; // sequence of Bayer filter colors: 2bit for each of 2x8 pixels grid indicate 0=Red,1=Green1,2=Blue,(3=green2)
  int colors; // Number of colors of bayer filter (3 or 4)

  int black_point; // Black offset taken from dslr info by dcraw
  int cblack[4]; // Black for each color.
  unsigned short white[8][8]; // square of white registered by camera
  float cam_mul[4]; // Camera color multiplier taken from exif by dcraw
  float pre_mul[4];
  int maximum; // White (maximum) point taken from dslr info
  int rotate_deg; // 0,90,180,270 degree of rotation: info taken by dcraw from exif
  int fuji_width;
  
  double defgain;
  double iso_speed;
  double shutter;
  double aperture;
  double focal_len;
  time_t timestamp;
  char *make, *model;

  int exifbase, prefilters, ciff_base, ciff_len;
  int thumbLength;
  int thumbOffset;
  int thumbType;
  int thumbWidth;
  int thumbHeight;

  unsigned short* allocation;
  unsigned short** data;             // holds pixel values, data[i][j] corresponds to the ith row and jth column

  float coeff[3][3];
  float icoeff[3][3];
  
  int profile_len;
  char* profile_data; // Embedded ICC color profile

  RawImage(  const Glib::ustring name):allocation(NULL),data(NULL),profile_data(NULL),fname(name)
  {
  }
  ~RawImage()
  {
	  if(allocation){ delete [] allocation; allocation=NULL;}
	  if(data){ delete [] data; data=NULL;}
	  if(profile_data){ delete [] profile_data; profile_data=NULL;}
  }

  int loadRaw (bool loadData=true);

  void allocData()
  {
	  if (filters) {
		if (!allocation) {
				allocation = new unsigned short[height * width];
				data = new unsigned short*[height];
				for (int i = 0; i < height; i++)
					data[i] = allocation + i * width;
		}
	  }else{
		if (!allocation) {
				allocation = new unsigned short[3 * height * width];
				data = new unsigned short*[height];
				for (int i = 0; i < height; i++)
					data[i] = allocation + 3 * i * width;
		}
	  }
	  if(profile_len)
		  profile_data = new char[profile_len];
  }

};

}

#endif
