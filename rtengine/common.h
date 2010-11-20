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

  unsigned filters; // \TODO problem in pre_interpolate()
  int prefilters;

  unsigned short** data;             // holds pixel values, data[i][j] corresponds to the ith row and jth column

protected:
  Glib::ustring filename; // complete filename
  int width;  // with of the image as reported by dcraw
  int height; // height of the image as reported by dcraw
  int fuji_width;
  int colors; // Number of colors of bayer filter (3 or 4)
  int black; // Black offset taken from dslr info by dcraw
  int cblack[4]; // Black for each color.
  int maximum; // White (maximum) point taken from dslr info
  int rotate_deg; // 0,90,180,270 degree of rotation: info taken by dcraw from exif
  double iso_speed;
  double shutter;
  time_t timestamp;
  char *make, *model;
  unsigned short white[8][8]; // square of white registered by camera
  float cam_mul[4]; // Camera color multiplier taken from exif by dcraw
  float pre_mul[4];
  float coeff[3][3]; // rgb_cam coeff.

  int exifbase, ciff_base, ciff_len;
  int profile_len;
  char* profile_data; // Embedded ICC color profile

  unsigned short* allocation; // pointer to allocated memory
public:
  RawImage(  const Glib::ustring name);
  ~RawImage();
  int loadRaw (bool loadData=true);
  void allocData();
  std::string get_filename() const { return filename;}
  int get_width()  const { return width; }
  int get_height() const { return height; }
  int get_FujiWidth() const { return fuji_width; }
  bool isBayer() const { return filters!=0; }
  int get_colors() const { return colors;}
  int get_black()  const { return black;}
  int get_cblack(int i) const {return cblack[i];}
  int get_white() const { return maximum;}
  unsigned short get_whiteSample( int r, int c ) const { return white[r][c];}

  double get_ISOspeed() const {return iso_speed;}
  double get_shutter()  const {return shutter; }
  time_t get_timestamp() const { return timestamp;}
  int get_rotateDegree() const { return rotate_deg;}
  const std::string get_maker() const { return std::string(make); }
  const std::string get_model() const { return std::string(model); }

  float get_cam_mul(int c )const {return cam_mul[c];}
  float get_pre_mul(int c )const {return pre_mul[c];}
  float get_rgb_cam( int r, int c) const { return coeff[r][c];}

  int get_exifBase()  const {return exifbase; }
  int get_ciffBase() const {return ciff_base; }
  int get_ciffLen()  const {return ciff_len; }

  int get_profileLen() const {return profile_len;}
  char* get_profile() const { return profile_data;}

public:
  bool ISRED  (unsigned row, unsigned col) const { return ((filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==0);}
  bool ISGREEN(unsigned row, unsigned col) const { return ((filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==1);}
  bool ISBLUE (unsigned row, unsigned col) const { return ((filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==2);}
};

}

#endif
