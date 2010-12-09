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
#ifndef __RAWIMAGE_H
#define __RAWIMAGE_H

#include <time.h>
#include <glibmm.h>
#include <dcraw.h>

namespace rtengine {

struct badPix
{
	int x;
	int y;
	badPix(	int xc, int yc ):x(xc),y(yc){}
};

class RawImage: public DCraw
{
public:
  RawImage(  const Glib::ustring name );
  ~RawImage();

  int loadRaw (bool loadData=true, bool closeFile=true);
  int get_colorsCoeff( float *pre_mul, float *scale_mul, int *cblack  );
  void set_prefilters(){
      if (isBayer() && get_colors() == 3) {
         prefilters = filters;
  	     filters &= ~((filters & 0x55555555) << 1);
      }
  }
  dcrawImage_t get_image() { return image; }
  unsigned short** compress_image(); // revert to compressed pixels format and release image data
  unsigned short** data;             // holds pixel values, data[i][j] corresponds to the ith row and jth column
  unsigned prefilters;               // original filters saved ( used for 4 color processing )
protected:
  Glib::ustring filename; // complete filename
  int rotate_deg; // 0,90,180,270 degree of rotation: info taken by dcraw from exif
  char* profile_data; // Embedded ICC color profile
  unsigned short* allocation; // pointer to allocated memory

public:

  std::string get_filename() const { return filename;}
  int get_width()  const { return width; }
  int get_height() const { return height; }
  int get_FujiWidth() const { return fuji_width; }
  bool isBayer() const { return filters!=0; }
  unsigned get_filters() const { return filters; }
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
  float get_rgb_cam( int r, int c) const { return rgb_cam[r][c];}

  int get_exifBase()  const {return exif_base; }
  int get_ciffBase() const {return ciff_base; }
  int get_ciffLen()  const {return ciff_len; }

  int get_profileLen() const {return profile_length;}
  char* get_profile() const { return profile_data;}
  IMFILE *get_file() { return ifp; }
  bool is_supportedThumb() const ;
  int get_thumbOffset(){ return int(thumb_offset);}
  int get_thumbWidth(){ return int(thumb_width);}
  int get_thumbHeight(){ return int(thumb_height);}
  int get_thumbBPS(){ return thumb_load_raw ? 16 : 8; }
  bool get_thumbSwap() const;
  unsigned get_thumbLength(){ return thumb_length;}
public:
  // dcraw functions
  void scale_colors(){ DCraw::scale_colors(); }
  void pre_interpolate() { DCraw::pre_interpolate(); }

public:
  bool ISRED  (unsigned row, unsigned col) const { return ((filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==0);}
  bool ISGREEN(unsigned row, unsigned col) const { return ((filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==1);}
  bool ISBLUE (unsigned row, unsigned col) const { return ((filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==2);}
  unsigned FC (unsigned row, unsigned col) const { return (filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3); }
};

}

#endif // __RAWIMAGE_H
