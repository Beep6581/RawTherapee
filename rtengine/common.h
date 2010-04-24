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

struct RawImage {

  int width;
  int height;

  unsigned filters;

  double red_multiplier;
  double green_multiplier;
  double blue_multiplier;

  double camwb_red;
  double camwb_green;
  double camwb_blue;

  int blackpoint;
  int rgb_max;
  int rotate_deg;
  int fuji_width;
  
  double defgain;

  char *make, *model;

  int exifbase, prefilters, ciff_base, ciff_len;

  unsigned short* allocation;
  unsigned short** data;             // holds pixel values, data[i][j] corresponds to the ith row and jth column

  float coeff[3][3];
  float icoeff[3][3];
  
  int profile_len;
  char* profile_data;
};

#endif
