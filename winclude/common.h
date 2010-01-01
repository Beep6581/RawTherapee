#ifndef _COMMON_
#define _COMMON_

#define ISRED(image,row,col) \
	((image->filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==0)
#define ISGREEN(image,row,col) \
	((image->filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==1)
#define ISBLUE(image,row,col) \
	((image->filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==2)


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

  struct tm* time;
  float iso_speed, aperture, focal_len, shutter;
  char *make, *model;

  int exifbase, exiflocation, exiforder;

  unsigned short** data;             // holds pixel values, data[i][j] corresponds to the ith row and jth column

  float coeff[3][4];
  float icoeff[3][4];

};

#endif
