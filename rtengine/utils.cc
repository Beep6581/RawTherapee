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
#include "../rtengine/utils.h"
#include <cmath>
#include <cstring>
#include <cstdio>

#undef MAX
#undef MIN

#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))

namespace rtengine {

void bilinearInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh) {

    int ix = 0;
    for (int i=0; i<dh; i++) {
        int sy = i*sh/dh;
        if (sy>=sh) sy = sh-1;
        double dy = (double)i*sh/dh - sy;
        int ny = sy+1;
        if (ny>=sh) ny = sy;
        int or1 = 3*sw*sy;
        int or2 = 3*sw*ny;
        for (int j=0; j<dw; j++) {
            int sx = j*sw/dw;
            if (sx>=sw) sx = sw;
            double dx = (double)j*sw/dw - sx;
            int nx = sx+1;
            if (nx>=sw) nx = sx;
            int ofs11 = or1 + 3*sx;
            int ofs12 = or1 + 3*nx;
            int ofs21 = or2 + 3*sx;
            int ofs22 = or2 + 3*nx;
            unsigned int val = src[ofs11]*(1-dx)*(1-dy) + src[ofs12]*dx*(1-dy) + src[ofs21]*(1-dx)*dy + src[ofs22]*dx*dy;
            dst[ix++] = val;
            ofs11++; ofs12++; ofs21++; ofs22++;
            val = src[ofs11]*(1-dx)*(1-dy) + src[ofs12]*dx*(1-dy) + src[ofs21]*(1-dx)*dy + src[ofs22]*dx*dy;
            dst[ix++] = val;
            ofs11++; ofs12++; ofs21++; ofs22++;
            val = src[ofs11]*(1-dx)*(1-dy) + src[ofs12]*dx*(1-dy) + src[ofs21]*(1-dx)*dy + src[ofs22]*dx*dy;
            dst[ix++] = val;
        }
    }
}

void nearestInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh) {

    int ix = 0;
    for (int i=0; i<dh; i++) {
        int rofs = sw * (i * sh / dh);
        for (int j=0; j<dw; j++) {
            int dx = rofs + j * sw / dw;
            dx *= 3;
            dst[ix++] = src[dx++];
            dst[ix++] = src[dx++];
            dst[ix++] = src[dx++];
        }
    }
}

void rotate (unsigned char* img, int& w, int& h, int deg) {

    if (deg==0) {
        return;
    }
    
    unsigned char* rotated = new unsigned char[3*w*h];
    int ix = 0;
    if (deg==90) {
        for (int i=0; i<h; i++) 
            for (int j=0; j<w; j++) {
                rotated[3*(j*h+h-i-1)+0] = img[ix++];
                rotated[3*(j*h+h-i-1)+1] = img[ix++];
                rotated[3*(j*h+h-i-1)+2] = img[ix++];
            }
        int tmp = w;
        w = h;
        h = tmp;
    }
    else if (deg==270) {
        for (int i=0; i<h; i++) 
            for (int j=0; j<w; j++) {
                rotated[3*(h*(w-j-1)+i)+0] = img[ix++];
                rotated[3*(h*(w-j-1)+i)+1] = img[ix++];
                rotated[3*(h*(w-j-1)+i)+2] = img[ix++];
            }
        int tmp = w;
        w = h;
        h = tmp;
    }
    else if (deg==180)
        for (int i=0; i<h; i++) 
            for (int j=0; j<w; j++) {
                rotated[3*(w*(h-i-1)+w-j-1)+0] = img[ix++];
                rotated[3*(w*(h-i-1)+w-j-1)+1] = img[ix++];
                rotated[3*(w*(h-i-1)+w-j-1)+2] = img[ix++];
            }
    memcpy (img, rotated, 3*w*h);
    delete [] rotated;
}

void hflip (unsigned char* img, int w, int h) {

    unsigned char* flipped = new unsigned char[3*w*h];
    int ix = 0;
    for (int i=0; i<h; i++)
        for (int j=0; j<w; j++) {
            flipped[3*(w*i+w-1-j)+0] = img[ix++];
            flipped[3*(w*i+w-1-j)+1] = img[ix++];
            flipped[3*(w*i+w-1-j)+2] = img[ix++];
        }
    memcpy (img, flipped, 3*w*h);
    delete [] flipped;
}

void vflip (unsigned char* img, int w, int h) {

    unsigned char* flipped = new unsigned char[3*w*h];
    int ix = 0;
    for (int i=0; i<h; i++)
        for (int j=0; j<w; j++) {
            flipped[3*(w*(h-1-i)+j)+0] = img[ix++];
            flipped[3*(w*(h-1-i)+j)+1] = img[ix++];
            flipped[3*(w*(h-1-i)+j)+2] = img[ix++];
        }
    memcpy (img, flipped, 3*w*h);
    delete [] flipped;
}

void rgb2hsv (int r, int g, int b, float &h, float &s, float &v) {

	double var_R = r / 65535.0;
	double var_G = g / 65535.0;
	double var_B = b / 65535.0;

	double var_Min = MIN(MIN(var_R,var_G),var_B);
	double var_Max = MAX(MAX(var_R,var_G),var_B);
	double del_Max = var_Max - var_Min;
	v = var_Max;
	if (fabs(del_Max)<0.00001) {
		h = 0;
		s = 0;
	}
	else {
		s = del_Max/var_Max;

		if      ( var_R == var_Max ) h = (var_G - var_B)/del_Max;
		else if ( var_G == var_Max ) h = 2.0 + (var_B - var_R)/del_Max;
		else if ( var_B == var_Max ) h = 4.0 + (var_R - var_G)/del_Max;
		h /= 6.0;

		if ( h < 0 )  h += 1;
		if ( h > 1 )  h -= 1;
	}
}

void hsv2rgb (float h, float s, float v, int &r, int &g, int &b) {

	float h1 = h*6; // sector 0 to 5
	int i = floor( h1 );
	float f = h1 - i; // fractional part of h

	float p = v * ( 1 - s );
	float q = v * ( 1 - s * f );
	float t = v * ( 1 - s * ( 1 - f ) );

	float r1,g1,b1;

	if (i==0) {r1 = v;  g1 = t;  b1 = p;}
	else if (i==1) {r1 = q;  g1 = v;  b1 = p;}
	else if (i==2) {r1 = p;  g1 = v;  b1 = t;}
	else if (i==3) {r1 = p;  g1 = q;  b1 = v;}
	else if (i==4) {r1 = t;  g1 = p;  b1 = v;}
	else if (i==5) {r1 = v;  g1 = p;  b1 = q;}

	r = (int)( r1 * 65535);
	g = (int)( g1 * 65535);
	b = (int)( b1 * 65535);
}

// The same function but set float values instead if int
// Function copied for speed concerns
void hsv2rgb (float h, float s, float v, float &r, float &g, float &b) {

	float h1 = h*6; // sector 0 to 5
	int i = floor( h1 );
	float f = h1 - i; // fractional part of h

	float p = v * ( 1 - s );
	float q = v * ( 1 - s * f );
	float t = v * ( 1 - s * ( 1 - f ) );

	if (i==0) {r = v;  g = t;  b = p;}
	else if (i==1) {r = q;  g = v;  b = p;}
	else if (i==2) {r = p;  g = v;  b = t;}
	else if (i==3) {r = p;  g = q;  b = v;}
	else if (i==4) {r = t;  g = p;  b = v;}
	else if (i==5) {r = v;  g = p;  b = q;}
}

}


