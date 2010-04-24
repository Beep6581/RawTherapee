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
#include <utils.h>
#include <string.h>
#include <stdio.h>

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

}


