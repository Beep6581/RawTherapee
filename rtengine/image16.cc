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
#include <image16.h>
#include <imagefloat.h>
#include <image8.h>
#include <string.h>
#include <stdio.h>
#include <rtengine.h>

using namespace rtengine;

Image16::Image16 () 
  : unaligned (NULL), data (NULL), r (NULL), g (NULL), b (NULL){
}

Image16::Image16 (int w, int h) 
  : unaligned (NULL), width(w), height (h), data (NULL), r (NULL), g (NULL), b (NULL) {

    allocate (w, h);
}

Image16::~Image16 () {
  
    if (data!=NULL) {
        delete [] data;
        delete [] r;
        delete [] g;
        delete [] b;
    }
}

void Image16::allocate (int W, int H) {

	width=W;
	height=H;
    if (data!=NULL) {
        delete [] data;
        delete [] r;
        delete [] g;
        delete [] b;
    }

    r = new unsigned short*[height];
    g = new unsigned short*[height];
    b = new unsigned short*[height];
    data = new unsigned short[W*H*3];

    rowstride = W;
    planestride = rowstride * height;

    unsigned short* redstart   = data + 0*planestride;
    unsigned short* greenstart = data + 1*planestride;
    unsigned short* bluestart  = data + 2*planestride;
       

    for (int i=0; i<height; i++) {
        r[i] = redstart   + i*rowstride;
        g[i] = greenstart + i*rowstride;
        b[i] = bluestart  + i*rowstride;
    }


}

void Image16::getScanline (int row, unsigned char* buffer, int bps) {

    if (data==NULL)
        return;
        
    if (bps==16) {
        int ix = 0;
        unsigned short* sbuffer = (unsigned short*) buffer;
        for (int i=0; i<width; i++) {
            sbuffer[ix++] = r[row][i];
            sbuffer[ix++] = g[row][i];
            sbuffer[ix++] = b[row][i];
        }
    }
    else if (bps==8) {
        int ix = 0;
        for (int i=0; i<width; i++) {
            buffer[ix++] = r[row][i] >> 8;
            buffer[ix++] = g[row][i] >> 8;
            buffer[ix++] = b[row][i] >> 8;
        }
    }        
}

void Image16::setScanline (int row, unsigned char* buffer, int bps) {

    if (data==NULL)
        return;
        
    if (bps==8) {
        int ix = 0;
        for (int i=0; i<width; i++) {
            r[row][i] = buffer[ix++] << 8;
            g[row][i] = buffer[ix++] << 8;
            b[row][i] = buffer[ix++] << 8;
        }
    }
    else if (bps==16) {
        unsigned short* sbuffer = (unsigned short*) buffer;
        int ix = 0;
        for (int i=0; i<width; i++) {
            r[row][i] = sbuffer[ix++];
            g[row][i] = sbuffer[ix++];
            b[row][i] = sbuffer[ix++];
        }
    }
}

Image16* Image16::copy () {

  Image16* cp = new Image16 (width, height);

  for (int i=0; i<height; i++) {
    memcpy (cp->r[i], r[i], width*sizeof(unsigned short));
    memcpy (cp->g[i], g[i], width*sizeof(unsigned short));
    memcpy (cp->b[i], b[i], width*sizeof(unsigned short));
  }

  return cp;
}

Image16* Image16::rotate (int deg) {

  if (deg==90) {
    Image16* result = new Image16 (height, width);
    for (int i=0; i<width; i++)
      for (int j=0; j<height; j++) {
        result->r[i][j] = r[height-1-j][i];
        result->g[i][j] = g[height-1-j][i];
        result->b[i][j] = b[height-1-j][i];
      }
    return result;
  }
  else if (deg==270) {
    Image16* result = new Image16 (height, width);
    for (int i=0; i<width; i++)
      for (int j=0; j<height; j++) {
        result->r[i][j] = r[j][width-1-i];
        result->g[i][j] = g[j][width-1-i];
        result->b[i][j] = b[j][width-1-i];
      }
    return result;
  }
  else if (deg==180) {
    Image16* result = new Image16 (width, height);
    for (int i=0; i<height; i++)
      for (int j=0; j<width; j++) {
        result->r[i][j] = r[height-1-i][width-1-j];
        result->g[i][j] = g[height-1-i][width-1-j];
        result->b[i][j] = b[height-1-i][width-1-j];
      }
    return result;
  }
  else
    return NULL;
}

Image16* Image16::hflip () {
 
  Image16* result = new Image16 (width, height);
  for (int i=0; i<height; i++)
    for (int j=0; j<width; j++) {
      result->r[i][j] = r[i][width-1-j];
      result->g[i][j] = g[i][width-1-j];
      result->b[i][j] = b[i][width-1-j];
    }
  return result;

}

Image16* Image16::vflip () {

  Image16* result = new Image16 (width, height);
  for (int i=0; i<height; i++)
    for (int j=0; j<width; j++) {
      result->r[i][j] = r[height-1-i][j];
      result->g[i][j] = g[height-1-i][j];
      result->b[i][j] = b[height-1-i][j];
    }
  return result;

}

Image16* Image16::resize (int nw, int nh, TypeInterpolation interp) {

    if (interp == TI_Nearest) {
        Image16* res = new Image16 (nw, nh);
        for (int i=0; i<nh; i++) {
            int ri = i*height/nh;
            for (int j=0; j<nw; j++) {
                int ci = j*width/nw;
                res->r[i][j] = r[ri][ci];
                res->g[i][j] = g[ri][ci];
                res->b[i][j] = b[ri][ci];
            }
        }
        return res;
    }
    else if (interp == TI_Bilinear) {
        Image16* res = new Image16 (nw, nh);
        for (int i=0; i<nh; i++) {
            int sy = i*height/nh;
            if (sy>=height) sy = height-1;
            double dy = (double)i*height/nh - sy;
            int ny = sy+1;
            if (ny>=height) ny = sy;
            for (int j=0; j<nw; j++) {
                int sx = j*width/nw;
                if (sx>=width) sx = width;
                double dx = (double)j*width/nw - sx;
                int nx = sx+1;
                if (nx>=width) nx = sx;
                res->r[i][j] = r[sy][sx]*(1-dx)*(1-dy) + r[sy][nx]*dx*(1-dy) + r[ny][sx]*(1-dx)*dy + r[ny][nx]*dx*dy;
                res->g[i][j] = g[sy][sx]*(1-dx)*(1-dy) + g[sy][nx]*dx*(1-dy) + g[ny][sx]*(1-dx)*dy + g[ny][nx]*dx*dy;
                res->b[i][j] = b[sy][sx]*(1-dx)*(1-dy) + b[sy][nx]*dx*(1-dy) + b[ny][sx]*(1-dx)*dy + b[ny][nx]*dx*dy;
            }
        }
        return res;
    }
    return NULL;
}

Image8* 
Image16::to8() const
{
	Image8* img8 = new Image8(width,height);
	for ( int h = 0; h < height; ++h )
	{
		for ( int w = 0; w < width; ++w )
		{
			img8->r(h,w,r[h][w] >> 8);
			img8->g(h,w,g[h][w] >> 8);
			img8->b(h,w,b[h][w] >> 8);
		}
	}
	return img8;
}

Imagefloat* 
Image16::tofloat() const
{
	Imagefloat* imgfloat = new Imagefloat(width,height);
	for ( int h = 0; h < height; ++h )
	{
		for ( int w = 0; w < width; ++w )
		{
			imgfloat->r[h][w] = ((float)r[h][w]) ;
			imgfloat->g[h][w] = ((float)g[h][w]) ;
			imgfloat->b[h][w] = ((float)b[h][w]) ;
		}
	}
	return imgfloat;
}

// Parallized transformation; create transform with cmsFLAGS_NOCACHE!
void Image16::ExecCMSTransform(cmsHTRANSFORM hTransform) {
    #pragma omp parallel for
    for (int i=0; i<height; i++)
        cmsDoTransform(hTransform, data + 3*i*rowstride, data + 3*i*rowstride, rowstride);
}