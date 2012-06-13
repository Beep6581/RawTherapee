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
#include "imagefloat.h"
#include "image16.h"
#include "image8.h"
#include <cstring>
#include "rtengine.h"
#include "mytime.h"
#include "iccstore.h"

using namespace rtengine;

Imagefloat::Imagefloat () 
  : unaligned (NULL), data (NULL), r (NULL), g (NULL), b (NULL){
}

Imagefloat::Imagefloat (int w, int h) 
  : unaligned (NULL), width(w), height (h), data (NULL), r (NULL), g (NULL), b (NULL) {

    allocate (w, h);
}

Imagefloat::~Imagefloat () {
  
    if (data!=NULL) {
        delete [] data;
        delete [] r;
        delete [] g;
        delete [] b;
    }
}

void Imagefloat::allocate (int W, int H) {

	width=W;
	height=H;

	if (data!=NULL) {
        delete [] data;
        delete [] r;
        delete [] g;
        delete [] b;
    }

    /*
    int lsize  = width + 8 - width % 8;
    unaligned = new unsigned char[16 + 3 * lsize * sizeof(float) * height];
    memset(unaligned, 0, (16 + 3 * lsize * sizeof(float) * height) * sizeof(unsigned char));

    uintptr_t poin = (uintptr_t)unaligned + 16 - (uintptr_t)unaligned % 16;
    data = (float*) (poin);
	*/
    r = new float*[height];
    g = new float*[height];
    b = new float*[height];

    data = new float[W*H*3];
    rowstride = W;
    planestride = rowstride * H;

    float * redstart   = data + 0*planestride;
    float * greenstart = data + 1*planestride;
    float * bluestart  = data + 2*planestride;


    for (int i=0; i<height; i++) {
        r[i] = (redstart   + i*rowstride);
        g[i] = (greenstart + i*rowstride);
        b[i] = (bluestart  + i*rowstride);
    }


}



Imagefloat* Imagefloat::copy () {

  Imagefloat* cp = new Imagefloat (width, height);

  for (int i=0; i<height; i++) {
    memcpy (cp->r[i], r[i], width*sizeof(float));
    memcpy (cp->g[i], g[i], width*sizeof(float));
    memcpy (cp->b[i], b[i], width*sizeof(float));
  }

  return cp;
}

Imagefloat* Imagefloat::rotate (int deg) {

  if (deg==90) {
    Imagefloat* result = new Imagefloat (height, width);
    for (int i=0; i<width; i++)
      for (int j=0; j<height; j++) {
        result->r[i][j] = r[height-1-j][i];
        result->g[i][j] = g[height-1-j][i];
        result->b[i][j] = b[height-1-j][i];
      }
    return result;
  }
  else if (deg==270) {
    Imagefloat* result = new Imagefloat (height, width);
    for (int i=0; i<width; i++)
      for (int j=0; j<height; j++) {
        result->r[i][j] = r[j][width-1-i];
        result->g[i][j] = g[j][width-1-i];
        result->b[i][j] = b[j][width-1-i];
      }
    return result;
  }
  else if (deg==180) {
    Imagefloat* result = new Imagefloat (width, height);
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

Imagefloat* Imagefloat::hflip () {
 
  Imagefloat* result = new Imagefloat (width, height);
  for (int i=0; i<height; i++)
    for (int j=0; j<width; j++) {
      result->r[i][j] = r[i][width-1-j];
      result->g[i][j] = g[i][width-1-j];
      result->b[i][j] = b[i][width-1-j];
    }
  return result;

}

Imagefloat* Imagefloat::vflip () {

  Imagefloat* result = new Imagefloat (width, height);
  for (int i=0; i<height; i++)
    for (int j=0; j<width; j++) {
      result->r[i][j] = r[height-1-i][j];
      result->g[i][j] = g[height-1-i][j];
      result->b[i][j] = b[height-1-i][j];
    }
  return result;

}

/*Imagefloat* Imagefloat::resize (int nw, int nh, TypeInterpolation interp) {

    if (interp == TI_Nearest) {
        Imagefloat* res = new Imagefloat (nw, nh);
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
        Imagefloat* res = new Imagefloat (nw, nh);
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
}*/

Image8* 
Imagefloat::to8() const
{
	Image8* img8 = new Image8(width,height);
	for ( int h = 0; h < height; ++h )
	{
		for ( int w = 0; w < width; ++w )
		{
			img8->r(h,w,((int)r[h][w]) >> 8);
			img8->g(h,w,((int)g[h][w]) >> 8);
			img8->b(h,w,((int)b[h][w]) >> 8);
		}
	}
	return img8;
}

Image16* 
Imagefloat::to16() const
{
	Image16* img16 = new Image16(width,height);
	for ( int h = 0; h < height; ++h )
	{
		for ( int w = 0; w < width; ++w )
		{
			img16->r[h][w] = ((int)r[h][w]) ;
			img16->g[h][w] = ((int)g[h][w]) ;
			img16->b[h][w] = ((int)b[h][w]) ;
		}
	}
	return img16;
}


void Imagefloat::calcCroppedHistogram(const ProcParams &params, float scale, LUTu & hist) {
    hist.clear();

    // Set up factors to calc the lightness
    TMatrix wprof = iccStore->workingSpaceMatrix (params.icm.working);

	float facRed   = wprof[1][0];
	float facGreen = wprof[1][1];
	float facBlue  = wprof[1][2];
	

    // calc pixel size
    int x1, x2, y1, y2;
    params.crop.mapToResized(width, height, scale, x1, x2, y1, y2);

    /*
     * WARNING: this parallelization is not 100% safe. Some of its values
     * are underestimated by 1 or 2, but in practical, it doesn't seem to be
     * annoying in any way, so i guess we may leave it as is...
     * (Hombre)
     */
    #pragma omp parallel for
    for (int y=y1; y<y2; y++) {
        int i;
        for (int x=x1; x<x2; x++) {
            i = (int)(facRed * r[y][x] + facGreen * g[y][x] + facBlue * b[y][x]);
            if (i<0) i=0; else if (i>65535) i=65535;
            hist[i]++;
        }
    }
}

// Parallized transformation; create transform with cmsFLAGS_NOCACHE!
void Imagefloat::ExecCMSTransform(cmsHTRANSFORM hTransform) {
    // LittleCMS cannot parallize planar setups
    // so build temporary buffers to allow multi processor execution
        #pragma omp parallel for
    for (int y=0; y<height; y++) {
        float buffer[width*3];
        float *p=buffer, *pR=r[y], *pG=g[y], *pB=b[y];

        for (int x=0; x<width; x++) {
            *(p++) = *(pR++); *(p++) = *(pG++); *(p++) = *(pB++);
    }

        cmsDoTransform (hTransform, buffer, buffer, width);

        p=buffer; pR=r[y]; pG=g[y]; pB=b[y];
        for (int x=0; x<width; x++) {
            *(pR++) = *(p++); *(pG++) = *(p++); *(pB++) = *(p++);
}
    }
}
