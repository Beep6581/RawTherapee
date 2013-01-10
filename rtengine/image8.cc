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
#include <cstring>
#include <cstdio>
#include "image8.h"
#include "rtengine.h"

using namespace rtengine;


Image8::Image8 () {
}

Image8::Image8 (int w, int h) {
    allocate (w, h);
}

Image8::~Image8 () {
}

void Image8::getScanline (int row, unsigned char* buffer, int bps) {

    if (data==NULL)
        return;

    if (bps==8)
        memcpy (buffer, data + row*width*3, width*3);
    else if (bps==16) {
        unsigned short* sbuffer = (unsigned short*) buffer;
        for (int i=0, ix = row*width*3; i<width*3; i++, ix++)
            sbuffer[i] = (unsigned short)(data[ix]) << 8;
    }
}

void Image8::setScanline (int row, unsigned char* buffer, int bps, float *minValue, float *maxValue) {

    if (data==NULL)
        return;

    // For optimization purpose, we're assuming that this class never have to provide min/max bound
   assert(!minValue);

   switch (sampleFormat) {
   case (IIOSF_UNSIGNED_CHAR):
       memcpy (data + row*width*3, buffer, width*3);
       break;
   case (IIOSF_UNSIGNED_SHORT):
   {
        unsigned short* sbuffer = (unsigned short*) buffer;
        for (int i=0, ix = row*width*3; i<width*3; i++, ix++)
            data[ix] = sbuffer[i] >> 8;
        break;
   }
   default:
       // Other type are ignored, but could be implemented if necessary
       break;
   }
}

Image8* Image8::copy () {

  Image8* cp = new Image8 (width, height);
  copyData(cp);
  return cp;
}

void Image8::getStdImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, procparams::HRecParams hrp)
{
    // compute channel multipliers
    double drm, dgm, dbm;
    ctemp.getMultipliers (drm, dgm, dbm);
    float rm=drm,gm=dgm,bm=dbm;

    rm = 1.0 / rm;
    gm = 1.0 / gm;
    bm = 1.0 / bm;
    float mul_lum = 0.299*rm + 0.587*gm + 0.114*bm;
    rm /= mul_lum;
    gm /= mul_lum;
    bm /= mul_lum;

    int sx1, sy1, sx2, sy2;

    transform (pp, tran, sx1, sy1, sx2, sy2);

    int imwidth=image->width,imheight=image->height;
    if (((tran & TR_ROT) == TR_R90)||((tran & TR_ROT) == TR_R270)) {
        int swap = imwidth;
        imwidth=imheight;
        imheight=swap;
    }
    int istart = sy1;
    int maxx=width,maxy=height;
    int mtran = tran & TR_ROT;
    int skip = pp.skip;

    //if ((sx1 + skip*imwidth)>maxx) imwidth -- ; // we have a boundary condition that can cause errors

    // improve speed by integrating the area division into the multipliers
    // switched to using ints for the red/green/blue channel buffer.
    // Incidentally this improves accuracy too.
    float area=skip*skip;
    rm/=area;
    gm/=area;
    bm/=area;

    #define GCLIP( x ) Color::gamma_srgb(CLIP(x))

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    AlignedBuffer<float> abR(imwidth);
    AlignedBuffer<float> abG(imwidth);
    AlignedBuffer<float> abB(imwidth);
    float *lineR  = abR.data;
    float *lineG  = abG.data;
    float *lineB =  abB.data;

#ifdef _OPENMP
#pragma omp for
#endif
    for (int ix=0;ix<imheight;ix++) {
        int i=istart+skip*ix;if (i>=maxy-skip) i=maxy-skip-1; // avoid trouble
        for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {
            if (jx>=maxx-skip) jx=maxx-skip-1; // avoid trouble

            float rtot,gtot,btot;
            rtot=gtot=btot=0;

            for (int m=0; m<skip; m++)
                for (int n=0; n<skip; n++) {
                    unsigned short r_, g_, b_;
                    convertTo(r(i+m, jx+n), r_);
                    convertTo(g(i+m, jx+n), g_);
                    convertTo(b(i+m, jx+n), b_);
                    rtot += Color::igamma_srgb(r_);
                    gtot += Color::igamma_srgb(g_);
                    btot += Color::igamma_srgb(b_);
                }
            // convert back to gamma and clip
            lineR[j] = GCLIP(rm*rtot);
            lineG[j] = GCLIP(gm*gtot);
            lineB[j] = GCLIP(bm*btot);
        }

        if      (mtran == TR_NONE)
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {
                image->r(ix, j) = lineR[j];
                image->g(ix, j) = lineG[j];
                image->b(ix, j) = lineB[j];
            }
        else if (mtran == TR_R180)
            for (int j=0; j<imwidth; j++) {
                image->r(imheight-1-ix, imwidth-1-j) = lineR[j];
                image->g(imheight-1-ix, imwidth-1-j) = lineG[j];
                image->b(imheight-1-ix, imwidth-1-j) = lineB[j];
            }
        else if (mtran == TR_R90)
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {
                image->r(j, imheight-1-ix) = lineR[j];
                image->g(j, imheight-1-ix) = lineG[j];
                image->b(j, imheight-1-ix) = lineB[j];
            }
        else if (mtran == TR_R270)
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {
                image->r(imwidth-1-j, ix) = lineR[j];
                image->g(imwidth-1-j, ix) = lineG[j];
                image->b(imwidth-1-j, ix) = lineB[j];
            }
    }
#ifdef _OPENMP
}
#endif
    #undef GCLIP
}
