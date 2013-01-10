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
#include <tiffio.h>
#include "imagefloat.h"
#include "image16.h"
#include "image8.h"
#include <cstring>
#include "rtengine.h"
#include "mytime.h"
#include "iccstore.h"
#include "alignedbuffer.h"
#include "rt_math.h"
#include "color.h"

using namespace rtengine;

Imagefloat::Imagefloat () {
}

Imagefloat::Imagefloat (int w, int h) {
    allocate (w, h);
}

Imagefloat::~Imagefloat () {
}

// Call this method to handle floating points input values of different size
void Imagefloat::setScanline (int row, unsigned char* buffer, int bps, float *minValue, float *maxValue) {

    if (data==NULL)
        return;

    switch (sampleFormat) {
    case (IIOSF_FLOAT):
    {
        int ix = 0;
        float* sbuffer = (float*) buffer;
        for (int i=0; i<width; i++) {
            r(row,i) = sbuffer[ix]; if (minValue) { if (sbuffer[ix]<minValue[0]) minValue[0] = sbuffer[ix]; else if (sbuffer[ix]>maxValue[0]) maxValue[0] = sbuffer[ix]; ++ix; }
            g(row,i) = sbuffer[ix]; if (minValue) { if (sbuffer[ix]<minValue[1]) minValue[1] = sbuffer[ix]; else if (sbuffer[ix]>maxValue[1]) maxValue[1] = sbuffer[ix]; ++ix; }
            b(row,i) = sbuffer[ix]; if (minValue) { if (sbuffer[ix]<minValue[2]) minValue[2] = sbuffer[ix]; else if (sbuffer[ix]>maxValue[2]) maxValue[2] = sbuffer[ix]; ++ix; }
        }
        break;
    }
    case (IIOSF_LOGLUV24):
    case (IIOSF_LOGLUV32):
    {
        int ix = 0;
        float* sbuffer = (float*) buffer;
        float xyzvalues[3], rgbvalues[3];
        for (int i=0; i<width; i++) {
            xyzvalues[0] = sbuffer[ix++];
            xyzvalues[1] = sbuffer[ix++];
            xyzvalues[2] = sbuffer[ix++];
            // TODO: we may have to handle other color space than sRGB!
            Color::xyz2srgb(xyzvalues[0], xyzvalues[1], xyzvalues[2], rgbvalues[0], rgbvalues[1], rgbvalues[2]);
            r(row,i) = rgbvalues[0]; if (minValue) { if (rgbvalues[0]<minValue[0]) minValue[0] = rgbvalues[0]; else if (rgbvalues[0]>maxValue[0]) maxValue[0] = rgbvalues[0]; }
            g(row,i) = rgbvalues[1]; if (minValue) { if (rgbvalues[1]<minValue[1]) minValue[1] = rgbvalues[1]; else if (rgbvalues[1]>maxValue[1]) maxValue[1] = rgbvalues[1]; }
            b(row,i) = rgbvalues[2]; if (minValue) { if (rgbvalues[2]<minValue[2]) minValue[2] = rgbvalues[2]; else if (rgbvalues[2]>maxValue[2]) maxValue[2] = rgbvalues[2]; }
        }
        break;
    }
    default:
        // Other type are ignored, but could be implemented if necessary
        break;
    }
}

void Imagefloat::getScanline (int row, unsigned char* buffer, int bps) {

    if (data==NULL)
        return;

    if (bps==32) {
        int ix = 0;
        float* sbuffer = (float*) buffer;
        for (int i=0; i<width; i++) {
            sbuffer[ix++] = r(row,i);
            sbuffer[ix++] = g(row,i);
            sbuffer[ix++] = b(row,i);
        }
    }
}

Imagefloat* Imagefloat::copy () {

  Imagefloat* cp = new Imagefloat (width, height);
  copyData(cp);
  return cp;
}

// This is called by the StdImageSource class. We assume that fp images from StdImageSource don't have to deal with gamma
void Imagefloat::getStdImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, procparams::HRecParams hrp)
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

    // improve speed by integrating the area division into the multipliers
    // switched to using ints for the red/green/blue channel buffer.
    // Incidentally this improves accuracy too.
    float area=skip*skip;
    rm/=area;
    gm/=area;
    bm/=area;

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
                    rtot += r(i+m, jx+n);
                    gtot += g(i+m, jx+n);
                    btot += b(i+m, jx+n);
                }
            lineR[j] = CLIP(rm*rtot);
            lineG[j] = CLIP(gm*gtot);
            lineB[j] = CLIP(bm*btot);
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
}

Image8* 
Imagefloat::to8()
{
    Image8* img8 = new Image8(width,height);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for ( int h=0; h < height; ++h )
    {
        for ( int w=0; w < width; ++w )
        {
            img8->r(h, w) = (unsigned char)( (unsigned short)(r(h,w)) >> 8);
            img8->g(h, w) = (unsigned char)( (unsigned short)(g(h,w)) >> 8);
            img8->b(h, w) = (unsigned char)( (unsigned short)(b(h,w)) >> 8);
        }
    }
    return img8;
}

Image16* 
Imagefloat::to16()
{
    Image16* img16 = new Image16(width,height);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for ( int h=0; h < height; ++h )
    {
        for ( int w=0; w < width; ++w )
        {
            img16->r( h,w) = (unsigned short)(r(h,w));
            img16->g( h,w) = (unsigned short)(g(h,w));
            img16->b( h,w) = (unsigned short)(b(h,w));
        }
    }
    return img16;
}

void Imagefloat::normalizeFloat(float srcMinVal, float srcMaxVal) {

    float scale = MAXVALD / (srcMaxVal-srcMinVal);
    int w = width;
    int h = height;

#ifdef _OPENMP
#pragma omp parallel for firstprivate(w, h, srcMinVal, scale) schedule(dynamic, 5)
#endif
    for (int y=0; y<h; y++) {
        for (int x=0; x<w; x++) {
             r(y,x) = (r(y,x)-srcMinVal)*scale;
             g(y,x) = (g(y,x)-srcMinVal)*scale;
             b(y,x) = (b(y,x)-srcMinVal)*scale;
        }
    }
}

// convert values's range to [0;1] ; this method assumes that the input values's range is [0;65535]
void Imagefloat::normalizeFloatTo1() {

    int w = width;
    int h = height;

#ifdef _OPENMP
#pragma omp parallel for firstprivate(w, h) schedule(dynamic, 5)
#endif
    for (int y=0; y<h; y++) {
        for (int x=0; x<w; x++) {
             r(y,x) /= 65535.f;
             g(y,x) /= 65535.f;
             b(y,x) /= 65535.f;
        }
    }
}

// convert values's range to [0;65535 ; this method assumes that the input values's range is [0;1]
void Imagefloat::normalizeFloatTo65535() {

    int w = width;
    int h = height;

#ifdef _OPENMP
#pragma omp parallel for firstprivate(w, h) schedule(dynamic, 5)
#endif
    for (int y=0; y<h; y++) {
        for (int x=0; x<w; x++) {
             r(y,x) *= 65535.f;
             g(y,x) *= 65535.f;
             b(y,x) *= 65535.f;
        }
    }
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int y=y1; y<y2; y++) {
        int i;
        for (int x=x1; x<x2; x++) {
            i = (int)(facRed * r(y,x) + facGreen * g(y,x) + facBlue * b(y,x));
            if (i<0) i=0; else if (i>65535) i=65535;
#ifdef _OPENMP
// Access to hist[] must be atomic. In this case, we may need to see if this parallelization is worth it
#pragma omp atomic
#endif
            hist[i]++;
        }
    }
}

// Parallelized transformation; create transform with cmsFLAGS_NOCACHE!
void Imagefloat::ExecCMSTransform(cmsHTRANSFORM hTransform) {

    // LittleCMS cannot parallelize planar setups -- Hombre: LCMS2.4 can! But it we use this new feature, memory allocation
	// have to be modified too to build temporary buffers that allow multi processor execution
#ifdef _OPENMP
#pragma omp parallel
#endif
{
        AlignedBuffer<float> pBuf(width*3);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int y=0; y<height; y++) {
        float *p=pBuf.data, *pR=r(y), *pG=g(y), *pB=b(y);

        for (int x=0; x<width; x++) {
            *(p++) = *(pR++); *(p++) = *(pG++); *(p++) = *(pB++);
        }

        cmsDoTransform (hTransform, pBuf.data, pBuf.data, width);

        p=pBuf.data; pR=r(y); pG=g(y); pB=b(y);
        for (int x=0; x<width; x++) {
            *(pR++) = *(p++); *(pG++) = *(p++); *(pB++) = *(p++);
        }
} // End of parallelization
    }
}
