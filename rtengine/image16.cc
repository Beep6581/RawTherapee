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
#include "image16.h"
#include "imagefloat.h"
#include "image8.h"
#include <cstring>
#include <cstdio>
#include "rtengine.h"

using namespace rtengine;

Image16::Image16 () {
}

Image16::Image16 (int w, int h) {
    allocate (w, h);
}

Image16::~Image16 () {
}

void Image16::getScanline (int row, unsigned char* buffer, int bps) {

    if (data==NULL)
        return;

    if (bps==16) {
        int ix = 0;
        unsigned short* sbuffer = (unsigned short*) buffer;
        for (int i=0; i<width; i++) {
            sbuffer[ix++] = r(row,i);
            sbuffer[ix++] = g(row,i);
            sbuffer[ix++] = b(row,i);
        }
    }
    else if (bps==8) {
        int ix = 0;
        for (int i=0; i<width; i++) {
            buffer[ix++] = r(row,i) >> 8;
            buffer[ix++] = g(row,i) >> 8;
            buffer[ix++] = b(row,i) >> 8;
        }
    }
}

/*
 * void Image16::setScanline (int row, unsigned char* buffer, int bps, int minValue[3], int maxValue[3]);
 * has not been implemented yet, because as of now, this method is called for IIOSF_FLOAT sample format only
 */
void Image16::setScanline (int row, unsigned char* buffer, int bps, float *minValue, float *maxValue) {

    if (data==NULL)
        return;

     // For optimization purpose, we're assuming that this class never have to provide min/max bound
    assert(!minValue);

    switch (sampleFormat) {
    case (IIOSF_UNSIGNED_CHAR):
    {
        int ix = 0;
        for (int i=0; i<width; i++) {
            r(row,i) = (unsigned short)(buffer[ix++]) << 8;
            g(row,i) = (unsigned short)(buffer[ix++]) << 8;
            b(row,i) = (unsigned short)(buffer[ix++]) << 8;
        }
        break;
    }
    case (IIOSF_UNSIGNED_SHORT):
    {
        unsigned short* sbuffer = (unsigned short*) buffer;
        int ix = 0;
        for (int i=0; i<width; i++) {
            r(row,i) = sbuffer[ix++];
            g(row,i) = sbuffer[ix++];
            b(row,i) = sbuffer[ix++];
        }
        break;
    }
    default:
        // Other type are ignored, but could be implemented if necessary
        break;
    }
    /*
     * Not used for now
     *
     */
}

Image16* Image16::copy () {

  Image16* cp = new Image16 (width, height);
  copyData(cp);
  return cp;
}

void Image16::getStdImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, procparams::HRecParams hrp)
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
                    rtot += Color::igamma_srgb(r(i+m, jx+n));
                    gtot += Color::igamma_srgb(g(i+m, jx+n));
                    btot += Color::igamma_srgb(b(i+m, jx+n));
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

Image8* 
Image16::to8()
{
    Image8* img8 = new Image8(width,height);
    for ( int h = 0; h < height; ++h )
    {
        for ( int w = 0; w < width; ++w )
        {
            img8->r(h, w) = (unsigned char)( r(h,w) >> 8);
            img8->g(h, w) = (unsigned char)( g(h,w) >> 8);
            img8->b(h, w) = (unsigned char)( b(h,w) >> 8);
        }
    }
    return img8;
}

Imagefloat* 
Image16::tofloat()
{
    Imagefloat* imgfloat = new Imagefloat(width,height);
    for ( int h = 0; h < height; ++h )
    {
        for ( int w = 0; w < width; ++w )
        {
            imgfloat->r(h,w) = (float)r(h,w);
            imgfloat->g(h,w) = (float)g(h,w);
            imgfloat->b(h,w) = (float)b(h,w);
        }
    }
    return imgfloat;
}

// Parallized transformation; create transform with cmsFLAGS_NOCACHE!
void Image16::ExecCMSTransform(cmsHTRANSFORM hTransform) {
    //cmsDoTransform(hTransform, data, data, planestride);

    // LittleCMS cannot parallelize planar setups -- Hombre: LCMS2.4 can! But it we use this new feature, memory allocation have to be modified too
    // so build temporary buffers to allow multi processor execution
#ifdef _OPENMP
#pragma omp parallel
#endif
{
        AlignedBuffer<unsigned short> buffer(width*3);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int y=0; y<height; y++) {
        unsigned short *p=buffer.data, *pR=r(y), *pG=g(y), *pB=b(y);

        for (int x=0; x<width; x++) {
            *(p++) = *(pR++); *(p++) = *(pG++); *(p++) = *(pB++);
        }

        cmsDoTransform (hTransform, buffer.data, buffer.data, width);

        p=buffer.data; pR=r(y); pG=g(y); pB=b(y);
        for (int x=0; x<width; x++) {
            *(pR++) = *(p++); *(pG++) = *(p++); *(pB++) = *(p++);
        }
} // End of parallelization
    }
}
