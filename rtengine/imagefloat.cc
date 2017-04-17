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

Imagefloat::Imagefloat ()
{
}

Imagefloat::Imagefloat (int w, int h)
{
    allocate (w, h);
}

Imagefloat::~Imagefloat ()
{
}

// Call this method to handle floating points input values of different size
void Imagefloat::setScanline (int row, unsigned char* buffer, int bps, float *minValue, float *maxValue)
{

    if (data == nullptr) {
        return;
    }

    switch (sampleFormat) {
    case (IIOSF_FLOAT): {
        int ix = 0;
        float* sbuffer = (float*) buffer;

        for (int i = 0; i < width; i++) {
            r(row, i) = sbuffer[ix];

            if (minValue) {
                if (sbuffer[ix] < minValue[0]) {
                    minValue[0] = sbuffer[ix];
                } else if (sbuffer[ix] > maxValue[0]) {
                    maxValue[0] = sbuffer[ix];
                } ++ix;
            }

            g(row, i) = sbuffer[ix];

            if (minValue) {
                if (sbuffer[ix] < minValue[1]) {
                    minValue[1] = sbuffer[ix];
                } else if (sbuffer[ix] > maxValue[1]) {
                    maxValue[1] = sbuffer[ix];
                } ++ix;
            }

            b(row, i) = sbuffer[ix];

            if (minValue) {
                if (sbuffer[ix] < minValue[2]) {
                    minValue[2] = sbuffer[ix];
                } else if (sbuffer[ix] > maxValue[2]) {
                    maxValue[2] = sbuffer[ix];
                } ++ix;
            }
        }

        break;
    }

    case (IIOSF_LOGLUV24):
    case (IIOSF_LOGLUV32): {
        int ix = 0;
        float* sbuffer = (float*) buffer;
        float xyzvalues[3], rgbvalues[3];

        for (int i = 0; i < width; i++) {
            xyzvalues[0] = sbuffer[ix++];
            xyzvalues[1] = sbuffer[ix++];
            xyzvalues[2] = sbuffer[ix++];
            // TODO: we may have to handle other color space than sRGB!
            Color::xyz2srgb(xyzvalues[0], xyzvalues[1], xyzvalues[2], rgbvalues[0], rgbvalues[1], rgbvalues[2]);
            r(row, i) = rgbvalues[0];

            if (minValue) {
                if (rgbvalues[0] < minValue[0]) {
                    minValue[0] = rgbvalues[0];
                } else if (rgbvalues[0] > maxValue[0]) {
                    maxValue[0] = rgbvalues[0];
                }
            }

            g(row, i) = rgbvalues[1];

            if (minValue) {
                if (rgbvalues[1] < minValue[1]) {
                    minValue[1] = rgbvalues[1];
                } else if (rgbvalues[1] > maxValue[1]) {
                    maxValue[1] = rgbvalues[1];
                }
            }

            b(row, i) = rgbvalues[2];

            if (minValue) {
                if (rgbvalues[2] < minValue[2]) {
                    minValue[2] = rgbvalues[2];
                } else if (rgbvalues[2] > maxValue[2]) {
                    maxValue[2] = rgbvalues[2];
                }
            }
        }

        break;
    }

    default:
        // Other type are ignored, but could be implemented if necessary
        break;
    }
}

void Imagefloat::getScanline (int row, unsigned char* buffer, int bps)
{

    if (data == nullptr) {
        return;
    }

    if (bps == 32) {
        int ix = 0;
        float* sbuffer = (float*) buffer;

        for (int i = 0; i < width; i++) {
            sbuffer[ix++] = r(row, i);
            sbuffer[ix++] = g(row, i);
            sbuffer[ix++] = b(row, i);
        }
    }
}

Imagefloat* Imagefloat::copy ()
{

    Imagefloat* cp = new Imagefloat (width, height);
    copyData(cp);
    return cp;
}

// This is called by the StdImageSource class. We assume that fp images from StdImageSource don't have to deal with gamma
void Imagefloat::getStdImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, procparams::ToneCurveParams hrp)
{

    // compute channel multipliers
    double drm, dgm, dbm;
    ctemp.getMultipliers (drm, dgm, dbm);
    float rm = drm, gm = dgm, bm = dbm;

    rm = 1.0 / rm;
    gm = 1.0 / gm;
    bm = 1.0 / bm;
    float mul_lum = 0.299 * rm + 0.587 * gm + 0.114 * bm;
    rm /= mul_lum;
    gm /= mul_lum;
    bm /= mul_lum;

    int sx1, sy1, sx2, sy2;

    transform (pp, tran, sx1, sy1, sx2, sy2);

    int imwidth = image->width; // Destination image
    int imheight = image->height; // Destination image

    if (((tran & TR_ROT) == TR_R90) || ((tran & TR_ROT) == TR_R270)) {
        int swap = imwidth;
        imwidth = imheight;
        imheight = swap;
    }

    int maxx = width; // Source image
    int maxy = height; // Source image
    int mtran = tran & TR_ROT;
    int skip = pp.getSkip();

    // improve speed by integrating the area division into the multipliers
    // switched to using ints for the red/green/blue channel buffer.
    // Incidentally this improves accuracy too.
    float area = skip * skip;
    float rm2 = rm;
    float gm2 = gm;
    float bm2 = bm;
    rm /= area;
    gm /= area;
    bm /= area;

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

        for (int iy = 0; iy < imheight; iy++) {
            if (skip == 1) {
                // special case (speedup for 1:1 scale)
                // i: source image, first line of the current destination row
                int src_y = sy1 + iy;

                // overflow security check, not sure that it's necessary
                if (src_y >= maxy) {
                    continue;
                }

                for (int dst_x = 0, src_x = sx1; dst_x < imwidth; dst_x++, src_x++) {
                    // overflow security check, not sure that it's necessary
                    if (src_x >= maxx) {
                        continue;
                    }

                    lineR[dst_x] = CLIP(rm2 * r(src_y, src_x));
                    lineG[dst_x] = CLIP(gm2 * g(src_y, src_x));
                    lineB[dst_x] = CLIP(bm2 * b(src_y, src_x));
                }
            } else {
                // source image, first line of the current destination row
                int src_y = sy1 + skip * iy;

                if (src_y >= maxy) {
                    continue;
                }

                for (int dst_x = 0, src_x = sx1; dst_x < imwidth; dst_x++, src_x += skip) {
                    if (src_x >= maxx) {
                        continue;
                    }

                    int src_sub_width = MIN(maxx - src_x, skip);
                    int src_sub_height = MIN(maxy - src_y, skip);

                    float rtot, gtot, btot; // RGB accumulators
                    rtot = gtot = btot = 0.;

                    for (int src_sub_y = 0; src_sub_y < src_sub_height; src_sub_y++)
                        for (int src_sub_x = 0; src_sub_x < src_sub_width; src_sub_x++) {
                            rtot += r(src_y + src_sub_y, src_x + src_sub_x);
                            gtot += g(src_y + src_sub_y, src_x + src_sub_x);
                            btot += b(src_y + src_sub_y, src_x + src_sub_x);
                        }

                    // convert back to gamma and clip
                    if (src_sub_width == skip && src_sub_height == skip) {
                        // Common case where the sub-region is complete
                        lineR[dst_x] = CLIP(rm * rtot);
                        lineG[dst_x] = CLIP(gm * gtot);
                        lineB[dst_x] = CLIP(bm * btot);
                    } else {
                        // computing a special factor for this incomplete sub-region
                        float area = src_sub_width * src_sub_height;
                        lineR[dst_x] = CLIP(rm2 * rtot / area);
                        lineG[dst_x] = CLIP(gm2 * gtot / area);
                        lineB[dst_x] = CLIP(bm2 * btot / area);
                    }
                }
            }

            if      (mtran == TR_NONE)
                for (int dst_x = 0, src_x = sx1; dst_x < imwidth; dst_x++, src_x += skip) {
                    image->r(iy, dst_x) = lineR[dst_x];
                    image->g(iy, dst_x) = lineG[dst_x];
                    image->b(iy, dst_x) = lineB[dst_x];
                }
            else if (mtran == TR_R180)
                for (int dst_x = 0; dst_x < imwidth; dst_x++) {
                    image->r(imheight - 1 - iy, imwidth - 1 - dst_x) = lineR[dst_x];
                    image->g(imheight - 1 - iy, imwidth - 1 - dst_x) = lineG[dst_x];
                    image->b(imheight - 1 - iy, imwidth - 1 - dst_x) = lineB[dst_x];
                }
            else if (mtran == TR_R90)
                for (int dst_x = 0, src_x = sx1; dst_x < imwidth; dst_x++, src_x += skip) {
                    image->r(dst_x, imheight - 1 - iy) = lineR[dst_x];
                    image->g(dst_x, imheight - 1 - iy) = lineG[dst_x];
                    image->b(dst_x, imheight - 1 - iy) = lineB[dst_x];
                }
            else if (mtran == TR_R270)
                for (int dst_x = 0, src_x = sx1; dst_x < imwidth; dst_x++, src_x += skip) {
                    image->r(imwidth - 1 - dst_x, iy) = lineR[dst_x];
                    image->g(imwidth - 1 - dst_x, iy) = lineG[dst_x];
                    image->b(imwidth - 1 - dst_x, iy) = lineB[dst_x];
                }
        }

#ifdef _OPENMP
    }
#endif
}

Image8*
Imagefloat::to8()
{
    Image8* img8 = new Image8(width, height);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif

    for (int h = 0; h < height; ++h) {
        for (int w = 0; w < width; ++w) {
            img8->r(h, w) = uint16ToUint8Rounded(r(h, w));
            img8->g(h, w) = uint16ToUint8Rounded(g(h, w));
            img8->b(h, w) = uint16ToUint8Rounded(b(h, w));
        }
    }

    return img8;
}

Image16*
Imagefloat::to16()
{
    Image16* img16 = new Image16(width, height);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif

    for (int h = 0; h < height; ++h) {
        for (int w = 0; w < width; ++w) {
            img16->r(h, w) = r(h, w);
            img16->g(h, w) = g(h, w);
            img16->b(h, w) = b(h, w);
        }
    }

    return img16;
}

void Imagefloat::normalizeFloat(float srcMinVal, float srcMaxVal)
{

    float scale = MAXVALD / (srcMaxVal - srcMinVal);
    int w = width;
    int h = height;

#ifdef _OPENMP
    #pragma omp parallel for firstprivate(w, h, srcMinVal, scale) schedule(dynamic, 5)
#endif

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            r(y, x) = (r(y, x) - srcMinVal) * scale;
            g(y, x) = (g(y, x) - srcMinVal) * scale;
            b(y, x) = (b(y, x) - srcMinVal) * scale;
        }
    }
}

// convert values's range to [0;1] ; this method assumes that the input values's range is [0;65535]
void Imagefloat::normalizeFloatTo1()
{

    int w = width;
    int h = height;

#ifdef _OPENMP
    #pragma omp parallel for firstprivate(w, h) schedule(dynamic, 5)
#endif

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            r(y, x) /= 65535.f;
            g(y, x) /= 65535.f;
            b(y, x) /= 65535.f;
        }
    }
}

// convert values's range to [0;65535 ; this method assumes that the input values's range is [0;1]
void Imagefloat::normalizeFloatTo65535()
{

    int w = width;
    int h = height;

#ifdef _OPENMP
    #pragma omp parallel for firstprivate(w, h) schedule(dynamic, 5)
#endif

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            r(y, x) *= 65535.f;
            g(y, x) *= 65535.f;
            b(y, x) *= 65535.f;
        }
    }
}

void Imagefloat::calcCroppedHistogram(const ProcParams &params, float scale, LUTu & hist)
{

    hist.clear();

    // Set up factors to calc the lightness
    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params.icm.working);

    float facRed   = wprof[1][0];
    float facGreen = wprof[1][1];
    float facBlue  = wprof[1][2];


    // calc pixel size
    int x1, x2, y1, y2;
    params.crop.mapToResized(width, height, scale, x1, x2, y1, y2);

    #pragma omp parallel
    {
        LUTu histThr(65536);
        histThr.clear();
        #pragma omp for nowait

        for (int y = y1; y < y2; y++) {
            for (int x = x1; x < x2; x++) {
                int i = (int)(facRed * r(y, x) + facGreen * g(y, x) + facBlue * b(y, x));

                if (i < 0) {
                    i = 0;
                } else if (i > 65535) {
                    i = 65535;
                }

                histThr[i]++;
            }
        }

        #pragma omp critical
        {
            for(int i = 0; i <= 0xffff; i++) {
                hist[i] += histThr[i];
            }
        }
    }

}

// Parallelized transformation; create transform with cmsFLAGS_NOCACHE!
void Imagefloat::ExecCMSTransform(cmsHTRANSFORM hTransform)
{

    // LittleCMS cannot parallelize planar setups -- Hombre: LCMS2.4 can! But it we use this new feature, memory allocation
    // have to be modified too to build temporary buffers that allow multi processor execution
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        AlignedBuffer<float> pBuf(width * 3);

#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif

        for (int y = 0; y < height; y++)
        {
            float *p = pBuf.data, *pR = r(y), *pG = g(y), *pB = b(y);

            for (int x = 0; x < width; x++) {
                *(p++) = *(pR++);
                *(p++) = *(pG++);
                *(p++) = *(pB++);
            }

            cmsDoTransform (hTransform, pBuf.data, pBuf.data, width);

            p = pBuf.data;
            pR = r(y);
            pG = g(y);
            pB = b(y);

            for (int x = 0; x < width; x++) {
                *(pR++) = *(p++);
                *(pG++) = *(p++);
                *(pB++) = *(p++);
            }
        } // End of parallelization
    }
}
