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

Image16::Image16 ()
{
}

Image16::Image16 (int w, int h)
{
    allocate (w, h);
}

Image16::~Image16 ()
{
}

void Image16::getScanline (int row, unsigned char* buffer, int bps)
{

    if (data == NULL) {
        return;
    }

    if (bps == 16) {
        int ix = 0;
        unsigned short* sbuffer = (unsigned short*) buffer;

        for (int i = 0; i < width; i++) {
            sbuffer[ix++] = r(row, i);
            sbuffer[ix++] = g(row, i);
            sbuffer[ix++] = b(row, i);
        }
    } else if (bps == 8) {
        int ix = 0;
        int i = 0;
#ifdef __SSSE3__
        // process 48 values using SSSE3. Looks like a lot of code, but it only needs about one instruction per value, whereas scalar version needs about five instructions per value
        vmask reduceWord2Bytev = _mm_set_epi8(0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 15, 13, 11, 9, 7, 5, 3, 1);
        // we need fivev and sixv to reduce the number of registers used for permutation masks from 9 to 6
        vint fivev = _mm_set1_epi8(5);
        vint sixv = _mm_set1_epi8(6);

        for (; i < width - 15; i += 16, ix += 48) {
            // generate initial shuffle masks. Gaps are set to 0xf0 to allow calculating subsequent masks from previous ones
            vint redmaskv = _mm_set_epi8(5, 0xf0, 0xf0, 4, 0xf0, 0xf0, 3, 0xf0, 0xf0, 2, 0xf0, 0xf0, 1, 0xf0, 0xf0, 0);
            vint greenmaskv = _mm_set_epi8(0xf0, 0xf0, 4, 0xf0, 0xf0, 3, 0xf0, 0xf0, 2, 0xf0, 0xf0, 1, 0xf0, 0xf0, 0, 0xf0);
            vint bluemaskv = _mm_set_epi8(0xf0, 4, 0xf0, 0xf0, 3, 0xf0, 0xf0, 2, 0xf0, 0xf0, 1, 0xf0, 0xf0, 0, 0xf0, 0xf0);

            // load first 8 values for each colour
            vint red1v = _mm_loadu_si128((__m128i*)&r(row, i));
            vint green1v = _mm_loadu_si128((__m128i*)&g(row, i));
            vint blue1v = _mm_loadu_si128((__m128i*)&b(row, i));

            // load second 8 values for each colour
            vint red2v = _mm_loadu_si128((__m128i*)&r(row, i + 8));
            vint green2v = _mm_loadu_si128((__m128i*)&g(row, i + 8));
            vint blue2v = _mm_loadu_si128((__m128i*)&b(row, i + 8));

            // shuffle the high bytes of the values to the lower 64 bit of the register
            red1v = _mm_shuffle_epi8(red1v, reduceWord2Bytev);
            green1v = _mm_shuffle_epi8(green1v, reduceWord2Bytev);
            blue1v = _mm_shuffle_epi8(blue1v, reduceWord2Bytev);

            // shuffle the high bytes of the values to the lower 64 bit of the register
            red2v = _mm_shuffle_epi8(red2v, reduceWord2Bytev);
            green2v = _mm_shuffle_epi8(green2v, reduceWord2Bytev);
            blue2v = _mm_shuffle_epi8(blue2v, reduceWord2Bytev);

            // mix first and second 8 values of each colour together
            red1v = (vint)_mm_shuffle_pd((__m128d)red1v, (__m128d)red2v, 0);
            green1v = (vint)_mm_shuffle_pd((__m128d)green1v, (__m128d)green2v, 0);
            blue1v = (vint)_mm_shuffle_pd((__m128d)blue1v, (__m128d)blue2v, 0);

            // now we have the input in registers => let's generate the output

            // first we need r0g0b0r1g1b1r2g2b2r3g3b3r4g4b4r5
            vint destv = _mm_shuffle_epi8(red1v, redmaskv);
            vint greenv = _mm_shuffle_epi8(green1v, greenmaskv);
            destv = _mm_or_si128(destv, greenv);
            vint bluev = _mm_shuffle_epi8(blue1v, bluemaskv);
            destv = _mm_or_si128(destv, bluev);
            _mm_storeu_si128((__m128i*) & (buffer[ix]), destv);

            // then we need g5b5r6g6b6r7g7b7r8g8b8r9g9b9raga
            // we can calculate the shuffle masks from previous ones => needs only 6 instead of 9 registers to handle the 9 different shuffle masks
            vint tempmaskv = _mm_add_epi8(redmaskv, fivev);
            redmaskv = _mm_add_epi8(bluemaskv, sixv);
            bluemaskv = _mm_add_epi8(greenmaskv, fivev);
            greenmaskv = tempmaskv;
            destv = _mm_shuffle_epi8(red1v, redmaskv);
            greenv = _mm_shuffle_epi8(green1v, greenmaskv);
            destv = _mm_or_si128(destv, greenv);
            bluev = _mm_shuffle_epi8(blue1v, bluemaskv);
            destv = _mm_or_si128(destv, bluev);
            _mm_storeu_si128((__m128i*) & (buffer[ix + 16]), destv);

            // and last one is barbgbbbrcgcbcrdgdbdregeberfgfbf
            // we can calculate the shuffle masks from previous ones => needs only 6 instead of 9 registers to handle the 9 different shuffle masks
            tempmaskv = _mm_add_epi8(greenmaskv, fivev);
            greenmaskv = _mm_add_epi8(redmaskv, fivev);
            redmaskv = _mm_add_epi8(bluemaskv, sixv);
            bluemaskv = tempmaskv;
            destv = _mm_shuffle_epi8(red1v, redmaskv);
            greenv = _mm_shuffle_epi8(green1v, greenmaskv);
            destv = _mm_or_si128(destv, greenv);
            bluev = _mm_shuffle_epi8(blue1v, bluemaskv);
            destv = _mm_or_si128(destv, bluev);
            _mm_storeu_si128((__m128i*) & (buffer[ix + 32]), destv);
        }

#endif

        for (; i < width; i++) {
            buffer[ix++] = r(row, i) >> 8;
            buffer[ix++] = g(row, i) >> 8;
            buffer[ix++] = b(row, i) >> 8;
        }
    }
}

/*
 * void Image16::setScanline (int row, unsigned char* buffer, int bps, int minValue[3], int maxValue[3]);
 * has not been implemented yet, because as of now, this method is called for IIOSF_FLOAT sample format only
 */
void Image16::setScanline (int row, unsigned char* buffer, int bps, float *minValue, float *maxValue)
{

    if (data == NULL) {
        return;
    }

    // For optimization purpose, we're assuming that this class never have to provide min/max bound
    assert(!minValue);

    switch (sampleFormat) {
        case (IIOSF_UNSIGNED_CHAR): {
            int ix = 0;

            for (int i = 0; i < width; i++) {
                r(row, i) = (unsigned short)(buffer[ix++]) << 8;
                g(row, i) = (unsigned short)(buffer[ix++]) << 8;
                b(row, i) = (unsigned short)(buffer[ix++]) << 8;
            }

            break;
        }

        case (IIOSF_UNSIGNED_SHORT): {
            unsigned short* sbuffer = (unsigned short*) buffer;
            int ix = 0;

            for (int i = 0; i < width; i++) {
                r(row, i) = sbuffer[ix++];
                g(row, i) = sbuffer[ix++];
                b(row, i) = sbuffer[ix++];
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

Image16* Image16::copy ()
{

    Image16* cp = new Image16 (width, height);
    copyData(cp);
    return cp;
}

void Image16::getStdImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, procparams::ToneCurveParams hrp)
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
    int skip = pp.skip;

    //if ((sx1 + skip*imwidth)>maxx) imwidth -- ; // we have a boundary condition that can cause errors

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

        // Iterating all the rows of the destination image
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
#undef GCLIP
}

Image8*
Image16::to8()
{
    Image8* img8 = new Image8(width, height);

    for ( int h = 0; h < height; ++h ) {
        for ( int w = 0; w < width; ++w ) {
            img8->r(h, w) = (unsigned char)( r(h, w) >> 8);
            img8->g(h, w) = (unsigned char)( g(h, w) >> 8);
            img8->b(h, w) = (unsigned char)( b(h, w) >> 8);
        }
    }

    return img8;
}

Imagefloat*
Image16::tofloat()
{
    Imagefloat* imgfloat = new Imagefloat(width, height);

    for ( int h = 0; h < height; ++h ) {
        for ( int w = 0; w < width; ++w ) {
            imgfloat->r(h, w) = (float)r(h, w);
            imgfloat->g(h, w) = (float)g(h, w);
            imgfloat->b(h, w) = (float)b(h, w);
        }
    }

    return imgfloat;
}

// Parallized transformation; create transform with cmsFLAGS_NOCACHE!
void Image16::ExecCMSTransform(cmsHTRANSFORM hTransform)
{
    //cmsDoTransform(hTransform, data, data, planestride);

    // LittleCMS cannot parallelize planar setups -- Hombre: LCMS2.4 can! But it we use this new feature, memory allocation have to be modified too
    // so build temporary buffers to allow multi processor execution
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        AlignedBuffer<unsigned short> buffer(width * 3);

#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif

        for (int y = 0; y < height; y++)
        {
            unsigned short *p = buffer.data, *pR = r(y), *pG = g(y), *pB = b(y);

            for (int x = 0; x < width; x++) {
                *(p++) = *(pR++);
                *(p++) = *(pG++);
                *(p++) = *(pB++);
            }

            cmsDoTransform (hTransform, buffer.data, buffer.data, width);

            p = buffer.data;
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
