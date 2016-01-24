////////////////////////////////////////////////////////////////
//
//          AMaZE demosaic algorithm
// (Aliasing Minimization and Zipper Elimination)
//
//  copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//  optimized for speed by Ingo Weyrich
//
// incorporating ideas of Luis Sanz Rodrigues and Paul Lee
//
// code dated: May 27, 2010
//
//  amaze_interpolate_RT.cc is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include "rtengine.h"
#include "rawimagesource.h"
#include "rt_math.h"
#include "../rtgui/multilangmgr.h"
#include "sleef.c"
#include "opthelper.h"
#define BENCHMARK
#include "StopWatch.h"

namespace rtengine
{

SSEFUNCTION void RawImageSource::amaze_demosaic_RT(int winx, int winy, int winw, int winh)
{
    BENCHFUN

    volatile double progress = 0.0;

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::amaze]));
        plistener->setProgress (0.0);
    }

    const int width = winw, height = winh;
    const float clip_pt = 1.0 / initialGain;
    const float clip_pt8 = 0.8 / initialGain;


#define TS 160   // Tile size; the image is processed in square tiles to lower memory requirements and facilitate multi-threading
#define TSH 80   // half of Tile size

    //offset of R pixel within a Bayer quartet
    int ex, ey;

    //determine GRBG coset; (ey,ex) is the offset of the R subarray
    if (FC(0, 0) == 1) { //first pixel is G
        if (FC(0, 1) == 0) {
            ey = 0;
            ex = 1;
        } else {
            ey = 1;
            ex = 0;
        }
    } else {//first pixel is R or B
        if (FC(0, 0) == 0) {
            ey = 0;
            ex = 0;
        } else {
            ey = 1;
            ex = 1;
        }
    }

    //shifts of pointer value to access pixels in vertical and diagonal directions
    static const int v1 = TS, v2 = 2 * TS, v3 = 3 * TS, p1 = -TS + 1, p2 = -2 * TS + 2, p3 = -3 * TS + 3, m1 = TS + 1, m2 = 2 * TS + 2, m3 = 3 * TS + 3;

    //tolerance to avoid dividing by zero
    static const float eps = 1e-5, epssq = 1e-10;       //tolerance to avoid dividing by zero

    //adaptive ratios threshold
    static const float arthresh = 0.75;

    //gaussian on 5x5 quincunx, sigma=1.2
    static const float gaussodd[4] = {0.14659727707323927f, 0.103592713382435f, 0.0732036125103057f, 0.0365543548389495f};
    //nyquist texture test threshold
    static const float nyqthresh = 0.5;
    //gaussian on 5x5, sigma=1.2, multiplied with nyqthresh to save some time later in loop
    // Is this really sigma=1.2????, seems more like sigma = 1.672
    static const float gaussgrad[6] = {nyqthresh * 0.07384411893421103f, nyqthresh * 0.06207511968171489f, nyqthresh * 0.0521818194747806f,
                                       nyqthresh * 0.03687419286733595f, nyqthresh * 0.03099732204057846f, nyqthresh * 0.018413194161458882f
                                      };
    //gaussian on 5x5 alt quincunx, sigma=1.5
    static const float gausseven[2] = {0.13719494435797422f, 0.05640252782101291f};
    //guassian on quincunx grid
    static const float gquinc[4] = {0.169917f, 0.108947f, 0.069855f, 0.0287182f};

    typedef struct {
        float h;
        float v;
    } s_hv;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int progresscounter = 0;

#define CLF 1
        // assign working space
        char *buffer = (char *) calloc(13 * sizeof(float) * TS * TS + sizeof(float) * TS * TSH + sizeof(char) * TS * TSH + 18 * CLF * 64 + 63, 1);
        // aligned to 64 byte boundary
        char *data = (char*)( ( uintptr_t(buffer) + uintptr_t(63)) / 64 * 64);

        // green values
        float *rgbgreen         = (float (*))         data;
        // sum of square of horizontal gradient and square of vertical gradient
        float *delhvsqsum       = (float (*))         ((char*)rgbgreen + sizeof(float) * TS * TS + CLF * 64);
        // gradient based directional weights for interpolation
        float *dirwts0          = (float (*))         ((char*)delhvsqsum + sizeof(float) * TS * TS + CLF * 64);
        float *dirwts1          = (float (*))         ((char*)dirwts0 + sizeof(float) * TS * TS + CLF * 64);
        // vertically interpolated color differences G-R, G-B
        float *vcd              = (float (*))         ((char*)dirwts1 + sizeof(float) * TS * TS + CLF * 64);
        // horizontally interpolated color differences
        float *hcd              = (float (*))         ((char*)vcd + sizeof(float) * TS * TS + CLF * 64);
        // alternative vertical interpolation
        float *vcdalt           = (float (*))         ((char*)hcd + sizeof(float) * TS * TS + CLF * 64);
        // alternative horizontal interpolation
        float *hcdalt           = (float (*))         ((char*)vcdalt + sizeof(float) * TS * TS + CLF * 64);
        // square of average color difference
        float *cddiffsq         = (float (*))         ((char*)hcdalt + sizeof(float) * TS * TS + CLF * 64);
        // weight to give horizontal vs vertical interpolation
        float *hvwt             = (float (*))         ((char*)cddiffsq + sizeof(float) * TS * TS + 2 * CLF * 64);
        // final interpolated color difference
        float (*Dgrb)[TS * TSH] = (float (*)[TS * TSH])vcdalt; // there is no overlap in buffer usage => share
        // gradient in plus (NE/SW) direction
        float *delp             = (float (*))cddiffsq; // there is no overlap in buffer usage => share
        // gradient in minus (NW/SE) direction
        float *delm             = (float (*))         ((char*)delp + sizeof(float) * TS * TSH + CLF * 64);
        // diagonal interpolation of R+B
        float *rbint            = (float (*))delm; // there is no overlap in buffer usage => share
        // horizontal and vertical curvature of interpolated G (used to refine interpolation in Nyquist texture regions)
        s_hv  *Dgrb2            = (s_hv  (*))         ((char*)hvwt + sizeof(float) * TS * TSH + CLF * 64);
        // difference between up/down interpolations of G
        float *dgintv           = (float (*))Dgrb2;   // there is no overlap in buffer usage => share
        // difference between left/right interpolations of G
        float *dginth           = (float (*))         ((char*)dgintv + sizeof(float) * TS * TS + CLF * 64);
        // square of diagonal colour differences
        float *Dgrbsq1m         = (float (*))         ((char*)dginth + sizeof(float) * TS * TS + CLF * 64);
        float *Dgrbsq1p         = (float (*))         ((char*)Dgrbsq1m + sizeof(float) * TS * TSH + CLF * 64);
        // tile raw data
        float *cfa              = (float (*))         ((char*)Dgrbsq1p + sizeof(float) * TS * TSH + CLF * 64);
        // relative weight for combining plus and minus diagonal interpolations
        float *pmwt             = (float (*))delhvsqsum;  // there is no overlap in buffer usage => share
        // interpolated color difference R-B in minus and plus direction
        float *rbm              = (float (*))vcd;  // there is no overlap in buffer usage => share
        float *rbp              = (float (*))         ((char*)rbm + sizeof(float) * TS * TSH + CLF * 64);
        // nyquist texture flag 1=nyquist, 0=not nyquist
        unsigned char *nyquist  = (unsigned char (*)) ((char*)cfa + sizeof(float) * TS * TS + CLF * 64);
        /*
                rgbgreen   = (float (*))         data; //pointers to array
                delhvsqsum = (float (*))         ((char*)rgbgreen + sizeof(float) * TS * TS + CLF * 64);
                dirwts0    = (float (*))         ((char*)delhvsqsum + sizeof(float) * TS * TS + CLF * 64);
                dirwts1    = (float (*))         ((char*)dirwts0 + sizeof(float) * TS * TS + CLF * 64);
                vcd        = (float (*))         ((char*)dirwts1 + sizeof(float) * TS * TS + CLF * 64);
                hcd        = (float (*))         ((char*)vcd + sizeof(float) * TS * TS + CLF * 64);
                vcdalt     = (float (*))         ((char*)hcd + sizeof(float) * TS * TS + CLF * 64);
                hcdalt     = (float (*))         ((char*)vcdalt + sizeof(float) * TS * TS + CLF * 64);
                cddiffsq   = (float (*))         ((char*)hcdalt + sizeof(float) * TS * TS + CLF * 64);
                hvwt       = (float (*))         ((char*)cddiffsq + sizeof(float) * TS * TS + CLF * 64);
                Dgrb       = (float (*)[TS * TSH]) ((char*)hvwt + sizeof(float) * TS * TSH + CLF * 64);
                delp       = (float (*))         ((char*)Dgrb + sizeof(float) * TS * TS + CLF * 64);
                delm       = (float (*))         ((char*)delp + sizeof(float) * TS * TSH + CLF * 64);
                rbint      = (float (*))         ((char*)delm + sizeof(float) * TS * TSH + CLF * 64);
                Dgrb2      = (s_hv  (*))         ((char*)rbint + sizeof(float) * TS * TSH + CLF * 64);
                dgintv     = (float (*))         ((char*)Dgrb2 + sizeof(float) * TS * TS + CLF * 64);
                dginth     = (float (*))         ((char*)dgintv + sizeof(float) * TS * TS + CLF * 64);
                Dgrbsq1m   = (float (*))         ((char*)dginth + sizeof(float) * TS * TS + CLF * 64);
                Dgrbsq1p   = (float (*))         ((char*)Dgrbsq1m + sizeof(float) * TS * TSH + CLF * 64);
                cfa        = (float (*))         ((char*)Dgrbsq1p + sizeof(float) * TS * TSH + CLF * 64);
                pmwt       = (float (*))         ((char*)cfa + sizeof(float) * TS * TS + CLF * 64);
                rbm        = (float (*))         ((char*)pmwt + sizeof(float) * TS * TSH + CLF * 64);
                rbp        = (float (*))         ((char*)rbm + sizeof(float) * TS * TSH + CLF * 64);

                nyquist    = (char (*))          ((char*)rbp + sizeof(float) * TS * TSH + CLF * 64);
        */
#undef CLF

        // Main algorithm: Tile loop

        // Issue 1676
        // use collapse(2) to collapse the 2 loops to one large loop, so there is better scaling
#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2) nowait
#endif

        for (int top = winy - 16; top < winy + height; top += TS - 32)
            for (int left = winx - 16; left < winx + width; left += TS - 32) {
#ifdef __SSE2__
                // Using SSE2 we can zero the memory without cache pollution
                vfloat zerov = ZEROV;

                for(int i = 3 * TSH; i < (TS - 6)*TSH; i += 16) {
                    _mm_stream_ps((float*)&nyquist[i], zerov);
                }

#else
                memset(&nyquist[3 * TSH], 0, sizeof(unsigned char) * (TS - 6) * TSH);
#endif
                //location of tile bottom edge
                const int bottom = min(top + TS, winy + height + 16);
                //location of tile right edge
                const int right  = min(left + TS, winx + width + 16);
                //tile width  (=TS except for right edge of image)
                const int rr1 = bottom - top;
                //tile height (=TS except for bottom edge of image)
                const int cc1 = right - left;
                // bookkeeping for borders
                // min and max row/column in the tile
                int rrmin = top < winy ? 16 : 0;
                int ccmin = left < winx ? 16 : 0;
                int rrmax = bottom > (winy + height) ? winy + height - top : rr1;
                int ccmax = right > (winx + width) ? winx + width - left : cc1;

                // rgb from input CFA data
                // rgb values should be floating point number between 0 and 1
                // after white balance multipliers are applied
                // a 16 pixel border is added to each side of the image
#ifdef __SSE2__
                const vfloat c65535v = F2V( 65535.0f );

                //fill upper border
                if (rrmin > 0) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = ccmin, row = 32 - rr + top; cc < ccmax; cc++) {
                            cfa[rr * TS + cc] = (rawData[row][cc + left]) / 65535.0f;
                            rgbgreen[rr * TS + cc] = cfa[rr * TS + cc];
                        }
                }

                // fill inner part
                for (int rr = rrmin; rr < rrmax; rr++) {
                    int row = rr + top;
                    int cc = ccmin;

                    for (; cc < ccmax - 3; cc += 4) {
                        int indx1 = rr * TS + cc;
                        vfloat tempv = LVFU(rawData[row][cc + left]) / c65535v;
                        STVF(cfa[indx1], tempv );
                        STVF(rgbgreen[indx1], tempv );
                    }

                    for (; cc < ccmax; cc++) {
                        int indx1 = rr * TS + cc;
                        cfa[indx1] = (rawData[row][cc + left]) / 65535.0f;
                        rgbgreen[indx1] = cfa[indx1];
                    }
                }

                //fill lower border
                if (rrmax < rr1) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = ccmin; cc < ccmax; cc += 4) {
                            int indx1 = (rrmax + rr) * TS + cc;
                            vfloat tempv = LVFU(rawData[(winy + height - rr - 2)][left + cc]) / c65535v;
                            STVF(cfa[indx1], tempv );
                            STVF(rgbgreen[indx1], tempv );
                        }
                }

                //fill left border
                if (ccmin > 0) {
                    for (int rr = rrmin; rr < rrmax; rr++)
                        for (int cc = 0, row = rr + top; cc < 16; cc++) {
                            cfa[rr * TS + cc] = (rawData[row][32 - cc + left]) / 65535.0f;
                            rgbgreen[rr * TS + cc] = cfa[rr * TS + cc];
                        }
                }

                //fill right border
                if (ccmax < cc1) {
                    for (int rr = rrmin; rr < rrmax; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[rr * TS + ccmax + cc] = (rawData[(top + rr)][(winx + width - cc - 2)]) / 65535.0f;
                            rgbgreen[rr * TS + ccmax + cc] = cfa[rr * TS + ccmax + cc];
                        }
                }

                //also, fill the image corners
                if (rrmin > 0 && ccmin > 0) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc += 4) {
                            int indx1 = (rr) * TS + cc;
                            vfloat tempv = LVFU(rawData[winy + 32 - rr][winx + 32 - cc]) / c65535v;
                            STVF(cfa[indx1], tempv );
                            STVF(rgbgreen[indx1], tempv );
                        }
                }

                if (rrmax < rr1 && ccmax < cc1) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc += 4) {
                            int indx1 = (rrmax + rr) * TS + ccmax + cc;
                            vfloat tempv = LVFU(rawData[(winy + height - rr - 2)][(winx + width - cc - 2)]) / c65535v;
                            STVFU(cfa[indx1], tempv );
                            STVFU(rgbgreen[indx1], tempv );
                        }
                }

                if (rrmin > 0 && ccmax < cc1) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[(rr)*TS + ccmax + cc] = (rawData[(winy + 32 - rr)][(winx + width - cc - 2)]) / 65535.0f;
                            rgbgreen[(rr)*TS + ccmax + cc] = cfa[(rr) * TS + ccmax + cc];
                        }
                }

                if (rrmax < rr1 && ccmin > 0) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[(rrmax + rr)*TS + cc] = (rawData[(winy + height - rr - 2)][(winx + 32 - cc)]) / 65535.0f;
                            rgbgreen[(rrmax + rr)*TS + cc] = cfa[(rrmax + rr) * TS + cc];
                        }
                }

#else

                for (int rr = rrmin; rr < rrmax; rr++)
                    for (int row = rr + top, cc = ccmin; cc < ccmax; cc++) {
                        int indx1 = rr * TS + cc;
                        cfa[indx1] = (rawData[row][cc + left]) / 65535.0f;
                        rgbgreen[indx1] = cfa[indx1];
                    }

                //fill borders
                if (rrmin > 0) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = ccmin, row = 32 - rr + top; cc < ccmax; cc++) {
                            cfa[rr * TS + cc] = (rawData[row][cc + left]) / 65535.0f;
                            rgbgreen[rr * TS + cc] = cfa[rr * TS + cc];
                        }
                }

                if (rrmax < rr1) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = ccmin; cc < ccmax; cc++) {
                            cfa[(rrmax + rr)*TS + cc] = (rawData[(winy + height - rr - 2)][left + cc]) / 65535.0f;
                            rgbgreen[(rrmax + rr)*TS + cc] = cfa[(rrmax + rr) * TS + cc];
                        }
                }

                if (ccmin > 0) {
                    for (int rr = rrmin; rr < rrmax; rr++)
                        for (int cc = 0, row = rr + top; cc < 16; cc++) {
                            cfa[rr * TS + cc] = (rawData[row][32 - cc + left]) / 65535.0f;
                            rgbgreen[rr * TS + cc] = cfa[rr * TS + cc];
                        }
                }

                if (ccmax < cc1) {
                    for (int rr = rrmin; rr < rrmax; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[rr * TS + ccmax + cc] = (rawData[(top + rr)][(winx + width - cc - 2)]) / 65535.0f;
                            rgbgreen[rr * TS + ccmax + cc] = cfa[rr * TS + ccmax + cc];
                        }
                }

                //also, fill the image corners
                if (rrmin > 0 && ccmin > 0) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[(rr)*TS + cc] = (rawData[winy + 32 - rr][winx + 32 - cc]) / 65535.0f;
                            rgbgreen[(rr)*TS + cc] = cfa[(rr) * TS + cc];
                        }
                }

                if (rrmax < rr1 && ccmax < cc1) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[(rrmax + rr)*TS + ccmax + cc] = (rawData[(winy + height - rr - 2)][(winx + width - cc - 2)]) / 65535.0f;
                            rgbgreen[(rrmax + rr)*TS + ccmax + cc] = cfa[(rrmax + rr) * TS + ccmax + cc];
                        }
                }

                if (rrmin > 0 && ccmax < cc1) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[(rr)*TS + ccmax + cc] = (rawData[(winy + 32 - rr)][(winx + width - cc - 2)]) / 65535.0f;
                            rgbgreen[(rr)*TS + ccmax + cc] = cfa[(rr) * TS + ccmax + cc];
                        }
                }

                if (rrmax < rr1 && ccmin > 0) {
                    for (int rr = 0; rr < 16; rr++)
                        for (int cc = 0; cc < 16; cc++) {
                            cfa[(rrmax + rr)*TS + cc] = (rawData[(winy + height - rr - 2)][(winx + 32 - cc)]) / 65535.0f;
                            rgbgreen[(rrmax + rr)*TS + cc] = cfa[(rrmax + rr) * TS + cc];
                        }
                }

#endif

                //end of border fill
#ifdef __SSE2__
                const vfloat epsv = F2V( eps );

                for (int rr = 2; rr < rr1 - 2; rr++) {
                    for (int indx = rr * TS; indx < rr * TS + cc1; indx += 4) {
                        vfloat delhv = vabsf( LVFU( cfa[indx + 1] ) -  LVFU( cfa[indx - 1] ) );
                        vfloat delvv = vabsf( LVF( cfa[indx + v1] ) -  LVF( cfa[indx - v1] ) );
                        STVF(dirwts1[indx], epsv + vabsf( LVFU( cfa[indx + 2] ) - LVF( cfa[indx] )) + vabsf( LVF( cfa[indx] ) - LVFU( cfa[indx - 2] )) + delhv );
                        STVF(dirwts0[indx], epsv + vabsf( LVF( cfa[indx + v2] ) - LVF( cfa[indx] )) + vabsf( LVF( cfa[indx] ) - LVF( cfa[indx - v2] )) + delvv );
                        STVF(delhvsqsum[indx], SQRV(delhv) + SQRV(delvv));
                    }
                }

#else

                for (int rr = 2; rr < rr1 - 2; rr++)
                    for (int cc = 2, indx = (rr) * TS + cc; cc < cc1 - 2; cc++, indx++) {
                        // horizontal and vedrtical gradient
                        float delh = fabsf(cfa[indx + 1] - cfa[indx - 1]);
                        float delv = fabsf(cfa[indx + v1] - cfa[indx - v1]);
                        dirwts0[indx] = eps + fabsf(cfa[indx + v2] - cfa[indx]) + fabsf(cfa[indx] - cfa[indx - v2]) + delv;
                        dirwts1[indx] = eps + fabsf(cfa[indx + 2] - cfa[indx]) + fabsf(cfa[indx] - cfa[indx - 2]) + delh; //+fabsf(cfa[indx+2]-cfa[indx-2]);
                        delhvsqsum[indx] = SQR(delh) + SQR(delv);
                    }

#endif

                //interpolate vertical and horizontal color differences
#ifdef __SSE2__
                vfloat sgnv;

                if( !(FC(4, 4) & 1) ) {
                    sgnv = _mm_set_ps( 1.0f, -1.0f, 1.0f, -1.0f );
                } else {
                    sgnv = _mm_set_ps( -1.0f, 1.0f, -1.0f, 1.0f );
                }

                vfloat  zd5v = F2V( 0.5f );
                vfloat  onev = F2V( 1.0f );
                vfloat  arthreshv = F2V( arthresh );
                vfloat  clip_pt8v = F2V( clip_pt8 );

                for (int rr = 4; rr < rr1 - 4; rr++) {
                    sgnv = -sgnv;

                    for (int indx = rr * TS + 4; indx < rr * TS + cc1 - 7; indx += 4) {
                        //color ratios in each cardinal direction
                        vfloat cfav = LVF(cfa[indx]);
                        vfloat cruv = LVF(cfa[indx - v1]) * (LVF(dirwts0[indx - v2]) + LVF(dirwts0[indx])) / (LVF(dirwts0[indx - v2]) * (epsv + cfav) + LVF(dirwts0[indx]) * (epsv + LVF(cfa[indx - v2])));
                        vfloat crdv = LVF(cfa[indx + v1]) * (LVF(dirwts0[indx + v2]) + LVF(dirwts0[indx])) / (LVF(dirwts0[indx + v2]) * (epsv + cfav) + LVF(dirwts0[indx]) * (epsv + LVF(cfa[indx + v2])));
                        vfloat crlv = LVFU(cfa[indx - 1]) * (LVFU(dirwts1[indx - 2]) + LVF(dirwts1[indx])) / (LVFU(dirwts1[indx - 2]) * (epsv + cfav) + LVF(dirwts1[indx]) * (epsv + LVFU(cfa[indx - 2])));
                        vfloat crrv = LVFU(cfa[indx + 1]) * (LVFU(dirwts1[indx + 2]) + LVF(dirwts1[indx])) / (LVFU(dirwts1[indx + 2]) * (epsv + cfav) + LVF(dirwts1[indx]) * (epsv + LVFU(cfa[indx + 2])));

                        vfloat guhav = LVF(cfa[indx - v1]) + zd5v * (cfav - LVF(cfa[indx - v2]));
                        vfloat gdhav = LVF(cfa[indx + v1]) + zd5v * (cfav - LVF(cfa[indx + v2]));
                        vfloat glhav = LVFU(cfa[indx - 1]) + zd5v * (cfav - LVFU(cfa[indx - 2]));
                        vfloat grhav = LVFU(cfa[indx + 1]) + zd5v * (cfav - LVFU(cfa[indx + 2]));

                        vfloat guarv = vself(vmaskf_lt(vabsf(onev - cruv), arthreshv), cfav * cruv, guhav);
                        vfloat gdarv = vself(vmaskf_lt(vabsf(onev - crdv), arthreshv), cfav * crdv, gdhav);
                        vfloat glarv = vself(vmaskf_lt(vabsf(onev - crlv), arthreshv), cfav * crlv, glhav);
                        vfloat grarv = vself(vmaskf_lt(vabsf(onev - crrv), arthreshv), cfav * crrv, grhav);

                        vfloat hwtv = LVFU(dirwts1[indx - 1]) / (LVFU(dirwts1[indx - 1]) + LVFU(dirwts1[indx + 1]));
                        vfloat vwtv = LVF(dirwts0[indx - v1]) / (LVF(dirwts0[indx + v1]) + LVF(dirwts0[indx - v1]));

                        //interpolated G via adaptive weights of cardinal evaluations
                        vfloat Ginthhav = vintpf(hwtv, grhav, glhav);
                        vfloat Gintvhav = vintpf(vwtv, gdhav, guhav);

                        //interpolated color differences
                        vfloat hcdaltv = sgnv * (Ginthhav - cfav);
                        vfloat vcdaltv = sgnv * (Gintvhav - cfav);
                        STVF(hcdalt[indx], hcdaltv);
                        STVF(vcdalt[indx], vcdaltv);

                        vmask clipmask = vorm( vorm( vmaskf_gt( cfav, clip_pt8v ), vmaskf_gt( Gintvhav, clip_pt8v ) ), vmaskf_gt( Ginthhav, clip_pt8v ));
                        guarv = vself( clipmask, guhav, guarv);
                        gdarv = vself( clipmask, gdhav, gdarv);
                        glarv = vself( clipmask, glhav, glarv);
                        grarv = vself( clipmask, grhav, grarv);
                        STVF(vcd[indx], vself( clipmask, vcdaltv, sgnv * (vintpf(vwtv, gdarv, guarv) - cfav)));
                        STVF(hcd[indx], vself( clipmask, hcdaltv, sgnv * (vintpf(hwtv, grarv, glarv) - cfav)));
                        //differences of interpolations in opposite directions

                        STVF(dgintv[indx], vminf(SQRV(guhav - gdhav), SQRV(guarv - gdarv)));
                        STVF(dginth[indx], vminf(SQRV(glhav - grhav), SQRV(glarv - grarv)));

                    }
                }

#else

                for (int rr = 4; rr < rr1 - 4; rr++) {
                    bool fcswitch = FC(rr, 4) & 1;

                    for (int cc = 4, indx = rr * TS + cc; cc < cc1 - 4; cc++, indx++) {

                        //color ratios in each cardinal direction
                        float cru = cfa[indx - v1] * (dirwts0[indx - v2] + dirwts0[indx]) / (dirwts0[indx - v2] * (eps + cfa[indx]) + dirwts0[indx] * (eps + cfa[indx - v2]));
                        float crd = cfa[indx + v1] * (dirwts0[indx + v2] + dirwts0[indx]) / (dirwts0[indx + v2] * (eps + cfa[indx]) + dirwts0[indx] * (eps + cfa[indx + v2]));
                        float crl = cfa[indx - 1] * (dirwts1[indx - 2] + dirwts1[indx]) / (dirwts1[indx - 2] * (eps + cfa[indx]) + dirwts1[indx] * (eps + cfa[indx - 2]));
                        float crr = cfa[indx + 1] * (dirwts1[indx + 2] + dirwts1[indx]) / (dirwts1[indx + 2] * (eps + cfa[indx]) + dirwts1[indx] * (eps + cfa[indx + 2]));

                        //G interpolated in vert/hor directions using Hamilton-Adams method
                        float guha = cfa[indx - v1] + xdiv2f(cfa[indx] - cfa[indx - v2]);
                        float gdha = cfa[indx + v1] + xdiv2f(cfa[indx] - cfa[indx + v2]);
                        float glha = cfa[indx - 1] + xdiv2f(cfa[indx] - cfa[indx - 2]);
                        float grha = cfa[indx + 1] + xdiv2f(cfa[indx] - cfa[indx + 2]);

                        //G interpolated in vert/hor directions using adaptive ratios
                        float guar, gdar, glar, grar;

                        if (fabsf(1.0f - cru) < arthresh) {
                            guar = cfa[indx] * cru;
                        } else {
                            guar = guha;
                        }

                        if (fabsf(1.0f - crd) < arthresh) {
                            gdar = cfa[indx] * crd;
                        } else {
                            gdar = gdha;
                        }

                        if (fabsf(1.0f - crl) < arthresh) {
                            glar = cfa[indx] * crl;
                        } else {
                            glar = glha;
                        }

                        if (fabsf(1.0f - crr) < arthresh) {
                            grar = cfa[indx] * crr;
                        } else {
                            grar = grha;
                        }

                        //adaptive weights for vertical/horizontal directions
                        float hwt = dirwts1[indx - 1] / (dirwts1[indx - 1] + dirwts1[indx + 1]);
                        float vwt = dirwts0[indx - v1] / (dirwts0[indx + v1] + dirwts0[indx - v1]);

                        //interpolated G via adaptive weights of cardinal evaluations
                        float Gintvha = vwt * gdha + (1.0f - vwt) * guha;
                        float Ginthha = hwt * grha + (1.0f - hwt) * glha;

                        //interpolated color differences
                        if (fcswitch) {
                            vcd[indx] = cfa[indx] - (vwt * gdar + (1.0f - vwt) * guar);
                            hcd[indx] = cfa[indx] - (hwt * grar + (1.0f - hwt) * glar);
                            vcdalt[indx] = cfa[indx] - Gintvha;
                            hcdalt[indx] = cfa[indx] - Ginthha;
                        } else {
                            //interpolated color differences
                            vcd[indx] = (vwt * gdar + (1.0f - vwt) * guar) - cfa[indx];
                            hcd[indx] = (hwt * grar + (1.0f - hwt) * glar) - cfa[indx];
                            vcdalt[indx] = Gintvha - cfa[indx];
                            hcdalt[indx] = Ginthha - cfa[indx];
                        }

                        fcswitch = !fcswitch;

                        if (cfa[indx] > clip_pt8 || Gintvha > clip_pt8 || Ginthha > clip_pt8) {
                            //use HA if highlights are (nearly) clipped
                            guar = guha;
                            gdar = gdha;
                            glar = glha;
                            grar = grha;
                            vcd[indx] = vcdalt[indx];
                            hcd[indx] = hcdalt[indx];
                        }

                        //differences of interpolations in opposite directions
                        dgintv[indx] = min(SQR(guha - gdha), SQR(guar - gdar));
                        dginth[indx] = min(SQR(glha - grha), SQR(glar - grar));

                    }


                }

#endif



#ifdef __SSE2__
                vfloat  clip_ptv = F2V( clip_pt );
                vfloat  sgn3v;

                if( !(FC(4, 4) & 1) ) {
                    sgnv = _mm_set_ps( 1.0f, -1.0f, 1.0f, -1.0f );
                } else {
                    sgnv = _mm_set_ps( -1.0f, 1.0f, -1.0f, 1.0f );
                }

                sgn3v = sgnv + sgnv + sgnv;

                for (int rr = 4; rr < rr1 - 4; rr++) {
                    vfloat nsgnv = sgnv;
                    sgnv = -sgnv;
                    sgn3v = -sgn3v;

                    for (int indx = rr * TS + 4; indx < rr * TS + cc1 - 4; indx += 4) {
                        vfloat hcdv = LVF( hcd[indx] );
                        vfloat hcdvarv = SQRV(LVFU(hcd[indx - 2]) - hcdv) + SQRV(LVFU(hcd[indx - 2]) - LVFU(hcd[indx + 2])) + SQRV(hcdv - LVFU(hcd[indx + 2]));
                        vfloat hcdaltv = LVF( hcdalt[indx] );
                        vfloat hcdaltvarv = SQRV(LVFU(hcdalt[indx - 2]) - hcdaltv) + SQRV(LVFU(hcdalt[indx - 2]) - LVFU(hcdalt[indx + 2])) + SQRV(hcdaltv - LVFU(hcdalt[indx + 2]));
                        vfloat vcdv = LVF( vcd[indx] );
                        vfloat vcdvarv = SQRV(LVF(vcd[indx - v2]) - vcdv) + SQRV(LVF(vcd[indx - v2]) - LVF(vcd[indx + v2])) + SQRV(vcdv - LVF(vcd[indx + v2]));
                        vfloat vcdaltv = LVF( vcdalt[indx] );
                        vfloat vcdaltvarv = SQRV(LVF(vcdalt[indx - v2]) - vcdaltv) + SQRV(LVF(vcdalt[indx - v2]) - LVF(vcdalt[indx + v2])) + SQRV(vcdaltv - LVF(vcdalt[indx + v2]));

                        //choose the smallest variance; this yields a smoother interpolation
                        hcdv = vself( vmaskf_lt( hcdaltvarv, hcdvarv ), hcdaltv, hcdv);
                        vcdv = vself( vmaskf_lt( vcdaltvarv, vcdvarv ), vcdaltv, vcdv);

                        vfloat Ginthv = sgnv * hcdv + LVF( cfa[indx] );
                        vfloat temp2v = sgn3v * hcdv;
                        vfloat hwtv = onev + temp2v / ( epsv + Ginthv + LVF( cfa[indx]));
                        vmask hcdmask = vmaskf_gt( nsgnv * hcdv, ZEROV );
                        vfloat hcdoldv = hcdv;
                        vfloat tempv = nsgnv * (LVF(cfa[indx]) - ULIMV( Ginthv, LVFU(cfa[indx - 1]), LVFU(cfa[indx + 1]) ));
                        hcdv = vself( vmaskf_lt( temp2v, -(LVF(cfa[indx]) + Ginthv)), tempv, vintpf(hwtv, hcdv, tempv));
                        hcdv = vself( hcdmask, hcdv, hcdoldv );
                        hcdv = vself( vmaskf_gt( Ginthv, clip_ptv), tempv, hcdv);
                        STVF(hcd[indx], hcdv);

                        vfloat Gintvv = sgnv * vcdv + LVF( cfa[indx] );
                        temp2v = sgn3v * vcdv;
                        vfloat vwtv = onev + temp2v / ( epsv + Gintvv + LVF( cfa[indx]));
                        vmask vcdmask = vmaskf_gt( nsgnv * vcdv, ZEROV );
                        vfloat vcdoldv = vcdv;
                        tempv = nsgnv * (LVF(cfa[indx]) - ULIMV( Gintvv, LVF(cfa[indx - v1]), LVF(cfa[indx + v1]) ));
                        vcdv = vself( vmaskf_lt( temp2v, -(LVF(cfa[indx]) + Gintvv)), tempv, vintpf(vwtv, vcdv, tempv));
                        vcdv = vself( vcdmask, vcdv, vcdoldv );
                        vcdv = vself( vmaskf_gt( Gintvv, clip_ptv), tempv, vcdv);
                        STVF(vcd[indx], vcdv);
                        STVFU(cddiffsq[indx], SQRV(vcdv - hcdv));
                    }

                }

#else

                for (int rr = 4; rr < rr1 - 4; rr++) {
                    for (int cc = 4, indx = rr * TS + cc, c = FC(rr, cc) & 1; cc < cc1 - 4; cc++, indx++) {
                        float hcdvar = 3.0f * (SQR(hcd[indx - 2]) + SQR(hcd[indx]) + SQR(hcd[indx + 2])) - SQR(hcd[indx - 2] + hcd[indx] + hcd[indx + 2]);
                        float hcdaltvar = 3.0f * (SQR(hcdalt[indx - 2]) + SQR(hcdalt[indx]) + SQR(hcdalt[indx + 2])) - SQR(hcdalt[indx - 2] + hcdalt[indx] + hcdalt[indx + 2]);
                        float vcdvar = 3.0f * (SQR(vcd[indx - v2]) + SQR(vcd[indx]) + SQR(vcd[indx + v2])) - SQR(vcd[indx - v2] + vcd[indx] + vcd[indx + v2]);
                        float vcdaltvar = 3.0f * (SQR(vcdalt[indx - v2]) + SQR(vcdalt[indx]) + SQR(vcdalt[indx + v2])) - SQR(vcdalt[indx - v2] + vcdalt[indx] + vcdalt[indx + v2]);

                        //choose the smallest variance; this yields a smoother interpolation
                        if (hcdaltvar < hcdvar) {
                            hcd[indx] = hcdalt[indx];
                        }

                        if (vcdaltvar < vcdvar) {
                            vcd[indx] = vcdalt[indx];
                        }

                        //bound the interpolation in regions of high saturation

                        //vertical and horizontal G interpolations
                        float Gintv, Ginth;

                        if (c) {//G site
                            Ginth = -hcd[indx] + cfa[indx]; //R or B
                            Gintv = -vcd[indx] + cfa[indx]; //B or R

                            if (hcd[indx] > 0) {
                                if (3.0f * hcd[indx] > (Ginth + cfa[indx])) {
                                    hcd[indx] = -ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]) + cfa[indx];
                                } else {
                                    float hwt = 1.0f - 3.0f * hcd[indx] / (eps + Ginth + cfa[indx]);
                                    hcd[indx] = hwt * hcd[indx] + (1.0f - hwt) * (-ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]) + cfa[indx]);
                                }
                            }

                            if (vcd[indx] > 0) {
                                if (3.0f * vcd[indx] > (Gintv + cfa[indx])) {
                                    vcd[indx] = -ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]) + cfa[indx];
                                } else {
                                    float vwt = 1.0f - 3.0f * vcd[indx] / (eps + Gintv + cfa[indx]);
                                    vcd[indx] = vwt * vcd[indx] + (1.0f - vwt) * (-ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]) + cfa[indx]);
                                }
                            }

                            if (Ginth > clip_pt) {
                                hcd[indx] = -ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]) + cfa[indx];    //for RT implementation
                            }

                            if (Gintv > clip_pt) {
                                vcd[indx] = -ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]) + cfa[indx];
                            }

                            //if (Ginth > pre_mul[c]) hcd[indx]=-ULIM(Ginth,cfa[indx-1],cfa[indx+1])+cfa[indx];//for dcraw implementation
                            //if (Gintv > pre_mul[c]) vcd[indx]=-ULIM(Gintv,cfa[indx-v1],cfa[indx+v1])+cfa[indx];

                        } else {//R or B site

                            Ginth = hcd[indx] + cfa[indx]; //interpolated G
                            Gintv = vcd[indx] + cfa[indx];

                            if (hcd[indx] < 0) {
                                if (3.0f * hcd[indx] < -(Ginth + cfa[indx])) {
                                    hcd[indx] = ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]) - cfa[indx];
                                } else {
                                    float hwt = 1.0f + 3.0f * hcd[indx] / (eps + Ginth + cfa[indx]);
                                    hcd[indx] = hwt * hcd[indx] + (1.0f - hwt) * (ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]) - cfa[indx]);
                                }
                            }

                            if (vcd[indx] < 0) {
                                if (3.0f * vcd[indx] < -(Gintv + cfa[indx])) {
                                    vcd[indx] = ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]) - cfa[indx];
                                } else {
                                    float vwt = 1.0f + 3.0f * vcd[indx] / (eps + Gintv + cfa[indx]);
                                    vcd[indx] = vwt * vcd[indx] + (1.0f - vwt) * (ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]) - cfa[indx]);
                                }
                            }

                            if (Ginth > clip_pt) {
                                hcd[indx] = ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]) - cfa[indx];    //for RT implementation
                            }

                            if (Gintv > clip_pt) {
                                vcd[indx] = ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]) - cfa[indx];
                            }

                            cddiffsq[indx] = SQR(vcd[indx] - hcd[indx]);
                        }

                        c = !c;
                    }
                }

#endif



#ifdef __SSE2__
                vfloat  epssqv = F2V( epssq );

                for (int rr = 6; rr < rr1 - 6; rr++) {
                    for (int indx = rr * TS + 6 + (FC(rr, 2) & 1); indx < rr * TS + cc1 - 6; indx += 8) {
                        //compute colour difference variances in cardinal directions
                        vfloat tempv = LC2VFU(vcd[indx]);
                        vfloat uavev = tempv + LC2VFU(vcd[indx - v1]) + LC2VFU(vcd[indx - v2]) + LC2VFU(vcd[indx - v3]);
                        vfloat davev = tempv + LC2VFU(vcd[indx + v1]) + LC2VFU(vcd[indx + v2]) + LC2VFU(vcd[indx + v3]);
                        vfloat Dgrbvvaruv = SQRV(tempv - uavev) + SQRV(LC2VFU(vcd[indx - v1]) - uavev) + SQRV(LC2VFU(vcd[indx - v2]) - uavev) + SQRV(LC2VFU(vcd[indx - v3]) - uavev);
                        vfloat Dgrbvvardv = SQRV(tempv - davev) + SQRV(LC2VFU(vcd[indx + v1]) - davev) + SQRV(LC2VFU(vcd[indx + v2]) - davev) + SQRV(LC2VFU(vcd[indx + v3]) - davev);

                        vfloat hwtv = LC2VFU(dirwts1[indx - 1]) / (LC2VFU(dirwts1[indx - 1]) + LC2VFU(dirwts1[indx + 1]));
                        vfloat vwtv = LC2VFU(dirwts0[indx - v1]) / (LC2VFU(dirwts0[indx + v1]) + LC2VFU(dirwts0[indx - v1]));

                        tempv = LC2VFU(hcd[indx]);
                        vfloat lavev = tempv + vaddc2vfu(hcd[indx - 3]) + LC2VFU(hcd[indx - 1]);
                        vfloat ravev = tempv + vaddc2vfu(hcd[indx + 1]) + LC2VFU(hcd[indx + 3]);

                        vfloat Dgrbhvarlv = SQRV(tempv - lavev) + SQRV(LC2VFU(hcd[indx - 1]) - lavev) + SQRV(LC2VFU(hcd[indx - 2]) - lavev) + SQRV(LC2VFU(hcd[indx - 3]) - lavev);
                        vfloat Dgrbhvarrv = SQRV(tempv - ravev) + SQRV(LC2VFU(hcd[indx + 1]) - ravev) + SQRV(LC2VFU(hcd[indx + 2]) - ravev) + SQRV(LC2VFU(hcd[indx + 3]) - ravev);


                        vfloat vcdvarv = epssqv + vintpf(vwtv, Dgrbvvardv, Dgrbvvaruv);
                        vfloat hcdvarv = epssqv + vintpf(hwtv, Dgrbhvarrv, Dgrbhvarlv);

                        //compute fluctuations in up/down and left/right interpolations of colors
                        Dgrbvvaruv = LC2VFU(dgintv[indx - v1]) + LC2VFU(dgintv[indx - v2]);
                        Dgrbvvardv = LC2VFU(dgintv[indx + v1]) + LC2VFU(dgintv[indx + v2]);

                        Dgrbhvarlv = vaddc2vfu(dginth[indx - 2]);
                        Dgrbhvarrv = vaddc2vfu(dginth[indx + 1]);

                        vfloat vcdvar1v = epssqv + LC2VFU(dgintv[indx]) + vintpf(vwtv, Dgrbvvardv, Dgrbvvaruv);
                        vfloat hcdvar1v = epssqv + LC2VFU(dginth[indx]) + vintpf(hwtv, Dgrbhvarrv, Dgrbhvarlv);

                        //determine adaptive weights for G interpolation
                        vfloat varwtv = hcdvarv / (vcdvarv + hcdvarv);
                        vfloat diffwtv = hcdvar1v / (vcdvar1v + hcdvar1v);

                        //if both agree on interpolation direction, choose the one with strongest directional discrimination;
                        //otherwise, choose the u/d and l/r difference fluctuation weights
                        vmask decmask = vandm( vmaskf_gt( (zd5v - varwtv) * (zd5v - diffwtv), ZEROV ), vmaskf_lt( vabsf( zd5v - diffwtv), vabsf( zd5v - varwtv) ) );
                        STVFU(hvwt[indx >> 1], vself( decmask, varwtv, diffwtv));
                    }
                }

#else

                for (int rr = 6; rr < rr1 - 6; rr++) {
                    for (int cc = 6 + (FC(rr, 2) & 1), indx = rr * TS + cc; cc < cc1 - 6; cc += 2, indx += 2) {

                        //compute color difference variances in cardinal directions

                        float uave = vcd[indx] + vcd[indx - v1] + vcd[indx - v2] + vcd[indx - v3];
                        float dave = vcd[indx] + vcd[indx + v1] + vcd[indx + v2] + vcd[indx + v3];
                        float lave = hcd[indx] + hcd[indx - 1] + hcd[indx - 2] + hcd[indx - 3];
                        float rave = hcd[indx] + hcd[indx + 1] + hcd[indx + 2] + hcd[indx + 3];

                        //color difference (G-R or G-B) variance in up/down/left/right directions
                        float Dgrbvvaru = SQR(vcd[indx] - uave) + SQR(vcd[indx - v1] - uave) + SQR(vcd[indx - v2] - uave) + SQR(vcd[indx - v3] - uave);
                        float Dgrbvvard = SQR(vcd[indx] - dave) + SQR(vcd[indx + v1] - dave) + SQR(vcd[indx + v2] - dave) + SQR(vcd[indx + v3] - dave);
                        float Dgrbhvarl = SQR(hcd[indx] - lave) + SQR(hcd[indx - 1] - lave) + SQR(hcd[indx - 2] - lave) + SQR(hcd[indx - 3] - lave);
                        float Dgrbhvarr = SQR(hcd[indx] - rave) + SQR(hcd[indx + 1] - rave) + SQR(hcd[indx + 2] - rave) + SQR(hcd[indx + 3] - rave);

                        float hwt = dirwts1[indx - 1] / (dirwts1[indx - 1] + dirwts1[indx + 1]);
                        float vwt = dirwts0[indx - v1] / (dirwts0[indx + v1] + dirwts0[indx - v1]);

                        float vcdvar = epssq + vwt * Dgrbvvard + (1.0f - vwt) * Dgrbvvaru;
                        float hcdvar = epssq + hwt * Dgrbhvarr + (1.0f - hwt) * Dgrbhvarl;

                        //compute fluctuations in up/down and left/right interpolations of colors
                        Dgrbvvaru = (dgintv[indx]) + (dgintv[indx - v1]) + (dgintv[indx - v2]);
                        Dgrbvvard = (dgintv[indx]) + (dgintv[indx + v1]) + (dgintv[indx + v2]);
                        Dgrbhvarl = (dginth[indx]) + (dginth[indx - 1]) + (dginth[indx - 2]);
                        Dgrbhvarr = (dginth[indx]) + (dginth[indx + 1]) + (dginth[indx + 2]);

                        float vcdvar1 = epssq + vwt * Dgrbvvard + (1.0f - vwt) * Dgrbvvaru;
                        float hcdvar1 = epssq + hwt * Dgrbhvarr + (1.0f - hwt) * Dgrbhvarl;

                        //determine adaptive weights for G interpolation
                        float varwt = hcdvar / (vcdvar + hcdvar);
                        float diffwt = hcdvar1 / (vcdvar1 + hcdvar1);

                        //if both agree on interpolation direction, choose the one with strongest directional discrimination;
                        //otherwise, choose the u/d and l/r difference fluctuation weights
                        if ((0.5 - varwt) * (0.5 - diffwt) > 0 && fabsf(0.5 - diffwt) < fabsf(0.5 - varwt)) {
                            hvwt[indx >> 1] = varwt;
                        } else {
                            hvwt[indx >> 1] = diffwt;
                        }

                    }
                }

#endif


                // Nyquist test
                int nystartrow = 0;
                int nyendrow = 0;
                int nystartcol = TS + 1;
                int nyendcol = 0;

                for (int rr = 6; rr < rr1 - 6; rr++) {
                    for (int cc = 6 + (FC(rr, 2) & 1), indx = rr * TS + cc; cc < cc1 - 6; cc += 2, indx += 2) {

                        //nyquist texture test: ask if difference of vcd compared to hcd is larger or smaller than RGGB gradients
                        // TODO_INGO: currently this part needs 10 float mults, 36 float adds, 4 int mults and 44 int adds for every second pixel
                        // it reads 304 bytes for every second pixel and writes <= 1 byte for every second pixel
                        // a precalculated vectorized version could do this with 1/4 of the operations
                        // but it would read 304 bytes for every second pixel and write 8 bytes for every second pixel for the precalculation
                        // (though the vectorized read should be faster than the scalar version)
                        // and read 8 bytes for every second pixel and write 1 byte for every second pixel for final calculation (maybe this last step can be avoided too)
                        float nyqtest1 = gaussodd[0] * cddiffsq[indx] +
                                         gaussodd[1] * (cddiffsq[(indx - m1)] + cddiffsq[(indx + p1)] +
                                                        cddiffsq[(indx - p1)] + cddiffsq[(indx + m1)]) +
                                         gaussodd[2] * (cddiffsq[(indx - v2)] + cddiffsq[(indx - 2)] +
                                                        cddiffsq[(indx + 2)] + cddiffsq[(indx + v2)]) +
                                         gaussodd[3] * (cddiffsq[(indx - m2)] + cddiffsq[(indx + p2)] +
                                                        cddiffsq[(indx - p2)] + cddiffsq[(indx + m2)]);
                        float nyqtest2 = gaussgrad[0] * delhvsqsum[indx] +
                                         gaussgrad[1] * (delhvsqsum[indx - v1] + delhvsqsum[indx + 1] +
                                                         delhvsqsum[indx - 1] + delhvsqsum[indx + v1]) +
                                         gaussgrad[2] * (delhvsqsum[indx - m1] + delhvsqsum[indx + p1] +
                                                         delhvsqsum[indx - p1] + delhvsqsum[indx + m1]) +
                                         gaussgrad[3] * (delhvsqsum[indx - v2] + delhvsqsum[indx - 2] +
                                                         delhvsqsum[indx + 2] + delhvsqsum[indx + v2]) +
                                         gaussgrad[4] * (delhvsqsum[indx - 2 * TS - 1] + delhvsqsum[indx - 2 * TS + 1] +
                                                         delhvsqsum[indx - TS - 2] + delhvsqsum[indx - TS + 2] +
                                                         delhvsqsum[indx + TS - 2] + delhvsqsum[indx + TS + 2] +
                                                         delhvsqsum[indx + 2 * TS - 1] + delhvsqsum[indx + 2 * TS + 1]) +
                                         gaussgrad[5] * (delhvsqsum[indx - m2] + delhvsqsum[indx + p2] +
                                                         delhvsqsum[indx - p2] + delhvsqsum[indx + m2]);


                        if(nyqtest1 > nyqtest2) {
                            nyquist[indx >> 1] = 1;    //nyquist=1 for nyquist region
                            nystartrow = nystartrow ? nystartrow : rr;
                            nyendrow = rr;
                            nystartcol = nystartcol > cc ? cc : nystartcol;
                            nyendcol = nyendcol < cc ? cc : nyendcol;
                        }
                    }
                }


                bool doNyquist = nystartrow != nyendrow && nystartcol != nyendcol;

                if(doNyquist) {
                    nyendrow ++; // because of < condition
                    nyendcol ++; // because of < condition
                    nystartcol -= (nystartcol & 1);
                    nystartrow = std::max(8, nystartrow);
                    nyendrow = std::min(rr1 - 8, nyendrow);
                    nystartcol = std::max(8, nystartcol);
                    nyendcol = std::min(cc1 - 8, nyendcol);

                    for (int rr = nystartrow; rr < nyendrow; rr++) {
                        for (int indx = rr * TS + nystartcol + (FC(rr, 2) & 1); indx < rr * TS + nyendcol; indx += 2) {
                            // TODO_INGO: if you look at the comments below, it does not seem to be correct to include nyquist[indx >> 1] into the summation
                            // Also this implementation has loop dependencies, which are not correct IMHO
                            // An implementation which uses a second buffer could avoid this dependencies and could be vectorized by factor 16 too (we're working with single bytes here)
                            // That would lead to differences in output compared to current code, but also would lead to more consistent output when changing TS
                            unsigned int nyquistneighbours = (nyquist[(indx - v2) >> 1] + nyquist[(indx - m1) >> 1] + nyquist[(indx + p1) >> 1] +
                                                              nyquist[(indx - 2) >> 1] + nyquist[indx >> 1] + nyquist[(indx + 2) >> 1] +
                                                              nyquist[(indx - p1) >> 1] + nyquist[(indx + m1) >> 1] + nyquist[(indx + v2) >> 1]);

                            //if most of your neighbours are named Nyquist, it's likely that you're one too
                            if (nyquistneighbours > 4) {
                                nyquist[indx >> 1] = 1;
                            }

                            //or not
                            if (nyquistneighbours < 4) {
                                nyquist[indx >> 1] = 0;
                            }
                        }
                    }

                    // end of Nyquist test

                    // in areas of Nyquist texture, do area interpolation
                    for (int rr = nystartrow; rr < nyendrow; rr++)
                        for (int indx = rr * TS + nystartcol + (FC(rr, 2) & 1); indx < rr * TS + nyendcol; indx += 2) {

                            if (nyquist[indx >> 1]) {
                                // area interpolation

                                float sumcfa = 0.f, sumh = 0.f, sumv = 0.f, sumsqh = 0.f, sumsqv = 0.f, areawt = 0.f;

                                for (int i = -6; i < 7; i += 2) {
                                    int indx1 = indx + (i * TS) - 6;

                                    for (int j = -6; j < 7; j += 2, indx1 += 2) {

                                        if (nyquist[indx1 >> 1]) {
                                            float cfatemp = cfa[indx1];
                                            sumcfa += cfatemp;
                                            sumh += (cfa[indx1 - 1] + cfa[indx1 + 1]);
                                            sumv += (cfa[indx1 - v1] + cfa[indx1 + v1]);
                                            sumsqh += SQR(cfatemp - cfa[indx1 - 1]) + SQR(cfatemp - cfa[indx1 + 1]);
                                            sumsqv += SQR(cfatemp - cfa[indx1 - v1]) + SQR(cfatemp - cfa[indx1 + v1]);
                                            areawt += 1;
                                        }
                                    }
                                }

                                //horizontal and vertical color differences, and adaptive weight
                                sumh = sumcfa - xdiv2f(sumh);
                                sumv = sumcfa - xdiv2f(sumv);
                                sumsqh = xdiv2f(sumsqh);
                                sumsqv = xdiv2f(sumsqv);
                                float hcdvar = epssq + fabsf(areawt * sumsqh - sumh * sumh);
                                float vcdvar = epssq + fabsf(areawt * sumsqv - sumv * sumv);
                                hvwt[indx >> 1] = hcdvar / (vcdvar + hcdvar);

                                // end of area interpolation

                            }
                        }
                }


                //populate G at R/B sites
                for (int rr = 8; rr < rr1 - 8; rr++)
                    for (int indx = rr * TS + 8 + (FC(rr, 2) & 1); indx < rr * TS + cc1 - 8; indx += 2) {

                        //first ask if one gets more directional discrimination from nearby B/R sites
                        float hvwtalt = xdivf(hvwt[(indx - m1) >> 1] + hvwt[(indx + p1) >> 1] + hvwt[(indx - p1) >> 1] + hvwt[(indx + m1) >> 1], 2);

                        hvwt[indx >> 1] = fabsf(0.5f - hvwt[indx >> 1]) < fabsf(0.5f - hvwtalt) ? hvwtalt : hvwt[indx >> 1];
                        //a better result was obtained from the neighbours

                        Dgrb[0][indx >> 1] = intp(hvwt[indx >> 1], vcd[indx], hcd[indx]); //evaluate color differences

                        rgbgreen[indx] = cfa[indx] + Dgrb[0][indx >> 1]; //evaluate G (finally!)

                        //local curvature in G (preparation for nyquist refinement step)
                        Dgrb2[indx >> 1].h = nyquist[indx >> 1] ? SQR(rgbgreen[indx] - xdiv2f(rgbgreen[indx - 1] + rgbgreen[indx + 1])) : 0.f;
                        Dgrb2[indx >> 1].v = nyquist[indx >> 1] ? SQR(rgbgreen[indx] - xdiv2f(rgbgreen[indx - v1] + rgbgreen[indx + v1])) : 0.f;
                    }


                //end of standard interpolation

                // refine Nyquist areas using G curvatures
                if(doNyquist) {
                    for (int rr = nystartrow; rr < nyendrow; rr++)
                        // TODO_INGO: maybe this part is also worth vectorizing using _mm_movemask_ps
                        for (int indx = rr * TS + nystartcol + (FC(rr, 2) & 1); indx < rr * TS + nyendcol; indx += 2) {

                            if (nyquist[indx >> 1]) {
                                //local averages (over Nyquist pixels only) of G curvature squared
                                float gvarh = epssq + (gquinc[0] * Dgrb2[indx >> 1].h +
                                                       gquinc[1] * (Dgrb2[(indx - m1) >> 1].h + Dgrb2[(indx + p1) >> 1].h + Dgrb2[(indx - p1) >> 1].h + Dgrb2[(indx + m1) >> 1].h) +
                                                       gquinc[2] * (Dgrb2[(indx - v2) >> 1].h + Dgrb2[(indx - 2) >> 1].h + Dgrb2[(indx + 2) >> 1].h + Dgrb2[(indx + v2) >> 1].h) +
                                                       gquinc[3] * (Dgrb2[(indx - m2) >> 1].h + Dgrb2[(indx + p2) >> 1].h + Dgrb2[(indx - p2) >> 1].h + Dgrb2[(indx + m2) >> 1].h));
                                float gvarv = epssq + (gquinc[0] * Dgrb2[indx >> 1].v +
                                                       gquinc[1] * (Dgrb2[(indx - m1) >> 1].v + Dgrb2[(indx + p1) >> 1].v + Dgrb2[(indx - p1) >> 1].v + Dgrb2[(indx + m1) >> 1].v) +
                                                       gquinc[2] * (Dgrb2[(indx - v2) >> 1].v + Dgrb2[(indx - 2) >> 1].v + Dgrb2[(indx + 2) >> 1].v + Dgrb2[(indx + v2) >> 1].v) +
                                                       gquinc[3] * (Dgrb2[(indx - m2) >> 1].v + Dgrb2[(indx + p2) >> 1].v + Dgrb2[(indx - p2) >> 1].v + Dgrb2[(indx + m2) >> 1].v));
                                //use the results as weights for refined G interpolation
                                Dgrb[0][indx >> 1] = (hcd[indx] * gvarv + vcd[indx] * gvarh) / (gvarv + gvarh);
                                rgbgreen[indx] = cfa[indx] + Dgrb[0][indx >> 1];
                            }
                        }
                }


#ifdef __SSE2__

                for (int rr = 6; rr < rr1 - 6; rr++) {
                    if((FC(rr, 2) & 1) == 0) {
                        for (int cc = 6, indx = (rr) * TS + cc; cc < cc1 - 6; cc += 8, indx += 8) {
                            vfloat tempv = LC2VFU(cfa[indx + 1]);
                            vfloat Dgrbsq1pv = (SQRV(tempv - LC2VFU(cfa[indx + 1 - p1])) + SQRV(tempv - LC2VFU(cfa[indx + 1 + p1])));
                            STVFU(delp[indx >> 1], vabsf(LC2VFU(cfa[indx + p1]) - LC2VFU(cfa[indx - p1])));
                            STVFU(delm[indx >> 1], vabsf(LC2VFU(cfa[indx + m1]) - LC2VFU(cfa[indx - m1])));
                            vfloat Dgrbsq1mv = (SQRV(tempv - LC2VFU(cfa[indx + 1 - m1])) + SQRV(tempv - LC2VFU(cfa[indx + 1 + m1])));
                            STVFU(Dgrbsq1m[indx >> 1], Dgrbsq1mv );
                            STVFU(Dgrbsq1p[indx >> 1], Dgrbsq1pv );
                        }
                    } else {
                        for (int cc = 6, indx = (rr) * TS + cc; cc < cc1 - 6; cc += 8, indx += 8) {
                            vfloat tempv = LC2VFU(cfa[indx]);
                            vfloat Dgrbsq1pv = (SQRV(tempv - LC2VFU(cfa[indx - p1])) + SQRV(tempv - LC2VFU(cfa[indx + p1])));
                            STVFU(delp[indx >> 1], vabsf(LC2VFU(cfa[indx + 1 + p1]) - LC2VFU(cfa[indx + 1 - p1])));
                            STVFU(delm[indx >> 1], vabsf(LC2VFU(cfa[indx + 1 + m1]) - LC2VFU(cfa[indx + 1 - m1])));
                            vfloat Dgrbsq1mv = (SQRV(tempv - LC2VFU(cfa[indx - m1])) + SQRV(tempv - LC2VFU(cfa[indx + m1])));
                            STVFU(Dgrbsq1m[indx >> 1], Dgrbsq1mv );
                            STVFU(Dgrbsq1p[indx >> 1], Dgrbsq1pv );
                        }
                    }
                }

#else

                for (int rr = 6; rr < rr1 - 6; rr++) {
                    if((FC(rr, 2) & 1) == 0) {
                        for (int cc = 6, indx = (rr) * TS + cc; cc < cc1 - 6; cc += 2, indx += 2) {
                            delp[indx >> 1] = fabsf(cfa[indx + p1] - cfa[indx - p1]);
                            delm[indx >> 1] = fabsf(cfa[indx + m1] - cfa[indx - m1]);
                            Dgrbsq1p[indx >> 1] = (SQR(cfa[indx + 1] - cfa[indx + 1 - p1]) + SQR(cfa[indx + 1] - cfa[indx + 1 + p1]));
                            Dgrbsq1m[indx >> 1] = (SQR(cfa[indx + 1] - cfa[indx + 1 - m1]) + SQR(cfa[indx + 1] - cfa[indx + 1 + m1]));
                        }
                    } else {
                        for (int cc = 6, indx = (rr) * TS + cc; cc < cc1 - 6; cc += 2, indx += 2) {
                            Dgrbsq1p[indx >> 1] = (SQR(cfa[indx] - cfa[indx - p1]) + SQR(cfa[indx] - cfa[indx + p1]));
                            Dgrbsq1m[indx >> 1] = (SQR(cfa[indx] - cfa[indx - m1]) + SQR(cfa[indx] - cfa[indx + m1]));
                            delp[indx >> 1] = fabsf(cfa[indx + 1 + p1] - cfa[indx + 1 - p1]);
                            delm[indx >> 1] = fabsf(cfa[indx + 1 + m1] - cfa[indx + 1 - m1]);
                        }
                    }
                }

#endif

                // diagonal interpolation correction

#ifdef __SSE2__
                vfloat gausseven0v = F2V(gausseven[0]);
                vfloat gausseven1v = F2V(gausseven[1]);
#endif

                for (int rr = 8; rr < rr1 - 8; rr++) {
#ifdef __SSE2__

                    for (int indx = rr * TS + 8 + (FC(rr, 2) & 1), indx1 = indx >> 1; indx < rr * TS + cc1 - 8; indx += 8, indx1 += 4) {

                        //diagonal color ratios
                        vfloat cfav = LC2VFU(cfa[indx]);

                        vfloat temp1v = LC2VFU(cfa[indx + m1]);
                        vfloat temp2v = LC2VFU(cfa[indx + m2]);
                        vfloat rbsev = vmul2f(temp1v) / (epsv + cfav + temp2v );
                        rbsev = vself(vmaskf_lt(vabsf(onev - rbsev), arthreshv), cfav * rbsev, temp1v + zd5v * (cfav - temp2v));

                        temp1v = LC2VFU(cfa[indx - m1]);
                        temp2v = LC2VFU(cfa[indx - m2]);
                        vfloat rbnwv = vmul2f(temp1v) / (epsv + cfav + temp2v );
                        rbnwv = vself(vmaskf_lt(vabsf(onev - rbnwv), arthreshv), cfav * rbnwv, temp1v + zd5v * (cfav - temp2v));

                        temp1v = epsv + LVFU(delm[indx1]);
                        vfloat wtsev = temp1v + LVFU(delm[(indx + m1) >> 1]) + LVFU(delm[(indx + m2) >> 1]); //same as for wtu,wtd,wtl,wtr
                        vfloat wtnwv = temp1v + LVFU(delm[(indx - m1) >> 1]) + LVFU(delm[(indx - m2) >> 1]);

                        vfloat rbmv = (wtsev * rbnwv + wtnwv * rbsev) / (wtsev + wtnwv);

                        temp1v = ULIMV(rbmv , LC2VFU(cfa[indx - m1]), LC2VFU(cfa[indx + m1]));
                        vfloat wtv = vmul2f(cfav - rbmv) / (epsv + rbmv + cfav);
                        temp2v = vintpf(wtv, rbmv, temp1v);

                        temp2v = vself(vmaskf_lt(rbmv + rbmv, cfav), temp1v, temp2v);
                        temp2v = vself(vmaskf_lt(rbmv, cfav), temp2v, rbmv);
                        STVFU(rbm[indx1], vself(vmaskf_gt(temp2v, clip_ptv), ULIMV(temp2v , LC2VFU(cfa[indx - m1]), LC2VFU(cfa[indx + m1])), temp2v ));


                        temp1v = LC2VFU(cfa[indx + p1]);
                        temp2v = LC2VFU(cfa[indx + p2]);
                        vfloat rbnev = vmul2f(temp1v) / (epsv + cfav + temp2v );
                        rbnev = vself(vmaskf_lt(vabsf(onev - rbnev), arthreshv), cfav * rbnev, temp1v + zd5v * (cfav - temp2v));

                        temp1v = LC2VFU(cfa[indx - p1]);
                        temp2v = LC2VFU(cfa[indx - p2]);
                        vfloat rbswv = vmul2f(temp1v) / (epsv + cfav + temp2v );
                        rbswv = vself(vmaskf_lt(vabsf(onev - rbswv), arthreshv), cfav * rbswv, temp1v + zd5v * (cfav - temp2v));

                        temp1v = epsv + LVFU(delp[indx1]);
                        vfloat wtnev = temp1v + LVFU(delp[(indx + p1) >> 1]) + LVFU(delp[(indx + p2) >> 1]);
                        vfloat wtswv = temp1v + LVFU(delp[(indx - p1) >> 1]) + LVFU(delp[(indx - p2) >> 1]);

                        vfloat rbpv = (wtnev * rbswv + wtswv * rbnev) / (wtnev + wtswv);

                        temp1v = ULIMV(rbpv , LC2VFU(cfa[indx - p1]), LC2VFU(cfa[indx + p1]));
                        wtv = vmul2f(cfav - rbpv) / (epsv + rbpv + cfav);
                        temp2v = vintpf(wtv, rbpv, temp1v);

                        temp2v = vself(vmaskf_lt(rbpv + rbpv, cfav), temp1v, temp2v);
                        temp2v = vself(vmaskf_lt(rbpv, cfav), temp2v, rbpv);
                        STVFU(rbp[indx1], vself(vmaskf_gt(temp2v, clip_ptv), ULIMV(temp2v , LC2VFU(cfa[indx - p1]), LC2VFU(cfa[indx + p1])), temp2v ));

                        vfloat rbvarmv = epssqv + (gausseven0v * (LVFU(Dgrbsq1m[(indx - v1) >> 1]) + LVFU(Dgrbsq1m[(indx - 1) >> 1]) + LVFU(Dgrbsq1m[(indx + 1) >> 1]) + LVFU(Dgrbsq1m[(indx + v1) >> 1])) +
                                                   gausseven1v * (LVFU(Dgrbsq1m[(indx - v2 - 1) >> 1]) + LVFU(Dgrbsq1m[(indx - v2 + 1) >> 1]) + LVFU(Dgrbsq1m[(indx - 2 - v1) >> 1]) + LVFU(Dgrbsq1m[(indx + 2 - v1) >> 1]) +
                                                           LVFU(Dgrbsq1m[(indx - 2 + v1) >> 1]) + LVFU(Dgrbsq1m[(indx + 2 + v1) >> 1]) + LVFU(Dgrbsq1m[(indx + v2 - 1) >> 1]) + LVFU(Dgrbsq1m[(indx + v2 + 1) >> 1])));
                        STVFU(pmwt[indx1] , rbvarmv / ((epssqv + (gausseven0v * (LVFU(Dgrbsq1p[(indx - v1) >> 1]) + LVFU(Dgrbsq1p[(indx - 1) >> 1]) + LVFU(Dgrbsq1p[(indx + 1) >> 1]) + LVFU(Dgrbsq1p[(indx + v1) >> 1])) +
                                                        gausseven1v * (LVFU(Dgrbsq1p[(indx - v2 - 1) >> 1]) + LVFU(Dgrbsq1p[(indx - v2 + 1) >> 1]) + LVFU(Dgrbsq1p[(indx - 2 - v1) >> 1]) + LVFU(Dgrbsq1p[(indx + 2 - v1) >> 1]) +
                                                                LVFU(Dgrbsq1p[(indx - 2 + v1) >> 1]) + LVFU(Dgrbsq1p[(indx + 2 + v1) >> 1]) + LVFU(Dgrbsq1p[(indx + v2 - 1) >> 1]) + LVFU(Dgrbsq1p[(indx + v2 + 1) >> 1])))) + rbvarmv));

                    }

#else

                    for (int cc = 8 + (FC(rr, 2) & 1), indx = rr * TS + cc, indx1 = indx >> 1; cc < cc1 - 8; cc += 2, indx += 2, indx1++) {

                        //diagonal color ratios
                        float crse = xmul2f(cfa[indx + m1]) / (eps + cfa[indx] + (cfa[indx + m2]));
                        float crnw = xmul2f(cfa[indx - m1]) / (eps + cfa[indx] + (cfa[indx - m2]));
                        float crne = xmul2f(cfa[indx + p1]) / (eps + cfa[indx] + (cfa[indx + p2]));
                        float crsw = xmul2f(cfa[indx - p1]) / (eps + cfa[indx] + (cfa[indx - p2]));
                        //color differences in diagonal directions
                        float rbse, rbnw, rbne, rbsw;

                        //assign B/R at R/B sites
                        if (fabsf(1.0f - crse) < arthresh) {
                            rbse = cfa[indx] * crse;    //use this if more precise diag interp is necessary
                        } else {
                            rbse = (cfa[indx + m1]) + xdiv2f(cfa[indx] - cfa[indx + m2]);
                        }

                        if (fabsf(1.0f - crnw) < arthresh) {
                            rbnw = cfa[indx] * crnw;
                        } else {
                            rbnw = (cfa[indx - m1]) + xdiv2f(cfa[indx] - cfa[indx - m2]);
                        }

                        if (fabsf(1.0f - crne) < arthresh) {
                            rbne = cfa[indx] * crne;
                        } else {
                            rbne = (cfa[indx + p1]) + xdiv2f(cfa[indx] - cfa[indx + p2]);
                        }

                        if (fabsf(1.0f - crsw) < arthresh) {
                            rbsw = cfa[indx] * crsw;
                        } else {
                            rbsw = (cfa[indx - p1]) + xdiv2f(cfa[indx] - cfa[indx - p2]);
                        }

                        float wtse = eps + delm[indx1] + delm[(indx + m1) >> 1] + delm[(indx + m2) >> 1]; //same as for wtu,wtd,wtl,wtr
                        float wtnw = eps + delm[indx1] + delm[(indx - m1) >> 1] + delm[(indx - m2) >> 1];
                        float wtne = eps + delp[indx1] + delp[(indx + p1) >> 1] + delp[(indx + p2) >> 1];
                        float wtsw = eps + delp[indx1] + delp[(indx - p1) >> 1] + delp[(indx - p2) >> 1];


                        rbm[indx1] = (wtse * rbnw + wtnw * rbse) / (wtse + wtnw);
                        rbp[indx1] = (wtne * rbsw + wtsw * rbne) / (wtne + wtsw);

                        //variance of R-B in plus/minus directions
                        float rbvarm = epssq + (gausseven[0] * (Dgrbsq1m[(indx - v1) >> 1] + Dgrbsq1m[(indx - 1) >> 1] + Dgrbsq1m[(indx + 1) >> 1] + Dgrbsq1m[(indx + v1) >> 1]) +
                                                gausseven[1] * (Dgrbsq1m[(indx - v2 - 1) >> 1] + Dgrbsq1m[(indx - v2 + 1) >> 1] + Dgrbsq1m[(indx - 2 - v1) >> 1] + Dgrbsq1m[(indx + 2 - v1) >> 1] +
                                                                Dgrbsq1m[(indx - 2 + v1) >> 1] + Dgrbsq1m[(indx + 2 + v1) >> 1] + Dgrbsq1m[(indx + v2 - 1) >> 1] + Dgrbsq1m[(indx + v2 + 1) >> 1]));
                        pmwt[indx1] = rbvarm / ((epssq + (gausseven[0] * (Dgrbsq1p[(indx - v1) >> 1] + Dgrbsq1p[(indx - 1) >> 1] + Dgrbsq1p[(indx + 1) >> 1] + Dgrbsq1p[(indx + v1) >> 1]) +
                                                          gausseven[1] * (Dgrbsq1p[(indx - v2 - 1) >> 1] + Dgrbsq1p[(indx - v2 + 1) >> 1] + Dgrbsq1p[(indx - 2 - v1) >> 1] + Dgrbsq1p[(indx + 2 - v1) >> 1] +
                                                                  Dgrbsq1p[(indx - 2 + v1) >> 1] + Dgrbsq1p[(indx + 2 + v1) >> 1] + Dgrbsq1p[(indx + v2 - 1) >> 1] + Dgrbsq1p[(indx + v2 + 1) >> 1]))) + rbvarm);

                        //bound the interpolation in regions of high saturation

                        if (rbp[indx1] < cfa[indx]) {
                            if (xmul2f(rbp[indx1]) < cfa[indx]) {
                                rbp[indx1] = ULIM(rbp[indx1] , cfa[indx - p1], cfa[indx + p1]);
                            } else {
                                float pwt = xmul2f(cfa[indx] - rbp[indx1]) / (eps + rbp[indx1] + cfa[indx]);
                                rbp[indx1] = pwt * rbp[indx1] + (1.0f - pwt) * ULIM(rbp[indx1], cfa[indx - p1], cfa[indx + p1]);
                            }
                        }

                        if (rbm[indx1] < cfa[indx]) {
                            if (xmul2f(rbm[indx1]) < cfa[indx]) {
                                rbm[indx1] = ULIM(rbm[indx1] , cfa[indx - m1], cfa[indx + m1]);
                            } else {
                                float mwt = xmul2f(cfa[indx] - rbm[indx1]) / (eps + rbm[indx1] + cfa[indx]);
                                rbm[indx1] = mwt * rbm[indx1] + (1.0f - mwt) * ULIM(rbm[indx1], cfa[indx - m1], cfa[indx + m1]);
                            }
                        }

                        if (rbp[indx1] > clip_pt) {
                            rbp[indx1] = ULIM(rbp[indx1], cfa[indx - p1], cfa[indx + p1]);
                        }

                        if (rbm[indx1] > clip_pt) {
                            rbm[indx1] = ULIM(rbm[indx1], cfa[indx - m1], cfa[indx + m1]);
                        }
                    }

#endif
                }

#ifdef __SSE2__
                vfloat zd25v = F2V(0.25f);
#endif

                for (int rr = 10; rr < rr1 - 10; rr++)
#ifdef __SSE2__
                    for (int indx = rr * TS + 10 + (FC(rr, 2) & 1), indx1 = indx >> 1; indx < rr * TS + cc1 - 10; indx += 8, indx1 += 4) {

                        //first ask if one gets more directional discrimination from nearby B/R sites
                        vfloat pmwtaltv = zd25v * (LVFU(pmwt[(indx - m1) >> 1]) + LVFU(pmwt[(indx + p1) >> 1]) + LVFU(pmwt[(indx - p1) >> 1]) + LVFU(pmwt[(indx + m1) >> 1]));
                        vfloat tempv = LVFU(pmwt[indx1]);
                        tempv = vself(vmaskf_lt(vabsf(zd5v - tempv), vabsf(zd5v - pmwtaltv)), pmwtaltv, tempv);
                        STVFU(pmwt[indx1], tempv);
                        STVFU(rbint[indx1], zd5v * (LC2VFU(cfa[indx]) + vintpf(tempv, LVFU(rbp[indx1]), LVFU(rbm[indx1]))));
                    }

#else

                    for (int cc = 10 + (FC(rr, 2) & 1), indx = rr * TS + cc, indx1 = indx >> 1; cc < cc1 - 10; cc += 2, indx += 2, indx1++) {

                        //first ask if one gets more directional discrimination from nearby B/R sites
                        float pmwtalt = xdivf(pmwt[(indx - m1) >> 1] + pmwt[(indx + p1) >> 1] + pmwt[(indx - p1) >> 1] + pmwt[(indx + m1) >> 1], 2);

                        if (fabsf(0.5 - pmwt[indx1]) < fabsf(0.5 - pmwtalt)) {
                            pmwt[indx1] = pmwtalt;   //a better result was obtained from the neighbours
                        }

                        rbint[indx1] = xdiv2f(cfa[indx] + rbm[indx1] * (1.0f - pmwt[indx1]) + rbp[indx1] * pmwt[indx1]); //this is R+B, interpolated
                    }

#endif

                for (int rr = 12; rr < rr1 - 12; rr++)
#ifdef __SSE2__
                    for (int indx = rr * TS + 12 + (FC(rr, 2) & 1), indx1 = indx >> 1; indx < rr * TS + cc1 - 12; indx += 8, indx1 += 4) {
                        vmask copymask = vmaskf_ge(vabsf(zd5v - LVFU(pmwt[indx1])), vabsf(zd5v - LVFU(hvwt[indx1])));

                        if(_mm_movemask_ps((vfloat)copymask)) { // if for any of the 4 pixels the condition is true, do the math for all 4 pixels and mask the unused out at the end
                            //now interpolate G vertically/horizontally using R+B values
                            //unfortunately, since G interpolation cannot be done diagonally this may lead to color shifts
                            //color ratios for G interpolation
                            vfloat rbintv = LVFU(rbint[indx1]);

                            //interpolated G via adaptive ratios or Hamilton-Adams in each cardinal direction
                            vfloat cruv = vmul2f(LC2VFU(cfa[indx - v1])) / (epsv + rbintv + LVFU(rbint[(indx1 - v1)]));
                            vfloat guv = rbintv * cruv;
                            vfloat gu2v = LC2VFU(cfa[indx - v1]) + zd5v * (rbintv - LVFU(rbint[(indx1 - v1)]));
                            guv = vself(vmaskf_lt(vabsf(onev - cruv), arthreshv), guv, gu2v);

                            vfloat crdv = vmul2f(LC2VFU(cfa[indx + v1])) / (epsv + rbintv + LVFU(rbint[(indx1 + v1)]));
                            vfloat gdv = rbintv * crdv;
                            vfloat gd2v = LC2VFU(cfa[indx + v1]) + zd5v * (rbintv - LVFU(rbint[(indx1 + v1)]));
                            gdv = vself(vmaskf_lt(vabsf(onev - crdv), arthreshv), gdv, gd2v);

                            vfloat Gintvv = (LC2VFU(dirwts0[indx - v1]) * gdv + LC2VFU(dirwts0[indx + v1]) * guv) / (LC2VFU(dirwts0[indx + v1]) + LC2VFU(dirwts0[indx - v1]));
                            vfloat Gint1v = ULIMV(Gintvv , LC2VFU(cfa[indx - v1]), LC2VFU(cfa[indx + v1]));
                            vfloat vwtv = vmul2f(rbintv - Gintvv) / (epsv + Gintvv + rbintv);
                            vfloat Gint2v = vintpf(vwtv, Gintvv, Gint1v);
                            Gint1v = vself(vmaskf_lt(vmul2f(Gintvv), rbintv), Gint1v, Gint2v);
                            Gintvv = vself(vmaskf_lt(Gintvv, rbintv), Gint1v, Gintvv);
                            Gintvv = vself(vmaskf_gt(Gintvv, clip_ptv), ULIMV(Gintvv, LC2VFU(cfa[indx - v1]), LC2VFU(cfa[indx + v1])), Gintvv);

                            vfloat crlv = vmul2f(LC2VFU(cfa[indx - 1])) / (epsv + rbintv + LVFU(rbint[(indx1 - 1)]));
                            vfloat glv = rbintv * crlv;
                            vfloat gl2v = LC2VFU(cfa[indx - 1]) + zd5v * (rbintv - LVFU(rbint[(indx1 - 1)]));
                            glv = vself(vmaskf_lt(vabsf(onev - crlv), arthreshv), glv, gl2v);

                            vfloat crrv = vmul2f(LC2VFU(cfa[indx + 1])) / (epsv + rbintv + LVFU(rbint[(indx1 + 1)]));
                            vfloat grv = rbintv * crrv;
                            vfloat gr2v = LC2VFU(cfa[indx + 1]) + zd5v * (rbintv - LVFU(rbint[(indx1 + 1)]));
                            grv = vself(vmaskf_lt(vabsf(onev - crrv), arthreshv), grv, gr2v);

                            vfloat Ginthv = (LC2VFU(dirwts1[indx - 1]) * grv + LC2VFU(dirwts1[indx + 1]) * glv) / (LC2VFU(dirwts1[indx - 1]) + LC2VFU(dirwts1[indx + 1]));
                            vfloat Gint1h = ULIMV(Ginthv , LC2VFU(cfa[indx - 1]), LC2VFU(cfa[indx + 1]));
                            vfloat hwtv = vmul2f(rbintv - Ginthv) / (epsv + Ginthv + rbintv);
                            vfloat Gint2h = vintpf(hwtv, Ginthv, Gint1h);
                            Gint1h = vself(vmaskf_lt(vmul2f(Ginthv), rbintv), Gint1h, Gint2h);
                            Ginthv = vself(vmaskf_lt(Ginthv, rbintv), Gint1h, Ginthv);
                            Ginthv = vself(vmaskf_gt(Ginthv, clip_ptv), ULIMV(Ginthv, LC2VFU(cfa[indx - 1]), LC2VFU(cfa[indx + 1])), Ginthv);

                            vfloat greenv = vself(copymask, vintpf(LVFU(hvwt[indx1]), Gintvv, Ginthv), LC2VFU(rgbgreen[indx]));
                            STC2VFU(rgbgreen[indx], greenv);

                            STVFU(Dgrb[0][indx1], vself(copymask, greenv - LC2VFU(cfa[indx]), LVFU(Dgrb[0][indx1])));
                        }
                    }

#else

                    for (int cc = 12 + (FC(rr, 2) & 1), indx = rr * TS + cc, indx1 = indx >> 1; cc < cc1 - 12; cc += 2, indx += 2, indx1++) {

                        if (fabsf(0.5 - pmwt[indx >> 1]) < fabsf(0.5 - hvwt[indx >> 1]) ) {
                            continue;
                        }

                        //now interpolate G vertically/horizontally using R+B values
                        //unfortunately, since G interpolation cannot be done diagonally this may lead to color shifts

                        //color ratios for G interpolation
                        float cru = cfa[indx - v1] * 2.0 / (eps + rbint[indx1] + rbint[(indx1 - v1)]);
                        float crd = cfa[indx + v1] * 2.0 / (eps + rbint[indx1] + rbint[(indx1 + v1)]);
                        float crl = cfa[indx - 1] * 2.0 / (eps + rbint[indx1] + rbint[(indx1 - 1)]);
                        float crr = cfa[indx + 1] * 2.0 / (eps + rbint[indx1] + rbint[(indx1 + 1)]);

                        //interpolation of G in four directions
                        float gu, gd, gl, gr;

                        //interpolated G via adaptive ratios or Hamilton-Adams in each cardinal direction
                        if (fabsf(1.f - cru) < arthresh) {
                            gu = rbint[indx1] * cru;
                        } else {
                            gu = cfa[indx - v1] + xdiv2f(rbint[indx1] - rbint[(indx1 - v1)]);
                        }

                        if (fabsf(1.f - crd) < arthresh) {
                            gd = rbint[indx1] * crd;
                        } else {
                            gd = cfa[indx + v1] + xdiv2f(rbint[indx1] - rbint[(indx1 + v1)]);
                        }

                        if (fabsf(1.f - crl) < arthresh) {
                            gl = rbint[indx1] * crl;
                        } else {
                            gl = cfa[indx - 1] + xdiv2f(rbint[indx1] - rbint[(indx1 - 1)]);
                        }

                        if (fabsf(1.f - crr) < arthresh) {
                            gr = rbint[indx1] * crr;
                        } else {
                            gr = cfa[indx + 1] + xdiv2f(rbint[indx1] - rbint[(indx1 + 1)]);
                        }

                        //interpolated G via adaptive weights of cardinal evaluations
                        float Gintv = (dirwts0[indx - v1] * gd + dirwts0[indx + v1] * gu) / (dirwts0[indx + v1] + dirwts0[indx - v1]);
                        float Ginth = (dirwts1[indx - 1] * gr + dirwts1[indx + 1] * gl) / (dirwts1[indx - 1] + dirwts1[indx + 1]);

                        //bound the interpolation in regions of high saturation
                        if (Gintv < rbint[indx1]) {
                            if (2 * Gintv < rbint[indx1]) {
                                Gintv = ULIM(Gintv , cfa[indx - v1], cfa[indx + v1]);
                            } else {
                                float vwt = 2.0 * (rbint[indx1] - Gintv) / (eps + Gintv + rbint[indx1]);
                                Gintv = vwt * Gintv + (1.f - vwt) * ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]);
                            }
                        }

                        if (Ginth < rbint[indx1]) {
                            if (2 * Ginth < rbint[indx1]) {
                                Ginth = ULIM(Ginth , cfa[indx - 1], cfa[indx + 1]);
                            } else {
                                float hwt = 2.0 * (rbint[indx1] - Ginth) / (eps + Ginth + rbint[indx1]);
                                Ginth = hwt * Ginth + (1.f - hwt) * ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]);
                            }
                        }

                        if (Ginth > clip_pt) {
                            Ginth = ULIM(Ginth, cfa[indx - 1], cfa[indx + 1]);
                        }

                        if (Gintv > clip_pt) {
                            Gintv = ULIM(Gintv, cfa[indx - v1], cfa[indx + v1]);
                        }

                        rgbgreen[indx] = Ginth * (1.f - hvwt[indx1]) + Gintv * hvwt[indx1];
                        Dgrb[0][indx >> 1] = rgbgreen[indx] - cfa[indx];
                    }

#endif

                //end of diagonal interpolation correction

                //fancy chrominance interpolation
                //(ey,ex) is location of R site
                for (int rr = 13 - ey; rr < rr1 - 12; rr += 2)
                    for (int indx1 = (rr * TS + 13 - ex) >> 1; indx1 < (rr * TS + cc1 - 12) >> 1; indx1++) { //B coset
                        Dgrb[1][indx1] = Dgrb[0][indx1]; //split out G-B from G-R
                        Dgrb[0][indx1] = 0;
                    }

#ifdef __SSE2__
                vfloat oned325v = F2V( 1.325f );
                vfloat zd175v = F2V( 0.175f );
                vfloat zd075v = F2V( 0.075f );
#endif

                for (int rr = 14; rr < rr1 - 14; rr++)
#ifdef __SSE2__
                    for (int cc = 14 + (FC(rr, 2) & 1), indx = rr * TS + cc, c = 1 - FC(rr, cc) / 2; cc < cc1 - 14; cc += 8, indx += 8) {
                        vfloat tempv = epsv + vabsf(LVFU(Dgrb[c][(indx - m1) >> 1]) - LVFU(Dgrb[c][(indx + m1) >> 1]));
                        vfloat temp2v = epsv + vabsf(LVFU(Dgrb[c][(indx + p1) >> 1]) - LVFU(Dgrb[c][(indx - p1) >> 1]));
                        vfloat wtnwv = onev / (tempv + vabsf(LVFU(Dgrb[c][(indx - m1) >> 1]) - LVFU(Dgrb[c][(indx - m3) >> 1])) + vabsf(LVFU(Dgrb[c][(indx + m1) >> 1]) - LVFU(Dgrb[c][(indx - m3) >> 1])));
                        vfloat wtnev = onev / (temp2v + vabsf(LVFU(Dgrb[c][(indx + p1) >> 1]) - LVFU(Dgrb[c][(indx + p3) >> 1])) + vabsf(LVFU(Dgrb[c][(indx - p1) >> 1]) - LVFU(Dgrb[c][(indx + p3) >> 1])));
                        vfloat wtswv = onev / (temp2v + vabsf(LVFU(Dgrb[c][(indx - p1) >> 1]) - LVFU(Dgrb[c][(indx + m3) >> 1])) + vabsf(LVFU(Dgrb[c][(indx + p1) >> 1]) - LVFU(Dgrb[c][(indx - p3) >> 1])));
                        vfloat wtsev = onev / (tempv + vabsf(LVFU(Dgrb[c][(indx + m1) >> 1]) - LVFU(Dgrb[c][(indx - p3) >> 1])) + vabsf(LVFU(Dgrb[c][(indx - m1) >> 1]) - LVFU(Dgrb[c][(indx + m3) >> 1])));

                        STVFU(Dgrb[c][indx >> 1], (wtnwv * (oned325v * LVFU(Dgrb[c][(indx - m1) >> 1]) - zd175v * LVFU(Dgrb[c][(indx - m3) >> 1]) - zd075v * (LVFU(Dgrb[c][(indx - m1 - 2) >> 1]) + LVFU(Dgrb[c][(indx - m1 - v2) >> 1])) ) +
                                                   wtnev * (oned325v * LVFU(Dgrb[c][(indx + p1) >> 1]) - zd175v * LVFU(Dgrb[c][(indx + p3) >> 1]) - zd075v * (LVFU(Dgrb[c][(indx + p1 + 2) >> 1]) + LVFU(Dgrb[c][(indx + p1 + v2) >> 1])) ) +
                                                   wtswv * (oned325v * LVFU(Dgrb[c][(indx - p1) >> 1]) - zd175v * LVFU(Dgrb[c][(indx - p3) >> 1]) - zd075v * (LVFU(Dgrb[c][(indx - p1 - 2) >> 1]) + LVFU(Dgrb[c][(indx - p1 - v2) >> 1])) ) +
                                                   wtsev * (oned325v * LVFU(Dgrb[c][(indx + m1) >> 1]) - zd175v * LVFU(Dgrb[c][(indx + m3) >> 1]) - zd075v * (LVFU(Dgrb[c][(indx + m1 + 2) >> 1]) + LVFU(Dgrb[c][(indx + m1 + v2) >> 1])) )) / (wtnwv + wtnev + wtswv + wtsev));
                    }

#else

                    for (int cc = 14 + (FC(rr, 2) & 1), indx = rr * TS + cc, c = 1 - FC(rr, cc) / 2; cc < cc1 - 14; cc += 2, indx += 2) {
                        float wtnw = 1.0f / (eps + fabsf(Dgrb[c][(indx - m1) >> 1] - Dgrb[c][(indx + m1) >> 1]) + fabsf(Dgrb[c][(indx - m1) >> 1] - Dgrb[c][(indx - m3) >> 1]) + fabsf(Dgrb[c][(indx + m1) >> 1] - Dgrb[c][(indx - m3) >> 1]));
                        float wtne = 1.0f / (eps + fabsf(Dgrb[c][(indx + p1) >> 1] - Dgrb[c][(indx - p1) >> 1]) + fabsf(Dgrb[c][(indx + p1) >> 1] - Dgrb[c][(indx + p3) >> 1]) + fabsf(Dgrb[c][(indx - p1) >> 1] - Dgrb[c][(indx + p3) >> 1]));
                        float wtsw = 1.0f / (eps + fabsf(Dgrb[c][(indx - p1) >> 1] - Dgrb[c][(indx + p1) >> 1]) + fabsf(Dgrb[c][(indx - p1) >> 1] - Dgrb[c][(indx + m3) >> 1]) + fabsf(Dgrb[c][(indx + p1) >> 1] - Dgrb[c][(indx - p3) >> 1]));
                        float wtse = 1.0f / (eps + fabsf(Dgrb[c][(indx + m1) >> 1] - Dgrb[c][(indx - m1) >> 1]) + fabsf(Dgrb[c][(indx + m1) >> 1] - Dgrb[c][(indx - p3) >> 1]) + fabsf(Dgrb[c][(indx - m1) >> 1] - Dgrb[c][(indx + m3) >> 1]));

                        Dgrb[c][indx >> 1] = (wtnw * (1.325f * Dgrb[c][(indx - m1) >> 1] - 0.175f * Dgrb[c][(indx - m3) >> 1] - 0.075f * Dgrb[c][(indx - m1 - 2) >> 1] - 0.075f * Dgrb[c][(indx - m1 - v2) >> 1] ) +
                                              wtne * (1.325f * Dgrb[c][(indx + p1) >> 1] - 0.175f * Dgrb[c][(indx + p3) >> 1] - 0.075f * Dgrb[c][(indx + p1 + 2) >> 1] - 0.075f * Dgrb[c][(indx + p1 + v2) >> 1] ) +
                                              wtsw * (1.325f * Dgrb[c][(indx - p1) >> 1] - 0.175f * Dgrb[c][(indx - p3) >> 1] - 0.075f * Dgrb[c][(indx - p1 - 2) >> 1] - 0.075f * Dgrb[c][(indx - p1 - v2) >> 1] ) +
                                              wtse * (1.325f * Dgrb[c][(indx + m1) >> 1] - 0.175f * Dgrb[c][(indx + m3) >> 1] - 0.075f * Dgrb[c][(indx + m1 + 2) >> 1] - 0.075f * Dgrb[c][(indx + m1 + v2) >> 1] )) / (wtnw + wtne + wtsw + wtse);
                    }

#endif
                //tile vars
                //counters for pixel location in the image
                int row, col;
                //counters for pixel location within the tile
                int cc;
                //pointer counters within the tile
                int indx;

                // end of tile initialization

#ifdef __SSE2__
                int offset;
                vfloat twov = F2V(2.f);
                vmask selmask;

                if((FC(16, 2) & 1) == 1) {
                    selmask = _mm_set_epi32(0xffffffff, 0, 0xffffffff, 0);
                    offset = 1;
                } else {
                    selmask = _mm_set_epi32(0, 0xffffffff, 0, 0xffffffff);
                    offset = 0;
                }

#endif

                for (int rr = 16; rr < rr1 - 16; rr++) {
#ifdef __SSE2__
                    offset = 1 - offset;
                    selmask = vnotm(selmask);

                    for (cc = 16, indx = rr * TS + cc, row = rr + top; cc < cc1 - 18 - (cc1 & 1); cc += 4, indx += 4) {
                        col = cc + left;
                        vfloat greenv = LVF(rgbgreen[indx]);
                        vfloat temp00v = vdup(LVF(hvwt[(indx - v1) >> 1]));
                        vfloat temp01v = vdup(LVF(hvwt[(indx + v1) >> 1]));
                        vfloat tempv =  onev / (temp00v + twov - vdup(LVFU(hvwt[(indx + 1 + offset) >> 1])) - vdup(LVFU(hvwt[(indx - 1 + offset) >> 1])) + temp01v);

                        vfloat redv1  = greenv - (temp00v * vdup(LVF(Dgrb[0][(indx - v1) >> 1])) + (onev - vdup(LVFU(hvwt[(indx + 1 + offset) >> 1]))) * vdup(LVFU(Dgrb[0][(indx + 1 + offset) >> 1])) + (onev - vdup(LVFU(hvwt[(indx - 1 + offset) >> 1]))) * vdup(LVFU(Dgrb[0][(indx - 1 + offset) >> 1])) + temp01v * vdup(LVF(Dgrb[0][(indx + v1) >> 1]))) * tempv;
                        vfloat bluev1 = greenv - (temp00v * vdup(LVF(Dgrb[1][(indx - v1) >> 1])) + (onev - vdup(LVFU(hvwt[(indx + 1 + offset) >> 1]))) * vdup(LVFU(Dgrb[1][(indx + 1 + offset) >> 1])) + (onev - vdup(LVFU(hvwt[(indx - 1 + offset) >> 1]))) * vdup(LVFU(Dgrb[1][(indx - 1 + offset) >> 1])) + temp01v * vdup(LVF(Dgrb[1][(indx + v1) >> 1]))) * tempv;
                        vfloat redv2  = greenv - vdup(LVF(Dgrb[0][indx >> 1]));
                        vfloat bluev2 = greenv - vdup(LVF(Dgrb[1][indx >> 1]));
                        STVFU(red[row][col], c65535v * vself(selmask, redv1, redv2));
                        STVFU(blue[row][col], c65535v * vself(selmask, bluev1, bluev2));
                    }

                    if(offset == 0) {
                        for (indx = rr * TS + cc; cc < cc1 - 16 - (cc1 & 1); cc += 2, indx++) {
                            col = cc + left;
                            float temp =  1.0f / (hvwt[(indx - v1) >> 1] + 2.0f - hvwt[(indx + 1) >> 1] - hvwt[(indx - 1) >> 1] + hvwt[(indx + v1) >> 1]);
                            red[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[0][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[0][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[0][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[0][(indx + v1) >> 1]) *
                                                        temp);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[1][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[1][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[1][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[1][(indx + v1) >> 1]) *
                                                         temp);

                            indx++;
                            col++;
                            red[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[0][indx >> 1]);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[1][indx >> 1]);
                        }

                        if(cc1 & 1) { // width of tile is odd
                            col = cc + left;
                            float temp =  1.0f / (hvwt[(indx - v1) >> 1] + 2.0f - hvwt[(indx + 1) >> 1] - hvwt[(indx - 1) >> 1] + hvwt[(indx + v1) >> 1]);
                            red[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[0][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[0][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[0][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[0][(indx + v1) >> 1]) *
                                                        temp);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[1][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[1][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[1][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[1][(indx + v1) >> 1]) *
                                                         temp);
                        }
                    } else {
                        for (indx = rr * TS + cc; cc < cc1 - 16 - (cc1 & 1); cc += 2, indx++) {
                            col = cc + left;
                            red[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[0][indx >> 1]);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[1][indx >> 1]);

                            indx++;
                            col++;
                            float temp =  1.0f / (hvwt[(indx - v1) >> 1] + 2.0f - hvwt[(indx + 1) >> 1] - hvwt[(indx - 1) >> 1] + hvwt[(indx + v1) >> 1]);
                            red[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[0][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[0][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[0][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[0][(indx + v1) >> 1]) *
                                                        temp);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[1][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[1][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[1][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[1][(indx + v1) >> 1]) *
                                                         temp);
                        }

                        if(cc1 & 1) { // width of tile is odd
                            col = cc + left;
                            red[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[0][indx >> 1]);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[1][indx >> 1]);
                        }
                    }

#else

                    if((FC(rr, 2) & 1) == 1) {
                        for (cc = 16, indx = rr * TS + cc, row = rr + top; cc < cc1 - 16 - (cc1 & 1); cc += 2, indx++) {
                            col = cc + left;
                            float temp =  1.0f / (hvwt[(indx - v1) >> 1] + 2.0f - hvwt[(indx + 1) >> 1] - hvwt[(indx - 1) >> 1] + hvwt[(indx + v1) >> 1]);
                            red[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[0][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[0][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[0][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[0][(indx + v1) >> 1]) *
                                                        temp);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[1][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[1][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[1][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[1][(indx + v1) >> 1]) *
                                                         temp);

                            indx++;
                            col++;
                            red[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[0][indx >> 1]);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[1][indx >> 1]);
                        }

                        if(cc1 & 1) { // width of tile is odd
                            col = cc + left;
                            float temp =  1.0f / (hvwt[(indx - v1) >> 1] + 2.0f - hvwt[(indx + 1) >> 1] - hvwt[(indx - 1) >> 1] + hvwt[(indx + v1) >> 1]);
                            red[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[0][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[0][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[0][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[0][(indx + v1) >> 1]) *
                                                        temp);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[1][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[1][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[1][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[1][(indx + v1) >> 1]) *
                                                         temp);
                        }
                    } else {
                        for (cc = 16, indx = rr * TS + cc, row = rr + top; cc < cc1 - 16 - (cc1 & 1); cc += 2, indx++) {
                            col = cc + left;
                            red[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[0][indx >> 1]);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[1][indx >> 1]);

                            indx++;
                            col++;
                            float temp =  1.0f / (hvwt[(indx - v1) >> 1] + 2.0f - hvwt[(indx + 1) >> 1] - hvwt[(indx - 1) >> 1] + hvwt[(indx + v1) >> 1]);
                            red[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[0][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[0][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[0][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[0][(indx + v1) >> 1]) *
                                                        temp);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - ((hvwt[(indx - v1) >> 1]) * Dgrb[1][(indx - v1) >> 1] + (1.0f - hvwt[(indx + 1) >> 1]) * Dgrb[1][(indx + 1) >> 1] + (1.0f - hvwt[(indx - 1) >> 1]) * Dgrb[1][(indx - 1) >> 1] + (hvwt[(indx + v1) >> 1]) * Dgrb[1][(indx + v1) >> 1]) *
                                                         temp);
                        }

                        if(cc1 & 1) { // width of tile is odd
                            col = cc + left;
                            red[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[0][indx >> 1]);
                            blue[row][col] = 65535.0f * (rgbgreen[indx] - Dgrb[1][indx >> 1]);
                        }
                    }

#endif
                }

                // copy smoothed results back to image matrix
                for (int rr = 16; rr < rr1 - 16; rr++) {
                    int row = rr + top;
                    int cc = 16;
#ifdef __SSE2__

                    for (; cc < cc1 - 19; cc += 4) {
                        STVFU(green[row][cc + left], LVF(rgbgreen[rr * TS + cc]) * c65535v);
                    }

#endif

                    for (; cc < cc1 - 16; cc++) {
                        green[row][cc + left] = 65535.0f * rgbgreen[rr * TS + cc];
                    }
                }

                //end of main loop

                if(plistener) {
                    progresscounter++;

                    if(progresscounter % 32 == 0) {
#ifdef _OPENMP
                        #pragma omp critical (amazeprogress)
#endif
                        {
                            progress += (double)32 * ((TS - 32) * (TS - 32)) / (height * width);
                            progress = progress > 1.0 ? 1.0 : progress;
                            plistener->setProgress(progress);
                        }
                    }
                }
            }

        // clean up
        free(buffer);
    }

    if(plistener) {
        plistener->setProgress(1.0);
    }

    // done

#undef TS

}
}
