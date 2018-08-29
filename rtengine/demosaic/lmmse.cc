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
#include <cmath>
#include <cassert>

#include "../rawimagesource.h"
#include "../rt_math.h"
#include "../color.h"
#include "../../rtgui/multilangmgr.h"
#include "../procparams.h"
#include "../sleef.c"
#include "../opthelper.h"
#include "../median.h"
#define BENCHMARK
#include "../StopWatch.h"

using namespace std;

namespace rtengine
{

#define CLIREF(x) LIM(x,-200000.0f,200000.0f) // avoid overflow : do not act directly on image[] or pix[]
#define x1125(a) (a + xdivf(a, 3))
#define x0875(a) (a - xdivf(a, 3))
#define x0250(a) xdivf(a, 2)
#define x00625(a) xdivf(a, 4)
#define x0125(a) xdivf(a, 3)

// LSMME demosaicing algorithm
// L. Zhang and X. Wu,
// Color demozaicing via directional Linear Minimum Mean Square-error Estimation,
// IEEE Trans. on Image Processing, vol. 14, pp. 2167-2178,
// Dec. 2005.
// Adapted to RawTherapee by Jacques Desmis 3/2013
// Improved speed and reduced memory consumption by Ingo Weyrich 2/2015
//TODO Tiles to reduce memory consumption
void RawImageSource::lmmse_demosaic(int winw, int winh, array2D<float> &rawData, array2D<float> &red, array2D<float> &green, array2D<float> &blue, int iterations)
{
    const int width = winw, height = winh;
    const int ba = 10;
    const int rr1 = height + 2 * ba;
    const int cc1 = width + 2 * ba;
    const int w1 = cc1;
    const int w2 = 2 * w1;
    const int w3 = 3 * w1;
    const int w4 = 4 * w1;
    float h0, h1, h2, h3, h4, hs;
    h0 = 1.0f;
    h1 = exp( -1.0f / 8.0f);
    h2 = exp( -4.0f / 8.0f);
    h3 = exp( -9.0f / 8.0f);
    h4 = exp(-16.0f / 8.0f);
    hs = h0 + 2.0f * (h1 + h2 + h3 + h4);
    h0 /= hs;
    h1 /= hs;
    h2 /= hs;
    h3 /= hs;
    h4 /= hs;
    int passref = 0;
    int iter = 0;

    if(iterations <= 4) {
        iter = iterations - 1;
        passref = 0;
    } else if (iterations <= 6) {
        iter = 3;
        passref = iterations - 4;
    } else if (iterations <= 8) {
        iter = 3;
        passref = iterations - 6;
    }

    bool applyGamma = true;

    if(iterations == 0) {
        applyGamma = false;
        iter = 0;
    } else {
        applyGamma = true;
    }

    float *rix[5];
    float *qix[5];
    float *buffer = (float *)calloc(rr1 * cc1 * 5 * sizeof(float), 1);

    if(buffer == nullptr) { // allocation of big block of memory failed, try to get 5 smaller ones
        printf("lmmse_interpolate_omp: allocation of big memory block failed, try to get 5 smaller ones now...\n");
        bool allocationFailed = false;

        for(int i = 0; i < 5; i++) {
            qix[i] = (float *)calloc(rr1 * cc1 * sizeof(float), 1);

            if(!qix[i]) { // allocation of at least one small block failed
                allocationFailed = true;
            }
        }

        if(allocationFailed) { // fall back to igv_interpolate
            printf("lmmse_interpolate_omp: allocation of 5 small memory blocks failed, falling back to igv_interpolate...\n");

            for(int i = 0; i < 5; i++) { // free the already allocated buffers
                if(qix[i]) {
                    free(qix[i]);
                }
            }

            igv_demosaic(winw, winh);
            return;
        }
    } else {
        qix[0] = buffer;

        for(int i = 1; i < 5; i++) {
            qix[i] = qix[i - 1] + rr1 * cc1;
        }
    }

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::LMMSE)));
        plistener->setProgress (0.0);
    }


    LUTf *gamtab;

    if(applyGamma) {
        gamtab = &(Color::gammatab_24_17a);
    } else {
        gamtab = new LUTf(65536, LUT_CLIP_ABOVE | LUT_CLIP_BELOW);
        gamtab->makeIdentity(65535.f);
    }


#ifdef _OPENMP
    #pragma omp parallel private(rix)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rrr = ba; rrr < rr1 - ba; rrr++) {
            for (int ccc = ba, row = rrr - ba; ccc < cc1 - ba; ccc++) {
                int col = ccc - ba;
                float *rix = qix[4] + rrr * cc1 + ccc;
                rix[0] = (*gamtab)[rawData[row][col]];
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.1);
            }
        }

        // G-R(B)
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int rr = 2; rr < rr1 - 2; rr++) {
            // G-R(B) at R(B) location
            for (int cc = 2 + (FC(rr, 2) & 1); cc < cc1 - 2; cc += 2) {
                rix[4] = qix[4] + rr * cc1 + cc;
                float v0 = x00625(rix[4][-w1 - 1] + rix[4][-w1 + 1] + rix[4][w1 - 1] + rix[4][w1 + 1]) + x0250(rix[4][0]);
                // horizontal
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[0][0] = - x0250(rix[4][ -2] + rix[4][ 2]) + xdiv2f(rix[4][ -1] + rix[4][0] + rix[4][ 1]);
                float Y = v0 + xdiv2f(rix[0][0]);

                if (rix[4][0] > 1.75f * Y) {
                    rix[0][0] = median(rix[0][0], rix[4][ -1], rix[4][ 1]);
                } else {
                    rix[0][0] = LIM(rix[0][0], 0.0f, 1.0f);
                }

                rix[0][0] -= rix[4][0];
                // vertical
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[1][0] = -x0250(rix[4][-w2] + rix[4][w2]) + xdiv2f(rix[4][-w1] + rix[4][0] + rix[4][w1]);
                Y = v0 + xdiv2f(rix[1][0]);

                if (rix[4][0] > 1.75f * Y) {
                    rix[1][0] = median(rix[1][0], rix[4][-w1], rix[4][w1]);
                } else {
                    rix[1][0] = LIM(rix[1][0], 0.0f, 1.0f);
                }

                rix[1][0] -= rix[4][0];
            }

            // G-R(B) at G location
            for (int ccc = 2 + (FC(rr, 3) & 1); ccc < cc1 - 2; ccc += 2) {
                rix[0] = qix[0] + rr * cc1 + ccc;
                rix[1] = qix[1] + rr * cc1 + ccc;
                rix[4] = qix[4] + rr * cc1 + ccc;
                rix[0][0] = x0250(rix[4][ -2] + rix[4][ 2]) - xdiv2f(rix[4][ -1] + rix[4][0] + rix[4][ 1]);
                rix[1][0] = x0250(rix[4][-w2] + rix[4][w2]) - xdiv2f(rix[4][-w1] + rix[4][0] + rix[4][w1]);
                rix[0][0] = LIM(rix[0][0], -1.0f, 0.0f) + rix[4][0];
                rix[1][0] = LIM(rix[1][0], -1.0f, 0.0f) + rix[4][0];
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.2);
            }
        }


        // apply low pass filter on differential colors
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 4; rr < rr1 - 4; rr++)
            for (int cc = 4; cc < cc1 - 4; cc++) {
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[2] = qix[2] + rr * cc1 + cc;
                rix[2][0] = h0 * rix[0][0] + h1 * (rix[0][ -1] + rix[0][ 1]) + h2 * (rix[0][ -2] + rix[0][ 2]) + h3 * (rix[0][ -3] + rix[0][ 3]) + h4 * (rix[0][ -4] + rix[0][ 4]);
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[3] = qix[3] + rr * cc1 + cc;
                rix[3][0] = h0 * rix[1][0] + h1 * (rix[1][-w1] + rix[1][w1]) + h2 * (rix[1][-w2] + rix[1][w2]) + h3 * (rix[1][-w3] + rix[1][w3]) + h4 * (rix[1][-w4] + rix[1][w4]);
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.3);
            }
        }

        // interpolate G-R(B) at R(B)
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 4; rr < rr1 - 4; rr++) {
            int cc = 4 + (FC(rr, 4) & 1);
#ifdef __SSE2__
            __m128 p1v, p2v, p3v, p4v, p5v, p6v, p7v, p8v, p9v, muv, vxv, vnv, xhv, vhv, xvv, vvv;
            __m128 epsv = _mm_set1_ps(1e-7);
            __m128 ninev = _mm_set1_ps(9.f);

            for (; cc < cc1 - 10; cc += 8) {
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[2] = qix[2] + rr * cc1 + cc;
                rix[3] = qix[3] + rr * cc1 + cc;
                rix[4] = qix[4] + rr * cc1 + cc;
                // horizontal
                p1v = LC2VFU(rix[2][-4]);
                p2v = LC2VFU(rix[2][-3]);
                p3v = LC2VFU(rix[2][-2]);
                p4v = LC2VFU(rix[2][-1]);
                p5v = LC2VFU(rix[2][ 0]);
                p6v = LC2VFU(rix[2][ 1]);
                p7v = LC2VFU(rix[2][ 2]);
                p8v = LC2VFU(rix[2][ 3]);
                p9v = LC2VFU(rix[2][ 4]);
                muv = (p1v + p2v + p3v + p4v + p5v + p6v + p7v + p8v + p9v) / ninev;
                vxv = epsv + SQRV(p1v - muv) + SQRV(p2v - muv) + SQRV(p3v - muv) + SQRV(p4v - muv) + SQRV(p5v - muv) + SQRV(p6v - muv) + SQRV(p7v - muv) + SQRV(p8v - muv) + SQRV(p9v - muv);
                p1v -= LC2VFU(rix[0][-4]);
                p2v -= LC2VFU(rix[0][-3]);
                p3v -= LC2VFU(rix[0][-2]);
                p4v -= LC2VFU(rix[0][-1]);
                p5v -= LC2VFU(rix[0][ 0]);
                p6v -= LC2VFU(rix[0][ 1]);
                p7v -= LC2VFU(rix[0][ 2]);
                p8v -= LC2VFU(rix[0][ 3]);
                p9v -= LC2VFU(rix[0][ 4]);
                vnv = epsv + SQRV(p1v) + SQRV(p2v) + SQRV(p3v) + SQRV(p4v) + SQRV(p5v) + SQRV(p6v) + SQRV(p7v) + SQRV(p8v) + SQRV(p9v);
                xhv = (LC2VFU(rix[0][0]) * vxv + LC2VFU(rix[2][0]) * vnv) / (vxv + vnv);
                vhv = vxv * vnv / (vxv + vnv);

                // vertical
                p1v = LC2VFU(rix[3][-w4]);
                p2v = LC2VFU(rix[3][-w3]);
                p3v = LC2VFU(rix[3][-w2]);
                p4v = LC2VFU(rix[3][-w1]);
                p5v = LC2VFU(rix[3][  0]);
                p6v = LC2VFU(rix[3][ w1]);
                p7v = LC2VFU(rix[3][ w2]);
                p8v = LC2VFU(rix[3][ w3]);
                p9v = LC2VFU(rix[3][ w4]);
                muv = (p1v + p2v + p3v + p4v + p5v + p6v + p7v + p8v + p9v) / ninev;
                vxv = epsv + SQRV(p1v - muv) + SQRV(p2v - muv) + SQRV(p3v - muv) + SQRV(p4v - muv) + SQRV(p5v - muv) + SQRV(p6v - muv) + SQRV(p7v - muv) + SQRV(p8v - muv) + SQRV(p9v - muv);
                p1v -= LC2VFU(rix[1][-w4]);
                p2v -= LC2VFU(rix[1][-w3]);
                p3v -= LC2VFU(rix[1][-w2]);
                p4v -= LC2VFU(rix[1][-w1]);
                p5v -= LC2VFU(rix[1][  0]);
                p6v -= LC2VFU(rix[1][ w1]);
                p7v -= LC2VFU(rix[1][ w2]);
                p8v -= LC2VFU(rix[1][ w3]);
                p9v -= LC2VFU(rix[1][ w4]);
                vnv = epsv + SQRV(p1v) + SQRV(p2v) + SQRV(p3v) + SQRV(p4v) + SQRV(p5v) + SQRV(p6v) + SQRV(p7v) + SQRV(p8v) + SQRV(p9v);
                xvv = (LC2VFU(rix[1][0]) * vxv + LC2VFU(rix[3][0]) * vnv) / (vxv + vnv);
                vvv = vxv * vnv / (vxv + vnv);
                // interpolated G-R(B)
                muv = (xhv * vvv + xvv * vhv) / (vhv + vvv);
                STC2VFU(rix[4][0], muv);
            }

#endif

            for (; cc < cc1 - 4; cc += 2) {
                rix[0] = qix[0] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[2] = qix[2] + rr * cc1 + cc;
                rix[3] = qix[3] + rr * cc1 + cc;
                rix[4] = qix[4] + rr * cc1 + cc;
                // horizontal
                float p1 = rix[2][-4];
                float p2 = rix[2][-3];
                float p3 = rix[2][-2];
                float p4 = rix[2][-1];
                float p5 = rix[2][ 0];
                float p6 = rix[2][ 1];
                float p7 = rix[2][ 2];
                float p8 = rix[2][ 3];
                float p9 = rix[2][ 4];
                float mu = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) / 9.f;
                float vx = 1e-7 + SQR(p1 - mu) + SQR(p2 - mu) + SQR(p3 - mu) + SQR(p4 - mu) + SQR(p5 - mu) + SQR(p6 - mu) + SQR(p7 - mu) + SQR(p8 - mu) + SQR(p9 - mu);
                p1 -= rix[0][-4];
                p2 -= rix[0][-3];
                p3 -= rix[0][-2];
                p4 -= rix[0][-1];
                p5 -= rix[0][ 0];
                p6 -= rix[0][ 1];
                p7 -= rix[0][ 2];
                p8 -= rix[0][ 3];
                p9 -= rix[0][ 4];
                float vn = 1e-7 + SQR(p1) + SQR(p2) + SQR(p3) + SQR(p4) + SQR(p5) + SQR(p6) + SQR(p7) + SQR(p8) + SQR(p9);
                float xh = (rix[0][0] * vx + rix[2][0] * vn) / (vx + vn);
                float vh = vx * vn / (vx + vn);

                // vertical
                p1 = rix[3][-w4];
                p2 = rix[3][-w3];
                p3 = rix[3][-w2];
                p4 = rix[3][-w1];
                p5 = rix[3][  0];
                p6 = rix[3][ w1];
                p7 = rix[3][ w2];
                p8 = rix[3][ w3];
                p9 = rix[3][ w4];
                mu = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) / 9.f;
                vx = 1e-7 + SQR(p1 - mu) + SQR(p2 - mu) + SQR(p3 - mu) + SQR(p4 - mu) + SQR(p5 - mu) + SQR(p6 - mu) + SQR(p7 - mu) + SQR(p8 - mu) + SQR(p9 - mu);
                p1 -= rix[1][-w4];
                p2 -= rix[1][-w3];
                p3 -= rix[1][-w2];
                p4 -= rix[1][-w1];
                p5 -= rix[1][  0];
                p6 -= rix[1][ w1];
                p7 -= rix[1][ w2];
                p8 -= rix[1][ w3];
                p9 -= rix[1][ w4];
                vn = 1e-7 + SQR(p1) + SQR(p2) + SQR(p3) + SQR(p4) + SQR(p5) + SQR(p6) + SQR(p7) + SQR(p8) + SQR(p9);
                float xv = (rix[1][0] * vx + rix[3][0] * vn) / (vx + vn);
                float vv = vx * vn / (vx + vn);
                // interpolated G-R(B)
                rix[4][0] = (xh * vv + xv * vh) / (vh + vv);
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.4);
            }
        }

        // copy CFA values
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 0; rr < rr1; rr++)
            for (int cc = 0, row = rr - ba; cc < cc1; cc++) {
                int col = cc - ba;
                int c = FC(rr, cc);
                rix[c] = qix[c] + rr * cc1 + cc;

                if ((row >= 0) & (row < height) & (col >= 0) & (col < width)) {
                    rix[c][0] = (*gamtab)[rawData[row][col]];
                } else {
                    rix[c][0] = 0.f;
                }

                if (c != 1) {
                    rix[1] = qix[1] + rr * cc1 + cc;
                    rix[4] = qix[4] + rr * cc1 + cc;
                    rix[1][0] = rix[c][0] + rix[4][0];
                }
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.5);
            }
        }

        // bilinear interpolation for R/B
        // interpolate R/B at G location
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 1; rr < rr1 - 1; rr++)
            for (int cc = 1 + (FC(rr, 2) & 1), c = FC(rr, cc + 1); cc < cc1 - 1; cc += 2) {
                rix[c] = qix[c] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[c][0] = rix[1][0] + xdiv2f(rix[c][ -1] - rix[1][ -1] + rix[c][ 1] - rix[1][ 1]);
                c = 2 - c;
                rix[c] = qix[c] + rr * cc1 + cc;
                rix[c][0] = rix[1][0] + xdiv2f(rix[c][-w1] - rix[1][-w1] + rix[c][w1] - rix[1][w1]);
                c = 2 - c;
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.6);
            }
        }

        // interpolate R/B at B/R location
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int rr = 1; rr < rr1 - 1; rr++)
            for (int cc = 1 + (FC(rr, 1) & 1), c = 2 - FC(rr, cc); cc < cc1 - 1; cc += 2) {
                rix[c] = qix[c] + rr * cc1 + cc;
                rix[1] = qix[1] + rr * cc1 + cc;
                rix[c][0] = rix[1][0] + x0250(rix[c][-w1] - rix[1][-w1] + rix[c][ -1] - rix[1][ -1] + rix[c][  1] - rix[1][  1] + rix[c][ w1] - rix[1][ w1]);
            }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if (plistener) {
                plistener->setProgress (0.7);
            }
        }

    }// End of parallelization 1

    // median filter/
    for (int pass = 0; pass < iter; pass++) {
        // Apply 3x3 median filter
        // Compute median(R-G) and median(B-G)

#ifdef _OPENMP
        #pragma omp parallel for private(rix)
#endif

        for (int rr = 1; rr < rr1 - 1; rr++) {
            for (int c = 0; c < 3; c += 2) {
                int d = c + 3 - (c == 0 ? 0 : 1);
                int cc = 1;
#ifdef __SSE2__

                for (; cc < cc1 - 4; cc += 4) {
                    rix[d] = qix[d] + rr * cc1 + cc;
                    rix[c] = qix[c] + rr * cc1 + cc;
                    rix[1] = qix[1] + rr * cc1 + cc;
                    // Assign 3x3 differential color values
                    const std::array<vfloat, 9> p = {
                        LVFU(rix[c][-w1 - 1]) - LVFU(rix[1][-w1 - 1]),
                        LVFU(rix[c][-w1]) - LVFU(rix[1][-w1]),
                        LVFU(rix[c][-w1 + 1]) - LVFU(rix[1][-w1 + 1]),
                        LVFU(rix[c][   -1]) - LVFU(rix[1][   -1]),
                        LVFU(rix[c][  0]) - LVFU(rix[1][  0]),
                        LVFU(rix[c][    1]) - LVFU(rix[1][    1]),
                        LVFU(rix[c][ w1 - 1]) - LVFU(rix[1][ w1 - 1]),
                        LVFU(rix[c][ w1]) - LVFU(rix[1][ w1]),
                        LVFU(rix[c][ w1 + 1]) - LVFU(rix[1][ w1 + 1])
                    };
                    _mm_storeu_ps(&rix[d][0], median(p));
                }

#endif

                for (; cc < cc1 - 1; cc++) {
                    rix[d] = qix[d] + rr * cc1 + cc;
                    rix[c] = qix[c] + rr * cc1 + cc;
                    rix[1] = qix[1] + rr * cc1 + cc;
                    // Assign 3x3 differential color values
                    const std::array<float, 9> p = {
                        rix[c][-w1 - 1] - rix[1][-w1 - 1],
                        rix[c][-w1] - rix[1][-w1],
                        rix[c][-w1 + 1] - rix[1][-w1 + 1],
                        rix[c][   -1] - rix[1][   -1],
                        rix[c][  0] - rix[1][  0],
                        rix[c][    1] - rix[1][    1],
                        rix[c][ w1 - 1] - rix[1][ w1 - 1],
                        rix[c][ w1] - rix[1][ w1],
                        rix[c][ w1 + 1] - rix[1][ w1 + 1]
                    };
                    rix[d][0] = median(p);
                }
            }
        }

        // red/blue at GREEN pixel locations & red/blue and green at BLUE/RED pixel locations
#ifdef _OPENMP
        #pragma omp parallel for private (rix)
#endif

        for (int rr = 0; rr < rr1; rr++) {
            rix[0] = qix[0] + rr * cc1;
            rix[1] = qix[1] + rr * cc1;
            rix[2] = qix[2] + rr * cc1;
            rix[3] = qix[3] + rr * cc1;
            rix[4] = qix[4] + rr * cc1;
            int c0 = FC(rr, 0);
            int c1 = FC(rr, 1);

            if(c0 == 1) {
                c1 = 2 - c1;
                int d = c1 + 3 - (c1 == 0 ? 0 : 1);
                int cc;

                for (cc = 0; cc < cc1 - 1; cc += 2) {
                    rix[0][0] = rix[1][0] + rix[3][0];
                    rix[2][0] = rix[1][0] + rix[4][0];
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                    rix[c1][0] = rix[1][0] + rix[d][0];
                    rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                }

                if(cc < cc1) { // remaining pixel, only if width is odd
                    rix[0][0] = rix[1][0] + rix[3][0];
                    rix[2][0] = rix[1][0] + rix[4][0];
                }
            } else {
                c0 = 2 - c0;
                int d = c0 + 3 - (c0 == 0 ? 0 : 1);
                int cc;

                for (cc = 0; cc < cc1 - 1; cc += 2) {
                    rix[c0][0] = rix[1][0] + rix[d][0];
                    rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                    rix[0][0] = rix[1][0] + rix[3][0];
                    rix[2][0] = rix[1][0] + rix[4][0];
                    rix[0]++;
                    rix[1]++;
                    rix[2]++;
                    rix[3]++;
                    rix[4]++;
                }

                if(cc < cc1) { // remaining pixel, only if width is odd
                    rix[c0][0] = rix[1][0] + rix[d][0];
                    rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
                }
            }
        }
    }

    if (plistener) {
        plistener->setProgress (0.8);
    }

    if(applyGamma) {
        gamtab = &(Color::igammatab_24_17);
    } else {
        gamtab->makeIdentity();
    }

    array2D<float>* rgb[3];
    rgb[0] = &red;
    rgb[1] = &green;
    rgb[2] = &blue;

    // copy result back to image matrix
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int row = 0; row < height; row++) {
        for (int col = 0, rr = row + ba; col < width; col++) {
            int cc = col + ba;
            int c = FC(row, col);

            for (int ii = 0; ii < 3; ii++)
                if (ii != c) {
                    float *rix = qix[ii] + rr * cc1 + cc;
                    (*(rgb[ii]))[row][col] = (*gamtab)[65535.f * rix[0]];
                } else {
                    (*(rgb[ii]))[row][col] = CLIP(rawData[row][col]);
                }
        }
    }

    if (plistener) {
        plistener->setProgress (1.0);
    }

    if(buffer) {
        free(buffer);
    } else
        for(int i = 0; i < 5; i++) {
            free(qix[i]);
        }

    if(!applyGamma) {
        delete gamtab;
    }

    if(iterations > 4 && iterations <= 6) {
        refinement(passref);
    } else if(iterations > 6) {
        refinement_lassus(passref);
    }

}

/*
   Refinement based on EECI demosaicing algorithm by L. Chang and Y.P. Tan
   Paul Lee
   Adapted for RawTherapee - Jacques Desmis 04/2013
*/

#ifdef __SSE2__
#define CLIPV(a) LIMV(a,ZEROV,c65535v)
#endif
void RawImageSource::refinement(int PassCount)
{

    int width = W;
    int height = H;
    int w1 = width;
    int w2 = 2 * w1;

    if (plistener) {
        plistener->setProgressStr (M("TP_RAW_DMETHOD_PROGRESSBAR_REFINE"));
    }

    array2D<float> *rgb[3];
    rgb[0] = &red;
    rgb[1] = &green;
    rgb[2] = &blue;

    for (int b = 0; b < PassCount; b++) {
        if (plistener) {
            plistener->setProgress ((float)b / PassCount);
        }


#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float *pix[3];

            /* Reinforce interpolated green pixels on RED/BLUE pixel locations */
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 2; row < height - 2; row++) {
                int col = 2 + (FC(row, 2) & 1);
                int c = FC(row, col);
#ifdef __SSE2__
                __m128 dLv, dRv, dUv, dDv, v0v;
                __m128 onev = _mm_set1_ps(1.f);
                __m128 zd5v = _mm_set1_ps(0.5f);
                __m128 c65535v = _mm_set1_ps(65535.f);

                for (; col < width - 8; col += 8) {
                    int indx = row * width + col;
                    pix[c] = (float*)(*rgb[c]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    dLv = onev / (onev + vabsf(LC2VFU(pix[c][ -2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dRv = onev / (onev + vabsf(LC2VFU(pix[c][  2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dUv = onev / (onev + vabsf(LC2VFU(pix[c][-w2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    dDv = onev / (onev + vabsf(LC2VFU(pix[c][ w2]) - LC2VFU(pix[c][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    v0v = CLIPV(LC2VFU(pix[c][0]) + zd5v + ((LC2VFU(pix[1][-1]) - LC2VFU(pix[c][-1])) * dLv + (LC2VFU(pix[1][1]) - LC2VFU(pix[c][1])) * dRv + (LC2VFU(pix[1][-w1]) - LC2VFU(pix[c][-w1])) * dUv + (LC2VFU(pix[1][w1]) - LC2VFU(pix[c][w1])) * dDv ) / (dLv + dRv + dUv + dDv));
                    STC2VFU(pix[1][0], v0v);
                }

#endif

                for (; col < width - 2; col += 2) {
                    int indx = row * width + col;
                    pix[c] = (float*)(*rgb[c]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    float dL = 1.f / (1.f + fabsf(pix[c][ -2] - pix[c][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dR = 1.f / (1.f + fabsf(pix[c][  2] - pix[c][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dU = 1.f / (1.f + fabsf(pix[c][-w2] - pix[c][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float dD = 1.f / (1.f + fabsf(pix[c][ w2] - pix[c][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float v0 = (pix[c][0] + 0.5f + ((pix[1][ -1] - pix[c][ -1]) * dL + (pix[1][  1] - pix[c][  1]) * dR + (pix[1][-w1] - pix[c][-w1]) * dU + (pix[1][ w1] - pix[c][ w1]) * dD ) / (dL + dR + dU + dD));
                    pix[1][0] = CLIP(v0);
                }
            }

            /* Reinforce interpolated red/blue pixels on GREEN pixel locations */
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 2; row < height - 2; row++) {
                int col = 2 + (FC(row, 3) & 1);
                int c = FC(row, col + 1);
#ifdef __SSE2__
                __m128 dLv, dRv, dUv, dDv, v0v;
                __m128 onev = _mm_set1_ps(1.f);
                __m128 zd5v = _mm_set1_ps(0.5f);
                __m128 c65535v = _mm_set1_ps(65535.f);

                for (; col < width - 8; col += 8) {
                    int indx = row * width + col;
                    pix[1] = (float*)(*rgb[1]) + indx;

                    for (int i = 0; i < 2; c = 2 - c, i++) {
                        pix[c] = (float*)(*rgb[c]) + indx;
                        dLv = onev / (onev + vabsf(LC2VFU(pix[1][ -2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][ 1]) - LC2VFU(pix[c][ -1])));
                        dRv = onev / (onev + vabsf(LC2VFU(pix[1][  2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][ 1]) - LC2VFU(pix[c][ -1])));
                        dUv = onev / (onev + vabsf(LC2VFU(pix[1][-w2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][w1]) - LC2VFU(pix[c][-w1])));
                        dDv = onev / (onev + vabsf(LC2VFU(pix[1][ w2]) - LC2VFU(pix[1][0])) + vabsf(LC2VFU(pix[c][w1]) - LC2VFU(pix[c][-w1])));
                        v0v = CLIPV(LC2VFU(pix[1][0]) + zd5v - ((LC2VFU(pix[1][-1]) - LC2VFU(pix[c][-1])) * dLv + (LC2VFU(pix[1][1]) - LC2VFU(pix[c][1])) * dRv + (LC2VFU(pix[1][-w1]) - LC2VFU(pix[c][-w1])) * dUv + (LC2VFU(pix[1][w1]) - LC2VFU(pix[c][w1])) * dDv ) / (dLv + dRv + dUv + dDv));
                        STC2VFU(pix[c][0], v0v);
                    }
                }

#endif

                for (; col < width - 2; col += 2) {
                    int indx = row * width + col;
                    pix[1] = (float*)(*rgb[1]) + indx;

                    for (int i = 0; i < 2; c = 2 - c, i++) {
                        pix[c] = (float*)(*rgb[c]) + indx;
                        float dL = 1.f / (1.f + fabsf(pix[1][ -2] - pix[1][0]) + fabsf(pix[c][ 1] - pix[c][ -1]));
                        float dR = 1.f / (1.f + fabsf(pix[1][  2] - pix[1][0]) + fabsf(pix[c][ 1] - pix[c][ -1]));
                        float dU = 1.f / (1.f + fabsf(pix[1][-w2] - pix[1][0]) + fabsf(pix[c][w1] - pix[c][-w1]));
                        float dD = 1.f / (1.f + fabsf(pix[1][ w2] - pix[1][0]) + fabsf(pix[c][w1] - pix[c][-w1]));
                        float v0 = (pix[1][0] + 0.5f - ((pix[1][ -1] - pix[c][ -1]) * dL + (pix[1][  1] - pix[c][  1]) * dR + (pix[1][-w1] - pix[c][-w1]) * dU + (pix[1][ w1] - pix[c][ w1]) * dD ) / (dL + dR + dU + dD));
                        pix[c][0] = CLIP(v0);
                    }
                }
            }

            /* Reinforce integrated red/blue pixels on BLUE/RED pixel locations */
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 2; row < height - 2; row++) {
                int col = 2 + (FC(row, 2) & 1);
                int c = 2 - FC(row, col);
#ifdef __SSE2__
                __m128 dLv, dRv, dUv, dDv, v0v;
                __m128 onev = _mm_set1_ps(1.f);
                __m128 zd5v = _mm_set1_ps(0.5f);
                __m128 c65535v = _mm_set1_ps(65535.f);

                for (; col < width - 8; col += 8) {
                    int indx = row * width + col;
                    pix[0] = (float*)(*rgb[0]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    pix[2] = (float*)(*rgb[2]) + indx;
                    int d = 2 - c;
                    dLv = onev / (onev + vabsf(LC2VFU(pix[d][ -2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dRv = onev / (onev + vabsf(LC2VFU(pix[d][  2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][ 1]) - LC2VFU(pix[1][ -1])));
                    dUv = onev / (onev + vabsf(LC2VFU(pix[d][-w2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    dDv = onev / (onev + vabsf(LC2VFU(pix[d][ w2]) - LC2VFU(pix[d][0])) + vabsf(LC2VFU(pix[1][w1]) - LC2VFU(pix[1][-w1])));
                    v0v = CLIPV(LC2VFU(pix[1][0]) + zd5v - ((LC2VFU(pix[1][-1]) - LC2VFU(pix[c][-1])) * dLv + (LC2VFU(pix[1][1]) - LC2VFU(pix[c][1])) * dRv + (LC2VFU(pix[1][-w1]) - LC2VFU(pix[c][-w1])) * dUv + (LC2VFU(pix[1][w1]) - LC2VFU(pix[c][w1])) * dDv ) / (dLv + dRv + dUv + dDv));
                    STC2VFU(pix[c][0], v0v);
                }

#endif

                for (; col < width - 2; col += 2) {
                    int indx = row * width + col;
                    pix[0] = (float*)(*rgb[0]) + indx;
                    pix[1] = (float*)(*rgb[1]) + indx;
                    pix[2] = (float*)(*rgb[2]) + indx;
                    int d = 2 - c;
                    float dL = 1.f / (1.f + fabsf(pix[d][ -2] - pix[d][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dR = 1.f / (1.f + fabsf(pix[d][  2] - pix[d][0]) + fabsf(pix[1][ 1] - pix[1][ -1]));
                    float dU = 1.f / (1.f + fabsf(pix[d][-w2] - pix[d][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float dD = 1.f / (1.f + fabsf(pix[d][ w2] - pix[d][0]) + fabsf(pix[1][w1] - pix[1][-w1]));
                    float v0 = (pix[1][0] + 0.5f - ((pix[1][ -1] - pix[c][ -1]) * dL + (pix[1][  1] - pix[c][  1]) * dR + (pix[1][-w1] - pix[c][-w1]) * dU + (pix[1][ w1] - pix[c][ w1]) * dD ) / (dL + dR + dU + dD));
                    pix[c][0] = CLIP(v0);
                }
            }
        } // end parallel
    }

}
#ifdef __SSE2__
#undef CLIPV
#endif


// Refinement based on EECI demozaicing algorithm by L. Chang and Y.P. Tan
// from "Lassus" : Luis Sanz Rodriguez, adapted by Jacques Desmis - JDC - and Oliver Duis for RawTherapee
// increases the signal to noise ratio (PSNR) # +1 to +2 dB : tested with Dcraw : eg: Lighthouse + AMaZE : whitout refinement:39.96dB, with refinement:41.86 dB
// reduce color artifacts, improves the interpolation
// but it's relatively slow
//
// Should be DISABLED if it decreases image quality by increases some image noise and generates blocky edges
void RawImageSource::refinement_lassus(int PassCount)
{
    // const int PassCount=1;

    // if (settings->verbose) printf("Refinement\n");

    int u = W, v = 2 * u, w = 3 * u, x = 4 * u, y = 5 * u;
    float (*image)[3];
    image = (float(*)[3]) calloc(W * H, sizeof * image);
#ifdef _OPENMP
    #pragma omp parallel shared(image)
#endif
    {
        // convert red, blue, green to image
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                image[i * W + j][0] = red  [i][j];
                image[i * W + j][1] = green[i][j];
                image[i * W + j][2] = blue [i][j];
            }
        }

        for (int b = 0; b < PassCount; b++) {
            if (plistener) {
                plistener->setProgressStr (M("TP_RAW_DMETHOD_PROGRESSBAR_REFINE"));
                plistener->setProgress ((float)b / PassCount);
            }

            // Reinforce interpolated green pixels on RED/BLUE pixel locations
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 6; row < H - 6; row++) {
                for (int col = 6 + (FC(row, 2) & 1), c = FC(row, col); col < W - 6; col += 2) {
                    float (*pix)[3] = image + row * W + col;

                    // Cubic Spline Interpolation by Li and Randhawa, modified by Luis Sanz Rodriguez

                    float f[4];
                    f[0] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-v][c]) - x0875(pix[0][c]) - x0250(pix[-x][c]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[-w][1])) + fabs(x0875(pix[-w][1]) - x1125(pix[-u][1]) + x0250(pix[-y][1])));
                    f[1] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+2][c]) - x0875(pix[0][c]) - x0250(pix[+4][c]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[+3][1])) + fabs(x0875(pix[+3][1]) - x1125(pix[+1][1]) + x0250(pix[+5][1])));
                    f[2] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-2][c]) - x0875(pix[0][c]) - x0250(pix[-4][c]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[-3][1])) + fabs(x0875(pix[-3][1]) - x1125(pix[-1][1]) + x0250(pix[-5][1])));
                    f[3] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+v][c]) - x0875(pix[0][c]) - x0250(pix[+x][c]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[+w][1])) + fabs(x0875(pix[+w][1]) - x1125(pix[+u][1]) + x0250(pix[+y][1])));

                    float g[4];//CLIREF avoid overflow
                    g[0] = pix[0][c] + (x0875(CLIREF(pix[-u][1] - pix[-u][c])) + x0125(CLIREF(pix[+u][1] - pix[+u][c])));
                    g[1] = pix[0][c] + (x0875(CLIREF(pix[+1][1] - pix[+1][c])) + x0125(CLIREF(pix[-1][1] - pix[-1][c])));
                    g[2] = pix[0][c] + (x0875(CLIREF(pix[-1][1] - pix[-1][c])) + x0125(CLIREF(pix[+1][1] - pix[+1][c])));
                    g[3] = pix[0][c] + (x0875(CLIREF(pix[+u][1] - pix[+u][c])) + x0125(CLIREF(pix[-u][1] - pix[-u][c])));

                    pix[0][1] = (f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]);

                }
            }

            // Reinforce interpolated red/blue pixels on GREEN pixel locations
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 6; row < H - 6; row++) {
                for (int col = 6 + (FC(row, 3) & 1), c = FC(row, col + 1); col < W - 6; col += 2) {
                    float (*pix)[3] = image + row * W + col;

                    for (int i = 0; i < 2; c = 2 - c, i++) {
                        float f[4];
                        f[0] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[-v][1]) - x1125(pix[0][1]) + x0250(pix[-x][1]))) + fabs(pix[u] [c] - pix[-u][c]) + fabs(pix[-w][c] - pix[-u][c]));
                        f[1] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[+2][1]) - x1125(pix[0][1]) + x0250(pix[+4][1]))) + fabs(pix[+1][c] - pix[-1][c]) + fabs(pix[+3][c] - pix[+1][c]));
                        f[2] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[-2][1]) - x1125(pix[0][1]) + x0250(pix[-4][1]))) + fabs(pix[+1][c] - pix[-1][c]) + fabs(pix[-3][c] - pix[-1][c]));
                        f[3] = 1.0f / (1.0f + xmul2f(fabs(x0875(pix[+v][1]) - x1125(pix[0][1]) + x0250(pix[+x][1]))) + fabs(pix[u] [c] - pix[-u][c]) + fabs(pix[+w][c] - pix[+u][c]));

                        float g[5];//CLIREF avoid overflow
                        g[0] = CLIREF(pix[-u][1] - pix[-u][c]);
                        g[1] = CLIREF(pix[+1][1] - pix[+1][c]);
                        g[2] = CLIREF(pix[-1][1] - pix[-1][c]);
                        g[3] = CLIREF(pix[+u][1] - pix[+u][c]);
                        g[4] = ((f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]));
                        pix[0][c] = pix[0][1] - (0.65f * g[4] + 0.35f * CLIREF(pix[0][1] - pix[0][c]));
                    }
                }
            }

            // Reinforce integrated red/blue pixels on BLUE/RED pixel locations
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 6; row < H - 6; row++) {
                for (int col = 6 + (FC(row, 2) & 1), c = 2 - FC(row, col), d = 2 - c; col < W - 6; col += 2) {
                    float (*pix)[3] = image + row * W + col;

                    float f[4];
                    f[0] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-v][d]) - x0875(pix[0][d]) - x0250(pix[-x][d]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[-w][1])) + fabs(x0875(pix[-w][1]) - x1125(pix[-u][1]) + x0250(pix[-y][1])));
                    f[1] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+2][d]) - x0875(pix[0][d]) - x0250(pix[+4][d]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[+3][1])) + fabs(x0875(pix[+3][1]) - x1125(pix[+1][1]) + x0250(pix[+5][1])));
                    f[2] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[-2][d]) - x0875(pix[0][d]) - x0250(pix[-4][d]))) + fabs(x0875(pix[1][1]) - x1125(pix[-1][1]) + x0250(pix[-3][1])) + fabs(x0875(pix[-3][1]) - x1125(pix[-1][1]) + x0250(pix[-5][1])));
                    f[3] = 1.0f / (1.0f + xmul2f(fabs(x1125(pix[+v][d]) - x0875(pix[0][d]) - x0250(pix[+x][d]))) + fabs(x0875(pix[u][1]) - x1125(pix[-u][1]) + x0250(pix[+w][1])) + fabs(x0875(pix[+w][1]) - x1125(pix[+u][1]) + x0250(pix[+y][1])));

                    float g[5];
                    g[0] = (x0875((pix[-u][1] - pix[-u][c])) + x0125((pix[-v][1] - pix[-v][c])));
                    g[1] = (x0875((pix[+1][1] - pix[+1][c])) + x0125((pix[+2][1] - pix[+2][c])));
                    g[2] = (x0875((pix[-1][1] - pix[-1][c])) + x0125((pix[-2][1] - pix[-2][c])));
                    g[3] = (x0875((pix[+u][1] - pix[+u][c])) + x0125((pix[+v][1] - pix[+v][c])));

                    g[4] = (f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3]) / (f[0] + f[1] + f[2] + f[3]);

                    const std::array<float, 9> p = {
                        pix[-u - 1][1] - pix[-u - 1][c],
                        pix[-u + 0][1] - pix[-u + 0][c],
                        pix[-u + 1][1] - pix[-u + 1][c],
                        pix[+0 - 1][1] - pix[+0 - 1][c],
                        pix[+0 + 0][1] - pix[+0 + 0][c],
                        pix[+0 + 1][1] - pix[+0 + 1][c],
                        pix[+u - 1][1] - pix[+u - 1][c],
                        pix[+u + 0][1] - pix[+u + 0][c],
                        pix[+u + 1][1] - pix[+u + 1][c]
                    };

                    const float med = median(p);

                    pix[0][c] = LIM(pix[0][1] - (1.30f * g[4] - 0.30f * (pix[0][1] - pix[0][c])), 0.99f * (pix[0][1] - med), 1.01f * (pix[0][1] - med));

                }
            }

        }

        // put modified values to red, green, blue
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                red  [i][j] = image[i * W + j][0];
                green[i][j] = image[i * W + j][1];
                blue [i][j] = image[i * W + j][2];
            }
        }
    }

    free(image);

}
} /* namespace */
