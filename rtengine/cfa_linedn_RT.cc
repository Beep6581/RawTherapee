////////////////////////////////////////////////////////////////
//
//          CFA line denoise by DCT filtering
//
//  copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//  parallelized 2013 by Ingo Weyrich
//
// code dated: June 7, 2010
//
//  cfa_linedn_RT.cc is free software: you can redistribute it and/or modify
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

#include <cmath>

#include "rtengine.h"
#include "rawimagesource.h"
#include "rt_math.h"

#define TS 224      // Tile size of 224 instead of 512 speeds up processing

#define CLASS

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace std;
using namespace rtengine;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::CLASS cfa_linedn(float noise)
{
    // local variables
    int height = H, width = W;

    const float clip_pt = 0.8 * initialGain * 65535.0;

    const float eps = 1e-5;       //tolerance to avoid dividing by zero

    const float gauss[5] = {0.20416368871516755, 0.18017382291138087, 0.1238315368057753, 0.0662822452863612, 0.02763055063889883};
    const float rolloff[8] = {0, 0.135335, 0.249352, 0.411112, 0.606531, 0.800737, 0.945959, 1}; //gaussian with sigma=3
    const float window[8] = {0, .25, .75, 1, 1, .75, .25, 0}; //sine squared

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (plistener) {
        plistener->setProgressStr ("Line Denoise...");
        plistener->setProgress (0.0);
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float noisevar = SQR(3 * noise * 65535); // _noise_ (as a fraction of saturation) is input to the algorithm
    float noisevarm4 = 4.0f * noisevar;
    volatile double progress = 0.0;
    float* RawDataTmp = (float*)malloc( width * height * sizeof(float));
    #pragma omp parallel
    {

        // allocate memory and assure the arrays don't have same 64 byte boundary to avoid L1 conflict misses
        float *cfain = (float*)malloc(4 * TS * TS * sizeof(float) + 3 * 16 * sizeof(float));
        float *cfablur = (cfain + (TS * TS) + 1 * 16);
        float *cfadiff = (cfain + (2 * TS * TS) + 2 * 16);
        float *cfadn = (cfain + (3 * TS * TS) + 3 * 16);


        float linehvar[4], linevvar[4], noisefactor[4][8][2], coeffsq;
        float dctblock[4][8][8];

        #pragma omp for

        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++) {
                RawDataTmp[i * width + j] = rawData[i][j];
            }

        // Main algorithm: Tile loop
        #pragma omp for schedule(dynamic) collapse(2)

        for (int top = 0; top < height - 16; top += TS - 32)
            for (int left = 0; left < width - 16; left += TS - 32) {

                int bottom = min(top + TS, height);
                int right  = min(left + TS, width);
                int numrows = bottom - top;
                int numcols = right - left;
                int indx1;

                // load CFA data; data should be in linear gamma space, before white balance multipliers are applied
                for (int rr = top; rr < top + numrows; rr++)
                    for (int cc = left, indx = (rr - top) * TS; cc < left + numcols; cc++, indx++) {
                        cfain[indx] = rawData[rr][cc];
                    }

                //pad the block to a multiple of 16 on both sides

                if (numcols < TS) {
                    indx1 = numcols % 16;

                    for (int i = 0; i < (16 - indx1); i++)
                        for (int rr = 0; rr < numrows; rr++) {
                            cfain[(rr)*TS + numcols + i] = cfain[(rr) * TS + numcols - i - 1];
                        }

                    numcols += 16 - indx1;
                }

                if (numrows < TS) {
                    indx1 = numrows % 16;

                    for (int i = 0; i < (16 - indx1); i++)
                        for (int cc = 0; cc < numcols; cc++) {
                            cfain[(numrows + i)*TS + cc] = cfain[(numrows - i - 1) * TS + cc];
                        }

                    numrows += 16 - indx1;
                }

                //The cleaning algorithm starts here

                //gaussian blur of CFA data
                for (int rr = 8; rr < numrows - 8; rr++) {
                    for (int indx = rr * TS; indx < rr * TS + numcols; indx++) {
                        cfablur[indx] = gauss[0] * cfain[indx];

                        for (int i = 1; i < 5; i++) {
                            cfablur[indx] += gauss[i] * (cfain[indx - (2 * i) * TS] + cfain[indx + (2 * i) * TS]);
                        }
                    }

                    for (int indx = rr * TS + 8; indx < rr * TS + numcols - 8; indx++) {
                        cfadn[indx] = gauss[0] * cfablur[indx];

                        for (int i = 1; i < 5; i++) {
                            cfadn[indx] += gauss[i] * (cfablur[indx - 2 * i] + cfablur[indx + 2 * i]);
                        }

                        cfadiff[indx] = cfain[indx] - cfadn[indx]; // hipass cfa data
                    }
                }

                //begin block DCT
                for (int rr = 8; rr < numrows - 22; rr += 8) // (rr,cc) shift by 8 to overlap blocks
                    for (int cc = 8; cc < numcols - 22; cc += 8) {
                        for (int ey = 0; ey < 2; ey++) // (ex,ey) specify RGGB subarray
                            for (int ex = 0; ex < 2; ex++) {
                                //grab an 8x8 block of a given RGGB channel
                                for (int i = 0; i < 8; i++)
                                    for (int j = 0; j < 8; j++) {
                                        dctblock[2 * ey + ex][i][j] = cfadiff[(rr + 2 * i + ey) * TS + cc + 2 * j + ex];
                                    }

                                ddct8x8s(-1, dctblock[2 * ey + ex]); //forward DCT
                            }

                        for (int ey = 0; ey < 2; ey++) // (ex,ey) specify RGGB subarray
                            for (int ex = 0; ex < 2; ex++) {
                                linehvar[2 * ey + ex] = linevvar[2 * ey + ex] = 0;

                                for (int i = 4; i < 8; i++) {
                                    linehvar[2 * ey + ex] += SQR(dctblock[2 * ey + ex][0][i]);
                                    linevvar[2 * ey + ex] += SQR(dctblock[2 * ey + ex][i][0]);
                                }

                                //Wiener filter for line denoising; roll off low frequencies
                                for (int i = 1; i < 8; i++) {
                                    coeffsq = SQR(dctblock[2 * ey + ex][i][0]); //vertical
                                    noisefactor[2 * ey + ex][i][0] = coeffsq / (coeffsq + rolloff[i] * noisevar + eps);
                                    coeffsq = SQR(dctblock[2 * ey + ex][0][i]); //horizontal
                                    noisefactor[2 * ey + ex][i][1] = coeffsq / (coeffsq + rolloff[i] * noisevar + eps);
                                    // noisefactor labels are [RGGB subarray][row/col position][0=vert,1=hor]
                                }
                            }

                        //horizontal lines
                        if (noisevarm4 > (linehvar[0] + linehvar[1])) { //horizontal lines
                            for (int i = 1; i < 8; i++) {
                                dctblock[0][0][i] *= 0.5f * (noisefactor[0][i][1] + noisefactor[1][i][1]); //or should we use MIN???
                                dctblock[1][0][i] *= 0.5f * (noisefactor[0][i][1] + noisefactor[1][i][1]); //or should we use MIN???
                            }
                        }

                        if (noisevarm4 > (linehvar[2] + linehvar[3])) { //horizontal lines
                            for (int i = 1; i < 8; i++) {
                                dctblock[2][0][i] *= 0.5f * (noisefactor[2][i][1] + noisefactor[3][i][1]); //or should we use MIN???
                                dctblock[3][0][i] *= 0.5f * (noisefactor[2][i][1] + noisefactor[3][i][1]); //or should we use MIN???
                            }
                        }

                        //vertical lines
                        if (noisevarm4 > (linevvar[0] + linevvar[2])) { //vertical lines
                            for (int i = 1; i < 8; i++) {
                                dctblock[0][i][0] *= 0.5f * (noisefactor[0][i][0] + noisefactor[2][i][0]); //or should we use MIN???
                                dctblock[2][i][0] *= 0.5f * (noisefactor[0][i][0] + noisefactor[2][i][0]); //or should we use MIN???
                            }
                        }

                        if (noisevarm4 > (linevvar[1] + linevvar[3])) { //vertical lines
                            for (int i = 1; i < 8; i++) {
                                dctblock[1][i][0] *= 0.5f * (noisefactor[1][i][0] + noisefactor[3][i][0]); //or should we use MIN???
                                dctblock[3][i][0] *= 0.5f * (noisefactor[1][i][0] + noisefactor[3][i][0]); //or should we use MIN???
                            }
                        }

                        for (int ey = 0; ey < 2; ey++) // (ex,ey) specify RGGB subarray
                            for (int ex = 0; ex < 2; ex++) {
                                ddct8x8s(1, dctblock[2 * ey + ex]); //inverse DCT

                                //multiply by window fn and add to output (cfadn)
                                for (int i = 0; i < 8; i++)
                                    for (int j = 0; j < 8; j++) {
                                        cfadn[(rr + 2 * i + ey)*TS + cc + 2 * j + ex] += window[i] * window[j] * dctblock[2 * ey + ex][i][j];
                                    }
                            }
                    }

                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                // copy smoothed results to temporary buffer
                for (int rr = 16; rr < numrows - 16; rr++) {
                    int row = rr + top;

                    for (int col = 16 + left, indx = rr * TS + 16; indx < rr * TS + numcols - 16; indx++, col++) {
                        if (rawData[row][col] < clip_pt && cfadn[indx] < clip_pt) {
                            RawDataTmp[row * width + col] = CLIP((int)(cfadn[indx] + 0.5f));
                        }
                    }
                }

                if(plistener) {
                    progress += (double)((TS - 32) * (TS - 32)) / (height * width);

                    if (progress > 1.0) {
                        progress = 1.0;
                    }

                    plistener->setProgress(progress);
                }

            }

        // clean up
        free(cfain);

// copy temporary buffer back to image matrix
        #pragma omp for

        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++) {
                rawData[i][j] = RawDataTmp[i * width + j];
            }

    } // end of parallel processing

    free(RawDataTmp);
}
#undef TS


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
 Discrete Cosine Transform Code

 Copyright(C) 1997 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
 You may use, copy, modify this code for any purpose and
 without fee. You may distribute this ORIGINAL package.
 */


/*
 Short Discrete Cosine Transform
 data length :8x8
 method      :row-column, radix 4 FFT
 functions
 ddct8x8s  : 8x8 DCT
 function prototypes
 void ddct8x8s(int isgn, float **a);
 */


/*
 -------- 8x8 DCT (Discrete Cosine Transform) / Inverse of DCT --------
 [definition]
 <case1> Normalized 8x8 IDCT
 C[k1][k2] = (1/4) * sum_j1=0^7 sum_j2=0^7
 a[j1][j2] * s[j1] * s[j2] *
 cos(pi*j1*(k1+1/2)/8) *
 cos(pi*j2*(k2+1/2)/8), 0<=k1<8, 0<=k2<8
 (s[0] = 1/sqrt(2), s[j] = 1, j > 0)
 <case2> Normalized 8x8 DCT
 C[k1][k2] = (1/4) * s[k1] * s[k2] * sum_j1=0^7 sum_j2=0^7
 a[j1][j2] *
 cos(pi*(j1+1/2)*k1/8) *
 cos(pi*(j2+1/2)*k2/8), 0<=k1<8, 0<=k2<8
 (s[0] = 1/sqrt(2), s[j] = 1, j > 0)
 [usage]
 <case1>
 ddct8x8s(1, a);
 <case2>
 ddct8x8s(-1, a);
 [parameters]
 a[0...7][0...7] :input/output data (double **)
 output data
 a[k1][k2] = C[k1][k2], 0<=k1<8, 0<=k2<8
 */


/* Cn_kR = sqrt(2.0/n) * cos(pi/2*k/n) */
/* Cn_kI = sqrt(2.0/n) * sin(pi/2*k/n) */
/* Wn_kR = cos(pi/2*k/n) */
/* Wn_kI = sin(pi/2*k/n) */
#define C8_1R   0.49039264020161522456
#define C8_1I   0.09754516100806413392
#define C8_2R   0.46193976625564337806
#define C8_2I   0.19134171618254488586
#define C8_3R   0.41573480615127261854
#define C8_3I   0.27778511650980111237
#define C8_4R   0.35355339059327376220
#define W8_4R   0.70710678118654752440


void RawImageSource::ddct8x8s(int isgn, float a[8][8])
{
    int j;
    float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    float xr, xi;

    if (isgn < 0) {
        for (j = 0; j <= 7; j++) {
            x0r = a[0][j] + a[7][j];
            x1r = a[0][j] - a[7][j];
            x0i = a[2][j] + a[5][j];
            x1i = a[2][j] - a[5][j];
            x2r = a[4][j] + a[3][j];
            x3r = a[4][j] - a[3][j];
            x2i = a[6][j] + a[1][j];
            x3i = a[6][j] - a[1][j];
            xr = x0r + x2r;
            xi = x0i + x2i;
            a[0][j] = C8_4R * (xr + xi);
            a[4][j] = C8_4R * (xr - xi);
            xr = x0r - x2r;
            xi = x0i - x2i;
            a[2][j] = C8_2R * xr - C8_2I * xi;
            a[6][j] = C8_2R * xi + C8_2I * xr;
            xr = W8_4R * (x1i - x3i);
            x1i = W8_4R * (x1i + x3i);
            x3i = x1i - x3r;
            x1i += x3r;
            x3r = x1r - xr;
            x1r += xr;
            a[1][j] = C8_1R * x1r - C8_1I * x1i;
            a[7][j] = C8_1R * x1i + C8_1I * x1r;
            a[3][j] = C8_3R * x3r - C8_3I * x3i;
            a[5][j] = C8_3R * x3i + C8_3I * x3r;
        }

        for (j = 0; j <= 7; j++) {
            x0r = a[j][0] + a[j][7];
            x1r = a[j][0] - a[j][7];
            x0i = a[j][2] + a[j][5];
            x1i = a[j][2] - a[j][5];
            x2r = a[j][4] + a[j][3];
            x3r = a[j][4] - a[j][3];
            x2i = a[j][6] + a[j][1];
            x3i = a[j][6] - a[j][1];
            xr = x0r + x2r;
            xi = x0i + x2i;
            a[j][0] = C8_4R * (xr + xi);
            a[j][4] = C8_4R * (xr - xi);
            xr = x0r - x2r;
            xi = x0i - x2i;
            a[j][2] = C8_2R * xr - C8_2I * xi;
            a[j][6] = C8_2R * xi + C8_2I * xr;
            xr = W8_4R * (x1i - x3i);
            x1i = W8_4R * (x1i + x3i);
            x3i = x1i - x3r;
            x1i += x3r;
            x3r = x1r - xr;
            x1r += xr;
            a[j][1] = C8_1R * x1r - C8_1I * x1i;
            a[j][7] = C8_1R * x1i + C8_1I * x1r;
            a[j][3] = C8_3R * x3r - C8_3I * x3i;
            a[j][5] = C8_3R * x3i + C8_3I * x3r;
        }
    } else {
        for (j = 0; j <= 7; j++) {
            x1r = C8_1R * a[1][j] + C8_1I * a[7][j];
            x1i = C8_1R * a[7][j] - C8_1I * a[1][j];
            x3r = C8_3R * a[3][j] + C8_3I * a[5][j];
            x3i = C8_3R * a[5][j] - C8_3I * a[3][j];
            xr = x1r - x3r;
            xi = x1i + x3i;
            x1r += x3r;
            x3i -= x1i;
            x1i = W8_4R * (xr + xi);
            x3r = W8_4R * (xr - xi);
            xr = C8_2R * a[2][j] + C8_2I * a[6][j];
            xi = C8_2R * a[6][j] - C8_2I * a[2][j];
            x0r = C8_4R * (a[0][j] + a[4][j]);
            x0i = C8_4R * (a[0][j] - a[4][j]);
            x2r = x0r - xr;
            x2i = x0i - xi;
            x0r += xr;
            x0i += xi;
            a[0][j] = x0r + x1r;
            a[7][j] = x0r - x1r;
            a[2][j] = x0i + x1i;
            a[5][j] = x0i - x1i;
            a[4][j] = x2r - x3i;
            a[3][j] = x2r + x3i;
            a[6][j] = x2i - x3r;
            a[1][j] = x2i + x3r;
        }

        for (j = 0; j <= 7; j++) {
            x1r = C8_1R * a[j][1] + C8_1I * a[j][7];
            x1i = C8_1R * a[j][7] - C8_1I * a[j][1];
            x3r = C8_3R * a[j][3] + C8_3I * a[j][5];
            x3i = C8_3R * a[j][5] - C8_3I * a[j][3];
            xr = x1r - x3r;
            xi = x1i + x3i;
            x1r += x3r;
            x3i -= x1i;
            x1i = W8_4R * (xr + xi);
            x3r = W8_4R * (xr - xi);
            xr = C8_2R * a[j][2] + C8_2I * a[j][6];
            xi = C8_2R * a[j][6] - C8_2I * a[j][2];
            x0r = C8_4R * (a[j][0] + a[j][4]);
            x0i = C8_4R * (a[j][0] - a[j][4]);
            x2r = x0r - xr;
            x2i = x0i - xi;
            x0r += xr;
            x0i += xi;
            a[j][0] = x0r + x1r;
            a[j][7] = x0r - x1r;
            a[j][2] = x0i + x1i;
            a[j][5] = x0i - x1i;
            a[j][4] = x2r - x3i;
            a[j][3] = x2r + x3i;
            a[j][6] = x2i - x3r;
            a[j][1] = x2i + x3r;
        }
    }
}

