////////////////////////////////////////////////////////////////
//
//			CFA line denoise by DCT filtering
//
//	copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: June 7, 2010
//
//	cfa_linedn_RT.cc is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////


#define TS 512		// Tile size

#define CLASS


/*#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <float.h>
#include <limits.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>*/
#include <math.h>


//#include "shrtdct_float.c"


#define SQR(x) ((x)*(x))
//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
//#define LIM(x,min,max) MAX(min,MIN(x,max))
//#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
//#define CLIP(x) LIM(x,0,65535)

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::CLASS cfa_linedn(float noise)
{  
	// local variables
	int height=H, width=W;
	int top, bottom, left, right, row, col;
	int rr, cc, rr1, cc1, c, indx, i, j;
	int ex, ey; 
	int verbose=1;
	
	float eps=1e-10;			//tolerance to avoid dividing by zero
	
	float gauss[5] = {0.20416368871516755, 0.18017382291138087, 0.1238315368057753, 0.0662822452863612, 0.02763055063889883};
	float rolloff[8] = {0, 0.135335, 0.249352, 0.411112, 0.606531, 0.800737, 0.945959, 1}; //gaussian with sigma=3
	float window[8] = {0, .25, .75, 1, 1, .75, .25, 0}; //sine squared
	float noisevar, linehvar, linevvar, coeffsq;
	
	float aarr[8][8], *dctblock[8];
    for (i = 0; i < 8; i++) dctblock[i] = aarr[i];
	
	char		*buffer;			// TS*TS*16
	float         (*cfain);			// TS*TS*4
	float         (*cfablur);		// TS*TS*4
	float         (*cfadiff);		// TS*TS*4
	float         (*cfadn);			// TS*TS*4
	
	double dt;
	clock_t t1, t2;
	
	//clock_t t1_main,       t2_main       = 0;
	
	// start
	//if (verbose) fprintf (stderr,_("CFA line denoise ...\n"));
	//t1 = clock();
	
	
	// assign working space
	buffer = (char *) malloc(16*TS*TS);
	//merror(buffer,"cfa_linedn()");
	memset(buffer,0,16*TS*TS);
	// rgb array
	cfain			= (float (*))		buffer; //pointers to rows of array
 	cfablur			= (float (*))			(buffer +  4*TS*TS);
	cfadiff			= (float (*))			(buffer +  8*TS*TS);
	cfadn			= (float (*))			(buffer +  12*TS*TS);
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if (plistener) {
		plistener->setProgressStr ("Line Denoise...");
		plistener->setProgress (0.0);
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	noisevar=SQR(3*noise*65535); // _noise_ (as a fraction of saturation) is input to the algorithm

	// Main algorithm: Tile loop
	for (top=0; top < height-16; top += TS-32)
		for (left=0; left < width-16; left += TS-32) {
			bottom = MIN( top+TS,height);
			right  = MIN(left+TS, width);
			rr1 = bottom - top;
			cc1 = right - left;
			// load CFA data; data should be in linear gamma space, before white balance multipliers are applied
			for (rr=0; rr < rr1; rr++)
				for (row=rr+top, cc=0, indx=rr*TS+cc; cc < cc1; cc++, indx++) {
					col = cc+left;
					c = FC(rr,cc);
					cfain[indx] = ri->data[row][col]; 
				}
			//pad the block to a multiple of 16 on both sides
			
			if (cc1 < TS) {
				indx=cc1 % 16;
				for (i=0; i<(16-indx); i++) 
					for (rr=0; rr<rr1; rr++) 
						cfain[(rr)*TS+cc1+i+1]=cfain[(rr)*TS+cc1-i];
				cc1 += 16-indx;
			}
			
			if (rr1 < TS) {
				indx=rr1 % 16;
				for (i=0; i<(16-indx); i++)
					for (cc=0; cc<cc1; cc++)
						cfain[(rr1+i+1)*TS+cc]=cfain[(rr1-i)*TS+cc];
				rr1 += 16-indx;
			}
			
			//The cleaning algorithm starts here
			
			
			//gaussian blur of CFA data
			for (rr=8; rr < rr1-8; rr++)
				for (cc=0, indx=rr*TS+cc; cc < cc1; cc++, indx++) {
					
					cfablur[indx]=gauss[0]*cfain[indx];
					for (i=1; i<5; i++) {
						cfablur[indx] += gauss[i]*(cfain[indx-(2*i)*TS]+cfain[indx+(2*i)*TS]);
					}
				}
			for (rr=8; rr < rr1-8; rr++)
				for (cc=8, indx=rr*TS+cc; cc < cc1-8; cc++, indx++) {
					
					cfadn[indx] = gauss[0]*cfablur[indx];
					for (i=1; i<5; i++) {
						cfadn[indx] += gauss[i]*(cfablur[indx-2*i]+cfablur[indx+2*i]);
					}
					cfadiff[indx]=cfain[indx]-cfadn[indx]; // hipass cfa data
				}
			
			//begin block DCT
			for (ey=0; ey<2; ey++) // (ex,ey) specify RGGB subarray
				for (ex=0; ex<2; ex++) 
					for (rr=8+ey; rr < rr1-22; rr+=8) // (rr,cc) shift by 8 to overlap blocks
						for (cc=8+ex; cc < cc1-22; cc+=8) {
							//grab an 8x8 block of a given RGGB channel
							for (i=0; i<8; i++) 
								for (j=0; j<8; j++) {
									dctblock[i][j]=cfadiff[(rr+2*i)*TS+cc+2*j];
								}
							
							 ddct8x8s(-1, dctblock); //forward DCT
							 
							linehvar=linevvar=0;
							for (i=4; i<8; i++) {
								linehvar += SQR(dctblock[0][i]);
								linevvar += SQR(dctblock[i][0]);
							}
							//Wiener filter for line denoising; roll off low frequencies
							if (noisevar>linehvar) {
								for (i=1; i<8; i++) {
									coeffsq=SQR(dctblock[0][i]);
									dctblock[0][i] *= coeffsq/(coeffsq+rolloff[i]*noisevar+eps);
								}
							}
							if (noisevar>linevvar) {
								for (i=1; i<8; i++) {
									coeffsq=SQR(dctblock[i][0]);
									dctblock[i][0] *= coeffsq/(coeffsq+rolloff[i]*noisevar+eps);
								}
							}
							
							ddct8x8s(1, dctblock); //inverse DCT
							
							//multiply by window fn and add to output (cfadn)
							for (i=0; i<8; i++) 
								for (j=0; j<8; j++) {
									cfadn[(rr+2*i)*TS+cc+2*j] += window[i]*window[j]*dctblock[i][j];
								}
						}	

			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// copy smoothed results back to image matrix
			for (rr=16; rr < rr1-16; rr++)
				for (row=rr+top, cc=16, indx=rr*TS+cc; cc < cc1-16; cc++, indx++) {
					col = cc + left;
					ri->data[row][col] = CLIP((int)(cfadn[indx]+ 0.5)); 
				}
			if(plistener) plistener->setProgress(fabs((float)top/height));
		}
	
	// clean up
	free(buffer);
	
	// done
	/*t2 = clock();
	dt = ((double)(t2-t1)) / CLOCKS_PER_SEC;
	if (verbose) {
		fprintf(stderr,_("elapsed time = %5.3fs\n"),dt);
	}*/
	
	
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


void RawImageSource::ddct8x8s(int isgn, float **a)
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

