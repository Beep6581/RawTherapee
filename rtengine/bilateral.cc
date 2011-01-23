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

#include <math.h>
#include "mytime.h"
#include <omp.h>
#include "buffer.h"
#include <float.h>
#include <string.h>

// This seems ugly, but way faster than any other solutions I tried

#define ELEM(a,b) (src[i - a][j - b] * ec[(int)((src[i - a][j - b]-src[i][j])*scale)+0x10000])
#define SULY(a,b) (ec[(int)((src[i - a][j - b]-src[i][j])*scale)+0x10000])

#define BL_BEGIN(a,b)   float min = FLT_MAX, max = FLT_MIN; \
						for (int i=0; i<H; i++) \
							for (int j=0; j<W; j++) { \
								if (src[i][j]<min) \
									min = src[i][j]; \
								if (src[i][j]>max) \
									max = src[i][j]; \
							} \
						float scale = 0x20000 / (max-min) / 2; \
                        float* ec = new float [0x20000]; \
                        for (int i=0; i<0x20000; i++) \
                            ec[i] = exp(-(double)(i-0x10000)*(double)(i-0x10000) / (2.0*sens*sens) / (scale*scale)); \
                        int rstart = b; \
                        int rend = H-b; \
                        int cstart = b; \
                        int cend = W-b;

#define BL_END(b)       buffer[i][j] = v; }} delete [] ec; \
                        for (int i=0; i<H; i++)  \
                            for (int j=0; j<W; j++)  \
                                if (i<rstart || j<cstart || i>=rend || j>=cend) \
                                    dst[i][j] = src[i][j]; \
                                else \
                                    dst[i][j] = buffer[i][j];

#define BL_OPER3(a11,a12,a21,a22) for (int i=rstart; i<rend; i++) { \
									for (int j=cstart; j<cend; j++) { \
								   float v = a11*ELEM(-1,-1) + a12*ELEM(-1,0) + a11*ELEM(-1,1) + \
                                        a21*ELEM(0,-1) + a22*ELEM(0,0) + a21*ELEM(0,1) + \
                                        a11*ELEM(1,-1) + a12*ELEM(1,0) + a11*ELEM(1,1); \
                                   v /= a11*SULY(-1,-1) + a12*SULY(-1,0) + a11*SULY(-1,1) + \
                                        a21*SULY(0,-1) + a22*SULY(0,0) + a21*SULY(0,1) + \
                                        a11*SULY(1,-1) + a12*SULY(1,0) + a11*SULY(1,1);


#define BL_OPER5(a11,a12,a13,a21,a22,a23,a31,a32,a33) for (int i=rstart; i<rend; i++) { \
														for (int j=cstart; j<cend; j++) { \
													   float v = a11*ELEM(-2,-2) + a12*ELEM(-2,-1) + a13*ELEM(-2,0) + a12*ELEM(-2,1) + a11*ELEM(-2,2) + \
                                                            a21*ELEM(-1,-2) + a22*ELEM(-1,-1) + a23*ELEM(-1,0) + a22*ELEM(-1,1) + a21*ELEM(-1,2) + \
                                                            a31*ELEM(0,-2) + a32*ELEM(0,-1) + a33*ELEM(0,0) + a32*ELEM(0,1) + a31*ELEM(0,2) + \
                                                            a21*ELEM(1,-2) + a22*ELEM(1,-1) + a23*ELEM(1,0) + a22*ELEM(1,1) + a21*ELEM(1,2) + \
                                                            a11*ELEM(2,-2) + a12*ELEM(2,-1) + a13*ELEM(2,0) + a12*ELEM(2,1) + a11*ELEM(2,2); \
                                                       v /= a11*SULY(-2,-2) + a12*SULY(-2,-1) + a13*SULY(-2,0) + a12*SULY(-2,1) + a11*SULY(-2,2) + \
                                                            a21*SULY(-1,-2) + a22*SULY(-1,-1) + a23*SULY(-1,0) + a22*SULY(-1,1) + a21*SULY(-1,2) + \
                                                            a31*SULY(0,-2) + a32*SULY(0,-1) + a33*SULY(0,0) + a32*SULY(0,1) + a31*SULY(0,2) + \
                                                            a21*SULY(1,-2) + a22*SULY(1,-1) + a23*SULY(1,0) + a22*SULY(1,1) + a21*SULY(1,2) + \
                                                            a11*SULY(2,-2) + a12*SULY(2,-1) + a13*SULY(2,0) + a12*SULY(2,1) + a11*SULY(2,2);

#define BL_OPER7(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44) \
                                                      for (int i=rstart; i<rend; i++) { \
                                                    	  for (int j=cstart; j<cend; j++) { \
													   float v = a11*ELEM(-3,-3) + a12*ELEM(-3,-2) + a13*ELEM(-3,-1) + a14*ELEM(-3,0) + a13*ELEM(-3,1) + a12*ELEM(-3,2) + a11*ELEM(-3,3) + \
                                                            a21*ELEM(-2,-3) + a22*ELEM(-2,-2) + a23*ELEM(-2,-1) + a24*ELEM(-2,0) + a23*ELEM(-2,1) + a22*ELEM(-2,2) + a21*ELEM(-2,3) + \
                                                            a31*ELEM(-1,-3) + a32*ELEM(-1,-2) + a33*ELEM(-1,-1) + a34*ELEM(-1,0) + a33*ELEM(-1,1) + a32*ELEM(-1,2) + a31*ELEM(-1,3) + \
                                                            a41*ELEM(0,-3) + a42*ELEM(0,-2) + a43*ELEM(0,-1) + a44*ELEM(0,0) + a43*ELEM(0,1) + a42*ELEM(0,2) + a41*ELEM(0,3) + \
                                                            a31*ELEM(1,-3) + a32*ELEM(1,-2) + a33*ELEM(1,-1) + a34*ELEM(1,0) + a33*ELEM(1,1) + a32*ELEM(1,2) + a31*ELEM(1,3) + \
                                                            a21*ELEM(2,-3) + a22*ELEM(2,-2) + a23*ELEM(2,-1) + a24*ELEM(2,0) + a23*ELEM(2,1) + a22*ELEM(2,2) + a21*ELEM(2,3) + \
                                                            a11*ELEM(3,-3) + a12*ELEM(3,-2) + a13*ELEM(3,-1) + a14*ELEM(3,0) + a13*ELEM(3,1) + a12*ELEM(3,2) + a11*ELEM(3,3); \
                                                       v /= a11*SULY(-3,-3) + a12*SULY(-3,-2) + a13*SULY(-3,-1) + a14*SULY(-3,0) + a13*SULY(-3,1) + a12*SULY(-3,2) + a11*SULY(-3,3) + \
                                                            a21*SULY(-2,-3) + a22*SULY(-2,-2) + a23*SULY(-2,-1) + a24*SULY(-2,0) + a23*SULY(-2,1) + a22*SULY(-2,2) + a21*SULY(-2,3) + \
                                                            a31*SULY(-1,-3) + a32*SULY(-1,-2) + a33*SULY(-1,-1) + a34*SULY(-1,0) + a33*SULY(-1,1) + a32*SULY(-1,2) + a31*SULY(-1,3) + \
                                                            a41*SULY(0,-3) + a42*SULY(0,-2) + a43*SULY(0,-1) + a44*SULY(0,0) + a43*SULY(0,1) + a42*SULY(0,2) + a41*SULY(0,3) + \
                                                            a31*SULY(1,-3) + a32*SULY(1,-2) + a33*SULY(1,-1) + a34*SULY(1,0) + a33*SULY(1,1) + a32*SULY(1,2) + a31*SULY(1,3) + \
                                                            a21*SULY(2,-3) + a22*SULY(2,-2) + a23*SULY(2,-1) + a24*SULY(2,0) + a23*SULY(2,1) + a22*SULY(2,2) + a21*SULY(2,3) + \
                                                            a11*SULY(3,-3) + a12*SULY(3,-2) + a13*SULY(3,-1) + a14*SULY(3,0) + a13*SULY(3,1) + a12*SULY(3,2) + a11*SULY(3,3);


#define BL_OPER9(a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55) \
                                                      for (int i=rstart; i<rend; i++) { \
                                                    	  for (int j=cstart; j<cend; j++) { \
                                                      float v = a11*ELEM(-4,-4) + a12*ELEM(-4,-3) + a13*ELEM(-4,-2) + a14*ELEM(-4,-1) + a15*ELEM(-4,0) + a14*ELEM(-4,1) + a13*ELEM(-4,2) + a12*ELEM(-4,3) + a11*ELEM(-4,4) + \
                                                            a21*ELEM(-3,-4) + a22*ELEM(-3,-3) + a23*ELEM(-3,-2) + a24*ELEM(-3,-1) + a25*ELEM(-3,0) + a24*ELEM(-3,1) + a23*ELEM(-3,2) + a22*ELEM(-3,3) + a21*ELEM(-3,4) + \
                                                            a31*ELEM(-2,-4) + a32*ELEM(-2,-3) + a33*ELEM(-2,-2) + a34*ELEM(-2,-1) + a35*ELEM(-2,0) + a34*ELEM(-2,1) + a33*ELEM(-2,2) + a32*ELEM(-2,3) + a31*ELEM(-2,4) + \
                                                            a41*ELEM(-1,-4) + a42*ELEM(-1,-3) + a43*ELEM(-1,-2) + a44*ELEM(-1,-1) + a45*ELEM(-1,0) + a44*ELEM(-1,1) + a43*ELEM(-1,2) + a42*ELEM(-1,3) + a41*ELEM(-1,4) + \
                                                            a51*ELEM(0,-4) + a52*ELEM(0,-3) + a53*ELEM(0,-2) + a54*ELEM(0,-1) + a55*ELEM(0,0) + a54*ELEM(0,1) + a53*ELEM(0,2) + a52*ELEM(0,3) + a51*ELEM(0,4) + \
                                                            a41*ELEM(1,-4) + a42*ELEM(1,-3) + a43*ELEM(1,-2) + a44*ELEM(1,-1) + a45*ELEM(1,0) + a44*ELEM(1,1) + a43*ELEM(1,2) + a42*ELEM(1,3) + a41*ELEM(1,4) + \
                                                            a31*ELEM(2,-4) + a32*ELEM(2,-3) + a33*ELEM(2,-2) + a34*ELEM(2,-1) + a35*ELEM(2,0) + a34*ELEM(2,1) + a33*ELEM(2,2) + a32*ELEM(2,3) + a31*ELEM(2,4) + \
                                                            a21*ELEM(3,-4) + a22*ELEM(3,-3) + a23*ELEM(3,-2) + a24*ELEM(3,-1) + a25*ELEM(3,0) + a24*ELEM(3,1) + a23*ELEM(3,2) + a22*ELEM(3,3) + a21*ELEM(3,4) + \
                                                            a11*ELEM(4,-4) + a12*ELEM(4,-3) + a13*ELEM(4,-2) + a14*ELEM(4,-1) + a15*ELEM(4,0) + a14*ELEM(4,1) + a13*ELEM(4,2) + a12*ELEM(4,3) + a11*ELEM(4,4); \
                                                       v /= a11*SULY(-4,-4) + a12*SULY(-4,-3) + a13*SULY(-4,-2) + a14*SULY(-4,-1) + a15*SULY(-4,0) + a14*SULY(-4,1) + a13*SULY(-4,2) + a12*SULY(-4,3) + a11*SULY(-4,4) + \
                                                            a21*SULY(-3,-4) + a22*SULY(-3,-3) + a23*SULY(-3,-2) + a24*SULY(-3,-1) + a25*SULY(-3,0) + a24*SULY(-3,1) + a23*SULY(-3,2) + a22*SULY(-3,3) + a21*SULY(-3,4) + \
                                                            a31*SULY(-2,-4) + a32*SULY(-2,-3) + a33*SULY(-2,-2) + a34*SULY(-2,-1) + a35*SULY(-2,0) + a34*SULY(-2,1) + a33*SULY(-2,2) + a32*SULY(-2,3) + a31*SULY(-2,4) + \
                                                            a41*SULY(-1,-4) + a42*SULY(-1,-3) + a43*SULY(-1,-2) + a44*SULY(-1,-1) + a45*SULY(-1,0) + a44*SULY(-1,1) + a43*SULY(-1,2) + a42*SULY(-1,3) + a41*SULY(-1,4) + \
                                                            a51*SULY(0,-4) + a52*SULY(0,-3) + a53*SULY(0,-2) + a54*SULY(0,-1) + a55*SULY(0,0) + a54*SULY(0,1) + a53*SULY(0,2) + a52*SULY(0,3) + a51*SULY(0,4) + \
                                                            a41*SULY(1,-4) + a42*SULY(1,-3) + a43*SULY(1,-2) + a44*SULY(1,-1) + a45*SULY(1,0) + a44*SULY(1,1) + a43*SULY(1,2) + a42*SULY(1,3) + a41*SULY(1,4) + \
                                                            a31*SULY(2,-4) + a32*SULY(2,-3) + a33*SULY(2,-2) + a34*SULY(2,-1) + a35*SULY(2,0) + a34*SULY(2,1) + a33*SULY(2,2) + a32*SULY(2,3) + a31*SULY(2,4) + \
                                                            a21*SULY(3,-4) + a22*SULY(3,-3) + a23*SULY(3,-2) + a24*SULY(3,-1) + a25*SULY(3,0) + a24*SULY(3,1) + a23*SULY(3,2) + a22*SULY(3,3) + a21*SULY(3,4) + \
                                                            a11*SULY(4,-4) + a12*SULY(4,-3) + a13*SULY(4,-2) + a14*SULY(4,-1) + a15*SULY(4,0) + a14*SULY(4,1) + a13*SULY(4,2) + a12*SULY(4,3) + a11*SULY(4,4);

#define BL_OPER11(a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,a31,a32,a33,a34,a35,a36,a41,a42,a43,a44,a45,a46,a51,a52,a53,a54,a55,a56,a61,a62,a63,a64,a65,a66) \
                                                      for (int i=rstart; i<rend; i++) { \
														for (int j=cstart; j<cend; j++) { \
                                                      float v = a11*ELEM(-5,-5) + a12*ELEM(-5,-4) + a13*ELEM(-5,-3) + a14*ELEM(-5,-2) + a15*ELEM(-5,-1) + a16*ELEM(-5,0) + a15*ELEM(-5,1) + a14*ELEM(-5,2) + a13*ELEM(-5,3) + a12*ELEM(-5,4) + a11*ELEM(-5,5) + \
                                                            a21*ELEM(-4,-5) + a22*ELEM(-4,-4) + a23*ELEM(-4,-3) + a24*ELEM(-4,-2) + a25*ELEM(-4,-1) + a26*ELEM(-4,0) + a25*ELEM(-4,1) + a24*ELEM(-4,2) + a23*ELEM(-4,3) + a22*ELEM(-4,4) + a21*ELEM(-4,5) + \
                                                            a31*ELEM(-3,-5) + a32*ELEM(-3,-4) + a33*ELEM(-3,-3) + a34*ELEM(-3,-2) + a35*ELEM(-3,-1) + a36*ELEM(-3,0) + a35*ELEM(-3,1) + a34*ELEM(-3,2) + a33*ELEM(-3,3) + a32*ELEM(-3,4) + a31*ELEM(-3,5) + \
                                                            a41*ELEM(-2,-5) + a42*ELEM(-2,-4) + a43*ELEM(-2,-3) + a44*ELEM(-2,-2) + a45*ELEM(-2,-1) + a46*ELEM(-2,0) + a45*ELEM(-2,1) + a44*ELEM(-2,2) + a43*ELEM(-2,3) + a42*ELEM(-2,4) + a41*ELEM(-4,5) + \
                                                            a51*ELEM(-1,-5) + a52*ELEM(-1,-4) + a53*ELEM(-1,-3) + a54*ELEM(-1,-2) + a55*ELEM(-1,-1) + a56*ELEM(-1,0) + a55*ELEM(-1,1) + a54*ELEM(-1,2) + a53*ELEM(-1,3) + a52*ELEM(-1,4) + a51*ELEM(-1,5) + \
                                                            a61*ELEM(0,-5) + a62*ELEM(0,-4) + a63*ELEM(0,-3) + a64*ELEM(0,-2) + a65*ELEM(0,-1) + a66*ELEM(0,0) + a65*ELEM(0,1) + a64*ELEM(0,2) + a63*ELEM(0,3) + a62*ELEM(0,4) + a61*ELEM(0,5) + \
                                                            a51*ELEM(1,-5) + a52*ELEM(1,-4) + a53*ELEM(1,-3) + a54*ELEM(1,-2) + a55*ELEM(1,-1) + a56*ELEM(1,0) + a55*ELEM(1,1) + a54*ELEM(1,2) + a53*ELEM(1,3) + a52*ELEM(1,4) + a51*ELEM(1,5) + \
                                                            a41*ELEM(2,-5) + a42*ELEM(2,-4) + a43*ELEM(2,-3) + a44*ELEM(2,-2) + a45*ELEM(2,-1) + a46*ELEM(2,0) + a45*ELEM(2,1) + a44*ELEM(2,2) + a43*ELEM(2,3) + a42*ELEM(2,4) + a41*ELEM(2,5) + \
                                                            a31*ELEM(3,-5) + a32*ELEM(3,-4) + a33*ELEM(3,-3) + a34*ELEM(3,-2) + a35*ELEM(3,-1) + a36*ELEM(3,0) + a35*ELEM(3,1) + a34*ELEM(3,2) + a33*ELEM(3,3) + a32*ELEM(3,4) + a31*ELEM(3,5) + \
                                                            a21*ELEM(4,-5) + a22*ELEM(4,-4) + a23*ELEM(4,-3) + a24*ELEM(4,-2) + a25*ELEM(4,-1) + a26*ELEM(4,0) + a25*ELEM(4,1) + a24*ELEM(4,2) + a23*ELEM(4,3) + a22*ELEM(4,4) + a21*ELEM(4,5) + \
                                                            a11*ELEM(5,-5) + a12*ELEM(5,-4) + a13*ELEM(5,-3) + a14*ELEM(5,-2) + a15*ELEM(5,-1) + a16*ELEM(5,0) + a15*ELEM(5,1) + a14*ELEM(5,2) + a13*ELEM(5,3) + a12*ELEM(5,4) + a11*ELEM(5,5); \
                                                       v /= a11*SULY(-5,-5) + a12*SULY(-5,-4) + a13*SULY(-5,-3) + a14*SULY(-5,-2) + a15*SULY(-5,-1) + a16*SULY(-5,0) + a15*SULY(-5,1) + a14*SULY(-5,2) + a13*SULY(-5,3) + a12*SULY(-5,4) + a11*SULY(-5,5) + \
                                                            a21*SULY(-4,-5) + a22*SULY(-4,-4) + a23*SULY(-4,-3) + a24*SULY(-4,-2) + a25*SULY(-4,-1) + a26*SULY(-4,0) + a25*SULY(-4,1) + a24*SULY(-4,2) + a23*SULY(-4,3) + a22*SULY(-4,4) + a21*SULY(-4,5) + \
                                                            a31*SULY(-3,-5) + a32*SULY(-3,-4) + a33*SULY(-3,-3) + a34*SULY(-3,-2) + a35*SULY(-3,-1) + a36*SULY(-3,0) + a35*SULY(-3,1) + a34*SULY(-3,2) + a33*SULY(-3,3) + a32*SULY(-3,4) + a31*SULY(-3,5) + \
                                                            a41*SULY(-2,-5) + a42*SULY(-2,-4) + a43*SULY(-2,-3) + a44*SULY(-2,-2) + a45*SULY(-2,-1) + a46*SULY(-2,0) + a45*SULY(-2,1) + a44*SULY(-2,2) + a43*SULY(-2,3) + a42*SULY(-2,4) + a41*SULY(-4,5) + \
                                                            a51*SULY(-1,-5) + a52*SULY(-1,-4) + a53*SULY(-1,-3) + a54*SULY(-1,-2) + a55*SULY(-1,-1) + a56*SULY(-1,0) + a55*SULY(-1,1) + a54*SULY(-1,2) + a53*SULY(-1,3) + a52*SULY(-1,4) + a51*SULY(-1,5) + \
                                                            a61*SULY(0,-5) + a62*SULY(0,-4) + a63*SULY(0,-3) + a64*SULY(0,-2) + a65*SULY(0,-1) + a66*SULY(0,0) + a65*SULY(0,1) + a64*SULY(0,2) + a63*SULY(0,3) + a62*SULY(0,4) + a61*SULY(0,5) + \
                                                            a51*SULY(1,-5) + a52*SULY(1,-4) + a53*SULY(1,-3) + a54*SULY(1,-2) + a55*SULY(1,-1) + a56*SULY(1,0) + a55*SULY(1,1) + a54*SULY(1,2) + a53*SULY(1,3) + a52*SULY(1,4) + a51*SULY(1,5) + \
                                                            a41*SULY(2,-5) + a42*SULY(2,-4) + a43*SULY(2,-3) + a44*SULY(2,-2) + a45*SULY(2,-1) + a46*SULY(2,0) + a45*SULY(2,1) + a44*SULY(2,2) + a43*SULY(2,3) + a42*SULY(2,4) + a41*SULY(2,5) + \
                                                            a31*SULY(3,-5) + a32*SULY(3,-4) + a33*SULY(3,-3) + a34*SULY(3,-2) + a35*SULY(3,-1) + a36*SULY(3,0) + a35*SULY(3,1) + a34*SULY(3,2) + a33*SULY(3,3) + a32*SULY(3,4) + a31*SULY(3,5) + \
                                                            a21*SULY(4,-5) + a22*SULY(4,-4) + a23*SULY(4,-3) + a24*SULY(4,-2) + a25*SULY(4,-1) + a26*SULY(4,0) + a25*SULY(4,1) + a24*SULY(4,2) + a23*SULY(4,3) + a22*SULY(4,4) + a21*SULY(4,5) + \
                                                            a11*SULY(5,-5) + a12*SULY(5,-4) + a13*SULY(5,-3) + a14*SULY(5,-2) + a15*SULY(5,-1) + a16*SULY(5,0) + a15*SULY(5,1) + a14*SULY(5,2) + a13*SULY(5,3) + a12*SULY(5,4) + a11*SULY(5,5); \

namespace rtengine {

// sigma = 0.5
void bilateral05 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(318,1)
	#pragma omp parallel for if (multiThread)
    BL_OPER3(1,7,7,55)
    BL_END(1)
}

// sigma = 0.6
void bilateral06 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(768,1)
	#pragma omp parallel for if (multiThread)
    BL_OPER3(1,4,4,16)
    BL_END(1)
}

// sigma = 0.7
void bilateral07 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(366,2)
	#pragma omp parallel for if (multiThread)
    BL_OPER5(0,0,1,0,8,21,1,21,59)
    BL_END(2)
}

// sigma = 0.8
void bilateral08 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(753,2)
	#pragma omp parallel for if (multiThread)
    BL_OPER5(0,0,1,0,5,10,1,10,23)
    BL_END(2)
}

// sigma = 0.9
void bilateral09 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(595,2)
	#pragma omp parallel for if (multiThread)
    BL_OPER5(0,1,2,1,6,12,2,12,22)
    BL_END(2)
}

// sigma = 1.0
void bilateral10 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(910,2)
	#pragma omp parallel for if (multiThread)
    BL_OPER5(0,1,2,1,4,7,2,7,12)
    BL_END(2)
}

// sigma = 1.1
void bilateral11 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(209,3)
	#pragma omp parallel for if (multiThread)
    BL_OPER7(0,0,1,1,0,2,5,8,1,5,18,27,1,8,27,41)
    BL_END(3)
}

// sigma = 1.2
void bilateral12 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(322,3)
	#pragma omp parallel for if (multiThread)
    BL_OPER7(0,0,1,1,0,1,4,6,1,4,11,16,1,6,16,23)
    BL_END(3)
}

// sigma = 1.3
void bilateral13 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(336,3)
	#pragma omp parallel for if (multiThread)
    BL_OPER7(0,0,1,1,0,2,4,6,1,4,11,14,1,6,14,19)
    BL_END(3)
}

// sigma = 1.4
void bilateral14 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(195,3)
	#pragma omp parallel for if (multiThread)
    BL_OPER7(0,1,2,3,1,4,8,10,2,8,17,21,3,10,21,28)
    BL_END(3)
}

// sigma = 1.5
void bilateral15 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(132,4)
	#pragma omp parallel for if (multiThread)
    BL_OPER9(0,0,0,1,1,0,1,2,4,5,0,2,6,12,14,1,4,12,22,28,1,5,14,28,35)
    BL_END(4)
}

// sigma = 1.6
void bilateral16 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(180,4)
	#pragma omp parallel for if (multiThread)
    BL_OPER9(0,0,0,1,1,0,1,2,3,4,0,2,5,9,10,1,3,9,15,19,1,4,10,19,23)
    BL_END(4)
}

// sigma = 1.7
void bilateral17 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(195,4)
	#pragma omp parallel for if (multiThread)
    BL_OPER9(0,0,1,1,1,0,1,2,3,4,1,2,5,8,9,1,3,8,13,16,1,4,9,16,19)
    BL_END(4)
}

// sigma = 1.8
void bilateral18 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(151,4)
	#pragma omp parallel for if (multiThread)
    BL_OPER9(0,0,1,2,2,0,1,3,5,5,1,3,6,10,12,2,5,10,16,19,2,5,12,19,22)
    BL_END(4)
}

// sigma = 1.9
void bilateral19 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(151,4)
	#pragma omp parallel for if (multiThread)
    BL_OPER9(0,0,1,2,2,0,1,3,4,5,1,3,5,8,9,2,4,8,12,14,2,5,9,14,16)
    BL_END(4)
}

// sigma = 2
void bilateral20 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(116,5)
	#pragma omp parallel for if (multiThread)
    BL_OPER11(0,0,0,1,1,1,0,0,1,2,3,3,0,1,2,4,7,7,1,2,4,8,12,14,1,3,7,12,18,20,1,3,7,14,20,23)
    BL_END(5)
}

// sigma = 2.1
void bilateral21 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(127,5)
	#pragma omp parallel for if (multiThread)
    BL_OPER11(0,0,0,1,1,1,0,0,1,2,3,3,0,1,2,4,6,7,1,2,4,8,11,12,1,3,6,11,15,17,1,3,7,12,17,19)
    BL_END(5)
}

// sigma = 2.2
void bilateral22 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(109,5)
	#pragma omp parallel for if (multiThread)
    BL_OPER11(0,0,0,1,1,2,0,1,2,3,3,4,1,2,3,5,7,8,1,3,5,9,12,13,1,3,7,12,16,18,2,4,8,13,18,20)
    BL_END(5)
}

// sigma = 2.3
void bilateral23 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(132,5)
	#pragma omp parallel for if (multiThread)
    BL_OPER11(0,0,1,1,1,1,0,1,1,2,3,3,1,1,3,5,6,7,1,2,5,7,10,11,1,3,6,10,13,14,1,3,7,11,14,16)
    BL_END(5)
}

// sigma = 2.4
void bilateral24 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

    BL_BEGIN(156,5)
	#pragma omp parallel for if (multiThread)
    BL_OPER11(0,0,1,1,1,1,0,1,1,2,3,3,1,1,3,4,5,6,1,2,4,6,8,9,1,3,5,8,10,11,1,3,6,9,11,12)
    BL_END(5)
}

// sigma = 2.5
void bilateral25 (float** src, float** dst, float** buffer, int W, int H, float sens, bool multiThread) {

	BL_BEGIN(173,5)
	#pragma omp parallel for if (multiThread)
    BL_OPER11(0,0,1,1,1,1,0,1,1,2,3,3,1,1,2,4,5,5,1,2,4,5,7,7,1,3,5,7,9,9,1,3,5,7,9,10)
    BL_END(5)
}

// main bilateral filter
void bilateral (float** src, float** dst, float** buffer, int W, int H, float sigma, float sens, bool multiThread) {

    if (sigma<0.45)
		#pragma omp parallel for if (multiThread)
    	for (int i=0; i<H; i++) {
            memcpy (buffer[i], src[i], W*sizeof(float));
            memcpy (dst[i], buffer[i], W*sizeof(float));
        }
    else if (sigma<0.55)
        bilateral05 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.65)
        bilateral06 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.75)
        bilateral07 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.85)
        bilateral08 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.95)
        bilateral09 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.05)
        bilateral10 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.15)
        bilateral11 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.25)
        bilateral12 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.35)
        bilateral13 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.45)
        bilateral14 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.55)
        bilateral15 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.65)
        bilateral16 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.75)
        bilateral17 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.85)
        bilateral18 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.95)
        bilateral19 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.05)
        bilateral20 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.15)
        bilateral21 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.25)
        bilateral22 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.35)
        bilateral23 (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.45)
        bilateral24 (src, dst, buffer, W, H, sens, multiThread);
    else
        bilateral25 (src, dst, buffer, W, H, sens, multiThread);
}

}
