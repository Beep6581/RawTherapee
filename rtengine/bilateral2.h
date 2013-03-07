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
#ifndef _BILATERAL2_
#define _BILATERAL2_

#include <cmath>
#include <cstring>
#include <cstdio>
#include <glibmm.h>

#include "rtengine.h"
#include "rt_math.h"
#include "alignedbuffer.h"
#include "mytime.h"
#include "gauss.h"

#include "array2D.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace rtengine;

// This seems ugly, but way faster than any other solutions I tried
#define ELEM(a,b) (src[i - a][j - b] * ec[src[i - a][j - b]-src[i][j]+65536.0f])
#define SULY(a,b) (ec[src[i - a][j - b]-src[i][j]+65536.0f])
//#define SULY(a,b) (ec[((int)(src[i - a][j - b]-src[i][j]+0x10000))])
#define BL_BEGIN(a,b)   double scale = (a); \
                        LUTf ec (0x20000); \
                        for (int i=0; i<0x20000; i++) \
                            ec[i] = (exp(-(double)(i-0x10000)*(double)(i-0x10000) / (2.0*sens*sens))*scale); \
                        int rstart = b; \
                        int rend = H-b; \
                        int cstart = b; \
                        int cend = W-b;

#define BL_FREE 		buffer[i][j] = v; }};
#define BL_END(b)       for (int i=0; i<H; i++)  \
                            for (int j=0; j<W; j++)  \
                                if (i<rstart || j<cstart || i>=rend || j>=cend) \
                                    dst[i][j] = src[i][j]; \
                                else \
                                    dst[i][j] = buffer[i][j];

#define BL_OPER3(a11,a12,a21,a22) for (int i=rstart; i<rend; i++) { \
									for (int j=cstart; j<cend; j++) { \
								   A v = a11*ELEM(-1,-1) + a12*ELEM(-1,0) + a11*ELEM(-1,1) + \
                                        a21*ELEM(0,-1) + a22*ELEM(0,0) + a21*ELEM(0,1) + \
                                        a11*ELEM(1,-1) + a12*ELEM(1,0) + a11*ELEM(1,1); \
                                   v /= a11*SULY(-1,-1) + a12*SULY(-1,0) + a11*SULY(-1,1) + \
                                        a21*SULY(0,-1) + a22*SULY(0,0) + a21*SULY(0,1) + \
                                        a11*SULY(1,-1) + a12*SULY(1,0) + a11*SULY(1,1); 


#define BL_OPER5(a11,a12,a13,a21,a22,a23,a31,a32,a33) for (int i=rstart; i<rend; i++) { \
														for (int j=cstart; j<cend; j++) { \
													   A v = a11*ELEM(-2,-2) + a12*ELEM(-2,-1) + a13*ELEM(-2,0) + a12*ELEM(-2,1) + a11*ELEM(-2,2) + \
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
													   A v = a11*ELEM(-3,-3) + a12*ELEM(-3,-2) + a13*ELEM(-3,-1) + a14*ELEM(-3,0) + a13*ELEM(-3,1) + a12*ELEM(-3,2) + a11*ELEM(-3,3) + \
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
                                                      A v = a11*ELEM(-4,-4) + a12*ELEM(-4,-3) + a13*ELEM(-4,-2) + a14*ELEM(-4,-1) + a15*ELEM(-4,0) + a14*ELEM(-4,1) + a13*ELEM(-4,2) + a12*ELEM(-4,3) + a11*ELEM(-4,4) + \
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
                                                      A v = a11*ELEM(-5,-5) + a12*ELEM(-5,-4) + a13*ELEM(-5,-3) + a14*ELEM(-5,-2) + a15*ELEM(-5,-1) + a16*ELEM(-5,0) + a15*ELEM(-5,1) + a14*ELEM(-5,2) + a13*ELEM(-5,3) + a12*ELEM(-5,4) + a11*ELEM(-5,5) + \
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


// sigma = 0.5
template<class T, class A> void bilateral05 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
    
    BL_BEGIN(318,1)
	#pragma omp for
    BL_OPER3(1,7,7,55)
    BL_FREE
#pragma omp for
    BL_END(1)
}

// sigma = 0.6
template<class T, class A> void bilateral06 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
    
    BL_BEGIN(768,1)
	#pragma omp for
    BL_OPER3(1,4,4,16)
    BL_FREE
#pragma omp for
    BL_END(1)
}

// sigma = 0.7
template<class T, class A> void bilateral07 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
    
    BL_BEGIN(366,2)
	#pragma omp for
    BL_OPER5(0,0,1,0,8,21,1,21,59)
    BL_FREE
#pragma omp for
    BL_END(2)
}

// sigma = 0.8
template<class T, class A> void bilateral08 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
    
    BL_BEGIN(753,2)
	#pragma omp for
    BL_OPER5(0,0,1,0,5,10,1,10,23)
    BL_FREE
#pragma omp for
    BL_END(2)
}

// sigma = 0.9
template<class T, class A> void bilateral09 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(595,2)
	#pragma omp for
    BL_OPER5(0,1,2,1,6,12,2,12,22)
    BL_FREE
#pragma omp for
    BL_END(2)
}

// sigma = 1.0
template<class T, class A> void bilateral10 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {

    BL_BEGIN(910,2)
	#pragma omp for
    BL_OPER5(0,1,2,1,4,7,2,7,12)
    BL_FREE
#pragma omp for
    BL_END(2)
}

// sigma = 1.1
template<class T, class A> void bilateral11 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(209,3)
	#pragma omp for
    BL_OPER7(0,0,1,1,0,2,5,8,1,5,18,27,1,8,27,41)
    BL_FREE
#pragma omp for
    BL_END(3)
}

// sigma = 1.2
template<class T, class A> void bilateral12 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(322,3)
	#pragma omp for
    BL_OPER7(0,0,1,1,0,1,4,6,1,4,11,16,1,6,16,23)
    BL_FREE
#pragma omp for
    BL_END(3)
}

// sigma = 1.3
template<class T, class A> void bilateral13 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(336,3)
	#pragma omp for
    BL_OPER7(0,0,1,1,0,2,4,6,1,4,11,14,1,6,14,19)
    BL_FREE
#pragma omp for
    BL_END(3)
}

// sigma = 1.4
template<class T, class A> void bilateral14 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(195,3)
	#pragma omp for
    BL_OPER7(0,1,2,3,1,4,8,10,2,8,17,21,3,10,21,28)
    BL_FREE
#pragma omp for
    BL_END(3)   
}

// sigma = 1.5
template<class T, class A> void bilateral15 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(132,4)
	#pragma omp for
    BL_OPER9(0,0,0,1,1,0,1,2,4,5,0,2,6,12,14,1,4,12,22,28,1,5,14,28,35)
    BL_FREE
#pragma omp for
    BL_END(4)   
}

// sigma = 1.6
template<class T, class A> void bilateral16 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(180,4)
	#pragma omp for
    BL_OPER9(0,0,0,1,1,0,1,2,3,4,0,2,5,9,10,1,3,9,15,19,1,4,10,19,23)
    BL_FREE
#pragma omp for
    BL_END(4)   
}

// sigma = 1.7
template<class T, class A> void bilateral17 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(195,4)
	#pragma omp for
    BL_OPER9(0,0,1,1,1,0,1,2,3,4,1,2,5,8,9,1,3,8,13,16,1,4,9,16,19)
    BL_FREE
#pragma omp for
    BL_END(4)   
}

// sigma = 1.8
template<class T, class A> void bilateral18 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(151,4)
	#pragma omp for
    BL_OPER9(0,0,1,2,2,0,1,3,5,5,1,3,6,10,12,2,5,10,16,19,2,5,12,19,22)
    BL_FREE
#pragma omp for
    BL_END(4)   
}

// sigma = 1.9
template<class T, class A> void bilateral19 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {
   
    BL_BEGIN(151,4)
	#pragma omp for
    BL_OPER9(0,0,1,2,2,0,1,3,4,5,1,3,5,8,9,2,4,8,12,14,2,5,9,14,16)
    BL_FREE
#pragma omp for
    BL_END(4)   
}

// sigma = 2
template<class T, class A> void bilateral20 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {

    BL_BEGIN(116,5)
	#pragma omp for
    BL_OPER11(0,0,0,1,1,1,0,0,1,2,3,3,0,1,2,4,7,7,1,2,4,8,12,14,1,3,7,12,18,20,1,3,7,14,20,23)
    BL_FREE
#pragma omp for
    BL_END(5)   
}

// sigma = 2.1
template<class T, class A> void bilateral21 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {

    BL_BEGIN(127,5)
	#pragma omp for
    BL_OPER11(0,0,0,1,1,1,0,0,1,2,3,3,0,1,2,4,6,7,1,2,4,8,11,12,1,3,6,11,15,17,1,3,7,12,17,19)
    BL_FREE
#pragma omp for
    BL_END(5)
}

// sigma = 2.2
template<class T, class A> void bilateral22 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {

    BL_BEGIN(109,5)
	#pragma omp for
    BL_OPER11(0,0,0,1,1,2,0,1,2,3,3,4,1,2,3,5,7,8,1,3,5,9,12,13,1,3,7,12,16,18,2,4,8,13,18,20)
    BL_FREE
#pragma omp for
    BL_END(5)
}

// sigma = 2.3
template<class T, class A> void bilateral23 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {

    BL_BEGIN(132,5)
	#pragma omp for
    BL_OPER11(0,0,1,1,1,1,0,1,1,2,3,3,1,1,3,5,6,7,1,2,5,7,10,11,1,3,6,10,13,14,1,3,7,11,14,16)
    BL_FREE
#pragma omp for
    BL_END(5)
}

// sigma = 2.4
template<class T, class A> void bilateral24 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {

    BL_BEGIN(156,5)
	#pragma omp for
    BL_OPER11(0,0,1,1,1,1,0,1,1,2,3,3,1,1,3,4,5,6,1,2,4,6,8,9,1,3,5,8,10,11,1,3,6,9,11,12)
    BL_FREE
#pragma omp for
    BL_END(5)
}

// sigma = 2.5
template<class T, class A> void bilateral25 (T** src, T** dst, T** buffer, int W, int H, double sens, bool multiThread) {

	BL_BEGIN(173,5)
	#pragma omp for
    BL_OPER11(0,0,1,1,1,1,0,1,1,2,3,3,1,1,2,4,5,5,1,2,4,5,7,7,1,3,5,7,9,9,1,3,5,7,9,10)
    BL_FREE
#pragma omp for
    BL_END(5)
}

// main bilateral filter
template<class T, class A> void bilateral (T** src, T** dst, T** buffer, int W, int H, double sigma, double sens, bool multiThread) {
//parallel if (multiThread)
    if (sigma<0.45)
#ifdef _OPENMP
#pragma omp  for
#endif
    	for (int i=0; i<H; i++) {
            memcpy (buffer[i], src[i], W*sizeof(T));
            memcpy (dst[i], buffer[i], W*sizeof(T));
        }
    else if (sigma<0.55)
        bilateral05<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.65)
        bilateral06<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.75)
        bilateral07<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.85)
        bilateral08<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<0.95)
        bilateral09<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.05)
        bilateral10<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.15)
        bilateral11<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.25)
        bilateral12<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.35)
        bilateral13<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.45)
        bilateral14<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.55)
        bilateral15<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.65)
        bilateral16<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.75)
        bilateral17<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.85)
        bilateral18<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<1.95)
        bilateral19<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.05)
        bilateral20<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.15)
        bilateral21<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.25)
        bilateral22<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.35)
        bilateral23<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else if (sigma<2.45)
        bilateral24<T,A> (src, dst, buffer, W, H, sens, multiThread);
    else 
        bilateral25<T,A> (src, dst, buffer, W, H, sens, multiThread);
}

// START OF EXPERIMENTAL CODE: O(1) bilateral box filter
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define BINBIT 7 //bit depth of histogram -- there are 2^BINBIT discrete levels
#define TRANSBIT 9 //bit shift = 16-BINBIT taking short ints to bit depth BINBIT

template<class T> void bilateral (T** src, T** dst, int W, int H, int sigmar, double sigmas, int row_from, int row_to) {

    // range weights
    /*LUTf ec(0x20000);
    for (int i=0; i<0x20000; i++) 
        ec[i] = exp(-(double)(i-0x10000)*(double)(i-0x10000) / (2.0*sigmar*sigmar)); */

    // histogram
    LUTu hist (1<<BINBIT);
    LUTu rhist(1<<BINBIT);

    // buffer for the final image
    /*float** buff_final = new float*[H];
    float * real_buff_final = new float[W*H];
    for (int i=0; i<H; i++) {
        buff_final[i] = real_buff_final + i*W;
        memset (buff_final[i], 0, W*sizeof(float));
    }*/
	
	array2D<float> buff_final(W,H,ARRAY2D_CLEAR_DATA);
    
    int r = sigmas;
   
    // calculate histogram at the beginning of the row
    rhist.clear();
    for (int x = MAX(0,row_from-r); x<=MIN(H,row_from+r); x++)
        for (int y = 0; y<r+1; y++)
            rhist[((int)src[x][y])>>TRANSBIT]++;

    sigmar*=2;

    for (int i=row_from; i<row_to; i++) {
        
        // calculate histogram at the beginning of the row
        if (i>r)
            for (int x = 0; x<=MIN(H,r); x++) 
                rhist[((int)src[i-r-1][x])>>TRANSBIT]--;
        if (i<H-r)
            for (int x = 0; x<=MIN(H,r); x++) 
                rhist[((int)src[i+r][x])>>TRANSBIT]++;
            
        hist=rhist;
        for (int j=0; j<W; j++) {
        
            // subtract pixels at the left and add pixels at the right
            if (j>r)
                for (int x=MAX(0,i-r); x<=MIN(i+r,H-1); x++) 
                    hist[(int)(src[x][j-r-1])>>TRANSBIT]--;
            if (j<W-r)
                for (int x=MAX(0,i-r); x<=MIN(i+r,H-1); x++) 
                    hist[((int)src[x][j+r])>>TRANSBIT]++;

            // calculate pixel value
            float totwt = 0.0, weight;
            for (int k=0; k<=(sigmar>>TRANSBIT); k++) {               
                float w = 1.0 - (double)k/(sigmar>>TRANSBIT);
                int v = (int)(src[i][j])>>TRANSBIT;
				//float frac = ((float)(src[i][j]-(v<<TRANSBIT)))/(1<<TRANSBIT);
                if (v-k >= 0) {
                    weight = hist [v-k/*+frac*/] * w;
					totwt += weight;
                    buff_final[i][j] += weight * (src[i][j]-(k<<TRANSBIT));
                }
                if (v+k < (1<<BINBIT)) {
                    weight = hist [v+k/*+frac*/] * w;
					totwt += weight;
                    buff_final[i][j] += weight * (src[i][j]+(k<<TRANSBIT));
                }
            }
            buff_final[i][j] /= totwt;
        }
    }
    for (int i=row_from; i<row_to; i++) 
        for (int j=0; j<W; j++) 
            dst[i][j] = (T)CLIP(buff_final[i][j]);

    //delete [] real_buff_final;
    //delete [] buff_final;
}
#undef BINBIT
#undef TRANSBIT

#endif
