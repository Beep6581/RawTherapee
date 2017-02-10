/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2017 RawTherapee development team
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR float PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bilateral.h"
#include "rtengine.h"
#include "LUT.h"

namespace
{

template<typename T>
constexpr T suly(T** source, const LUTf& ec, int i, int j, int a, int b)
{
    return ec[source[i - a][j - b] - source[i][j] + 65536.0f];
}

template<typename T>
constexpr T elem(T** source, const LUTf& ec, int i, int j, int a, int b)
{
    return source[i - a][j - b] * suly(source, ec, i, j, a, b);
}

#define BL_BEGIN(a,b)   int rstart = b; \
                        int rend = height-b; \
                        int cstart = b; \
                        int cend = width-b;

#define BL_FREE         buffer[i][j] = v; }};
#define BL_END(b)       for (int i=0; i<height; i++)  \
                            for (int j=0; j<width; j++)  \
                                if (i<rstart || j<cstart || i>=rend || j>=cend) \
                                    destination[i][j] = source[i][j]; \
                                else \
                                    destination[i][j] = buffer[i][j];

#define BL_OPER7(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44) \
                                                      for (int i=rstart; i<rend; i++) { \
                                                          for (int j=cstart; j<cend; j++) { \
                                                       float v = a11 * elem(source, ec, i, j, -3, -3) + a12 * elem(source, ec, i, j, -3, -2) + a13 * elem(source, ec, i, j, -3, -1) + a14 * elem(source, ec, i, j, -3, 0) + a13 * elem(source, ec, i, j, -3, 1) + a12 * elem(source, ec, i, j, -3, 2) + a11 * elem(source, ec, i, j, -3, 3) + \
                                                            a21 * elem(source, ec, i, j, -2, -3) + a22 * elem(source, ec, i, j, -2, -2) + a23 * elem(source, ec, i, j, -2, -1) + a24 * elem(source, ec, i, j, -2, 0) + a23 * elem(source, ec, i, j, -2, 1) + a22 * elem(source, ec, i, j, -2, 2) + a21 * elem(source, ec, i, j, -2, 3) + \
                                                            a31 * elem(source, ec, i, j, -1, -3) + a32 * elem(source, ec, i, j, -1, -2) + a33 * elem(source, ec, i, j, -1, -1) + a34 * elem(source, ec, i, j, -1, 0) + a33 * elem(source, ec, i, j, -1, 1) + a32 * elem(source, ec, i, j, -1, 2) + a31 * elem(source, ec, i, j, -1, 3) + \
                                                            a41 * elem(source, ec, i, j, 0, -3) + a42 * elem(source, ec, i, j, 0, -2) + a43 * elem(source, ec, i, j, 0, -1) + a44 * elem(source, ec, i, j, 0, 0) + a43 * elem(source, ec, i, j, 0, 1) + a42 * elem(source, ec, i, j, 0, 2) + a41 * elem(source, ec, i, j, 0, 3) + \
                                                            a31 * elem(source, ec, i, j, 1, -3) + a32 * elem(source, ec, i, j, 1, -2) + a33 * elem(source, ec, i, j, 1, -1) + a34 * elem(source, ec, i, j, 1, 0) + a33 * elem(source, ec, i, j, 1, 1) + a32 * elem(source, ec, i, j, 1, 2) + a31 * elem(source, ec, i, j, 1, 3) + \
                                                            a21 * elem(source, ec, i, j, 2, -3) + a22 * elem(source, ec, i, j, 2, -2) + a23 * elem(source, ec, i, j, 2, -1) + a24 * elem(source, ec, i, j, 2, 0) + a23 * elem(source, ec, i, j, 2, 1) + a22 * elem(source, ec, i, j, 2, 2) + a21 * elem(source, ec, i, j, 2, 3) + \
                                                            a11 * elem(source, ec, i, j, 3, -3) + a12 * elem(source, ec, i, j, 3, -2) + a13 * elem(source, ec, i, j, 3, -1) + a14 * elem(source, ec, i, j, 3, 0) + a13 * elem(source, ec, i, j, 3, 1) + a12 * elem(source, ec, i, j, 3, 2) + a11 * elem(source, ec, i, j, 3, 3); \
                                                       v /= a11 * suly(source, ec, i, j, -3, -3) + a12 * suly(source, ec, i, j, -3, -2) + a13 * suly(source, ec, i, j, -3, -1) + a14 * suly(source, ec, i, j, -3, 0) + a13 * suly(source, ec, i, j, -3, 1) + a12 * suly(source, ec, i, j, -3, 2) + a11 * suly(source, ec, i, j, -3, 3) + \
                                                            a21 * suly(source, ec, i, j, -2, -3) + a22 * suly(source, ec, i, j, -2, -2) + a23 * suly(source, ec, i, j, -2, -1) + a24 * suly(source, ec, i, j, -2, 0) + a23 * suly(source, ec, i, j, -2, 1) + a22 * suly(source, ec, i, j, -2, 2) + a21 * suly(source, ec, i, j, -2, 3) + \
                                                            a31 * suly(source, ec, i, j, -1, -3) + a32 * suly(source, ec, i, j, -1, -2) + a33 * suly(source, ec, i, j, -1, -1) + a34 * suly(source, ec, i, j, -1, 0) + a33 * suly(source, ec, i, j, -1, 1) + a32 * suly(source, ec, i, j, -1, 2) + a31 * suly(source, ec, i, j, -1, 3) + \
                                                            a41 * suly(source, ec, i, j, 0, -3) + a42 * suly(source, ec, i, j, 0, -2) + a43 * suly(source, ec, i, j, 0, -1) + a44 * suly(source, ec, i, j, 0, 0) + a43 * suly(source, ec, i, j, 0, 1) + a42 * suly(source, ec, i, j, 0, 2) + a41 * suly(source, ec, i, j, 0, 3) + \
                                                            a31 * suly(source, ec, i, j, 1, -3) + a32 * suly(source, ec, i, j, 1, -2) + a33 * suly(source, ec, i, j, 1, -1) + a34 * suly(source, ec, i, j, 1, 0) + a33 * suly(source, ec, i, j, 1, 1) + a32 * suly(source, ec, i, j, 1, 2) + a31 * suly(source, ec, i, j, 1, 3) + \
                                                            a21 * suly(source, ec, i, j, 2, -3) + a22 * suly(source, ec, i, j, 2, -2) + a23 * suly(source, ec, i, j, 2, -1) + a24 * suly(source, ec, i, j, 2, 0) + a23 * suly(source, ec, i, j, 2, 1) + a22 * suly(source, ec, i, j, 2, 2) + a21 * suly(source, ec, i, j, 2, 3) + \
                                                            a11 * suly(source, ec, i, j, 3, -3) + a12 * suly(source, ec, i, j, 3, -2) + a13 * suly(source, ec, i, j, 3, -1) + a14 * suly(source, ec, i, j, 3, 0) + a13 * suly(source, ec, i, j, 3, 1) + a12 * suly(source, ec, i, j, 3, 2) + a11 * suly(source, ec, i, j, 3, 3);


#define BL_OPER9(a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55) \
                                                      for (int i=rstart; i<rend; i++) { \
                                                          for (int j=cstart; j<cend; j++) { \
                                                      float v = a11 * elem(source, ec, i, j, -4, -4) + a12 * elem(source, ec, i, j, -4, -3) + a13 * elem(source, ec, i, j, -4, -2) + a14 * elem(source, ec, i, j, -4, -1) + a15 * elem(source, ec, i, j, -4, 0) + a14 * elem(source, ec, i, j, -4, 1) + a13 * elem(source, ec, i, j, -4, 2) + a12 * elem(source, ec, i, j, -4, 3) + a11 * elem(source, ec, i, j, -4, 4) + \
                                                            a21 * elem(source, ec, i, j, -3, -4) + a22 * elem(source, ec, i, j, -3, -3) + a23 * elem(source, ec, i, j, -3, -2) + a24 * elem(source, ec, i, j, -3, -1) + a25 * elem(source, ec, i, j, -3, 0) + a24 * elem(source, ec, i, j, -3, 1) + a23 * elem(source, ec, i, j, -3, 2) + a22 * elem(source, ec, i, j, -3, 3) + a21 * elem(source, ec, i, j, -3, 4) + \
                                                            a31 * elem(source, ec, i, j, -2, -4) + a32 * elem(source, ec, i, j, -2, -3) + a33 * elem(source, ec, i, j, -2, -2) + a34 * elem(source, ec, i, j, -2, -1) + a35 * elem(source, ec, i, j, -2, 0) + a34 * elem(source, ec, i, j, -2, 1) + a33 * elem(source, ec, i, j, -2, 2) + a32 * elem(source, ec, i, j, -2, 3) + a31 * elem(source, ec, i, j, -2, 4) + \
                                                            a41 * elem(source, ec, i, j, -1, -4) + a42 * elem(source, ec, i, j, -1, -3) + a43 * elem(source, ec, i, j, -1, -2) + a44 * elem(source, ec, i, j, -1, -1) + a45 * elem(source, ec, i, j, -1, 0) + a44 * elem(source, ec, i, j, -1, 1) + a43 * elem(source, ec, i, j, -1, 2) + a42 * elem(source, ec, i, j, -1, 3) + a41 * elem(source, ec, i, j, -1, 4) + \
                                                            a51 * elem(source, ec, i, j, 0, -4) + a52 * elem(source, ec, i, j, 0, -3) + a53 * elem(source, ec, i, j, 0, -2) + a54 * elem(source, ec, i, j, 0, -1) + a55 * elem(source, ec, i, j, 0, 0) + a54 * elem(source, ec, i, j, 0, 1) + a53 * elem(source, ec, i, j, 0, 2) + a52 * elem(source, ec, i, j, 0, 3) + a51 * elem(source, ec, i, j, 0, 4) + \
                                                            a41 * elem(source, ec, i, j, 1, -4) + a42 * elem(source, ec, i, j, 1, -3) + a43 * elem(source, ec, i, j, 1, -2) + a44 * elem(source, ec, i, j, 1, -1) + a45 * elem(source, ec, i, j, 1, 0) + a44 * elem(source, ec, i, j, 1, 1) + a43 * elem(source, ec, i, j, 1, 2) + a42 * elem(source, ec, i, j, 1, 3) + a41 * elem(source, ec, i, j, 1, 4) + \
                                                            a31 * elem(source, ec, i, j, 2, -4) + a32 * elem(source, ec, i, j, 2, -3) + a33 * elem(source, ec, i, j, 2, -2) + a34 * elem(source, ec, i, j, 2, -1) + a35 * elem(source, ec, i, j, 2, 0) + a34 * elem(source, ec, i, j, 2, 1) + a33 * elem(source, ec, i, j, 2, 2) + a32 * elem(source, ec, i, j, 2, 3) + a31 * elem(source, ec, i, j, 2, 4) + \
                                                            a21 * elem(source, ec, i, j, 3, -4) + a22 * elem(source, ec, i, j, 3, -3) + a23 * elem(source, ec, i, j, 3, -2) + a24 * elem(source, ec, i, j, 3, -1) + a25 * elem(source, ec, i, j, 3, 0) + a24 * elem(source, ec, i, j, 3, 1) + a23 * elem(source, ec, i, j, 3, 2) + a22 * elem(source, ec, i, j, 3, 3) + a21 * elem(source, ec, i, j, 3, 4) + \
                                                            a11 * elem(source, ec, i, j, 4, -4) + a12 * elem(source, ec, i, j, 4, -3) + a13 * elem(source, ec, i, j, 4, -2) + a14 * elem(source, ec, i, j, 4, -1) + a15 * elem(source, ec, i, j, 4, 0) + a14 * elem(source, ec, i, j, 4, 1) + a13 * elem(source, ec, i, j, 4, 2) + a12 * elem(source, ec, i, j, 4, 3) + a11 * elem(source, ec, i, j, 4, 4); \
                                                       v /= a11 * suly(source, ec, i, j, -4, -4) + a12 * suly(source, ec, i, j, -4, -3) + a13 * suly(source, ec, i, j, -4, -2) + a14 * suly(source, ec, i, j, -4, -1) + a15 * suly(source, ec, i, j, -4, 0) + a14 * suly(source, ec, i, j, -4, 1) + a13 * suly(source, ec, i, j, -4, 2) + a12 * suly(source, ec, i, j, -4, 3) + a11 * suly(source, ec, i, j, -4, 4) + \
                                                            a21 * suly(source, ec, i, j, -3, -4) + a22 * suly(source, ec, i, j, -3, -3) + a23 * suly(source, ec, i, j, -3, -2) + a24 * suly(source, ec, i, j, -3, -1) + a25 * suly(source, ec, i, j, -3, 0) + a24 * suly(source, ec, i, j, -3, 1) + a23 * suly(source, ec, i, j, -3, 2) + a22 * suly(source, ec, i, j, -3, 3) + a21 * suly(source, ec, i, j, -3, 4) + \
                                                            a31 * suly(source, ec, i, j, -2, -4) + a32 * suly(source, ec, i, j, -2, -3) + a33 * suly(source, ec, i, j, -2, -2) + a34 * suly(source, ec, i, j, -2, -1) + a35 * suly(source, ec, i, j, -2, 0) + a34 * suly(source, ec, i, j, -2, 1) + a33 * suly(source, ec, i, j, -2, 2) + a32 * suly(source, ec, i, j, -2, 3) + a31 * suly(source, ec, i, j, -2, 4) + \
                                                            a41 * suly(source, ec, i, j, -1, -4) + a42 * suly(source, ec, i, j, -1, -3) + a43 * suly(source, ec, i, j, -1, -2) + a44 * suly(source, ec, i, j, -1, -1) + a45 * suly(source, ec, i, j, -1, 0) + a44 * suly(source, ec, i, j, -1, 1) + a43 * suly(source, ec, i, j, -1, 2) + a42 * suly(source, ec, i, j, -1, 3) + a41 * suly(source, ec, i, j, -1, 4) + \
                                                            a51 * suly(source, ec, i, j, 0, -4) + a52 * suly(source, ec, i, j, 0, -3) + a53 * suly(source, ec, i, j, 0, -2) + a54 * suly(source, ec, i, j, 0, -1) + a55 * suly(source, ec, i, j, 0, 0) + a54 * suly(source, ec, i, j, 0, 1) + a53 * suly(source, ec, i, j, 0, 2) + a52 * suly(source, ec, i, j, 0, 3) + a51 * suly(source, ec, i, j, 0, 4) + \
                                                            a41 * suly(source, ec, i, j, 1, -4) + a42 * suly(source, ec, i, j, 1, -3) + a43 * suly(source, ec, i, j, 1, -2) + a44 * suly(source, ec, i, j, 1, -1) + a45 * suly(source, ec, i, j, 1, 0) + a44 * suly(source, ec, i, j, 1, 1) + a43 * suly(source, ec, i, j, 1, 2) + a42 * suly(source, ec, i, j, 1, 3) + a41 * suly(source, ec, i, j, 1, 4) + \
                                                            a31 * suly(source, ec, i, j, 2, -4) + a32 * suly(source, ec, i, j, 2, -3) + a33 * suly(source, ec, i, j, 2, -2) + a34 * suly(source, ec, i, j, 2, -1) + a35 * suly(source, ec, i, j, 2, 0) + a34 * suly(source, ec, i, j, 2, 1) + a33 * suly(source, ec, i, j, 2, 2) + a32 * suly(source, ec, i, j, 2, 3) + a31 * suly(source, ec, i, j, 2, 4) + \
                                                            a21 * suly(source, ec, i, j, 3, -4) + a22 * suly(source, ec, i, j, 3, -3) + a23 * suly(source, ec, i, j, 3, -2) + a24 * suly(source, ec, i, j, 3, -1) + a25 * suly(source, ec, i, j, 3, 0) + a24 * suly(source, ec, i, j, 3, 1) + a23 * suly(source, ec, i, j, 3, 2) + a22 * suly(source, ec, i, j, 3, 3) + a21 * suly(source, ec, i, j, 3, 4) + \
                                                            a11 * suly(source, ec, i, j, 4, -4) + a12 * suly(source, ec, i, j, 4, -3) + a13 * suly(source, ec, i, j, 4, -2) + a14 * suly(source, ec, i, j, 4, -1) + a15 * suly(source, ec, i, j, 4, 0) + a14 * suly(source, ec, i, j, 4, 1) + a13 * suly(source, ec, i, j, 4, 2) + a12 * suly(source, ec, i, j, 4, 3) + a11 * suly(source, ec, i, j, 4, 4);

#define BL_OPER11(a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,a31,a32,a33,a34,a35,a36,a41,a42,a43,a44,a45,a46,a51,a52,a53,a54,a55,a56,a61,a62,a63,a64,a65,a66) \
                                                      for (int i=rstart; i<rend; i++) { \
                                                        for (int j=cstart; j<cend; j++) { \
                                                      float v = a11 * elem(source, ec, i, j, -5, -5) + a12 * elem(source, ec, i, j, -5, -4) + a13 * elem(source, ec, i, j, -5, -3) + a14 * elem(source, ec, i, j, -5, -2) + a15 * elem(source, ec, i, j, -5, -1) + a16 * elem(source, ec, i, j, -5, 0) + a15 * elem(source, ec, i, j, -5, 1) + a14 * elem(source, ec, i, j, -5, 2) + a13 * elem(source, ec, i, j, -5, 3) + a12 * elem(source, ec, i, j, -5, 4) + a11 * elem(source, ec, i, j, -5, 5) + \
                                                            a21 * elem(source, ec, i, j, -4, -5) + a22 * elem(source, ec, i, j, -4, -4) + a23 * elem(source, ec, i, j, -4, -3) + a24 * elem(source, ec, i, j, -4, -2) + a25 * elem(source, ec, i, j, -4, -1) + a26 * elem(source, ec, i, j, -4, 0) + a25 * elem(source, ec, i, j, -4, 1) + a24 * elem(source, ec, i, j, -4, 2) + a23 * elem(source, ec, i, j, -4, 3) + a22 * elem(source, ec, i, j, -4, 4) + a21 * elem(source, ec, i, j, -4, 5) + \
                                                            a31 * elem(source, ec, i, j, -3, -5) + a32 * elem(source, ec, i, j, -3, -4) + a33 * elem(source, ec, i, j, -3, -3) + a34 * elem(source, ec, i, j, -3, -2) + a35 * elem(source, ec, i, j, -3, -1) + a36 * elem(source, ec, i, j, -3, 0) + a35 * elem(source, ec, i, j, -3, 1) + a34 * elem(source, ec, i, j, -3, 2) + a33 * elem(source, ec, i, j, -3, 3) + a32 * elem(source, ec, i, j, -3, 4) + a31 * elem(source, ec, i, j, -3, 5) + \
                                                            a41 * elem(source, ec, i, j, -2, -5) + a42 * elem(source, ec, i, j, -2, -4) + a43 * elem(source, ec, i, j, -2, -3) + a44 * elem(source, ec, i, j, -2, -2) + a45 * elem(source, ec, i, j, -2, -1) + a46 * elem(source, ec, i, j, -2, 0) + a45 * elem(source, ec, i, j, -2, 1) + a44 * elem(source, ec, i, j, -2, 2) + a43 * elem(source, ec, i, j, -2, 3) + a42 * elem(source, ec, i, j, -2, 4) + a41 * elem(source, ec, i, j, -2, 5) + \
                                                            a51 * elem(source, ec, i, j, -1, -5) + a52 * elem(source, ec, i, j, -1, -4) + a53 * elem(source, ec, i, j, -1, -3) + a54 * elem(source, ec, i, j, -1, -2) + a55 * elem(source, ec, i, j, -1, -1) + a56 * elem(source, ec, i, j, -1, 0) + a55 * elem(source, ec, i, j, -1, 1) + a54 * elem(source, ec, i, j, -1, 2) + a53 * elem(source, ec, i, j, -1, 3) + a52 * elem(source, ec, i, j, -1, 4) + a51 * elem(source, ec, i, j, -1, 5) + \
                                                            a61 * elem(source, ec, i, j, 0, -5) + a62 * elem(source, ec, i, j, 0, -4) + a63 * elem(source, ec, i, j, 0, -3) + a64 * elem(source, ec, i, j, 0, -2) + a65 * elem(source, ec, i, j, 0, -1) + a66 * elem(source, ec, i, j, 0, 0) + a65 * elem(source, ec, i, j, 0, 1) + a64 * elem(source, ec, i, j, 0, 2) + a63 * elem(source, ec, i, j, 0, 3) + a62 * elem(source, ec, i, j, 0, 4) + a61 * elem(source, ec, i, j, 0, 5) + \
                                                            a51 * elem(source, ec, i, j, 1, -5) + a52 * elem(source, ec, i, j, 1, -4) + a53 * elem(source, ec, i, j, 1, -3) + a54 * elem(source, ec, i, j, 1, -2) + a55 * elem(source, ec, i, j, 1, -1) + a56 * elem(source, ec, i, j, 1, 0) + a55 * elem(source, ec, i, j, 1, 1) + a54 * elem(source, ec, i, j, 1, 2) + a53 * elem(source, ec, i, j, 1, 3) + a52 * elem(source, ec, i, j, 1, 4) + a51 * elem(source, ec, i, j, 1, 5) + \
                                                            a41 * elem(source, ec, i, j, 2, -5) + a42 * elem(source, ec, i, j, 2, -4) + a43 * elem(source, ec, i, j, 2, -3) + a44 * elem(source, ec, i, j, 2, -2) + a45 * elem(source, ec, i, j, 2, -1) + a46 * elem(source, ec, i, j, 2, 0) + a45 * elem(source, ec, i, j, 2, 1) + a44 * elem(source, ec, i, j, 2, 2) + a43 * elem(source, ec, i, j, 2, 3) + a42 * elem(source, ec, i, j, 2, 4) + a41 * elem(source, ec, i, j, 2, 5) + \
                                                            a31 * elem(source, ec, i, j, 3, -5) + a32 * elem(source, ec, i, j, 3, -4) + a33 * elem(source, ec, i, j, 3, -3) + a34 * elem(source, ec, i, j, 3, -2) + a35 * elem(source, ec, i, j, 3, -1) + a36 * elem(source, ec, i, j, 3, 0) + a35 * elem(source, ec, i, j, 3, 1) + a34 * elem(source, ec, i, j, 3, 2) + a33 * elem(source, ec, i, j, 3, 3) + a32 * elem(source, ec, i, j, 3, 4) + a31 * elem(source, ec, i, j, 3, 5) + \
                                                            a21 * elem(source, ec, i, j, 4, -5) + a22 * elem(source, ec, i, j, 4, -4) + a23 * elem(source, ec, i, j, 4, -3) + a24 * elem(source, ec, i, j, 4, -2) + a25 * elem(source, ec, i, j, 4, -1) + a26 * elem(source, ec, i, j, 4, 0) + a25 * elem(source, ec, i, j, 4, 1) + a24 * elem(source, ec, i, j, 4, 2) + a23 * elem(source, ec, i, j, 4, 3) + a22 * elem(source, ec, i, j, 4, 4) + a21 * elem(source, ec, i, j, 4, 5) + \
                                                            a11 * elem(source, ec, i, j, 5, -5) + a12 * elem(source, ec, i, j, 5, -4) + a13 * elem(source, ec, i, j, 5, -3) + a14 * elem(source, ec, i, j, 5, -2) + a15 * elem(source, ec, i, j, 5, -1) + a16 * elem(source, ec, i, j, 5, 0) + a15 * elem(source, ec, i, j, 5, 1) + a14 * elem(source, ec, i, j, 5, 2) + a13 * elem(source, ec, i, j, 5, 3) + a12 * elem(source, ec, i, j, 5, 4) + a11 * elem(source, ec, i, j, 5, 5); \
                                                       v /= a11 * suly(source, ec, i, j, -5, -5) + a12 * suly(source, ec, i, j, -5, -4) + a13 * suly(source, ec, i, j, -5, -3) + a14 * suly(source, ec, i, j, -5, -2) + a15 * suly(source, ec, i, j, -5, -1) + a16 * suly(source, ec, i, j, -5, 0) + a15 * suly(source, ec, i, j, -5, 1) + a14 * suly(source, ec, i, j, -5, 2) + a13 * suly(source, ec, i, j, -5, 3) + a12 * suly(source, ec, i, j, -5, 4) + a11 * suly(source, ec, i, j, -5, 5) + \
                                                            a21 * suly(source, ec, i, j, -4, -5) + a22 * suly(source, ec, i, j, -4, -4) + a23 * suly(source, ec, i, j, -4, -3) + a24 * suly(source, ec, i, j, -4, -2) + a25 * suly(source, ec, i, j, -4, -1) + a26 * suly(source, ec, i, j, -4, 0) + a25 * suly(source, ec, i, j, -4, 1) + a24 * suly(source, ec, i, j, -4, 2) + a23 * suly(source, ec, i, j, -4, 3) + a22 * suly(source, ec, i, j, -4, 4) + a21 * suly(source, ec, i, j, -4, 5) + \
                                                            a31 * suly(source, ec, i, j, -3, -5) + a32 * suly(source, ec, i, j, -3, -4) + a33 * suly(source, ec, i, j, -3, -3) + a34 * suly(source, ec, i, j, -3, -2) + a35 * suly(source, ec, i, j, -3, -1) + a36 * suly(source, ec, i, j, -3, 0) + a35 * suly(source, ec, i, j, -3, 1) + a34 * suly(source, ec, i, j, -3, 2) + a33 * suly(source, ec, i, j, -3, 3) + a32 * suly(source, ec, i, j, -3, 4) + a31 * suly(source, ec, i, j, -3, 5) + \
                                                            a41 * suly(source, ec, i, j, -2, -5) + a42 * suly(source, ec, i, j, -2, -4) + a43 * suly(source, ec, i, j, -2, -3) + a44 * suly(source, ec, i, j, -2, -2) + a45 * suly(source, ec, i, j, -2, -1) + a46 * suly(source, ec, i, j, -2, 0) + a45 * suly(source, ec, i, j, -2, 1) + a44 * suly(source, ec, i, j, -2, 2) + a43 * suly(source, ec, i, j, -2, 3) + a42 * suly(source, ec, i, j, -2, 4) + a41 * suly(source, ec, i, j, -2, 5) + \
                                                            a51 * suly(source, ec, i, j, -1, -5) + a52 * suly(source, ec, i, j, -1, -4) + a53 * suly(source, ec, i, j, -1, -3) + a54 * suly(source, ec, i, j, -1, -2) + a55 * suly(source, ec, i, j, -1, -1) + a56 * suly(source, ec, i, j, -1, 0) + a55 * suly(source, ec, i, j, -1, 1) + a54 * suly(source, ec, i, j, -1, 2) + a53 * suly(source, ec, i, j, -1, 3) + a52 * suly(source, ec, i, j, -1, 4) + a51 * suly(source, ec, i, j, -1, 5) + \
                                                            a61 * suly(source, ec, i, j, 0, -5) + a62 * suly(source, ec, i, j, 0, -4) + a63 * suly(source, ec, i, j, 0, -3) + a64 * suly(source, ec, i, j, 0, -2) + a65 * suly(source, ec, i, j, 0, -1) + a66 * suly(source, ec, i, j, 0, 0) + a65 * suly(source, ec, i, j, 0, 1) + a64 * suly(source, ec, i, j, 0, 2) + a63 * suly(source, ec, i, j, 0, 3) + a62 * suly(source, ec, i, j, 0, 4) + a61 * suly(source, ec, i, j, 0, 5) + \
                                                            a51 * suly(source, ec, i, j, 1, -5) + a52 * suly(source, ec, i, j, 1, -4) + a53 * suly(source, ec, i, j, 1, -3) + a54 * suly(source, ec, i, j, 1, -2) + a55 * suly(source, ec, i, j, 1, -1) + a56 * suly(source, ec, i, j, 1, 0) + a55 * suly(source, ec, i, j, 1, 1) + a54 * suly(source, ec, i, j, 1, 2) + a53 * suly(source, ec, i, j, 1, 3) + a52 * suly(source, ec, i, j, 1, 4) + a51 * suly(source, ec, i, j, 1, 5) + \
                                                            a41 * suly(source, ec, i, j, 2, -5) + a42 * suly(source, ec, i, j, 2, -4) + a43 * suly(source, ec, i, j, 2, -3) + a44 * suly(source, ec, i, j, 2, -2) + a45 * suly(source, ec, i, j, 2, -1) + a46 * suly(source, ec, i, j, 2, 0) + a45 * suly(source, ec, i, j, 2, 1) + a44 * suly(source, ec, i, j, 2, 2) + a43 * suly(source, ec, i, j, 2, 3) + a42 * suly(source, ec, i, j, 2, 4) + a41 * suly(source, ec, i, j, 2, 5) + \
                                                            a31 * suly(source, ec, i, j, 3, -5) + a32 * suly(source, ec, i, j, 3, -4) + a33 * suly(source, ec, i, j, 3, -3) + a34 * suly(source, ec, i, j, 3, -2) + a35 * suly(source, ec, i, j, 3, -1) + a36 * suly(source, ec, i, j, 3, 0) + a35 * suly(source, ec, i, j, 3, 1) + a34 * suly(source, ec, i, j, 3, 2) + a33 * suly(source, ec, i, j, 3, 3) + a32 * suly(source, ec, i, j, 3, 4) + a31 * suly(source, ec, i, j, 3, 5) + \
                                                            a21 * suly(source, ec, i, j, 4, -5) + a22 * suly(source, ec, i, j, 4, -4) + a23 * suly(source, ec, i, j, 4, -3) + a24 * suly(source, ec, i, j, 4, -2) + a25 * suly(source, ec, i, j, 4, -1) + a26 * suly(source, ec, i, j, 4, 0) + a25 * suly(source, ec, i, j, 4, 1) + a24 * suly(source, ec, i, j, 4, 2) + a23 * suly(source, ec, i, j, 4, 3) + a22 * suly(source, ec, i, j, 4, 4) + a21 * suly(source, ec, i, j, 4, 5) + \
                                                            a11 * suly(source, ec, i, j, 5, -5) + a12 * suly(source, ec, i, j, 5, -4) + a13 * suly(source, ec, i, j, 5, -3) + a14 * suly(source, ec, i, j, 5, -2) + a15 * suly(source, ec, i, j, 5, -1) + a16 * suly(source, ec, i, j, 5, 0) + a15 * suly(source, ec, i, j, 5, 1) + a14 * suly(source, ec, i, j, 5, 2) + a13 * suly(source, ec, i, j, 5, 3) + a12 * suly(source, ec, i, j, 5, 4) + a11 * suly(source, ec, i, j, 5, 5); \

}

class rtengine::Bilateral::Implementation final
{
public:
    Implementation(
        float** _source,
        float** _destination,
        float** _buffer,
        int _width,
        int _height,
        double _sigma,
        double _sensitivity
    ) :
        source(_source),
        destination(_destination),
        buffer(_buffer),
        width(_width),
        height(_height),
        sigma(_sigma),
        ec(0x20000)
    {
        const double sensitivity_square_times_two = _sensitivity * _sensitivity * 2.0;

        for (int i=0; i<0x20000; i++) {
            const double i_biased = i - 0x10000;
            ec[i] = std::exp(-i_biased * i_biased / sensitivity_square_times_two);
        }
    }

    void compute()
    {
        if (sigma < 0.45)
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < height; i++) {
                memcpy(buffer[i], source[i], width * sizeof(float));
                memcpy(destination[i], buffer[i], width * sizeof(float));
            }
        else if (sigma < 0.55) {
            bilateral05();
        } else if (sigma < 0.65) {
            bilateral06();
        } else if (sigma < 0.75) {
            bilateral07();
        } else if (sigma < 0.85) {
            bilateral08();
        } else if (sigma < 0.95) {
            bilateral09();
        } else if (sigma < 1.05) {
            bilateral10();
        } else if (sigma < 1.15) {
            bilateral11();
        } else if (sigma < 1.25) {
            bilateral12();
        } else if (sigma < 1.35) {
            bilateral13();
        } else if (sigma < 1.45) {
            bilateral14();
        } else if (sigma < 1.55) {
            bilateral15();
        } else if (sigma < 1.65) {
            bilateral16();
        } else if (sigma < 1.75) {
            bilateral17();
        } else if (sigma < 1.85) {
            bilateral18();
        } else if (sigma < 1.95) {
            bilateral19();
        } else if (sigma < 2.05) {
            bilateral20();
        } else if (sigma < 2.15) {
            bilateral21();
        } else if (sigma < 2.25) {
            bilateral22();
        } else if (sigma < 2.35) {
            bilateral23();
        } else if (sigma < 2.45) {
            bilateral24();
        } else {
            bilateral25();
        }
    }

private:
    void blOper3(int border, float c00, float c01, float c11)
    {
        const auto suly = [this](int i, int j, int a, int b) {
            return ec[source[i - a][j - b] - source[i][j] + 65536.0f];
        };

        const auto elem = [this, &suly](int i, int j, int a, int b) {
            return source[i - a][j - b] * suly(i, j, a, b);
        };

        const int rstart = border;
        const int rend = height - border;
        const int cstart = border;
        const int cend = width - border;

        #pragma omp for
        for (int i = rstart; i < rend; ++i) {
            for (int j = cstart; j < cend; ++j) {
                float v = c11 * elem(i, j, -1, -1) + c01 * elem(i, j, -1, 0) + c11 * elem(i, j, -1, 1) +
                          c01 * elem(i, j, 0, -1)  + c00 * elem(i, j, 0, 0)  + c01 * elem(i, j, 0, 1) +
                          c11 * elem(i, j, 1, -1)  + c01 * elem(i, j, 1, 0)  + c11 * elem(i, j, 1, 1);
                v /= c11 * suly(i, j, -1, -1) + c01 * suly(i, j, -1, 0) + c11 * suly(i, j, -1, 1) +
                     c01 * suly(i, j, 0, -1)  + c00 * suly(i, j, 0, 0)  + c01 * suly(i, j, 0, 1) +
                     c11 * suly(i, j, 1, -1)  + c01 * suly(i, j, 1, 0)  + c11 * suly(i, j, 1, 1);
                buffer[i][j] = v;
            }
        }

        #pragma omp for
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                if (i<rstart || j<cstart || i>=rend || j>=cend) {
                    destination[i][j] = source[i][j];
                } else {
                    destination[i][j] = buffer[i][j];
                }
            }
        }
    }

    void blOper5(int border, float c00, float c01, float c02, float c11, float c12, float c22)
    {
        const auto suly = [this](int i, int j, int a, int b) {
            return ec[source[i - a][j - b] - source[i][j] + 65536.0f];
        };

        const auto elem = [this, &suly](int i, int j, int a, int b) {
            return source[i - a][j - b] * suly(i, j, a, b);
        };

        const int rstart = border;
        const int rend = height - border;
        const int cstart = border;
        const int cend = width - border;

        #pragma omp for
        for (int i = rstart; i < rend; ++i) {
            for (int j = cstart; j < cend; ++j) {
                float v = c00 * elem(i, j, -2, -2) + c01 * elem(i, j, -2, -1) + c02 * elem(i, j, -2, 0) + c01 * elem(i, j, -2, 1) + c00 * elem(i, j, -2, 2) +
                     c01 * elem(i, j, -1, -2) + c11 * elem(i, j, -1, -1) + c12 * elem(i, j, -1, 0) + c11 * elem(i, j, -1, 1) + c01 * elem(i, j, -1, 2) +
                     c02 * elem(i, j, 0, -2) + c12 * elem(i, j, 0, -1) + c22 * elem(i, j, 0, 0) + c12 * elem(i, j, 0, 1) + c02 * elem(i, j, 0, 2) +
                     c01 * elem(i, j, 1, -2) + c11 * elem(i, j, 1, -1) + c12 * elem(i, j, 1, 0) + c11 * elem(i, j, 1, 1) + c01 * elem(i, j, 1, 2) +
                     c00 * elem(i, j, 2, -2) + c01 * elem(i, j, 2, -1) + c02 * elem(i, j, 2, 0) + c01 * elem(i, j, 2, 1) + c00 * elem(i, j, 2, 2);
                v /= c00 * suly(i, j, -2, -2) + c01 * suly(i, j, -2, -1) + c02 * suly(i, j, -2, 0) + c01 * suly(i, j, -2, 1) + c00 * suly(i, j, -2, 2) +
                     c01 * suly(i, j, -1, -2) + c11 * suly(i, j, -1, -1) + c12 * suly(i, j, -1, 0) + c11 * suly(i, j, -1, 1) + c01 * suly(i, j, -1, 2) +
                     c02 * suly(i, j, 0, -2) + c12 * suly(i, j, 0, -1) + c22 * suly(i, j, 0, 0) + c12 * suly(i, j, 0, 1) + c02 * suly(i, j, 0, 2) +
                     c01 * suly(i, j, 1, -2) + c11 * suly(i, j, 1, -1) + c12 * suly(i, j, 1, 0) + c11 * suly(i, j, 1, 1) + c01 * suly(i, j, 1, 2) +
                     c00 * suly(i, j, 2, -2) + c01 * suly(i, j, 2, -1) + c02 * suly(i, j, 2, 0) + c01 * suly(i, j, 2, 1) + c00 * suly(i, j, 2, 2);
                buffer[i][j] = v;
            }
        }

        #pragma omp for
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                if (i<rstart || j<cstart || i>=rend || j>=cend) {
                    destination[i][j] = source[i][j];
                } else {
                    destination[i][j] = buffer[i][j];
                }
            }
        }
    }

    // sigma = 0.5
    void bilateral05()
    {
        blOper3(1, 55, 7, 1);
    }

    // sigma = 0.6
    void bilateral06()
    {
        blOper3(1, 16, 4, 1);
    }

    // sigma = 0.7
    void bilateral07()
    {
        blOper5(2, 0, 1, 2, 6, 12, 22);
    }

    // sigma = 0.8
    void bilateral08()
    {
        blOper5(2, 0, 0, 1, 5, 10, 23);
    }

    // sigma = 0.9
    void bilateral09()
    {
        blOper5(2, 0, 1, 2, 6, 12, 22);
    }

    // sigma = 1.0
    void bilateral10()
    {
        blOper5(2, 0, 1, 2, 4, 7, 12);
    }

    // sigma = 1.1
    void bilateral11()
    {
        BL_BEGIN(209, 3)
        #pragma omp for
        BL_OPER7(0, 0, 1, 1, 0, 2, 5, 8, 1, 5, 18, 27, 1, 8, 27, 41)
        BL_FREE
        #pragma omp for
        BL_END(3)
    }

    // sigma = 1.2
    void bilateral12()
    {
        BL_BEGIN(322, 3)
        #pragma omp for
        BL_OPER7(0, 0, 1, 1, 0, 1, 4, 6, 1, 4, 11, 16, 1, 6, 16, 23)
        BL_FREE
        #pragma omp for
        BL_END(3)
    }

    // sigma = 1.3
    void bilateral13()
    {
        BL_BEGIN(336, 3)
        #pragma omp for
        BL_OPER7(0, 0, 1, 1, 0, 2, 4, 6, 1, 4, 11, 14, 1, 6, 14, 19)
        BL_FREE
        #pragma omp for
        BL_END(3)
    }

    // sigma = 1.4
    void bilateral14()
    {
        BL_BEGIN(195, 3)
        #pragma omp for
        BL_OPER7(0, 1, 2, 3, 1, 4, 8, 10, 2, 8, 17, 21, 3, 10, 21, 28)
        BL_FREE
        #pragma omp for
        BL_END(3)
    }

    // sigma = 1.5
    void bilateral15()
    {
        BL_BEGIN(132, 4)
        #pragma omp for
        BL_OPER9(0, 0, 0, 1, 1, 0, 1, 2, 4, 5, 0, 2, 6, 12, 14, 1, 4, 12, 22, 28, 1, 5, 14, 28, 35)
        BL_FREE
        #pragma omp for
        BL_END(4)
    }

    // sigma = 1.6
    void bilateral16()
    {
        BL_BEGIN(180, 4)
        #pragma omp for
        BL_OPER9(0, 0, 0, 1, 1, 0, 1, 2, 3, 4, 0, 2, 5, 9, 10, 1, 3, 9, 15, 19, 1, 4, 10, 19, 23)
        BL_FREE
        #pragma omp for
        BL_END(4)
    }

    // sigma = 1.7
    void bilateral17()
    {
        BL_BEGIN(195, 4)
        #pragma omp for
        BL_OPER9(0, 0, 1, 1, 1, 0, 1, 2, 3, 4, 1, 2, 5, 8, 9, 1, 3, 8, 13, 16, 1, 4, 9, 16, 19)
        BL_FREE
        #pragma omp for
        BL_END(4)
    }

    // sigma = 1.8
    void bilateral18()
    {
        BL_BEGIN(151, 4)
        #pragma omp for
        BL_OPER9(0, 0, 1, 2, 2, 0, 1, 3, 5, 5, 1, 3, 6, 10, 12, 2, 5, 10, 16, 19, 2, 5, 12, 19, 22)
        BL_FREE
        #pragma omp for
        BL_END(4)
    }

    // sigma = 1.9
    void bilateral19()
    {
        BL_BEGIN(151, 4)
        #pragma omp for
        BL_OPER9(0, 0, 1, 2, 2, 0, 1, 3, 4, 5, 1, 3, 5, 8, 9, 2, 4, 8, 12, 14, 2, 5, 9, 14, 16)
        BL_FREE
        #pragma omp for
        BL_END(4)
    }

    // sigma = 2
    void bilateral20()
    {
        BL_BEGIN(116, 5)
        #pragma omp for
        BL_OPER11(0, 0, 0, 1, 1, 1, 0, 0, 1, 2, 3, 3, 0, 1, 2, 4, 7, 7, 1, 2, 4, 8, 12, 14, 1, 3, 7, 12, 18, 20, 1, 3, 7, 14, 20, 23)
        BL_FREE
        #pragma omp for
        BL_END(5)
    }

    // sigma = 2.1
    void bilateral21()
    {
        BL_BEGIN(127, 5)
        #pragma omp for
        BL_OPER11(0, 0, 0, 1, 1, 1, 0, 0, 1, 2, 3, 3, 0, 1, 2, 4, 6, 7, 1, 2, 4, 8, 11, 12, 1, 3, 6, 11, 15, 17, 1, 3, 7, 12, 17, 19)
        BL_FREE
        #pragma omp for
        BL_END(5)
    }

    // sigma = 2.2
    void bilateral22()
    {
        BL_BEGIN(109, 5)
        #pragma omp for
        BL_OPER11(0, 0, 0, 1, 1, 2, 0, 1, 2, 3, 3, 4, 1, 2, 3, 5, 7, 8, 1, 3, 5, 9, 12, 13, 1, 3, 7, 12, 16, 18, 2, 4, 8, 13, 18, 20)
        BL_FREE
        #pragma omp for
        BL_END(5)
    }

    // sigma = 2.3
    void bilateral23()
    {
        BL_BEGIN(132, 5)
        #pragma omp for
        BL_OPER11(0, 0, 1, 1, 1, 1, 0, 1, 1, 2, 3, 3, 1, 1, 3, 5, 6, 7, 1, 2, 5, 7, 10, 11, 1, 3, 6, 10, 13, 14, 1, 3, 7, 11, 14, 16)
        BL_FREE
        #pragma omp for
        BL_END(5)
    }

    // sigma = 2.4
    void bilateral24()
    {
        BL_BEGIN(156, 5)
        #pragma omp for
        BL_OPER11(0, 0, 1, 1, 1, 1, 0, 1, 1, 2, 3, 3, 1, 1, 3, 4, 5, 6, 1, 2, 4, 6, 8, 9, 1, 3, 5, 8, 10, 11, 1, 3, 6, 9, 11, 12)
        BL_FREE
        #pragma omp for
        BL_END(5)
    }

    // sigma = 2.5
    void bilateral25()
    {
        BL_BEGIN(173, 5)
        #pragma omp for
        BL_OPER11(0, 0, 1, 1, 1, 1, 0, 1, 1, 2, 3, 3, 1, 1, 2, 4, 5, 5, 1, 2, 4, 5, 7, 7, 1, 3, 5, 7, 9, 9, 1, 3, 5, 7, 9, 10)
        BL_FREE
        #pragma omp for
        BL_END(5)
    }

    float** const source;
    float** const destination;
    float** const buffer;
    const int width;
    const int height;
    const double sigma;
    const LUTf ec;
};

rtengine::Bilateral::Bilateral(
    float** source,
    float** destination,
    float** buffer,
    int width,
    int height,
    double sigma,
    double sensitivity
) :
    implementation(
        new Implementation(
            source,
            destination,
            buffer,
            width,
            height,
            sigma,
            sensitivity
        )
    )
{
};

rtengine::Bilateral::~Bilateral() = default;

void rtengine::Bilateral::compute()
{
    implementation->compute();
}
