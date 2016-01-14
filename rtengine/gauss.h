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
#ifndef _GAUSS_H_
#define _GAUSS_H_

#include <string.h>
enum eGaussType {GAUSS_STANDARD, GAUSS_MULT, GAUSS_DIV};
#define GAUSS_SKIP 0.25
#define GAUSS_3X3_LIMIT 0.6
#define GAUSS_DOUBLE 70.0

void gaussianBlur(float** src, float** dst, const int W, const int H, const double sigma, float *buffer = NULL, eGaussType gausstype = GAUSS_STANDARD, float** buffer2 = NULL);

#endif