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
#include <bilateral2.h>

void bilateral_unsigned (unsigned short** src, unsigned short** dst, unsigned short** buffer, Dim dim, double sigma, double sens) {

    bilateral<unsigned short, unsigned int> (src, dst, buffer, dim, sigma, sens);
}

void bilateral_signed (short** src, short** dst, short** buffer, Dim dim, double sigma, double sens) {

    bilateral<short, int> (src, dst, buffer, dim, sigma, sens);
}

void bilateral_box_unsigned (unsigned short** src, unsigned short** dst, int W, int H, int sigmar, double sigmas, bilateralparams row) {

    bilateral<unsigned short> (src, dst, W, H, sigmar, sigmas, row.row_from, row.row_to);
}

