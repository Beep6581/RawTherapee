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
#include <gauss.h>

void gaussHorizontal_unsigned (unsigned short** src, unsigned short** dst, AlignedBuffer<double>* buffer, int W, int row_from, int row_to, double sigma) {

    gaussHorizontal<unsigned short> (src, dst, buffer, W, row_from, row_to, sigma);
}

void gaussVertical_unsigned (unsigned short** src, unsigned short** dst, AlignedBuffer<double>* buffer, int H, int col_from, int col_to, double sigma) {

    gaussVertical<unsigned short> (src, dst, buffer, H, col_from, col_to, sigma);
}

void gaussHorizontal_signed (short** src, short** dst, AlignedBuffer<double>* buffer, int W, int row_from, int row_to, double sigma) {

    gaussHorizontal<short> (src, dst, buffer, W, row_from, row_to, sigma);
}

void gaussVertical_signed (short** src, short** dst, AlignedBuffer<double>* buffer, int H, int col_from, int col_to, double sigma) {

    gaussVertical<short> (src, dst, buffer, H, col_from, col_to, sigma);
}

void gaussHorizontal_float (float** src, float** dst, AlignedBuffer<double>* buffer, int W, int row_from, int row_to, double sigma) {

    gaussHorizontal<float> (src, dst, buffer, W, row_from, row_to, sigma);
}

void gaussVertical_float (float** src, float** dst, AlignedBuffer<double>* buffer, int H, int col_from, int col_to, double sigma) {

    gaussVertical<float> (src, dst, buffer, H, col_from, col_to, sigma);
}
