/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (C) 2019 Ingo Weyrich <heckflosse67@gmx.de>
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

namespace rtengine
{
void compute13x13kernel2(float sigma, float kernel[13][13]);
void compute9x9kernel2(float sigma, float kernel[9][9]);
void compute7x7kernel2(float sigma, float kernel[7][7]);
void compute5x5kernel2(float sigma, float kernel[5][5]);
void compute3x3kernel2(float sigma, float kernel[3][3]);

void boxblur(float** src, float** dst, int radius, int W, int H, bool multiThread);
void boxblur(float* src, float* dst, int radius, int W, int H, bool multiThread);
void boxabsblur(float** src, float** dst, int radius, int W, int H, bool multiThread);
void boxabsblur(float* src, float* dst, int radius, int W, int H, bool multiThread);

}
