/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

#pragma once

template<typename T> class array2D;

namespace rtengine
{


void guidedFilter(const array2D<float> &guide, const array2D<float> &src, array2D<float> &dst, int r, float epsilon, bool multithread, int subsampling=0);

void guidedFilterLog(float base, array2D<float> &chan, int r, float eps, bool multithread, int subsampling=0);

void guidedFilterLog(const array2D<float> &guide, float base, array2D<float> &chan, int r, float eps, bool multithread, int subsampling=0);

} // namespace rtengine
