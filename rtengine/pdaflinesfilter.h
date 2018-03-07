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

#include "rawimagesource.h"
#include <vector>

namespace rtengine {

class PDAFLinesFilter {
public:
    PDAFLinesFilter(RawImage *ri);
    int mark(array2D<float> &rawData, PixelsMap &bpMap);

private:
    int markLine(array2D<float> &rawData, PixelsMap &bpMap, int y);
    
    RawImage *ri_;
    int W_;
    int H_;
    std::vector<int> pattern_;
    int offset_;
    std::vector<bool> rowmap_;
};

} // namespace rtengine
