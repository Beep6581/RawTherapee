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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <vector>

#include "noncopyable.h"
#include "rawimagesource.h"

namespace rtengine
{

class RawImage;

class PDAFLinesFilter: public rtengine::NonCopyable
{
public:
    explicit PDAFLinesFilter(RawImage *ri);
    ~PDAFLinesFilter();
    
    int mark(const array2D<float> &rawData, PixelsMap &bpMap);
    RawImageSource::GreenEqulibrateThreshold &greenEqThreshold();
    std::unique_ptr<RawImageSource::CFALineDenoiseRowBlender> lineDenoiseRowBlender();

private:
    int markLine(const array2D<float>& rawData, PixelsMap &bpMap, int y);

    RawImage *ri_;
    int W_;
    int H_;
    std::vector<int> pattern_;
    int offset_;
    std::vector<bool> rowmap_;
    RawImageSource::GreenEqulibrateThreshold *gthresh_;
};

} // namespace rtengine
