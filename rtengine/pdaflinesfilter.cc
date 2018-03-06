/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "pdaflinesfilter.h"
#include "settings.h"
#include <iostream>

namespace rtengine {

extern const Settings *settings;


PDAFLinesFilter::PDAFLinesFilter(RawImage *ri):
    ri_(ri),
    W_(ri->get_width()),
    H_(ri->get_height())
{
    if (ri_->get_maker() == "Sony") {
        if (ri_->get_model() == "ILCE-7M3") {
            // A7III, from https://www.dpreview.com/forums/post/60843139
            // in the original post:
            //   P 5 P 17 P 11 P 11 P 17 P 11 P 5 P 11 P 11 P 11 P 17 P 11 P 5 P 11 P 11 P 17 P 5 P 11 P 17 P 5 P 17 P 5 P 11 P 11 P 11 P 17 P 5 P 11 P 11 P 11 P 5 P 17 P 5 P 17 P 11
            //   
            // rotated to match the start of the frame
            //   P 11 P 11 P 11 P 17 P 11 P 5 P 11 P 11 P 17 P 5 P 11 P 17 P 5 P 17 P 5 P 11 P 11 P 11 P 17 P 5 P 11 P 11 P 11 P 5 P 17 P 5 P 17 P 11 P 5 P 17 P 11 P 11 P 17 P 11 P 5
            pattern_ = {
                0, 12, 24, 36, 54, 66, 72, 84, 96, 114, 120, 132, 150, 156, 174, 180, 192, 204, 216, 234, 240, 252, 264, 276, 282, 300, 306, 324, 336, 342, 360, 372, 384, 402, 414, 420
            };
            offset_ = 9;
        } else if (ri_->get_model() == "ILCE-6000") {
            // detected by hand, using the picture from https://www.dpreview.com/forums/thread/3923513
            // P 11 P 23 P 17 P 17 P 17 P 23 P 11 P 17 P 17 P 17 P 23 P 11 P 23 P 11 P 17 P 23 P 11 P 17 P 17 P 23 P 17 P 11 P 17 P 17 P 17 P 23 P 17 P 11 P 17 P 17 P 23 P 11 P 17 P 11 P 23
            pattern_ = {
                0, 12, 36, 54, 72, 90, 114, 126, 144, 162, 180, 204, 216, 240, 252, 270, 294, 306, 324, 342, 366, 384, 396, 414, 432, 450, 474, 492, 504, 522, 540, 564, 576, 594, 606, 630
            };
            offset_ = 3;
        } else if (ri_->get_model() == "ILCE-9") {
            // the A9 is the same as the A7III, rotated of 1 position
            // source: https://www.dpreview.com/forums/post/60857788
            // P 11 P 11 P 11 P 17 P 11 P 5 P 11 P 11 P 17 P 5 P 11 P 17 P 5 P 17 P 5 P 11 P 11 P 11 P 17 P 5 P 11 P 11 P 11 P 5 P 17 P 5 P 17 P 11 P 5 P 17 P 11 P 11 P 17 P 11 P 5
            pattern_ = {
                0, 12, 24, 36, 54, 66, 72, 84, 96, 114, 120, 132, 150, 156, 174, 180, 192, 204, 216, 234, 240, 252, 264, 276, 282, 300, 306, 324, 336, 342, 360, 372, 384, 402, 414, 420
            };
            offset_ = -7;
        }
    }
}


int PDAFLinesFilter::mark(array2D<float> &rawData, PixelsMap &bpMap)
{
    if (pattern_.empty()) {
        if (settings->verbose) {
            std::cout << "no PDAF pattern known for " << ri_->get_maker() << " " << ri_->get_model() << std::endl;
        }
        return 0;
    }
        
    size_t idx = 0;
    int off = offset_;
        
    int found = 0;
    for (int y = 1; y < H_-1; ++y) {
        int yy = pattern_[idx] + off;
        if (y == yy) {
            int n = 0;
            for (int x = 1; x < W_-1; ++x) {
                if (ri_->FC(y, x) == 1) {
                    const float
                        g0 = rawData[y][x],
                        g1 = rawData[y-1][x+1],
                        g2 = rawData[y+1][x+1],
                        g3 = rawData[y-1][x-1],
                        g4 = rawData[y+1][x-1];
                    const float g = max(g0, g1, g2);
                    if (g > max(g3, g4)) {
                        int row = y;
                        if (g == g1) {
                            --row;
                        } else if (g == g2) {
                            ++row;
                        }
                        bpMap.set(x, row);
                        bpMap.set(x-1, row);
                        bpMap.set(x+1, row);
                        n += 2;
                    }
                }
            }
            found += n;
            if (n && settings->verbose) {
                std::cout << "marked " << n << " pixels in PDAF line at " << y << std::endl;
            }
        } else if (y > yy) {
            ++idx;
            if (idx >= pattern_.size()) {
                idx = 0;
                off += pattern_.back();
            }
        }
    }
    return found;
}

} // namespace rtengine
