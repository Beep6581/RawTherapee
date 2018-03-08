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

#include "pdaflinesfilter.h"
#include "settings.h"
#include <iostream>

namespace rtengine {

extern const Settings *settings;


namespace {

class PDAFGreenEqulibrateThreshold: public RawImageSource::GreenEqulibrateThreshold {
    static constexpr float base_threshold = 0.5f;
public:
    PDAFGreenEqulibrateThreshold(int w, int h, int ntiles=20):
        RawImageSource::GreenEqulibrateThreshold(base_threshold),
        w_(w),
        h_(h)
    {
        tw_ = (w_ + 1) / ntiles;
        th_ = (h_ + 1) / ntiles;
        area_ = tw_ * th_;
        tiles_.resize(ntiles+1, std::vector<int>(ntiles+1));
    }

    void increment(int row, int col)
    {
        auto &r = tiles_[row / th_];
        ++r[col / tw_];
    }
    
    float operator()(int row, int col) const
    {
        int y = row / th_;
        int x = col / tw_;

        float f = tile_factor(y, x);
        int cy = y * th_ + th_/2;
        int cx = x * tw_ + tw_/2;

        if (std::abs(y - cy) > std::abs(x - cx)) {
            int y1 = y > cy ? y+1 : y-1;
            if (y1 >= 0 && size_t(y1) < tiles_.size()) {
                float f2 = tile_factor(y1, x);
                int d = std::abs(cy - row);
                f = f * float(th_ - d)/float(th_) + f2 * float(d)/float(th_);
            }
        } else {
            int x1 = x > cx ? x+1 : x-1;
            if (x1 >= 0 && size_t(x1) < tiles_[y].size()) {
                float f2 = tile_factor(y, x1);
                int d = std::abs(cx - col);
                f = f * float(tw_ - d)/float(tw_) + f2 * float(d)/float(tw_);
            }
        }
        return thresh_ * f;
    }

    void print() const
    {
        std::cout << "PDAFGreenEqulibrateThreshold:\n";
        for (size_t row = 0; row < tiles_.size(); ++row) {
            for (size_t col = 0; col < tiles_.size(); ++col) {
                std::cout << " " << tile_factor(row, col);
            }
            std::cout << std::endl;
        }
    }

private:
    float tile_factor(int y, int x) const
    {
        return float(tiles_[y][x] * 12) / area_;
    }
    
    int w_;
    int h_;
    int ntiles_;
    int tw_;
    int th_;
    float area_;
    std::vector<std::vector<int>> tiles_;
};

} // namespace



PDAFLinesFilter::PDAFLinesFilter(RawImage *ri):
    ri_(ri),
    W_(ri->get_width()),
    H_(ri->get_height())
{
    gthresh_ = new PDAFGreenEqulibrateThreshold(W_, H_);
    
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


PDAFLinesFilter::~PDAFLinesFilter()
{
    delete gthresh_;
}


RawImageSource::GreenEqulibrateThreshold &PDAFLinesFilter::greenEqThreshold()
{
    return *gthresh_;
}


int PDAFLinesFilter::markLine(array2D<float> &rawData, PixelsMap &bpMap, int y)
{
    rowmap_.clear();
    rowmap_.resize((W_+1)/2, false);
    int marked = 0;
    
    for (int x = 1; x < W_-1; ++x) {
        if (ri_->FC(y, x) == 1) {
            const float
                g0 = rawData[y][x],
                g1 = rawData[y-1][x+1],
                g2 = rawData[y+1][x+1],
                g3 = rawData[y-1][x-1],
                g4 = rawData[y+1][x-1];
            if (g0 > max(g1, g2, g3, g4)) {
                const float gu = (g2 + g4) / 2.f;
                const float gd = (g1 + g3) / 2.f;
                const float gM = max(gu, gd);
                const float gm = min(gu, gd);
                const float d = (gM - gm) / gM;
                if (d < 0.2f && (g0 - (gm + gM)/2.f) / g0 > std::min(d, 0.1f)) {
                    rowmap_[x/2] = true;
                }
            }
        }
    }

    PDAFGreenEqulibrateThreshold *m = static_cast<PDAFGreenEqulibrateThreshold *>(gthresh_);

    for (int x = 2; x < W_-2; ++x) {
        if (ri_->FC(y, x) == 1) {
            const int i = x/2;
            if (rowmap_[i-1] && rowmap_[i] && rowmap_[i+1]) {
                for (int xx = x-2; xx <= x+2; ++xx) {
                    if (!bpMap.get(xx, y)) {
                        bpMap.set(xx, y);
                        m->increment(y, xx);
                        ++marked;
                    }
                }
            }
        }
    }

    return marked;
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
            int n = markLine(rawData, bpMap, y) + markLine(rawData, bpMap, y-1) + markLine(rawData, bpMap, y+1);
            if (n) {
                found += n;
                if (settings->verbose) {
                    std::cout << "marked " << n << " pixels in PDAF line at " << y << std::endl;
                }
            }
        } else if (y > yy) {
            ++idx;
            if (idx >= pattern_.size()) {
                idx = 0;
                off += pattern_.back();
            }
        }
    }

    if (settings->verbose) {
        static_cast<PDAFGreenEqulibrateThreshold *>(gthresh_)->print();
    }
    
    return found;
}

} // namespace rtengine
