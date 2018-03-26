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
#include "camconst.h"

namespace rtengine {

extern const Settings *settings;


namespace {

class PDAFGreenEqulibrateThreshold: public RawImageSource::GreenEqulibrateThreshold {
    static constexpr float BASE_THRESHOLD = 0.5f;
    static constexpr int TILE_SIZE = 200;
    static constexpr float AREA = TILE_SIZE * TILE_SIZE;
    static constexpr int PIXEL_COUNT_FACTOR = 12;
    
public:
    PDAFGreenEqulibrateThreshold(int w, int h):
        RawImageSource::GreenEqulibrateThreshold(BASE_THRESHOLD),
        w_(w),
        h_(h)
    {
        int ctiles = w_ / TILE_SIZE;
        int rtiles = h_ / TILE_SIZE;
        tiles_.resize(rtiles+1, std::vector<float>(ctiles+1));
    }

    void processTiles()
    {
        for(size_t i = 0; i < tiles_.size(); ++i) {
            for(size_t j = 0; j < tiles_[i].size(); ++j) {
                tiles_[i][j] = tiles_[i][j] * PIXEL_COUNT_FACTOR / (AREA * AREA);
            }
        }
    }

    void increment(int row, int col)
    {
        auto &r = tiles_[row / TILE_SIZE];
        ++r[col / TILE_SIZE];
    }
    
    float operator()(int row, int col) const
    {
        int y = row / TILE_SIZE;
        int x = col / TILE_SIZE;

        int cy = y * TILE_SIZE + TILE_SIZE/2;
        int cx = x * TILE_SIZE + TILE_SIZE/2;

        int x1 = col > cx ? x+1 : x-1;
        int y1 = row > cy ? y+1 : y-1;

        float fxy = tile_factor(y, x);
        float f = 0.f;

        if (x1 >= 0 && size_t(x1) < tiles_[y].size()) {
            if (y1 >= 0 && size_t(y1) < tiles_.size()) {
                // bilinear interpolation
                float fx1y = tile_factor(y, x1);
                float fx1y1 = tile_factor(y1, x1);
                float fxy1 = tile_factor(y1, x);

                // x direction
                int d = std::abs(cx - col);
                float f1 = fxy * (TILE_SIZE - d) + fx1y * float(d);
                float f2 = fxy1 * (TILE_SIZE - d) + fx1y1 * float(d);
                // y direction
                d = std::abs(cy - row);
                f = (f1 * (TILE_SIZE - d) + f2 * float(d));
            } else {
                float f2 = tile_factor(y, x1);
                int d = std::abs(cx - col);
                f = (fxy * float(TILE_SIZE - d) + f2 * float(d)) * TILE_SIZE;
            }
        } else if (y1 >= 0 && size_t(y1) < tiles_.size()) {
            float f2 = tile_factor(y1, x);
            int d = std::abs(cy - row);
            f = (fxy * float(TILE_SIZE - d) + f2 * float(d)) * TILE_SIZE;
        } else {
            f = fxy * AREA;
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
        return tiles_[y][x];
    }
    
    int w_;
    int h_;
    std::vector<std::vector<float>> tiles_;
};


class PDAFLineDenoiseRowFilter: public RawImageSource::CFALineDenoiseRowBlender {
public:
    PDAFLineDenoiseRowFilter(const std::vector<int> &pattern, int offset):
        pattern_(pattern),
        offset_(offset)
    {}

    float operator()(int row) const
    {
        static constexpr float BORDER[] = { 1.f, 1.f, 0.8f, 0.5f, 0.2f };
        static constexpr int BORDER_WIDTH = sizeof(BORDER)/sizeof(float) - 1;
        
        if (!pattern_.empty()) {
            int key = (row - offset_) % pattern_.back();
            auto it = std::lower_bound(pattern_.begin(), pattern_.end(), key);

            int b = *it;
            int d = b - key;

            if (it > pattern_.begin()) {
                int b2 = *(it-1);
                int d2 = key - b2;
                float f = BORDER[std::min(std::min(d, d2), BORDER_WIDTH)];
                return f;
            } else {
                float f = BORDER[std::min(d, BORDER_WIDTH)];
                return f;
            }
        }
        return 0.f;
    }

private:
    std::vector<int> pattern_;
    int offset_;
};

} // namespace



PDAFLinesFilter::PDAFLinesFilter(RawImage *ri):
    ri_(ri),
    W_(ri->get_width()),
    H_(ri->get_height())
{
    gthresh_ = new PDAFGreenEqulibrateThreshold(W_, H_);

    CameraConstantsStore* ccs = CameraConstantsStore::getInstance();
    CameraConst *cc = ccs->get(ri_->get_maker().c_str(), ri_->get_model().c_str());

    if (cc) {
        pattern_ = cc->get_pdafPattern();
        if(!pattern_.empty()) {
            offset_ = cc->get_pdafOffset();
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


std::unique_ptr<RawImageSource::CFALineDenoiseRowBlender> PDAFLinesFilter::lineDenoiseRowBlender()
{
    return std::unique_ptr<RawImageSource::CFALineDenoiseRowBlender>(new PDAFLineDenoiseRowFilter(pattern_, offset_));
}


int PDAFLinesFilter::markLine(array2D<float> &rawData, PixelsMap &bpMap, int y)
{
    rowmap_.clear();
    rowmap_.resize((W_+1)/2, false);
    int marked = 0;
    
    for (int x = 1 + (ri_->FC(y, 0) & 1); x < W_-1; x += 2) {
        const float
            g0 = rawData[y][x],
            g1 = rawData[y-1][x+1],
            g2 = rawData[y+1][x+1],
            g3 = rawData[y-1][x-1],
            g4 = rawData[y+1][x-1];
        if (g0 > max(g1, g2, g3, g4)) {
            const float gu = g2 + g4;
            const float gd = g1 + g3;
            const float gM = max(gu, gd);
            const float gm = min(gu, gd);
            const float d = (gM - gm) / gM;
            if (d < 0.2f && (1.f - (gm + gM)/(4.f * g0)) > std::min(d, 0.1f)) {
                rowmap_[x/2] = true;
            }
        }
    }

    PDAFGreenEqulibrateThreshold *m = static_cast<PDAFGreenEqulibrateThreshold *>(gthresh_);

    for (int x = 2 + (ri_->FC(y, 1) & 1); x < W_-2; x += 2) {
        const int i = x/2;
        if (rowmap_[i+1]) {
            if (rowmap_[i]) {
                if (rowmap_[i-1]) {
                    for (int xx = x-2; xx <= x+2; ++xx) {
                        if (!bpMap.get(xx, y)) {
                            bpMap.set(xx, y);
                            m->increment(y, xx);
                            ++marked;
                        }
                    }
                }
            } else {
                x += 2;
            }
        } else {
            x += 4;
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
    static_cast<PDAFGreenEqulibrateThreshold *>(gthresh_)->processTiles();
    return found;
}

} // namespace rtengine
