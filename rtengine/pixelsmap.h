#pragma once

/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2019 Gabor Horvath <hgabor@rawtherapee.com>
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

#include <cstdint>
#include <cstring>
#include <vector>

#include "noncopyable.h"

namespace rtengine
{

struct badPix {
    uint16_t x;
    uint16_t y;
    badPix(uint16_t xc, uint16_t yc): x(xc), y(yc) {}
};

class PixelsMap :
    public NonCopyable
{
    int w; // line width in base_t units
    int h; // height
    typedef unsigned long base_t;
    static constexpr size_t base_t_size = sizeof(base_t);
    base_t *pm;

public:
    PixelsMap(int width, int height)
        : w((width / (base_t_size * 8)) + 1), h(height), pm(new base_t[h * w])
    {
        clear();
    }

    ~PixelsMap()
    {
        delete [] pm;
    }
    int width() const
    {
        return w;
    }
    int height() const
    {
        return h;
    }

    // if a pixel is set returns true
    bool get(int x, int y) const
    {
        return (pm[y * w + x / (base_t_size * 8)] & (base_t)1 << (x % (base_t_size * 8))) != 0;
    }

    // set a pixel
    void set(int x, int y)
    {
        pm[y * w + x / (base_t_size * 8)] |= (base_t)1 << (x % (base_t_size * 8)) ;
    }

    // set pixels from a list
    int set(const std::vector<badPix>& bp)
    {
        for (const auto& bad_pix : bp) {
            set(bad_pix.x, bad_pix.y);
        }

        return bp.size();
    }

    void clear()
    {
        std::memset(pm, 0, h * w * base_t_size);
    }
    // return 0 if at least one pixel in the word(base_t) is set, otherwise return the number of pixels to skip to the next word base_t
    int skipIfZero(int x, int y) const
    {
        return pm[y * w + x / (base_t_size * 8)] == 0 ? base_t_size * 8 - x % (base_t_size * 8) : 0;
    }
};

}
