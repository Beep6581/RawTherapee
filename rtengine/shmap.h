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
#ifndef __SHMAP__
#define __SHMAP__

#include <image16.h>

namespace rtengine {

class SHMap {

    public: 
        int W, H;
        unsigned short** map;
        unsigned short   max, min, avg;
        
     SHMap (int w, int h);
    ~SHMap ();

    void update (Image16* img, unsigned short** buffer, double radius, double lumi[3], bool hq);
    void forceStat (unsigned short max_, unsigned short min_, unsigned short avg_);
};
};
#endif
