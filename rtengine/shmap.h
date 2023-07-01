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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "noncopyable.h"

template<typename T>
class LUT;

using LUTf = LUT<float>;

namespace rtengine
{

class Imagefloat;
class LabImage;

class SHMap :
    public NonCopyable
{

public:
    float** map;
    float   max_f, min_f, avg;

    SHMap (int w, int h);
    ~SHMap ();
    void updateLab (LabImage* img, double radius, bool hq, int skip);

    void update (Imagefloat* img, double radius, double lumi[3], bool hq, int skip);
    void updateL (float** L, double radius, bool hq, int skip);
    void forceStat (float max_, float min_, float avg_);

private:
    int W, H;
    void fillLuminanceLab( LabImage * img, float **luminance);

    void dirpyr_shmap(float ** data_fine, float ** data_coarse, int width, int height, const LUTf& rangefn, int level, int scale);

};

}
