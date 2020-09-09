/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "coord2d.h"
#include "procparams.h"
#include "imagesource.h"

namespace rtengine {

class PerspectiveCorrection {
public:
    struct Params
    {
        double angle;
        double pitch;
        double yaw;
    };

    static Params autocompute(ImageSource *src, bool corr_pitch, bool corr_yaw, const procparams::ProcParams *pparams, const FramesMetaData *metadata);

    //static void autocrop(int width, int height, bool fixratio, const procparams::PerspectiveParams &params, const FramesMetaData *metadata, int &x, int &y, int &w, int &h);
};

} // namespace rtengine
