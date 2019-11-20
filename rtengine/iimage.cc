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

#include "color.h"
#include "procparams.h"
#include "rtengine.h"

const char rtengine::sImage8[] =     "Image8";
const char rtengine::sImage16[] =    "Image16";
const char rtengine::sImagefloat[] = "Imagefloat";
int rtengine::getCoarseBitMask( const procparams::CoarseTransformParams &coarse)
{
    int tr = TR_NONE;

    if (coarse.rotate == 90) {
        tr |= TR_R90;
    } else if (coarse.rotate == 180) {
        tr |= TR_R180;
    } else if (coarse.rotate == 270) {
        tr |= TR_R270;
    }

    if (coarse.hflip) {
        tr |= TR_HFLIP;
    }

    if (coarse.vflip) {
        tr |= TR_VFLIP;
    }

    return tr;
}

const LUTf& rtengine::getigammatab() {
    return Color::igammatab_srgb;
}
