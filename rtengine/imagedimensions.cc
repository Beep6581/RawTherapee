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

#include "imagedimensions.h"
#include "rtengine.h"

void ImageDimensions::transform (const PreviewProps & pp, int tran, int &sx1, int &sy1, int &sx2, int &sy2)
{

    int sw = width, sh = height;

    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = height;
        sh = width;
    }

    int ppx = pp.x, ppy = pp.y;

    if (tran & TR_HFLIP) {
        ppx = sw - pp.x - pp.w;
    }

    if (tran & TR_VFLIP) {
        ppy = sh - pp.y - pp.h;
    }

    sx1 = ppx;
    sy1 = ppy;
    sx2 = ppx + pp.w;
    sy2 = ppy + pp.h;

    if ((tran & TR_ROT) == TR_R180) {
        sx1 = width - ppx - pp.w;
        sy1 = height - ppy - pp.h;
        sx2 = sx1 + pp.w;
        sy2 = sy1 + pp.h;
    } else if ((tran & TR_ROT) == TR_R90) {
        sx1 = ppy;
        sy1 = height - ppx - pp.w;
        sx2 = sx1 + pp.h;
        sy2 = sy1 + pp.w;
    } else if ((tran & TR_ROT) == TR_R270) {
        sx1 = width - ppy - pp.h;
        sy1 = ppx;
        sx2 = sx1 + pp.h;
        sy2 = sy1 + pp.w;
    }

    //printf ("ppx %d ppy %d ppw %d pph %d s: %d %d %d %d\n",pp.x, pp.y,pp.w,pp.h,sx1,sy1,sx2,sy2);
    if (sx1 < 0) {
        sx1 = 0;
    }

    if (sy1 < 0) {
        sy1 = 0;
    }
}

