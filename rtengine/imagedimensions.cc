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

PreviewProps::PreviewProps(int _x, int _y, int _width, int _height, int _skip) :
    x(_x),
    y(_y),
    width(_width),
    height(_height),
    skip(_skip)
{
}

int PreviewProps::getX() const
{
    return x;
}

int PreviewProps::getY() const
{
    return y;
}

int PreviewProps::getWidth() const
{
    return width;
}

int PreviewProps::getHeight() const
{
    return height;
}

int PreviewProps::getSkip() const
{
    return skip;
}

ImageDimensions::ImageDimensions() :
    width(-1),
    height(-1)
{
}

void ImageDimensions::transform(const PreviewProps& pp, int tran, int& sx1, int& sy1, int& sx2, int& sy2) const
{
    int sw = width;
    int sh = height;

    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        std::swap(sw, sh);
    }

    int ppx = pp.getX();
    int ppy = pp.getY();

    if (tran & TR_HFLIP) {
        ppx = sw - pp.getX() - pp.getWidth();
    }

    if (tran & TR_VFLIP) {
        ppy = sh - pp.getY() - pp.getHeight();
    }

    sx1 = ppx;
    sy1 = ppy;
    sx2 = ppx + pp.getWidth();
    sy2 = ppy + pp.getHeight();

    if ((tran & TR_ROT) == TR_R180) {
        sx1 = width - ppx - pp.getWidth();
        sy1 = height - ppy - pp.getHeight();
        sx2 = sx1 + pp.getWidth();
        sy2 = sy1 + pp.getHeight();
    } else if ((tran & TR_ROT) == TR_R90) {
        sx1 = ppy;
        sy1 = height - ppx - pp.getWidth();
        sx2 = sx1 + pp.getHeight();
        sy2 = sy1 + pp.getWidth();
    } else if ((tran & TR_ROT) == TR_R270) {
        sx1 = width - ppy - pp.getHeight();
        sy1 = ppx;
        sx2 = sx1 + pp.getHeight();
        sy2 = sy1 + pp.getWidth();
    }

    if (sx1 < 0) {
        sx1 = 0;
    }

    if (sy1 < 0) {
        sy1 = 0;
    }
}

