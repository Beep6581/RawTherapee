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

#ifndef _IMAGEDIMENSIONS_
#define _IMAGEDIMENSIONS_

class PreviewProps
{
public:
    int x, y, w, h, skip;
    PreviewProps (int x, int y, int w, int h, int skip);

    void set (int x, int y, int w, int h, int skip);
};

/*
 * Description of an image dimension, with getter and setter
 */
class ImageDimensions
{

public:
    int width;
    int height;

public:
    ImageDimensions ();
    int getW ();
    int getH ();
    int getWidth () const;
    int getHeight () const;
    void transform (const PreviewProps & pp, int tran, int &sx1, int &sy1, int &sx2, int &sy2);
};

inline PreviewProps::PreviewProps (int x, int y, int w, int h, int skip) :
        x (x), y (y), w (w), h (h), skip (skip) {
}

inline void PreviewProps::set (int x, int y, int w, int h, int skip) {
    this->x = x;
    this->y = y;
    this->w = w;
    this->h = h;
    this->skip = skip;
}

inline ImageDimensions::ImageDimensions () :
        width (-1), height (-1) {
}

inline int ImageDimensions::getW () {
    return width;
}

inline int ImageDimensions::getH () {
    return height;
}

inline int ImageDimensions::getWidth () const {
    return width;
}

inline int ImageDimensions::getHeight () const {
    return height;
}

#endif
