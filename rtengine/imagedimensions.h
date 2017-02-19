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

#pragma once

class PreviewProps
{
public:
    PreviewProps(int _x, int _y, int _width, int _height, int _skip);

    int getX() const;
    int getY() const;
    int getWidth() const;
    int getHeight() const;
    int getSkip() const;

private:
    int x;
    int y;
    int width;
    int height;
    int skip;
};

/*
 * Description of an image dimension, with getter
 */
class ImageDimensions
{
public:
    ImageDimensions();

    int getWidth() const
    {
        return width;
    }
    int getHeight() const
    {
        return height;
    }

    void transform(const PreviewProps& pp, int tran, int& sx1, int& sy1, int& sx2, int& sy2) const;

protected:
    int width;
    int height;
};
