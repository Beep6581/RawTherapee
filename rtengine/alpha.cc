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

#include "alpha.h"

namespace rtengine
{

Alpha::Alpha () {}

Alpha::Alpha (int width, int height)
{
    if (width > 0 && height > 0) {
        surface = Cairo::ImageSurface::create(Cairo::FORMAT_A8, width, height);
    }
}

/*
Alpha::~Alpha () {
    surface->unreference();
}
*/

void Alpha::setSize(int width, int height)
{
    if (width > 0 && height > 0) {
        if (surface) {
            if (width != getWidth() && height != getHeight()) {
                surface.clear();  // does this delete the referenced object? Unreferencing doesn't work, since Cairo expect to have a non null refCount in the destructor!
            } else {
                return;
            }
        }

        surface = Cairo::ImageSurface::create(Cairo::FORMAT_A8, width, height);
    } else if (surface) {
        surface.clear();
    }
}

int Alpha::getWidth()
{
    if (surface) {
        return surface->get_width();
    }

    return -1;
}

int Alpha::getHeight()
{
    if (surface) {
        return surface->get_height();
    }

    return -1;
}

}
