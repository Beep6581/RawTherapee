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
#ifndef _ALPHA_H_
#define _ALPHA_H_

#include <gtkmm.h>
#include <assert.h>

#define CHECK_BOUNDS 0

namespace rtengine
{

/// Alpha channel class (8 bits)
class Alpha
{
protected:
    Cairo::RefPtr<Cairo::ImageSurface> surface;

public:
    Alpha ();
    Alpha (int width, int height);
    //~Alpha ();

    void setSize(int width, int height);
    int getWidth();
    int getHeight();

    const Cairo::RefPtr<Cairo::ImageSurface> getSurface ();

    // TODO: to make the editing faster, we should add an iterator class

    // Will send back the start of a row
    unsigned char* operator () (unsigned row) const;
    // Will send back a value at a given row, col position
    unsigned char& operator () (unsigned row, unsigned col);
    const unsigned char operator () (unsigned row, unsigned col) const;
};



inline const Cairo::RefPtr<Cairo::ImageSurface> Alpha::getSurface () {
    return surface; // to be used in bitmap edition
}

inline const unsigned char Alpha::operator () (unsigned row,
        unsigned col) const {
    return *(surface->get_data () + row * surface->get_width () + col);
}

inline unsigned char& Alpha::operator () (unsigned row, unsigned col) {
    return *(surface->get_data () + row * surface->get_width () + col);
}

inline unsigned char* Alpha::operator () (unsigned row) const {
    return surface->get_data () + row * surface->get_width ();
}

}

#endif
