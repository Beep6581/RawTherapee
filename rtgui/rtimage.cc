/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2018 Jean-Christophe FRISCH <natureh.510@gmail.com>
 *  Copyright (c) 2022 Pierre CABRERA <pierre.cab@gmail.com>
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

#include "rtimage.h"

RTImage::RTImage () {}

RTImage::RTImage (const Glib::ustring& iconName, const Gtk::IconSize iconSize) :
    Gtk::Image(iconName, iconSize),
    size(iconSize)
{
}

void RTImage::set_from_icon_name(const Glib::ustring& iconName)
{
    Gtk::Image::set_from_icon_name(iconName, this->size);
}

void RTImage::set_from_icon_name(const Glib::ustring& iconName, const Gtk::IconSize iconSize)
{
    this->size = iconSize;
    Gtk::Image::set_from_icon_name(iconName, this->size);
}
