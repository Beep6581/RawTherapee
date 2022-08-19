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
#pragma once

#include <gtkmm/image.h>

/**
 * @brief A derived class of Gtk::Image in order to handle theme-related icon sets.
 */
class RTImage final : public Gtk::Image
{
private:
    Gtk::IconSize size;

public:
    RTImage ();
    explicit RTImage (const Glib::ustring& iconName, const Gtk::IconSize iconSize = Gtk::ICON_SIZE_SMALL_TOOLBAR);

    void set_from_icon_name(const Glib::ustring& iconName);
    void set_from_icon_name(const Glib::ustring& iconName, const Gtk::IconSize iconSize);
};
