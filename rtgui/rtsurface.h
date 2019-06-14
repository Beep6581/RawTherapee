/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

#include <gtkmm/image.h>
#include "rtscalable.h"

/**
 * @brief A derived class of Gtk::Image in order to handle theme-related icon sets.
 */
class RTSurface : public RTScalable
{

private:

    static double dpiBack; // used to keep track of master dpi change
    static int scaleBack;  // used to keep track of master scale change
    Cairo::RefPtr<Cairo::ImageSurface> surface;
    void changeImage (Glib::ustring imageName);

public:

    RTSurface ();
    RTSurface (const RTSurface&  other);
    RTSurface (Glib::ustring fileName, Glib::ustring rtlFileName = Glib::ustring());

    void setImage (Glib::ustring fileName, Glib::ustring rtlFileName = Glib::ustring());
    int getWidth() const;
    int getHeight() const;
    bool hasSurface() const;

    const Cairo::RefPtr<Cairo::ImageSurface>& get() const;

    static void init();
    static void updateImages ();
    static void setDPInScale (const double newDPI, const int newScale);
    static void setScale (const int newScale);

    void from(Glib::RefPtr<RTSurface> other);
};
