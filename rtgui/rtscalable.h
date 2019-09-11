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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <gtkmm/image.h>

/**
 * @brief A master class for derived class of Gtk::Image in order to handle theme-related icon sets.
 */
class RTScalable
{
    static double dpi;
    static int scale;
    static Gtk::TextDirection direction;  // cached value for text-direction
    static void deleteDir(const Glib::ustring& path);

protected:
    static void setDPInScale (const double newDPI, const int newScale);
    static Cairo::RefPtr<Cairo::ImageSurface> loadImage(const Glib::ustring &fname, double dpi);
    static Gtk::TextDirection getDirection();


public:

#ifdef __APPLE__
    static constexpr double baseDPI = 72.;
    static constexpr double baseHiDPI = 144.;
    static constexpr int    baseFontSize = 12;
#else
    static constexpr double baseDPI = 96.;
    static constexpr double baseHiDPI = 192.;
    static constexpr int    baseFontSize = 9;
#endif

    static void init(Gtk::Window *window);
    static void cleanup(bool all = false);
    static double getDPI ();
    static double getTweakedDPI ();   // The returned value is tweaked DPI to adapt to main the font size. Maybe not an ideal solution.
    static int getScale ();
};
