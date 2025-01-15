/*
 *  This file is part of RawTherapee.
 *
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

#include <gtkmm.h>

/**
 * A static class in order to handle Hi-DPI.
 *
 * About Cairo size convention (for surface):
 *     Cairo size is expressed in "px" with a 96 DPI (i.e. px per inch) default resolution
 *
 * About Pango size convention (for font):
 *     Pango size can be expressed in two different units:
 *         - Absolute size (i.e. "px")
 *         - Non-absolute size (i.e. "pt"): The default resolution is 72 DPI (i.e. pt per inch). To
 *             convert the size to "px", use the following formula:
 *                 "size in px" = "size in pt" * ("device resolution" / 72)
 *     Notes:
 *         - By default, size is expressed in non-absolute size (i.e. "pt"). Conversion between absolute
 *            and non-absolute size is ensured by Pango.
 *         - On MacOS, font is already scaled by the System library (i.e. "size in px" = "size in pt" * 1.).
 *            Refer to https://gitlab.gnome.org/GNOME/gtk/-/blob/gtk-3-24/gdk/quartz/gdkscreen-quartz.c
 *
 * Hi-DPI implementation according to the OS (source: GDK code):
 *     - Windows: A default DPI of 96 is considered. Current DPI parameter is provided by the OS.
 *        Scale is calculated by (int)("current DPI" / 96). If Scale is greater than 1, DPI is
 *        forced to 96.
 *     - MacOS: Scale is calculated from OS parameters (= "Retina screen width" / "Virtual width").
 *        DPI is forced to 72.
 *     - Linux: DPI is calculated from OS parameter (= 96 * "text-scaling-factor"). Note: "text-scaling-factor"
 *        is different from "device factor".
 */
class RTScalable
{
private:
    static double dpi;
    static int scale;
    static void getDPInScale(const Gtk::Window* window, double &newDPI, int &newScale);

protected:
    static Cairo::RefPtr<Cairo::ImageSurface> loadSurfaceFromIcon(const Glib::ustring &icon_name, const Gtk::IconSize iconSize = Gtk::ICON_SIZE_SMALL_TOOLBAR);
    static Cairo::RefPtr<Cairo::ImageSurface> loadSurfaceFromPNG(const Glib::ustring &fname, const bool is_path = false);
    static Cairo::RefPtr<Cairo::ImageSurface> loadSurfaceFromSVG(const Glib::ustring &fname, const int width = -1, const int height = -1, const bool is_path = false);

public:
    static constexpr double pangoDPI = 72.; // Pango default DPI for "pt" size
    static constexpr double baseDPI = 96.; // Cairo default DPI
    static void init(const Gtk::Window* window);
    static void setDPInScale(const Gtk::Window* window);
    static void setDPInScale(const double newDPI, const int newScale);
    static double getDPI();
    static int getScale();
    static double getGlobalScale();
    static int scalePixelSize(const int pixel_size);
    static double scalePixelSize(const double pixel_size);
};
