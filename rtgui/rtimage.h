/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include "rtscalable.h"

/**
 * @brief A derived class of Gtk::Image in order to handle theme-related icon sets.
 */
class RTImage final : public Gtk::Image, public RTScalable
{
    static double dpiBack; // used to keep track of master dpi change
    static int scaleBack;  // used to keep track of master scale change
    //bool on_configure_event(GdkEventConfigure* configure_event);

protected:
    Cairo::RefPtr<Cairo::ImageSurface> surface;
    Glib::RefPtr<Gdk::Pixbuf> pixbuf;

public:
    RTImage ();
    RTImage (RTImage &other);
    explicit RTImage (Glib::RefPtr<Gdk::Pixbuf> &pixbuf);
    explicit RTImage (Cairo::RefPtr<Cairo::ImageSurface> &surf);
    explicit RTImage(Cairo::RefPtr<Cairo::ImageSurface> other);
    explicit RTImage (Glib::RefPtr<RTImage> &other);
    explicit RTImage (const Glib::ustring& fileName, const Glib::ustring& rtlFileName = Glib::ustring());

    void setImage (const Glib::ustring& fileName, const Glib::ustring& rtlFileName = Glib::ustring());
    void changeImage (const Glib::ustring& imageName);
    Cairo::RefPtr<Cairo::ImageSurface> get_surface();
    int get_width();
    int get_height();


    static void init();
    static void cleanup(bool all = false);
    static void updateImages ();
    static void setDPInScale (const double newDPI, const int newScale);
    static void setScale (const int newScale);

    static Glib::RefPtr<Gdk::Pixbuf> createPixbufFromFile (const Glib::ustring& fileName);
    static Cairo::RefPtr<Cairo::ImageSurface> createImgSurfFromFile (const Glib::ustring& fileName);

};
