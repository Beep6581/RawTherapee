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

#include "rtscalable.h"

/**
 * @brief A derived class of Gtk::Image in order to handle theme-related icon sets.
 */
class RTSurface :
    public RTScalable
{
public:
    RTSurface();
    explicit RTSurface(const Glib::ustring& fileName, const Glib::ustring& rtlFileName = {});

    void setImage(const Glib::ustring& fileName, const Glib::ustring& rtlFileName = {});

    int getWidth() const;
    int getHeight() const;
    bool hasSurface() const;

    Cairo::RefPtr<const Cairo::ImageSurface> get() const;
    const Cairo::RefPtr<Cairo::ImageSurface>& get();

    static void init();
    static void updateImages();
    static void setDPInScale(double newDPI, int newScale);

private:
    void changeImage(const Glib::ustring& imageName);

    static double dpiBack; // used to keep track of master dpi change
    static int scaleBack;  // used to keep track of master scale change
    Cairo::RefPtr<Cairo::ImageSurface> surface;
};
