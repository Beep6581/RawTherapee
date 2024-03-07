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

#include "rtscalable.h"

/**
 * @brief A custom class in order to handle Hi-DPI surface.
 */
class RTSurface final : public RTScalable
{
public:
    enum class RTSurfaceType {
        InvalidType,
        IconType,
        PNGType,
        SVGType
    };

private:
    double dpiBack; // Used to identify dpi change
    int scaleBack;  // Used to identify scale change
    RTSurfaceType type;
    Glib::ustring name;
    Gtk::IconSize icon_size;
    Cairo::RefPtr<Cairo::ImageSurface> surface;

public:
    RTSurface();
    explicit RTSurface(const Glib::ustring &icon_name, const Gtk::IconSize iconSize);
    explicit RTSurface(const Glib::ustring &fname);

    int getWidth();
    int getHeight();
    bool hasSurface();

    Cairo::RefPtr<Cairo::ImageSurface> get();
    void updateSurface();
};
