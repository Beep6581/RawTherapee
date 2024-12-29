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

#include <iostream>

#include "rtsurface.h"

#include "options.h"

RTSurface::RTSurface() :
    surface(new Cairo::ImageSurface(nullptr, false))
{
    // Initialize "back" parameters from RTScalable
    dpiBack = RTScalable::getDPI();
    scaleBack = RTScalable::getScale();

    // Initialize other private parameters
    type = RTSurfaceType::InvalidType;
    name = "";
    icon_size = Gtk::ICON_SIZE_INVALID;
}

RTSurface::RTSurface(const Glib::ustring &icon_name, const Gtk::IconSize iconSize) :
    RTSurface()
{
    // Create surface
    surface = RTScalable::loadSurfaceFromIcon(icon_name, iconSize);

    if (surface) {
        // Save private parameters
        type = RTSurfaceType::IconType;
        name = icon_name;
        icon_size = iconSize;
    }
}

RTSurface::RTSurface(const Glib::ustring &fname) :
    RTSurface()
{
    // Create surface based on file extension
    const auto pos = fname.find_last_of('.');

    if (pos < fname.length()) {
        const auto fext = fname.substr(pos + 1, fname.length()).lowercase();

        // Case where fname is a PNG file
        if (fext == "png") {
            // Create surface from PNG file
            surface = RTScalable::loadSurfaceFromPNG(fname);

            if (surface) {
                // Save private parameter
                type = RTSurfaceType::PNGType;
                name = fname;
            }
        }

        // Case where fname is a SVG file
        if (fext == "svg") {
            // Create surface from SVG file
            surface = RTScalable::loadSurfaceFromSVG(fname);

            if (surface) {
                // Save private parameter
                type = RTSurfaceType::SVGType;
                name = fname;
            }
        }
    }
}

int RTSurface::getWidth()
{
    int w, h;

    if (hasSurface()) {
        switch (type) {
            case RTSurfaceType::IconType:
                // Get width from Gtk::IconSize
                if (!Gtk::IconSize::lookup(icon_size, w, h)) { // Size in invalid
                    w = h = -1; // Invalid case
                }

                return w;

            case RTSurfaceType::PNGType:
                // Directly return surface width
                return surface->get_width();

            case RTSurfaceType::SVGType:
                // Returned size shall consider the scaling
                return (surface->get_width() / RTScalable::getScale());

            case RTSurfaceType::InvalidType:
            default:
                // Invalid case
                return -1;
        }
    } else {
        // Invalid case
        return -1;
    }
}

int RTSurface::getHeight()
{
    int w, h;

    if (hasSurface()) {
        switch (type) {
            case RTSurfaceType::IconType:
                // Get width from Gtk::IconSize
                if (!Gtk::IconSize::lookup(icon_size, w, h)) { // Size in invalid
                    w = h = -1; // Invalid case
                }

                return h;

            case RTSurfaceType::PNGType:
                // Directly return surface width
                return surface->get_height();

            case RTSurfaceType::SVGType:
                // Returned size shall consider the scaling
                return (surface->get_height() / RTScalable::getScale());

            case RTSurfaceType::InvalidType:
            default:
                // Invalid case
                return -1;
        }
    } else {
        // Invalid case
        return -1;
    }
}

bool RTSurface::hasSurface()
{
    return static_cast<bool>(surface);
}

Cairo::RefPtr<Cairo::ImageSurface> RTSurface::get()
{
    if (dpiBack != RTScalable::getDPI() ||
        scaleBack != RTScalable::getScale()) {
            updateSurface();

            // Save new DPI and scale
            dpiBack = RTScalable::getDPI();
            scaleBack = RTScalable::getScale();
        }

    return surface;
}

void RTSurface::updateSurface()
{
    // Update surface based on the scale
    switch (type) {
        case RTSurfaceType::IconType :
            surface = RTScalable::loadSurfaceFromIcon(name, icon_size);
            break;
        case RTSurfaceType::PNGType :
            surface = RTScalable::loadSurfaceFromPNG(name);
            break;
        case RTSurfaceType::SVGType :
            surface = RTScalable::loadSurfaceFromSVG(name);
            break;
        default :
            break;
    }
}
