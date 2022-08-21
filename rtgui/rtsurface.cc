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
    type = RTSurface::InvalidType;
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
        type = RTSurface::IconType;
        name = icon_name;
        icon_size = iconSize;
    }
}

RTSurface::RTSurface(const Glib::ustring &fname) :
    RTSurface()
{
    // Create surface based on file extension
    const auto pos = fname.find_last_of('.');

    if (pos >= 0 && pos < fname.length()) {
        const auto fext = fname.substr(pos + 1, fname.length()).lowercase();

        // Case where fname is a PNG file
        if (fext == "png") {
            // Create surface from PNG file
            surface = RTScalable::loadSurfaceFromPNG(fname);

            if (surface) {
                // Save private parameter
                type = RTSurface::PNGType;
                name = fname;
            }
        }

        // Case where fname is a SVG file
        if (fext == "svg") {
            // Create surface from SVG file
            surface = RTScalable::loadSurfaceFromSVG(fname);

            if (surface) {
                // Save private parameter
                type = RTSurface::SVGType;
                name = fname;
            }
        }
    }
}

int RTSurface::getWidth()
{
    return
        surface
            ? (surface->get_width() / RTScalable::getScale())
            : -1;
}

int RTSurface::getHeight()
{
    return
        surface
            ? (surface->get_height() / RTScalable::getScale())
            : -1;
}

bool RTSurface::hasSurface()
{
    return static_cast<bool>(surface);
}

Cairo::RefPtr<Cairo::ImageSurface> RTSurface::get()
{
    if (dpiBack != RTScalable::getDPI() ||
        scaleBack != RTScalable::getScale()) {
            switch (type) {
                case RTSurface::IconType :
                    surface = RTScalable::loadSurfaceFromIcon(name, icon_size);
                    break;
                case RTSurface::PNGType :
                    surface = RTScalable::loadSurfaceFromPNG(name);
                    break;
                case RTSurface::SVGType :
                    surface = RTScalable::loadSurfaceFromSVG(name);
                    break;
                default :
                    break;
            }

            // Save new DPI and scale
            dpiBack = RTScalable::getDPI();
            scaleBack = RTScalable::getScale();
        }

    return surface;
}

RTPixbuf::RTPixbuf()
{
    // Initialize "back" parameters from RTScalable
    dpiBack = RTPixbuf::getDPI();
    scaleBack = RTPixbuf::getScale();

    // Initialize other private parameters
    type = RTPixbuf::InvalidType;
    name = "";
    icon_size = Gtk::ICON_SIZE_INVALID;
}

RTPixbuf::RTPixbuf(const Glib::ustring &icon_name, const Gtk::IconSize iconSize) :
    RTPixbuf()
{
    // Create surface
    const Cairo::RefPtr<Cairo::ImageSurface> surface = RTScalable::loadSurfaceFromIcon(icon_name, iconSize);

    if (surface) {
        // Downscale surface before creating pixbuf
        const Cairo::RefPtr<Cairo::Surface> _surface = Cairo::Surface::create(surface,
            Cairo::CONTENT_COLOR_ALPHA,
            surface->get_width() / RTScalable::getScale(),
            surface->get_height() / RTScalable::getScale());
        const auto c = Cairo::Context::create(_surface);
        c->scale(1. / static_cast<double>(RTScalable::getScale()), 1. / static_cast<double>(RTScalable::getScale()));
        c->set_source(surface, 0., 0.);
        c->paint();

        // Create pixbuf from surface
        pixbuf = Gdk::Pixbuf::create(_surface,
            0,
            0,
            surface->get_width() / RTScalable::getScale(),
            surface->get_height() / RTScalable::getScale());

        // Save private parameters
        type = RTPixbuf::IconType;
        name = icon_name;
        icon_size = iconSize;
    }
}

RTPixbuf::RTPixbuf(const Glib::ustring &fname) :
    RTPixbuf()
{
    // Create surface based on file extension
    const auto pos = fname.find_last_of('.');

    if (pos >= 0 && pos < fname.length()) {
        const auto fext = fname.substr(pos + 1, fname.length()).lowercase();

        // Case where fname is a PNG file
        if (fext == "png") {
            // Create surface from PNG file
            const Cairo::RefPtr<Cairo::ImageSurface> surface = RTScalable::loadSurfaceFromPNG(fname);

            if (surface) {
                // Downscale surface before creating pixbuf
                const Cairo::RefPtr<Cairo::Surface> _surface = Cairo::Surface::create(surface,
                    Cairo::CONTENT_COLOR_ALPHA,
                    surface->get_width() / RTScalable::getScale(),
                    surface->get_height() / RTScalable::getScale());
                const auto c = Cairo::Context::create(_surface);
                c->scale(1. / static_cast<double>(RTScalable::getScale()), 1. / static_cast<double>(RTScalable::getScale()));
                c->set_source(surface, 0., 0.);
                c->paint();

                // Create pixbuf from surface
                pixbuf = Gdk::Pixbuf::create(_surface,
                    0,
                    0,
                    surface->get_width() / RTScalable::getScale(),
                    surface->get_height() / RTScalable::getScale());

                // Save private parameter
                type = RTPixbuf::PNGType;
                name = fname;
            }
        }

        // Case where fname is a SVG file
        if (fext == "svg") {
            // Create surface from SVG file
            const Cairo::RefPtr<Cairo::ImageSurface> surface = RTScalable::loadSurfaceFromSVG(fname);

            if (surface) {
                // Downscale surface before creating pixbuf
                const Cairo::RefPtr<Cairo::Surface> _surface = Cairo::Surface::create(surface,
                    Cairo::CONTENT_COLOR_ALPHA,
                    surface->get_width() / RTScalable::getScale(),
                    surface->get_height() / RTScalable::getScale());
                const auto c = Cairo::Context::create(_surface);
                c->scale(1. / static_cast<double>(RTScalable::getScale()), 1. / static_cast<double>(RTScalable::getScale()));
                c->set_source(surface, 0., 0.);
                c->paint();

                // Create pixbuf from surface
                pixbuf = Gdk::Pixbuf::create(_surface,
                    0,
                    0,
                    surface->get_width() / RTScalable::getScale(),
                    surface->get_height() / RTScalable::getScale());

                // Save private parameter
                type = RTPixbuf::SVGType;
                name = fname;
            }
        }
    }
}

int RTPixbuf::getWidth()
{
    return
        pixbuf
            ? pixbuf->get_width()
            : -1;
}

int RTPixbuf::getHeight()
{
    return
        pixbuf
            ? pixbuf->get_height()
            : -1;
}

bool RTPixbuf::hasPixbuf()
{
    return static_cast<bool>(pixbuf);
}

Glib::RefPtr<Gdk::Pixbuf> RTPixbuf::get()
{
    if (dpiBack != RTScalable::getDPI() ||
        scaleBack != RTScalable::getScale()) {
            Cairo::RefPtr<Cairo::ImageSurface> surface;

            // Surface needs to be regenerated
            switch (type) {
                case RTPixbuf::IconType :
                    surface = RTScalable::loadSurfaceFromIcon(name, icon_size);
                    pixbuf = Gdk::Pixbuf::create(surface, 0, 0, surface->get_width(), surface->get_height());
                    break;
                case RTPixbuf::PNGType :
                    surface = RTScalable::loadSurfaceFromPNG(name);
                    pixbuf = Gdk::Pixbuf::create(surface, 0, 0, surface->get_width(), surface->get_height());
                    break;
                case RTPixbuf::SVGType :
                    surface = RTScalable::loadSurfaceFromSVG(name);
                    pixbuf = Gdk::Pixbuf::create(surface, 0, 0, surface->get_width(), surface->get_height());
                    break;
                default :
                    break;
            }

            // Save new DPI and scale
            dpiBack = RTScalable::getDPI();
            scaleBack = RTScalable::getScale();
        }

    return pixbuf;
}
