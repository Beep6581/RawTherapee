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

#include <iostream>

#include "rtsurface.h"

#include "options.h"

namespace
{

using SurfaceCache = std::map<std::string, Cairo::RefPtr<Cairo::ImageSurface>>;

SurfaceCache surfaceCache;

}

RTSurface::RTSurface() :
    surface(new Cairo::ImageSurface(nullptr, false))
{
}

RTSurface::RTSurface(const Glib::ustring& fileName, const Glib::ustring& rtlFileName) :
    RTSurface()
{
    setImage(fileName, rtlFileName);
}

void RTSurface::setImage(const Glib::ustring& fileName, const Glib::ustring& rtlFileName)
{
    const Glib::ustring& imageName =
        !rtlFileName.empty() && getDirection() == Gtk::TEXT_DIR_RTL
            ? rtlFileName
            : fileName;

    changeImage (imageName);
}

int RTSurface::getWidth() const
{
    return
        surface
            ? surface->get_width()
            : -1;
}

int RTSurface::getHeight() const
{
    return
        surface
            ? surface->get_height()
            : -1;
}

bool RTSurface::hasSurface() const
{
    return static_cast<bool>(surface);
}

Cairo::RefPtr<const Cairo::ImageSurface> RTSurface::get() const
{
    return surface;
}

const Cairo::RefPtr<Cairo::ImageSurface>& RTSurface::get()
{
    return surface;
}

void RTSurface::init()
{
    dpiBack = getDPI();
    scaleBack = getScale();
}

void RTSurface::updateImages()
{
    const double tweakedDpi = getTweakedDPI();

    for (auto& entry : surfaceCache) {
        entry.second = loadImage(entry.first, tweakedDpi);
    }
}

void RTSurface::setDPInScale(const double newDPI, const int newScale)
{
    if (
        getScale() != newScale
        || (
            getScale() == 1
            && getDPI() != newDPI
        )
    ) {
        setDPInScale(newDPI, newScale);
        dpiBack = getDPI();
        scaleBack = getScale();

        updateImages();
    }
}

void RTSurface::changeImage(const Glib::ustring& imageName)
{
    const SurfaceCache::const_iterator iterator = surfaceCache.find(imageName);

    if (iterator != surfaceCache.end()) {
        surface = iterator->second;
    } else {
        surface = loadImage(imageName, getTweakedDPI());

        // HOMBRE: As of now, GDK_SCALE is forced to 1, so setting the Cairo::ImageSurface scale is not required
        //         Anyway, this might be of use one day
        /*
        double x=0., y=0.;
        cairo_surface_get_device_scale(surface->cobj(), &x, &y);
        if (getScale() == 2) {
            cairo_surface_set_device_scale(surface->cobj(), 0.5, 0.5); // Not sure if it should be 0.5 or 2.0 here !
            surface->flush();
        }
        */

        surfaceCache.emplace(imageName, surface);
    }
}

double RTSurface::dpiBack = 0.;

int RTSurface::scaleBack = 0;
