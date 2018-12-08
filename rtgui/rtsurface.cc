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

#include <iostream>

#include "options.h"
#include "../rtengine/icons.h"
#include "rtsurface.h"

namespace
{

std::map<std::string, Cairo::RefPtr<Cairo::ImageSurface>> surfaceCache;

}

double RTSurface::dpiBack = 0.;
int RTSurface::scaleBack = 0;

RTSurface::RTSurface () : RTScalable()
{
    Cairo::RefPtr<Cairo::ImageSurface> imgSurf(new Cairo::ImageSurface(nullptr, false));
    surface = imgSurf;
}

RTSurface::RTSurface(const RTSurface& other) : RTScalable()
{
    surface = other.surface;
}

RTSurface::RTSurface (Glib::ustring fileName, Glib::ustring rtlFileName) : RTScalable()
{
    Cairo::RefPtr<Cairo::ImageSurface> imgSurf(new Cairo::ImageSurface(nullptr, false));
    surface = imgSurf;
    setImage (fileName, rtlFileName);
}

void RTSurface::setImage (Glib::ustring fileName, Glib::ustring rtlFileName)
{
    Glib::ustring imageName;

    if (!rtlFileName.empty() && getDirection() == Gtk::TEXT_DIR_RTL) {
        imageName = rtlFileName;
    } else {
        imageName = fileName;
    }

    changeImage (imageName);
}

void RTSurface::setDPInScale (const double newDPI, const int newScale)
{
    if (getScale() != newScale || (getScale() == 1 && getDPI() != newDPI)) {
        RTScalable::setDPInScale(newDPI, newScale);
        dpiBack = getDPI();
        scaleBack = getScale();
        printf("RTSurface::setDPInScale : New scale = %d & new DPI = %.3f (%.3f asked) -> Reloading all RTSurface\n", scaleBack, dpiBack, newDPI);
        updateImages();
    }
}

void RTSurface::changeImage (Glib::ustring imageName)
{
    auto iterator = surfaceCache.find (imageName);

    if (iterator == surfaceCache.end ()) {
        double requestedDPI = getDPI();
        const auto imagePath = rtengine::findIconAbsolutePath (imageName, requestedDPI);
        surface = Cairo::ImageSurface::create_from_png(imagePath);

        if (getDPI() > 0. && requestedDPI != -1 && requestedDPI != getDPI()) {
            // scale the bitmap
            printf("Resizing from %.1f to %.1f DPI\n", requestedDPI, getDPI());
            resizeImage(surface, getDPI() / requestedDPI);
        }

        // HOMBRE: As of now, GDK_SCALE is forced to 1, so setting the Cairo::ImageSurface scale is not required
        /*
        double x=0., y=0.;
        cairo_surface_get_device_scale(surface->cobj(), &x, &y);
        printf("   -> Cairo::ImageSurface is now %dx%d (scale: %.1f)\n", surface->get_width(), surface->get_height(), (float)x);
        if (getScale() == 2) {
            cairo_surface_set_device_scale(surface->cobj(), 0.5, 0.5);
            cairo_surface_get_device_scale(surface->cobj(), &x, &y);
            surface->flush();
            printf("      Cairo::ImageSurface is now %dx%d (scale: %.1f)\n", surface->get_width(), surface->get_height(), (float)x);
        }
        */

        iterator = surfaceCache.emplace (imageName, surface).first;
    }

    surface = iterator->second;
}

int RTSurface::getWidth() const
{
    return surface ? surface->get_width() : -1;
}

int RTSurface::getHeight() const
{
    return surface ? surface->get_height() : -1;
}

void RTSurface::init()
{
    dpiBack = RTScalable::getDPI();
    scaleBack = RTScalable::getScale();
}

void RTSurface::updateImages()
{
    for (auto& entry : surfaceCache) {
        double imgDPI = getDPI();
        const auto imagePath = rtengine::findIconAbsolutePath (entry.first, imgDPI);
        entry.second = Cairo::ImageSurface::create_from_png(imagePath);
        printf("RTSurface::updateImages : %s\n", imagePath.c_str());
    }
}

void RTSurface::from(Glib::RefPtr<RTSurface> other)
{
    surface = other->surface;
}

bool RTSurface::hasSurface() const
{
    return surface ? true : false;
}
