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

#include "rtimage.h"

#include <cassert>
#include <iostream>

#include "../rtengine/settings.h"

namespace
{

std::map<std::string, Glib::RefPtr<Gdk::Pixbuf> > pixbufCache;
std::map<std::string, Cairo::RefPtr<Cairo::ImageSurface> > surfaceCache;

}

double RTImage::dpiBack = 0.;
int RTImage::scaleBack = 0;

RTImage::RTImage () {}

RTImage::RTImage (RTImage &other) : surface(other.surface), pixbuf(other.pixbuf)
{
    if (pixbuf) {
        set(pixbuf);
    } else if (surface) {
        set(surface);
    }
}

RTImage::RTImage (const Glib::ustring& fileName, const Glib::ustring& rtlFileName) : Gtk::Image()
{
    setImage (fileName, rtlFileName);
}

RTImage::RTImage (Glib::RefPtr<Gdk::Pixbuf> &pbuf)
{
    if (surface) {
        surface.clear();
    }
    if (pbuf) {
        set(pbuf);
        this->pixbuf = pbuf;
    }
}

RTImage::RTImage (Cairo::RefPtr<Cairo::ImageSurface> &surf)
{
    if (pixbuf) {
        pixbuf.clear();
    }
    if (surf) {
        set(surf);
        surface = surf;
    }
}

RTImage::RTImage (Glib::RefPtr<RTImage> &other)
{
    if (other) {
        if (other->get_surface()) {
            surface = other->get_surface();
            set(surface);
        } else {
            pixbuf = other->get_pixbuf();
            set(pixbuf);
        }
    }
}

void RTImage::setImage (const Glib::ustring& fileName, const Glib::ustring& rtlFileName)
{
    Glib::ustring imageName;

    if (!rtlFileName.empty() && getDirection() == Gtk::TEXT_DIR_RTL) {
        imageName = rtlFileName;
    } else {
        imageName = fileName;
    }

    changeImage (imageName);
}

/*
 * On windows, if scale = 2, the dpi is non significant, i.e. should be considered = 192
 */
void RTImage::setDPInScale (const double newDPI, const int newScale)
{
    if (scaleBack != newScale || (scaleBack == 1 && dpiBack != newDPI)) {
        RTScalable::setDPInScale(newDPI, newScale);
        dpiBack = getDPI();
        scaleBack = getScale();
        updateImages();
    }
}

void RTImage::changeImage (const Glib::ustring& imageName)
{
    clear ();

    if (imageName.empty()) {
        return;
    }

    if (pixbuf) {
        auto iterator = pixbufCache.find (imageName);
        assert(iterator != pixbufCache.end ());
        pixbuf = iterator->second;
        set(iterator->second);
    } else  {  // if no Pixbuf is set, we update or create a Cairo::ImageSurface
        auto iterator = surfaceCache.find (imageName);
        if (iterator == surfaceCache.end ()) {
            auto surf = createImgSurfFromFile(imageName);
            iterator = surfaceCache.emplace (imageName, surf).first;
        }
        surface = iterator->second;
        set(iterator->second);
    }
}

Cairo::RefPtr<Cairo::ImageSurface> RTImage::get_surface()
{
    return surface;
}

int RTImage::get_width()
{
    if (surface) {
        return surface->get_width();
    }
    if (pixbuf) {
        return pixbuf->get_width();
    }
    return -1;
}

int RTImage::get_height()
{
    if (surface) {
        return surface->get_height();
    }
    if (pixbuf) {
        return pixbuf->get_height();
    }
    return -1;
}

void RTImage::init()
{
    dpiBack = RTScalable::getDPI();
    scaleBack = RTScalable::getScale();
}

void RTImage::cleanup(bool all)
{
    for (auto& entry : pixbufCache) {
        entry.second.reset();
    }
    for (auto& entry : surfaceCache) {
        entry.second.clear();
    }
    RTScalable::cleanup(all);
}

void RTImage::updateImages()
{
    for (auto& entry : pixbufCache) {
        entry.second = createPixbufFromFile(entry.first);
    }
    for (auto& entry : surfaceCache) {
        entry.second = createImgSurfFromFile(entry.first);
    }
}

Glib::RefPtr<Gdk::Pixbuf> RTImage::createPixbufFromFile (const Glib::ustring& fileName)
{
    Cairo::RefPtr<Cairo::ImageSurface> imgSurf = createImgSurfFromFile(fileName);
    return Gdk::Pixbuf::create(imgSurf, 0, 0, imgSurf->get_width(), imgSurf->get_height());
}

Cairo::RefPtr<Cairo::ImageSurface> RTImage::createImgSurfFromFile (const Glib::ustring& fileName)
{
    Cairo::RefPtr<Cairo::ImageSurface> surf;

    try {
        surf = loadImage(fileName, getTweakedDPI());

        // HOMBRE: As of now, GDK_SCALE is forced to 1, so setting the Cairo::ImageSurface scale is not required
        /*
        double x=0., y=0.;
        cairo_surface_get_device_scale(surf->cobj(), &x, &y);
        if (getScale() == 2) {
            cairo_surface_set_device_scale(surf->cobj(), 0.5, 0.5);
            cairo_surface_get_device_scale(surf->cobj(), &x, &y);
            surf->flush();
        }
        */
    } catch (const Glib::Exception& exception) {
        if (rtengine::settings->verbose) {
            std::cerr << "Failed to load image \"" << fileName << "\": " << exception.what() << std::endl;
        }
    }

    return surf;
}
