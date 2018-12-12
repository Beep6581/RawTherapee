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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "rtimage.h"

#include <iostream>

#include "options.h"

namespace
{

std::map<std::string, Glib::RefPtr<Gdk::Pixbuf> > pixbufCache;
std::map<std::string, Cairo::RefPtr<Cairo::ImageSurface> > surfaceCache;

}

double RTImage::dpiBack = 0.;
int RTImage::scaleBack = 0;

RTImage::RTImage () {}

RTImage::RTImage (RTImage &other)
{
    dpiBack = other.dpiBack;
    scaleBack = other.scaleBack;
    pixbuf = other.pixbuf;
    surface = other.surface;
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

RTImage::RTImage (Glib::RefPtr<Gdk::Pixbuf> &pixbuf)
{
    if (pixbuf) {
        set(pixbuf);
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
        //printf("RTImage::setDPInScale : New scale = %d & new DPI = %.3f (%.3f asked) -> Reloading all RTImage\n", scaleBack, dpiBack, newDPI);
        updateImages();
    }
}

void RTImage::changeImage (const Glib::ustring& imageName)
{
    clear ();

    if (pixbuf) {
        auto iterator = pixbufCache.find (imageName);
        //printf("changeImage / pixbufCache[%d] : \"%s\" %s!\n", (int)(pixbufCache.size()), imageName.c_str(), iterator == pixbufCache.end () ? "not found" : "found");
        assert(iterator != pixbufCache.end ());
        pixbuf = iterator->second;
        set(iterator->second);
    } else  {  // if no Pixbuf is set, we update or create a Cairo::ImageSurface
        auto iterator = surfaceCache.find (imageName);
        //printf("changeImage / surfaceCache[%d] : \"%s\" %s!\n", (int)(surfaceCache.size()), imageName.c_str(), iterator == surfaceCache.end () ? "not found" : "found");
        if (iterator == surfaceCache.end ()) {
            auto surf = createImgSurfFromFile(imageName);
            iterator = surfaceCache.emplace (imageName, surf).first;
        }
        surface = iterator->second;
        set(iterator->second);
    }
}

void RTImage::init()
{
    dpiBack = RTScalable::getDPI();
    scaleBack = RTScalable::getScale();
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
    Glib::RefPtr<Gdk::Pixbuf> pixbuf = Gdk::Pixbuf::create(imgSurf, 0, 0, imgSurf->get_width(), imgSurf->get_height());
    return pixbuf;
}

Cairo::RefPtr<Cairo::ImageSurface> RTImage::createImgSurfFromFile (const Glib::ustring& fileName)
{
    Cairo::RefPtr<Cairo::ImageSurface> surf;

    //printf("Creating \"%s\"\n", fileName.c_str());
    try {
        surf = loadImage(fileName, getTweakedDPI());

        // HOMBRE: As of now, GDK_SCALE is forced to 1, so setting the Cairo::ImageSurface scale is not required
        /*
        double x=0., y=0.;
        cairo_surface_get_device_scale(surf->cobj(), &x, &y);
        //printf("   -> Cairo::ImageSurface is now %dx%d (scale: %.1f)\n", surf->get_width(), surf->get_height(), (float)x);
        if (getScale() == 2) {
            cairo_surface_set_device_scale(surf->cobj(), 0.5, 0.5);
            cairo_surface_get_device_scale(surf->cobj(), &x, &y);
            surf->flush();
            //printf("      Cairo::ImageSurface is now %dx%d (scale: %.1f)\n", surf->get_width(), surf->get_height(), (float)x);
        }
        */
    } catch (const Glib::Exception& exception) {
        if (options.rtSettings.verbose) {
            std::cerr << "Failed to load image \"" << fileName << "\": " << exception.what() << std::endl;
        }
    }

    return surf;
}
