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
#include "../rtengine/icons.h"

namespace
{

std::map<std::string, Cairo::RefPtr<Cairo::ImageSurface>> surfaceCache;

}

double RTImage::dpiBack = 0.;
int RTImage::scaleBack = 0;

RTImage::RTImage () {}

RTImage::RTImage (const Glib::ustring& fileName, const Glib::ustring& rtlFileName)
{
    setImage (fileName, rtlFileName);
}

RTImage::RTImage (RTImage &other)
{
    from(&other);
}

RTImage::RTImage (Glib::RefPtr<Gdk::Pixbuf> &pixbuf)
{
    if (pixbuf) {
        set(pixbuf);
    }
}

RTImage::RTImage (Glib::RefPtr<RTImage> &other)
{
    if (other) {
        from(other.get());
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
        printf("RTImage::setDPInScale : New scale = %d & new DPI = %.3f (%.3f asked) -> Reloading all RTImage\n", scaleBack, dpiBack, newDPI);
        updateImages();
    }
}

void RTImage::changeImage (const Glib::ustring& imageName)
{
    clear ();

    auto iterator = surfaceCache.find (imageName);

    if (iterator == surfaceCache.end ()) {
        const auto pixbuf = createFromFile(imageName);
        iterator = surfaceCache.emplace (imageName, pixbuf).first;
    }

    surface = iterator->second;
    set(iterator->second);
}

void RTImage::init()
{
    dpiBack = RTScalable::getDPI();
    scaleBack = RTScalable::getScale();
}

void RTImage::updateImages()
{
    for (auto& entry : surfaceCache) {
        entry.second = createFromFile(entry.first);
    }
}

void RTImage::from(RTImage* other)
{
    if (!other) {
        return;
    }

    if (other->get_pixbuf()) {
        set(other->get_pixbuf());
    } else {
        surface = other->surface;
        set(surface);
    }
}

void RTImage::from(Glib::RefPtr<RTImage> other)
{
    if (other) {
        from (other.get());
    }
}

Cairo::RefPtr<Cairo::ImageSurface> RTImage::createFromFile (const Glib::ustring& fileName)
{
    Cairo::RefPtr<Cairo::ImageSurface> surf;

    try {

        double requestedDPI = getDPI();
        const auto filePath = rtengine::findIconAbsolutePath (fileName, requestedDPI);
        // requestedDPI is now the icon's DPI, not necessarily the one we asked

        if (!filePath.empty ()) {
            surf = Cairo::ImageSurface::create_from_png(filePath);
        }

        if (!surf) {
            return surf; // nullptr
        }

        if (getDPI() > 0. && requestedDPI != -1 && requestedDPI != getDPI()) {
            // scale the bitmap
            printf("Resizing from %.1f to %.1f DPI\n", requestedDPI, getDPI());
            resizeImage(surf, getDPI() / requestedDPI);
        }

        double x=0., y=0.;
        cairo_surface_get_device_scale(surf->cobj(), &x, &y);
        printf("   -> Cairo::ImageSurface is now %dx%d (scale: %.1f)\n", surf->get_width(), surf->get_height(), (float)x);
        if (getScale() == 2) {
            cairo_surface_set_device_scale(surf->cobj(), 0.5, 0.5);
            cairo_surface_get_device_scale(surf->cobj(), &x, &y);
            surf->flush();
            printf("      Cairo::ImageSurface is now %dx%d (scale: %.1f)\n", surf->get_width(), surf->get_height(), (float)x);
        }
    } catch (const Glib::Exception& exception) {
        if (options.rtSettings.verbose) {
            std::cerr << "Failed to load image \"" << fileName << "\": " << exception.what() << std::endl;
        }
    }

    return surf;
}

/*
bool RTImage::on_configure_event(GdkEventConfigure* configure_event)
{
}
*/
