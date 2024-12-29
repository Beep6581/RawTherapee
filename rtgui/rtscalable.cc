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

#include "rtscalable.h"

#include <iostream>
#include <librsvg/rsvg.h>

#include "rtengine/settings.h"
#include "guiutils.h"

extern Glib::ustring argv0;

// Default static parameter values
double RTScalable::dpi = 96.;
int RTScalable::scale = 1;

void RTScalable::getDPInScale(const Gtk::Window* window, double &newDPI, int &newScale)
{
    if (window) {
        const auto screen = window->get_screen();
        newDPI = screen->get_resolution(); // Get DPI retrieved from the OS

        if (window->get_scale_factor() > 0) {
             // Get scale factor associated to the window
            newScale = window->get_scale_factor();
        } else {
            newScale = 1; // Default minimum value of 1 as scale is used to scale surface
        }
    }
}

Cairo::RefPtr<Cairo::ImageSurface> RTScalable::loadSurfaceFromIcon(const Glib::ustring &iconName, const Gtk::IconSize iconSize)
{
    GThreadLock lock; // All icon theme access or image access on separate thread HAVE to be protected

    Cairo::RefPtr<Cairo::ImageSurface> surf; // Create Cairo::RefPtr<Cairo::ImageSurface> nullptr

    // Get icon theme
    const auto theme = Gtk::IconTheme::get_default();

    // Get pixel size from Gtk::IconSize
    int wSize, hSize;

    if (!Gtk::IconSize::lookup(iconSize, wSize, hSize)) { // Size in invalid
        wSize = hSize = 16; // Set to a default size of 16px (i.e. Gtk::ICON_SIZE_SMALL_TOOLBAR one)
    }

    // Get scale based on DPI and scale
    // Note: hSize not used because icon are considered squared
    const int size = wSize;

    // Looking for corresponding icon (if existing)
    const auto iconInfo = theme->lookup_icon(iconName, size);

    if (!iconInfo) {
        std::cerr << "Failed to load icon \"" << iconName << "\" for size " << size << "px" << std::endl;

        return surf;
    }

    const auto iconPath = iconInfo.get_filename();

    if (iconPath.empty()) {
        std::cerr << "Failed to load icon \"" << iconName << "\" for size " << size << "px" << std::endl;

        return surf;
    }

    // Create surface from corresponding icon
    const auto pos = iconPath.find_last_of('.');

    if (pos < iconPath.length()) {
        const auto fext = iconPath.substr(pos + 1, iconPath.length()).lowercase();

        // Case where iconPath is a PNG file
        if (fext == "png") {
            // Create surface from PNG file
            surf = RTScalable::loadSurfaceFromPNG(iconPath, true);
        }

        // Case where iconPath is a SVG file
        if (fext == "svg") {
            // Create surface from SVG file
            surf = RTScalable::loadSurfaceFromSVG(iconPath, size, size, true);
        }
    }

    return surf;
}

Cairo::RefPtr<Cairo::ImageSurface> RTScalable::loadSurfaceFromPNG(const Glib::ustring &fname, const bool is_path)
{
    GThreadLock lock; // All icon theme access or image access on separate thread HAVE to be protected

    Cairo::RefPtr<Cairo::ImageSurface> surf; // Create Cairo::RefPtr<Cairo::ImageSurface> nullptr

    Glib::ustring path;

    if (is_path) {
        // Directly use fname as a path
        path = fname;
    } else {
        // Look for PNG file in "images" folder
        Glib::ustring imagesFolder = Glib::build_filename(argv0, "images");
        path = Glib::build_filename(imagesFolder, fname);
    }

    // Create surface from PNG file if file exist
    if (Glib::file_test(path.c_str(), Glib::FILE_TEST_EXISTS)) {
        surf = Cairo::ImageSurface::create_from_png(path);
    } else {
        std::cerr << "Failed to load PNG file \"" << fname << "\"" << std::endl;
    }

    return surf;
}

Cairo::RefPtr<Cairo::ImageSurface> RTScalable::loadSurfaceFromSVG(const Glib::ustring &fname, const int width, const int height, const bool is_path)
{
    GThreadLock lock; // All icon theme access or image access on separate thread HAVE to be protected

    Cairo::RefPtr<Cairo::ImageSurface> surf; // Create Cairo::RefPtr<Cairo::ImageSurface> nullptr

    Glib::ustring path;

    if (is_path) {
        // Directly use fname as a path
        path = fname;
    } else {
        // Look for SVG file in "images" folder
        Glib::ustring imagesFolder = Glib::build_filename(argv0, "images");
        path = Glib::build_filename(imagesFolder, fname);
    }

    // Create surface from SVG file if file exist
    if (Glib::file_test(path.c_str(), Glib::FILE_TEST_EXISTS)) {
        // Read content of SVG file
        std::string svgFile;
        try {
            svgFile = Glib::file_get_contents(path);
        }
        catch (Glib::FileError &err) {
            std::cerr << "Failed to load SVG file \"" << fname << "\": " << err.what() << std::endl;
            return surf;
        }

        // Create surface with librsvg library
        GError* error = nullptr;
        RsvgHandle* handle = rsvg_handle_new_from_data((unsigned const char*)svgFile.c_str(), svgFile.length(), &error);

        if (error) {
            std::cerr << "Failed to load SVG file \"" << fname << "\": " << std::endl
                      << Glib::ustring(error->message) << std::endl;
            free(error);
            return surf;
        }

        int w, h;

        if (width == -1 || height == -1) {
            // Use SVG image natural width and height
            double _w, _h;
            const bool has_dim = rsvg_handle_get_intrinsic_size_in_pixels(handle, &_w, &_h); // Get SVG image dimensions
            if (has_dim) {
                w = std::ceil(_w);
                h = std::ceil(_h);
            } else {
                w = h = 16; // Set to a default size of 16px (i.e. Gtk::ICON_SIZE_SMALL_TOOLBAR one)
            }
        } else {
            // Use given width and height
            w = width;
            h = height;
        }

        // Create an upscaled surface to avoid blur effect
        surf = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32,
            w * RTScalable::getScale(),
            h * RTScalable::getScale());

        // Render (and erase with) default surface background
        Cairo::RefPtr<Cairo::Context> c = Cairo::Context::create(surf);
        c->set_source_rgba (0., 0., 0., 0.);
        c->set_operator (Cairo::OPERATOR_CLEAR);
        c->paint();

        // Render upscaled surface based on SVG image
        error = nullptr;
        RsvgRectangle rect = {
            .x = 0.,
            .y = 0.,
            .width = static_cast<double>(w * RTScalable::getScale()),
            .height = static_cast<double>(h * RTScalable::getScale())
        };
        c->set_operator (Cairo::OPERATOR_OVER);
        const bool success = rsvg_handle_render_document(handle, c->cobj(), &rect, &error);

        if (!success && error) {
            std::cerr << "Failed to load SVG file \"" << fname << "\": " << std::endl
                      << Glib::ustring(error->message) << std::endl;
            free(error);
            return surf;
        }

        g_object_unref(handle);

        // Set device scale to avoid blur effect
        cairo_surface_set_device_scale(surf->cobj(),
            static_cast<double>(RTScalable::getScale()),
            static_cast<double>(RTScalable::getScale()));
    } else {
        std::cerr << "Failed to load SVG file \"" << fname << "\"" << std::endl;
    }

    return surf;
}

void RTScalable::init(const Gtk::Window* window)
{
    // Retrieve DPI and Scale paremeters from OS
    getDPInScale(window, dpi, scale);
}

void RTScalable::setDPInScale (const Gtk::Window* window)
{
    getDPInScale(window, dpi, scale);
}

void RTScalable::setDPInScale (const double newDPI, const int newScale)
{
    dpi = newDPI;
    scale = newScale;
}

double RTScalable::getDPI ()
{
    return dpi;
}

int RTScalable::getScale ()
{
    return scale;
}

double RTScalable::getGlobalScale()
{
    return (RTScalable::getDPI() / RTScalable::baseDPI);
}

int RTScalable::scalePixelSize(const int pixel_size)
{
    const double s = getGlobalScale();
    return static_cast<int>(pixel_size * s + 0.5); // Rounded scaled size
}

double RTScalable::scalePixelSize(const double pixel_size)
{
    const double s = getGlobalScale();
    return (pixel_size * s);
}
