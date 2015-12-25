/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2011 Jean-Christophe FRISCH <natureh@free.fr>
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

std::vector<Glib::ustring> imagePaths;
std::map<Glib::ustring, Glib::RefPtr<Gdk::Pixbuf>> pixbufCache;

bool loadIconSet(const Glib::ustring& iconSet)
{
    try {

        Glib::KeyFile keyFile;
        keyFile.load_from_file (iconSet);

        auto iconSetDir = keyFile.get_string ("General", "Iconset");

        if (!iconSetDir.empty ()) {
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir, "actions"));
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir));
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir, "devices"));
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir, "places"));
        }

        iconSetDir = keyFile.get_string ("General", "FallbackIconset");

        if (!iconSetDir.empty ()) {
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir, "actions"));
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir));
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir, "devices"));
            imagePaths.push_back (Glib::build_filename (argv0, "images", iconSetDir, "places"));
        }

        return true;

    } catch (const Glib::Exception& exception) {

        if (options.rtSettings.verbose) {
            std::cerr << "Failed to load icon set \"" << iconSet << "\": " << exception.what() << std::endl;
        }

        return false;

    }
}

}

RTImage::RTImage (const Glib::ustring& fileName, const Glib::ustring& rtlFileName) : Gtk::Image()
{
    Glib::ustring imageName;

    if (!rtlFileName.empty () && get_direction () == Gtk::TEXT_DIR_RTL) {
        imageName = rtlFileName;
    } else {
        imageName = fileName;
    }

    changeImage (imageName);
}

void RTImage::changeImage (const Glib::ustring& imageName)
{
    clear ();

    auto iterator = pixbufCache.find (imageName);

    if (iterator == pixbufCache.end ()) {
        const auto imagePath = findIconAbsolutePath (imageName);
        const auto pixbuf = Gdk::Pixbuf::create_from_file (imagePath);

        iterator = pixbufCache.emplace (imageName, pixbuf).first;
    }

    set(iterator->second);
}

void RTImage::updateImages()
{
    for (auto& entry : pixbufCache) {
        const auto imagePath = findIconAbsolutePath (entry.first);
        entry.second = Gdk::Pixbuf::create_from_file (imagePath);
    }
}

Glib::ustring RTImage::findIconAbsolutePath (const Glib::ustring& iconName)
{
    try {

        for (const auto& imagePath : imagePaths) {
            const auto iconPath = Glib::build_filename(imagePath, iconName);

            if (Glib::file_test(iconPath, Glib::FILE_TEST_IS_REGULAR)) {
                return iconPath;
            }
        }

    } catch(const Glib::Exception&) {}

    if (options.rtSettings.verbose) {
        std::cerr << "Icon \"" << iconName << "\" could not be found!" << std::endl;
    }

    return Glib::ustring();
}

void RTImage::setPaths (const Options& options)
{
    Glib::ustring iconSet;

    // Either use the system icon set or the one specified in the options.
    if (options.useSystemTheme) {
        iconSet = Glib::build_filename (argv0, "themes", "system.iconset");
    } else {
        iconSet = Glib::build_filename (argv0, "themes", options.theme + ".iconset");
    }

    imagePaths.clear ();

    if (!loadIconSet (iconSet)) {
        // If the preferred icon set is unavailable, fall back to the default icon set.
        loadIconSet (Glib::build_filename (argv0, "themes", "Default.iconset"));
    }

    // The images folder is the second fallback solution.
    imagePaths.push_back (Glib::build_filename(argv0, "images"));
}

Glib::RefPtr<Gdk::Pixbuf> RTImage::createFromFile (const Glib::ustring& fileName)
{
    Glib::RefPtr<Gdk::Pixbuf> pixbuf;

    try {

        const auto filePath = findIconAbsolutePath (fileName);

        if (!filePath.empty ()) {
            pixbuf = Gdk::Pixbuf::create_from_file (filePath);
        }

    } catch (const Glib::Exception& exception) {

        if (options.rtSettings.verbose) {
            std::cerr << "Failed to load image \"" << fileName << "\": " << exception.what() << std::endl;
        }

    }

    return pixbuf;
}

Cairo::RefPtr<Cairo::ImageSurface> RTImage::createFromPng (const Glib::ustring& fileName)
{
    Cairo::RefPtr<Cairo::ImageSurface> surface;

    try {

        const auto filePath = findIconAbsolutePath (fileName);

        if (!filePath.empty()) {
            surface = Cairo::ImageSurface::create_from_png (Glib::locale_from_utf8 (filePath));
        }

    } catch (const Glib::Exception& exception) {

        if (options.rtSettings.verbose) {
            std::cerr << "Failed to load PNG \"" << fileName << "\": " << exception.what() << std::endl;
        }

    }

    return surface;
}


