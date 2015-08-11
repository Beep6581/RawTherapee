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
#include "../rtengine/safegtk.h"

extern Glib::ustring argv0;
extern Options options;

std::vector<Glib::ustring> imagesPaths;
std::map<Glib::ustring, Glib::RefPtr<Gdk::Pixbuf> > pixBufMap; // List of image buffers in order to live update them on theme switch and to avoid a lot of file accesses

/*
 * RTImage is a derived class of Gtk::Image, in order to handle theme related iconsets
 */
RTImage::RTImage(Glib::ustring fileName, Glib::ustring rtlFileName) : Gtk::Image()
{
    Glib::ustring mapKey;

    if (rtlFileName.length()) {
        if (get_direction() == Gtk::TEXT_DIR_RTL) {
            mapKey = rtlFileName;
        } else {
            mapKey = fileName;
        }
    } else {
        mapKey = fileName;
    }

    std::map<Glib::ustring, Glib::RefPtr<Gdk::Pixbuf> >::iterator it;
    it = pixBufMap.find(mapKey);

    if (it != pixBufMap.end()) {
        set(it->second);
    } else {
        Glib::RefPtr<Gdk::Pixbuf> tempPixPuf = Gdk::Pixbuf::create_from_file(findIconAbsolutePath(mapKey));
        pixBufMap.insert(std::pair<Glib::ustring, Glib::RefPtr<Gdk::Pixbuf> >(mapKey, tempPixPuf));
        set(tempPixPuf);
    }
}

void RTImage::updateImages()
{
    std::map<Glib::ustring, Glib::RefPtr<Gdk::Pixbuf> >::iterator it;

    for (it = pixBufMap.begin(); it != pixBufMap.end(); ++it) {
        Glib::ustring fullPath = findIconAbsolutePath(it->first);
        it->second = Gdk::Pixbuf::create_from_file(fullPath);
    }
}

// DONE (was TODO: Maybe this could be optimized: in order to avoid looking up for an icon file in the filesystem on each popupmenu selection, maybe we could find a way to copy the image data from another RTImage)
void RTImage::changeImage(Glib::ustring &newImage)
{
    clear();
    std::map<Glib::ustring, Glib::RefPtr<Gdk::Pixbuf> >::iterator it;
    it = pixBufMap.find(newImage);

    if (it != pixBufMap.end()) {
        set(it->second);
    } else {
        Glib::ustring fullPath = findIconAbsolutePath(newImage);
        Glib::RefPtr<Gdk::Pixbuf> tempPixPuf = Gdk::Pixbuf::create_from_file(fullPath);
        pixBufMap.insert(std::pair<Glib::ustring, Glib::RefPtr<Gdk::Pixbuf> >(newImage, tempPixPuf));
        set(tempPixPuf);
    }
}

Glib::ustring RTImage::findIconAbsolutePath(const Glib::ustring &iconFName)
{
    Glib::ustring path;

    for (unsigned int i = 0; i < imagesPaths.size(); i++) {
        path = Glib::build_filename(imagesPaths[i], iconFName);

        if (safe_file_test(path, Glib::FILE_TEST_EXISTS)) {
            return path;
        }
    }

    printf("\"%s\" not found!\n", iconFName.c_str());
    return "";
}

void RTImage::setPaths(Options &opt)
{
    Glib::ustring configFilename;
    rtengine::SafeKeyFile keyFile;
    bool hasKeyFile = true;

    imagesPaths.clear();

    // system theme will use the theme set in system.iconset
    if (opt.useSystemTheme) {
        configFilename = Glib::build_filename(argv0, Glib::build_filename("themes", "system.iconset"));
    }
    // Gtk theme will use the theme set in it's *.iconset fiel, if it exists
    else {
        configFilename = Glib::build_filename(argv0, Glib::build_filename("themes", Glib::ustring::format(opt.theme, ".iconset")));
    }

    try {
        if (!safe_file_test(configFilename, Glib::FILE_TEST_EXISTS) || !keyFile.load_from_file (configFilename)) {
            // ...otherwise fallback to the iconset set in default.iconset
            configFilename = Glib::build_filename(argv0, Glib::build_filename("themes", "Default.iconset"));

            if (!keyFile.load_from_file (configFilename)) {
                hasKeyFile = false;
            }
        }
    } catch (Glib::Error &err) {
        if (options.rtSettings.verbose) {
            printf("RTImage::setPaths / Error code %d while reading values from \"%s\":\n%s\n", err.code(), configFilename.c_str(), err.what().c_str());
        }
    } catch (...) {
        if (options.rtSettings.verbose) {
            printf("RTImage::setPaths / Unknown exception while trying to load \"%s\"!\n", configFilename.c_str());
        }
    }

    if (hasKeyFile && keyFile.has_group ("General")) {
        Glib::ustring iSet;

        if (keyFile.has_key ("General", "Iconset")) {
            iSet = keyFile.get_string ("General", "Iconset");
        }

        if (iSet.length()) {
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "actions"))));
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", iSet)));
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "devices"))));
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "places"))));
        }

        iSet.clear();

        if (keyFile.has_key ("General", "FallbackIconset")) {
            iSet = keyFile.get_string ("General", "FallbackIconset");
        }

        if (iSet.length()) {
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "actions"))));
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", iSet)));
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "devices"))));
            imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "places"))));
        }
    }

    // The images/ folder is the second fallback solution
    imagesPaths.push_back (Glib::build_filename(argv0, "images"));
}

