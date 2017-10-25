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

#include "icons.h"

#include <iostream>
#pragma GCC diagnostic warning "-Wextra"

namespace rtengine
{

std::vector<Glib::ustring> imagePaths;

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

Glib::ustring findIconAbsolutePath (const Glib::ustring& iconName)
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

void setPaths ()
{
    // TODO: Forcing the Dark theme, so reading the icon set files is useless for now...

    /*Glib::ustring iconSet;

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
    }*/

    imagePaths.clear ();

    imagePaths.push_back (Glib::build_filename(argv0, "images", "Dark"));
    imagePaths.push_back (Glib::build_filename(argv0, "images", "Dark", "actions"));
    imagePaths.push_back (Glib::build_filename(argv0, "images", "Dark", "devices"));
    imagePaths.push_back (Glib::build_filename(argv0, "images", "Dark", "places"));

    // The images folder is the second fallback solution.
    imagePaths.push_back (Glib::build_filename(argv0, "images"));
}

}
