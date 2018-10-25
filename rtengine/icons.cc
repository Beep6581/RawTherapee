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

#include "icons.h"

#include <iostream>

namespace rtengine
{

Glib::ustring imagePath;
Glib::ustring imagePath1;
Glib::ustring imagePath2;

Glib::ustring findIconAbsolutePath (const Glib::ustring& iconName, double &dpi)
{
    try {
        double fallBackDPI = 0.;
        Glib::ustring path;
        Glib::ustring fallBackPath;

        if (dpi == 96.) {
            path = imagePath1;
            fallBackPath = imagePath2;
            fallBackDPI = 192.;
        } else {
            path = imagePath2;
            fallBackPath = imagePath1;
            dpi = 192.;
            fallBackDPI = 96.;
        }

        auto iconPath = Glib::build_filename(path, iconName);

        if (Glib::file_test(iconPath, Glib::FILE_TEST_IS_REGULAR)) {
            return iconPath;
        } else {
            // fallback solution
            iconPath = Glib::build_filename(fallBackPath, iconName);
            if (Glib::file_test(iconPath, Glib::FILE_TEST_IS_REGULAR)) {
                dpi = fallBackDPI;
                return iconPath;
            } else {
                // second fallback solution: icon not resolution dependent
                iconPath = Glib::build_filename(imagePath, iconName);
                if (Glib::file_test(iconPath, Glib::FILE_TEST_IS_REGULAR)) {
                    dpi = -1;
                    return iconPath;
                }
            }
        }

    } catch(const Glib::Exception&) {}

    if (options.rtSettings.verbose) {
        std::cerr << "Icon \"" << iconName << "\" (" << dpi << " dpi) could not be found!" << std::endl;
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

    imagePath1 = Glib::build_filename(argv0, "images", "1", "dark");
    imagePath2 = Glib::build_filename(argv0, "images", "2", "dark");

    // The images folder is the second fallback solution.
    imagePath = Glib::build_filename(argv0, "images");
}

}
