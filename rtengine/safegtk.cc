/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2010 Sasha Vasko <sasha@aftercode.net>
 *  Copyright (c) 2010 Oliver Duis <www.oliverduis.de>
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

#include "safegtk.h"

#include <fcntl.h>
#include <cstdio>
#include <memory>

#include <glib/gstdio.h>

#ifdef WIN32
#include <windows.h>
#endif

/*
 * For an unknown reason, Glib::filename_to_utf8 doesn't work on Windows, so we're using
 * Glib::filename_to_utf8 for Linux/Apple and Glib::locale_to_utf8 for Windows
 */
Glib::ustring safe_filename_to_utf8 (const std::string& src)
{
    Glib::ustring utf8_str;

#ifdef WIN32

    try {
        utf8_str = Glib::locale_to_utf8(src);
    } catch (const Glib::Error& e) {
        utf8_str = Glib::convert_with_fallback(src, "UTF-8", "ISO-8859-1", "?");
    }

#else

    utf8_str = Glib::filename_to_utf8(src);

#endif

    return utf8_str;
}

Glib::ustring safe_locale_to_utf8 (const std::string& src)
{
    Glib::ustring utf8_str;

    try {
        utf8_str = Glib::locale_to_utf8(src);
    } catch (const Glib::Error& e) {
        utf8_str = Glib::convert_with_fallback(src, "UTF-8", "ISO-8859-1", "?");
    }

    return utf8_str;
}

std::string safe_locale_from_utf8 (const Glib::ustring& utf8_str)
{
    std::string str;

    try {
        str = Glib::locale_from_utf8(utf8_str);
    } catch (Glib::Error&) {}

    return str;
}

FILE * safe_g_fopen(const Glib::ustring& src, const gchar *mode)
{
    return g_fopen(src.c_str(), mode);
}

bool safe_file_test (const Glib::ustring& filename, Glib::FileTest test)
{
    return Glib::file_test (filename, test);
}
