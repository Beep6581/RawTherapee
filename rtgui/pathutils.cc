/*
 *  This file is part of RawTherapee.
 *
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

#include <glibmm/convert.h>
#include <glibmm/miscutils.h>

#include "pathutils.h"


Glib::ustring removeExtension (const Glib::ustring& filename)
{

    Glib::ustring bname = Glib::path_get_basename(filename);
    size_t lastdot = bname.find_last_of ('.');
    size_t lastwhitespace = bname.find_last_of (" \t\f\v\n\r");

    if (lastdot != bname.npos && (lastwhitespace == bname.npos || lastdot > lastwhitespace)) {
        return filename.substr (0, filename.size() - (bname.size() - lastdot));
    } else {
        return filename;
    }
}

Glib::ustring getExtension (const Glib::ustring& filename)
{

    Glib::ustring bname = Glib::path_get_basename(filename);
    size_t lastdot = bname.find_last_of ('.');
    size_t lastwhitespace = bname.find_last_of (" \t\f\v\n\r");

    if (lastdot != bname.npos && (lastwhitespace == bname.npos || lastdot > lastwhitespace)) {
        return filename.substr (filename.size() - (bname.size() - lastdot) + 1, filename.npos);
    } else {
        return "";
    }
}


// For an unknown reason, Glib::filename_to_utf8 doesn't work on reliably Windows,
// so we're using Glib::filename_to_utf8 for Linux/Apple and Glib::locale_to_utf8 for Windows.
Glib::ustring fname_to_utf8(const std::string &fname)
{
#ifdef _WIN32

    try {
        return Glib::locale_to_utf8(fname);
    } catch (Glib::Error&) {
        return Glib::convert_with_fallback(fname, "UTF-8", "ISO-8859-1", "?");
    }

#else

    return Glib::filename_to_utf8(fname);

#endif
}
