/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-201 Gabor Horvath <hgabor@rawtherapee.com>
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

#include <memory>
#include <iostream>

#include <glib/gstdio.h>
#include <giomm.h>

#ifdef _WIN32
#include <windows.h>
#endif

namespace {
std::string getMD5 (const Glib::ustring& fname)
{

    auto file = Gio::File::create_for_path (fname);

    if (file && file->query_exists ())   {

#ifdef _WIN32

        std::unique_ptr<wchar_t, GFreeFunc> wfname (reinterpret_cast<wchar_t*> (g_utf8_to_utf16 (fname.c_str (), -1, NULL, NULL, NULL)), g_free);

        WIN32_FILE_ATTRIBUTE_DATA fileAttr;

        if (GetFileAttributesExW (wfname.get (), GetFileExInfoStandard, &fileAttr)) {
            // We use name, size and creation time to identify a file.
            const auto identifier = Glib::ustring::compose ("%1-%2-%3-%4", fileAttr.nFileSizeLow, fileAttr.ftCreationTime.dwHighDateTime, fileAttr.ftCreationTime.dwLowDateTime, fname);
            return Glib::Checksum::compute_checksum (Glib::Checksum::CHECKSUM_MD5, identifier);
        }

#else

        try {

            if (auto info = file->query_info ()) {
                // We only use name and size to identify a file.
                const auto identifier = Glib::ustring::compose ("%1%2", fname, info->get_size ());
                return Glib::Checksum::compute_checksum (Glib::Checksum::CHECKSUM_MD5, identifier);
            }

        } catch (Gio::Error&) {}

#endif

    }

    return {};
}

}
