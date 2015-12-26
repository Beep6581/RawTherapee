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

// Opens a file for binary writing and request exclusive lock (cases were you need "wb" mode plus locking)
// (Important on Windows to prevent Explorer to crash RT when parallel scanning e.g. a currently written image file)
FILE * safe_g_fopen_WriteBinLock(const Glib::ustring& fname)
{
    FILE* f = NULL;

#ifdef WIN32
    // g_fopen just uses _wfopen internally on Windows, does not lock access and has no options to set this
    // so use a native function to work around this problem
    wchar_t *wFname = (wchar_t*)g_utf8_to_utf16 (fname.c_str(), -1, NULL, NULL, NULL);
    HANDLE hFile = CreateFileW(wFname, GENERIC_READ | GENERIC_WRITE, 0 /* no sharing allowed */, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    g_free(wFname);

    if (hFile == INVALID_HANDLE_VALUE) {
        f = NULL;
    } else {
        f = _fdopen( _open_osfhandle((intptr_t)hFile, 0) , "wb");
    }

#else
    f = safe_g_fopen(fname, "wb");
#endif

    return f;
}

// Covers old UNIX ::open, which expects ANSI instead of UTF8 on Windows
int safe_open_ReadOnly(const char *fname)
{
    int fd = -1;

#ifdef WIN32
    // First convert UTF8 to UTF16, then use Windows function to open
    wchar_t *wFname = (wchar_t*)g_utf8_to_utf16 (fname, -1, NULL, NULL, NULL);
    HANDLE hFile = CreateFileW(wFname, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    g_free(wFname);

    // convert back to old file descriptor format
    if (hFile != INVALID_HANDLE_VALUE) {
        fd = _open_osfhandle((intptr_t)hFile, 0);
    }

#else
    fd = ::open(fname, O_RDONLY);
#endif

    return fd;
}


FILE * safe_g_fopen(const Glib::ustring& src, const gchar *mode)
{
    return g_fopen(src.c_str(), mode);
}

bool safe_file_test (const Glib::ustring& filename, Glib::FileTest test)
{
    return Glib::file_test (filename, test);
}
