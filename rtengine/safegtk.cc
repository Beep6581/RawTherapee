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

FILE * safe_g_fopen(const Glib::ustring& src, const gchar *mode)
{
    return g_fopen(src.c_str(), mode);
}

bool safe_file_test (const Glib::ustring& filename, Glib::FileTest test)
{
    return Glib::file_test (filename, test);
}
