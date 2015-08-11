/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#ifndef _BROWSERFILTER_
#define _BROWSERFILTER_

#include "exiffiltersettings.h"
#include <glibmm.h>

class BrowserFilter
{

public:
    bool showRanked[6];
    bool showCLabeled[6];
    bool showTrash;
    bool showNotTrash;
    bool showEdited[2];
    bool showRecentlySaved[2];
    bool multiselect;

    Glib::ustring queryString;
    Glib::ustring queryFileName;

    bool exifFilterEnabled;
    ExifFilterSettings exifFilter;

    BrowserFilter ();
};

#endif
