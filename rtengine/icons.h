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
#ifndef _ICONS_
#define _ICONS_

#include <glibmm.h>
#include "../rtgui/options.h"

namespace rtengine
{

/**
 * @brief Find the absolute path for an icon name of the desired DPI.
 *
 * @return the absolute path to the icon file
 * @return deiredDPI is updated to another DPI if the desired one wasn't found but a fallback solution exist, or -1 if the icon is not resolution dependent
 */
Glib::ustring findIconAbsolutePath (const Glib::ustring& iconName, double &dpi);
void setPaths ();

}

#endif
