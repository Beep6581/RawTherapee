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
#ifndef __UTILS_
#define __UTILS_

#include <gtkmm.h>
#include <rtengine.h>

bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference=true);
void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh);
Glib::ustring removeExtension (const Glib::ustring& filename);
Glib::ustring getExtension (const Glib::ustring& filename);
void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams);

#endif
