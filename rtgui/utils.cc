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
#include <guiutils.h>
#include <options.h>

void removeIfThere (Gtk::Container* cont, Gtk::Widget* w) {

    Glib::ListHandle<Gtk::Widget*> list = cont->get_children ();
    Glib::ListHandle<Gtk::Widget*>::iterator i = list.begin ();
    for (; i!=list.end() && *i!=w; i++);
    if (i!=list.end()) {
        w->reference ();
        cont->remove (*w);
    }
}

void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh) {

    if (options.thumbInterp==0)
        rtengine::nearestInterp (src, sw, sh, dst, dw, dh);
    else if (options.thumbInterp==1)
        rtengine::bilinearInterp (src, sw, sh, dst, dw, dh);
}
