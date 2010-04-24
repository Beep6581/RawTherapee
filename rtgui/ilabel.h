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
#ifndef _ILABEL_
#define _ILABEL_

#include <gtkmm.h>

class ILabel : public Gtk::DrawingArea {

    Glib::ustring label;

    public:
        ILabel (Glib::ustring lab);
        bool on_expose_event(GdkEventExpose* event);
        void on_realize();
        void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
};

#endif

