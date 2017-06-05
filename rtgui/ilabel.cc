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
#include "ilabel.h"

ILabel::ILabel (const Glib::ustring &lab) : label(lab) {}

void ILabel::on_realize()
{

    Gtk::DrawingArea::on_realize();
    add_events(Gdk::EXPOSURE_MASK);

    Glib::RefPtr<Pango::Layout> fn = create_pango_layout(label);
    fn->set_markup (label);
    int labw, labh;
    fn->get_pixel_size (labw, labh);
    set_size_request (2 + labw, 2 + labh);
}

bool ILabel::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    Glib::RefPtr<Gtk::StyleContext> style = get_style_context ();

    Gtk::StateFlags state = get_state_flags();
    Gdk::RGBA c = style->get_background_color(state);
    cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
    cr->rectangle (0, 0, get_width (), get_height());
    cr->fill ();

    Glib::RefPtr<Pango::Layout> fn = create_pango_layout (label);
    fn->set_markup (label);
    cr->move_to(1., 1.);
    fn->add_to_cairo_context(cr);
    c = style->get_color (state);
    cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
    cr->fill();

    return true;
}

void ILabel::on_style_updated ()
{

    Glib::RefPtr<Pango::Layout> fn = create_pango_layout(label);
    fn->set_markup (label);
    int labw, labh;
    fn->get_pixel_size (labw, labh);
    set_size_request (2 + labw, 2 + labh);
}
