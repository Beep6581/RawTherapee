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

#include "shcselector.h"
#include "multilangmgr.h"
#include "mycurve.h"

SHCSelector::SHCSelector() : ColoredBar(RTO_Left2Right), movingPosition(-1), tmpX(0.0), tmpPos(0.0), wslider(0.0), cl(nullptr)
{

    positions[0] = defaults[0] = 0.25;
    positions[1] = defaults[1] = 0.5;
    positions[2] = defaults[2] = 0.75;
    leftMargin = RADIUS;
    rightMargin = RADIUS;

    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->add_class(GTK_STYLE_CLASS_DEFAULT);
    style->add_class(GTK_STYLE_CLASS_SCALE);
    style->add_class(GTK_STYLE_CLASS_SLIDER);

    // TODO: This is a hack :) ; change this name to a specific one and create a new entry in all gtkrc theme files
    set_name("ThresholdSelector");
    set_can_focus(false);
    set_tooltip_text(M("SHCSELECTOR_TOOLTIP"));
}

Gtk::SizeRequestMode SHCSelector::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_HEIGHT_FOR_WIDTH;
}

void SHCSelector::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    int minimumWidth = 0;
    int naturalWidth = 0;
    get_preferred_width_vfunc (minimumWidth, naturalWidth);
    get_preferred_height_for_width_vfunc (minimumWidth, minimum_height, natural_height);
}

void SHCSelector::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = 100;
    natural_width = 150;
}

void SHCSelector::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    natural_height = minimum_height = 14;
}

void SHCSelector::get_preferred_width_for_height_vfunc (int width, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}

void SHCSelector::setMargins(int left, int right)
{
    leftMargin = left;
    rightMargin = right;
}

void SHCSelector::setDefaults (double spos, double cpos, double hpos)
{
    defaults[0] = spos;
    defaults[1] = cpos;
    defaults[2] = hpos;
}

void SHCSelector::setPositions (double spos, double cpos, double hpos)
{

    positions[0] = spos;
    positions[1] = cpos;
    positions[2] = hpos;

    queue_draw ();
}

void SHCSelector::getPositions (double& spos, double& cpos, double& hpos)
{

    spos = positions[0];
    cpos = positions[1];
    hpos = positions[2];
}

void SHCSelector::on_realize()
{

    Gtk::DrawingArea::on_realize();

    add_events(Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
}

bool SHCSelector::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    Gdk::RGBA c;

    int w = get_width () - leftMargin - rightMargin;
    int h = get_height ();

    wslider = std::max(int(h / 5), 10);
    double hwslider = double(wslider) / 2.;

    Gtk::StateFlags state = !is_sensitive() ? Gtk::STATE_FLAG_INSENSITIVE : Gtk::STATE_FLAG_NORMAL;
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();




    //style->render_background(cr, leftMargin, 0, w, h);
    //return true;




    // clear bg

    // set the box's colors
    cr->set_line_width (1.0);
    cr->set_line_cap(Cairo::LINE_CAP_BUTT);

    if (is_sensitive() && canGetColors()) {
        // gradient background
        Glib::RefPtr<Gdk::Window> win = get_window();
        // this will eventually create/update the off-screen pixmap
        setDrawRectangle(win, leftMargin + 1, 1, w - 2, int(float(h) * 5.5f / 7.f + 0.5f));
        // that we're displaying here
        ColoredBar::expose(cr);
    }

    /*  useless
    else {
        // solid background
        // draw the box's background
        style->render_background(cr, leftMargin+1, 1, w-2, int(float(h)*5.5f/7.f+0.5f));
    }
    */

    // draw the box's borders
    cr->set_line_width (1.);
    c = style->get_border_color(state);
    cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
    cr->rectangle (leftMargin + 0.5, 0.5, w - 1, int(float(h) * 5.5f / 7.f + 0.5f) + 1);
    cr->stroke ();

    // draw sliders
    //cr->set_line_width (1.);
    for (int i = 0; i < 3; i++) {
        if (i == movingPosition) {
            style->set_state(Gtk::STATE_FLAG_ACTIVE);
        }
        /*
        else if (i==litCursor)
            style->set_state(Gtk::STATE_FLAG_PRELIGHT);
        */
        else if (!is_sensitive()) {
            style->set_state(Gtk::STATE_FLAG_INSENSITIVE);
        } else {
            style->set_state(Gtk::STATE_FLAG_NORMAL);
        }

        style->render_slider(cr, leftMargin + 0.5 + (w - 1)*positions[i] - hwslider, vb, wslider, h - vb, Gtk::ORIENTATION_HORIZONTAL);
        style->set_state(Gtk::STATE_FLAG_NORMAL);
    }

    /*
    for (int i=0; i<3; i++) {
        double posX = leftMargin+0.5+(w-1)*positions[i];
        double arrowY = h-(h*3.5/7.-0.5)-vb;
        double baseY = h-0.5-vb;
        double centerY = (arrowY+baseY)/2.;
        cr->move_to (posX, arrowY);
        cr->line_to (posX+hwslider, centerY);
        cr->line_to (posX+hwslider, baseY);
        cr->line_to (posX-hwslider, baseY);
        cr->line_to (posX-hwslider, centerY);
        cr->close_path();
        if (i==movingPosition) {
            // moved (selected)
            c = style->get_background_color(Gtk::STATE_FLAG_SELECTED);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->fill_preserve ();
            c = style->get_border_color (Gtk::STATE_FLAG_SELECTED);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->stroke ();
        }
        */

    /*
    else if (i==litCursor) {
        // prelight
        c = style->get_background_color(Gtk::STATE_FLAG_PRELIGHT);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->fill_preserve ();
        c = style->get_border_color (Gtk::STATE_FLAG_PRELIGHT);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->stroke ();
    }
    */

    /*
    else {
        // normal
        c = style->get_background_color(is_sensitive() ? Gtk::STATE_FLAG_ACTIVE : Gtk::STATE_FLAG_INSENSITIVE);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->fill_preserve ();
        c = style->get_border_color (is_sensitive() ? Gtk::STATE_FLAG_ACTIVE : Gtk::STATE_FLAG_INSENSITIVE);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->stroke ();
    }
    }
    */

    // draw text for the slider that is being moved
    cr->set_line_width (0.5);

    if (movingPosition >= 0) {
        int i = movingPosition;
        int offset;
        int layout_width, layout_height;
        Glib::RefPtr<Pango::Layout> layout = create_pango_layout(Glib::ustring::format(std::setprecision(2), positions[i]));
        layout->get_pixel_size(layout_width, layout_height);
        offset = positions[i] > 0.5 ? -layout_width - 1 - hwslider : 1 + hwslider;
        c = style->get_background_color(state);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());

        cr->set_line_width(3.);
        cr->set_line_join(Cairo::LINE_JOIN_ROUND);
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);

        cr->move_to (leftMargin + w * positions[i] + offset, 0.);
        layout->add_to_cairo_context (cr);
        cr->stroke_preserve();
        c = style->get_color(Gtk::STATE_FLAG_PRELIGHT);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->fill ();
    }

    return true;
}

bool SHCSelector::on_button_press_event (GdkEventButton* event)
{

    // check if a slider is under the cursor
    double w = double(get_width () - leftMargin - rightMargin);
    movingPosition = -1;

    for (int i = 0; i < 3; i++)
        if (event->x > double(leftMargin) + w * positions[i] - wslider / 2. && event->x < double(leftMargin) + w * positions[i] + wslider / 2) {
            movingPosition = i;
            tmpX = event->x;
            tmpPos = positions[i];
            break;
        }

    queue_draw ();
    return true;
}

bool SHCSelector::on_button_release_event (GdkEventButton* event)
{

    if (event->button == 1) {
        if (movingPosition >= 0) {
            movingPosition = -1;
            queue_draw ();
        }
    } else if (event->button == 3) {
        if (movingPosition >= 0) {
            movingPosition = -1;
        }

        // right mouse button reset the selector to the stored default values
        if (reset()) {
            // rest has modified the values
            if (cl) {
                cl->shcChanged ();
            }
        }
    }

    return true;
}

bool SHCSelector::on_motion_notify_event (GdkEventMotion* event)
{

    if (movingPosition >= 0) {
        int w = get_width ();
        positions[movingPosition] = tmpPos + (event->x - tmpX) / w;

        if (positions[movingPosition] < 0) {
            positions[movingPosition] = 0.0;
        }

        if (movingPosition > 0 && positions[movingPosition] < positions[movingPosition - 1] + wslider / w) {
            positions[movingPosition] = positions[movingPosition - 1] + wslider / w;
        }

        if (positions[movingPosition] > 1.0) {
            positions[movingPosition] = 1.0;
        }

        if (movingPosition < 2 && positions[movingPosition] > positions[movingPosition + 1] - wslider / w) {
            positions[movingPosition] = positions[movingPosition + 1] - wslider / w;
        }

        if (cl) {
            cl->shcChanged ();
        }

        queue_draw ();
    }

    return true;
}

void SHCSelector::styleChanged (const Glib::RefPtr<Gtk::Style>& style)
{

    queue_draw ();
}

bool SHCSelector::reset ()      //  : movingPosition(-1), cl(NULL) {
{
    if ( positions[0] != defaults[0] ||
            positions[1] != defaults[1] ||
            positions[2] != defaults[2]
       ) {

        positions[0] = defaults[0];
        positions[1] = defaults[1];
        positions[2] = defaults[2];
        queue_draw ();
        return true;
    }

    return false;
}

void SHCSelector::refresh()
{
    setDirty(true);
    Glib::RefPtr<Gdk::Window> win = get_window();

    if (win) {
        win->invalidate(true);
    }
}
