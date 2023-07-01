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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iomanip>

#include "shcselector.h"
#include "multilangmgr.h"
#include "mycurve.h"
#include "rtscalable.h"

SHCSelector::SHCSelector() : movingPosition(-1), tmpX(0.0), tmpPos(0.0), wslider(0.0), cl(nullptr), coloredBar(RTO_Left2Right)
{

    int s = RTScalable::getScale();

    positions[0] = defaults[0] = 0.25;
    positions[1] = defaults[1] = 0.5;
    positions[2] = defaults[2] = 0.75;
    leftMargin = (RADIUS - 1.5) * s;
    rightMargin = (RADIUS - 1.5) * s;

    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->add_class("drawingarea");
    style->add_class(GTK_STYLE_CLASS_TROUGH);
    //style->add_class(GTK_STYLE_CLASS_SCALE);
    //style->add_class(GTK_STYLE_CLASS_SLIDER);

    // TODO: This is a hack :) ; change this name to a specific one and create a new entry in all CSS theme files
    set_name("ThresholdSelector");
    set_can_focus(false);
    set_tooltip_text(M("SHCSELECTOR_TOOLTIP"));
}

Gtk::SizeRequestMode SHCSelector::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_CONSTANT_SIZE;
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
    int s = RTScalable::getScale();
    minimum_width = 100 * s;
    natural_width = 150 * s;
}

void SHCSelector::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    natural_height = minimum_height = 14 * RTScalable::getScale();
}

void SHCSelector::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
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

void SHCSelector::updateBackBuffer()
{

    if (!get_realized() || !isDirty() || !get_width() || !get_height())  {
        return;
    }

    // This will create or update the size of the BackBuffer::surface
    setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, get_width(), get_height(), true);

    if (!surface) {
        return;
    }

    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

    cr->set_source_rgba (0., 0., 0., 0.);
    cr->set_operator (Cairo::OPERATOR_CLEAR);
    cr->paint ();
    cr->set_operator (Cairo::OPERATOR_OVER);

    int w = get_width () - leftMargin - rightMargin;
    int h = get_height ();

    double s = RTScalable::getScale();

    wslider = (double)std::max(h / 5, 10) * s;
    double hwslider = wslider / 2.;

    // clear bg
    cr->set_source_rgba (0., 0., 0., 0.);
    cr->set_operator (Cairo::OPERATOR_CLEAR);
    cr->paint ();
    cr->set_operator (Cairo::OPERATOR_OVER);


    // set the box's colors
    cr->set_line_width (1.0 * s);
    cr->set_antialias(Cairo::ANTIALIAS_SUBPIXEL);
    cr->set_line_cap(Cairo::LINE_CAP_BUTT);

    int coloredBarHeight = (int)((double)h * 5.5 / 7. + 0.5);
    if (is_sensitive() && coloredBar.canGetColors()) {
        // gradient background

        // this will eventually create/update the off-screen BackBuffer
        coloredBar.setDrawRectangle(leftMargin + 1 * (int)s, 1 * (int)s, w - 2 * (int)s, coloredBarHeight - 2 * (int)s);
        // that we're displaying here
        coloredBar.expose(*this, cr);
    } else {
        style->render_background(cr, leftMargin + 1 * (int)s, 1 * (int)s, w - 2 * (int)s, coloredBarHeight - 2 * (int)s);
    }
    // draw the box's borders
    style->render_frame(cr, leftMargin, 0, w, coloredBarHeight);


    // draw sliders
    for (int i = 0; i < 3; i++) {
        if (i == movingPosition) {
            style->set_state(Gtk::STATE_FLAG_ACTIVE);
        }
        else if (!is_sensitive()) {
            style->set_state(Gtk::STATE_FLAG_INSENSITIVE);
        } else {
            style->set_state(Gtk::STATE_FLAG_NORMAL);
        }

        style->render_slider(cr, (double)leftMargin + 1. * s + ((double)w - 2. * s) * positions[i] - (double)hwslider, (double)vb * s, wslider, (double)h - (double)vb * s, Gtk::ORIENTATION_VERTICAL);
        style->set_state(Gtk::STATE_FLAG_NORMAL);
    }

    // draw text for the slider that is being moved
    if (movingPosition >= 0) {
        int i = movingPosition;
        int offset;
        int layout_width = 0, layout_height = 0;

        Glib::RefPtr<Pango::Context> context = get_pango_context () ;
        Pango::FontDescription fontd(get_style_context()->get_font());

        // update font
        fontd.set_weight (Pango::WEIGHT_NORMAL);
        fontd.set_absolute_size((double)h * 0.8 * (double)Pango::SCALE);
        context->set_font_description (fontd);

        Glib::RefPtr<Pango::Layout> layout = create_pango_layout(Glib::ustring::format(std::setprecision(2), positions[i]));
        layout->get_pixel_size(layout_width, layout_height);
        offset = positions[i] > 0.5 ? -layout_width - 1 * (int)s - hwslider : 1 * (int)s + hwslider;
        cr->set_source_rgb (0., 0., 0.);

        cr->set_line_width(3. * s);
        cr->set_line_join(Cairo::LINE_JOIN_ROUND);
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);

        cr->move_to ((double)leftMargin + (double)w * positions[i] + (double)offset, 0.);
        layout->add_to_cairo_context (cr);
        cr->stroke_preserve();
        cr->set_line_width(0.5 * s);
        cr->set_source_rgb (1., 1., 1.);
        cr->fill ();
    }
}

bool SHCSelector::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    // on_realize & updateBackBuffer have to be called before
    if (get_realized() && get_width() && get_height()) {
        if (isDirty()) {
            updateBackBuffer();
        }

        if (surface) {
            copySurface(cr);
        }
    }

    return true;
}

bool SHCSelector::on_button_press_event (GdkEventButton* event)
{

    // check if a slider is under the cursor
    double w = double(get_width () - leftMargin - rightMargin);
    movingPosition = -1;

    double s = RTScalable::getScale();

    for (int i = 0; i < 3; i++) {
        double currPos = double(leftMargin) + 1. * s + (w - 2. * s) * positions[i];
        double hwslider = wslider / 2.;
        if (event->x >= currPos - hwslider && event->x <= currPos + hwslider) {
            movingPosition = i;
            tmpX = event->x;
            tmpPos = positions[i];
            break;
        }
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
        double s = RTScalable::getScale();
        double innerw = double(get_width () - leftMargin - rightMargin) - 2. * s;
        positions[movingPosition] = tmpPos + (event->x - tmpX) / innerw;

        if (positions[movingPosition] < 0) {
            positions[movingPosition] = 0.0;
        }

        if (movingPosition > 0 && positions[movingPosition] < positions[movingPosition - 1] + wslider / innerw) {
            positions[movingPosition] = positions[movingPosition - 1] + wslider / innerw;
        }

        if (positions[movingPosition] > 1.0) {
            positions[movingPosition] = 1.0;
        }

        if (movingPosition < 2 && positions[movingPosition] > positions[movingPosition + 1] - wslider / innerw) {
            positions[movingPosition] = positions[movingPosition + 1] - wslider / innerw;
        }

        if (cl) {
            cl->shcChanged ();
        }

        setDirty(true);
        queue_draw ();
    }

    return true;
}

void SHCSelector::styleChanged (const Glib::RefPtr<Gtk::StyleContext>& style)
{

    setDirty(true);
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
        setDirty(true);
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
