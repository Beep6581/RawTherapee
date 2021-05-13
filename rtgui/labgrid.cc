/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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

// adapted from the "color correction" module of Darktable. Original copyright follows
/*
    copyright (c) 2009--2010 johannes hanika.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "labgrid.h"

#include "../rtengine/color.h"
#include "options.h"
#include "rtimage.h"

using rtengine::Color;


//-----------------------------------------------------------------------------
// LabGridArea
//-----------------------------------------------------------------------------

bool LabGridArea::notifyListener()
{
    if (listener) {
        const auto round =
            [](float v) -> float
            {
                return int(v * 1000) / 1000.f;
            };
        listener->panelChanged(evt, Glib::ustring::compose(evtMsg, round(high_a), round(high_b), round(low_a), round(low_b)));
    }
    return false;
}


LabGridArea::LabGridArea(rtengine::ProcEvent evt, const Glib::ustring &msg, bool enable_low):
    Gtk::DrawingArea(),
    evt(evt), evtMsg(msg),
    litPoint(NONE),
    low_a(0.f), high_a(0.f), low_b(0.f), high_b(0.f),
    defaultLow_a(0.f), defaultHigh_a(0.f), defaultLow_b(0.f), defaultHigh_b(0.f),
    listener(nullptr),
    edited(false),
    isDragged(false),
    low_enabled(enable_low)
{
    set_can_focus(false); // prevent moving the grid while you're moving a point
    add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::POINTER_MOTION_MASK);
    set_name("LabGrid");
    get_style_context()->add_class("drawingarea");
}

void LabGridArea::getParams(double &la, double &lb, double &ha, double &hb) const
{
    la = low_a;
    ha = high_a;
    lb = low_b;
    hb = high_b;
}


void LabGridArea::setParams(double la, double lb, double ha, double hb, bool notify)
{
    const double lo = -1.0;
    const double hi = 1.0;
    low_a = rtengine::LIM(la, lo, hi);
    low_b = rtengine::LIM(lb, lo, hi);
    high_a = rtengine::LIM(ha, lo, hi);
    high_b = rtengine::LIM(hb, lo, hi);
    queue_draw();
    if (notify) {
        notifyListener();
    }
}

void LabGridArea::setDefault (double la, double lb, double ha, double hb)
{
    defaultLow_a = la;
    defaultLow_b = lb;
    defaultHigh_a = ha;
    defaultHigh_b = hb;
}


void LabGridArea::reset(bool toInitial)
{
    if (toInitial) {
        setParams(defaultLow_a, defaultLow_b, defaultHigh_a, defaultHigh_b, true);
    } else {
        setParams(0., 0., 0., 0., true);
    }
}


void LabGridArea::setEdited(bool yes)
{
    edited = yes;
}


bool LabGridArea::getEdited() const
{
    return edited;
}


void LabGridArea::setListener(ToolPanelListener *l)
{
    listener = l;
}


void LabGridArea::on_style_updated ()
{
    setDirty(true);
    queue_draw ();
}


bool LabGridArea::on_draw(const ::Cairo::RefPtr<Cairo::Context> &crf)
{
    Gtk::Allocation allocation = get_allocation();
    allocation.set_x(0);
    allocation.set_y(0);

    // setDrawRectangle will allocate the backbuffer Surface
    if (setDrawRectangle(Cairo::FORMAT_ARGB32, allocation)) {
        setDirty(true);
    }

    if (!isDirty() || !surfaceCreated()) {
        return true;
    }

    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    Gtk::Border padding = getPadding(style);  // already scaled
    Cairo::RefPtr<Cairo::Context> cr = getContext();

    if (isDirty()) {
        int width = allocation.get_width();
        int height = allocation.get_height();

        int s = RTScalable::getScale();

        cr->set_line_cap(Cairo::LINE_CAP_SQUARE);

        // clear background
        cr->set_source_rgba (0., 0., 0., 0.);
        cr->set_operator (Cairo::OPERATOR_CLEAR);
        cr->paint ();
        cr->set_operator (Cairo::OPERATOR_OVER);
        style->render_background(cr,
                inset * s + padding.get_left() - s,
                inset * s + padding.get_top() - s,
                width - 2 * inset * s - padding.get_right() - padding.get_left() + 2 * s,
                height - 2 * inset * s - padding.get_top() - padding.get_bottom() + 2 * s
                );

        // drawing the cells
        cr->translate(inset * s + padding.get_left(), inset * s + padding.get_top());
        cr->set_antialias(Cairo::ANTIALIAS_NONE);
        width -= 2 * inset * s + padding.get_right() + padding.get_left();
        height -= 2 * inset * s + padding.get_top() + padding.get_bottom();

        // flip y:
        cr->translate(0, height);
        cr->scale(1., -1.);
        const int cells = 8;
        float step = 12000.f / float(cells/2);
        double cellW = double(width) / double(cells);
        double cellH = double(height) / double(cells);
        double cellYMin = 0.;
        double cellYMax = std::floor(cellH);
        for (int j = 0; j < cells; j++) {
            double cellXMin = 0.;
            double cellXMax = std::floor(cellW);
            for (int i = 0; i < cells; i++) {
                float R, G, B;
                float x, y, z;
                int ii = i - cells/2;
                int jj = j - cells/2;
                float a = step * (ii + 0.5f);
                float b = step * (jj + 0.5f);
                Color::Lab2XYZ(25000.f, a, b, x, y, z);
                Color::xyz2srgb(x, y, z, R, G, B);
                cr->set_source_rgb(R / 65535.f, G / 65535.f, B / 65535.f);
                cr->rectangle(
                        cellXMin,
                        cellYMin,
                        cellXMax - cellXMin - (i == cells-1 ? 0. : double(s)),
                        cellYMax - cellYMin - (j == cells-1 ? 0. : double(s))
                        );
                cellXMin = cellXMax;
                cellXMax = std::floor(cellW * double(i+2) + 0.01);
                cr->fill();
            }
            cellYMin = cellYMax;
            cellYMax = std::floor(cellH * double(j+2) + 0.01);
        }

        // drawing the connection line
        cr->set_antialias(Cairo::ANTIALIAS_DEFAULT);
        float loa, hia, lob, hib;
        loa = .5 * (width + width * low_a);
        hia = .5 * (width + width * high_a);
        lob = .5 * (height + height * low_b);
        hib = .5 * (height + height * high_b);
        cr->set_line_width(2. * double(s));
        cr->set_source_rgb(0.6, 0.6, 0.6);
        cr->move_to(loa, lob);
        cr->line_to(hia, hib);
        cr->stroke();

        // drawing points
        if (low_enabled) {
            cr->set_source_rgb(0.1, 0.1, 0.1);
            if (litPoint == LOW) {
                cr->arc(loa, lob, 5 * s, 0, 2. * rtengine::RT_PI);
            } else {
                cr->arc(loa, lob, 3 * s, 0, 2. * rtengine::RT_PI);
            }
            cr->fill();
        }

        cr->set_source_rgb(0.9, 0.9, 0.9);
        if (litPoint == HIGH) {
            cr->arc(hia, hib, 5 * s, 0, 2. * rtengine::RT_PI);
        } else {
            cr->arc(hia, hib, 3 * s, 0, 2. * rtengine::RT_PI);
        }
        cr->fill();
    }

    copySurface(crf);
    return false;
}


bool LabGridArea::on_button_press_event(GdkEventButton *event)
{
    if (event->button == 1) {
        if (event->type == GDK_2BUTTON_PRESS) {
            switch (litPoint) {
            case NONE:
                low_a = low_b = high_a = high_b = 0.f;
                break;
            case LOW:
                low_a = low_b = 0.f;
                break;
            case HIGH:
                high_a = high_b = 0.f;
                break;
            }
            edited = true;
            notifyListener();
            queue_draw();
        } else if (event->type == GDK_BUTTON_PRESS && litPoint != NONE) {
            isDragged = true;
        }
        return false;
    }
    return true;
}


bool LabGridArea::on_button_release_event(GdkEventButton *event)
{
    if (event->button == 1) {
        isDragged = false;
        return false;
    }
    return true;
}


bool LabGridArea::on_motion_notify_event(GdkEventMotion *event)
{
    if (isDragged && delayconn.connected()) {
        delayconn.disconnect();
    }

    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    Gtk::Border padding = getPadding(style);  // already scaled

    State oldLitPoint = litPoint;

    int s = RTScalable::getScale();
    int width = get_allocated_width() - 2 * inset * s - padding.get_right() - padding.get_left();
    int height = get_allocated_height() - 2 * inset * s - padding.get_top() - padding.get_bottom();
    const float mouse_x = std::min(double(std::max(event->x - inset * s - padding.get_right(), 0.)), double(width));
    const float mouse_y = std::min(double(std::max(get_allocated_height() - 1 - event->y - inset * s - padding.get_bottom(), 0.)), double(height));
    const float ma = (2.f * mouse_x - width) / width;
    const float mb = (2.f * mouse_y - height) / height;
    if (isDragged) {
        if (litPoint == LOW) {
            low_a = ma;
            low_b = mb;
        } else if (litPoint == HIGH) {
            high_a = ma;
            high_b = mb;
        }
        edited = true;
        grab_focus();
        if (options.adjusterMinDelay == 0) {
            notifyListener();
        } else {
            delayconn = Glib::signal_timeout().connect(sigc::mem_fun(*this, &LabGridArea::notifyListener), options.adjusterMinDelay);
        }
        queue_draw();
    } else {
        litPoint = NONE;
        float la = low_a;
        float lb = low_b;
        float ha = high_a;
        float hb = high_b;
        const float thrs = 0.05f;
        const float distlo = (la - ma) * (la - ma) + (lb - mb) * (lb - mb);
        const float disthi = (ha - ma) * (ha - ma) + (hb - mb) * (hb - mb);
        if (low_enabled && distlo < thrs * thrs && distlo < disthi) {
            litPoint = LOW;
        } else if (disthi < thrs * thrs && disthi <= distlo) {
            litPoint = HIGH;
        }
        if ((oldLitPoint == NONE && litPoint != NONE) || (oldLitPoint != NONE && litPoint == NONE)) {
            queue_draw();
        }
    }
    return true;
}


Gtk::SizeRequestMode LabGridArea::get_request_mode_vfunc() const
{
    return Gtk::SIZE_REQUEST_HEIGHT_FOR_WIDTH;
}


void LabGridArea::get_preferred_width_vfunc(int &minimum_width, int &natural_width) const
{
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    Gtk::Border padding = getPadding(style);  // already scaled
    int s = RTScalable::getScale();
    int p = padding.get_left() + padding.get_right();

    minimum_width = 50 * s + p;
    natural_width = 150 * s + p;  // same as GRAPH_SIZE from mycurve.h
}


void LabGridArea::get_preferred_height_for_width_vfunc(int width, int &minimum_height, int &natural_height) const
{
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    Gtk::Border padding = getPadding(style);  // already scaled

    minimum_height = natural_height = width - padding.get_left() - padding.get_right() + padding.get_top() + padding.get_bottom();
}


bool LabGridArea::lowEnabled() const
{
    return low_enabled;
}


void LabGridArea::setLowEnabled(bool yes)
{
    if (low_enabled != yes) {
        low_enabled = yes;
        queue_draw();
    }
}


//-----------------------------------------------------------------------------
// LabGrid
//-----------------------------------------------------------------------------

LabGrid::LabGrid(rtengine::ProcEvent evt, const Glib::ustring &msg, bool enable_low):
    grid(evt, msg, enable_low)
{
    Gtk::Button *reset = Gtk::manage(new Gtk::Button());
    reset->set_tooltip_markup(M("ADJUSTER_RESET_TO_DEFAULT"));
    reset->add(*Gtk::manage(new RTImage("undo-small.png", "redo-small.png")));
    reset->signal_button_release_event().connect(sigc::mem_fun(*this, &LabGrid::resetPressed));

    setExpandAlignProperties(reset, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
    reset->set_relief(Gtk::RELIEF_NONE);
    reset->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
    reset->set_can_focus(false);
    reset->set_size_request(-1, 20);

    pack_start(grid, true, true);
    pack_start(*reset, false, false);
    show_all_children();
}


bool LabGrid::resetPressed(GdkEventButton *event)
{
    grid.reset(event->state & GDK_CONTROL_MASK);
    return false;    
}
