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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
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
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "labgrid.h"

using rtengine::Color;


bool LabGrid::notifyListener()
{
    if (listener) {
        listener->panelChanged(evt, Glib::ustring::compose(M("TP_COLORTONING_LABGRID_VALUES"), int(low_a), int(low_b), int(high_a), int(high_b)));
    }
    return false;
}


LabGrid::LabGrid(rtengine::ProcEvent evt):
    Gtk::DrawingArea(),
    evt(evt), litPoint(NONE),
    low_a(0.f), high_a(0.f), low_b(0.f), high_b(0.f),
    defaultLow_a(0.f), defaultHigh_a(0.f), defaultLow_b(0.f), defaultHigh_b(0.f),
    listener(nullptr),
    edited(false),
    isDragged(false)
{
    set_can_focus(true);
    add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::POINTER_MOTION_MASK);
}

void LabGrid::getParams(double &la, double &lb, double &ha, double &hb) const
{
    la = low_a;
    ha = high_a;
    lb = low_b;
    hb = high_b;
}


void LabGrid::setParams(double la, double lb, double ha, double hb, bool notify)
{
    low_a = la;
    low_b = lb;
    high_a = ha;
    high_b = hb;
    queue_draw();
    if (notify) {
        notifyListener();
    }
}

void LabGrid::setDefault (double la, double lb, double ha, double hb)
{
    defaultLow_a = la;
    defaultLow_b = lb;
    defaultHigh_a = ha;
    defaultHigh_b = hb;
}


void LabGrid::reset(bool toInitial)
{
    if (toInitial) {
        setParams(defaultLow_a, defaultLow_b, defaultHigh_a, defaultHigh_b, true);
    } else {
        setParams(0., 0., 0., 0., true);
    }
}


void LabGrid::setEdited(bool yes)
{
    edited = yes;
}


bool LabGrid::getEdited() const
{
    return edited;
}


void LabGrid::setListener(ToolPanelListener *l)
{
    listener = l;
}


void LabGrid::on_style_updated ()
{
    setDirty(true);
    queue_draw ();
}


bool LabGrid::on_draw(const ::Cairo::RefPtr<Cairo::Context> &crf)
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
    Cairo::RefPtr<Cairo::Context> cr = getContext();

    if (isDirty()) {
        int width = allocation.get_width();
        int height = allocation.get_height();

        cr->set_line_cap(Cairo::LINE_CAP_SQUARE);

        // clear background
        cr->set_source_rgba (0., 0., 0., 0.);
        cr->set_operator (Cairo::OPERATOR_CLEAR);
        cr->paint ();
        cr->set_operator (Cairo::OPERATOR_OVER);
        style->render_background(cr, 0, 0, width, height);

        // drawing the cells
        cr->translate(inset, inset);
        cr->set_antialias(Cairo::ANTIALIAS_NONE);
        width -= 2 * inset;
        height -= 2 * inset;
        // flip y:
        cr->translate(0, height);
        cr->scale(1., -1.);
        const int cells = 8;
        float step = rtengine::ColorToningParams::LABGRID_CORR_MAX / float(cells/2);
        for (int j = 0; j < cells; j++) {
            for (int i = 0; i < cells; i++) {
                float R, G, B;
                float x, y, z;
                int ii = i - cells/2;
                int jj = j - cells/2;
                float a = step * (ii + 0.5);
                float b = step * (jj + 0.5);
                Color::Lab2XYZ(25000.f, a, b, x, y, z);
                Color::xyz2srgb(x, y, z, R, G, B);
                cr->set_source_rgb(R / 65535.f, G / 65535.f, B / 65535.f);
                cr->rectangle(width * i / float(cells), height * j / float(cells), width / float(cells) - 1, height / float(cells) - 1);
                cr->fill();
            }
        }

        // drawing the connection line
        cr->set_antialias(Cairo::ANTIALIAS_DEFAULT);
        float loa, hia, lob, hib;
        loa = .5f * (width + width * low_a / rtengine::ColorToningParams::LABGRID_CORR_MAX);
        hia = .5f * (width + width * high_a / rtengine::ColorToningParams::LABGRID_CORR_MAX);
        lob = .5f * (height + height * low_b / rtengine::ColorToningParams::LABGRID_CORR_MAX);
        hib = .5f * (height + height * high_b / rtengine::ColorToningParams::LABGRID_CORR_MAX);
        cr->set_line_width(2.);
        cr->set_source_rgb(0.6, 0.6, 0.6);
        cr->move_to(loa, lob);
        cr->line_to(hia, hib);
        cr->stroke();

        // drawing points
        cr->set_source_rgb(0.1, 0.1, 0.1);
        if (litPoint == LOW) {
            cr->arc(loa, lob, 5, 0, 2. * rtengine::RT_PI);
        } else {
            cr->arc(loa, lob, 3, 0, 2. * rtengine::RT_PI);
        }
        cr->fill();

        cr->set_source_rgb(0.9, 0.9, 0.9);
        if (litPoint == HIGH) {
            cr->arc(hia, hib, 5, 0, 2. * rtengine::RT_PI);
        } else {
            cr->arc(hia, hib, 3, 0, 2. * rtengine::RT_PI);
        }
        cr->fill();
    }

    copySurface(crf);
    return false;
}


bool LabGrid::on_button_press_event(GdkEventButton *event)
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


bool LabGrid::on_button_release_event(GdkEventButton *event)
{
    if (event->button == 1) {
        isDragged = false;
        return false;
    }
    return true;
}


bool LabGrid::on_motion_notify_event(GdkEventMotion *event)
{
    if (isDragged && delayconn.connected()) {
        delayconn.disconnect();
    }

    State oldLitPoint = litPoint;

    int width = get_allocated_width() - 2 * inset, height = get_allocated_height() - 2 * inset;
    const float mouse_x = std::min(std::max(event->x - inset, 0.), double(width));
    const float mouse_y = std::min(std::max(height - 1 - event->y + inset, 0.), double(height));
    const float ma = (2.0 * mouse_x - width) / (float)width;
    const float mb = (2.0 * mouse_y - height) / (float)height;
    if (isDragged) {
        if (litPoint == LOW) {
            low_a = ma * rtengine::ColorToningParams::LABGRID_CORR_MAX;
            low_b = mb * rtengine::ColorToningParams::LABGRID_CORR_MAX;
        } else if (litPoint == HIGH) {
            high_a = ma * rtengine::ColorToningParams::LABGRID_CORR_MAX;
            high_b = mb * rtengine::ColorToningParams::LABGRID_CORR_MAX;
        }
        edited = true;
        grab_focus();
        if (options.adjusterMinDelay == 0) {
            notifyListener();
        } else {
            delayconn = Glib::signal_timeout().connect(sigc::mem_fun(*this, &LabGrid::notifyListener), options.adjusterMinDelay);
        }
        queue_draw();
    } else {
        litPoint = NONE;
        float la = low_a / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        float lb = low_b / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        float ha = high_a / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        float hb = high_b / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        const float thrs = 0.05f;
        const float distlo = (la - ma) * (la - ma) + (lb - mb) * (lb - mb);
        const float disthi = (ha - ma) * (ha - ma) + (hb - mb) * (hb - mb);
        if (distlo < thrs * thrs && distlo < disthi) {
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


Gtk::SizeRequestMode LabGrid::get_request_mode_vfunc() const
{
    return Gtk::SIZE_REQUEST_HEIGHT_FOR_WIDTH;
}


void LabGrid::get_preferred_width_vfunc(int &minimum_width, int &natural_width) const
{
    minimum_width = 50;
    natural_width = 150;  // same as GRAPH_SIZE from mycurve.h
}


void LabGrid::get_preferred_height_for_width_vfunc(int width, int &minimum_height, int &natural_height) const
{
    minimum_height = natural_height = width;
}
