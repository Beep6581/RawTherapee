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


bool LabGrid::notify_listener()
{
    if (listener_) {
        listener_->panelChanged(evt_, Glib::ustring::compose("%1 %2 %3 %4", int(low_a_), int(low_b_), int(high_a_), int(high_b_)));
    }
    return false;
}
    

LabGrid::LabGrid(rtengine::ProcEvent evt):
    Gtk::DrawingArea(),
    evt_(evt), selected_(OFF),
    low_a_(0.f), high_a_(0.f), low_b_(0.f), high_b_(0.f),
    listener_(nullptr),
    edited_(false)
{
    set_can_focus(true);
    add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::POINTER_MOTION_MASK);
}


void LabGrid::get_params(double &la, double &lb, double &ha, double &hb) const
{
    la = low_a_;
    ha = high_a_;
    lb = low_b_;
    hb = high_b_;
}


void LabGrid::set_params(double la, double lb, double ha, double hb, bool notify)
{
    low_a_ = la;
    low_b_ = lb;
    high_a_ = ha;
    high_b_ = hb;
    queue_draw();
    if (notify) {
        notify_listener();
    }
}


void LabGrid::set_edited(bool yes)
{
    edited_ = yes;
}


bool LabGrid::get_edited() const
{
    return edited_;
}


void LabGrid::set_listener(ToolPanelListener *l)
{
    listener_ = l;
}


bool LabGrid::on_draw(const ::Cairo::RefPtr<Cairo::Context> &crf)
{
    int width = get_allocated_width();
    int height = get_allocated_height();
    Cairo::RefPtr<Cairo::ImageSurface> cst =
        Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, width, height);
    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(cst);
    // clear bg
    cr->set_source_rgb(.2, .2, .2);
    cr->paint();
    
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
    cr->set_antialias(Cairo::ANTIALIAS_DEFAULT);
    float loa, hia, lob, hib;
    loa = .5f * (width + width * low_a_ / rtengine::ColorToningParams::LABGRID_CORR_MAX);
    hia = .5f * (width + width * high_a_ / rtengine::ColorToningParams::LABGRID_CORR_MAX);
    lob = .5f * (height + height * low_b_ / rtengine::ColorToningParams::LABGRID_CORR_MAX);
    hib = .5f * (height + height * high_b_ / rtengine::ColorToningParams::LABGRID_CORR_MAX);
    cr->set_line_width(2.);
    cr->set_source_rgb(0.6, 0.6, 0.6);
    cr->move_to(loa, lob);
    cr->line_to(hia, hib);
    cr->stroke();

    cr->set_source_rgb(0.1, 0.1, 0.1);
    if (selected_ == LOW) {
        cr->arc(loa, lob, 5, 0, 2. * rtengine::RT_PI);
    } else {
        cr->arc(loa, lob, 3, 0, 2. * rtengine::RT_PI);
    }
    cr->fill();

    cr->set_source_rgb(0.9, 0.9, 0.9);
    if (selected_ == HIGH) {
        cr->arc(hia, hib, 5, 0, 2. * rtengine::RT_PI);
    } else {
        cr->arc(hia, hib, 3, 0, 2. * rtengine::RT_PI);
    }
    cr->fill();

    crf->set_source(cst, 0, 0);
    crf->paint();
        
    return true;
}


bool LabGrid::on_button_press_event(GdkEventButton *event)
{
    if (event->button == 1 && event->type == GDK_2BUTTON_PRESS) {
        switch (selected_) {
        case OFF:
            low_a_ = low_b_ = high_a_ = high_b_ = 0.f;
            break;
        case LOW:
            low_a_ = low_b_ = 0.f;
            break;
        case HIGH:
            high_a_ = high_b_ = 0.f;
            break;
        }
        edited_ = true;
        queue_draw();
        notify_listener();
    }
    return true;
}
    

bool LabGrid::on_motion_notify_event(GdkEventMotion *event)
{
    if (delay_conn_.connected()) {
        delay_conn_.disconnect();
    }
        
    int width = get_allocated_width() - 2 * inset, height = get_allocated_height() - 2 * inset;
    const float mouse_x = std::min(std::max(event->x - inset, 0.), double(width));
    const float mouse_y = std::min(std::max(height - 1 - event->y + inset, 0.), double(height));
    const float ma = (2.0 * mouse_x - width) / (float)width;
    const float mb = (2.0 * mouse_y - height) / (float)height;
    bool refresh = selected_ != OFF;
    if (event->state & GDK_BUTTON1_MASK) {
        if (selected_ == LOW) {
            low_a_ = ma * rtengine::ColorToningParams::LABGRID_CORR_MAX;
            low_b_ = mb * rtengine::ColorToningParams::LABGRID_CORR_MAX;
        } else if (selected_ == HIGH) {
            high_a_ = ma * rtengine::ColorToningParams::LABGRID_CORR_MAX;
            high_b_ = mb * rtengine::ColorToningParams::LABGRID_CORR_MAX;
        }
    } else {
        selected_ = OFF;
        float la = low_a_ / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        float lb = low_b_ / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        float ha = high_a_ / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        float hb = high_b_ / rtengine::ColorToningParams::LABGRID_CORR_MAX;
        const float thrs = 0.05f;
        const float distlo = (la - ma) * (la - ma) + (lb - mb) * (lb - mb);
        const float disthi = (ha - ma) * (ha - ma) + (hb - mb) * (hb - mb);
        if (distlo < thrs * thrs && distlo < disthi) {
            selected_ = LOW;
        } else if (disthi < thrs * thrs && disthi <= distlo) {
            selected_ = HIGH;
        }
    }
    if (selected_ != OFF) {
        grab_focus();
        if (options.adjusterMinDelay == 0) {
            notify_listener();
        } else {
            delay_conn_ = Glib::signal_timeout().connect(sigc::mem_fun(*this, &LabGrid::notify_listener), options.adjusterMinDelay);
        }
    }
    if (refresh) {
        queue_draw();
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
    natural_width = 200;
}


void LabGrid::get_preferred_height_for_width_vfunc(int width, int &minimum_height, int &natural_height) const
{
    minimum_height = natural_height = width;
}
