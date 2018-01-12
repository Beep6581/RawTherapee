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

#pragma once

#include <gtkmm.h>
#include "eventmapper.h"
#include "toolpanel.h"


class LabGrid: public Gtk::DrawingArea {
private:
    rtengine::ProcEvent evt_;
    enum State { OFF, HIGH, LOW };
    State selected_;
    float low_a_;
    float high_a_;
    float low_b_;
    float high_b_;
    ToolPanelListener *listener_;
    bool edited_;
    sigc::connection delay_conn_;
    static const int inset = 2;

    bool notify_listener();
    
public:
    LabGrid(rtengine::ProcEvent evt);

    void get_params(double &la, double &lb, double &ha, double &hb) const;
    void set_params(double la, double lb, double ha, double hb, bool notify);
    void set_edited(bool yes);
    bool get_edited() const;
    void set_listener(ToolPanelListener *l);

    bool on_draw(const ::Cairo::RefPtr<Cairo::Context> &crf);
    bool on_button_press_event(GdkEventButton *event);
    bool on_motion_notify_event(GdkEventMotion *event);
    Gtk::SizeRequestMode get_request_mode_vfunc() const;
    void get_preferred_width_vfunc(int &minimum_width, int &natural_width) const;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const;
};

