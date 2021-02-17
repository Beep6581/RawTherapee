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

#pragma once

#include <gtkmm.h>
#include "eventmapper.h"
#include "toolpanel.h"


class LabGridArea final : public Gtk::DrawingArea, public BackBuffer {
private:
    rtengine::ProcEvent evt;
    Glib::ustring evtMsg;
    
    enum State { NONE, HIGH, LOW };
    State litPoint;
    double low_a;
    double high_a;
    double low_b;
    double high_b;

    double defaultLow_a;
    double defaultHigh_a;
    double defaultLow_b;
    double defaultHigh_b;

    ToolPanelListener *listener;
    bool edited;
    bool isDragged;
    sigc::connection delayconn;
    static const int inset = 5;

    bool low_enabled;

    bool notifyListener();
    void getLitPoint();

public:
    LabGridArea(rtengine::ProcEvent evt, const Glib::ustring &msg, bool enable_low=true);

    void getParams(double &la, double &lb, double &ha, double &hb) const;
    void setParams(double la, double lb, double ha, double hb, bool notify);
    void setDefault (double la, double lb, double ha, double hb);
    void setEdited(bool yes);
    bool getEdited() const;
    void reset(bool toInitial);
    void setListener(ToolPanelListener *l);

    bool lowEnabled() const;
    void setLowEnabled(bool yes);

    bool on_draw(const ::Cairo::RefPtr<Cairo::Context> &crf) override;
    void on_style_updated () override;
    bool on_button_press_event(GdkEventButton *event) override;
    bool on_button_release_event(GdkEventButton *event) override;
    bool on_motion_notify_event(GdkEventMotion *event) override;
    Gtk::SizeRequestMode get_request_mode_vfunc() const override;
    void get_preferred_width_vfunc(int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
};


class LabGrid: public Gtk::Box {
private:
    LabGridArea grid;

    bool resetPressed(GdkEventButton *event);
    
public:
    LabGrid(rtengine::ProcEvent evt, const Glib::ustring &msg, bool enable_low=true);

    void getParams(double &la, double &lb, double &ha, double &hb) const { return grid.getParams(la, lb, ha, hb); }
    void setParams(double la, double lb, double ha, double hb, bool notify) { grid.setParams(la, lb, ha, hb, notify); }
    void setDefault (double la, double lb, double ha, double hb) { grid.setDefault(la, lb, ha, hb); }
    void setEdited(bool yes) { grid.setEdited(yes); }
    bool getEdited() const { return grid.getEdited(); }
    void reset(bool toInitial) { grid.reset(toInitial); }
    void setListener(ToolPanelListener *l) { grid.setListener(l); }
    bool lowEnabled() const { return grid.lowEnabled(); }
    void setLowEnabled(bool yes) { grid.setLowEnabled(yes); }
};

