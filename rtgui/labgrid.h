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


class LabGridArea final : public Gtk::DrawingArea {
private:
    rtengine::ProcEvent evt;
    Glib::ustring evtMsg;

    enum State { NONE, HIGH, LOW, GRE};
    State litPoint;
    double low_a;
    double high_a;
    double low_b;
    double high_b;
    double gre_x;
    double gre_y;
    double whi_x;
    double whi_y;
    double me_x;
    double me_y;
    double ghs_x6;
    double ghs_y6;
    double ghs_x7;
    double ghs_y7;
    double ghs_x8;
    double ghs_y8;
    double ghs_x9;
    double ghs_y9;
    double ghs_x10;
    double ghs_y10;
    double ghs_x11;
    double ghs_y11;//+4 12 11
    double defaultLow_a;
    double defaultHigh_a;
    double defaultLow_b;
    double defaultHigh_b;
    double defaultgre_x;
    double defaultgre_y;
    double defaultwhi_x;
    double defaultwhi_y;
    double defaultme_x;
    double defaultme_y;
    double default_gsx6;//added for GHS 
    double default_gsy6;
    double default_gsx7;
    double default_gsy7;
    double default_gsx8;
    double default_gsy8;
    double default_gsx9;
    double default_gsy9;
    double default_gsx10;//+4 12 11
    double default_gsy10;
    double default_gsx11;
    double default_gsy11;

    ToolPanelListener *listener;
    bool edited;
    bool isDragged;
    sigc::connection delayconn;
    static const int inset = 5;

    bool low_enabled;
    bool ciexy_enabled;
    bool ghs_enabled;
    bool mous_enabled;

    bool notifyListener();
    void getLitPoint();

public:
    LabGridArea(rtengine::ProcEvent evt, const Glib::ustring &msg, bool enable_low=true, bool ciexy=false, bool ghs=false, bool mous=false);

    void getParams(double &la, double &lb, double &ha, double &hb, double &gx, double &gy, double &wx, double &wy, double &mx, double &my, 
        double &gx6, double &gy6, double &gx7, double &gy7, double &gx8, double &gy8, double &gx9, double &gy9, double &gx10, double &gy10, double &gx11, double &gy11) const;//+4 12 11
    void setParams(double la, double lb, double ha, double hb, double gx, double gy, double wx, double wy, double mx, double my,  double gx6, double gy6, double gx7, double gy7, double gx8, double gy8, double gx9, double gy9, double gx10, double gy10, double gx11, double gy11, bool notify);//+4 12 11
    void setDefault (double la, double lb, double ha, double hb, double gx, double gy, double wx, double wy, double mx, double my, double gx6, double gy6, double gx7, double gy7, double gx8, double gy8, double gx9, double gy9, double gx10, double gy10, double gx11, double gy11);//+4 12 11
    void setEdited(bool yes);
    bool getEdited() const;
    void reset(bool toInitial);
    void setListener(ToolPanelListener *l);

    bool lowEnabled() const;
    void setLowEnabled(bool yes);
    bool ciexyEnabled() const;
    void setciexyEnabled(bool yes);
    bool ghsEnabled() const;
    void setghsEnabled(bool yes);
    bool mousEnabled() const;
    void setmousEnabled(bool yes);

    bool on_draw(const ::Cairo::RefPtr<Cairo::Context> &cr) override;
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
    LabGrid(rtengine::ProcEvent evt, const Glib::ustring &msg, bool enable_low=true, bool ciexy=false, bool ghs=false, bool mous=true);

    void getParams(double &la, double &lb, double &ha, double &hb, double &gx, double &gy, double &wx, double &wy, double &mx, double &my, double &gx6, double &gy6, double &gx7, double &gy7, double &gx8, double &gy8, double &gx9, double &gy9, double &gx10, double &gy10, double &gx11, double &gy11) 
                const { return grid.getParams(la, lb, ha, hb, gx, gy, wx, wy, mx, my, gx6, gy6, gx7, gy7, gx8, gy8, gx9, gy9, gx10, gy10, gx11, gy11); }//+4 12 11
    void setParams(double la, double lb, double ha, double hb, double gx, double gy, double wx, double wy, double mx, double my, double gx6, double gy6, double gx7, double gy7, double gx8, double gy8, double gx9, double gy9, double gx10, double gy10, double gx11, double gy11, bool notify) 
                { grid.setParams(la, lb, ha, hb, gx, gy, wx, wy, mx, my, gx6, gy6, gx7, gy7, gx8, gy8, gx9, gy9, gx10, gy10, gx11, gy11, notify); }//+4 12 11
    void setDefault (double la, double lb, double ha, double hb, double gx, double gy, double wx, double wy, double mx, double my, double gx6, double gy6, double gx7, double gy7, double gx8, double gy8, double gx9, double gy9, double gx10, double gy10, double gx11, double gy11)
                { grid.setDefault(la, lb, ha, hb, gx, gy, wx, wy, mx, my, gx6, gy6, gx7, gy7, gx8, gy8, gx9, gy9, gx10, gy10, gx11, gy11); }//+4 12 11
    void setEdited(bool yes) { grid.setEdited(yes); }
    bool getEdited() const { return grid.getEdited(); }
    void reset(bool toInitial) { grid.reset(toInitial); }
    void setListener(ToolPanelListener *l) { grid.setListener(l); }
    bool lowEnabled() const { return grid.lowEnabled(); }
    void setLowEnabled(bool yes) { grid.setLowEnabled(yes); }
    bool ciexyEnabled() const { return grid.ciexyEnabled(); }
    void setciexyEnabled(bool yes) { grid.setciexyEnabled(yes); }
    bool ghsEnabled() const { return grid.ghsEnabled(); }
    void setghsEnabled(bool yes) { grid.setghsEnabled(yes); }
    bool mousEnabled() const { return grid.mousEnabled(); }
    void setmousEnabled(bool yes) { grid.setmousEnabled(yes); }

};

