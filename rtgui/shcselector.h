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
#ifndef _SHCSELECTOR_
#define _SHCSELECTOR_

#include <gtkmm.h>
#include "coloredbar.h"

class SHCListener
{
public:
    virtual ~SHCListener() {}
    virtual void shcChanged () {}
};

class SHCSelector : public Gtk::DrawingArea, public ColoredBar
{

protected:

    int movingPosition;
    double tmpX, tmpPos;

    double defaults[3];
    double positions[3];
    double wslider;

    // left margin, essentially a workaround to take care of an eventual right colored bar (e.g. for curves)
    int leftMargin;
    // right margin, essentially a workaround to take care of an eventual right colored bar
    int rightMargin;

    const static int hb = 3;  // horizontal border
    const static int vb = 2;  // vertical border

    SHCListener* cl;

public:

    SHCSelector();

    void setSHCListener (SHCListener* l)
    {
        cl = l;;
    }

    void setMargins(int left, int right);
    void setDefaults (double spos, double cpos, double hpos);
    void setPositions (double spos, double cpos, double hpos);
    void getPositions (double& spos, double& cpos, double& hpos);
    void on_realize();
    bool on_expose_event(GdkEventExpose* event);
    bool on_button_press_event (GdkEventButton* event);
    bool on_button_release_event (GdkEventButton* event);
    bool on_motion_notify_event (GdkEventMotion* event);
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
    bool reset ();
    void refresh();
};

#endif

