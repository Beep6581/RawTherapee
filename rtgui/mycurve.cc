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
#include "mycurve.h"
#include "../rtengine/curves.h"
#include <cstring>
#include <gdkmm/types.h>

MyCurve::MyCurve () : pipetteR(-1.f), pipetteG(-1.f), pipetteB(-1.f), pipetteVal(-1.f), listener(NULL), cursor_type( CSArrow)
{

    graphX = get_allocation().get_width() - RADIUS * 2;
    graphY = get_allocation().get_height() - RADIUS * 2;
    prevGraphW = graphW;
    prevGraphH = graphH;
    buttonPressed = false;
    snapTo = ST_None;
    leftBar = NULL;
    bottomBar = NULL;
    colorProvider = NULL;
    sized = RS_Pending;
    snapToElmt = -100;
    curveIsDirty = true;
    edited_point = -1;

    set_extension_events(Gdk::EXTENSION_EVENTS_ALL);
#if defined (__APPLE__)
    // Workaround: disabling POINTER_MOTION_HINT_MASK as for gtk 2.24.22 the get_pointer() function is buggy for quartz and modifier mask is not updated correctly.
    // This workaround should be removed when bug is fixed in GTK2 or when migrating to GTK3
    add_events(Gdk::EXPOSURE_MASK | Gdk::POINTER_MOTION_MASK |  Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::BUTTON1_MOTION_MASK);
#else
    add_events(Gdk::EXPOSURE_MASK | Gdk::POINTER_MOTION_MASK |  Gdk::POINTER_MOTION_HINT_MASK | Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::BUTTON1_MOTION_MASK);
#endif

    mcih = new MyCurveIdleHelper;
    mcih->myCurve = this;
    mcih->destroyed = false;
    mcih->pending = 0;
}

MyCurve::~MyCurve ()
{

    if (mcih->pending) {
        mcih->destroyed = true;
    } else {
        delete mcih;
    }
}

int MyCurve::calcDimensions ()
{
    int newRequestedW, newRequestedH;

    newRequestedW = newRequestedH = get_allocation().get_width();

    if (leftBar && !bottomBar) {
        newRequestedH -= getBarWidth() + CBAR_MARGIN - RADIUS;
    }

    if (!leftBar && bottomBar) {
        newRequestedH += getBarWidth() + CBAR_MARGIN - RADIUS;
    }

    graphW = newRequestedW - RADIUS - (leftBar   ? (getBarWidth() + CBAR_MARGIN) : RADIUS);
    graphH = newRequestedH - RADIUS - (bottomBar ? (getBarWidth() + CBAR_MARGIN) : RADIUS);
    graphX = newRequestedW - RADIUS - graphW;
    graphY = RADIUS + graphH;

    return newRequestedH;
}

void MyCurve::setColoredBar (ColoredBar *left, ColoredBar *bottom)
{
    leftBar = left;
    bottomBar = bottom;
}

void MyCurve::notifyListener ()
{

    if (listener) {
        listener->curveChanged ();
    }
}

bool MyCurve::snapCoordinateX(double testedVal, double realVal)
{

    double dist = realVal - testedVal;

    if (dist < 0.) {
        dist = -dist;
    }

    if (dist < snapToMinDistX) {
        snapToMinDistX = dist;
        snapToValX = testedVal;
        return true;
    }

    return false;
}

bool MyCurve::snapCoordinateY(double testedVal, double realVal)
{

    double dist = realVal - testedVal;

    if (dist < 0.) {
        dist = -dist;
    }

    if (dist < snapToMinDistY) {
        snapToMinDistY = dist;
        snapToValY = testedVal;
        return true;
    }

    return false;
}

float MyCurve::getVal(LUTf &curve, int x)
{
    if ((graphW - 2) == curve.getSize()) {
        return curve[x];
    } else {
        return curve.getVal01(float(x) / (graphW - 3));
    }
}

void MyCurve::on_style_changed (const Glib::RefPtr<Gtk::Style>& style)
{
    setDirty(true);
    queue_draw ();
}

void MyCurve::refresh()
{
    if (leftBar != NULL) {
        leftBar->setDirty(true);
    }

    if (bottomBar != NULL) {
        bottomBar->setDirty(true);
    }

    setDirty(true);

    Glib::RefPtr<Gdk::Window> win = get_window();

    if (win) {
        win->invalidate(true);
    }
}
