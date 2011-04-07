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
#include <mycurve.h>
#include <curves.h>
#include <string.h>
#include <gdkmm/types.h>

MyCurve::MyCurve () : listener(NULL) {

    cursor_type = CSArrow;
    innerWidth = get_allocation().get_width() - RADIUS * 2;
    innerHeight = get_allocation().get_height() - RADIUS * 2;
    prevInnerHeight = innerHeight;
    buttonPressed = false;
    snapTo = ST_None;
    colorProvider = NULL;
    sized = RS_Pending;

    set_extension_events(Gdk::EXTENSION_EVENTS_ALL);
    add_events(Gdk::EXPOSURE_MASK |	Gdk::POINTER_MOTION_MASK |	Gdk::POINTER_MOTION_HINT_MASK |	Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::BUTTON1_MOTION_MASK);

    mcih = new MyCurveIdleHelper;
    mcih->myCurve = this;
    mcih->destroyed = false;
    mcih->pending = 0;
}

MyCurve::~MyCurve () {

    if (mcih->pending)
        mcih->destroyed = true;
    else
        delete mcih;
}

void MyCurve::notifyListener () {

    if (listener)
        listener->curveChanged ();
}
