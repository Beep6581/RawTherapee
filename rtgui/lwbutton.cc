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
#include "lwbutton.h"
#include "guiutils.h"
#include "rtsurface.h"

LWButton::LWButton (Cairo::RefPtr<RTSurface> i, int aCode, void* aData, Alignment ha, Alignment va, Glib::ustring* tooltip)
    : xpos(0), ypos(0), halign(ha), valign(va), icon(i), bgr(0.0), bgg(0.0), bgb(0.0), fgr(0.0), fgg(0.0), fgb(0.0), state(Normal), listener(nullptr), actionCode(aCode), actionData(aData), toolTip(tooltip)
{

    if (i)  {
        w = i->getWidth () + 2;
        h = i->getHeight () + 2;
    } else {
        w = h = 2;
    }
}

void LWButton::getSize (int& minw, int& minh) const 
{

    minw = w;
    minh = h;
}

void LWButton::setPosition (int x, int y)
{

    xpos = x;
    ypos = y;
}

void LWButton::addPosition (int x, int y)
{
    xpos += x;
    ypos += y;
}

void LWButton::getPosition (int& x, int& y) const
{

    x = xpos;
    y = ypos;
}

void LWButton::setIcon (Cairo::RefPtr<RTSurface> i)
{

    icon = i;

    if (i)  {
        w = i->getWidth () + 2;
        h = i->getHeight () + 2;
    } else {
        w = h = 2;
    }
}

Cairo::RefPtr<RTSurface> LWButton::getIcon () const
{

    return icon;
}

void LWButton::setColors (const Gdk::RGBA& bg, const Gdk::RGBA& fg)
{

    bgr = bg.get_red ();
    bgg = bg.get_green ();
    bgb = bg.get_blue ();
    fgr = fg.get_red ();
    fgg = fg.get_green ();
    fgb = fg.get_blue ();
}

bool LWButton::inside (int x, int y) const
{

    return x > xpos && x < xpos + w && y > ypos && y < ypos + h;
}

bool LWButton::motionNotify  (int x, int y)
{

    bool in = inside (x, y);
    State nstate = state;

    if (state == Normal && in) {
        nstate = Over;
    } else if (state == Over && !in) {
        nstate = Normal;
    } else if (state == Pressed_In && !in) {
        nstate = Pressed_Out;
    } else if (state == Pressed_Out && in) {
        nstate = Pressed_In;
    }

    if (state != nstate) {
        state = nstate;

        if (listener) {
            listener->redrawNeeded (this);
        }

        return true;
    }

    return in;
}

bool LWButton::pressNotify   (int x, int y)
{

    bool in = inside (x, y);
    State nstate = state;

    if (in && (state == Normal || state == Over || state == Pressed_Out)) {
        nstate = Pressed_In;
    } else if (!in && state == Pressed_In) {
        nstate = Normal;
    }

    if (state != nstate) {
        state = nstate;

        if (listener) {
            listener->redrawNeeded (this);
        }

        return true;
    }

    return in;
}

bool LWButton::releaseNotify (int x, int y)
{

    bool in = inside (x, y);
    State nstate;
    bool action = false;

    if (in && (state == Pressed_In || state == Pressed_Out)) {
        nstate = Over;
        action = true;
    } else {
        nstate = Normal;
    }

    bool ret = action;

    if (state != nstate) {
        state = nstate;

        if (listener) {
            listener->redrawNeeded (this);
        }

        ret = true;
    }

    if (action && listener) {
        listener->buttonPressed (this, actionCode, actionData);
    }

    return ret;
}

void LWButton::redraw (Cairo::RefPtr<Cairo::Context> context)
{

    GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected
    context->set_line_width (1.0);
    context->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    context->rectangle (xpos + 0.5, ypos + 0.5, w - 1.0, h - 1.0);

    if (state == Pressed_In) {
        context->set_source_rgb (fgr, fgg, fgb);
    } else {
        context->set_source_rgba (bgr, bgg, bgb, 0);
    }

    context->fill_preserve ();

    if (state == Over) {
        context->set_source_rgb (fgr, fgg, fgb);
    } else {
        context->set_source_rgba (bgr, bgg, bgb, 0);
    }

    context->stroke ();
    int dilat = 1;

    if (state == Pressed_In) {
        dilat++;
    }

    if (icon) {
        context->set_source (icon->get(), xpos + dilat, ypos + dilat);
        context->paint ();
    }
}

void LWButton::getAlignment (Alignment& ha, Alignment& va) const
{

    ha = halign;
    va = valign;
}

Glib::ustring LWButton::getToolTip (int x, int y) const
{
    if (inside(x, y) && toolTip) {
        return *toolTip;
    } else {
        return {};
    }
}

void LWButton::setToolTip (Glib::ustring* tooltip)
{

    toolTip = tooltip;
}

