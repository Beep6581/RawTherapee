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
#include "lwbutton.h"
#include "guiutils.h"

LWButton::LWButton (Cairo::RefPtr<Cairo::ImageSurface> i, int aCode, void* aData, Alignment ha, Alignment va, Glib::ustring tooltip)
    : halign(ha), valign(va), icon(i), state(Normal), listener(NULL), actionCode(aCode), actionData(aData), toolTip(tooltip)
{

    if (i)  {
        w = i->get_width () + 2;
        h = i->get_height () + 2;
    } else {
        w = h = 2;
    }
}

void LWButton::getSize (int& minw, int& minh)
{

    minw = w;
    minh = h;
}

void LWButton::setPosition (int x, int y)
{

    xpos = x;
    ypos = y;
}

void LWButton::getPosition (int& x, int& y)
{

    x = xpos;
    y = ypos;
}

void LWButton::setIcon (Cairo::RefPtr<Cairo::ImageSurface> i)
{

    icon = i;

    if (i)  {
        w = i->get_width () + 2;
        h = i->get_height () + 2;
    } else {
        w = h = 2;
    }
}

Cairo::RefPtr<Cairo::ImageSurface> LWButton::getIcon ()
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

bool LWButton::inside (int x, int y)
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
    State nstate = state;
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

    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
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
        context->set_source (icon, xpos + dilat, ypos + dilat);
        context->paint ();
    }
}

void LWButton::getAlignment (Alignment& ha, Alignment& va)
{

    ha = halign;
    va = valign;
}

Glib::ustring LWButton::getToolTip (int x, int y)
{

    if (inside (x, y)) {
        return toolTip;
    } else {
        return "";
    }
}

void LWButton::setToolTip (const Glib::ustring& tooltip)
{

    toolTip = tooltip;
}

