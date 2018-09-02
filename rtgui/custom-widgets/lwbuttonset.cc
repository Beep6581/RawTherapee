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
#include "lwbuttonset.h"

LWButtonSet::LWButtonSet () : aw(0), ah(0), ax(0), ay(0)
{
}

LWButtonSet::~LWButtonSet ()
{

    for (size_t i = 0; i < buttons.size(); i++) {
        delete buttons[i];
    }
}

void LWButtonSet::add (LWButton* b)
{

    buttons.push_back (b);
}

void LWButtonSet::getMinimalDimensions (int& w, int& h)
{

    w = 0;
    h = 0;

    for (size_t i = 0; i < buttons.size(); i++) {
        int bw, bh;
        buttons[i]->getSize (bw, bh);
        w += bw;

        if (bh > h) {
            h = bh;
        }
    }
}

void LWButtonSet::arrangeButtons (int x, int y, int w, int h)
{

    int mw, mh;
    getMinimalDimensions (mw, mh);

    if (w < 0) {
        w = mw;
    }

    if (h < 0) {
        h = mh;
    }

    int begx = x;
    int endx = x + w - 1;

    for (size_t i = 0; i < buttons.size(); i++) {
        LWButton::Alignment halign, valign;
        int bx = 0, by = 0, bw = 0, bh = 0;
        buttons[i]->getSize (bw, bh);
        buttons[i]->getAlignment (halign, valign);

        if (halign == LWButton::Left) {
            bx = begx;
            begx += bw;
        } else if (halign == LWButton::Right) {
            bx = endx - bw;
            endx -= bw;
        }

        if (valign == LWButton::Top) {
            by = y;
        } else if (valign == LWButton::Bottom) {
            by = y + h - bh - 1;
        } else if (valign == LWButton::Center) {
            by = y + (h - bh) / 2;
        }

        buttons[i]->setPosition (bx, by);
    }

    aw = w;
    ah = h;
    ax = x;
    ay = y;
}

void LWButtonSet::move (int nx, int ny)
{

    for (size_t i = 0; i < buttons.size(); i++) {
        int x, y;
        buttons[i]->getPosition (x, y);
        buttons[i]->setPosition (x + nx - ax, y + ny - ay);
    }

    ax = nx;
    ay = ny;
}

void LWButtonSet::redraw (Cairo::RefPtr<Cairo::Context> context)
{

    for (size_t i = 0; i < buttons.size(); i++) {
        buttons[i]->redraw (context);
    }
}

bool LWButtonSet::motionNotify (int x, int y)
{

    bool res = false;

    for (size_t i = 0; i < buttons.size(); i++) {
        bool handled = buttons[i]->motionNotify (x, y);
        res = res || handled;
    }

    return res;
}

bool LWButtonSet::pressNotify (int x, int y)
{

    bool res = false;

    for (size_t i = 0; i < buttons.size(); i++) {
        bool handled = buttons[i]->pressNotify (x, y);
        res = res || handled;
    }

    return res;
}

bool LWButtonSet::releaseNotify (int x, int y)
{

    bool res = false;

    for (size_t i = 0; i < buttons.size(); i++) {
        bool handled = buttons[i]->releaseNotify (x, y);
        res = res || handled;
    }

    return res;
}

bool LWButtonSet::inside (int x, int y)
{

    for (size_t i = 0; i < buttons.size(); i++)
        if (buttons[i]->inside (x, y)) {
            return true;
        }

    return false;
}

void LWButtonSet::setButtonListener (LWButtonListener* bl)
{

    for (size_t i = 0; i < buttons.size(); i++) {
        buttons[i]->setButtonListener (bl);
    }
}

void LWButtonSet::getAllocatedDimensions (int& w, int& h)
{

    w = aw;
    h = ah;
}

void LWButtonSet::setColors (const Gdk::RGBA& bg, const Gdk::RGBA& fg)
{

    for (size_t i = 0; i < buttons.size(); i++) {
        buttons[i]->setColors (bg, fg);
    }
}

Glib::ustring LWButtonSet::getToolTip (int x, int y)
{

    for (size_t i = 0; i < buttons.size(); i++) {
        Glib::ustring ttip = buttons[i]->getToolTip (x, y);

        if (ttip != "") {
            return ttip;
        }
    }

    return "";
}
