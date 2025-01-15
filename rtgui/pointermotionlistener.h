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
#pragma once

#include <glibmm/ustring.h>
#include <sigc++/sigc++.h>


namespace rtengine
{

namespace procparams
{

struct ColorManagementParams;

}

}


class PointerMotionListener
{
protected:
    sigc::signal<void> sig_cycle_rgb;
    sigc::signal<void> sig_cycle_hsv;

public:
    virtual ~PointerMotionListener() = default;
    virtual void pointerMoved(bool validPos, const rtengine::procparams::ColorManagementParams &cmp, int x, int y, int r, int g, int b, bool isRaw = false) = 0;
    virtual void getRGBText (int r, int g, int b, Glib::ustring &sR, Glib::ustring &sG, Glib::ustring &sB, bool isRaw = false) { sR = "--"; sG = "--"; sB = "--"; }
    virtual void getHSVText (float h, float s, float v, Glib::ustring &sH, Glib::ustring &sS, Glib::ustring &sV) { sH = "--"; sS = "--"; sV = "--"; }
    virtual void getLABText (float l, float a, float b, Glib::ustring &sL, Glib::ustring &sA, Glib::ustring &sB) { sL = "--"; sA = "--"; sB = "--"; }

    sigc::signal<void> signal_cycle_rgb()
    {
        return sig_cycle_rgb;
    }
    sigc::signal<void> signal_cycle_hsv()
    {
        return sig_cycle_hsv;
    }
};
