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
#ifndef _LWBUTTONSET_
#define _LWBUTTONSET_

#include <gtkmm.h>
#include "lwbutton.h"
#include <vector>

class LWButtonSet
{

protected:
    std::vector<LWButton*> buttons;
    int aw, ah, ax, ay;
public:
    LWButtonSet ();
    ~LWButtonSet ();

    void add (LWButton* b);

    void    getMinimalDimensions (int& w, int& h);
    void    getAllocatedDimensions (int& w, int& h);
    void    arrangeButtons (int x, int y, int w, int h);
    void    setColors     (const Gdk::RGBA& bg, const Gdk::RGBA& fg);
    bool    motionNotify  (int x, int y);
    bool    pressNotify   (int x, int y);
    bool    releaseNotify (int x, int y);
    void    move          (int nx, int ny);
    bool    inside        (int x, int y);

    Glib::ustring getToolTip (int x, int y);

    void    setButtonListener   (LWButtonListener* bl);
    void    redraw              (Cairo::RefPtr<Cairo::Context> context);
};

#endif
