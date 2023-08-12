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

#include <gtkmm.h>

#include "delayed.h"
#include "options.h"
#include "pointermotionlistener.h"

class PreviewWindow;

class Navigator final :
    public Gtk::Frame,
    public PointerMotionListener
{

    typedef const double (*TMatrix)[3];

private:
    DelayedCall<bool, const rtengine::procparams::ColorManagementParams *, int, int, int, int, int, bool> pointer_moved_delayed_call;

    Options::NavigatorUnit currentRGBUnit;
    Options::NavigatorUnit currentHSVUnit;
    void cycleUnitsRGB (GdkEventButton *event);
    void cycleUnitsHSV (GdkEventButton *event);

protected:
    Gtk::Label* dimension;
    Gtk::Label* position;
    Gtk::Label *R, *G, *B;
    Gtk::Label *H, *S, *V;
    Gtk::Label *LAB_A, *LAB_B, *LAB_L;

    Gtk::Label *lR, *lG, *lB;
    Gtk::Label *lH, *lS, *lV;
    Gtk::Label *lLAB_A, *lLAB_B, *lLAB_L;


public:
    PreviewWindow* previewWindow;

    Navigator();
    ~Navigator() override;

    // pointermotionlistener interface
    //  void pointerMoved (bool validPos, int x, int y, int r, int g, int b);
    void pointerMoved(bool validPos, const rtengine::procparams::ColorManagementParams &cmp, int x, int y, int r, int g, int b, bool raw = false) override;
    void setInvalid (int fullWidth = -1, int fullHeight = -1);

    void getRGBText (int r, int g, int b, Glib::ustring &sR, Glib::ustring &sG, Glib::ustring &sB, bool isRaw = false) override;
    void getHSVText (float h, float s, float v, Glib::ustring &sH, Glib::ustring &sS, Glib::ustring &sV) override;
    void getLABText (float l, float a, float b, Glib::ustring &sL, Glib::ustring &sA, Glib::ustring &sB) override;

};
