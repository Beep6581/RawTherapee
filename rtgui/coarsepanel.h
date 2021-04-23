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

#include "toolpanel.h"

class CoarsePanel final :
    public Gtk::Box,
    public ToolPanel
{

protected:
    Gtk::Button* rotate_left;
    Gtk::Button* rotate_right;
    Gtk::ToggleButton* hflip;
    Gtk::ToggleButton* vflip;
    int degree;
    bool oldhflip, oldvflip, degreechanged;

public:

    CoarsePanel ();

    void read               (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write              (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void initBatchBehavior  ();

    void rotateLeft     ();
    void rotateRight    ();
    void flipHorizontal ();
    void flipVertical   ();
};
