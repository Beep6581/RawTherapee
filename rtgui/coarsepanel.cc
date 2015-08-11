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
#include "coarsepanel.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

CoarsePanel::CoarsePanel () : ToolPanel ()
{

    degree = 0;
    degreechanged = true;

    Gtk::Image* rotateli = Gtk::manage (new RTImage ("stock-rotate-270.png"));
    rotate_left = Gtk::manage (new Gtk::Button ());
    rotate_left->add (*rotateli);
    rotate_left->set_relief(Gtk::RELIEF_NONE);
    pack_start (*rotate_left);

    Gtk::Image* rotateri = Gtk::manage (new RTImage ("stock-rotate-90.png"));
    rotate_right = Gtk::manage (new Gtk::Button ());
    rotate_right->add (*rotateri);
    rotate_right->set_relief(Gtk::RELIEF_NONE);
    pack_start (*rotate_right);

    Gtk::Image* fliphi = Gtk::manage (new RTImage ("stock-flip-horizontal.png"));
    hflip = Gtk::manage (new Gtk::ToggleButton ());
    hflip->add (*fliphi);
    hflip->set_relief(Gtk::RELIEF_NONE);
    pack_start (*hflip);

    Gtk::Image* flipvi = Gtk::manage (new RTImage ("stock-flip-vertical.png"));
    vflip = Gtk::manage (new Gtk::ToggleButton ());
    vflip->add (*flipvi);
    vflip->set_relief(Gtk::RELIEF_NONE);
    pack_start (*vflip);

    rotate_left->set_tooltip_markup (M("TP_COARSETRAF_TOOLTIP_ROTLEFT"));
    rotate_right->set_tooltip_markup (M("TP_COARSETRAF_TOOLTIP_ROTRIGHT"));
    vflip->set_tooltip_text (M("TP_COARSETRAF_TOOLTIP_VFLIP"));
    hflip->set_tooltip_text (M("TP_COARSETRAF_TOOLTIP_HFLIP"));

    rotate_left->signal_pressed().connect( sigc::mem_fun(*this, &CoarsePanel::rotateLeft) );
    rotate_right->signal_pressed().connect( sigc::mem_fun(*this, &CoarsePanel::rotateRight) );
    hflip->signal_toggled().connect( sigc::mem_fun(*this, &CoarsePanel::flipHorizontal) );
    vflip->signal_toggled().connect( sigc::mem_fun(*this, &CoarsePanel::flipVertical) );

    show_all_children ();
}

void CoarsePanel::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    degree = pp->coarse.rotate;

    if (pedited) {
        hflip->set_active (pedited->coarse.hflip ? pp->coarse.hflip : false);
        vflip->set_active (pedited->coarse.vflip ? pp->coarse.vflip : false);
        degreechanged = false;
        oldhflip = pp->coarse.hflip;
        oldvflip = pp->coarse.vflip;
    } else {
        hflip->set_active (pp->coarse.hflip);
        vflip->set_active (pp->coarse.vflip);
    }

    enableListener ();
}

void CoarsePanel::write (ProcParams* pp, ParamsEdited* pedited)
{

    if (pedited) {
        pedited->coarse.rotate = degreechanged;
        pedited->coarse.hflip = oldhflip != hflip->get_active ();
        pedited->coarse.vflip = oldvflip != vflip->get_active ();
    }

    pp->coarse.rotate = degree;
    pp->coarse.hflip = hflip->get_active ();
    pp->coarse.vflip = vflip->get_active ();
}

void CoarsePanel::initBatchBehavior ()
{

    disableListener ();

    degree = 0;
    hflip->set_active (false);
    vflip->set_active (false);

    enableListener ();
}

void CoarsePanel::rotateLeft ()
{

    //Rotate one way or the opposite depending if the image is already flipped or not
    if ( (vflip->get_active()) == (hflip->get_active ()) ) {
        degree = (degree + 270) % 360;
    } else {
        degree = (degree + 90) % 360;
    }

    degreechanged = true;

    if (listener) {
        listener->panelChanged (EvCTRotate, Glib::ustring::format (degree));
    }
}

void CoarsePanel::rotateRight ()
{

    //Rotate one way or the opposite depending if the image is already flipped or not
    if ( (vflip->get_active()) == (hflip->get_active ()) ) {
        degree = (degree + 90) % 360;
    } else {
        degree = (degree + 270) % 360;
    }

    degreechanged = true;

    if (listener) {
        listener->panelChanged (EvCTRotate, Glib::ustring::format (degree));
    }
}

void CoarsePanel::flipHorizontal ()
{

    if (listener) {
        if (hflip->get_active ()) {
            listener->panelChanged (EvCTHFlip, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvCTHFlip, M("GENERAL_DISABLED"));
        }
    }
}

void CoarsePanel::flipVertical   ()
{

    if (listener) {
        if (vflip->get_active ()) {
            listener->panelChanged (EvCTVFlip, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvCTVFlip, M("GENERAL_DISABLED"));
        }
    }
}


