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
#include <rotate.h>
#include <iomanip>
#include <guiutils.h>

extern Glib::ustring argv0;

using namespace rtengine;
using namespace rtengine::procparams;

Rotate::Rotate () : ToolPanel (), degAdd(false) {

  rlistener = NULL;

  degree = Gtk::manage (new Adjuster (M("TP_ROTATE_DEGREE"), -45, 45, 0.01, 0));
  degree->setAdjusterListener (this); 
  pack_start (*degree);

  fill = Gtk::manage (new Gtk::CheckButton (M("TP_ROTATE_FILL")));
  pack_start (*fill);

  selectStraight = Gtk::manage (new Gtk::Button (M("TP_ROTATE_SELECTLINE")));
  Gtk::Image* selimg = Gtk::manage (new Gtk::Image (argv0+"/images/straighten16.png"));
  selectStraight->set_image (*selimg);
  pack_start (*selectStraight, Gtk::PACK_SHRINK, 2);

  autoCrop = Gtk::manage (new Gtk::Button (M("TP_ROTATE_AUTOCROP")));
  pack_start (*autoCrop, Gtk::PACK_SHRINK, 2); 

  selectStraight->signal_pressed().connect( sigc::mem_fun(*this, &Rotate::selectStraightPressed) );
  autoCrop->signal_pressed().connect( sigc::mem_fun(*this, &Rotate::autoCropPressed) );
  fillConn = fill->signal_toggled().connect( sigc::mem_fun(*this, &Rotate::fillPressed) );

  fill->set_active (true);
  show_all ();
}

void Rotate::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        degree->setEditedState (pedited->rotate.degree ? Edited : UnEdited);
        fill->set_inconsistent (!pedited->rotate.fill);
    }

    degree->setValue (pp->rotate.degree);
    fillConn.block (true);
    fill->set_active (pp->rotate.fill);
    fillConn.block (false);
    
    lastFill = pp->rotate.fill;
    
    enableListener ();
}

void Rotate::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->rotate.degree = degree->getValue ();
    pp->rotate.fill   = fill->get_active ();

    if (pedited) {
        pedited->rotate.degree = degree->getEditedState ();
        pedited->rotate.fill   = !fill->get_inconsistent();
    }
}

void Rotate::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    degree->setDefault (defParams->rotate.degree);

    if (pedited) 
        degree->setDefaultEditedState (pedited->rotate.degree ? Edited : UnEdited);
    else 
        degree->setDefaultEditedState (Irrelevant);
}

void Rotate::adjusterChanged (Adjuster* a, double newval) {

    if (listener) 
        listener->panelChanged (EvROTDegree, Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), degree->getValue()));
}

void Rotate::straighten (double deg) {

    degree->setValue (degree->getValue()+deg);
    degree->setEditedState (Edited);
    if (listener) 
        listener->panelChanged (EvROTDegree, Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), degree->getValue()));
}

void Rotate::selectStraightPressed () {

    if (rlistener)
        rlistener->straightenRequested ();
}

void Rotate::autoCropPressed () {

    if (rlistener)
        rlistener->autoCropRequested ();
}

void Rotate::fillPressed () {

    if (batchMode) {
        if (fill->get_inconsistent()) {
            fill->set_inconsistent (false);
            fillConn.block (true);
            fill->set_active (false);
            fillConn.block (false);
        }
        else if (lastFill)
            fill->set_inconsistent (true);

        lastFill = fill->get_active ();
    }

    if (listener) {
        if (fill->get_active())
            listener->panelChanged (EvROTDegree, M("TP_ROTATE_FILL")+' '+M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvROTFill, M("TP_ROTATE_FILL")+' '+M("GENERAL_DISABLED"));
    }
}

void Rotate::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    degree->showEditedCB ();
    removeIfThere (this, autoCrop);
}

void Rotate::setAdjusterBehavior (bool brotadd) {

    if (!degAdd && brotadd || degAdd && !brotadd)
        degree->setLimits (-45, 45, 0.01, 0);
    
    degAdd = brotadd;
}
