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
#include "cacorrection.h"
#include <iomanip>
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

CACorrection::CACorrection () : FoldableToolPanel(this) {

    Gtk::Image* icaredL =   Gtk::manage (new RTImage ("ajd-ca-red1.png"));
    Gtk::Image* icaredR =   Gtk::manage (new RTImage ("ajd-ca-red2.png"));
    Gtk::Image* icablueL =  Gtk::manage (new RTImage ("ajd-ca-blue1.png"));
    Gtk::Image* icablueR =  Gtk::manage (new RTImage ("ajd-ca-blue2.png"));

    red = Gtk::manage (new Adjuster (M("TP_CACORRECTION_RED"), -0.005, 0.005, 0.0001, 0, icaredL, icaredR));
    red->setAdjusterListener (this); 

    blue = Gtk::manage (new Adjuster (M("TP_CACORRECTION_BLUE"), -0.005, 0.005, 0.0001, 0, icablueL, icablueR));
    blue->setAdjusterListener (this); 

    pack_start (*red);
    pack_start (*blue);

    show_all();
}

void CACorrection::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        red->setEditedState (pedited->cacorrection.red ? Edited : UnEdited);
        blue->setEditedState (pedited->cacorrection.blue ? Edited : UnEdited);
    }

    red->setValue (pp->cacorrection.red);
    blue->setValue (pp->cacorrection.blue);
    
    enableListener ();
}

void CACorrection::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->cacorrection.red  = red->getValue ();
    pp->cacorrection.blue = blue->getValue ();

    if (pedited) {
        pedited->cacorrection.red = red->getEditedState ();
        pedited->cacorrection.blue = blue->getEditedState ();
    }
}

void CACorrection::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    red->setDefault (defParams->cacorrection.red);
    blue->setDefault (defParams->cacorrection.blue);

    if (pedited) {
        red->setDefaultEditedState (pedited->cacorrection.red ? Edited : UnEdited);
        blue->setDefaultEditedState (pedited->cacorrection.blue ? Edited : UnEdited);
    }
    else {
        red->setDefaultEditedState (Irrelevant);
        blue->setDefaultEditedState (Irrelevant);
    }
}

void CACorrection::adjusterChanged (Adjuster* a, double newval) {

    if (listener) 
        listener->panelChanged (EvCACorr, Glib::ustring::compose ("%1=%3\n%2=%4", M("TP_CACORRECTION_RED"), M("TP_CACORRECTION_BLUE"), Glib::ustring::format (std::setw(5), std::fixed, std::setprecision(4), red->getValue()), Glib::ustring::format (std::setw(5), std::fixed, std::setprecision(4), blue->getValue())));
}

void CACorrection::setAdjusterBehavior (bool badd) {

	red->setAddMode(badd);
	blue->setAddMode(badd);
}

void CACorrection::trimValues (rtengine::procparams::ProcParams* pp) {

	red->trimValue(pp->cacorrection.red);
	blue->trimValue(pp->cacorrection.blue);
}

void CACorrection::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    red->showEditedCB ();
    blue->showEditedCB ();
}
