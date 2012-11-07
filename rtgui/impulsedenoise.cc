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
#include "impulsedenoise.h"
#include <cmath>
#include <iomanip>
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;

ImpulseDenoise::ImpulseDenoise () : Gtk::VBox(), FoldableToolPanel(this) {

	set_border_width(4);

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);

	thresh = Gtk::manage (new Adjuster (M("TP_IMPULSEDENOISE_THRESH"), 0, 100, 1, 50));

	pack_start (*enabled);
	pack_start (*Gtk::manage (new  Gtk::HSeparator()));
	pack_start (*thresh);

	enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &ImpulseDenoise::enabledChanged) );
	thresh->setAdjusterListener (this);

	show_all_children ();
}

void ImpulseDenoise::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        thresh->setEditedState    (pedited->impulseDenoise.thresh ? Edited : UnEdited);
        enabled->set_inconsistent (!pedited->impulseDenoise.enabled);
    }

    enaConn.block (true);
    enabled->set_active (pp->impulseDenoise.enabled);
    enaConn.block (false);
    
    lastEnabled = pp->impulseDenoise.enabled;

    thresh->setValue (pp->impulseDenoise.thresh);

    enableListener ();
}

void ImpulseDenoise::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->impulseDenoise.thresh    = thresh->getValue ();
    pp->impulseDenoise.enabled   = enabled->get_active();
	
    if (pedited) {
        pedited->impulseDenoise.thresh        = thresh->getEditedState ();
        pedited->impulseDenoise.enabled       = !enabled->get_inconsistent();
    }
}

void ImpulseDenoise::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    thresh->setDefault (defParams->impulseDenoise.thresh);

    if (pedited) 
        thresh->setDefaultEditedState (pedited->impulseDenoise.thresh ? Edited : UnEdited);
    else
        thresh->setDefaultEditedState (Irrelevant);
}

void ImpulseDenoise::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {

        listener->panelChanged (EvIDNThresh, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
    }
}

void ImpulseDenoise::enabledChanged () {

    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            enaConn.block (true);
            enabled->set_active (false);
            enaConn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);

        lastEnabled = enabled->get_active ();
    }

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvIDNEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvIDNEnabled, M("GENERAL_DISABLED"));
    }  
}

void ImpulseDenoise::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    thresh->showEditedCB ();
}

void ImpulseDenoise::setAdjusterBehavior (bool threshadd) {

	thresh->setAddMode(threshadd);
}

void ImpulseDenoise::trimValues (rtengine::procparams::ProcParams* pp) {

	thresh->trimValue(pp->impulseDenoise.thresh);
}
