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
#include <impulsedenoise.h>
#include <math.h>
#include <iomanip>
#include <guiutils.h>

using namespace rtengine;
using namespace rtengine::procparams;

ImpulseDenoise::ImpulseDenoise () : ToolPanel() {

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);

	//amount = Gtk::manage (new Adjuster (M("TP_DETAIL_AMOUNT"), 1, 100, 1, 30));

	pack_start (*enabled);
	pack_start (*Gtk::manage (new  Gtk::HSeparator()));
	//pack_start (*amount);

	enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &ImpulseDenoise::enabledChanged) );
	//amount->setAdjusterListener (this);

	show_all_children ();
}

void ImpulseDenoise::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        //amount->setEditedState    (pedited->impulseDenoise.amount ? Edited : UnEdited);
        enabled->set_inconsistent (!pedited->impulseDenoise.enabled);
    }

    enaConn.block (true);
    enabled->set_active (pp->impulseDenoise.enabled);
    enaConn.block (false);
    
    lastEnabled = pp->impulseDenoise.enabled;

    //amount->setValue (pp->impulseDenoise.amount);

    enableListener ();
}

void ImpulseDenoise::write (ProcParams* pp, ParamsEdited* pedited) {

    //pp->impulseDenoise.amount    = amount->getValue ();
    pp->impulseDenoise.enabled   = enabled->get_active();
	
    if (pedited) {
        //pedited->impulseDenoise.amount        = amount->getEditedState ();
        pedited->impulseDenoise.enabled       = !enabled->get_inconsistent();
    }
}

void ImpulseDenoise::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    //amount->setDefault (defParams->impulseDenoise.amount);

    /*if (pedited) 
        amount->setDefaultEditedState (pedited->impulseDenoise.amount ? Edited : UnEdited);
    else
        amount->setDefaultEditedState (Irrelevant);*/
}

/*void ImpulseDenoise::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {

        listener->panelChanged (EvCDNRadius, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
    }
}*/

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
            listener->panelChanged (EvCDNEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvCDNEnabled, M("GENERAL_DISABLED"));
    }  
}

void ImpulseDenoise::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    //amount->showEditedCB ();
}
