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
#include <colordenoise.h>
#include <math.h>
#include <iomanip>
#include <guiutils.h>

using namespace rtengine;
using namespace rtengine::procparams;

ColorDenoise::ColorDenoise () : ToolPanel() {

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);

	amount = Gtk::manage (new Adjuster ("Amount", 1, 100, 1, 30));

	pack_start (*enabled);
	pack_start (*Gtk::manage (new  Gtk::HSeparator()));
	pack_start (*amount);

	enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &ColorDenoise::enabledChanged) );
	amount->setAdjusterListener (this);

	show_all_children ();
}

void ColorDenoise::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        amount->setEditedState    (pedited->colorDenoise.amount ? Edited : UnEdited);
        enabled->set_inconsistent (!pedited->colorDenoise.enabled);
    }

    enaConn.block (true);
    enabled->set_active (pp->colorDenoise.enabled);
    enaConn.block (false);
    
    lastEnabled = pp->colorDenoise.enabled;

    amount->setValue (pp->colorDenoise.amount);

    enableListener ();
}

void ColorDenoise::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->colorDenoise.amount    = amount->getValue ();
    pp->colorDenoise.enabled   = enabled->get_active();
	
    if (pedited) {
        pedited->colorDenoise.amount        = amount->getEditedState ();
        pedited->colorDenoise.enabled       = !enabled->get_inconsistent();
    }
}

void ColorDenoise::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    amount->setDefault (defParams->colorDenoise.amount);

    if (pedited) 
        amount->setDefaultEditedState (pedited->colorDenoise.amount ? Edited : UnEdited);
    else
        amount->setDefaultEditedState (Irrelevant);
}

void ColorDenoise::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {

        listener->panelChanged (EvCDNRadius, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
    }
}

void ColorDenoise::enabledChanged () {

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

void ColorDenoise::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    amount->showEditedCB ();
}
