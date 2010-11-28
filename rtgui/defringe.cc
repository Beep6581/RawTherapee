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
#include <defringe.h>
#include <iomanip>
#include <math.h>

using namespace rtengine;
using namespace rtengine::procparams;

Defringe::Defringe () : ToolPanel ()  {

  enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
  enabled->set_active (false);
  enabled->show ();
  pack_start (*enabled);

  Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
  hsep1->show ();
  pack_start (*hsep1);

  enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Defringe::enabledChanged) );
  //edgConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Defringe::edgeChanged) );

  radius  = Gtk::manage (new Adjuster (M("TP_DEFRINGE_RADIUS"), 0.5, 5.0, 0.1, 2.0));
  threshold    = Gtk::manage (new Adjuster (M("TP_DEFRINGE_THRESHOLD"), 0, 100, 1, 25));
  radius->setAdjusterListener (this);
  threshold->setAdjusterListener (this); 
  radius->show();
  threshold->show();

  pack_start (*radius);
  pack_start (*threshold);
}

void Defringe::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        radius->setEditedState    (pedited->defringe.radius ? Edited : UnEdited);
        threshold->setEditedState      (pedited->defringe.threshold ? Edited : UnEdited);
        enabled->set_inconsistent (!pedited->defringe.enabled);
    }

    enaConn.block (true);
    enabled->set_active (pp->defringe.enabled);
    enaConn.block (false);
    
    lastEnabled = pp->defringe.enabled;

    radius->setValue    (pp->defringe.radius);
    threshold->setValue      (pp->defringe.threshold);

    enableListener ();
}

void Defringe::write (ProcParams* pp, ParamsEdited* pedited) {

  pp->defringe.radius        = radius->getValue ();
  pp->defringe.threshold = (int)threshold->getValue ();
  pp->defringe.enabled       = enabled->get_active();

    if (pedited) {
        pedited->defringe.radius        = radius->getEditedState ();
        pedited->defringe.threshold = threshold->getEditedState ();
        pedited->defringe.enabled       = !enabled->get_inconsistent();
    }
}

void Defringe::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    radius->setDefault (defParams->defringe.radius);
    threshold->setDefault (defParams->defringe.threshold);

    if (pedited) {
        radius->setDefaultEditedState (pedited->defringe.radius ? Edited : UnEdited);
        threshold->setDefaultEditedState   (pedited->defringe.threshold ? Edited : UnEdited);
    }
    else {
        radius->setDefaultEditedState (Irrelevant);
        threshold->setDefaultEditedState   (Irrelevant);
    }
}

void Defringe::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {
        
        if (a==radius)
            listener->panelChanged (EvLDNRadius, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
        else if (a==threshold) 
            listener->panelChanged (EvLDNEdgeTolerance, Glib::ustring::format ((int)a->getValue()));
    }
}

void Defringe::enabledChanged () {

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
            listener->panelChanged (EvLDNEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvLDNEnabled, M("GENERAL_DISABLED"));
    }  
}

void Defringe::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    radius->showEditedCB ();
    threshold->showEditedCB ();
}

/*void Defringe::setAdjusterBehavior (bool bthresholdtoladd) {

    if (!thresholdtolAdd && bthresholdtoladd)
		threshold->setLimits (-10000, 10000, 100, 0);
	else if (thresholdtolAdd && !bthresholdtoladd)
		threshold->setLimits (100, 10000, 100, 1000);
    
    thresholdtolAdd = bthresholdtoladd;
}*/
