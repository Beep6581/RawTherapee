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
#include <lumadenoise.h>
#include <iomanip>
#include <math.h>

using namespace rtengine;
using namespace rtengine::procparams;

LumaDenoise::LumaDenoise () : ToolPanel ()  {

  enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
  enabled->set_active (false);
  enabled->show ();
  pack_start (*enabled);

  Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
  hsep1->show ();
  pack_start (*hsep1);

  enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &LumaDenoise::enabledChanged) );

  radius  = Gtk::manage (new Adjuster (M("TP_LUMADENOISE_RADIUS"), 0.5, 50, 0.1, 1.9));
  edge    = Gtk::manage (new Adjuster (M("TP_LUMADENOISE_EDGETOLERANCE"), 10, 30000, 100, 1500));
  radius->setAdjusterListener (this);
  edge->setAdjusterListener (this); 
  radius->show();
  edge->show();

  pack_start (*radius);
  pack_start (*edge);
}

void LumaDenoise::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        radius->setEditedState    (pedited->lumaDenoise.radius ? Edited : UnEdited);
        edge->setEditedState      (pedited->lumaDenoise.edgetolerance ? Edited : UnEdited);
        enabled->set_inconsistent (!pedited->lumaDenoise.enabled);
    }

    enaConn.block (true);
    enabled->set_active (pp->lumaDenoise.enabled);
    enaConn.block (false);
    
    lastEnabled = pp->lumaDenoise.enabled;

    radius->setValue    (pp->lumaDenoise.radius);
    edge->setValue      (pp->lumaDenoise.edgetolerance);

    enableListener ();
}

void LumaDenoise::write (ProcParams* pp, ParamsEdited* pedited) {

  pp->lumaDenoise.radius        = radius->getValue ();
  pp->lumaDenoise.edgetolerance = (int)edge->getValue ();
  pp->lumaDenoise.enabled       = enabled->get_active();

    if (pedited) {
        pedited->lumaDenoise.radius        = radius->getEditedState ();
        pedited->lumaDenoise.edgetolerance = edge->getEditedState ();
        pedited->lumaDenoise.enabled       = !enabled->get_inconsistent();
    }
}

void LumaDenoise::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    radius->setDefault (defParams->lumaDenoise.radius);
    edge->setDefault (defParams->lumaDenoise.edgetolerance);

    if (pedited) {
        radius->setDefaultEditedState (pedited->lumaDenoise.radius ? Edited : UnEdited);
        edge->setDefaultEditedState   (pedited->lumaDenoise.edgetolerance ? Edited : UnEdited);
    }
    else {
        radius->setDefaultEditedState (Irrelevant);
        edge->setDefaultEditedState   (Irrelevant);
    }
}

void LumaDenoise::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {
        
        if (a==radius)
            listener->panelChanged (EvLDNRadius, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
        else if (a==edge) 
            listener->panelChanged (EvLDNEdgeTolerance, Glib::ustring::format ((int)a->getValue()));
    }
}

void LumaDenoise::enabledChanged () {

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

void LumaDenoise::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    radius->showEditedCB ();
    edge->showEditedCB ();
}

void LumaDenoise::setAdjusterBehavior (bool bedgetoladd) {

    if (!edgetolAdd && bedgetoladd)
		edge->setLimits (-10000, 10000, 100, 0);
	else if (edgetolAdd && !bedgetoladd)
		edge->setLimits (10, 30000, 100, 1500);
    
    edgetolAdd = bedgetoladd;
}
