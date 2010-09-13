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
#include <dirpyrdenoise.h>
#include <iomanip>
#include <math.h>

using namespace rtengine;
using namespace rtengine::procparams;

DirPyrDenoise::DirPyrDenoise () : ToolPanel ()  {

  enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
  enabled->set_active (false);
  enabled->show ();
  pack_start (*enabled);

  Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
  hsep1->show ();
  pack_start (*hsep1);

  enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::enabledChanged) );

  luma  = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_LUMA"), 0, 100, 1, 10));
  chroma    = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_CHROMA"), 0, 100, 1, 10));
  luma->setAdjusterListener (this);
  chroma->setAdjusterListener (this); 
  luma->show();
  chroma->show();

  pack_start (*luma);
  pack_start (*chroma);
}

void DirPyrDenoise::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        luma->setEditedState    (pedited->dirpyrDenoise.luma ? Edited : UnEdited);
        chroma->setEditedState      (pedited->dirpyrDenoise.chroma ? Edited : UnEdited);
        enabled->set_inconsistent (!pedited->dirpyrDenoise.enabled);
    }

    enaConn.block (true);
    enabled->set_active (pp->dirpyrDenoise.enabled);
    enaConn.block (false);
    
    lastEnabled = pp->dirpyrDenoise.enabled;

    luma->setValue    (pp->dirpyrDenoise.luma);
    chroma->setValue      (pp->dirpyrDenoise.chroma);

    enableListener ();
}

void DirPyrDenoise::write (ProcParams* pp, ParamsEdited* pedited) {

  pp->dirpyrDenoise.luma        = luma->getValue ();
  pp->dirpyrDenoise.chroma = (int)chroma->getValue ();
  pp->dirpyrDenoise.enabled       = enabled->get_active();

    if (pedited) {
        pedited->dirpyrDenoise.luma     = luma->getEditedState ();
        pedited->dirpyrDenoise.chroma	= chroma->getEditedState ();
        pedited->dirpyrDenoise.enabled  = !enabled->get_inconsistent();
    }
}

void DirPyrDenoise::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    luma->setDefault (defParams->dirpyrDenoise.luma);
    chroma->setDefault (defParams->dirpyrDenoise.chroma);

    if (pedited) {
        luma->setDefaultEditedState (pedited->dirpyrDenoise.luma ? Edited : UnEdited);
        chroma->setDefaultEditedState   (pedited->dirpyrDenoise.chroma ? Edited : UnEdited);
    }
    else {
        luma->setDefaultEditedState (Irrelevant);
        chroma->setDefaultEditedState   (Irrelevant);
    }
}

void DirPyrDenoise::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {
        
        if (a==luma)
            listener->panelChanged (EvLDNRadius, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
        else if (a==chroma) 
            listener->panelChanged (EvLDNEdgeTolerance, Glib::ustring::format ((int)a->getValue()));
    }
}

void DirPyrDenoise::enabledChanged () {

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

void DirPyrDenoise::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    luma->showEditedCB ();
    chroma->showEditedCB ();
}

/*void DirPyrDenoise::setAdjusterBehavior (bool bedgetoladd) {

    if (!edgetolAdd && bedgetoladd)
		edge->setLimits (-10000, 10000, 100, 0);
	else if (edgetolAdd && !bedgetoladd)
		edge->setLimits (10, 30000, 100, 1500);
    
    edgetolAdd = bedgetoladd;
}*/
