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
#include <colorboost.h>
#include <iomanip>
#include <guiutils.h>

using namespace rtengine;
using namespace rtengine::procparams;

ColorBoost::ColorBoost () : ToolPanel(), cbAdd(false) {

  colorboost = new Adjuster (M("TP_COLORBOOST_AMOUNT"), -100, 300, 1, 0);

  pack_start (*colorboost);
  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  avoidclip = Gtk::manage (new Gtk::CheckButton (M("TP_COLORBOOST_AVOIDCOLORCLIP")));

  pack_start (*avoidclip);
  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  enablelimiter = Gtk::manage (new Gtk::CheckButton (M("TP_COLORBOOST_ENABLESATLIMITER")));
  pack_start (*enablelimiter);

  saturationlimiter = new Adjuster (M("TP_COLORBOOST_SATLIMIT"), 0, 200, 0.1, 100);
  saturationlimiter->show ();
  saturationlimiter->reference ();  

  colorboost->setAdjusterListener (this);
  saturationlimiter->setAdjusterListener (this);
  acconn = avoidclip->signal_toggled().connect( sigc::mem_fun(*this, &ColorBoost::avoidclip_toggled) );
  elconn = enablelimiter->signal_toggled().connect( sigc::mem_fun(*this, &ColorBoost::enablelimiter_toggled) );
  
  show_all_children ();
}

ColorBoost::~ColorBoost () {

    delete saturationlimiter;
}

void ColorBoost::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        colorboost->setEditedState (pedited->colorBoost.amount ? Edited : UnEdited);
        saturationlimiter->setEditedState (pedited->colorBoost.saturationlimit ? Edited : UnEdited);
        avoidclip->set_inconsistent (!pedited->colorBoost.avoidclip);
        enablelimiter->set_inconsistent (!pedited->colorBoost.enable_saturationlimiter);
    }

    colorboost->setValue (pp->colorBoost.amount);
    saturationlimiter->setValue (pp->colorBoost.saturationlimit);
    acconn.block (true);
    avoidclip->set_active (pp->colorBoost.avoidclip);
    acconn.block (false);
    elconn.block (true);
    enablelimiter->set_active (pp->colorBoost.enable_saturationlimiter);
    elconn.block (false);

    removeIfThere (this, saturationlimiter, false);
    if (enablelimiter->get_active () || enablelimiter->get_inconsistent())
        pack_start (*saturationlimiter);

    lastACVal = pp->colorBoost.avoidclip;
    lastELVal = pp->colorBoost.enable_saturationlimiter;
    
    enableListener ();
}

void ColorBoost::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->colorBoost.amount = (int)colorboost->getValue();
    pp->colorBoost.avoidclip = avoidclip->get_active ();
    pp->colorBoost.enable_saturationlimiter = enablelimiter->get_active ();
    pp->colorBoost.saturationlimit = saturationlimiter->getValue ();

    if (pedited) {
        pedited->colorBoost.amount = colorboost->getEditedState ();
        pedited->colorBoost.avoidclip = !avoidclip->get_inconsistent();
        pedited->colorBoost.enable_saturationlimiter = !enablelimiter->get_inconsistent();
        pedited->colorBoost.saturationlimit = saturationlimiter->getEditedState ();
    }
}

void ColorBoost::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    colorboost->setDefault (defParams->colorBoost.amount);
    saturationlimiter->setDefault (defParams->colorBoost.saturationlimit);

    if (pedited) {
        colorboost->setDefaultEditedState (pedited->colorBoost.amount ? Edited : UnEdited);
        saturationlimiter->setDefaultEditedState (pedited->colorBoost.saturationlimit ? Edited : UnEdited);
    }
    else {
        colorboost->setDefaultEditedState (Irrelevant);
        saturationlimiter->setDefaultEditedState (Irrelevant);
    }
}

void ColorBoost::avoidclip_toggled () {

    if (batchMode) {
        if (avoidclip->get_inconsistent()) {
            avoidclip->set_inconsistent (false);
            acconn.block (true);
            avoidclip->set_active (false);
            acconn.block (false);
        }
        else if (lastACVal)
            avoidclip->set_inconsistent (true);

        lastACVal = avoidclip->get_active ();
    }

    if (listener) {
        if (avoidclip->get_active ()) 
            listener->panelChanged (EvCBAvoidClip, M("GENERAL_ENABLED"));
        else            
            listener->panelChanged (EvCBAvoidClip, M("GENERAL_DISABLED"));
    }
}

void ColorBoost::enablelimiter_toggled () {

    if (batchMode) {
        if (enablelimiter->get_inconsistent()) {
            enablelimiter->set_inconsistent (false);
            elconn.block (true);
            enablelimiter->set_active (false);
            elconn.block (false);
        }
        else if (lastELVal)
            enablelimiter->set_inconsistent (true);

        lastELVal = enablelimiter->get_active ();
    }

    removeIfThere (this, saturationlimiter, false);
    if (enablelimiter->get_active () || enablelimiter->get_inconsistent())
        pack_start (*saturationlimiter);

    if (listener) {
        if (enablelimiter->get_active ()) 
            listener->panelChanged (EvCBSatLimiter, M("GENERAL_ENABLED"));
        else            
            listener->panelChanged (EvCBSatLimiter, M("GENERAL_DISABLED"));
    }
}

void ColorBoost::adjusterChanged (Adjuster* a, double newval) {

    if (listener) {
        if (a!=saturationlimiter) 
            listener->panelChanged (EvCBBoost, Glib::ustring::format ((int)a->getValue()));
        else
            listener->panelChanged (EvCBSatLimit, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
    }
}

void ColorBoost::setAdjusterBehavior (bool bcbadd) {

    if (!cbAdd && bcbadd || cbAdd && !bcbadd)
        colorboost->setLimits (-100, 100, 1, 0);
    
    cbAdd = bcbadd;
}

void ColorBoost::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    colorboost->showEditedCB ();
    saturationlimiter->showEditedCB ();
}
