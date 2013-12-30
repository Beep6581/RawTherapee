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
#include "shadowshighlights.h"

using namespace rtengine;
using namespace rtengine::procparams;

ShadowsHighlights::ShadowsHighlights () : Gtk::VBox(), FoldableToolPanel(this) {

  set_border_width(4);

  enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
  enabled->set_active (false);
  pack_start (*enabled);
  enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &ShadowsHighlights::enabledChanged) );

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  hq = Gtk::manage (new Gtk::CheckButton (M("TP_SHADOWSHLIGHTS_SHARPMASK")));
  hq->set_active (false);
  pack_start (*hq);
  hqConn = hq->signal_toggled().connect( sigc::mem_fun(*this, &ShadowsHighlights::hqChanged) );

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));
  highlights   = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), 0, 100, 1, 0));
  h_tonalwidth = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_HLTONALW"), 10, 100, 1, 80));
  pack_start (*highlights);
  pack_start (*h_tonalwidth);

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  shadows      = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_SHADOWS"), 0, 100, 1, 0));
  s_tonalwidth = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_SHTONALW"), 10, 100, 1, 80));
  pack_start (*shadows);
  pack_start (*s_tonalwidth);

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  lcontrast = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_LOCALCONTR"), 0, 100, 1, 0));
  pack_start (*lcontrast);

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  radius = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_RADIUS"), 5, 100, 1, 30));
  pack_start (*radius);

  radius->setAdjusterListener (this);
  highlights->setAdjusterListener (this); 
  h_tonalwidth->setAdjusterListener (this); 
  shadows->setAdjusterListener (this); 
  s_tonalwidth->setAdjusterListener (this); 
  lcontrast->setAdjusterListener (this); 
  
  show_all_children ();
}

void ShadowsHighlights::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        radius->setEditedState      (pedited->sh.radius ? Edited : UnEdited);
        lcontrast->setEditedState   (pedited->sh.localcontrast ? Edited : UnEdited);
        highlights->setEditedState  (pedited->sh.highlights ? Edited : UnEdited);
        h_tonalwidth->setEditedState (pedited->sh.htonalwidth ? Edited : UnEdited);
        shadows->setEditedState     (pedited->sh.shadows ? Edited : UnEdited);
        s_tonalwidth->setEditedState (pedited->sh.stonalwidth ? Edited : UnEdited);
        enabled->set_inconsistent   (!pedited->sh.enabled);
        hq->set_inconsistent        (!pedited->sh.hq);
    }

    enaConn.block (true);
    enabled->set_active (pp->sh.enabled);
    enaConn.block (false);
    hqConn.block (true);
    hq->set_active (pp->sh.hq);
    hqConn.block (false);
    
    lastEnabled = pp->sh.enabled;
    lastHQ = pp->sh.hq;

    radius->setValue        (pp->sh.radius);
    lcontrast->setValue     (pp->sh.localcontrast);
    highlights->setValue    (pp->sh.highlights);
    h_tonalwidth->setValue  (pp->sh.htonalwidth);
    shadows->setValue       (pp->sh.shadows);
    s_tonalwidth->setValue  (pp->sh.stonalwidth);

    enableListener ();
}

void ShadowsHighlights::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->sh.radius        = (int)radius->getValue ();
    pp->sh.localcontrast = (int)lcontrast->getValue ();
    pp->sh.highlights    = (int)highlights->getValue ();
    pp->sh.htonalwidth   = (int)h_tonalwidth->getValue ();
    pp->sh.shadows       = (int)shadows->getValue ();
    pp->sh.stonalwidth   = (int)s_tonalwidth->getValue ();
    pp->sh.enabled       = enabled->get_active();
    pp->sh.hq            = hq->get_active();

    if (pedited) {
        pedited->sh.radius          = radius->getEditedState ();
        pedited->sh.localcontrast   = lcontrast->getEditedState ();
        pedited->sh.highlights      = highlights->getEditedState ();
        pedited->sh.htonalwidth     = h_tonalwidth->getEditedState ();
        pedited->sh.shadows         = shadows->getEditedState ();
        pedited->sh.stonalwidth     = s_tonalwidth->getEditedState ();
        pedited->sh.enabled         = !enabled->get_inconsistent();
        pedited->sh.hq              = !hq->get_inconsistent();
    }
}

void ShadowsHighlights::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    radius->setDefault (defParams->sh.radius);
    lcontrast->setDefault (defParams->sh.localcontrast);
    highlights->setDefault (defParams->sh.highlights);
    h_tonalwidth->setDefault (defParams->sh.htonalwidth);
    shadows->setDefault (defParams->sh.shadows);
    s_tonalwidth->setDefault (defParams->sh.stonalwidth);

    if (pedited) {
        radius->setDefaultEditedState       (pedited->sh.radius ? Edited : UnEdited);
        lcontrast->setDefaultEditedState    (pedited->sh.localcontrast ? Edited : UnEdited);
        highlights->setDefaultEditedState   (pedited->sh.highlights ? Edited : UnEdited);
        h_tonalwidth->setDefaultEditedState (pedited->sh.htonalwidth ? Edited : UnEdited);
        shadows->setDefaultEditedState      (pedited->sh.shadows ? Edited : UnEdited);
        s_tonalwidth->setDefaultEditedState (pedited->sh.stonalwidth ? Edited : UnEdited);
    }
    else {
        radius->setDefaultEditedState       (Irrelevant);
        lcontrast->setDefaultEditedState    (Irrelevant);
        highlights->setDefaultEditedState   (Irrelevant);
        h_tonalwidth->setDefaultEditedState (Irrelevant);
        shadows->setDefaultEditedState      (Irrelevant);
        s_tonalwidth->setDefaultEditedState (Irrelevant);
    }
}

void ShadowsHighlights::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {

        Glib::ustring costr = Glib::ustring::format ((int)a->getValue());

        if (a==highlights) 
            listener->panelChanged (EvSHHighlights, costr);
        else if (a==h_tonalwidth) 
            listener->panelChanged (EvSHHLTonalW, costr);
        else if (a==shadows) 
            listener->panelChanged (EvSHShadows, costr);
        else if (a==s_tonalwidth) 
            listener->panelChanged (EvSHSHTonalW, costr);
        else if (a==radius) 
            listener->panelChanged (EvSHRadius, costr);
        else if (a==lcontrast) 
            listener->panelChanged (EvSHLContrast, costr);
    }
}

void ShadowsHighlights::enabledChanged () {

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
        if (enabled->get_active())
            listener->panelChanged (EvSHEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvSHEnabled, M("GENERAL_DISABLED"));
    }
}

void ShadowsHighlights::hqChanged () {

    if (batchMode) {
        if (hq->get_inconsistent()) {
            hq->set_inconsistent (false);
            hqConn.block (true);
            hq->set_active (false);
            hqConn.block (false);
        }
        else if (lastHQ)
            hq->set_inconsistent (true);

        lastHQ = hq->get_active ();
    }
    
    if (listener) {
        if (hq->get_active())
            listener->panelChanged (EvSHHighQuality, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvSHHighQuality, M("GENERAL_DISABLED"));
    }
}

void ShadowsHighlights::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    radius->showEditedCB ();
    lcontrast->showEditedCB ();
    highlights->showEditedCB ();
    h_tonalwidth->showEditedCB ();
    shadows->showEditedCB ();
    s_tonalwidth->showEditedCB ();
}

void ShadowsHighlights::setAdjusterBehavior (bool hadd, bool sadd, bool lcadd) {

	highlights->setAddMode(hadd);
	shadows->setAddMode(sadd);
	lcontrast->setAddMode(lcadd);
}

void ShadowsHighlights::trimValues (rtengine::procparams::ProcParams* pp) {

	highlights->trimValue(pp->sh.highlights);
	shadows->trimValue(pp->sh.shadows);
	lcontrast->trimValue(pp->sh.localcontrast);
}
