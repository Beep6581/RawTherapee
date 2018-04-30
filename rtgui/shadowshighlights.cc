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

ShadowsHighlights::ShadowsHighlights () : FoldableToolPanel(this, "shadowshighlights", M("TP_SHADOWSHLIGHTS_LABEL"), false, true)
{
    highlights   = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), 0, 100, 1, 0));
    h_tonalwidth = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_HLTONALW"), 10, 100, 1, 70));
    pack_start (*highlights);
    pack_start (*h_tonalwidth);

    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    shadows      = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_SHADOWS"), 0, 100, 1, 0));
    s_tonalwidth = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_SHTONALW"), 10, 100, 1, 30));
    pack_start (*shadows);
    pack_start (*s_tonalwidth);

    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    radius = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_RADIUS"), 1, 100, 1, 40));
    pack_start (*radius);

    radius->setAdjusterListener (this);
    highlights->setAdjusterListener (this);
    h_tonalwidth->setAdjusterListener (this);
    shadows->setAdjusterListener (this);
    s_tonalwidth->setAdjusterListener (this);

    show_all_children ();
}

void ShadowsHighlights::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        radius->setEditedState       (pedited->sh.radius ? Edited : UnEdited);
        highlights->setEditedState   (pedited->sh.highlights ? Edited : UnEdited);
        h_tonalwidth->setEditedState (pedited->sh.htonalwidth ? Edited : UnEdited);
        shadows->setEditedState      (pedited->sh.shadows ? Edited : UnEdited);
        s_tonalwidth->setEditedState (pedited->sh.stonalwidth ? Edited : UnEdited);
        set_inconsistent             (multiImage && !pedited->sh.enabled);
    }

    setEnabled (pp->sh.enabled);

    radius->setValue        (pp->sh.radius);
    highlights->setValue    (pp->sh.highlights);
    h_tonalwidth->setValue  (pp->sh.htonalwidth);
    shadows->setValue       (pp->sh.shadows);
    s_tonalwidth->setValue  (pp->sh.stonalwidth);

    enableListener ();
}

void ShadowsHighlights::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->sh.radius        = (int)radius->getValue ();
    pp->sh.highlights    = (int)highlights->getValue ();
    pp->sh.htonalwidth   = (int)h_tonalwidth->getValue ();
    pp->sh.shadows       = (int)shadows->getValue ();
    pp->sh.stonalwidth   = (int)s_tonalwidth->getValue ();
    pp->sh.enabled       = getEnabled();

    if (pedited) {
        pedited->sh.radius          = radius->getEditedState ();
        pedited->sh.highlights      = highlights->getEditedState ();
        pedited->sh.htonalwidth     = h_tonalwidth->getEditedState ();
        pedited->sh.shadows         = shadows->getEditedState ();
        pedited->sh.stonalwidth     = s_tonalwidth->getEditedState ();
        pedited->sh.enabled         = !get_inconsistent();
    }
}

void ShadowsHighlights::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    radius->setDefault (defParams->sh.radius);
    highlights->setDefault (defParams->sh.highlights);
    h_tonalwidth->setDefault (defParams->sh.htonalwidth);
    shadows->setDefault (defParams->sh.shadows);
    s_tonalwidth->setDefault (defParams->sh.stonalwidth);

    if (pedited) {
        radius->setDefaultEditedState       (pedited->sh.radius ? Edited : UnEdited);
        highlights->setDefaultEditedState   (pedited->sh.highlights ? Edited : UnEdited);
        h_tonalwidth->setDefaultEditedState (pedited->sh.htonalwidth ? Edited : UnEdited);
        shadows->setDefaultEditedState      (pedited->sh.shadows ? Edited : UnEdited);
        s_tonalwidth->setDefaultEditedState (pedited->sh.stonalwidth ? Edited : UnEdited);
    } else {
        radius->setDefaultEditedState       (Irrelevant);
        highlights->setDefaultEditedState   (Irrelevant);
        h_tonalwidth->setDefaultEditedState (Irrelevant);
        shadows->setDefaultEditedState      (Irrelevant);
        s_tonalwidth->setDefaultEditedState (Irrelevant);
    }
}

void ShadowsHighlights::adjusterChanged (Adjuster* a, double newval)
{

    if (listener && getEnabled()) {

        Glib::ustring costr = Glib::ustring::format ((int)a->getValue());

        if (a == highlights) {
            listener->panelChanged (EvSHHighlights, costr);
        } else if (a == h_tonalwidth) {
            listener->panelChanged (EvSHHLTonalW, costr);
        } else if (a == shadows) {
            listener->panelChanged (EvSHShadows, costr);
        } else if (a == s_tonalwidth) {
            listener->panelChanged (EvSHSHTonalW, costr);
        } else if (a == radius) {
            listener->panelChanged (EvSHRadius, costr);
        }
    }
}

void ShadowsHighlights::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvSHEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvSHEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvSHEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void ShadowsHighlights::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    radius->showEditedCB ();
    highlights->showEditedCB ();
    h_tonalwidth->showEditedCB ();
    shadows->showEditedCB ();
    s_tonalwidth->showEditedCB ();
}

void ShadowsHighlights::setAdjusterBehavior (bool hadd, bool sadd)
{

    highlights->setAddMode(hadd);
    shadows->setAddMode(sadd);
}

void ShadowsHighlights::trimValues (rtengine::procparams::ProcParams* pp)
{

    highlights->trimValue(pp->sh.highlights);
    shadows->trimValue(pp->sh.shadows);
}
