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
#include "rawcacorrection.h"
#include "guiutils.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

RAWCACorr::RAWCACorr () : FoldableToolPanel(this, "rawcacorrection", M("TP_CHROMATABERR_LABEL"))
{
    Gtk::Image* icaredL =   Gtk::manage (new RTImage ("circle-red-cyan-small.png"));
    Gtk::Image* icaredR =   Gtk::manage (new RTImage ("circle-cyan-red-small.png"));
    Gtk::Image* icablueL =  Gtk::manage (new RTImage ("circle-blue-yellow-small.png"));
    Gtk::Image* icablueR =  Gtk::manage (new RTImage ("circle-yellow-blue-small.png"));

    caAutocorrect = Gtk::manage (new CheckBox(M("TP_RAWCACORR_AUTO"), multiImage));
    caAutocorrect->setCheckBoxListener (this);

    caRed = Gtk::manage(new Adjuster (M("TP_RAWCACORR_CARED"), -4.0, 4.0, 0.1, 0, icaredL, icaredR));
    caRed->setAdjusterListener (this);

    if (caRed->delay < options.adjusterMaxDelay) {
        caRed->delay = options.adjusterMaxDelay;
    }

    caRed->show();
    caBlue = Gtk::manage(new Adjuster (M("TP_RAWCACORR_CABLUE"), -4.0, 4.0, 0.1, 0, icablueL, icablueR));
    caBlue->setAdjusterListener (this);

    if (caBlue->delay < options.adjusterMaxDelay) {
        caBlue->delay = options.adjusterMaxDelay;
    }

    caBlue->show();

    pack_start( *caAutocorrect, Gtk::PACK_SHRINK, 4);
    pack_start( *caRed, Gtk::PACK_SHRINK, 4);
    pack_start( *caBlue, Gtk::PACK_SHRINK, 4);

}

void RAWCACorr::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited ) {
        caAutocorrect->setEdited(pedited->raw.ca_autocorrect);
        caRed->setEditedState( pedited->raw.cared ? Edited : UnEdited );
        caBlue->setEditedState( pedited->raw.cablue ? Edited : UnEdited );
    }

    // disable Red and Blue sliders when caAutocorrect is enabled
    caRed->set_sensitive(!pp->raw.ca_autocorrect);
    caBlue->set_sensitive(!pp->raw.ca_autocorrect);

    caAutocorrect->setValue(pp->raw.ca_autocorrect);
    caRed->setValue (pp->raw.cared);
    caBlue->setValue (pp->raw.cablue);

    enableListener ();
}

void RAWCACorr::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.ca_autocorrect = caAutocorrect->getLastActive();
    pp->raw.cared = caRed->getValue();
    pp->raw.cablue = caBlue->getValue();

    if (pedited) {
        pedited->raw.ca_autocorrect = !caAutocorrect->get_inconsistent();
        pedited->raw.cared = caRed->getEditedState ();
        pedited->raw.cablue = caBlue->getEditedState ();
    }

}

void RAWCACorr::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {

        Glib::ustring value = a->getTextValue();

        if (a == caRed) {
            listener->panelChanged (EvPreProcessCARed,  value );
        } else if (a == caBlue) {
            listener->panelChanged (EvPreProcessCABlue,  value );
        }
    }
}

void RAWCACorr::checkBoxToggled (CheckBox* c, CheckValue newval)
{
    if (c == caAutocorrect) {
        if (!batchMode) {
            // disable Red and Blue sliders when caAutocorrect is enabled
            caRed->set_sensitive(!caAutocorrect->getLastActive ());
            caBlue->set_sensitive(!caAutocorrect->getLastActive ());
        }
        if (listener) {
            listener->panelChanged (EvPreProcessAutoCA, caAutocorrect->getLastActive() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
        }
    }
}

void RAWCACorr::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    caRed->showEditedCB ();
    caBlue->showEditedCB ();
}

void RAWCACorr::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    caRed->setDefault( defParams->raw.cared);
    caBlue->setDefault( defParams->raw.cablue);

    if (pedited) {
        caRed->setDefaultEditedState( pedited->raw.cared ? Edited : UnEdited);
        caBlue->setDefaultEditedState( pedited->raw.cablue ? Edited : UnEdited);
    } else {
        caRed->setDefaultEditedState( Irrelevant );
        caBlue->setDefaultEditedState( Irrelevant );
    }
}

void RAWCACorr::setAdjusterBehavior (bool caadd)
{

    caRed->setAddMode(caadd);
    caBlue->setAddMode(caadd);
}

void RAWCACorr::trimValues (rtengine::procparams::ProcParams* pp)
{

    caRed->trimValue(pp->raw.cared);
    caBlue->trimValue(pp->raw.cablue);
}
