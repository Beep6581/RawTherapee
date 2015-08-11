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
#include "rawexposure.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

RAWExposure::RAWExposure () : FoldableToolPanel(this, "rawexposure", M("TP_EXPOS_WHITEPOINT_LABEL"))
{
    PexPos = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_LINEAR"), 0.1, 16.0, 0.01, 1));
    PexPos->setAdjusterListener (this);

    if (PexPos->delay < 1000) {
        PexPos->delay = 1000;
    }

    PexPos->show();
    PexPreser = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_PRESER"), 0, 2.5, 0.1, 0));
    PexPreser->setAdjusterListener (this);

    if (PexPreser->delay < 1000) {
        PexPreser->delay = 1000;
    }

    PexPreser->show();

    pack_start( *PexPos, Gtk::PACK_SHRINK, 4);//exposi
    // raw highlight exposure setting is obsolete, removing from GUI
    //pack_start( *PexPreser, Gtk::PACK_SHRINK, 4);
}

void RAWExposure::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited ) {
        PexPos->setEditedState( pedited->raw.exPos ? Edited : UnEdited );
        PexPreser->setEditedState( pedited->raw.exPreser ? Edited : UnEdited );
    }

    PexPos->setValue (pp->raw.expos);
    PexPreser->setValue (pp->raw.preser);//exposi

    enableListener ();
}

void RAWExposure::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.expos = PexPos->getValue();
    pp->raw.preser = PexPreser->getValue();//exposi

    if (pedited) {
        pedited->raw.exPos = PexPos->getEditedState ();
        pedited->raw.exPreser = PexPreser->getEditedState ();//exposi
    }

}

void RAWExposure::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {
        Glib::ustring value = a->getTextValue();

        if (a == PexPos ) {
            listener->panelChanged (EvPreProcessExpCorrLinear,  value );
        } else if (a == PexPreser && ABS(PexPos->getValue() - 1.0) > 0.0001) { // update takes long, only do it if it would have an effect
            listener->panelChanged (EvPreProcessExpCorrPH,  value );
        }
    }
}

void RAWExposure::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    PexPos->showEditedCB ();
    PexPreser->showEditedCB ();//exposure
}

void RAWExposure::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    PexPos->setDefault( defParams->raw.expos);
    PexPreser->setDefault( defParams->raw.preser);

    if (pedited) {
        PexPos->setDefaultEditedState( pedited->raw.exPos ? Edited : UnEdited);
        PexPreser->setDefaultEditedState( pedited->raw.exPreser ? Edited : UnEdited);
    } else {
        PexPos->setDefaultEditedState( Irrelevant );
        PexPreser->setDefaultEditedState( Irrelevant );
    }
}

void RAWExposure::setAdjusterBehavior (bool pexposadd, bool pexpreseradd)
{

    PexPos->setAddMode(pexposadd);
    PexPreser->setAddMode(pexpreseradd);
}

void RAWExposure::trimValues (rtengine::procparams::ProcParams* pp)
{

    PexPos->trimValue(pp->raw.expos);
    PexPreser->trimValue(pp->raw.preser);
}
