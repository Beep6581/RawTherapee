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
#include "xtransrawexposure.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

XTransRAWExposure::XTransRAWExposure () : FoldableToolPanel(this, "xtransrawexposure", M("TP_EXPOS_BLACKPOINT_LABEL"))
{
    PexBlackRed = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACK_RED"), -2048, 2048, 0.1, 0)); //black level
    PexBlackRed->setAdjusterListener (this);

    if (PexBlackRed->delay < options.adjusterMaxDelay) {
        PexBlackRed->delay = options.adjusterMaxDelay;
    }

    PexBlackRed->show();
    PexBlackGreen = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACK_GREEN"), -2048, 2048, 0.1, 0)); //black level
    PexBlackGreen->setAdjusterListener (this);

    if (PexBlackGreen->delay < options.adjusterMaxDelay) {
        PexBlackGreen->delay = options.adjusterMaxDelay;
    }

    PexBlackGreen->show();
    PexBlackBlue = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACK_BLUE"), -2048, 2048, 0.1, 0)); //black level
    PexBlackBlue->setAdjusterListener (this);

    if (PexBlackBlue->delay < options.adjusterMaxDelay) {
        PexBlackBlue->delay = options.adjusterMaxDelay;
    }

    PexBlackBlue->show();

    pack_start( *PexBlackRed, Gtk::PACK_SHRINK, 0);//black
    pack_start( *PexBlackGreen, Gtk::PACK_SHRINK, 0);//black
    pack_start( *PexBlackBlue, Gtk::PACK_SHRINK, 0);//black
}

void XTransRAWExposure::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited ) {
        PexBlackRed->setEditedState( pedited->raw.xtranssensor.exBlackRed ? Edited : UnEdited );
        PexBlackGreen->setEditedState( pedited->raw.xtranssensor.exBlackGreen ? Edited : UnEdited );
        PexBlackBlue->setEditedState( pedited->raw.xtranssensor.exBlackBlue ? Edited : UnEdited );
    }

    PexBlackRed->setValue (pp->raw.xtranssensor.blackred);//black
    PexBlackGreen->setValue (pp->raw.xtranssensor.blackgreen);//black
    PexBlackBlue->setValue (pp->raw.xtranssensor.blackblue);//black

    enableListener ();
}

void XTransRAWExposure::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.xtranssensor.blackred   = PexBlackRed->getValue();// black
    pp->raw.xtranssensor.blackgreen = PexBlackGreen->getValue();// black
    pp->raw.xtranssensor.blackblue  = PexBlackBlue->getValue();// black

    if (pedited) {
        pedited->raw.xtranssensor.exBlackRed = PexBlackRed->getEditedState ();//black
        pedited->raw.xtranssensor.exBlackGreen = PexBlackGreen->getEditedState ();//black
        pedited->raw.xtranssensor.exBlackBlue = PexBlackBlue->getEditedState ();//black
    }

}

void XTransRAWExposure::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {
        Glib::ustring value = a->getTextValue();

        if      (a == PexBlackRed) {
            listener->panelChanged (EvPreProcessExpBlackRed, value);
        } else if (a == PexBlackGreen) {
            listener->panelChanged (EvPreProcessExpBlackGreen, value);
        } else if (a == PexBlackBlue) {
            listener->panelChanged (EvPreProcessExpBlackBlue, value);
        }
    }
}

void XTransRAWExposure::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    PexBlackRed->showEditedCB ();//black
    PexBlackGreen->showEditedCB ();//black
    PexBlackBlue->showEditedCB ();//black

}

void XTransRAWExposure::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    PexBlackRed->setDefault( defParams->raw.xtranssensor.blackred);
    PexBlackGreen->setDefault( defParams->raw.xtranssensor.blackgreen);
    PexBlackBlue->setDefault( defParams->raw.xtranssensor.blackblue);

    if (pedited) {
        PexBlackRed->setDefaultEditedState( pedited->raw.xtranssensor.exBlackRed ? Edited : UnEdited);
        PexBlackGreen->setDefaultEditedState( pedited->raw.xtranssensor.exBlackGreen ? Edited : UnEdited);
        PexBlackBlue->setDefaultEditedState( pedited->raw.xtranssensor.exBlackBlue ? Edited : UnEdited);

    } else {
        PexBlackRed->setDefaultEditedState( Irrelevant );
        PexBlackGreen->setDefaultEditedState( Irrelevant );
        PexBlackBlue->setDefaultEditedState( Irrelevant );

    }
}

void XTransRAWExposure::setAdjusterBehavior (bool pexblackadd)
{

    PexBlackRed->setAddMode(pexblackadd);
    PexBlackGreen->setAddMode(pexblackadd);
    PexBlackBlue->setAddMode(pexblackadd);
}

void XTransRAWExposure::trimValues (rtengine::procparams::ProcParams* pp)
{

    PexBlackRed->trimValue(pp->raw.xtranssensor.blackred);
    PexBlackGreen->trimValue(pp->raw.xtranssensor.blackgreen);
    PexBlackBlue->trimValue(pp->raw.xtranssensor.blackblue);
}
