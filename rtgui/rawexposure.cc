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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <sstream>

#include "rawexposure.h"

#include "guiutils.h"
#include "options.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

RAWExposure::RAWExposure () : FoldableToolPanel(this, "rawexposure", M("TP_EXPOS_WHITEPOINT_LABEL"))
{
    PexPos = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_LINEAR"), 0.1, 16.0, 0.01, 1));
    PexPos->setAdjusterListener (this);

    PexPos->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    PexPos->show();
    pack_start( *PexPos, Gtk::PACK_SHRINK, 4);//exposi
    PexPos->setLogScale(100, 0);
}

void RAWExposure::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited ) {
        PexPos->setEditedState( pedited->raw.exPos ? Edited : UnEdited );
    }

    PexPos->setValue (pp->raw.expos);

    enableListener ();
}

void RAWExposure::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.expos = PexPos->getValue();

    if (pedited) {
        pedited->raw.exPos = PexPos->getEditedState ();
    }

}

void RAWExposure::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        Glib::ustring value = a->getTextValue();

        if (a == PexPos ) {
            listener->panelChanged (EvPreProcessExpCorrLinear,  value );
        }
    }
}

void RAWExposure::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    PexPos->showEditedCB ();
}

void RAWExposure::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    PexPos->setDefault( defParams->raw.expos);

    if (pedited) {
        PexPos->setDefaultEditedState( pedited->raw.exPos ? Edited : UnEdited);
    } else {
        PexPos->setDefaultEditedState( Irrelevant );
    }
}

void RAWExposure::setAdjusterBehavior (bool pexposadd)
{

    PexPos->setAddMode(pexposadd);
}

void RAWExposure::trimValues (rtengine::procparams::ProcParams* pp)
{

    PexPos->trimValue(pp->raw.expos);
}
