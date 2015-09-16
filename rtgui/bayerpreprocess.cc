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
#include "bayerpreprocess.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

BayerPreProcess::BayerPreProcess () : FoldableToolPanel(this, "bayerpreprocess", M("TP_PREPROCESS_LABEL"), true)
{
    lineDenoise = Gtk::manage(new Adjuster (M("TP_PREPROCESS_LINEDENOISE"), 0, 1000, 1, 0));
    lineDenoise->setAdjusterListener (this);

    if (lineDenoise->delay < options.adjusterMaxDelay) {
        lineDenoise->delay = options.adjusterMaxDelay;
    }

    lineDenoise->show();

    greenEqThreshold = Gtk::manage(new Adjuster (M("TP_PREPROCESS_GREENEQUIL"), 0, 100, 1, 0));
    greenEqThreshold->setAdjusterListener (this);

    if (greenEqThreshold->delay < options.adjusterMaxDelay) {
        greenEqThreshold->delay = options.adjusterMaxDelay;
    }

    greenEqThreshold->show();

    pack_start( *lineDenoise, Gtk::PACK_SHRINK, 4);

    pack_start( *Gtk::manage (new  Gtk::HSeparator()));

    pack_start( *greenEqThreshold, Gtk::PACK_SHRINK, 4);

}

void BayerPreProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited ) {
        lineDenoise->setEditedState( pedited->raw.bayersensor.linenoise ? Edited : UnEdited );
        greenEqThreshold->setEditedState( pedited->raw.bayersensor.greenEq ? Edited : UnEdited );
    }

    lineDenoise->setValue (pp->raw.bayersensor.linenoise);
    greenEqThreshold->setValue (pp->raw.bayersensor.greenthresh);

    enableListener ();
}

void BayerPreProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.bayersensor.linenoise = lineDenoise->getIntValue();
    pp->raw.bayersensor.greenthresh = greenEqThreshold->getIntValue();

    if (pedited) {
        pedited->raw.bayersensor.linenoise = lineDenoise->getEditedState ();
        pedited->raw.bayersensor.greenEq = greenEqThreshold->getEditedState ();
    }
}

void BayerPreProcess::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {

        Glib::ustring value = a->getTextValue();

        if (a == greenEqThreshold) {
            listener->panelChanged (EvPreProcessGEquilThresh,  value );
        } else if (a == lineDenoise) {
            listener->panelChanged (EvPreProcessLineDenoise,  value );
        }
    }
}

void BayerPreProcess::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    lineDenoise->showEditedCB ();
    greenEqThreshold->showEditedCB ();
}

void BayerPreProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    lineDenoise->setDefault( defParams->raw.bayersensor.linenoise);
    greenEqThreshold->setDefault (defParams->raw.bayersensor.greenthresh);

    if (pedited) {
        lineDenoise->setDefaultEditedState( pedited->raw.bayersensor.linenoise ? Edited : UnEdited);
        greenEqThreshold->setDefaultEditedState(pedited->raw.bayersensor.greenEq ? Edited : UnEdited);
    } else {
        lineDenoise->setDefaultEditedState( Irrelevant );
        greenEqThreshold->setDefaultEditedState(Irrelevant );
    }
}

void BayerPreProcess::setAdjusterBehavior (bool linedenoiseadd, bool greenequiladd)
{

    lineDenoise->setAddMode(linedenoiseadd);
    greenEqThreshold->setAddMode(greenequiladd);
}

void BayerPreProcess::trimValues (rtengine::procparams::ProcParams* pp)
{

    lineDenoise->trimValue(pp->raw.bayersensor.linenoise);
    greenEqThreshold->trimValue(pp->raw.bayersensor.greenthresh);
}
