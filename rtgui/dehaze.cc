/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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
#include "dehaze.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

Dehaze::Dehaze(): FoldableToolPanel(this, "dehaze", M("TP_DEHAZE_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvDehazeEnabled = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_ENABLED");
    EvDehazeStrength = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_STRENGTH");
    
    strength = Gtk::manage(new Adjuster(M("TP_DEHAZE_STRENGTH"), 0., 100., 1., 50.));
    strength->setAdjusterListener(this);
    strength->show();

    pack_start(*strength);
}


void Dehaze::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        strength->setEditedState(pedited->dehaze.strength ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->dehaze.enabled);
    }

    setEnabled(pp->dehaze.enabled);
    strength->setValue(pp->dehaze.strength);

    enableListener();
}


void Dehaze::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->dehaze.strength = strength->getValue();
    pp->dehaze.enabled = getEnabled();

    if (pedited) {
        pedited->dehaze.strength = strength->getEditedState();
        pedited->dehaze.enabled = !get_inconsistent();
    }
}

void Dehaze::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    strength->setDefault(defParams->dehaze.strength);

    if (pedited) {
        strength->setDefaultEditedState(pedited->dehaze.strength ? Edited : UnEdited);
    } else {
        strength->setDefaultEditedState(Irrelevant);
    }
}


void Dehaze::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvDehazeStrength, a->getTextValue());
    }
}


void Dehaze::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Dehaze::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    strength->showEditedCB();
}


void Dehaze::setAdjusterBehavior(bool strengthAdd)
{
    strength->setAddMode(strengthAdd);
}

