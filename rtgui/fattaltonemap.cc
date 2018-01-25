/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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
#include "fattaltonemap.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

FattalToneMapping::FattalToneMapping(): FoldableToolPanel(this, "fattal", M("TP_TM_FATTAL_LABEL"), true, true)
{
    amount = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_AMOUNT"), 1., 100., 1., 30.));
    threshold = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_THRESHOLD"), -100., 100., 1., 0.0));

    amount->setAdjusterListener(this);
    threshold->setAdjusterListener(this);

    amount->show();
    threshold->show();

    pack_start(*amount);
    pack_start(*threshold);
}

void FattalToneMapping::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if(pedited) {
        threshold->setEditedState(pedited->fattal.threshold ? Edited : UnEdited);
        amount->setEditedState(pedited->fattal.amount ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->fattal.enabled);
    }

    setEnabled(pp->fattal.enabled);
    threshold->setValue(pp->fattal.threshold);
    amount->setValue(pp->fattal.amount);

    enableListener();
}

void FattalToneMapping::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->fattal.threshold = threshold->getValue();
    pp->fattal.amount = amount->getValue();
    pp->fattal.enabled = getEnabled();

    if(pedited) {
        pedited->fattal.threshold = threshold->getEditedState();
        pedited->fattal.amount = amount->getEditedState();
        pedited->fattal.enabled = !get_inconsistent();
    }
}

void FattalToneMapping::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    threshold->setDefault(defParams->fattal.threshold);
    amount->setDefault(defParams->fattal.amount);

    if(pedited) {
        threshold->setDefaultEditedState(pedited->fattal.threshold ? Edited : UnEdited);
        amount->setDefaultEditedState(pedited->fattal.amount ? Edited : UnEdited);
    } else {
        threshold->setDefaultEditedState(Irrelevant);
        amount->setDefaultEditedState(Irrelevant);
    }
}

void FattalToneMapping::adjusterChanged(Adjuster* a, double newval)
{
    if(listener && getEnabled()) {
        if(a == threshold) {
            listener->panelChanged(EvTMFattalThreshold, a->getTextValue());
        } else if(a == amount) {
            listener->panelChanged(EvTMFattalAmount, a->getTextValue());
        }
    }
}

void FattalToneMapping::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvTMFattalEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvTMFattalEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvTMFattalEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void FattalToneMapping::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    threshold->showEditedCB();
    amount->showEditedCB();
}

void FattalToneMapping::setAdjusterBehavior (bool alphaAdd, bool betaAdd)
{
    threshold->setAddMode(alphaAdd);
    amount->setAddMode(betaAdd);
}

