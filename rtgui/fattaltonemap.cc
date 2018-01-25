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
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

FattalToneMapping::FattalToneMapping(): FoldableToolPanel(this, "fattal", M("TP_TM_FATTAL_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvTMFattalAnchor = m->newEvent(HDR, "HISTORY_MSG_TM_FATTAL_ANCHOR");

    amount = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_AMOUNT"), 1., 100., 1., 30.));
    threshold = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_THRESHOLD"), -100., 100., 1., 0.0));
    Gtk::Image *al = Gtk::manage(new RTImage("adj-black.png"));
    Gtk::Image *ar = Gtk::manage(new RTImage("adj-white.png"));
    anchor = Gtk::manage(new Adjuster(M("TP_TM_FATTAL_ANCHOR"), 1, 100, 1, 50, al, ar));

    amount->setAdjusterListener(this);
    threshold->setAdjusterListener(this);
    anchor->setAdjusterListener(this);

    amount->show();
    threshold->show();
    anchor->show();

    pack_start(*amount);
    pack_start(*threshold);
    pack_start(*anchor);
}

void FattalToneMapping::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        threshold->setEditedState(pedited->fattal.threshold ? Edited : UnEdited);
        amount->setEditedState(pedited->fattal.amount ? Edited : UnEdited);
        anchor->setEditedState(pedited->fattal.anchor ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->fattal.enabled);
    }

    setEnabled(pp->fattal.enabled);
    threshold->setValue(pp->fattal.threshold);
    amount->setValue(pp->fattal.amount);
    anchor->setValue(pp->fattal.anchor);

    enableListener();
}

void FattalToneMapping::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->fattal.threshold = threshold->getValue();
    pp->fattal.amount = amount->getValue();
    pp->fattal.anchor = anchor->getValue();
    pp->fattal.enabled = getEnabled();

    if(pedited) {
        pedited->fattal.threshold = threshold->getEditedState();
        pedited->fattal.amount = amount->getEditedState();
        pedited->fattal.anchor = anchor->getEditedState();
        pedited->fattal.enabled = !get_inconsistent();
    }
}

void FattalToneMapping::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    threshold->setDefault(defParams->fattal.threshold);
    amount->setDefault(defParams->fattal.amount);
    anchor->setDefault(defParams->fattal.anchor);

    if(pedited) {
        threshold->setDefaultEditedState(pedited->fattal.threshold ? Edited : UnEdited);
        amount->setDefaultEditedState(pedited->fattal.amount ? Edited : UnEdited);
        anchor->setDefaultEditedState(pedited->fattal.anchor ? Edited : UnEdited);
    } else {
        threshold->setDefaultEditedState(Irrelevant);
        amount->setDefaultEditedState(Irrelevant);
        anchor->setDefaultEditedState(Irrelevant);
    }
}

void FattalToneMapping::adjusterChanged(Adjuster* a, double newval)
{
    if(listener && getEnabled()) {
        if(a == threshold) {
            listener->panelChanged(EvTMFattalThreshold, a->getTextValue());
        } else if(a == amount) {
            listener->panelChanged(EvTMFattalAmount, a->getTextValue());
        } else if(a == anchor) {
            listener->panelChanged(EvTMFattalAnchor, a->getTextValue());
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
    anchor->showEditedCB();
}

void FattalToneMapping::setAdjusterBehavior(bool amountAdd, bool thresholdAdd, bool anchorAdd)
{
    amount->setAddMode(amountAdd);
    threshold->setAddMode(thresholdAdd);
    anchor->setAddMode(anchorAdd);
}

