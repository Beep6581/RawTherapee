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

//    setEnabledTooltipMarkup(M("TP_EPD_TOOLTIP"));
    
    alpha = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_ALPHA"), 0.0, 2.0, 0.01, 1.0));
    beta = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_BETA"), 0.0, 2.0, 0.01, 1.0));

    alpha->setAdjusterListener(this);
    beta->setAdjusterListener(this);

    alpha->show();
    beta->show();

    pack_start(*alpha);
    pack_start(*beta);
}

void FattalToneMapping::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if(pedited) {
        alpha->setEditedState(pedited->fattal.alpha ? Edited : UnEdited);
        beta->setEditedState(pedited->fattal.beta ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->fattal.enabled);
    }

    setEnabled(pp->fattal.enabled);
    alpha->setValue(pp->fattal.alpha);
    beta->setValue(pp->fattal.beta);

    enableListener();
}

void FattalToneMapping::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->fattal.alpha = alpha->getValue();
    pp->fattal.beta = beta->getValue();
    pp->fattal.enabled = getEnabled();

    if(pedited) {
        pedited->fattal.alpha = alpha->getEditedState();
        pedited->fattal.beta = beta->getEditedState();
        pedited->fattal.enabled = !get_inconsistent();
    }
}

void FattalToneMapping::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    alpha->setDefault(defParams->fattal.alpha);
    beta->setDefault(defParams->fattal.beta);

    if(pedited) {
        alpha->setDefaultEditedState(pedited->fattal.alpha ? Edited : UnEdited);
        beta->setDefaultEditedState(pedited->fattal.beta ? Edited : UnEdited);
    } else {
        alpha->setDefaultEditedState(Irrelevant);
        beta->setDefaultEditedState(Irrelevant);
    }
}

void FattalToneMapping::adjusterChanged(Adjuster* a, double newval)
{
    if(listener && getEnabled()) {
        if(a == alpha) {
            listener->panelChanged(EvTMFattalAlpha, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        } else if(a == beta) {
            listener->panelChanged(EvTMFattalBeta, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
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

    alpha->showEditedCB();
    beta->showEditedCB();
}

