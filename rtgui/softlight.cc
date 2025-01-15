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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <cmath>
#include <iomanip>

#include "softlight.h"

#include "eventmapper.h"

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring SoftLight::TOOL_NAME = "softlight";

SoftLight::SoftLight(): FoldableToolPanel(this, TOOL_NAME, M("TP_SOFTLIGHT_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvSoftLightEnabled = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_SOFTLIGHT_ENABLED");
    EvSoftLightStrength = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_SOFTLIGHT_STRENGTH");
    
    strength = Gtk::manage(new Adjuster(M("TP_SOFTLIGHT_STRENGTH"), 0., 100., 1., 30.));
    strength->setAdjusterListener(this);
    strength->show();

    pack_start(*strength);
}


void SoftLight::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        strength->setEditedState(pedited->softlight.strength ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->softlight.enabled);
    }

    setEnabled(pp->softlight.enabled);
    strength->setValue(pp->softlight.strength);

    enableListener();
}


void SoftLight::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->softlight.strength = strength->getValue();
    pp->softlight.enabled = getEnabled();

    if (pedited) {
        pedited->softlight.strength = strength->getEditedState();
        pedited->softlight.enabled = !get_inconsistent();
    }
}

void SoftLight::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    strength->setDefault(defParams->softlight.strength);

    if (pedited) {
        strength->setDefaultEditedState(pedited->softlight.strength ? Edited : UnEdited);
    } else {
        strength->setDefaultEditedState(Irrelevant);
    }
}


void SoftLight::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvSoftLightStrength, a->getTextValue());
    }
}

void SoftLight::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvSoftLightEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvSoftLightEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvSoftLightEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void SoftLight::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    strength->showEditedCB();
}


void SoftLight::setAdjusterBehavior(bool strengthAdd)
{
    strength->setAddMode(strengthAdd);
}

