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
#include "localcontrast.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

LocalContrast::LocalContrast(): FoldableToolPanel(this, "localcontrast", M("TP_LOCALCONTRAST_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvLocalContrastEnabled = m->newEvent(RGBCURVE, "HISTORY_MSG_LOCALCONTRAST_ENABLED");
    EvLocalContrastRadius = m->newEvent(RGBCURVE, "HISTORY_MSG_LOCALCONTRAST_RADIUS");
    EvLocalContrastAmount = m->newEvent(RGBCURVE, "HISTORY_MSG_LOCALCONTRAST_AMOUNT");
    EvLocalContrastDarkness = m->newEvent(RGBCURVE, "HISTORY_MSG_LOCALCONTRAST_DARKNESS");
    EvLocalContrastLightness = m->newEvent(RGBCURVE, "HISTORY_MSG_LOCALCONTRAST_LIGHTNESS");
    
    radius = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_RADIUS"), 3., 200., 1., 80.));
    amount = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_AMOUNT"), 0., 1., 0.01, 0.2));
    darkness = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_DARKNESS"), 0., 3., 0.01, 1.));
    lightness = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_LIGHTNESS"), 0., 3., 0.01, 1.));

    radius->setAdjusterListener(this);
    amount->setAdjusterListener(this);
    darkness->setAdjusterListener(this);
    lightness->setAdjusterListener(this);

    radius->show();
    amount->show();
    darkness->show();
    lightness->show();

    pack_start(*radius);
    pack_start(*amount);
    pack_start(*darkness);
    pack_start(*lightness);
}


void LocalContrast::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        radius->setEditedState(pedited->localContrast.radius ? Edited : UnEdited);
        amount->setEditedState(pedited->localContrast.amount ? Edited : UnEdited);
        darkness->setEditedState(pedited->localContrast.darkness ? Edited : UnEdited);
        lightness->setEditedState(pedited->localContrast.lightness ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->localContrast.enabled);
    }

    setEnabled(pp->localContrast.enabled);
    radius->setValue(pp->localContrast.radius);
    amount->setValue(pp->localContrast.amount);
    darkness->setValue(pp->localContrast.darkness);
    lightness->setValue(pp->localContrast.lightness);

    enableListener();
}


void LocalContrast::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->localContrast.radius = radius->getValue();
    pp->localContrast.amount = amount->getValue();
    pp->localContrast.darkness = darkness->getValue();
    pp->localContrast.lightness = lightness->getValue();
    pp->localContrast.enabled = getEnabled();

    if (pedited) {
        pedited->localContrast.radius = radius->getEditedState();
        pedited->localContrast.amount = amount->getEditedState();
        pedited->localContrast.darkness = darkness->getEditedState();
        pedited->localContrast.lightness = lightness->getEditedState();
        pedited->localContrast.enabled = !get_inconsistent();
    }
}

void LocalContrast::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    radius->setDefault(defParams->localContrast.radius);
    amount->setDefault(defParams->localContrast.amount);
    darkness->setDefault(defParams->localContrast.darkness);
    lightness->setDefault(defParams->localContrast.lightness);

    if (pedited) {
        radius->setDefaultEditedState(pedited->localContrast.radius ? Edited : UnEdited);
        amount->setDefaultEditedState(pedited->localContrast.amount ? Edited : UnEdited);
        darkness->setDefaultEditedState(pedited->localContrast.darkness ? Edited : UnEdited);
        lightness->setDefaultEditedState(pedited->localContrast.lightness ? Edited : UnEdited);
    } else {
        radius->setDefaultEditedState(Irrelevant);
        amount->setDefaultEditedState(Irrelevant);
        darkness->setDefaultEditedState(Irrelevant);
        lightness->setDefaultEditedState(Irrelevant);
    }
}


void LocalContrast::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == radius) {
            listener->panelChanged(EvLocalContrastRadius, a->getTextValue());
        } else if (a == amount) {
            listener->panelChanged(EvLocalContrastAmount, a->getTextValue());
        } else if (a == darkness) {
            listener->panelChanged(EvLocalContrastDarkness, a->getTextValue());
        } else if (a == lightness) {
            listener->panelChanged(EvLocalContrastLightness, a->getTextValue());
        }
    }
}


void LocalContrast::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void LocalContrast::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    radius->showEditedCB();
    amount->showEditedCB();
    darkness->showEditedCB();
    lightness->showEditedCB();
}


void LocalContrast::setAdjusterBehavior(bool radiusAdd, bool amountAdd, bool darknessAdd, bool lightnessAdd)
{
    radius->setAddMode(radiusAdd);
    amount->setAddMode(amountAdd);
    darkness->setAddMode(darknessAdd);
    lightness->setAddMode(lightnessAdd);
}

