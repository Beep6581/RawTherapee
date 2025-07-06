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
#include "vignetting.h"

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Vignetting::TOOL_NAME = "vignetting";

Vignetting::Vignetting () : FoldableToolPanel(this, TOOL_NAME, M("TP_VIGNETTING_LABEL"))
{

    amount = Gtk::manage (new Adjuster (M("TP_VIGNETTING_AMOUNT"), -100, 100, 1, 0));
    amount->setAdjusterListener (this);

    radius = Gtk::manage (new Adjuster (M("TP_VIGNETTING_RADIUS"), 0, 100, 1, 50));
    radius->setAdjusterListener (this);

    strength = Gtk::manage (new Adjuster (M("TP_VIGNETTING_STRENGTH"), 1, 100, 1, 1));
    strength->setAdjusterListener (this);

    centerX = Gtk::manage (new Adjuster (M("TP_VIGNETTING_CENTER_X"), -100, 100, 1, 0));
    centerX->setAdjusterListener (this);

    centerY = Gtk::manage (new Adjuster (M("TP_VIGNETTING_CENTER_Y"), -100, 100, 1, 0));
    centerY->setAdjusterListener (this);

    pack_start (*amount);
    pack_start (*radius);
    pack_start (*strength);
    pack_start (*centerX);
    pack_start (*centerY);

    show_all();
}

void Vignetting::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        amount->setEditedState (pedited->vignetting.amount ? Edited : UnEdited);
        radius->setEditedState (pedited->vignetting.radius ? Edited : UnEdited);
        strength->setEditedState (pedited->vignetting.strength ? Edited : UnEdited);
        centerX->setEditedState (pedited->vignetting.centerX ? Edited : UnEdited);
        centerY->setEditedState (pedited->vignetting.centerY ? Edited : UnEdited);
    }

    amount->setValue (pp->vignetting.amount);
    radius->setValue (pp->vignetting.radius);
    strength->setValue (pp->vignetting.strength);
    centerX->setValue (pp->vignetting.centerX);
    centerY->setValue (pp->vignetting.centerY);

    enableListener ();
}

void Vignetting::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->vignetting.amount = amount->getIntValue ();
    pp->vignetting.radius = radius->getIntValue ();
    pp->vignetting.strength = strength->getIntValue ();
    pp->vignetting.centerX = centerX->getIntValue ();
    pp->vignetting.centerY = centerY->getIntValue ();

    if (pedited) {
        pedited->vignetting.amount = amount->getEditedState ();
        pedited->vignetting.radius = radius->getEditedState ();
        pedited->vignetting.strength = strength->getEditedState ();
        pedited->vignetting.centerX = centerX->getEditedState ();
        pedited->vignetting.centerY = centerY->getEditedState ();
    }
}

void Vignetting::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    amount->setDefault (defParams->vignetting.amount);
    radius->setDefault (defParams->vignetting.radius);
    strength->setDefault (defParams->vignetting.strength);
    centerX->setDefault (defParams->vignetting.centerX);
    centerY->setDefault (defParams->vignetting.centerY);

    if (pedited) {
        amount->setDefaultEditedState (pedited->vignetting.amount ? Edited : UnEdited);
        radius->setDefaultEditedState (pedited->vignetting.radius ? Edited : UnEdited);
        strength->setDefaultEditedState (pedited->vignetting.strength ? Edited : UnEdited);
        centerX->setDefaultEditedState (pedited->vignetting.centerX ? Edited : UnEdited);
        centerY->setDefaultEditedState (pedited->vignetting.centerY ? Edited : UnEdited);
    } else {
        amount->setDefaultEditedState (Irrelevant);
        radius->setDefaultEditedState (Irrelevant);
        strength->setDefaultEditedState (Irrelevant);
        centerX->setDefaultEditedState (Irrelevant);
        centerY->setDefaultEditedState (Irrelevant);
    }
}

void Vignetting::adjusterChanged(Adjuster* a, double newval)
{
    if (listener)  {
        if (a == amount) {
            listener->panelChanged (EvVignettingAmount, amount->getTextValue());
        } else if (a == radius) {
            listener->panelChanged (EvVignettingRadius, radius->getTextValue());
        } else if (a == strength) {
            listener->panelChanged (EvVignettingStrenght, strength->getTextValue());
        } else if (a == centerX || a == centerY) {
            listener->panelChanged (EvVignettingCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
        }
    }
}

void Vignetting::setAdjusterBehavior (bool amountadd, bool radiusadd, bool strengthadd, bool centeradd)
{

    amount->setAddMode(amountadd);
    radius->setAddMode(radiusadd);
    strength->setAddMode(strengthadd);
    centerX->setAddMode(centeradd);
    centerY->setAddMode(centeradd);
}

void Vignetting::trimValues (rtengine::procparams::ProcParams* pp)
{

    amount->trimValue(pp->vignetting.amount);
    radius->trimValue(pp->vignetting.radius);
    strength->trimValue(pp->vignetting.strength);
    centerX->trimValue(pp->vignetting.centerX);
    centerY->trimValue(pp->vignetting.centerY);
}

void Vignetting::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    amount->showEditedCB ();
    radius->showEditedCB ();
    strength->showEditedCB ();
    centerX->showEditedCB ();
    centerY->showEditedCB ();
}
