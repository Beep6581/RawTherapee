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
#include <vignetting.h>
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

Vignetting::Vignetting () : Gtk::VBox(), FoldableToolPanel(this) {

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

void Vignetting::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) 
        amount->setEditedState (pedited->vignetting.amount ? Edited : UnEdited);

    amount->setValue (pp->vignetting.amount);
    radius->setValue (pp->vignetting.radius);
    strength->setValue (pp->vignetting.strength);
    centerX->setValue (pp->vignetting.centerX);
    centerY->setValue (pp->vignetting.centerY);

    enableListener ();
}

void Vignetting::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->vignetting.amount = (int)amount->getValue ();
    pp->vignetting.radius = (int)radius->getValue ();
    pp->vignetting.strength = (int)strength->getValue ();
    pp->vignetting.centerX = (int)centerX->getValue ();
    pp->vignetting.centerY = (int)centerY->getValue ();

    if (pedited) 
        pedited->vignetting.amount = amount->getEditedState ();
}

void Vignetting::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    amount->setDefault (defParams->vignetting.amount);
    radius->setDefault (defParams->vignetting.radius);
    strength->setDefault (defParams->vignetting.strength);
    centerX->setDefault (defParams->vignetting.centerX);
    centerY->setDefault (defParams->vignetting.centerY);

    if (pedited) 
        amount->setDefaultEditedState (pedited->vignetting.amount ? Edited : UnEdited);
    else 
        amount->setDefaultEditedState (Irrelevant);
}

void Vignetting::adjusterChanged (Adjuster* a, double newval) {

    if (listener) 
        listener->panelChanged (EvVignetting, Glib::ustring::compose ("%1=%5\n%2=%6\n%3=%7\n%4=%8 %9", M("TP_VIGNETTING_AMOUNT"), M("TP_VIGNETTING_RADIUS"), M("TP_VIGNETTING_STRENGTH"), M("TP_VIGNETTING_CENTER"), (int)amount->getValue(), (int)radius->getValue(), (int)strength->getValue(), (int)centerX->getValue(), (int)centerY->getValue()));
}

void Vignetting::setAdjusterBehavior (bool vadd) {

	amount->setAddMode(vadd);
}

void Vignetting::trimValues (rtengine::procparams::ProcParams* pp) {

	amount->trimValue(pp->vignetting.amount);
}

void Vignetting::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    amount->showEditedCB ();
}
