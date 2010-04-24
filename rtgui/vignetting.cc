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

Vignetting::Vignetting () : vigAdd(false) {

    amount = Gtk::manage (new Adjuster (M("TP_VIGNETTING_AMOUNT"), -100, 100, 1, 0));
    amount->setAdjusterListener (this); 

    radius = Gtk::manage (new Adjuster (M("TP_VIGNETTING_RADIUS"), 0, 100, 1, 50));
    radius->setAdjusterListener (this); 

    pack_start (*amount);
    pack_start (*radius);

    show_all();
}

void Vignetting::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) 
        amount->setEditedState (pedited->vignetting.amount ? Edited : UnEdited);

    amount->setValue (pp->vignetting.amount);
    radius->setValue (pp->vignetting.radius);

    enableListener ();
}

void Vignetting::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->vignetting.amount = (int)amount->getValue ();
    pp->vignetting.radius = (int)radius->getValue ();

    if (pedited) 
        pedited->vignetting.amount = amount->getEditedState ();
}

void Vignetting::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    amount->setDefault (defParams->vignetting.amount);
    radius->setDefault (defParams->vignetting.radius);

    if (pedited) 
        amount->setDefaultEditedState (pedited->vignetting.amount ? Edited : UnEdited);
    else 
        amount->setDefaultEditedState (Irrelevant);
}

void Vignetting::adjusterChanged (Adjuster* a, double newval) {

    if (listener) 
        listener->panelChanged (EvVignetting, Glib::ustring::compose ("%1=%3\n%2=%4", M("TP_VIGNETTING_AMOUNT"), M("TP_VIGNETTING_RADIUS"), (int)amount->getValue(), (int)radius->getValue()));
}

void Vignetting::setAdjusterBehavior (bool bvadd) {

    if (!vigAdd && bvadd || vigAdd && !bvadd)
        amount->setLimits (-100, 100, 1, 0);
    
    vigAdd = bvadd;
}

void Vignetting::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    amount->showEditedCB ();
}
