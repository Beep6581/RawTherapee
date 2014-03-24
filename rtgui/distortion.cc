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
#include "distortion.h"
#include <iomanip>
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

Distortion::Distortion (): FoldableToolPanel(this) {

    rlistener = NULL;
    autoDistor = Gtk::manage (new Gtk::Button (M("TP_DISTORTION_AUTO")));
    autoDistor->set_image (*Gtk::manage (new RTImage ("distortion-auto.png")));
    autoDistor->set_tooltip_text (M("TP_DISTORTION_AUTO_TIP"));
    idConn = autoDistor->signal_pressed().connect( sigc::mem_fun(*this, &Distortion::idPressed) );
    autoDistor->show();
    pack_start (*autoDistor);

    Gtk::Image* idistL =   Gtk::manage (new RTImage ("distortion-pincushion.png"));
    Gtk::Image* idistR =   Gtk::manage (new RTImage ("distortion-barrel.png"));

    distor = Gtk::manage (new Adjuster (M("TP_DISTORTION_AMOUNT"), -0.5, 0.5, 0.001, 0, idistL, idistR));
    distor->setAdjusterListener (this); 
    distor->show();
    pack_start (*distor);
}

void Distortion::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        distor->setEditedState (pedited->distortion.amount ? Edited : UnEdited);
    }

    distor->setValue (pp->distortion.amount);
    
    enableListener ();
}

void Distortion::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->distortion.amount = distor->getValue ();

    if (pedited) {
        pedited->distortion.amount = distor->getEditedState ();
    }
}

void Distortion::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    distor->setDefault (defParams->distortion.amount);

    if (pedited) {
        distor->setDefaultEditedState (pedited->distortion.amount ? Edited : UnEdited);
    }
    else  {
        distor->setDefaultEditedState (Irrelevant);
    }
}

void Distortion::adjusterChanged (Adjuster* a, double newval) {

    if (listener) 
        listener->panelChanged (EvDISTAmount, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), a->getValue()));
}

void Distortion::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    if (batchMode) {
        autoDistor->set_sensitive(false);
    }
    distor->showEditedCB ();
}

void Distortion::idPressed () {
    if (!batchMode) {
        if (rlistener) {
            double new_amount = rlistener->autoDistorRequested();
            distor->setValue(new_amount);
            adjusterChanged (distor, new_amount);
        }
    }
}

void Distortion::setAdjusterBehavior (bool vadd) {

	distor->setAddMode(vadd);
}

void Distortion::trimValues (rtengine::procparams::ProcParams* pp) {

	distor->trimValue(pp->distortion.amount);
}
