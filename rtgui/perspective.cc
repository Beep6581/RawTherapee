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
#include "perspective.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

PerspCorrection::PerspCorrection () : FoldableToolPanel(this, "perspective", M("TP_PERSPECTIVE_LABEL")) {

    Gtk::Image* ipersHL =   Gtk::manage (new RTImage ("perspective-h1.png"));
    Gtk::Image* ipersHR =   Gtk::manage (new RTImage ("perspective-h2.png"));
    Gtk::Image* ipersVL =   Gtk::manage (new RTImage ("perspective-v1.png"));
    Gtk::Image* ipersVR =   Gtk::manage (new RTImage ("perspective-v2.png"));

	horiz = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_HORIZONTAL"), -100, 100, 0.1, 0, ipersHL, ipersHR));
	horiz->setAdjusterListener (this);

	vert = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_VERTICAL"), -100, 100, 0.1, 0, ipersVL, ipersVR));
	vert->setAdjusterListener (this);

    pack_start (*horiz);
    pack_start (*vert);

    show_all();
}

void PerspCorrection::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
    	horiz->setEditedState (pedited->perspective.horizontal ? Edited : UnEdited);
    	vert->setEditedState (pedited->perspective.vertical ? Edited : UnEdited);
    }

    horiz->setValue (pp->perspective.horizontal);
    vert->setValue (pp->perspective.vertical);

    enableListener ();
}

void PerspCorrection::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->perspective.horizontal  = horiz->getValue ();
    pp->perspective.vertical = vert->getValue ();

    if (pedited) {
        pedited->perspective.horizontal = horiz->getEditedState ();
        pedited->perspective.vertical = vert->getEditedState ();
    }
}

void PerspCorrection::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

	horiz->setDefault (defParams->perspective.horizontal);
	vert->setDefault (defParams->perspective.vertical);

    if (pedited) {
    	horiz->setDefaultEditedState (pedited->perspective.horizontal ? Edited : UnEdited);
    	vert->setDefaultEditedState (pedited->perspective.vertical ? Edited : UnEdited);
    }
    else {
    	horiz->setDefaultEditedState (Irrelevant);
    	vert->setDefaultEditedState (Irrelevant);
    }
}

void PerspCorrection::adjusterChanged (Adjuster* a, double newval) {

    if (listener)
        listener->panelChanged (EvPerspCorr, Glib::ustring::compose ("%1=%3\n%2=%4", M("TP_PERSPECTIVE_HORIZONTAL"), M("TP_PERSPECTIVE_VERTICAL"), horiz->getValue(), vert->getValue()));
}

void PerspCorrection::setAdjusterBehavior (bool badd) {

	horiz->setAddMode(badd);
	vert->setAddMode(badd);
}

void PerspCorrection::trimValues (rtengine::procparams::ProcParams* pp) {

	horiz->trimValue(pp->perspective.horizontal);
	vert->trimValue(pp->perspective.vertical);
}

void PerspCorrection::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    horiz->showEditedCB ();
    vert->showEditedCB ();
}
