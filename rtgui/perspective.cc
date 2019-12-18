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
#include "perspective.h"

#include "rtimage.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

PerspCorrection::PerspCorrection () : FoldableToolPanel(this, "perspective", M("TP_PERSPECTIVE_LABEL"))
{

    Gtk::Image* ipersHL =   Gtk::manage (new RTImage ("perspective-horizontal-left-small.png"));
    Gtk::Image* ipersHR =   Gtk::manage (new RTImage ("perspective-horizontal-right-small.png"));
    Gtk::Image* ipersVL =   Gtk::manage (new RTImage ("perspective-vertical-bottom-small.png"));
    Gtk::Image* ipersVR =   Gtk::manage (new RTImage ("perspective-vertical-top-small.png"));
    Gtk::Image* ipersBL =   Gtk::manage (new RTImage ("perspective-vertical-bottom-small.png"));
    Gtk::Image* ipersBR =   Gtk::manage (new RTImage ("perspective-vertical-top-small.png"));

    vert = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_VERTICAL"), -85, 85, 0.1, 0, ipersVL, ipersVR));
    vert->setAdjusterListener (this);

    horiz = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_HORIZONTAL"), -85, 85, 0.1, 0, ipersHL, ipersHR));
    horiz->setAdjusterListener (this);

    vBias = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_VERTICAL_BIAS"), -85, 85, 0.1, 0, ipersBL, ipersBR));
    vBias->setAdjusterListener (this);

    fov = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_FOV"), 0.1, 150, 0.1, 65));
    fov->setAdjusterListener (this);

    pack_start (*vert);
    pack_start (*horiz);
    pack_start (*vBias);
    pack_start (*fov);

    horiz->setLogScale(2, 0);
    vert->setLogScale(2, 0);
    vBias->setLogScale(2, 0);
    fov->setLogScale(2, 0);

    show_all();
}

void PerspCorrection::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        horiz->setEditedState (pedited->perspective.horizontal ? Edited : UnEdited);
        vert->setEditedState (pedited->perspective.vertical ? Edited : UnEdited);
        vBias->setEditedState (pedited->perspective.vBias ? Edited : UnEdited);
        fov->setEditedState (pedited->perspective.fov ? Edited : UnEdited);
    }

    horiz->setValue (pp->perspective.horizontal);
    vert->setValue (pp->perspective.vertical);
    vBias->setValue (pp->perspective.vBias);
    fov->setValue (pp->perspective.fov);

    enableListener ();
}

void PerspCorrection::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->perspective.horizontal  = horiz->getValue ();
    pp->perspective.vertical = vert->getValue ();
    pp->perspective.vBias = vBias->getValue ();
    pp->perspective.fov = fov->getValue ();

    if (pedited) {
        pedited->perspective.horizontal = horiz->getEditedState ();
        pedited->perspective.vertical = vert->getEditedState ();
        pedited->perspective.vBias = vBias->getEditedState ();
        pedited->perspective.fov = fov->getEditedState ();
    }
}

void PerspCorrection::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    horiz->setDefault (defParams->perspective.horizontal);
    vert->setDefault (defParams->perspective.vertical);
    vBias->setDefault (defParams->perspective.vBias);
    fov->setDefault (defParams->perspective.fov);

    if (pedited) {
        horiz->setDefaultEditedState (pedited->perspective.horizontal ? Edited : UnEdited);
        vert->setDefaultEditedState (pedited->perspective.vertical ? Edited : UnEdited);
        vBias->setDefaultEditedState (pedited->perspective.vBias ? Edited : UnEdited);
        fov->setDefaultEditedState (pedited->perspective.fov ? Edited : UnEdited);
    } else {
        horiz->setDefaultEditedState (Irrelevant);
        vert->setDefaultEditedState (Irrelevant);
        vBias->setDefaultEditedState (Irrelevant);
        fov->setDefaultEditedState (Irrelevant);
    }
}

void PerspCorrection::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        listener->panelChanged (EvPerspCorr, Glib::ustring::compose ("%1=%5\n%2=%6\n%3=%7\n%4=%8", M("TP_PERSPECTIVE_HORIZONTAL"), M("TP_PERSPECTIVE_VERTICAL"), M("TP_PERSPECTIVE_VERTICAL_BIAS"), M("TP_PERSPECTIVE_FOV"), horiz->getValue(), vert->getValue(), vBias->getValue(), fov->getValue()));
    }
}

void PerspCorrection::setAdjusterBehavior (bool badd)
{

    horiz->setAddMode(badd);
    vert->setAddMode(badd);
    vBias->setAddMode(badd);
    fov->setAddMode(badd);
}

void PerspCorrection::trimValues (rtengine::procparams::ProcParams* pp)
{

    horiz->trimValue(pp->perspective.horizontal);
    vert->trimValue(pp->perspective.vertical);
    vBias->trimValue(pp->perspective.vBias);
    fov->trimValue(pp->perspective.fov);
}

void PerspCorrection::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    horiz->showEditedCB ();
    vert->showEditedCB ();
    vBias->showEditedCB ();
    fov->showEditedCB ();
}
