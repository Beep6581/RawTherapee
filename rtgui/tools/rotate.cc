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
#include <iomanip>

#include "rotate.h"

#include "guiutils.h"
#include "lensgeomlistener.h"
#include "rtimage.h"

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Rotate::TOOL_NAME = "rotate";

Rotate::Rotate () : FoldableToolPanel(this, TOOL_NAME, M("TP_ROTATE_LABEL"))
{

    rlistener = nullptr;

    //TODO the action of the rotation slider is counter-intuitive
    Gtk::Image* irotateL =   Gtk::manage (new RTImage ("rotate-right-small"));
    Gtk::Image* irotateR =   Gtk::manage (new RTImage ("rotate-left-small"));

    degree = Gtk::manage (new Adjuster (M("TP_ROTATE_DEGREE"), -45, 45, 0.01, 0, irotateL, irotateR));
    degree->setAdjusterListener (this);
    pack_start (*degree);

    selectStraight = Gtk::manage (new Gtk::Button (M("TP_ROTATE_SELECTLINE")));
    selectStraight->set_image (*Gtk::manage (new RTImage ("rotate-straighten-small")));
    selectStraight->get_style_context()->add_class("independent");
    pack_start (*selectStraight, Gtk::PACK_SHRINK, 2);

    selectStraight->signal_pressed().connect( sigc::mem_fun(*this, &Rotate::selectStraightPressed) );

    degree->setLogScale(2, 0);

    show_all ();
}

void Rotate::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        degree->setEditedState (pedited->rotate.degree ? Edited : UnEdited);
    }

    degree->setValue (pp->rotate.degree);

    enableListener ();
}

void Rotate::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->rotate.degree = degree->getValue ();

    if (pedited) {
        pedited->rotate.degree = degree->getEditedState ();
    }
}

void Rotate::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    degree->setDefault (defParams->rotate.degree);

    if (pedited) {
        degree->setDefaultEditedState (pedited->rotate.degree ? Edited : UnEdited);
    } else {
        degree->setDefaultEditedState (Irrelevant);
    }
}

void Rotate::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        listener->panelChanged(EvROTDegree, Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), degree->getValue()));
    }
}

void Rotate::straighten (double deg)
{

    degree->setValue (degree->getValue() + deg);
    degree->setEditedState (Edited);

    if (listener) {
        listener->panelChanged (EvROTDegree, Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), degree->getValue()));
    }
}

void Rotate::selectStraightPressed ()
{

    if (rlistener) {
        rlistener->straightenRequested ();
    }
}

void Rotate::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    degree->showEditedCB ();
}

void Rotate::setAdjusterBehavior (bool rotadd)
{

    degree->setAddMode(rotadd);
}

void Rotate::trimValues (rtengine::procparams::ProcParams* pp)
{

    degree->trimValue(pp->rotate.degree);
}
