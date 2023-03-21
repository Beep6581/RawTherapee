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

#include "distortion.h"

#include "rtimage.h"

#include "eventmapper.h"

#include "../rtengine/procparams.h"
#include "../rtengine/refreshmap.h"

using namespace rtengine;
using namespace rtengine::procparams;

Distortion::Distortion (): FoldableToolPanel(this, "distortion", M("TP_DISTORTION_LABEL"))
{
    auto mapper = ProcEventMapper::getInstance();
    EvDistortionDefish = mapper->newEvent(TRANSFORM, "HISTORY_MSG_DISTORTION_DEFISH");
    EvDistortionDefishVoid = mapper->newEvent(M_VOID, "HISTORY_MSG_DISTORTION_DEFISH");

    rlistener = nullptr;
    autoDistor = Gtk::manage (new Gtk::Button (M("GENERAL_AUTO")));
    autoDistor->set_image (*Gtk::manage (new RTImage ("distortion-auto-small.png")));
    autoDistor->get_style_context()->add_class("independent");
    autoDistor->set_alignment(0.5f, 0.5f);
    autoDistor->set_tooltip_text (M("TP_DISTORTION_AUTO_TOOLTIP"));
    idConn = autoDistor->signal_pressed().connect( sigc::mem_fun(*this, &Distortion::idPressed) );
    autoDistor->show();
    pack_start (*autoDistor);

    Gtk::Image* idistL =   Gtk::manage (new RTImage ("distortion-pincushion-small.png"));
    Gtk::Image* idistR =   Gtk::manage (new RTImage ("distortion-barrel-small.png"));

    distor = Gtk::manage (new Adjuster (M("TP_DISTORTION_AMOUNT"), -0.5, 0.5, 0.001, 0, idistL, idistR));
    distor->setAdjusterListener (this);
    distor->setLogScale(2, 0);
    distor->show();
    pack_start (*distor);

    defish = Gtk::manage(new Gtk::CheckButton(M("TP_DISTORTION_DEFISH")));
    defish->signal_toggled().connect(sigc::mem_fun(*this, &Distortion::defishChanged));
    defish->show();
    pack_start (*defish);

    focal_length = Gtk::manage (new Adjuster (M("TP_DISTORTION_FOCAL_LENGTH"), 0.5, 25, 0.01, DistortionParams::DEFAULT_FOCAL_LENGTH));
    focal_length->setAdjusterListener (this);
    focal_length->show();
    pack_start(*focal_length);
}

void Distortion::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        distor->setEditedState (pedited->distortion.amount ? Edited : UnEdited);
        focal_length->setEditedState (pedited->distortion.focal_length != DistortionParams::DEFAULT_FOCAL_LENGTH ? Edited : UnEdited);
    }

    distor->setValue (pp->distortion.amount);
    defish->set_active(pp->distortion.defish);
    focal_length->setValue(pp->distortion.focal_length);

    enableListener ();
}

void Distortion::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->distortion.amount = distor->getValue ();
    pp->distortion.defish = defish->get_active ();
    pp->distortion.focal_length = focal_length->getValue ();

    if (pedited) {
        pedited->distortion.amount = distor->getEditedState ();
        pedited->distortion.focal_length = focal_length->getEditedState ();
        pedited->distortion.defish = true;
    }
}

void Distortion::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    distor->setDefault (defParams->distortion.amount);
    focal_length->setDefault (defParams->distortion.focal_length);

    if (pedited) {
        distor->setDefaultEditedState (pedited->distortion.amount ? Edited : UnEdited);
    } else  {
        distor->setDefaultEditedState (Irrelevant);
    }
}

void Distortion::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        listener->panelChanged (EvDISTAmount, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), a->getValue()));
    }
}

void Distortion::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);

    if (batchMode) {
        autoDistor->set_sensitive(false);
    }

    distor->showEditedCB ();
}

void Distortion::idPressed ()
{
    if (!batchMode) {
        if (rlistener) {
            double new_amount = rlistener->autoDistorRequested();
            distor->setValue(new_amount);
            adjusterChanged (distor, new_amount);
        }
    }
}

void Distortion::defishChanged()
{
    if (listener) {
        listener->panelChanged(EvDistortionDefish, defish->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void Distortion::focalLengthChanged(Adjuster *a, const double newval)
{
}

void Distortion::setAdjusterBehavior (bool vadd)
{
    distor->setAddMode(vadd);
}

void Distortion::trimValues (rtengine::procparams::ProcParams* pp)
{

    distor->trimValue(pp->distortion.amount);
    focal_length->trimValue(pp->distortion.focal_length);
}
