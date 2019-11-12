/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Ingo Weyrich (heckflosse67@gmx.de)
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

#include <cmath>
#include <iomanip>

#include "pdsharpening.h"

#include "eventmapper.h"
#include "options.h"
#include "rtimage.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

PdSharpening::PdSharpening() :
    FoldableToolPanel(this, "capturesharpening", M("TP_PDSHARPENING_LABEL"), false, true),
    lastAutoContrast(true),
    lastAutoRadius(true)
{
    auto m = ProcEventMapper::getInstance();
    EvPdShrContrast = m->newEvent(CAPTURESHARPEN, "HISTORY_MSG_PDSHARPEN_CONTRAST");
    EvPdShrRange = m->newEvent(CAPTURESHARPEN, "HISTORY_MSG_PDSHARPEN_RANGE");
    EvPdShrDRadius = m->newEvent(CAPTURESHARPEN, "HISTORY_MSG_PDSHARPEN_RADIUS");
    EvPdShrDRadiusOffset = m->newEvent(CAPTURESHARPEN, "HISTORY_MSG_PDSHARPEN_RADIUS_BOOST");
    EvPdShrDIterations = m->newEvent(CAPTURESHARPEN, "HISTORY_MSG_PDSHARPEN_ITERATIONS");
    EvPdShrAutoContrast = m->newEvent(CAPTURESHARPEN, "HISTORY_MSG_PDSHARPEN_AUTO_CONTRAST");
    EvPdShrAutoRadius = m->newEvent(CAPTURESHARPEN, "HISTORY_MSG_PDSHARPEN_AUTO_RADIUS");

    Gtk::HBox* hb = Gtk::manage(new Gtk::HBox());
    hb->show();
    contrast = Gtk::manage(new Adjuster(M("TP_SHARPENING_CONTRAST"), 0, 200, 1, 10));
    contrast->setAdjusterListener(this);
    contrast->addAutoButton();
    contrast->setAutoValue(true);

    pack_start(*contrast);
    contrast->show();

    Gtk::Image *al = Gtk::manage(new RTImage("circle-black-small.png"));
    Gtk::Image *ar = Gtk::manage(new RTImage("circle-white-small.png"));
    range = Gtk::manage(new Adjuster(M("TP_SHARPENING_RANGE"), 0, 100, 1, 95, al, ar));
    range->setAdjusterListener(this);

    pack_start(*range);
    range->show();

    pack_start(*hb);

    Gtk::VBox* rld = Gtk::manage(new Gtk::VBox());
    dradius = Gtk::manage(new Adjuster(M("TP_SHARPENING_RADIUS"), 0.4, 1.15, 0.01, 0.75));
    dradius->addAutoButton();
    dradius->setAutoValue(true);
    dradiusOffset = Gtk::manage(new Adjuster(M("TP_SHARPENING_RADIUS_BOOST"), -0.5, 0.5, 0.01, 0.0));
    diter = Gtk::manage(new Adjuster(M("TP_SHARPENING_RLD_ITERATIONS"), 1, 100, 1, 20));
    rld->pack_start(*dradius);
    rld->pack_start(*dradiusOffset);
    rld->pack_start(*diter);
    dradius->show();
    dradiusOffset->show();
    diter->show();
    rld->show();
    pack_start(*rld);

    dradius->setAdjusterListener(this);
    dradiusOffset->setAdjusterListener(this);
    diter->setAdjusterListener(this);
    range->setAdjusterListener(this);

    contrast->delay = std::max(contrast->delay, options.adjusterMaxDelay);
    range->delay = std::max(range->delay, options.adjusterMaxDelay);
    dradius->delay = std::max(dradius->delay, options.adjusterMaxDelay);
    dradiusOffset->delay = std::max(dradiusOffset->delay, options.adjusterMaxDelay);
    diter->delay = std::max(diter->delay, options.adjusterMaxDelay);
}

PdSharpening::~PdSharpening()
{
    idle_register.destroy();
}


void PdSharpening::read(const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener();

    if (pedited) {
        contrast->setEditedState(pedited->pdsharpening.contrast ? Edited : UnEdited);
        contrast->setAutoInconsistent(multiImage && !pedited->pdsharpening.autoContrast);
        dradius->setAutoInconsistent(multiImage && !pedited->pdsharpening.autoRadius);
        dradius->setEditedState(pedited->pdsharpening.deconvradius ? Edited : UnEdited);
        dradiusOffset->setEditedState(pedited->pdsharpening.deconvradiusOffset ? Edited : UnEdited);
        diter->setEditedState(pedited->pdsharpening.deconviter ? Edited : UnEdited);
        range->setEditedState(pedited->pdsharpening.deconvrange ? Edited : UnEdited);

        set_inconsistent(multiImage && !pedited->pdsharpening.enabled);
    }

    setEnabled(pp->pdsharpening.enabled);

    contrast->setValue(pp->pdsharpening.contrast);
    contrast->setAutoValue(pp->pdsharpening.autoContrast);
    dradius->setValue(pp->pdsharpening.deconvradius);
    dradius->setAutoValue(pp->pdsharpening.autoRadius);
    dradiusOffset->setValue(pp->pdsharpening.deconvradiusOffset);
    diter->setValue(pp->pdsharpening.deconviter);
    range->setValue(pp->pdsharpening.deconvrange);
    lastAutoContrast = pp->pdsharpening.autoContrast;
    lastAutoRadius = pp->pdsharpening.autoRadius;

    enableListener();
}

void PdSharpening::write(ProcParams* pp, ParamsEdited* pedited)
{

    pp->pdsharpening.contrast = contrast->getValue();
    pp->pdsharpening.autoContrast = contrast->getAutoValue();
    pp->pdsharpening.enabled = getEnabled();
    pp->pdsharpening.deconvradius = dradius->getValue();
    pp->pdsharpening.autoRadius = dradius->getAutoValue();
    pp->pdsharpening.deconvradiusOffset = dradiusOffset->getValue();
    pp->pdsharpening.deconviter =(int)diter->getValue();
    pp->pdsharpening.deconvrange =(int)range->getValue();

    if (pedited) {
        pedited->pdsharpening.contrast = contrast->getEditedState();
        pedited->pdsharpening.autoContrast = !contrast->getAutoInconsistent();
        pedited->pdsharpening.deconvradius = dradius->getEditedState();
        pedited->pdsharpening.autoRadius = !dradius->getAutoInconsistent();
        pedited->pdsharpening.deconvradiusOffset = dradiusOffset->getEditedState();
        pedited->pdsharpening.deconviter = diter->getEditedState();
        pedited->pdsharpening.deconvrange = range->getEditedState();
        pedited->pdsharpening.enabled = !get_inconsistent();
    }
}

void PdSharpening::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited)
{

    contrast->setDefault(defParams->pdsharpening.contrast);
    dradius->setDefault(defParams->pdsharpening.deconvradius);
    dradiusOffset->setDefault(defParams->pdsharpening.deconvradiusOffset);
    diter->setDefault(defParams->pdsharpening.deconviter);
    range->setDefault(defParams->pdsharpening.deconvrange);

    if (pedited) {
        contrast->setDefaultEditedState(pedited->pdsharpening.contrast ? Edited : UnEdited);
        dradius->setDefaultEditedState(pedited->pdsharpening.deconvradius ? Edited : UnEdited);
        dradiusOffset->setDefaultEditedState(pedited->pdsharpening.deconvradiusOffset ? Edited : UnEdited);
        diter->setDefaultEditedState(pedited->pdsharpening.deconviter ? Edited : UnEdited);
        range->setDefaultEditedState(pedited->pdsharpening.deconvrange ? Edited : UnEdited);
    } else {
        contrast->setDefaultEditedState(Irrelevant);
        dradius->setDefaultEditedState(Irrelevant);
        dradiusOffset->setDefaultEditedState(Irrelevant);
        diter->setDefaultEditedState(Irrelevant);
        range->setDefaultEditedState(Irrelevant);
    }
}

void PdSharpening::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && (multiImage || getEnabled())) {

        Glib::ustring costr;

        if (a == dradius || a == dradiusOffset) {
            costr = Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), a->getValue());
        } else {
            costr = Glib::ustring::format((int)a->getValue());
        }

        if (a == contrast) {
            listener->panelChanged(EvPdShrContrast, costr);
        } else if (a == dradius) {
            listener->panelChanged(EvPdShrDRadius, costr);
        } else if (a == dradiusOffset) {
            listener->panelChanged(EvPdShrDRadiusOffset, costr);
        } else if (a == diter) {
            listener->panelChanged(EvPdShrDIterations, costr);
        } else if (a == range) {
            listener->panelChanged(EvPdShrRange, costr);
        }
    }
}

void PdSharpening::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvPdShrEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvPdShrEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvPdShrEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void PdSharpening::setBatchMode(bool batchMode)
{

    ToolPanel::setBatchMode(batchMode);

    contrast->showEditedCB();
    dradius->showEditedCB();
    dradiusOffset->showEditedCB();
    diter->showEditedCB();
    range->showEditedCB();
}

void PdSharpening::setAdjusterBehavior(bool contrastadd, bool radiusadd, bool iteradd, bool rangeadd)
{

    contrast->setAddMode(contrastadd);
    dradius->setAddMode(radiusadd);
    dradiusOffset->setAddMode(radiusadd);
    diter->setAddMode(iteradd);
    range->setAddMode(rangeadd);
}

void PdSharpening::trimValues(rtengine::procparams::ProcParams* pp)
{

    contrast->trimValue(pp->pdsharpening.contrast);
    dradius->trimValue(pp->pdsharpening.deconvradius);
    dradiusOffset->trimValue(pp->pdsharpening.deconvradiusOffset);
    diter->trimValue(pp->pdsharpening.deconviter);
    range->trimValue(pp->pdsharpening.deconvrange);
}

void PdSharpening::autoContrastChanged(double autoContrast)
{
    idle_register.add(
        [this, autoContrast]() -> bool
        {
            disableListener();
            contrast->setValue(autoContrast);
            enableListener();
            return false;
        }
    );
}

void PdSharpening::autoRadiusChanged(double autoRadius)
{
    idle_register.add(
        [this, autoRadius]() -> bool
        {
            disableListener();
            dradius->setValue(autoRadius);
            enableListener();
            return false;
        }
    );
}

void PdSharpening::adjusterAutoToggled(Adjuster* a, bool newval)
{
    if (multiImage) {
        if (a->getAutoInconsistent()) {
            a->setAutoInconsistent(false);
            a->setAutoValue(false);
        } else if (lastAutoRadius) {
            a->setAutoInconsistent(true);
        }

        (a == contrast ? lastAutoContrast : lastAutoRadius) = a->getAutoValue();
    }

    if (listener) {
        const auto event = (a == contrast ? EvPdShrAutoContrast : EvPdShrAutoRadius);
        if (a->getAutoInconsistent()) {
            listener->panelChanged(event, M("GENERAL_UNCHANGED"));
        } else if (a->getAutoValue()) {
            listener->panelChanged(event, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(event, M("GENERAL_DISABLED"));
        }
    }
}
