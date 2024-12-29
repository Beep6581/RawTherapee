/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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
#include <cmath>
#include <iomanip>

#include "dehaze.h"

#include "eventmapper.h"

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Dehaze::TOOL_NAME = "dehaze";

Dehaze::Dehaze(): FoldableToolPanel(this, TOOL_NAME, M("TP_DEHAZE_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvDehazeEnabled = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_ENABLED");
    EvDehazeStrength = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_STRENGTH");
    EvDehazeShowDepthMap = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_SHOW_DEPTH_MAP");
    EvDehazeDepth = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_DEPTH");
    EvDehazeSaturation = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_SATURATION");
    
    strength = Gtk::manage(new Adjuster(M("TP_DEHAZE_STRENGTH"), 0., 100., 1., 50.));
    strength->setAdjusterListener(this);
    strength->show();

    depth = Gtk::manage(new Adjuster(M("TP_DEHAZE_DEPTH"), 0., 100., 1., 25.));
    depth->setAdjusterListener(this);
    depth->show();

    saturation = Gtk::manage(new Adjuster(M("TP_DEHAZE_SATURATION"), 0., 100., 1., 50.));
    saturation->setAdjusterListener(this);
    saturation->show();

    showDepthMap = Gtk::manage(new Gtk::CheckButton(M("TP_DEHAZE_SHOW_DEPTH_MAP")));
    showDepthMap->signal_toggled().connect(sigc::mem_fun(*this, &Dehaze::showDepthMapChanged));
    showDepthMap->show();
    
    pack_start(*strength);
    pack_start(*depth);
    pack_start(*saturation);
    pack_start(*showDepthMap);
}

void Dehaze::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        strength->setEditedState(pedited->dehaze.strength ? Edited : UnEdited);
        saturation->setEditedState(pedited->dehaze.saturation ? Edited : UnEdited);
        depth->setEditedState(pedited->dehaze.depth ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->dehaze.enabled);
        showDepthMap->set_inconsistent(!pedited->dehaze.showDepthMap);
    }

    setEnabled(pp->dehaze.enabled);
    strength->setValue(pp->dehaze.strength);
    saturation->setValue(pp->dehaze.saturation);
    depth->setValue(pp->dehaze.depth);
    showDepthMap->set_active(pp->dehaze.showDepthMap);

    enableListener();
}

void Dehaze::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->dehaze.strength = strength->getValue();
    pp->dehaze.saturation = saturation->getValue();
    pp->dehaze.depth = depth->getValue();
    pp->dehaze.enabled = getEnabled();
    pp->dehaze.showDepthMap = showDepthMap->get_active();

    if (pedited) {
        pedited->dehaze.strength = strength->getEditedState();
        pedited->dehaze.saturation = saturation->getEditedState();
        pedited->dehaze.depth = depth->getEditedState();
        pedited->dehaze.enabled = !get_inconsistent();
        pedited->dehaze.showDepthMap = !showDepthMap->get_inconsistent();
    }
}

void Dehaze::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    strength->setDefault(defParams->dehaze.strength);
    saturation->setDefault(defParams->dehaze.saturation);
    depth->setDefault(defParams->dehaze.depth);

    if (pedited) {
        saturation->setDefaultEditedState(pedited->dehaze.saturation ? Edited : UnEdited);
        strength->setDefaultEditedState(pedited->dehaze.strength ? Edited : UnEdited);
        depth->setDefaultEditedState(pedited->dehaze.depth ? Edited : UnEdited);
    } else {
        saturation->setDefaultEditedState(Irrelevant);
        strength->setDefaultEditedState(Irrelevant);
        depth->setDefaultEditedState(Irrelevant);
    }
}

void Dehaze::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == strength) {
            listener->panelChanged(EvDehazeStrength, a->getTextValue());
        } else if (a == depth) {
            listener->panelChanged(EvDehazeDepth, a->getTextValue());
        } else if (a == saturation) {
            listener->panelChanged(EvDehazeSaturation, a->getTextValue());
        }
    }
}

void Dehaze::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Dehaze::showDepthMapChanged()
{
    if (listener) {
        listener->panelChanged(EvDehazeShowDepthMap, showDepthMap->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void Dehaze::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    strength->showEditedCB();
    depth->showEditedCB();
}

void Dehaze::setAdjusterBehavior(bool strengthAdd)
{
    strength->setAddMode(strengthAdd);
}
