/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2020 Alberto Romei <aldrop8@gmail.com>
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
#include <sstream>

#include "preprocesswb.h"
#include "eventmapper.h"

#include "guiutils.h"
#include "options.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

PreprocessWB::PreprocessWB() :
    FoldableToolPanel(this, "preprocesswb", M("TP_PREPROCWB_LABEL")),
    evPreprocessWBMode(ProcEventMapper::getInstance()->newEvent(FIRST, "HISTORY_MSG_PREPROCWB_MODE")),
    evPreprocessWBMults(ProcEventMapper::getInstance()->newEvent(FIRST, "HISTORY_MSG_PREPROCWB_MULTS")),
    mode(Gtk::manage(new MyComboBoxText())),
    red(Gtk::manage(new Adjuster(M("TP_PREPROCWB_RED"), 0.05, 20.0, 0.01, 1))),
    blue(Gtk::manage(new Adjuster(M("TP_PREPROCWB_BLUE"), 0.05, 20.0, 0.01, 1)))
{
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_PREPROCWB_MODE") + ": ")), Gtk::PACK_SHRINK, 0);

    mode->append(M("TP_PREPROCWB_MODE_CAMERA"));
    mode->append(M("TP_PREPROCWB_MODE_AUTO"));

    hb->pack_start(*mode);

    mode->set_active(0);
    mode->signal_changed().connect(sigc::mem_fun(*this, &PreprocessWB::modeChanged));

    red->setAdjusterListener(this);
    blue->setAdjusterListener(this);

    red->setLogScale(8, 1, true);
    blue->setLogScale(8, 1, true);

    if (red->delay < options.adjusterMaxDelay) {
        red->delay = options.adjusterMaxDelay;
    }

    if (blue->delay < options.adjusterMaxDelay) {
        blue->delay = options.adjusterMaxDelay;
    }

    mode->show();
    red->show();
    blue->show();

    pack_start(*hb, Gtk::PACK_SHRINK, 4);
    pack_start(*red, Gtk::PACK_SHRINK, 4);
    pack_start(*blue, Gtk::PACK_SHRINK, 4);
}

void PreprocessWB::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener();

    if (pedited) {
        red->setEditedState(pedited->raw.preprocessWB.red ? Edited : UnEdited);
        blue->setEditedState(pedited->raw.preprocessWB.blue ? Edited : UnEdited);
    }

    mode->set_active(int(pp->raw.preprocessWB.mode));
    red->setValue(pp->raw.preprocessWB.red);
    blue->setValue(pp->raw.preprocessWB.blue);

    enableListener();
}

void PreprocessWB::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    if (mode->get_active_row_number() != 2) {
        pp->raw.preprocessWB.mode = RAWParams::PreprocessWB::Mode(mode->get_active_row_number());
    }

    pp->raw.preprocessWB.red = red->getValue();
    pp->raw.preprocessWB.blue = blue->getValue();

    if (pedited) {
        pedited->raw.preprocessWB.mode = mode->get_active_row_number() != 2; // UNCHANGED entry, see setBatchMode
        pedited->raw.preprocessWB.red = red->getEditedState();
        pedited->raw.preprocessWB.blue = blue->getEditedState();
    }

}

void PreprocessWB::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        listener->panelChanged(evPreprocessWBMults,  Glib::ustring::compose(
                                   "R/G=%1 ; B/G=%2", red->getValue(), blue->getValue()));
    }
}

void PreprocessWB::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    if (batchMode) {
        mode->append(M("GENERAL_UNCHANGED"));
        mode->set_active_text(M("GENERAL_UNCHANGED"));
        red->showEditedCB();
        blue->showEditedCB();
    }
}

void PreprocessWB::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    red->setDefault(defParams->raw.preprocessWB.red);
    blue->setDefault(defParams->raw.preprocessWB.blue);

    if (pedited) {
        red->setDefaultEditedState(pedited->raw.preprocessWB.red ? Edited : UnEdited);
        blue->setDefaultEditedState(pedited->raw.preprocessWB.blue ? Edited : UnEdited);
    } else {
        red->setDefaultEditedState(Irrelevant);
        blue->setDefaultEditedState(Irrelevant);
    }
}

void PreprocessWB::setAdjusterBehavior(bool add)
{
    red->setAddMode(add);
    blue->setAddMode(add);
}

void PreprocessWB::trimValues(rtengine::procparams::ProcParams* pp)
{
    red->trimValue(pp->raw.preprocessWB.red);
    blue->trimValue(pp->raw.preprocessWB.blue);
}

void PreprocessWB::modeChanged()
{
    if (listener) {
        listener->panelChanged(evPreprocessWBMode, mode->get_active_text());
    }
}
