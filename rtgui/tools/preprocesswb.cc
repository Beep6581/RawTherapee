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

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring PreprocessWB::TOOL_NAME = "preprocesswb";

PreprocessWB::PreprocessWB() :
    FoldableToolPanel(this, TOOL_NAME, M("TP_PREPROCWB_LABEL")),
    evPreprocessWBMode(ProcEventMapper::getInstance()->newEvent(FIRST, "HISTORY_MSG_PREPROCWB_MODE")),
    mode(Gtk::manage(new MyComboBoxText()))
{
    Gtk::Box *hb = Gtk::manage(new Gtk::Box());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_PREPROCWB_MODE") + ": ")), Gtk::PACK_SHRINK, 0);

    mode->append(M("TP_PREPROCWB_MODE_CAMERA"));
    mode->append(M("TP_PREPROCWB_MODE_AUTO"));

    hb->pack_start(*mode);

    mode->set_active(0);
    mode->signal_changed().connect(sigc::mem_fun(*this, &PreprocessWB::modeChanged));

    mode->show();

    pack_start(*hb, Gtk::PACK_SHRINK, 4);
}

void PreprocessWB::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener();

    mode->set_active(int(pp->raw.preprocessWB.mode));

    enableListener();
}

void PreprocessWB::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    if (mode->get_active_row_number() != 2) {
        pp->raw.preprocessWB.mode = RAWParams::PreprocessWB::Mode(mode->get_active_row_number());
    }

    if (pedited) {
        pedited->raw.preprocessWB.mode = mode->get_active_row_number() != 2; // UNCHANGED entry, see setBatchMode
    }

}

void PreprocessWB::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    if (batchMode) {
        mode->append(M("GENERAL_UNCHANGED"));
        mode->set_active_text(M("GENERAL_UNCHANGED"));
    }
}

void PreprocessWB::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
}

void PreprocessWB::trimValues(rtengine::procparams::ProcParams* pp)
{
}

void PreprocessWB::modeChanged()
{
    if (listener) {
        listener->panelChanged(evPreprocessWBMode, mode->get_active_text());
    }
}
