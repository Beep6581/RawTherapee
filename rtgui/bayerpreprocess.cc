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
#include "bayerpreprocess.h"
#include "guiutils.h"
#include "eventmapper.h"
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

BayerPreProcess::BayerPreProcess () : FoldableToolPanel(this, "bayerpreprocess", M("TP_PREPROCESS_LABEL"), true)
{
    auto m = ProcEventMapper::getInstance();
    EvLineDenoiseDirection = m->newEvent(DARKFRAME, "HISTORY_MSG_PREPROCESS_LINEDENOISE_DIRECTION");
    EvPDAFLinesFilter = m->newEvent(DARKFRAME, "HISTORY_MSG_PREPROCESS_PDAFLINESFILTER");

    lineDenoise = Gtk::manage(new Adjuster (M("TP_PREPROCESS_LINEDENOISE"), 0, 1000, 1, 0));
    lineDenoise->setAdjusterListener (this);

    if (lineDenoise->delay < options.adjusterMaxDelay) {
        lineDenoise->delay = options.adjusterMaxDelay;
    }

    lineDenoise->show();

    greenEqThreshold = Gtk::manage(new Adjuster (M("TP_PREPROCESS_GREENEQUIL"), 0, 100, 1, 0));
    greenEqThreshold->setAdjusterListener (this);

    if (greenEqThreshold->delay < options.adjusterMaxDelay) {
        greenEqThreshold->delay = options.adjusterMaxDelay;
    }

    greenEqThreshold->show();

    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_PREPROCESS_LINEDENOISE_DIRECTION") + ": ")), Gtk::PACK_SHRINK, 0);
    lineDenoiseDirection = Gtk::manage(new MyComboBoxText());
    lineDenoiseDirection->append(M("TP_PREPROCESS_LINEDENOISE_DIRECTION_HORIZONTAL"));
    lineDenoiseDirection->append(M("TP_PREPROCESS_LINEDENOISE_DIRECTION_VERTICAL"));
    lineDenoiseDirection->append(M("TP_PREPROCESS_LINEDENOISE_DIRECTION_BOTH"));
    lineDenoiseDirection->append(M("TP_PREPROCESS_LINEDENOISE_DIRECTION_PDAF_LINES"));
    lineDenoiseDirection->show();
    lineDenoiseDirection->signal_changed().connect(sigc::mem_fun(*this, &BayerPreProcess::lineDenoiseDirectionChanged));
    
    hb->pack_start(*lineDenoiseDirection);

    pack_start( *lineDenoise, Gtk::PACK_SHRINK, 4);
    pack_start(*hb, Gtk::PACK_SHRINK, 4);

    pack_start( *Gtk::manage (new  Gtk::HSeparator()));

    pack_start( *greenEqThreshold, Gtk::PACK_SHRINK, 4);

    pdafLinesFilter = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_PDAFLINESFILTER"))));
    pdafLinesFilter->set_tooltip_markup(M("TP_PREPROCESS_PDAFLINESFILTER_TOOLTIP"));
    pdafLinesFilter->show();
    pdafLinesFilter->signal_toggled().connect(sigc::mem_fun(*this, &BayerPreProcess::pdafLinesFilterChanged), true);

    pack_start(*Gtk::manage(new Gtk::HSeparator()));
    pack_start(*pdafLinesFilter, Gtk::PACK_SHRINK, 4);
}

void BayerPreProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited ) {
        lineDenoise->setEditedState( pedited->raw.bayersensor.linenoise ? Edited : UnEdited );
        greenEqThreshold->setEditedState( pedited->raw.bayersensor.greenEq ? Edited : UnEdited );
        if (!pedited->raw.bayersensor.linenoiseDirection) {
            lineDenoiseDirection->set_active(3);
        }
        pdafLinesFilter->set_inconsistent(!pedited->raw.bayersensor.pdafLinesFilter);
    }

    lineDenoise->setValue (pp->raw.bayersensor.linenoise);
    int d = int(pp->raw.bayersensor.linenoiseDirection)-1;
    if (d == 4) {
        --d;
    }
    lineDenoiseDirection->set_active(d);
    greenEqThreshold->setValue (pp->raw.bayersensor.greenthresh);
    pdafLinesFilter->set_active(pp->raw.bayersensor.pdafLinesFilter);

    enableListener ();
}

void BayerPreProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.bayersensor.linenoise = lineDenoise->getIntValue();
    int d = lineDenoiseDirection->get_active_row_number() + 1;
    if (d == 4) {
        ++d;
    }
    pp->raw.bayersensor.linenoiseDirection = RAWParams::BayerSensor::LineNoiseDirection(d);
    pp->raw.bayersensor.greenthresh = greenEqThreshold->getIntValue();
    pp->raw.bayersensor.pdafLinesFilter = pdafLinesFilter->get_active();

    if (pedited) {
        pedited->raw.bayersensor.linenoise = lineDenoise->getEditedState ();
        pedited->raw.bayersensor.greenEq = greenEqThreshold->getEditedState ();
        pedited->raw.bayersensor.linenoise = lineDenoiseDirection->get_active_row_number() != 3;
        pedited->raw.bayersensor.pdafLinesFilter = !pdafLinesFilter->get_inconsistent();
    }
}

void BayerPreProcess::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {

        Glib::ustring value = a->getTextValue();

        if (a == greenEqThreshold) {
            listener->panelChanged (EvPreProcessGEquilThresh,  value );
        } else if (a == lineDenoise) {
            listener->panelChanged (EvPreProcessLineDenoise,  value );
        }
    }
}

void BayerPreProcess::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    lineDenoise->showEditedCB ();
    greenEqThreshold->showEditedCB ();
    if (batchMode) {
        lineDenoiseDirection->append(M("GENERAL_UNCHANGED"));
    }
}

void BayerPreProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    lineDenoise->setDefault( defParams->raw.bayersensor.linenoise);
    greenEqThreshold->setDefault (defParams->raw.bayersensor.greenthresh);

    if (pedited) {
        lineDenoise->setDefaultEditedState( pedited->raw.bayersensor.linenoise ? Edited : UnEdited);
        greenEqThreshold->setDefaultEditedState(pedited->raw.bayersensor.greenEq ? Edited : UnEdited);
    } else {
        lineDenoise->setDefaultEditedState( Irrelevant );
        greenEqThreshold->setDefaultEditedState(Irrelevant );
    }
}

void BayerPreProcess::setAdjusterBehavior (bool linedenoiseadd, bool greenequiladd)
{

    lineDenoise->setAddMode(linedenoiseadd);
    greenEqThreshold->setAddMode(greenequiladd);
}

void BayerPreProcess::trimValues (rtengine::procparams::ProcParams* pp)
{

    lineDenoise->trimValue(pp->raw.bayersensor.linenoise);
    greenEqThreshold->trimValue(pp->raw.bayersensor.greenthresh);
}


void BayerPreProcess::lineDenoiseDirectionChanged()
{
    if (listener) {
        listener->panelChanged(EvLineDenoiseDirection, lineDenoiseDirection->get_active_text());
    }
}

void BayerPreProcess::pdafLinesFilterChanged()
{
    if (listener) {
        listener->panelChanged(EvPDAFLinesFilter, pdafLinesFilter->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
