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
#include "xtransprocess.h"
#include "options.h"
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;

XTransProcess::XTransProcess () : FoldableToolPanel(this, "xtransprocess", M("TP_RAW_LABEL"), true)
{
    Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
    hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_DMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    method = Gtk::manage (new MyComboBoxText ());

    for (const auto method_string : RAWParams::XTransSensor::getMethodStrings()) {
        const std::string langKey =
            [method_string]() -> std::string
            {
                const std::string str(method_string);

                std::string res;
                for (const auto& c : str) {
                    switch (c) {
                        case '(':
                        case ')':
                        case ' ':
                        case '-': {
                            continue;
                        }

                        default: {
                            res += c;
                            break;
                        }
                    }
                }

                return res;
            }();
        method->append(M("TP_RAW_" + Glib::ustring(langKey).uppercase()));
    }

    method->set_active(0);
    hb1->set_tooltip_markup (M("TP_RAW_SENSOR_XTRANS_DMETHOD_TOOLTIP"));

    hb1->pack_end (*method, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *hb1, Gtk::PACK_SHRINK, 4);

    pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
    ccSteps = Gtk::manage (new Adjuster (M("TP_RAW_FALSECOLOR"), 0, 5, 1, 0 ));
    ccSteps->setAdjusterListener (this);

    if (ccSteps->delay < options.adjusterMaxDelay) {
        ccSteps->delay = options.adjusterMaxDelay;
    }

    ccSteps->show();
    pack_start( *ccSteps, Gtk::PACK_SHRINK, 4);

    methodconn = method->signal_changed().connect( sigc::mem_fun(*this, &XTransProcess::methodChanged) );
}


void XTransProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    methodconn.block (true);

    method->set_active(std::numeric_limits<int>::max());

    for (size_t i = 0; i < RAWParams::XTransSensor::getMethodStrings().size(); ++i)
        if( pp->raw.xtranssensor.method == RAWParams::XTransSensor::getMethodStrings()[i]) {
            method->set_active(i);
            oldSelection = i;
            break;
        }

    if(pedited ) {
        ccSteps->setEditedState (pedited->raw.xtranssensor.ccSteps ? Edited : UnEdited);

        if( !pedited->raw.xtranssensor.method ) {
            method->set_active(std::numeric_limits<int>::max()); // No name
        }
    }

    ccSteps->setValue (pp->raw.xtranssensor.ccSteps);

    methodconn.block (false);

    enableListener ();
}

void XTransProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.xtranssensor.ccSteps = ccSteps->getIntValue();

    int currentRow = method->get_active_row_number();

    if (currentRow >= 0 && currentRow < std::numeric_limits<int>::max()) {
        pp->raw.xtranssensor.method = procparams::RAWParams::XTransSensor::getMethodStrings()[currentRow];
    }

    if (pedited) {
        pedited->raw.xtranssensor.method = method->get_active_row_number() != std::numeric_limits<int>::max();
        pedited->raw.xtranssensor.ccSteps = ccSteps->getEditedState ();
    }
}

void XTransProcess::setBatchMode(bool batchMode)
{
    method->append (M("GENERAL_UNCHANGED"));
    method->set_active(std::numeric_limits<int>::max()); // No name
    ToolPanel::setBatchMode (batchMode);
    ccSteps->showEditedCB ();
}

void XTransProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    ccSteps->setDefault (defParams->raw.xtranssensor.ccSteps);

    if (pedited) {
        ccSteps->setDefaultEditedState(pedited->raw.xtranssensor.ccSteps ? Edited : UnEdited);
    } else {
        ccSteps->setDefaultEditedState(Irrelevant );
    }
}

void XTransProcess::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {
        if (a == ccSteps) {
            listener->panelChanged (EvDemosaicFalseColorIter, a->getTextValue() );
        }
    }
}

void XTransProcess::methodChanged ()
{
    const int curSelection = method->get_active_row_number();
    const RAWParams::XTransSensor::Method method = RAWParams::XTransSensor::Method(curSelection);

    Glib::ustring methodName;
    bool ppreq = false;

    if (curSelection >= 0 && curSelection < std::numeric_limits<int>::max()) {
        methodName = RAWParams::XTransSensor::getMethodStrings()[curSelection];

        if (method == RAWParams::XTransSensor::Method::MONO || RAWParams::XTransSensor::Method(oldSelection) == RAWParams::XTransSensor::Method::MONO) {
            ppreq = true;
        }
    }

    oldSelection = curSelection;

    if (listener) {
        listener->panelChanged (ppreq ? EvDemosaicMethodPreProc : EvDemosaicMethod, methodName);
    }
}
