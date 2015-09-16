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
#include "bayerprocess.h"
#include "options.h"
#include "guiutils.h"
using namespace rtengine;
using namespace rtengine::procparams;

BayerProcess::BayerProcess () : FoldableToolPanel(this, "bayerprocess", M("TP_RAW_LABEL"), true)
{
    Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
    hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_DMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    method = Gtk::manage (new MyComboBoxText ());

    for( size_t i = 0; i < procparams::RAWParams::BayerSensor::numMethods; i++) {
        method->append(procparams::RAWParams::BayerSensor::methodstring[i]);
    }

    method->set_active(0);
    hb1->set_tooltip_markup (M("TP_RAW_DMETHOD_TOOLTIP"));

    hb1->pack_end (*method, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *hb1, Gtk::PACK_SHRINK, 4);

    dcbOptions = Gtk::manage (new Gtk::VBox ());
    dcbOptions->set_border_width(4);

    dcbIterations = Gtk::manage (new Adjuster (M("TP_RAW_DCBITERATIONS"), 0, 5, 1, 2));
    dcbIterations->setAdjusterListener (this);

    if (dcbIterations->delay < options.adjusterMaxDelay) {
        dcbIterations->delay = options.adjusterMaxDelay;
    }

    dcbIterations->show();
    dcbEnhance = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_DCBENHANCE")));
    dcbOptions->pack_start(*dcbIterations);
    dcbOptions->pack_start(*dcbEnhance);
    pack_start( *dcbOptions, Gtk::PACK_SHRINK, 4);

    lmmseOptions = Gtk::manage (new Gtk::VBox ());
    lmmseOptions->set_border_width(4);

    lmmseIterations = Gtk::manage (new Adjuster (M("TP_RAW_LMMSEITERATIONS"), 0, 6, 1, 2));
    lmmseIterations->setAdjusterListener (this);
    lmmseIterations->set_tooltip_markup (M("TP_RAW_LMMSE_TOOLTIP"));

    if (lmmseIterations->delay < options.adjusterMaxDelay) {
        lmmseIterations->delay = options.adjusterMaxDelay;
    }

    lmmseIterations->show();
    lmmseOptions->pack_start(*lmmseIterations);
    pack_start( *lmmseOptions, Gtk::PACK_SHRINK, 4);

    pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
    ccSteps = Gtk::manage (new Adjuster (M("TP_RAW_FALSECOLOR"), 0, 5, 1, 0 ));
    ccSteps->setAdjusterListener (this);

    if (ccSteps->delay < options.adjusterMaxDelay) {
        ccSteps->delay = options.adjusterMaxDelay;
    }

    ccSteps->show();
    pack_start( *ccSteps, Gtk::PACK_SHRINK, 4);

    //pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
    //allOptions = Gtk::manage (new Gtk::VBox ());
    //allOptions->set_border_width(2);
    //allEnhance = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_ALLENHANCE")));
    //allOptions->pack_start(*allEnhance);
    //pack_start( *allOptions, Gtk::PACK_SHRINK, 4);

    methodconn = method->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::methodChanged) );
    dcbEnhconn = dcbEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::dcbEnhanceChanged), true);
    //allEnhconn = allEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::allEnhanceChanged), true);
}


void BayerProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    methodconn.block (true);
    dcbEnhconn.block (true);
    //allEnhconn.block (true);

    method->set_active(procparams::RAWParams::BayerSensor::numMethods);

    for( size_t i = 0; i < procparams::RAWParams::BayerSensor::numMethods; i++)
        if( pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[i]) {
            method->set_active(i);
            oldSelection = i;
            break;
        }

    if(pedited ) {
        ccSteps->setEditedState (pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
        dcbIterations->setEditedState ( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        dcbEnhance->set_inconsistent(!pedited->raw.bayersensor.dcbEnhance);
        //allEnhance->set_inconsistent(!pedited->raw.bayersensor.allEnhance);
        lmmseIterations->setEditedState ( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);

        if( !pedited->raw.bayersensor.method ) {
            method->set_active(procparams::RAWParams::BayerSensor::numMethods);    // No name
        }
    }

    //allEnhance->set_active(pp->raw.bayersensor.all_enhance);

    dcbIterations->setValue (pp->raw.bayersensor.dcb_iterations);
    dcbEnhance->set_active(pp->raw.bayersensor.dcb_enhance);
    ccSteps->setValue (pp->raw.bayersensor.ccSteps);

    if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::dcb] ||
            method->get_active_row_number() == procparams::RAWParams::BayerSensor::numMethods) {
        dcbOptions->show();
    } else {
        dcbOptions->hide();
    }

    lmmseIterations->setValue (pp->raw.bayersensor.lmmse_iterations);

    if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::lmmse] ||
            method->get_active_row_number() == procparams::RAWParams::BayerSensor::numMethods) {
        lmmseOptions->show();
    } else {
        lmmseOptions->hide();
    }

    // Flase color suppression is applied to all demozaicing method, so don't hide anything
    /*if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::eahd] ||
          pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::hphd] ||
          pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::vng4])
        ccSteps->show();
    else
        ccSteps->hide();*/

    lastDCBen = pp->raw.bayersensor.dcb_enhance;
    //lastALLen = pp->raw.bayersensor.all_enhance;

    methodconn.block (false);
    dcbEnhconn.block (false);
    //allEnhconn.block (false);

    enableListener ();
}

void BayerProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.bayersensor.ccSteps = ccSteps->getIntValue();
    pp->raw.bayersensor.dcb_iterations = dcbIterations->getIntValue();
    pp->raw.bayersensor.dcb_enhance = dcbEnhance->get_active();
    //pp->raw.bayersensor.all_enhance = allEnhance->get_active();
    pp->raw.bayersensor.lmmse_iterations = lmmseIterations->getIntValue();

    int currentRow = method->get_active_row_number();

    if( currentRow >= 0 && currentRow < procparams::RAWParams::BayerSensor::numMethods) {
        pp->raw.bayersensor.method = procparams::RAWParams::BayerSensor::methodstring[currentRow];
    }

    if (pedited) {
        pedited->raw.bayersensor.ccSteps = ccSteps->getEditedState ();
        pedited->raw.bayersensor.method = method->get_active_row_number() != procparams::RAWParams::BayerSensor::numMethods;
        pedited->raw.bayersensor.dcbIterations = dcbIterations->getEditedState ();
        pedited->raw.bayersensor.dcbEnhance = !dcbEnhance->get_inconsistent();
        //pedited->raw.bayersensor.allEnhance = !allEnhance->get_inconsistent();
        pedited->raw.bayersensor.lmmseIterations = lmmseIterations->getEditedState ();

    }
}

void BayerProcess::setBatchMode(bool batchMode)
{
    method->append (M("GENERAL_UNCHANGED"));
    method->set_active(procparams::RAWParams::BayerSensor::numMethods); // No name
    dcbOptions->hide();
    lmmseOptions->hide();
    ToolPanel::setBatchMode (batchMode);
    ccSteps->showEditedCB ();
    dcbIterations->showEditedCB ();
    lmmseIterations->showEditedCB ();
}

void BayerProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    dcbIterations->setDefault( defParams->raw.bayersensor.dcb_iterations);
    lmmseIterations->setDefault( defParams->raw.bayersensor.lmmse_iterations);
    ccSteps->setDefault (defParams->raw.bayersensor.ccSteps);

    if (pedited) {
        dcbIterations->setDefaultEditedState( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        lmmseIterations->setDefaultEditedState( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        ccSteps->setDefaultEditedState(pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
    } else {
        dcbIterations->setDefaultEditedState( Irrelevant );
        lmmseIterations->setDefaultEditedState( Irrelevant );
        ccSteps->setDefaultEditedState(Irrelevant );
    }
}

void BayerProcess::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {
        if (a == dcbIterations) {
            listener->panelChanged (EvDemosaicDCBIter, a->getTextValue() );
        } else if (a == ccSteps) {
            listener->panelChanged (EvDemosaicFalseColorIter, a->getTextValue() );
        } else if (a == lmmseIterations) {
            listener->panelChanged (EvDemosaicLMMSEIter, a->getTextValue() );
        }
    }
}

void BayerProcess::methodChanged ()
{
    int  curSelection = method->get_active_row_number();

    if ( curSelection == procparams::RAWParams::BayerSensor::dcb) {
        dcbOptions->show();
    } else {
        dcbOptions->hide();
    }

    if ( curSelection == procparams::RAWParams::BayerSensor::lmmse) {
        lmmseOptions->show();
    } else {
        lmmseOptions->hide();
    }

    Glib::ustring methodName = "";
    bool ppreq = false;

    if( curSelection >= 0 && curSelection < procparams::RAWParams::BayerSensor::numMethods) {
        methodName = procparams::RAWParams::BayerSensor::methodstring[curSelection];

        if (curSelection == procparams::RAWParams::BayerSensor::mono || oldSelection == procparams::RAWParams::BayerSensor::mono) {
            ppreq = true;
        }
    }

    oldSelection = curSelection;

    if (listener) {
        listener->panelChanged (ppreq ? EvDemosaicMethodPreProc : EvDemosaicMethod, methodName);
    }
}

void BayerProcess::dcbEnhanceChanged ()
{
    if (batchMode) {
        if (dcbEnhance->get_inconsistent()) {
            dcbEnhance->set_inconsistent (false);
            dcbEnhconn.block (true);
            dcbEnhance->set_active (false);
            dcbEnhconn.block (false);
        } else if (lastDCBen) {
            dcbEnhance->set_inconsistent (true);
        }

        lastDCBen = dcbEnhance->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvDemosaicDCBEnhanced, dcbEnhance->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

/*void BayerProcess::allEnhanceChanged ()
{
    if (batchMode) {
        if (allEnhance->get_inconsistent()) {
            allEnhance->set_inconsistent (false);
            allEnhconn.block (true);
            allEnhance->set_active (false);
            allEnhconn.block (false);
        }
        else if (lastALLen)
            allEnhance->set_inconsistent (true);

        lastALLen = allEnhance->get_active ();
    }
    if (listener)
        listener->panelChanged (EvDemosaicALLEnhanced, allEnhance->get_active()?M("GENERAL_ENABLED"):M("GENERAL_DISABLED"));
}*/
