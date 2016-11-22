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
        method->append_text(procparams::RAWParams::BayerSensor::methodstring[i]);
    }

    method->set_active(0);
    hb1->set_tooltip_markup (M("TP_RAW_DMETHOD_TOOLTIP"));

    hb1->pack_end (*method, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *hb1, Gtk::PACK_SHRINK, 4);

    imageNumberBox = Gtk::manage (new Gtk::HBox ());
    imageNumberBox->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_IMAGENUM") + ": ")), Gtk::PACK_SHRINK, 4);
    imageNumber = Gtk::manage (new MyComboBoxText ());
    imageNumber->append_text("1");
    imageNumber->append_text("2");
    imageNumber->append_text("3");
    imageNumber->append_text("4");
    imageNumber->set_active(0);
    imageNumberBox->set_tooltip_text(M("TP_RAW_IMAGENUM_TOOLTIP"));
    imageNumberBox->pack_end (*imageNumber, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *imageNumberBox, Gtk::PACK_SHRINK, 4);

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

    pixelShiftOptions = Gtk::manage (new Gtk::VBox ());
    pixelShiftOptions->set_border_width(4);
    pixelShiftShowMotion = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTSHOWMOTION")));
    pixelShiftOptions->pack_start(*pixelShiftShowMotion);

    pixelShiftAutomatic = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTADAPTIVE")));
    pixelShiftOptions->pack_start(*pixelShiftAutomatic);

    pixelShiftMotion = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTMOTION"), 0, 100, 1, 70));
    pixelShiftMotion->setAdjusterListener (this);
    pixelShiftMotion->set_tooltip_markup (M("TP_RAW_PIXELSHIFTMOTION_TOOLTIP"));

    if (pixelShiftMotion->delay < options.adjusterMaxDelay) {
        pixelShiftMotion->delay = options.adjusterMaxDelay;
    }
    pixelShiftMotion->show();
    pixelShiftOptions->pack_start(*pixelShiftMotion);

    pixelShiftMotionCorrection = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTMOTIONCORRECTION"), 1, 5, 2, 3));
    pixelShiftMotionCorrection->setAdjusterListener (this);
    pixelShiftMotionCorrection->set_tooltip_markup (M("TP_RAW_PIXELSHIFTMOTIONCORRECTION_TOOLTIP"));

    if (pixelShiftMotionCorrection->delay < options.adjusterMaxDelay) {
        pixelShiftMotionCorrection->delay = options.adjusterMaxDelay;
    }

    pixelShiftMotionCorrection->show();
    pixelShiftOptions->pack_start(*pixelShiftMotionCorrection);

    pixelShiftStddevFactor = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSTDDEVFACTOR"), 2, 8, 0.1, 5));
    pixelShiftStddevFactor->setAdjusterListener (this);
//    pixelShiftStddevFactor->set_tooltip_markup (M("TP_RAW_PIXELSHIFTSTDDEVFACTOR_TOOLTIP"));

    if (pixelShiftStddevFactor->delay < options.adjusterMaxDelay) {
        pixelShiftStddevFactor->delay = options.adjusterMaxDelay;
    }

    pixelShiftStddevFactor->show();
    pixelShiftOptions->pack_start(*pixelShiftStddevFactor);

    pixelShiftEperIso = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTEPERISO"), 0.3, 1.0, 0.05, 0.75));
    pixelShiftEperIso->setAdjusterListener (this);
//    pixelShiftStddevFactor->set_tooltip_markup (M("TP_RAW_PIXELSHIFTSTDDEVFACTOR_TOOLTIP"));

    if (pixelShiftEperIso->delay < options.adjusterMaxDelay) {
        pixelShiftEperIso->delay = options.adjusterMaxDelay;
    }

    pixelShiftEperIso->show();
    pixelShiftOptions->pack_start(*pixelShiftEperIso);

    pixelShiftNreadIso = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTNREADISO"), 1.0, 10.0, 0.05, 3.45));
    pixelShiftNreadIso->setAdjusterListener (this);
//    pixelShiftStddevFactor->set_tooltip_markup (M("TP_RAW_PIXELSHIFTSTDDEVFACTOR_TOOLTIP"));

    if (pixelShiftNreadIso->delay < options.adjusterMaxDelay) {
        pixelShiftNreadIso->delay = options.adjusterMaxDelay;
    }

    pixelShiftNreadIso->show();
    pixelShiftOptions->pack_start(*pixelShiftNreadIso);


    pixelShiftPrnu = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTPRNU"), 0.3, 2.0, 0.1, 1.0));
    pixelShiftPrnu->setAdjusterListener (this);
//    pixelShiftStddevFactor->set_tooltip_markup (M("TP_RAW_PIXELSHIFTSTDDEVFACTOR_TOOLTIP"));

    if (pixelShiftPrnu->delay < options.adjusterMaxDelay) {
        pixelShiftPrnu->delay = options.adjusterMaxDelay;
    }

    pixelShiftPrnu->show();
    pixelShiftOptions->pack_start(*pixelShiftPrnu);

    pack_start( *pixelShiftOptions, Gtk::PACK_SHRINK, 4);


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
    imagenumberconn = imageNumber->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::imageNumberChanged) );
    dcbEnhconn = dcbEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::dcbEnhanceChanged), true);
    pixelShiftShowMotionconn = pixelShiftShowMotion->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftShowMotionChanged), true);
    pixelShiftAutomaticconn = pixelShiftAutomatic->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftAutomaticChanged), true);
    //allEnhconn = allEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::allEnhanceChanged), true);
}


void BayerProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    methodconn.block (true);
    dcbEnhconn.block (true);
    imagenumberconn.block (true);
    //allEnhconn.block (true);

    method->set_active(procparams::RAWParams::BayerSensor::numMethods);
    imageNumber->set_active(pp->raw.bayersensor.imageNum);

    for( size_t i = 0; i < procparams::RAWParams::BayerSensor::numMethods; i++) {
        if( pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[i]) {
            method->set_active(i);
            oldMethod = i;
            break;
        }
    }

    if(pedited) {
        ccSteps->setEditedState (pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
        dcbIterations->setEditedState ( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        dcbEnhance->set_inconsistent(!pedited->raw.bayersensor.dcbEnhance);
        pixelShiftShowMotion->set_inconsistent(!pedited->raw.bayersensor.pixelshiftShowMotion);
        pixelShiftAutomatic->set_inconsistent(!pedited->raw.bayersensor.pixelShiftAutomatic);
        //allEnhance->set_inconsistent(!pedited->raw.bayersensor.allEnhance);
        lmmseIterations->setEditedState ( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        pixelShiftMotion->setEditedState ( pedited->raw.bayersensor.pixelshiftMotion ? Edited : UnEdited);
        pixelShiftMotionCorrection->setEditedState ( pedited->raw.bayersensor.pixelshiftMotionCorrection ? Edited : UnEdited);
        pixelShiftStddevFactor->setEditedState ( pedited->raw.bayersensor.pixelShiftStddevFactor ? Edited : UnEdited);
        pixelShiftEperIso->setEditedState ( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftNreadIso->setEditedState ( pedited->raw.bayersensor.pixelShiftNreadIso ? Edited : UnEdited);
        pixelShiftPrnu->setEditedState ( pedited->raw.bayersensor.pixelShiftPrnu ? Edited : UnEdited);

        if(!pedited->raw.bayersensor.method) {
            method->set_active(procparams::RAWParams::BayerSensor::numMethods);    // No name
        }
        if(!pedited->raw.bayersensor.imageNum) {
            imageNumber->set_active(4);
        }
    }

    //allEnhance->set_active(pp->raw.bayersensor.all_enhance);

    dcbIterations->setValue (pp->raw.bayersensor.dcb_iterations);
    dcbEnhance->set_active(pp->raw.bayersensor.dcb_enhance);
    pixelShiftShowMotion->set_active(pp->raw.bayersensor.pixelshiftShowMotion);
    pixelShiftAutomatic->set_active(pp->raw.bayersensor.pixelShiftAutomatic);
    ccSteps->setValue (pp->raw.bayersensor.ccSteps);
    lmmseIterations->setValue (pp->raw.bayersensor.lmmse_iterations);
    pixelShiftMotion->setValue (pp->raw.bayersensor.pixelshiftMotion);
    pixelShiftMotionCorrection->setValue (pp->raw.bayersensor.pixelshiftMotionCorrection);
    pixelShiftStddevFactor->setValue (pp->raw.bayersensor.pixelShiftStddevFactor);
    pixelShiftEperIso->setValue (pp->raw.bayersensor.pixelShiftEperIso);
    pixelShiftNreadIso->setValue (pp->raw.bayersensor.pixelShiftNreadIso);
    pixelShiftPrnu->setValue (pp->raw.bayersensor.pixelShiftPrnu);

    if (!batchMode) {
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::dcb] ||
                method->get_active_row_number() == procparams::RAWParams::BayerSensor::numMethods) {
            dcbOptions->show();
        } else {
            dcbOptions->hide();
        }
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::lmmse] ||
                method->get_active_row_number() == procparams::RAWParams::BayerSensor::numMethods) {
            lmmseOptions->show();
        } else {
            lmmseOptions->hide();
        }
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::pixelshift_simple] ||
                method->get_active_row_number() == procparams::RAWParams::BayerSensor::numMethods) {
            pixelShiftOptions->show();
        } else {
            pixelShiftOptions->hide();
        }

        // Flase color suppression is applied to all demozaicing method, so don't hide anything
        /*if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::eahd] ||
              pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::hphd] ||
              pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::vng4])
            ccSteps->show();
        else
            ccSteps->hide();*/
    }

    lastDCBen = pp->raw.bayersensor.dcb_enhance;
    //lastALLen = pp->raw.bayersensor.all_enhance;

    methodconn.block (false);
    imagenumberconn.block (false);
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
    pp->raw.bayersensor.pixelshiftMotion = pixelShiftMotion->getIntValue();
    pp->raw.bayersensor.pixelshiftMotionCorrection = pixelShiftMotionCorrection->getIntValue();
    pp->raw.bayersensor.pixelShiftStddevFactor = pixelShiftStddevFactor->getValue();
    pp->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getValue();
    pp->raw.bayersensor.pixelShiftNreadIso = pixelShiftNreadIso->getValue();
    pp->raw.bayersensor.pixelShiftPrnu = pixelShiftPrnu->getValue();
    pp->raw.bayersensor.pixelshiftShowMotion = pixelShiftShowMotion->get_active();
    pp->raw.bayersensor.pixelShiftAutomatic = pixelShiftAutomatic->get_active();

    int currentRow = method->get_active_row_number();
    if( currentRow >= 0 && currentRow < procparams::RAWParams::BayerSensor::numMethods) {
        pp->raw.bayersensor.method = procparams::RAWParams::BayerSensor::methodstring[currentRow];
    }

    currentRow = imageNumber->get_active_row_number();
    if (currentRow < 4) {
        pp->raw.bayersensor.imageNum = currentRow;
    }


    if (pedited) {
        pedited->raw.bayersensor.ccSteps = ccSteps->getEditedState ();
        pedited->raw.bayersensor.method = method->get_active_row_number() != procparams::RAWParams::BayerSensor::numMethods;
        pedited->raw.bayersensor.imageNum = imageNumber->get_active_row_number() < 4;
        pedited->raw.bayersensor.dcbIterations = dcbIterations->getEditedState ();
        pedited->raw.bayersensor.dcbEnhance = !dcbEnhance->get_inconsistent();
        //pedited->raw.bayersensor.allEnhance = !allEnhance->get_inconsistent();
        pedited->raw.bayersensor.lmmseIterations = lmmseIterations->getEditedState ();
        pedited->raw.bayersensor.pixelshiftMotion = pixelShiftMotion->getEditedState ();
        pedited->raw.bayersensor.pixelshiftMotionCorrection = pixelShiftMotionCorrection->getEditedState ();
        pedited->raw.bayersensor.pixelShiftStddevFactor = pixelShiftStddevFactor->getEditedState ();
        pedited->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftNreadIso = pixelShiftNreadIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftPrnu = pixelShiftPrnu->getEditedState ();
        pedited->raw.bayersensor.pixelshiftShowMotion = !pixelShiftShowMotion->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftAutomatic = !pixelShiftAutomatic->get_inconsistent();
    }
}

void BayerProcess::setBatchMode(bool batchMode)
{
    method->append_text (M("GENERAL_UNCHANGED"));
    method->set_active(procparams::RAWParams::BayerSensor::numMethods); // No name
    imageNumber->append_text (M("GENERAL_UNCHANGED"));
    imageNumber->set_active(4);
    dcbOptions->hide();
    lmmseOptions->hide();
    pixelShiftOptions->hide();
    ToolPanel::setBatchMode (batchMode);
    ccSteps->showEditedCB ();
    dcbIterations->showEditedCB ();
    lmmseIterations->showEditedCB ();
    pixelShiftMotion->showEditedCB ();
    pixelShiftMotionCorrection->showEditedCB ();
    pixelShiftStddevFactor->showEditedCB ();
    pixelShiftEperIso->showEditedCB ();
    pixelShiftNreadIso->showEditedCB ();
    pixelShiftPrnu->showEditedCB ();
}

void BayerProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    dcbIterations->setDefault( defParams->raw.bayersensor.dcb_iterations);
    lmmseIterations->setDefault( defParams->raw.bayersensor.lmmse_iterations);
    pixelShiftMotion->setDefault( defParams->raw.bayersensor.pixelshiftMotion);
    pixelShiftMotionCorrection->setDefault( defParams->raw.bayersensor.pixelshiftMotionCorrection);
    pixelShiftStddevFactor->setDefault( defParams->raw.bayersensor.pixelShiftStddevFactor);
    pixelShiftEperIso->setDefault( defParams->raw.bayersensor.pixelShiftEperIso);
    pixelShiftNreadIso->setDefault( defParams->raw.bayersensor.pixelShiftNreadIso);
    pixelShiftPrnu->setDefault( defParams->raw.bayersensor.pixelShiftPrnu);
    ccSteps->setDefault (defParams->raw.bayersensor.ccSteps);

    if (pedited) {
        dcbIterations->setDefaultEditedState( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        lmmseIterations->setDefaultEditedState( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        pixelShiftMotion->setDefaultEditedState( pedited->raw.bayersensor.pixelshiftMotion ? Edited : UnEdited);
        pixelShiftMotionCorrection->setDefaultEditedState( pedited->raw.bayersensor.pixelshiftMotionCorrection ? Edited : UnEdited);
        pixelShiftStddevFactor->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftStddevFactor ? Edited : UnEdited);
        pixelShiftEperIso->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftNreadIso->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftNreadIso ? Edited : UnEdited);
        pixelShiftPrnu->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftPrnu ? Edited : UnEdited);
        ccSteps->setDefaultEditedState(pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
    } else {
        dcbIterations->setDefaultEditedState( Irrelevant );
        lmmseIterations->setDefaultEditedState( Irrelevant );
        pixelShiftMotion->setDefaultEditedState( Irrelevant );
        pixelShiftMotionCorrection->setDefaultEditedState( Irrelevant );
        pixelShiftStddevFactor->setDefaultEditedState( Irrelevant );
        pixelShiftEperIso->setDefaultEditedState( Irrelevant );
        pixelShiftNreadIso->setDefaultEditedState( Irrelevant );
        pixelShiftPrnu->setDefaultEditedState( Irrelevant );
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
        } else if (a == pixelShiftMotion) {
            listener->panelChanged (EvDemosaicPixelshiftMotion, a->getTextValue() );
        } else if (a == pixelShiftMotionCorrection) {
            listener->panelChanged (EvDemosaicPixelshiftMotion, a->getTextValue() );
        } else if (a == pixelShiftStddevFactor) {
            listener->panelChanged (EvDemosaicPixelshiftMotion, a->getTextValue() );
        } else if (a == pixelShiftEperIso) {
            listener->panelChanged (EvDemosaicPixelshiftMotion, a->getTextValue() );
        } else if (a == pixelShiftNreadIso) {
            listener->panelChanged (EvDemosaicPixelshiftMotion, a->getTextValue() );
        } else if (a == pixelShiftPrnu) {
            listener->panelChanged (EvDemosaicPixelshiftMotion, a->getTextValue() );
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

    if ( curSelection == procparams::RAWParams::BayerSensor::pixelshift_simple) {
        pixelShiftOptions->show();
    } else {
        pixelShiftOptions->hide();
    }

    Glib::ustring methodName = "";
    bool ppreq = false;

    if( curSelection >= 0 && curSelection < procparams::RAWParams::BayerSensor::numMethods) {
        methodName = procparams::RAWParams::BayerSensor::methodstring[curSelection];

        if (curSelection == procparams::RAWParams::BayerSensor::mono || oldMethod == procparams::RAWParams::BayerSensor::mono) {
            ppreq = true;
        }
    }

    oldMethod = curSelection;

    if (listener) {
        listener->panelChanged (ppreq ? EvDemosaicMethodPreProc : EvDemosaicMethod, methodName);
    }
}

void BayerProcess::imageNumberChanged ()
{
    if (listener) {
        listener->panelChanged (EvRawImageNum, imageNumber->get_active_text());
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

void BayerProcess::pixelShiftShowMotionChanged ()
{
    if (batchMode) {
        if (pixelShiftShowMotion->get_inconsistent()) {
            pixelShiftShowMotion->set_inconsistent (false);
            pixelShiftShowMotionconn.block (true);
            pixelShiftShowMotion->set_active (false);
            pixelShiftShowMotionconn.block (false);
        } else if (lastDCBen) {
            pixelShiftShowMotion->set_inconsistent (true);
        }

        lastDCBen = pixelShiftShowMotion->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvDemosaicPixelshiftMotion, pixelShiftShowMotion->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftAutomaticChanged ()
{
    if (batchMode) {
        if (pixelShiftAutomatic->get_inconsistent()) {
            pixelShiftAutomatic->set_inconsistent (false);
            pixelShiftAutomaticconn.block (true);
            pixelShiftAutomatic->set_active (false);
            pixelShiftAutomaticconn.block (false);
        } else if (lastDCBen) {
            pixelShiftAutomatic->set_inconsistent (true);
        }

        lastDCBen = pixelShiftAutomatic->get_active ();
    }
    pixelShiftMotion->set_sensitive(!pixelShiftAutomatic->get_active ());
    pixelShiftEperIso->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftNreadIso->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftPrnu->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftStddevFactor->set_sensitive(pixelShiftAutomatic->get_active ());

    if (listener) {
        listener->panelChanged (EvDemosaicPixelshiftMotion, pixelShiftAutomatic->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
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
