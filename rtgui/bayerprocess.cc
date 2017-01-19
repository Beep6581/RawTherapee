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

    pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
    ccSteps = Gtk::manage (new Adjuster (M("TP_RAW_FALSECOLOR"), 0, 5, 1, 0 ));
    ccSteps->setAdjusterListener (this);

    if (ccSteps->delay < options.adjusterMaxDelay) {
        ccSteps->delay = options.adjusterMaxDelay;
    }

    ccSteps->show();
    pack_start( *ccSteps, Gtk::PACK_SHRINK, 4);


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

    pixelShiftFrame = Gtk::manage (new Gtk::VBox ());
    pixelShiftFrame->set_border_width(0);

    Gtk::HBox* hb3 = Gtk::manage (new Gtk::HBox ());
    hb3->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_PIXELSHIFTMOTIONMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    pixelShiftMotionMethod = Gtk::manage (new MyComboBoxText ());
    pixelShiftMotionMethod->append_text("Off");
    pixelShiftMotionMethod->append_text("Automatic");
    pixelShiftMotionMethod->append_text("Custom");
    pixelShiftMotionMethod->set_active(1);
    pixelShiftMotionMethod->show();
    hb3->pack_start(*pixelShiftMotionMethod);
    pixelShiftFrame->pack_start(*hb3);

    pixelShiftOptions = Gtk::manage (new Gtk::VBox ());
    pixelShiftOptions->set_border_width(0);

    pixelShiftShowMotion = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTSHOWMOTION")));
    pixelShiftShowMotion->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTION_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftShowMotion);

    pixelShiftShowMotionMaskOnly = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY")));
    pixelShiftShowMotionMaskOnly->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftShowMotionMaskOnly);

    pixelShiftAutomatic = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTADAPTIVE")));
    pixelShiftOptions->pack_start(*pixelShiftAutomatic);

    pixelShiftGreen = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTGREEN")));
    pixelShiftOptions->pack_start(*pixelShiftGreen);

    pixelShiftBlur = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTBLUR")));
    pixelShiftBlur->set_tooltip_text (M("TP_RAW_PIXELSHIFTBLUR_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftBlur);

    pixelShiftSmooth = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTSMOOTH")));
    pixelShiftSmooth->set_tooltip_text (M("TP_RAW_PIXELSHIFTSMOOTH_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftSmooth);

    pixelShiftHoleFill = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTHOLEFILL")));
    pixelShiftHoleFill->set_tooltip_text (M("TP_RAW_PIXELSHIFTHOLEFILL_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftHoleFill);

    pixelShiftMedian = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTMEDIAN")));
    pixelShiftMedian->set_tooltip_text (M("TP_RAW_PIXELSHIFTMEDIAN_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftMedian);

    pixelShiftMedian3 = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTMEDIAN3")));
    pixelShiftMedian3->set_tooltip_text (M("TP_RAW_PIXELSHIFTMEDIAN3_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftMedian3);

    pixelShiftNonGreenCross = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTNONGREENCROSS")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenCross);

    pixelShiftNonGreenCross2 = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTNONGREENCROSS2")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenCross2);

    pixelShiftNonGreenAmaze = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTNONGREENAMAZE")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenAmaze);

    pixelShiftNonGreenHorizontal = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTNONGREENHORIZONTAL")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenHorizontal);

    pixelShiftNonGreenVertical = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTNONGREENVERTICAL")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenVertical);

    pixelShiftExp0 = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_PIXELSHIFTEXP0")));
    pixelShiftOptions->pack_start(*pixelShiftExp0);

    pixelShiftMotion = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTMOTION"), 0, 100, 1, 70));
    pixelShiftMotion->setAdjusterListener (this);
    pixelShiftMotion->set_tooltip_text (M("TP_RAW_PIXELSHIFTMOTION_TOOLTIP"));

    if (pixelShiftMotion->delay < options.adjusterMaxDelay) {
        pixelShiftMotion->delay = options.adjusterMaxDelay;
    }
    pixelShiftMotion->show();
    pixelShiftOptions->pack_start(*pixelShiftMotion);


    Gtk::HBox* hb2 = Gtk::manage (new Gtk::HBox ());
    hb2->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_PIXELSHIFTMOTIONCORRECTION") + ": ")), Gtk::PACK_SHRINK, 0);
    pixelShiftMotionCorrection = Gtk::manage (new MyComboBoxText ());
    pixelShiftMotionCorrection->append_text("1x1");
    pixelShiftMotionCorrection->append_text("1x2");
    pixelShiftMotionCorrection->append_text("3x3");
    pixelShiftMotionCorrection->append_text("5x5");
    pixelShiftMotionCorrection->append_text("7x7");
    pixelShiftMotionCorrection->append_text("3x3 new");
    pixelShiftMotionCorrection->set_active(0);
    pixelShiftMotionCorrection->show();
    hb2->pack_start(*pixelShiftMotionCorrection);
    pixelShiftOptions->pack_start(*hb2);

    pixelShiftStddevFactorGreen = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSTDDEVFACTORGREEN"), 2, 8, 0.1, 5));
    pixelShiftStddevFactorGreen->setAdjusterListener (this);

    if (pixelShiftStddevFactorGreen->delay < options.adjusterMaxDelay) {
        pixelShiftStddevFactorGreen->delay = options.adjusterMaxDelay;
    }

    pixelShiftStddevFactorGreen->show();
    pixelShiftOptions->pack_start(*pixelShiftStddevFactorGreen);

    pixelShiftStddevFactorRed = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSTDDEVFACTORRED"), 1, 8, 0.1, 5));
    pixelShiftStddevFactorRed->setAdjusterListener (this);

    if (pixelShiftStddevFactorRed->delay < options.adjusterMaxDelay) {
        pixelShiftStddevFactorRed->delay = options.adjusterMaxDelay;
    }

    pixelShiftStddevFactorRed->show();
    pixelShiftOptions->pack_start(*pixelShiftStddevFactorRed);

    pixelShiftStddevFactorBlue = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSTDDEVFACTORBLUE"), 1, 8, 0.1, 5));
    pixelShiftStddevFactorBlue->setAdjusterListener (this);

    if (pixelShiftStddevFactorBlue->delay < options.adjusterMaxDelay) {
        pixelShiftStddevFactorBlue->delay = options.adjusterMaxDelay;
    }

    pixelShiftStddevFactorBlue->show();
    pixelShiftOptions->pack_start(*pixelShiftStddevFactorBlue);

    pixelShiftEperIso = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTEPERISO"), -2.0, 2.0, 0.05, 0.0));
    pixelShiftEperIso->setAdjusterListener (this);

    if (pixelShiftEperIso->delay < options.adjusterMaxDelay) {
        pixelShiftEperIso->delay = options.adjusterMaxDelay;
    }

    pixelShiftEperIso->show();
    pixelShiftOptions->pack_start(*pixelShiftEperIso);

    pixelShiftNreadIso = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTNREADISO"), -2.0, 2.0, 0.05, 0.0));
    pixelShiftNreadIso->setAdjusterListener (this);

    if (pixelShiftNreadIso->delay < options.adjusterMaxDelay) {
        pixelShiftNreadIso->delay = options.adjusterMaxDelay;
    }

    pixelShiftNreadIso->show();
    pixelShiftOptions->pack_start(*pixelShiftNreadIso);


    pixelShiftPrnu = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTPRNU"), 0.3, 2.0, 0.1, 1.0));
    pixelShiftPrnu->setAdjusterListener (this);

    if (pixelShiftPrnu->delay < options.adjusterMaxDelay) {
        pixelShiftPrnu->delay = options.adjusterMaxDelay;
    }

    pixelShiftPrnu->show();
    pixelShiftOptions->pack_start(*pixelShiftPrnu);

    pixelShiftSigma = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSIGMA"), 0.5, 25, 0.1, 1.0));
    pixelShiftSigma->setAdjusterListener (this);

    if (pixelShiftSigma->delay < options.adjusterMaxDelay) {
        pixelShiftSigma->delay = options.adjusterMaxDelay;
    }

    pixelShiftSigma->show();
    pixelShiftOptions->pack_start(*pixelShiftSigma);

    pixelShiftSum = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTMASKTHRESHOLD"), 1.0, 8.0, 0.1, 3.0));
    pixelShiftSum->setAdjusterListener (this);

    if (pixelShiftSum->delay < options.adjusterMaxDelay) {
        pixelShiftSum->delay = options.adjusterMaxDelay;
    }

    pixelShiftSum->show();
    pixelShiftOptions->pack_start(*pixelShiftSum);

    pixelShiftRedBlueWeight = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTREDBLUEWEIGHT"), 0.1, 1.0, 0.1, 0.7));
    pixelShiftRedBlueWeight->setAdjusterListener (this);

    if (pixelShiftRedBlueWeight->delay < options.adjusterMaxDelay) {
        pixelShiftRedBlueWeight->delay = options.adjusterMaxDelay;
    }

    pixelShiftRedBlueWeight->show();
    pixelShiftOptions->pack_start(*pixelShiftRedBlueWeight);

    pixelShiftFrame->pack_start(*pixelShiftOptions);
    pixelShiftOptions->hide();

    pack_start( *pixelShiftFrame, Gtk::PACK_SHRINK, 4);

    //pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
    //allOptions = Gtk::manage (new Gtk::VBox ());
    //allOptions->set_border_width(2);
    //allEnhance = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_ALLENHANCE")));
    //allOptions->pack_start(*allEnhance);
    //pack_start( *allOptions, Gtk::PACK_SHRINK, 4);

    methodconn = method->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::methodChanged) );
    psmcconn = pixelShiftMotionCorrection->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::psMotionCorrectionChanged) );
    imagenumberconn = imageNumber->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::imageNumberChanged) );
    dcbEnhconn = dcbEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::dcbEnhanceChanged), true);
    pixelShiftMotionMethodConn = pixelShiftMotionMethod->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::pixelShiftMotionMethodChanged) );
    pixelShiftShowMotionconn = pixelShiftShowMotion->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftShowMotionChanged), true);
    pixelShiftShowMotionMaskOnlyconn = pixelShiftShowMotionMaskOnly->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftShowMotionMaskOnlyChanged), true);
    pixelShiftAutomaticconn = pixelShiftAutomatic->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftAutomaticChanged), true);
    pixelShiftNonGreenHorizontalconn = pixelShiftNonGreenHorizontal->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenHorizontalChanged), true);
    pixelShiftNonGreenVerticalconn = pixelShiftNonGreenVertical->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenVerticalChanged), true);
    pixelShiftHoleFillconn = pixelShiftHoleFill->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftHoleFillChanged), true);
    pixelShiftMedianconn = pixelShiftMedian->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftMedianChanged), true);
    pixelShiftMedian3conn = pixelShiftMedian3->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftMedian3Changed), true);
    pixelShiftGreenconn = pixelShiftGreen->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftGreenChanged), true);
    pixelShiftBlurconn = pixelShiftBlur->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftBlurChanged), true);
    pixelShiftSmoothconn = pixelShiftSmooth->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftSmoothChanged), true);
    pixelShiftExp0conn = pixelShiftExp0->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftExp0Changed), true);
    pixelShiftNonGreenCrossconn = pixelShiftNonGreenCross->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenCrossChanged), true);
    pixelShiftNonGreenCross2conn = pixelShiftNonGreenCross2->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenCross2Changed), true);
    pixelShiftNonGreenAmazeconn = pixelShiftNonGreenAmaze->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenAmazeChanged), true);
}


void BayerProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    methodconn.block (true);
    dcbEnhconn.block (true);
    imagenumberconn.block (true);
    psmcconn.block (true);
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
        pixelShiftShowMotionMaskOnly->set_inconsistent(!pedited->raw.bayersensor.pixelshiftShowMotionMaskOnly);
        pixelShiftAutomatic->set_inconsistent(!pedited->raw.bayersensor.pixelShiftAutomatic);
        pixelShiftNonGreenHorizontal->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenHorizontal);
        pixelShiftNonGreenVertical->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenVertical);
        pixelShiftHoleFill->set_inconsistent(!pedited->raw.bayersensor.pixelShiftHoleFill);
        pixelShiftMedian->set_inconsistent(!pedited->raw.bayersensor.pixelShiftMedian);
        pixelShiftMedian3->set_inconsistent(!pedited->raw.bayersensor.pixelShiftMedian3);
        pixelShiftGreen->set_inconsistent(!pedited->raw.bayersensor.pixelShiftGreen);
        pixelShiftBlur->set_inconsistent(!pedited->raw.bayersensor.pixelShiftBlur);
        pixelShiftSmooth->set_inconsistent(!pedited->raw.bayersensor.pixelShiftSmooth);
        pixelShiftExp0->set_inconsistent(!pedited->raw.bayersensor.pixelShiftExp0);
        pixelShiftNonGreenCross->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenCross);
        pixelShiftNonGreenCross2->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenCross2);
        pixelShiftNonGreenAmaze->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenAmaze);
        lmmseIterations->setEditedState ( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        pixelShiftMotion->setEditedState ( pedited->raw.bayersensor.pixelShiftMotion ? Edited : UnEdited);
        pixelShiftStddevFactorGreen->setEditedState ( pedited->raw.bayersensor.pixelShiftStddevFactorGreen ? Edited : UnEdited);
        pixelShiftStddevFactorRed->setEditedState ( pedited->raw.bayersensor.pixelShiftStddevFactorRed ? Edited : UnEdited);
        pixelShiftStddevFactorBlue->setEditedState ( pedited->raw.bayersensor.pixelShiftStddevFactorBlue ? Edited : UnEdited);
        pixelShiftEperIso->setEditedState ( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftNreadIso->setEditedState ( pedited->raw.bayersensor.pixelShiftNreadIso ? Edited : UnEdited);
        pixelShiftPrnu->setEditedState ( pedited->raw.bayersensor.pixelShiftPrnu ? Edited : UnEdited);
        pixelShiftSigma->setEditedState ( pedited->raw.bayersensor.pixelShiftSigma ? Edited : UnEdited);
        pixelShiftSum->setEditedState ( pedited->raw.bayersensor.pixelShiftSum ? Edited : UnEdited);
        pixelShiftRedBlueWeight->setEditedState ( pedited->raw.bayersensor.pixelShiftRedBlueWeight ? Edited : UnEdited);

        if(!pedited->raw.bayersensor.method) {
            method->set_active(procparams::RAWParams::BayerSensor::numMethods);    // No name
        }
        if(!pedited->raw.bayersensor.imageNum) {
            imageNumber->set_active_text(M("GENERAL_UNCHANGED"));
        }
        if(!pedited->raw.bayersensor.pixelShiftMotionCorrection) {
            pixelShiftMotionCorrection->set_active_text(M("GENERAL_UNCHANGED"));
        }
        if(!pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod) {
            pixelShiftMotionMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }
    }

    //allEnhance->set_active(pp->raw.bayersensor.all_enhance);

    dcbIterations->setValue (pp->raw.bayersensor.dcb_iterations);
    dcbEnhance->set_active(pp->raw.bayersensor.dcb_enhance);
    pixelShiftShowMotion->set_active(pp->raw.bayersensor.pixelshiftShowMotion);
    pixelShiftShowMotionMaskOnly->set_active(pp->raw.bayersensor.pixelshiftShowMotionMaskOnly);
    pixelShiftAutomatic->set_active(pp->raw.bayersensor.pixelShiftAutomatic);
    pixelShiftNonGreenHorizontal->set_active(pp->raw.bayersensor.pixelShiftNonGreenHorizontal);
    pixelShiftNonGreenVertical->set_active(pp->raw.bayersensor.pixelShiftNonGreenVertical);
    pixelShiftHoleFill->set_active(pp->raw.bayersensor.pixelShiftHoleFill);
    pixelShiftMedian->set_active(pp->raw.bayersensor.pixelShiftMedian);
    pixelShiftMedian3->set_active(pp->raw.bayersensor.pixelShiftMedian3);
    pixelShiftGreen->set_active(pp->raw.bayersensor.pixelShiftGreen);
    pixelShiftBlur->set_active(pp->raw.bayersensor.pixelShiftBlur);
    pixelShiftSmooth->set_active(pp->raw.bayersensor.pixelShiftSmooth);
    pixelShiftExp0->set_active(pp->raw.bayersensor.pixelShiftExp0);
    pixelShiftNonGreenCross->set_active(pp->raw.bayersensor.pixelShiftNonGreenCross);
    pixelShiftNonGreenCross2->set_active(pp->raw.bayersensor.pixelShiftNonGreenCross2);
    pixelShiftNonGreenAmaze->set_active(pp->raw.bayersensor.pixelShiftNonGreenAmaze);
    ccSteps->setValue (pp->raw.bayersensor.ccSteps);
    lmmseIterations->setValue (pp->raw.bayersensor.lmmse_iterations);
    pixelShiftMotion->setValue (pp->raw.bayersensor.pixelShiftMotion);
    pixelShiftMotionCorrection->set_active ((int)pp->raw.bayersensor.pixelShiftMotionCorrection);
    pixelShiftMotionMethod->set_active ((int)pp->raw.bayersensor.pixelShiftMotionCorrectionMethod);
    pixelShiftStddevFactorGreen->setValue (pp->raw.bayersensor.pixelShiftStddevFactorGreen);
    pixelShiftStddevFactorRed->setValue (pp->raw.bayersensor.pixelShiftStddevFactorRed);
    pixelShiftStddevFactorBlue->setValue (pp->raw.bayersensor.pixelShiftStddevFactorBlue);
    pixelShiftEperIso->setValue (pp->raw.bayersensor.pixelShiftEperIso);
    pixelShiftNreadIso->setValue (pp->raw.bayersensor.pixelShiftNreadIso);
    pixelShiftPrnu->setValue (pp->raw.bayersensor.pixelShiftPrnu);
    pixelShiftSigma->setValue (pp->raw.bayersensor.pixelShiftSigma);
    pixelShiftSum->setValue (pp->raw.bayersensor.pixelShiftSum);
    pixelShiftRedBlueWeight->setValue (pp->raw.bayersensor.pixelShiftRedBlueWeight);

    pixelShiftHoleFill->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5);
    pixelShiftBlur->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5);
    pixelShiftSmooth->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5 && pixelShiftBlur->get_active());

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
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::pixelshift] ||
                method->get_active_row_number() == procparams::RAWParams::BayerSensor::numMethods) {
            if(pp->raw.bayersensor.pixelShiftMotionCorrectionMethod == RAWParams::BayerSensor::Custom) {
                pixelShiftOptions->show();
            } else {
                pixelShiftOptions->hide();
            }
            pixelShiftFrame->show();
        } else {
            pixelShiftFrame->hide();
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
    psmcconn.block (false);
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
    pp->raw.bayersensor.pixelShiftMotion = pixelShiftMotion->getIntValue();
    pp->raw.bayersensor.pixelShiftMotionCorrection = (RAWParams::BayerSensor::ePSMotionCorrection)pixelShiftMotionCorrection->get_active_row_number();
    pp->raw.bayersensor.pixelShiftMotionCorrectionMethod = (RAWParams::BayerSensor::ePSMotionCorrectionMethod)pixelShiftMotionMethod->get_active_row_number();
    pp->raw.bayersensor.pixelShiftStddevFactorGreen = pixelShiftStddevFactorGreen->getValue();
    pp->raw.bayersensor.pixelShiftStddevFactorRed = pixelShiftStddevFactorRed->getValue();
    pp->raw.bayersensor.pixelShiftStddevFactorBlue = pixelShiftStddevFactorBlue->getValue();
    pp->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getValue();
    pp->raw.bayersensor.pixelShiftNreadIso = pixelShiftNreadIso->getValue();
    pp->raw.bayersensor.pixelShiftPrnu = pixelShiftPrnu->getValue();
    pp->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getValue();
    pp->raw.bayersensor.pixelShiftSum = pixelShiftSum->getValue();
    pp->raw.bayersensor.pixelShiftRedBlueWeight = pixelShiftRedBlueWeight->getValue();
    pp->raw.bayersensor.pixelshiftShowMotion = pixelShiftShowMotion->get_active();
    pp->raw.bayersensor.pixelshiftShowMotionMaskOnly = pixelShiftShowMotionMaskOnly->get_active();
    pp->raw.bayersensor.pixelShiftAutomatic = pixelShiftAutomatic->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenHorizontal = pixelShiftNonGreenHorizontal->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenVertical = pixelShiftNonGreenVertical->get_active();
    pp->raw.bayersensor.pixelShiftHoleFill = pixelShiftHoleFill->get_active();
    pp->raw.bayersensor.pixelShiftMedian = pixelShiftMedian->get_active();
    pp->raw.bayersensor.pixelShiftMedian3 = pixelShiftMedian3->get_active();
    pp->raw.bayersensor.pixelShiftGreen = pixelShiftGreen->get_active();
    pp->raw.bayersensor.pixelShiftBlur = pixelShiftBlur->get_active();
    pp->raw.bayersensor.pixelShiftSmooth = pixelShiftSmooth->get_active();
    pp->raw.bayersensor.pixelShiftExp0 = pixelShiftExp0->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenCross = pixelShiftNonGreenCross->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenCross2 = pixelShiftNonGreenCross2->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenAmaze = pixelShiftNonGreenAmaze->get_active();

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
        pedited->raw.bayersensor.imageNum = imageNumber->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.dcbIterations = dcbIterations->getEditedState ();
        pedited->raw.bayersensor.dcbEnhance = !dcbEnhance->get_inconsistent();
        //pedited->raw.bayersensor.allEnhance = !allEnhance->get_inconsistent();
        pedited->raw.bayersensor.lmmseIterations = lmmseIterations->getEditedState ();
        pedited->raw.bayersensor.pixelShiftMotion = pixelShiftMotion->getEditedState ();
        pedited->raw.bayersensor.pixelShiftMotionCorrection = pixelShiftMotionCorrection->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod = pixelShiftMotionMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.pixelShiftStddevFactorGreen = pixelShiftStddevFactorGreen->getEditedState ();
        pedited->raw.bayersensor.pixelShiftStddevFactorRed = pixelShiftStddevFactorRed->getEditedState ();
        pedited->raw.bayersensor.pixelShiftStddevFactorBlue = pixelShiftStddevFactorBlue->getEditedState ();
        pedited->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftNreadIso = pixelShiftNreadIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftPrnu = pixelShiftPrnu->getEditedState ();
        pedited->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getEditedState ();
        pedited->raw.bayersensor.pixelShiftSum = pixelShiftSum->getEditedState ();
        pedited->raw.bayersensor.pixelShiftRedBlueWeight = pixelShiftRedBlueWeight->getEditedState ();
        pedited->raw.bayersensor.pixelshiftShowMotion = !pixelShiftShowMotion->get_inconsistent();
        pedited->raw.bayersensor.pixelshiftShowMotionMaskOnly = !pixelShiftShowMotionMaskOnly->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftAutomatic = !pixelShiftAutomatic->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenHorizontal = !pixelShiftNonGreenHorizontal->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenVertical = !pixelShiftNonGreenVertical->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftHoleFill = !pixelShiftHoleFill->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftMedian = !pixelShiftMedian->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftMedian3 = !pixelShiftMedian3->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftGreen = !pixelShiftGreen->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftBlur = !pixelShiftBlur->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftSmooth = !pixelShiftSmooth->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftExp0 = !pixelShiftExp0->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenCross = !pixelShiftNonGreenCross->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenCross2 = !pixelShiftNonGreenCross2->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenAmaze = !pixelShiftNonGreenAmaze->get_inconsistent();
    }
}

void BayerProcess::setBatchMode(bool batchMode)
{
    method->append_text (M("GENERAL_UNCHANGED"));
    method->set_active(procparams::RAWParams::BayerSensor::numMethods); // No name
    pixelShiftMotionCorrection->append_text (M("GENERAL_UNCHANGED"));
    pixelShiftMotionCorrection->set_active (4);
    pixelShiftMotionMethod->append_text (M("GENERAL_UNCHANGED"));
    pixelShiftMotionMethod->set_active (4);
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
    pixelShiftStddevFactorGreen->showEditedCB ();
    pixelShiftStddevFactorRed->showEditedCB ();
    pixelShiftStddevFactorBlue->showEditedCB ();
    pixelShiftEperIso->showEditedCB ();
    pixelShiftNreadIso->showEditedCB ();
    pixelShiftPrnu->showEditedCB ();
    pixelShiftSigma->showEditedCB ();
    pixelShiftSum->showEditedCB ();
    pixelShiftRedBlueWeight->showEditedCB ();
}

void BayerProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    dcbIterations->setDefault( defParams->raw.bayersensor.dcb_iterations);
    lmmseIterations->setDefault( defParams->raw.bayersensor.lmmse_iterations);
    pixelShiftMotion->setDefault( defParams->raw.bayersensor.pixelShiftMotion);
    pixelShiftStddevFactorGreen->setDefault( defParams->raw.bayersensor.pixelShiftStddevFactorGreen);
    pixelShiftStddevFactorRed->setDefault( defParams->raw.bayersensor.pixelShiftStddevFactorRed);
    pixelShiftStddevFactorBlue->setDefault( defParams->raw.bayersensor.pixelShiftStddevFactorBlue);
    pixelShiftEperIso->setDefault( defParams->raw.bayersensor.pixelShiftEperIso);
    pixelShiftNreadIso->setDefault( defParams->raw.bayersensor.pixelShiftNreadIso);
    pixelShiftPrnu->setDefault( defParams->raw.bayersensor.pixelShiftPrnu);
    pixelShiftSigma->setDefault( defParams->raw.bayersensor.pixelShiftSigma);
    pixelShiftSum->setDefault( defParams->raw.bayersensor.pixelShiftSum);
    pixelShiftRedBlueWeight->setDefault( defParams->raw.bayersensor.pixelShiftRedBlueWeight);
    ccSteps->setDefault (defParams->raw.bayersensor.ccSteps);

    if (pedited) {
        dcbIterations->setDefaultEditedState( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        lmmseIterations->setDefaultEditedState( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        pixelShiftMotion->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftMotion ? Edited : UnEdited);
        pixelShiftStddevFactorGreen->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftStddevFactorGreen ? Edited : UnEdited);
        pixelShiftStddevFactorRed->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftStddevFactorRed ? Edited : UnEdited);
        pixelShiftStddevFactorBlue->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftStddevFactorBlue ? Edited : UnEdited);
        pixelShiftEperIso->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftNreadIso->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftNreadIso ? Edited : UnEdited);
        pixelShiftPrnu->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftPrnu ? Edited : UnEdited);
        pixelShiftSigma->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftSigma ? Edited : UnEdited);
        pixelShiftSum->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftSum ? Edited : UnEdited);
        pixelShiftRedBlueWeight->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftRedBlueWeight ? Edited : UnEdited);
        ccSteps->setDefaultEditedState(pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
    } else {
        dcbIterations->setDefaultEditedState( Irrelevant );
        lmmseIterations->setDefaultEditedState( Irrelevant );
        pixelShiftMotion->setDefaultEditedState( Irrelevant );
        pixelShiftStddevFactorGreen->setDefaultEditedState( Irrelevant );
        pixelShiftStddevFactorRed->setDefaultEditedState( Irrelevant );
        pixelShiftStddevFactorBlue->setDefaultEditedState( Irrelevant );
        pixelShiftEperIso->setDefaultEditedState( Irrelevant );
        pixelShiftNreadIso->setDefaultEditedState( Irrelevant );
        pixelShiftPrnu->setDefaultEditedState( Irrelevant );
        pixelShiftSigma->setDefaultEditedState( Irrelevant );
        pixelShiftSum->setDefaultEditedState( Irrelevant );
        pixelShiftRedBlueWeight->setDefaultEditedState( Irrelevant );
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
            listener->panelChanged (EvPixelShiftMotion, a->getTextValue() );
        } else if (a == pixelShiftStddevFactorGreen) {
            listener->panelChanged (EvPixelShiftStddevFactorGreen, a->getTextValue() );
        } else if (a == pixelShiftStddevFactorRed) {
            listener->panelChanged (EvPixelShiftStddevFactorRed, a->getTextValue() );
        } else if (a == pixelShiftStddevFactorBlue) {
            listener->panelChanged (EvPixelShiftStddevFactorBlue, a->getTextValue() );
        } else if (a == pixelShiftEperIso) {
            listener->panelChanged (EvPixelShiftEperIso, a->getTextValue() );
        } else if (a == pixelShiftNreadIso) {
            listener->panelChanged (EvPixelShiftNreadIso, a->getTextValue() );
        } else if (a == pixelShiftPrnu) {
            listener->panelChanged (EvPixelShiftPrnu, a->getTextValue() );
        } else if (a == pixelShiftSigma) {
            listener->panelChanged (EvPixelShiftSigma, a->getTextValue() );
        } else if (a == pixelShiftSum) {
            listener->panelChanged (EvPixelShiftSum, a->getTextValue() );
        } else if (a == pixelShiftRedBlueWeight) {
            listener->panelChanged (EvPixelShiftRedBlueWeight, a->getTextValue() );
        }
    }
}

void BayerProcess::psMotionCorrectionChanged ()
{
    if(pixelShiftMotionCorrection->get_active_row_number() == 5) {
        pixelShiftBlur->set_sensitive(true);
        pixelShiftHoleFill->set_sensitive(true);
        pixelShiftSmooth->set_sensitive(pixelShiftBlur->get_active());
    } else {
        pixelShiftBlur->set_sensitive(false);
        pixelShiftHoleFill->set_sensitive(false);
        pixelShiftSmooth->set_sensitive(false);
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftMotionCorrection, pixelShiftMotionCorrection->get_active_text());
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

    if ( curSelection == procparams::RAWParams::BayerSensor::pixelshift) {
        if(pixelShiftMotionMethod->get_active_row_number() == 2) {
            pixelShiftOptions->show();
        } else {
            pixelShiftOptions->hide();
        }
        pixelShiftFrame->show();
    } else {
        pixelShiftFrame->hide();
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

void BayerProcess::pixelShiftMotionMethodChanged ()
{
    if(pixelShiftMotionMethod->get_active_row_number() == 0) {
        pixelShiftOptions->hide();
        pixelShiftShowMotion->hide();
        pixelShiftShowMotionMaskOnly->hide();
    } else if(pixelShiftMotionMethod->get_active_row_number() == 2) {
        pixelShiftOptions->show();
        pixelShiftShowMotion->show();
        pixelShiftShowMotionMaskOnly->show();
    } else {
        pixelShiftOptions->hide();
        pixelShiftShowMotion->show();
        pixelShiftShowMotionMaskOnly->show();
    }
    if (listener) {
        listener->panelChanged (EvPixelShiftMotionMethod, pixelShiftMotionMethod->get_active_text());
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
    pixelShiftShowMotionMaskOnly->set_sensitive(pixelShiftShowMotion->get_active ());
    if (listener) {
        listener->panelChanged (EvPixelshiftShowMotion, pixelShiftShowMotion->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftShowMotionMaskOnlyChanged ()
{
    if (batchMode) {
        if (pixelShiftShowMotionMaskOnly->get_inconsistent()) {
            pixelShiftShowMotionMaskOnly->set_inconsistent (false);
            pixelShiftShowMotionMaskOnlyconn.block (true);
            pixelShiftShowMotionMaskOnly->set_active (false);
            pixelShiftShowMotionMaskOnlyconn.block (false);
        } else if (lastDCBen) {
            pixelShiftShowMotionMaskOnly->set_inconsistent (true);
        }

        lastDCBen = pixelShiftShowMotionMaskOnly->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelshiftShowMotionMaskOnly, pixelShiftShowMotionMaskOnly->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
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
    pixelShiftSigma->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftSum->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftRedBlueWeight->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftStddevFactorGreen->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftStddevFactorRed->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftStddevFactorBlue->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftNonGreenHorizontal->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftNonGreenVertical->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftHoleFill->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5);
    pixelShiftMedian->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftMedian3->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMedian->get_active());
    pixelShiftGreen->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftBlur->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5);
    pixelShiftSmooth->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5 && pixelShiftBlur->get_active());
    pixelShiftExp0->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftNonGreenCross->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftNonGreenCross2->set_sensitive(pixelShiftAutomatic->get_active ());
    pixelShiftNonGreenAmaze->set_sensitive(pixelShiftAutomatic->get_active ());

    if (listener) {
        listener->panelChanged (EvPixelShiftAutomatic, pixelShiftAutomatic->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void BayerProcess::pixelShiftNonGreenHorizontalChanged ()
{
    if (batchMode) {
        if (pixelShiftNonGreenHorizontal->get_inconsistent()) {
            pixelShiftNonGreenHorizontal->set_inconsistent (false);
            pixelShiftNonGreenHorizontalconn.block (true);
            pixelShiftNonGreenHorizontal->set_active (false);
            pixelShiftNonGreenHorizontalconn.block (false);
        } else if (lastDCBen) {
            pixelShiftNonGreenHorizontal->set_inconsistent (true);
        }

        lastDCBen = pixelShiftNonGreenHorizontal->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftNonGreenHorizontal, pixelShiftNonGreenHorizontal->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftNonGreenVerticalChanged ()
{
    if (batchMode) {
        if (pixelShiftNonGreenVertical->get_inconsistent()) {
            pixelShiftNonGreenVertical->set_inconsistent (false);
            pixelShiftNonGreenVerticalconn.block (true);
            pixelShiftNonGreenVertical->set_active (false);
            pixelShiftNonGreenVerticalconn.block (false);
        } else if (lastDCBen) {
            pixelShiftNonGreenVertical->set_inconsistent (true);
        }

        lastDCBen = pixelShiftNonGreenVertical->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftNonGreenVertical, pixelShiftNonGreenVertical->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftHoleFillChanged ()
{
    if (batchMode) {
        if (pixelShiftHoleFill->get_inconsistent()) {
            pixelShiftHoleFill->set_inconsistent (false);
            pixelShiftHoleFillconn.block (true);
            pixelShiftHoleFill->set_active (false);
            pixelShiftHoleFillconn.block (false);
        } else if (lastDCBen) {
            pixelShiftHoleFill->set_inconsistent (true);
        }

        lastDCBen = pixelShiftHoleFill->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftHoleFill, pixelShiftHoleFill->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftMedianChanged ()
{
    if (batchMode) {
        if (pixelShiftMedian->get_inconsistent()) {
            pixelShiftMedian->set_inconsistent (false);
            pixelShiftMedianconn.block (true);
            pixelShiftMedian->set_active (false);
            pixelShiftMedianconn.block (false);
        } else if (lastDCBen) {
            pixelShiftMedian->set_inconsistent (true);
        }

        lastDCBen = pixelShiftMedian->get_active ();
    }

    pixelShiftMedian3->set_sensitive(pixelShiftMedian->get_active ());

    if (listener) {
        listener->panelChanged (EvPixelShiftMedian, pixelShiftMedian->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftMedian3Changed ()
{
    if (batchMode) {
        if (pixelShiftMedian3->get_inconsistent()) {
            pixelShiftMedian3->set_inconsistent (false);
            pixelShiftMedian3conn.block (true);
            pixelShiftMedian3->set_active (false);
            pixelShiftMedian3conn.block (false);
        } else if (lastDCBen) {
            pixelShiftMedian3->set_inconsistent (true);
        }

        lastDCBen = pixelShiftMedian3->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftMedian3, pixelShiftMedian3->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftGreenChanged ()
{
    if (batchMode) {
        if (pixelShiftGreen->get_inconsistent()) {
            pixelShiftGreen->set_inconsistent (false);
            pixelShiftGreenconn.block (true);
            pixelShiftGreen->set_active (false);
            pixelShiftGreenconn.block (false);
        } else if (lastDCBen) {
            pixelShiftGreen->set_inconsistent (true);
        }

        lastDCBen = pixelShiftGreen->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftGreen, pixelShiftGreen->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftBlurChanged ()
{
    if (batchMode) {
        if (pixelShiftBlur->get_inconsistent()) {
            pixelShiftBlur->set_inconsistent (false);
            pixelShiftBlurconn.block (true);
            pixelShiftBlur->set_active (false);
            pixelShiftBlurconn.block (false);
        } else if (lastDCBen) {
            pixelShiftBlur->set_inconsistent (true);
        }

        lastDCBen = pixelShiftBlur->get_active ();
    }

    pixelShiftSmooth->set_sensitive(pixelShiftBlur->get_active ());

    if (listener) {
        listener->panelChanged (EvPixelShiftBlur, pixelShiftBlur->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftSmoothChanged ()
{
    if (batchMode) {
        if (pixelShiftSmooth->get_inconsistent()) {
            pixelShiftSmooth->set_inconsistent (false);
            pixelShiftSmoothconn.block (true);
            pixelShiftSmooth->set_active (false);
            pixelShiftSmoothconn.block (false);
        } else if (lastDCBen) {
            pixelShiftSmooth->set_inconsistent (true);
        }

        lastDCBen = pixelShiftSmooth->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftSmooth, pixelShiftSmooth->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftExp0Changed ()
{
    if (batchMode) {
        if (pixelShiftExp0->get_inconsistent()) {
            pixelShiftExp0->set_inconsistent (false);
            pixelShiftExp0conn.block (true);
            pixelShiftExp0->set_active (false);
            pixelShiftExp0conn.block (false);
        } else if (lastDCBen) {
            pixelShiftExp0->set_inconsistent (true);
        }

        lastDCBen = pixelShiftExp0->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftExp0, pixelShiftExp0->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftNonGreenCrossChanged ()
{
    if (batchMode) {
        if (pixelShiftNonGreenCross->get_inconsistent()) {
            pixelShiftNonGreenCross->set_inconsistent (false);
            pixelShiftNonGreenCrossconn.block (true);
            pixelShiftNonGreenCross->set_active (false);
            pixelShiftNonGreenCrossconn.block (false);
        } else if (lastDCBen) {
            pixelShiftNonGreenCross->set_inconsistent (true);
        }

        lastDCBen = pixelShiftNonGreenCross->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftNonGreenCross, pixelShiftNonGreenCross->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftNonGreenCross2Changed ()
{
    if (batchMode) {
        if (pixelShiftNonGreenCross2->get_inconsistent()) {
            pixelShiftNonGreenCross2->set_inconsistent (false);
            pixelShiftNonGreenCross2conn.block (true);
            pixelShiftNonGreenCross2->set_active (false);
            pixelShiftNonGreenCross2conn.block (false);
        } else if (lastDCBen) {
            pixelShiftNonGreenCross2->set_inconsistent (true);
        }

        lastDCBen = pixelShiftNonGreenCross2->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftGreenAmaze, pixelShiftNonGreenCross2->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftNonGreenAmazeChanged ()
{
    if (batchMode) {
        if (pixelShiftNonGreenAmaze->get_inconsistent()) {
            pixelShiftNonGreenAmaze->set_inconsistent (false);
            pixelShiftNonGreenAmazeconn.block (true);
            pixelShiftNonGreenAmaze->set_active (false);
            pixelShiftNonGreenAmazeconn.block (false);
        } else if (lastDCBen) {
            pixelShiftNonGreenAmaze->set_inconsistent (true);
        }

        lastDCBen = pixelShiftNonGreenAmaze->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftNonGreenAmaze, pixelShiftNonGreenAmaze->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
