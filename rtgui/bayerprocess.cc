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

    imageNumberBox = Gtk::manage (new Gtk::HBox ());
    imageNumberBox->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_IMAGENUM") + ": ")), Gtk::PACK_SHRINK, 4);
    imageNumber = Gtk::manage (new MyComboBoxText ());
    imageNumber->append("1");
    imageNumber->append("2");
    imageNumber->append("3");
    imageNumber->append("4");
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

    dcbIterations = Gtk::manage (new Adjuster (M("TP_RAW_DCBITERATIONS"), 0, 5, 1, 2));
    dcbIterations->setAdjusterListener (this);

    if (dcbIterations->delay < options.adjusterMaxDelay) {
        dcbIterations->delay = options.adjusterMaxDelay;
    }

    dcbIterations->show();
    dcbEnhance = Gtk::manage (new MyCheckButton(M("TP_RAW_DCBENHANCE")));
    dcbOptions->pack_start(*dcbIterations);
    dcbOptions->pack_start(*dcbEnhance);
    pack_start( *dcbOptions, Gtk::PACK_SHRINK, 4);

    lmmseOptions = Gtk::manage (new Gtk::VBox ());

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

    pixelShiftEqualBright = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTEQUALBRIGHT")));
    pixelShiftEqualBright->set_tooltip_text (M("TP_RAW_PIXELSHIFTEQUALBRIGHT_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftEqualBright);

    Gtk::HBox* hb3 = Gtk::manage (new Gtk::HBox ());
    hb3->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_PIXELSHIFTMOTIONMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    pixelShiftMotionMethod = Gtk::manage (new MyComboBoxText ());
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_OFF"));
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_AUTO"));
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_CUSTOM"));
    pixelShiftMotionMethod->set_active(RAWParams::BayerSensor::ePSMotionCorrectionMethod::Automatic);
    pixelShiftMotionMethod->show();
    hb3->pack_start(*pixelShiftMotionMethod);
    pixelShiftFrame->pack_start(*hb3);

    pixelShiftOptions = Gtk::manage (new Gtk::VBox ());
    pixelShiftOptions->set_border_width(0);

    pixelShiftShowMotion = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTSHOWMOTION")));
    pixelShiftShowMotion->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTION_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftShowMotion);

    pixelShiftShowMotionMaskOnly = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY")));
    pixelShiftShowMotionMaskOnly->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftShowMotionMaskOnly);

#ifdef PIXELSHIFTDEV
    pixelShiftAutomatic = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTADAPTIVE")));
    pixelShiftOptions->pack_start(*pixelShiftAutomatic);
#endif
    pixelShiftGreen = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTGREEN")));
    pixelShiftOptions->pack_start(*pixelShiftGreen);

    pixelShiftNonGreenCross = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTNONGREENCROSS")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenCross);

    pixelShiftHoleFill = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTHOLEFILL")));
    pixelShiftHoleFill->set_tooltip_text (M("TP_RAW_PIXELSHIFTHOLEFILL_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftHoleFill);

    pixelShiftBlur = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTBLUR")));
    pixelShiftBlur->set_tooltip_text (M("TP_RAW_PIXELSHIFTSIGMA_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftBlur);

    pixelShiftSigma = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSIGMA"), 0.5, 25, 0.1, 1.0));
    pixelShiftSigma->set_tooltip_text (M("TP_RAW_PIXELSHIFTSIGMA_TOOLTIP"));
    pixelShiftSigma->setAdjusterListener (this);

    if (pixelShiftSigma->delay < options.adjusterMaxDelay) {
        pixelShiftSigma->delay = options.adjusterMaxDelay;
    }

    pixelShiftSigma->show();
    pixelShiftOptions->pack_start(*pixelShiftSigma);


    pixelShiftSmooth = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSMOOTH"), 0, 1, 0.05, 0.7));
    pixelShiftSmooth->set_tooltip_text (M("TP_RAW_PIXELSHIFTSMOOTH_TOOLTIP"));
    pixelShiftSmooth->setAdjusterListener (this);

    if (pixelShiftSmooth->delay < options.adjusterMaxDelay) {
        pixelShiftSmooth->delay = options.adjusterMaxDelay;
    }

    pixelShiftSmooth->show();
    pixelShiftOptions->pack_start(*pixelShiftSmooth);

    pixelShiftEperIso = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTEPERISO"), -5.0, 5.0, 0.05, 0.0));
    pixelShiftEperIso->set_tooltip_text(M("TP_RAW_PIXELSHIFTEPERISO_TOOLTIP"));
    pixelShiftEperIso->setAdjusterListener (this);

    if (pixelShiftEperIso->delay < options.adjusterMaxDelay) {
        pixelShiftEperIso->delay = options.adjusterMaxDelay;
    }

    pixelShiftEperIso->show();
    pixelShiftOptions->pack_start(*pixelShiftEperIso);


    pixelShiftMedian = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTMEDIAN")));
    pixelShiftMedian->set_tooltip_text (M("TP_RAW_PIXELSHIFTMEDIAN_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftMedian);


#ifdef PIXELSHIFTDEV
    pixelShiftMedian3 = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTMEDIAN3")));
    pixelShiftMedian3->set_tooltip_text (M("TP_RAW_PIXELSHIFTMEDIAN3_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftMedian3);

    pixelShiftNonGreenCross2 = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTNONGREENCROSS2")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenCross2);

    pixelShiftNonGreenAmaze = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTNONGREENAMAZE")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenAmaze);

    pixelShiftNonGreenHorizontal = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTNONGREENHORIZONTAL")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenHorizontal);

    pixelShiftNonGreenVertical = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTNONGREENVERTICAL")));
    pixelShiftOptions->pack_start(*pixelShiftNonGreenVertical);

    pixelShiftExp0 = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTEXP0")));
    pixelShiftOptions->pack_start(*pixelShiftExp0);
#endif
    pixelShiftLmmse = Gtk::manage (new MyCheckButton(M("TP_RAW_PIXELSHIFTLMMSE")));
    pixelShiftLmmse->set_tooltip_text (M("TP_RAW_PIXELSHIFTLMMSE_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftLmmse);

#ifdef PIXELSHIFTDEV
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
    pixelShiftMotionCorrection->append("1x1");
    pixelShiftMotionCorrection->append("1x2");
    pixelShiftMotionCorrection->append("3x3");
    pixelShiftMotionCorrection->append("5x5");
    pixelShiftMotionCorrection->append("7x7");
    pixelShiftMotionCorrection->append("3x3 new");
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
#endif

#ifdef PIXELSHIFTDEV
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
#endif

    pixelShiftFrame->pack_start(*pixelShiftOptions);
    pixelShiftOptions->hide();

    pack_start( *pixelShiftFrame, Gtk::PACK_SHRINK, 4);

    method->connect(method->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::methodChanged) ));
    imageNumber->connect(imageNumber->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::imageNumberChanged) ));
    dcbEnhance->connect ( dcbEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::dcbEnhanceChanged), true));
    pixelShiftMotionMethod->connect(pixelShiftMotionMethod->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::pixelShiftMotionMethodChanged) ));
    pixelShiftShowMotion->connect(pixelShiftShowMotion->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftShowMotionChanged), true));
    pixelShiftShowMotionMaskOnly->connect(pixelShiftShowMotionMaskOnly->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftShowMotionMaskOnlyChanged), true));
    pixelShiftHoleFill->connect(pixelShiftHoleFill->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftHoleFillChanged), true));
    pixelShiftMedian->connect(pixelShiftMedian->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftMedianChanged), true));
    pixelShiftGreen->connect(pixelShiftGreen->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftGreenChanged), true));
    pixelShiftBlur->connect(pixelShiftBlur->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftBlurChanged), true));
    pixelShiftLmmse->connect(pixelShiftLmmse->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftLmmseChanged), true));
    pixelShiftEqualBright->connect(pixelShiftEqualBright->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftEqualBrightChanged), true));
    pixelShiftNonGreenCross->connect(pixelShiftNonGreenCross->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenCrossChanged), true));
#ifdef PIXELSHIFTDEV
    pixelShiftMotionCorrection->connect(pixelShiftMotionCorrection->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::psMotionCorrectionChanged) ));
    pixelShiftAutomatic->connect(pixelShiftAutomatic->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftAutomaticChanged), true));
    pixelShiftNonGreenHorizontal->connect(pixelShiftNonGreenHorizontal->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenHorizontalChanged), true));
    pixelShiftNonGreenVertical->connect(pixelShiftNonGreenVertical->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenVerticalChanged), true));
    pixelShiftMedian3->connect(pixelShiftMedian3->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftMedian3Changed), true));
    pixelShiftExp0->connect(pixelShiftExp0->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftExp0Changed), true));
    pixelShiftNonGreenCross2->connect(pixelShiftNonGreenCross2->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenCross2Changed), true));
    pixelShiftNonGreenAmaze->connect(pixelShiftNonGreenAmaze->signal_toggled().connect ( sigc::mem_fun(*this, &BayerProcess::pixelShiftNonGreenAmazeChanged), true));
#endif
}


void BayerProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    method->block (true);
    dcbEnhance->block (true);
    imageNumber->block (true);
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
        pixelShiftHoleFill->set_inconsistent(!pedited->raw.bayersensor.pixelShiftHoleFill);
        pixelShiftMedian->set_inconsistent(!pedited->raw.bayersensor.pixelShiftMedian);
        pixelShiftGreen->set_inconsistent(!pedited->raw.bayersensor.pixelShiftGreen);
        pixelShiftBlur->set_inconsistent(!pedited->raw.bayersensor.pixelShiftBlur);
        pixelShiftSmooth->setEditedState ( pedited->raw.bayersensor.pixelShiftSmooth ? Edited : UnEdited);
        pixelShiftLmmse->set_inconsistent(!pedited->raw.bayersensor.pixelShiftLmmse);
        pixelShiftEqualBright->set_inconsistent(!pedited->raw.bayersensor.pixelShiftEqualBright);
        pixelShiftNonGreenCross->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenCross);
        lmmseIterations->setEditedState ( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        pixelShiftEperIso->setEditedState ( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftSigma->setEditedState ( pedited->raw.bayersensor.pixelShiftSigma ? Edited : UnEdited);
#ifdef PIXELSHIFTDEV
        pixelShiftNreadIso->setEditedState ( pedited->raw.bayersensor.pixelShiftNreadIso ? Edited : UnEdited);
        pixelShiftPrnu->setEditedState ( pedited->raw.bayersensor.pixelShiftPrnu ? Edited : UnEdited);
        pixelShiftStddevFactorGreen->setEditedState ( pedited->raw.bayersensor.pixelShiftStddevFactorGreen ? Edited : UnEdited);
        pixelShiftStddevFactorRed->setEditedState ( pedited->raw.bayersensor.pixelShiftStddevFactorRed ? Edited : UnEdited);
        pixelShiftStddevFactorBlue->setEditedState ( pedited->raw.bayersensor.pixelShiftStddevFactorBlue ? Edited : UnEdited);
        pixelShiftSum->setEditedState ( pedited->raw.bayersensor.pixelShiftSum ? Edited : UnEdited);
        pixelShiftAutomatic->set_inconsistent(!pedited->raw.bayersensor.pixelShiftAutomatic);
        pixelShiftNonGreenHorizontal->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenHorizontal);
        pixelShiftNonGreenVertical->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenVertical);
        pixelShiftMedian3->set_inconsistent(!pedited->raw.bayersensor.pixelShiftMedian3);
        pixelShiftExp0->set_inconsistent(!pedited->raw.bayersensor.pixelShiftExp0);
        pixelShiftNonGreenCross2->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenCross2);
        pixelShiftNonGreenAmaze->set_inconsistent(!pedited->raw.bayersensor.pixelShiftNonGreenAmaze);
        pixelShiftMotion->setEditedState ( pedited->raw.bayersensor.pixelShiftMotion ? Edited : UnEdited);
        pixelShiftRedBlueWeight->setEditedState ( pedited->raw.bayersensor.pixelShiftRedBlueWeight ? Edited : UnEdited);
#endif

        if(!pedited->raw.bayersensor.method) {
            method->set_active(procparams::RAWParams::BayerSensor::numMethods);    // No name
        }
        if(!pedited->raw.bayersensor.imageNum) {
            imageNumber->set_active_text(M("GENERAL_UNCHANGED"));
        }
#ifdef PIXELSHIFTDEV
        if(!pedited->raw.bayersensor.pixelShiftMotionCorrection) {
            pixelShiftMotionCorrection->set_active_text(M("GENERAL_UNCHANGED"));
        }
#endif
        if(!pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod) {
            pixelShiftMotionMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }
    }

    //allEnhance->set_active(pp->raw.bayersensor.all_enhance);

    dcbIterations->setValue (pp->raw.bayersensor.dcb_iterations);
    dcbEnhance->set_active(pp->raw.bayersensor.dcb_enhance);
    pixelShiftShowMotion->set_active(pp->raw.bayersensor.pixelshiftShowMotion);
    pixelShiftShowMotionMaskOnly->set_sensitive(pp->raw.bayersensor.pixelshiftShowMotion);
    pixelShiftShowMotionMaskOnly->set_active(pp->raw.bayersensor.pixelshiftShowMotionMaskOnly);
    pixelShiftHoleFill->set_active(pp->raw.bayersensor.pixelShiftHoleFill);
    pixelShiftMedian->set_active(pp->raw.bayersensor.pixelShiftMedian);
    pixelShiftGreen->set_active(pp->raw.bayersensor.pixelShiftGreen);
    pixelShiftBlur->set_active(pp->raw.bayersensor.pixelShiftBlur);
    pixelShiftSmooth->set_sensitive (pp->raw.bayersensor.pixelShiftBlur);
    pixelShiftSmooth->setValue (pp->raw.bayersensor.pixelShiftSmoothFactor);
    pixelShiftLmmse->set_active(pp->raw.bayersensor.pixelShiftLmmse);
    pixelShiftEqualBright->set_active(pp->raw.bayersensor.pixelShiftEqualBright);
    pixelShiftNonGreenCross->set_active(pp->raw.bayersensor.pixelShiftNonGreenCross);
    ccSteps->setValue (pp->raw.bayersensor.ccSteps);
    lmmseIterations->setValue (pp->raw.bayersensor.lmmse_iterations);
    pixelShiftMotionMethod->set_active ((int)pp->raw.bayersensor.pixelShiftMotionCorrectionMethod);
    pixelShiftEperIso->setValue (pp->raw.bayersensor.pixelShiftEperIso);
    pixelShiftSigma->setValue (pp->raw.bayersensor.pixelShiftSigma);
    pixelShiftSigma->set_sensitive (pp->raw.bayersensor.pixelShiftBlur);
#ifdef PIXELSHIFTDEV
    pixelShiftStddevFactorGreen->setValue (pp->raw.bayersensor.pixelShiftStddevFactorGreen);
    pixelShiftStddevFactorRed->setValue (pp->raw.bayersensor.pixelShiftStddevFactorRed);
    pixelShiftStddevFactorBlue->setValue (pp->raw.bayersensor.pixelShiftStddevFactorBlue);
    pixelShiftSum->setValue (pp->raw.bayersensor.pixelShiftSum);
    pixelShiftMedian3->set_active(pp->raw.bayersensor.pixelShiftMedian3);
    pixelShiftMedian3->set_sensitive(pixelShiftMedian->get_active());
    pixelShiftAutomatic->set_active(pp->raw.bayersensor.pixelShiftAutomatic);
    pixelShiftNonGreenHorizontal->set_active(pp->raw.bayersensor.pixelShiftNonGreenHorizontal);
    pixelShiftNonGreenVertical->set_active(pp->raw.bayersensor.pixelShiftNonGreenVertical);
    pixelShiftExp0->set_active(pp->raw.bayersensor.pixelShiftExp0);
    pixelShiftNonGreenCross2->set_active(pp->raw.bayersensor.pixelShiftNonGreenCross2);
    pixelShiftNonGreenAmaze->set_active(pp->raw.bayersensor.pixelShiftNonGreenAmaze);
    pixelShiftMotion->setValue (pp->raw.bayersensor.pixelShiftMotion);
    pixelShiftMotionCorrection->set_active ((int)pp->raw.bayersensor.pixelShiftMotionCorrection);
    pixelShiftHoleFill->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5);
    pixelShiftBlur->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5);
    pixelShiftSmooth->set_sensitive(pixelShiftAutomatic->get_active () && pixelShiftMotionCorrection->get_active_row_number() == 5 && pixelShiftBlur->get_active());
    pixelShiftNreadIso->setValue (pp->raw.bayersensor.pixelShiftNreadIso);
    pixelShiftPrnu->setValue (pp->raw.bayersensor.pixelShiftPrnu);
    pixelShiftRedBlueWeight->setValue (pp->raw.bayersensor.pixelShiftRedBlueWeight);
#endif
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

    dcbEnhance->setLastActive(pp->raw.bayersensor.dcb_enhance);
    pixelShiftShowMotion->setLastActive(pp->raw.bayersensor.pixelshiftShowMotion);
    pixelShiftShowMotionMaskOnly->setLastActive(pp->raw.bayersensor.pixelshiftShowMotionMaskOnly);
    pixelShiftNonGreenCross->setLastActive(pp->raw.bayersensor.pixelShiftNonGreenCross);
    pixelShiftGreen->setLastActive(pp->raw.bayersensor.pixelShiftGreen);
    pixelShiftBlur->setLastActive(pp->raw.bayersensor.pixelShiftBlur);
    pixelShiftHoleFill->setLastActive(pp->raw.bayersensor.pixelShiftHoleFill);
    pixelShiftMedian->setLastActive(pp->raw.bayersensor.pixelShiftMedian);
    pixelShiftLmmse->setLastActive(pp->raw.bayersensor.pixelShiftLmmse);
    pixelShiftEqualBright->setLastActive(pp->raw.bayersensor.pixelShiftEqualBright);
#ifdef PIXELSHIFTDEV
    pixelShiftMedian3->setLastActive(pp->raw.bayersensor.pixelShiftMedian3);
    pixelShiftAutomatic->setLastActive(pp->raw.bayersensor.pixelShiftAutomatic);
    pixelShiftNonGreenHorizontal->setLastActive(pp->raw.bayersensor.pixelShiftNonGreenHorizontal);
    pixelShiftNonGreenVertical->setLastActive(pp->raw.bayersensor.pixelShiftNonGreenVertical);
    pixelShiftNonGreenCross2->setLastActive(pp->raw.bayersensor.pixelShiftNonGreenCross2);
    pixelShiftNonGreenAmaze->setLastActive(pp->raw.bayersensor.pixelShiftNonGreenAmaze);
    pixelShiftExp0->setLastActive(pp->raw.bayersensor.pixelShiftExp0);
#endif
    
    //lastALLen = pp->raw.bayersensor.all_enhance;

    method->block (false);
#ifdef PIXELSHIFTDEV
    pixelShiftMotionCorrection->block (false);
#endif
    imageNumber->block (false);
    dcbEnhance->block (false);
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
    pp->raw.bayersensor.pixelShiftMotionCorrectionMethod = (RAWParams::BayerSensor::ePSMotionCorrectionMethod)pixelShiftMotionMethod->get_active_row_number();
    pp->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getValue();
    pp->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getValue();
    pp->raw.bayersensor.pixelshiftShowMotion = pixelShiftShowMotion->get_active();
    pp->raw.bayersensor.pixelshiftShowMotionMaskOnly = pixelShiftShowMotionMaskOnly->get_active();
    pp->raw.bayersensor.pixelShiftHoleFill = pixelShiftHoleFill->get_active();
    pp->raw.bayersensor.pixelShiftMedian = pixelShiftMedian->get_active();
    pp->raw.bayersensor.pixelShiftGreen = pixelShiftGreen->get_active();
    pp->raw.bayersensor.pixelShiftBlur = pixelShiftBlur->get_active();
    pp->raw.bayersensor.pixelShiftSmoothFactor = pixelShiftSmooth->getValue();
    pp->raw.bayersensor.pixelShiftLmmse = pixelShiftLmmse->get_active();
    pp->raw.bayersensor.pixelShiftEqualBright = pixelShiftEqualBright->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenCross = pixelShiftNonGreenCross->get_active();
#ifdef PIXELSHIFTDEV
    pp->raw.bayersensor.pixelShiftStddevFactorGreen = pixelShiftStddevFactorGreen->getValue();
    pp->raw.bayersensor.pixelShiftStddevFactorRed = pixelShiftStddevFactorRed->getValue();
    pp->raw.bayersensor.pixelShiftStddevFactorBlue = pixelShiftStddevFactorBlue->getValue();
    pp->raw.bayersensor.pixelShiftSum = pixelShiftSum->getValue();
    pp->raw.bayersensor.pixelShiftMedian3 = pixelShiftMedian3->get_active();
    pp->raw.bayersensor.pixelShiftMotion = pixelShiftMotion->getIntValue();
    pp->raw.bayersensor.pixelShiftMotionCorrection = (RAWParams::BayerSensor::ePSMotionCorrection)pixelShiftMotionCorrection->get_active_row_number();
    pp->raw.bayersensor.pixelShiftAutomatic = pixelShiftAutomatic->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenHorizontal = pixelShiftNonGreenHorizontal->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenVertical = pixelShiftNonGreenVertical->get_active();
    pp->raw.bayersensor.pixelShiftExp0 = pixelShiftExp0->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenCross2 = pixelShiftNonGreenCross2->get_active();
    pp->raw.bayersensor.pixelShiftNonGreenAmaze = pixelShiftNonGreenAmaze->get_active();
    pp->raw.bayersensor.pixelShiftNreadIso = pixelShiftNreadIso->getValue();
    pp->raw.bayersensor.pixelShiftPrnu = pixelShiftPrnu->getValue();
    pp->raw.bayersensor.pixelShiftRedBlueWeight = pixelShiftRedBlueWeight->getValue();
#endif

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
        pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod = pixelShiftMotionMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getEditedState ();
        pedited->raw.bayersensor.pixelshiftShowMotion = !pixelShiftShowMotion->get_inconsistent();
        pedited->raw.bayersensor.pixelshiftShowMotionMaskOnly = !pixelShiftShowMotionMaskOnly->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftHoleFill = !pixelShiftHoleFill->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftMedian = !pixelShiftMedian->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftGreen = !pixelShiftGreen->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftBlur = !pixelShiftBlur->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftSmooth = pixelShiftSmooth->getEditedState();
        pedited->raw.bayersensor.pixelShiftLmmse = !pixelShiftLmmse->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftEqualBright = !pixelShiftEqualBright->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenCross = !pixelShiftNonGreenCross->get_inconsistent();
#ifdef PIXELSHIFTDEV
        pedited->raw.bayersensor.pixelShiftStddevFactorGreen = pixelShiftStddevFactorGreen->getEditedState ();
        pedited->raw.bayersensor.pixelShiftStddevFactorRed = pixelShiftStddevFactorRed->getEditedState ();
        pedited->raw.bayersensor.pixelShiftStddevFactorBlue = pixelShiftStddevFactorBlue->getEditedState ();
        pedited->raw.bayersensor.pixelShiftSum = pixelShiftSum->getEditedState ();
        pedited->raw.bayersensor.pixelShiftMedian3 = !pixelShiftMedian3->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftMotion = pixelShiftMotion->getEditedState ();
        pedited->raw.bayersensor.pixelShiftMotionCorrection = pixelShiftMotionCorrection->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.pixelShiftAutomatic = !pixelShiftAutomatic->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenHorizontal = !pixelShiftNonGreenHorizontal->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenVertical = !pixelShiftNonGreenVertical->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftExp0 = !pixelShiftExp0->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenCross2 = !pixelShiftNonGreenCross2->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenAmaze = !pixelShiftNonGreenAmaze->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNreadIso = pixelShiftNreadIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftPrnu = pixelShiftPrnu->getEditedState ();
        pedited->raw.bayersensor.pixelShiftRedBlueWeight = pixelShiftRedBlueWeight->getEditedState ();
#endif
    }
}

void BayerProcess::setBatchMode(bool batchMode)
{
    method->append (M("GENERAL_UNCHANGED"));
    method->set_active(procparams::RAWParams::BayerSensor::numMethods); // No name
#ifdef PIXELSHIFTDEV
    pixelShiftMotionCorrection->append (M("GENERAL_UNCHANGED"));
    pixelShiftMotionCorrection->set_active (4);
#endif
    pixelShiftMotionMethod->append (M("GENERAL_UNCHANGED"));
    pixelShiftMotionMethod->set_active (4);
    imageNumber->append (M("GENERAL_UNCHANGED"));
    imageNumber->set_active(4);
    dcbOptions->hide();
    lmmseOptions->hide();
    pixelShiftOptions->hide();
    ToolPanel::setBatchMode (batchMode);
    ccSteps->showEditedCB ();
    dcbIterations->showEditedCB ();
    lmmseIterations->showEditedCB ();
#ifdef PIXELSHIFTDEV
    pixelShiftMotion->showEditedCB ();
    pixelShiftSum->showEditedCB ();
    pixelShiftStddevFactorGreen->showEditedCB ();
    pixelShiftStddevFactorRed->showEditedCB ();
    pixelShiftStddevFactorBlue->showEditedCB ();
    pixelShiftNreadIso->showEditedCB ();
    pixelShiftPrnu->showEditedCB ();
    pixelShiftRedBlueWeight->showEditedCB ();
#endif
    pixelShiftEperIso->showEditedCB ();
    pixelShiftSigma->showEditedCB ();
}

void BayerProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    dcbIterations->setDefault( defParams->raw.bayersensor.dcb_iterations);
    lmmseIterations->setDefault( defParams->raw.bayersensor.lmmse_iterations);
#ifdef PIXELSHIFTDEV
    pixelShiftMotion->setDefault( defParams->raw.bayersensor.pixelShiftMotion);
    pixelShiftSum->setDefault( defParams->raw.bayersensor.pixelShiftSum);
    pixelShiftStddevFactorGreen->setDefault( defParams->raw.bayersensor.pixelShiftStddevFactorGreen);
    pixelShiftStddevFactorRed->setDefault( defParams->raw.bayersensor.pixelShiftStddevFactorRed);
    pixelShiftStddevFactorBlue->setDefault( defParams->raw.bayersensor.pixelShiftStddevFactorBlue);
    pixelShiftNreadIso->setDefault( defParams->raw.bayersensor.pixelShiftNreadIso);
    pixelShiftPrnu->setDefault( defParams->raw.bayersensor.pixelShiftPrnu);
    pixelShiftRedBlueWeight->setDefault( defParams->raw.bayersensor.pixelShiftRedBlueWeight);
#endif
    pixelShiftEperIso->setDefault( defParams->raw.bayersensor.pixelShiftEperIso);
    pixelShiftSigma->setDefault( defParams->raw.bayersensor.pixelShiftSigma);
    ccSteps->setDefault (defParams->raw.bayersensor.ccSteps);

    if (pedited) {
        dcbIterations->setDefaultEditedState( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        lmmseIterations->setDefaultEditedState( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
#ifdef PIXELSHIFTDEV
        pixelShiftMotion->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftMotion ? Edited : UnEdited);
        pixelShiftSum->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftSum ? Edited : UnEdited);
        pixelShiftStddevFactorGreen->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftStddevFactorGreen ? Edited : UnEdited);
        pixelShiftStddevFactorRed->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftStddevFactorRed ? Edited : UnEdited);
        pixelShiftStddevFactorBlue->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftStddevFactorBlue ? Edited : UnEdited);
        pixelShiftNreadIso->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftNreadIso ? Edited : UnEdited);
        pixelShiftPrnu->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftPrnu ? Edited : UnEdited);
        pixelShiftRedBlueWeight->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftRedBlueWeight ? Edited : UnEdited);
#endif
        pixelShiftEperIso->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftSigma->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftSigma ? Edited : UnEdited);
        ccSteps->setDefaultEditedState(pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
    } else {
        dcbIterations->setDefaultEditedState( Irrelevant );
        lmmseIterations->setDefaultEditedState( Irrelevant );
#ifdef PIXELSHIFTDEV
        pixelShiftMotion->setDefaultEditedState( Irrelevant );
        pixelShiftSum->setDefaultEditedState( Irrelevant );
        pixelShiftStddevFactorGreen->setDefaultEditedState( Irrelevant );
        pixelShiftStddevFactorRed->setDefaultEditedState( Irrelevant );
        pixelShiftStddevFactorBlue->setDefaultEditedState( Irrelevant );
        pixelShiftNreadIso->setDefaultEditedState( Irrelevant );
        pixelShiftPrnu->setDefaultEditedState( Irrelevant );
        pixelShiftRedBlueWeight->setDefaultEditedState( Irrelevant );
#endif
        pixelShiftEperIso->setDefaultEditedState( Irrelevant );
        pixelShiftSigma->setDefaultEditedState( Irrelevant );
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
#ifdef PIXELSHIFTDEV
        } else if (a == pixelShiftMotion) {
            listener->panelChanged (EvPixelShiftMotion, a->getTextValue() );
        } else if (a == pixelShiftSum) {
            listener->panelChanged (EvPixelShiftSum, a->getTextValue() );
        } else if (a == pixelShiftStddevFactorGreen) {
            listener->panelChanged (EvPixelShiftStddevFactorGreen, a->getTextValue() );
        } else if (a == pixelShiftStddevFactorRed) {
            listener->panelChanged (EvPixelShiftStddevFactorRed, a->getTextValue() );
        } else if (a == pixelShiftStddevFactorBlue) {
            listener->panelChanged (EvPixelShiftStddevFactorBlue, a->getTextValue() );
        } else if (a == pixelShiftNreadIso) {
            listener->panelChanged (EvPixelShiftNreadIso, a->getTextValue() );
        } else if (a == pixelShiftPrnu) {
            listener->panelChanged (EvPixelShiftPrnu, a->getTextValue() );
        } else if (a == pixelShiftRedBlueWeight) {
            listener->panelChanged (EvPixelShiftRedBlueWeight, a->getTextValue() );
#endif
        } else if (a == pixelShiftEperIso) {
            listener->panelChanged (EvPixelShiftEperIso, a->getTextValue() );
        } else if (a == pixelShiftSigma) {
            listener->panelChanged (EvPixelShiftSigma, a->getTextValue() );
        } else if (a == pixelShiftSmooth) {
            listener->panelChanged (EvPixelShiftSmooth, a->getTextValue() );
        }
    }
}

#ifdef PIXELSHIFTDEV
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
#endif
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
        if(pixelShiftMotionMethod->get_active_row_number() == RAWParams::BayerSensor::ePSMotionCorrectionMethod::Custom) {
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
            dcbEnhance->block (true);
            dcbEnhance->set_active (false);
            dcbEnhance->block (false);
        } else if (dcbEnhance->getLastActive()) {
            dcbEnhance->set_inconsistent (true);
        }

        dcbEnhance->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvDemosaicDCBEnhanced, dcbEnhance->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftMotionMethodChanged ()
{
    if(pixelShiftMotionMethod->get_active_row_number() == RAWParams::BayerSensor::ePSMotionCorrectionMethod::Off) {
        pixelShiftOptions->hide();
        pixelShiftShowMotion->hide();
        pixelShiftShowMotionMaskOnly->hide();
    } else if(pixelShiftMotionMethod->get_active_row_number() == RAWParams::BayerSensor::ePSMotionCorrectionMethod::Custom) {
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
            pixelShiftShowMotion->block (true);
            pixelShiftShowMotion->set_active (false);
            pixelShiftShowMotion->block (false);
        } else if (pixelShiftShowMotion->getLastActive()) {
            pixelShiftShowMotion->set_inconsistent (true);
        }

        pixelShiftShowMotion->setLastActive();
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
            pixelShiftShowMotionMaskOnly->block (true);
            pixelShiftShowMotionMaskOnly->set_active (false);
            pixelShiftShowMotionMaskOnly->block (false);
        } else if (pixelShiftShowMotionMaskOnly->getLastActive()) {
            pixelShiftShowMotionMaskOnly->set_inconsistent (true);
        }

        pixelShiftShowMotionMaskOnly->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelshiftShowMotionMaskOnly, pixelShiftShowMotionMaskOnly->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

#ifdef PIXELSHIFTDEV
void BayerProcess::pixelShiftAutomaticChanged ()
{
    if (batchMode) {
        if (pixelShiftAutomatic->get_inconsistent()) {
            pixelShiftAutomatic->set_inconsistent (false);
            pixelShiftAutomatic->block (true);
            pixelShiftAutomatic->set_active (false);
            pixelShiftAutomatic->block (false);
        } else if (pixelShiftAutomatic->getLastActive()) {
            pixelShiftAutomatic->set_inconsistent (true);
        }

        pixelShiftAutomatic->setLastActive();
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
            pixelShiftNonGreenHorizontal->block (true);
            pixelShiftNonGreenHorizontal->set_active (false);
            pixelShiftNonGreenHorizontal->block (false);
        } else if (pixelShiftNonGreenHorizontal->getLastActive()) {
            pixelShiftNonGreenHorizontal->set_inconsistent (true);
        }

        pixelShiftNonGreenHorizontal->setLastActive();
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
            pixelShiftNonGreenVertical->block (true);
            pixelShiftNonGreenVertical->set_active (false);
            pixelShiftNonGreenVertical->block (false);
        } else if (pixelShiftNonGreenVertical->getLastActive()) {
            pixelShiftNonGreenVertical->set_inconsistent (true);
        }

        pixelShiftNonGreenVertical->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftNonGreenVertical, pixelShiftNonGreenVertical->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
#endif

void BayerProcess::pixelShiftHoleFillChanged ()
{
    if (batchMode) {
        if (pixelShiftHoleFill->get_inconsistent()) {
            pixelShiftHoleFill->set_inconsistent (false);
            pixelShiftHoleFill->block (true);
            pixelShiftHoleFill->set_active (false);
            pixelShiftHoleFill->block (false);
        } else if (pixelShiftHoleFill->getLastActive()) {
            pixelShiftHoleFill->set_inconsistent (true);
        }

        pixelShiftHoleFill->setLastActive();
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
            pixelShiftMedian->block (true);
            pixelShiftMedian->set_active (false);
            pixelShiftMedian->block (false);
        } else if (pixelShiftMedian->getLastActive()) {
            pixelShiftMedian->set_inconsistent (true);
        }

        pixelShiftMedian->setLastActive();
    }
#ifdef PIXELSHIFTDEV
    pixelShiftMedian3->set_sensitive(pixelShiftMedian->get_active ());
#endif
    if (listener) {
        listener->panelChanged (EvPixelShiftMedian, pixelShiftMedian->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
#ifdef PIXELSHIFTDEV
void BayerProcess::pixelShiftMedian3Changed ()
{
    if (batchMode) {
        if (pixelShiftMedian3->get_inconsistent()) {
            pixelShiftMedian3->set_inconsistent (false);
            pixelShiftMedian3->block (true);
            pixelShiftMedian3->set_active (false);
            pixelShiftMedian3->block (false);
        } else if (pixelShiftMedian3->getLastActive()) {
            pixelShiftMedian3->set_inconsistent (true);
        }

        pixelShiftMedian3->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftMedian3, pixelShiftMedian3->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
#endif
void BayerProcess::pixelShiftGreenChanged ()
{
    if (batchMode) {
        if (pixelShiftGreen->get_inconsistent()) {
            pixelShiftGreen->set_inconsistent (false);
            pixelShiftGreen->block (true);
            pixelShiftGreen->set_active (false);
            pixelShiftGreen->block (false);
        } else if (pixelShiftGreen->getLastActive()) {
            pixelShiftGreen->set_inconsistent (true);
        }

        pixelShiftGreen->setLastActive();
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
            pixelShiftBlur->block (true);
            pixelShiftBlur->set_active (false);
            pixelShiftBlur->block (false);
        } else if (pixelShiftBlur->getLastActive()) {
            pixelShiftBlur->set_inconsistent (true);
        }

        pixelShiftBlur->setLastActive();
    }

    pixelShiftSmooth->set_sensitive(pixelShiftBlur->get_active ());
    pixelShiftSigma->set_sensitive(pixelShiftBlur->get_active ());
    if (listener) {
        listener->panelChanged (EvPixelShiftBlur, pixelShiftBlur->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
#ifdef PIXELSHIFTDEV
void BayerProcess::pixelShiftExp0Changed ()
{
    if (batchMode) {
        if (pixelShiftExp0->get_inconsistent()) {
            pixelShiftExp0->set_inconsistent (false);
            pixelShiftExp0->block (true);
            pixelShiftExp0->set_active (false);
            pixelShiftExp0->block (false);
        } else if (pixelShiftExp0->getLastActive()) {
            pixelShiftExp0->set_inconsistent (true);
        }

        pixelShiftExp0->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftExp0, pixelShiftExp0->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
#endif
void BayerProcess::pixelShiftLmmseChanged ()
{
    if (batchMode) {
        if (pixelShiftLmmse->get_inconsistent()) {
            pixelShiftLmmse->set_inconsistent (false);
            pixelShiftLmmse->block (true);
            pixelShiftLmmse->set_active (false);
            pixelShiftLmmse->block (false);
        } else if (pixelShiftLmmse->getLastActive()) {
            pixelShiftLmmse->set_inconsistent (true);
        }

        pixelShiftLmmse->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftLmmse, pixelShiftLmmse->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftEqualBrightChanged ()
{
    if (batchMode) {
        if (pixelShiftEqualBright->get_inconsistent()) {
            pixelShiftEqualBright->set_inconsistent (false);
            pixelShiftEqualBright->block (true);
            pixelShiftEqualBright->set_active (false);
            pixelShiftEqualBright->block (false);
        } else if (pixelShiftEqualBright->getLastActive()) {
            pixelShiftEqualBright->set_inconsistent (true);
        }

        pixelShiftEqualBright->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftEqualBright, pixelShiftEqualBright->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void BayerProcess::pixelShiftNonGreenCrossChanged ()
{
    if (batchMode) {
        if (pixelShiftNonGreenCross->get_inconsistent()) {
            pixelShiftNonGreenCross->set_inconsistent (false);
            pixelShiftNonGreenCross->block (true);
            pixelShiftNonGreenCross->set_active (false);
            pixelShiftNonGreenCross->block (false);
        } else if (pixelShiftNonGreenCross->getLastActive()) {
            pixelShiftNonGreenCross->set_inconsistent (true);
        }

        pixelShiftNonGreenCross->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftNonGreenCross, pixelShiftNonGreenCross->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
#ifdef PIXELSHIFTDEV
void BayerProcess::pixelShiftNonGreenCross2Changed ()
{
    if (batchMode) {
        if (pixelShiftNonGreenCross2->get_inconsistent()) {
            pixelShiftNonGreenCross2->set_inconsistent (false);
            pixelShiftNonGreenCross2->block (true);
            pixelShiftNonGreenCross2->set_active (false);
            pixelShiftNonGreenCross2->block (false);
        } else if (pixelShiftNonGreenCross2->getLastActive()) {
            pixelShiftNonGreenCross2->set_inconsistent (true);
        }

        pixelShiftNonGreenCross2->setLastActive();
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
            pixelShiftNonGreenAmaze->block (true);
            pixelShiftNonGreenAmaze->set_active (false);
            pixelShiftNonGreenAmaze->block (false);
        } else if (pixelShiftNonGreenAmaze->getLastActive()) {
            pixelShiftNonGreenAmaze->set_inconsistent (true);
        }

        pixelShiftNonGreenAmaze->setLastActive();
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftNonGreenAmaze, pixelShiftNonGreenAmaze->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
#endif