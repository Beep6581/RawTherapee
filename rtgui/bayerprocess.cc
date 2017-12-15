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

    for(const auto method_string : procparams::RAWParams::BayerSensor::getMethodStrings()) {
        method->append(M("TP_RAW_" + Glib::ustring(method_string).uppercase()));
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
    dcbEnhance = Gtk::manage (new CheckBox(M("TP_RAW_DCBENHANCE"), multiImage));
    dcbEnhance->setCheckBoxListener (this);
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

    pixelShiftEqualBright = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTEQUALBRIGHT"), multiImage));
    pixelShiftEqualBright->setCheckBoxListener (this);
    pixelShiftEqualBright->set_tooltip_text (M("TP_RAW_PIXELSHIFTEQUALBRIGHT_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftEqualBright);

    pixelShiftEqualBrightChannel = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTEQUALBRIGHTCHANNEL"), multiImage));
    pixelShiftEqualBrightChannel->setCheckBoxListener (this);
    pixelShiftEqualBrightChannel->set_tooltip_text (M("TP_RAW_PIXELSHIFTEQUALBRIGHTCHANNEL_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftEqualBrightChannel);

    Gtk::HBox* hb3 = Gtk::manage (new Gtk::HBox ());
    hb3->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_PIXELSHIFTMOTIONMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    pixelShiftMotionMethod = Gtk::manage (new MyComboBoxText ());
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_OFF"));
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_AUTO"));
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_CUSTOM"));
    pixelShiftMotionMethod->set_active(toUnderlying(RAWParams::BayerSensor::PSMotionCorrectionMethod::AUTO));
    pixelShiftMotionMethod->show();
    hb3->pack_start(*pixelShiftMotionMethod);
    pixelShiftFrame->pack_start(*hb3);

    pixelShiftOptions = Gtk::manage (new Gtk::VBox ());
    pixelShiftOptions->set_border_width(0);

    pixelShiftShowMotion = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTSHOWMOTION"), multiImage));
    pixelShiftShowMotion->setCheckBoxListener (this);
    pixelShiftShowMotion->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTION_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftShowMotion);

    pixelShiftShowMotionMaskOnly = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY"), multiImage));
    pixelShiftShowMotionMaskOnly->setCheckBoxListener (this);
    pixelShiftShowMotionMaskOnly->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY_TOOLTIP"));
    pixelShiftFrame->pack_start(*pixelShiftShowMotionMaskOnly);

#ifdef PIXELSHIFTDEV
    pixelShiftAutomatic = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTADAPTIVE"), multiImage));
    pixelShiftAutomatic->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftAutomatic);
#endif
    pixelShiftGreen = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTGREEN"), multiImage));
    pixelShiftGreen->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftGreen);

    pixelShiftNonGreenCross = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTNONGREENCROSS"), multiImage));
    pixelShiftNonGreenCross->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftNonGreenCross);

    pixelShiftHoleFill = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTHOLEFILL"), multiImage));
    pixelShiftHoleFill->setCheckBoxListener (this);
    pixelShiftHoleFill->set_tooltip_text (M("TP_RAW_PIXELSHIFTHOLEFILL_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftHoleFill);

    pixelShiftBlur = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTBLUR"), multiImage));
    pixelShiftBlur->setCheckBoxListener (this);
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


    pixelShiftMedian = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTMEDIAN"), multiImage));
    pixelShiftMedian->setCheckBoxListener (this);
    pixelShiftMedian->set_tooltip_text (M("TP_RAW_PIXELSHIFTMEDIAN_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftMedian);


#ifdef PIXELSHIFTDEV
    pixelShiftMedian3 = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTMEDIAN3"), multiImage));
    pixelShiftMedian3->setCheckBoxListener (this);
    pixelShiftMedian3->set_tooltip_text (M("TP_RAW_PIXELSHIFTMEDIAN3_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftMedian3);

    pixelShiftNonGreenCross2 = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTNONGREENCROSS2"), multiImage));
    pixelShiftNonGreenCross2->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftNonGreenCross2);

    pixelShiftNonGreenAmaze = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTNONGREENAMAZE"), multiImage));
    pixelShiftNonGreenAmaze->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftNonGreenAmaze);

    pixelShiftNonGreenHorizontal = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTNONGREENHORIZONTAL"), multiImage));
    pixelShiftNonGreenHorizontal->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftNonGreenHorizontal);

    pixelShiftNonGreenVertical = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTNONGREENVERTICAL"), multiImage));
    pixelShiftNonGreenVertical->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftNonGreenVertical);

    pixelShiftExp0 = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTEXP0"), multiImage));
    pixelShiftExp0->setCheckBoxListener (this);
    pixelShiftOptions->pack_start(*pixelShiftExp0);
#endif
    pixelShiftLmmse = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTLMMSE"), multiImage));
    pixelShiftLmmse->setCheckBoxListener (this);
    pixelShiftLmmse->set_tooltip_text (M("TP_RAW_PIXELSHIFTLMMSE_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftLmmse);

    pixelShiftOneGreen = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTONEGREEN"), multiImage));
    pixelShiftOneGreen->setCheckBoxListener (this);
    pixelShiftOneGreen->set_tooltip_text (M("TP_RAW_PIXELSHIFTONEGREEN_TOOLTIP"));
    pixelShiftOptions->pack_start(*pixelShiftOneGreen);

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
    pixelShiftMotionMethod->connect(pixelShiftMotionMethod->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::pixelShiftMotionMethodChanged) ));
#ifdef PIXELSHIFTDEV
    pixelShiftMotionCorrection->connect(pixelShiftMotionCorrection->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::psMotionCorrectionChanged) ));
#endif
}


void BayerProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    method->block (true);
    imageNumber->block (true);
    //allEnhconn.block (true);
#ifdef PIXELSHIFTDEV
    pixelShiftMotionCorrection->block (true);
#endif

    method->set_active(std::numeric_limits<int>::max());
    imageNumber->set_active(pp->raw.bayersensor.imageNum);

    for (size_t i = 0; i < procparams::RAWParams::BayerSensor::getMethodStrings().size(); ++i) {
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodStrings()[i]) {
            method->set_active(i);
            oldMethod = i;
            break;
        }
    }

    //allEnhance->set_active(pp->raw.bayersensor.all_enhance);

    dcbIterations->setValue (pp->raw.bayersensor.dcb_iterations);
    dcbEnhance->setValue (pp->raw.bayersensor.dcb_enhance);
    pixelShiftShowMotion->setValue (pp->raw.bayersensor.pixelShiftShowMotion);
    if (!batchMode) {
        pixelShiftShowMotionMaskOnly->set_sensitive (pp->raw.bayersensor.pixelShiftShowMotion);
    }
    pixelShiftShowMotionMaskOnly->setValue (pp->raw.bayersensor.pixelShiftShowMotionMaskOnly);
    pixelShiftHoleFill->setValue (pp->raw.bayersensor.pixelShiftHoleFill);
    pixelShiftMedian->setValue (pp->raw.bayersensor.pixelShiftMedian);
    pixelShiftGreen->setValue (pp->raw.bayersensor.pixelShiftGreen);
    pixelShiftBlur->setValue (pp->raw.bayersensor.pixelShiftBlur);
    if (!batchMode) {
        pixelShiftSmooth->set_sensitive (pp->raw.bayersensor.pixelShiftBlur);
    }
    pixelShiftSmooth->setValue (pp->raw.bayersensor.pixelShiftSmoothFactor);
    pixelShiftLmmse->setValue (pp->raw.bayersensor.pixelShiftLmmse);
    pixelShiftOneGreen->setValue (pp->raw.bayersensor.pixelShiftOneGreen);
    pixelShiftEqualBright->setValue (pp->raw.bayersensor.pixelShiftEqualBright);
    pixelShiftEqualBrightChannel->set_sensitive (pp->raw.bayersensor.pixelShiftEqualBright);
    pixelShiftEqualBrightChannel->setValue (pp->raw.bayersensor.pixelShiftEqualBrightChannel);
    pixelShiftNonGreenCross->setValue (pp->raw.bayersensor.pixelShiftNonGreenCross);
    ccSteps->setValue (pp->raw.bayersensor.ccSteps);
    lmmseIterations->setValue (pp->raw.bayersensor.lmmse_iterations);
    pixelShiftMotionMethod->set_active ((int)pp->raw.bayersensor.pixelShiftMotionCorrectionMethod);
    pixelShiftEperIso->setValue (pp->raw.bayersensor.pixelShiftEperIso);
    pixelShiftSigma->setValue (pp->raw.bayersensor.pixelShiftSigma);
    if (!batchMode) {
        pixelShiftSigma->set_sensitive (pp->raw.bayersensor.pixelShiftBlur);
    }
#ifdef PIXELSHIFTDEV
    pixelShiftStddevFactorGreen->setValue (pp->raw.bayersensor.pixelShiftStddevFactorGreen);
    pixelShiftStddevFactorRed->setValue (pp->raw.bayersensor.pixelShiftStddevFactorRed);
    pixelShiftStddevFactorBlue->setValue (pp->raw.bayersensor.pixelShiftStddevFactorBlue);
    pixelShiftSum->setValue (pp->raw.bayersensor.pixelShiftSum);
    pixelShiftMedian3->setValue (pp->raw.bayersensor.pixelShiftMedian3);
    if (!batchMode) {
        pixelShiftMedian3->set_sensitive (pixelShiftMedian->getValue() != CheckValue::off);
    }
    pixelShiftAutomatic->setValue (pp->raw.bayersensor.pixelShiftAutomatic);
    pixelShiftNonGreenHorizontal->setValue (pp->raw.bayersensor.pixelShiftNonGreenHorizontal);
    pixelShiftNonGreenVertical->setValue (pp->raw.bayersensor.pixelShiftNonGreenVertical);
    pixelShiftExp0->setValue (pp->raw.bayersensor.pixelShiftExp0);
    pixelShiftNonGreenCross2->setValue (pp->raw.bayersensor.pixelShiftNonGreenCross2);
    pixelShiftNonGreenAmaze->setValue (pp->raw.bayersensor.pixelShiftNonGreenAmaze);
    pixelShiftMotion->setValue (pp->raw.bayersensor.pixelShiftMotion);
    pixelShiftMotionCorrection->set_active  ((int)pp->raw.bayersensor.pixelShiftMotionCorrection);
    if (!batchMode) {
        pixelShiftHoleFill->set_sensitive (pixelShiftAutomatic->getValue () != CheckValue::off && pixelShiftMotionCorrection->get_active_row_number() == 5);
        pixelShiftBlur->set_sensitive(pixelShiftAutomatic->getValue () != CheckValue::off && pixelShiftMotionCorrection->get_active_row_number() == 5);
        pixelShiftSmooth->set_sensitive(pixelShiftAutomatic->getValue () != CheckValue::off && pixelShiftMotionCorrection->get_active_row_number() == 5 && pixelShiftBlur->getValue() != CheckValue::off);
    }
    pixelShiftNreadIso->setValue (pp->raw.bayersensor.pixelShiftNreadIso);
    pixelShiftPrnu->setValue (pp->raw.bayersensor.pixelShiftPrnu);
    pixelShiftRedBlueWeight->setValue (pp->raw.bayersensor.pixelShiftRedBlueWeight);
#endif

    if(pedited) {
        ccSteps->setEditedState (pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
        dcbIterations->setEditedState ( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        dcbEnhance->setEdited (pedited->raw.bayersensor.dcbEnhance);
        pixelShiftShowMotion->setEdited (pedited->raw.bayersensor.pixelShiftShowMotion);
        pixelShiftShowMotionMaskOnly->setEdited (pedited->raw.bayersensor.pixelShiftShowMotionMaskOnly);
        pixelShiftHoleFill->setEdited (pedited->raw.bayersensor.pixelShiftHoleFill);
        pixelShiftMedian->setEdited(pedited->raw.bayersensor.pixelShiftMedian);
        pixelShiftGreen->setEdited (pedited->raw.bayersensor.pixelShiftGreen);
        pixelShiftBlur->setEdited (pedited->raw.bayersensor.pixelShiftBlur);
        pixelShiftSmooth->setEditedState ( pedited->raw.bayersensor.pixelShiftSmooth ? Edited : UnEdited);
        pixelShiftLmmse->setEdited (pedited->raw.bayersensor.pixelShiftLmmse);
        pixelShiftOneGreen->setEdited (pedited->raw.bayersensor.pixelShiftOneGreen);
        pixelShiftEqualBright->setEdited (pedited->raw.bayersensor.pixelShiftEqualBright);
        pixelShiftEqualBrightChannel->setEdited (pedited->raw.bayersensor.pixelShiftEqualBrightChannel);
        pixelShiftNonGreenCross->setEdited (pedited->raw.bayersensor.pixelShiftNonGreenCross);
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
        pixelShiftAutomatic->setEdited (pedited->raw.bayersensor.pixelShiftAutomatic);
        pixelShiftNonGreenHorizontal->setEdited (pedited->raw.bayersensor.pixelShiftNonGreenHorizontal);
        pixelShiftNonGreenVertical->setEdited (pedited->raw.bayersensor.pixelShiftNonGreenVertical);
        pixelShiftMedian3->setEdited (pedited->raw.bayersensor.pixelShiftMedian3);
        pixelShiftExp0->setEdited (pedited->raw.bayersensor.pixelShiftExp0);
        pixelShiftNonGreenCross2->setEdited (pedited->raw.bayersensor.pixelShiftNonGreenCross2);
        pixelShiftNonGreenAmaze->setEdited (pedited->raw.bayersensor.pixelShiftNonGreenAmaze);
        pixelShiftMotion->setEditedState ( pedited->raw.bayersensor.pixelShiftMotion ? Edited : UnEdited);
        pixelShiftRedBlueWeight->setEditedState ( pedited->raw.bayersensor.pixelShiftRedBlueWeight ? Edited : UnEdited);
#endif

        if(!pedited->raw.bayersensor.method) {
            method->set_active(std::numeric_limits<int>::max()); // No name
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

    if (!batchMode) {
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCB) ||
                method->get_active_row_number() == std::numeric_limits<int>::max()) {
            dcbOptions->show();
        } else {
            dcbOptions->hide();
        }
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::LMMSE) ||
                method->get_active_row_number() == std::numeric_limits<int>::max()) {
            lmmseOptions->show();
        } else {
            lmmseOptions->hide();
        }
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::PIXELSHIFT) ||
                method->get_active_row_number() == std::numeric_limits<int>::max()) {
            if(pp->raw.bayersensor.pixelShiftMotionCorrectionMethod == RAWParams::BayerSensor::PSMotionCorrectionMethod::CUSTOM) {
                pixelShiftOptions->show();
            } else {
                pixelShiftOptions->hide();
            }
            pixelShiftFrame->show();
        } else {
            pixelShiftFrame->hide();
        }

        // Flase color suppression is applied to all demozaicing method, so don't hide anything
        /*if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::EAHD) ||
              pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::HPHD) ||
              pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::VNG4))
            ccSteps->show();
        else
            ccSteps->hide();*/
    }

    //lastALLen = pp->raw.bayersensor.all_enhance;

    method->block (false);
#ifdef PIXELSHIFTDEV
    pixelShiftMotionCorrection->block (false);
#endif
    imageNumber->block (false);
    //allEnhconn.block (false);

    enableListener ();
}

void BayerProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.bayersensor.ccSteps = ccSteps->getIntValue();
    pp->raw.bayersensor.dcb_iterations = dcbIterations->getIntValue();
    pp->raw.bayersensor.dcb_enhance = dcbEnhance->getLastActive ();
    //pp->raw.bayersensor.all_enhance = allEnhance->getLastActive ();
    pp->raw.bayersensor.lmmse_iterations = lmmseIterations->getIntValue();
    pp->raw.bayersensor.pixelShiftMotionCorrectionMethod = (RAWParams::BayerSensor::PSMotionCorrectionMethod)pixelShiftMotionMethod->get_active_row_number();
    pp->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getValue();
    pp->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getValue();
    pp->raw.bayersensor.pixelShiftShowMotion = pixelShiftShowMotion->getLastActive ();
    pp->raw.bayersensor.pixelShiftShowMotionMaskOnly = pixelShiftShowMotionMaskOnly->getLastActive ();
    pp->raw.bayersensor.pixelShiftHoleFill = pixelShiftHoleFill->getLastActive ();
    pp->raw.bayersensor.pixelShiftMedian = pixelShiftMedian->getLastActive ();
    pp->raw.bayersensor.pixelShiftGreen = pixelShiftGreen->getLastActive ();
    pp->raw.bayersensor.pixelShiftBlur = pixelShiftBlur->getLastActive ();
    pp->raw.bayersensor.pixelShiftSmoothFactor = pixelShiftSmooth->getValue();
    pp->raw.bayersensor.pixelShiftLmmse = pixelShiftLmmse->getLastActive ();
    pp->raw.bayersensor.pixelShiftOneGreen = pixelShiftOneGreen->getLastActive ();
    pp->raw.bayersensor.pixelShiftEqualBright = pixelShiftEqualBright->getLastActive ();
    pp->raw.bayersensor.pixelShiftEqualBrightChannel = pixelShiftEqualBrightChannel->getLastActive ();
    pp->raw.bayersensor.pixelShiftNonGreenCross = pixelShiftNonGreenCross->getLastActive ();
#ifdef PIXELSHIFTDEV
    pp->raw.bayersensor.pixelShiftStddevFactorGreen = pixelShiftStddevFactorGreen->getValue();
    pp->raw.bayersensor.pixelShiftStddevFactorRed = pixelShiftStddevFactorRed->getValue();
    pp->raw.bayersensor.pixelShiftStddevFactorBlue = pixelShiftStddevFactorBlue->getValue();
    pp->raw.bayersensor.pixelShiftSum = pixelShiftSum->getValue();
    pp->raw.bayersensor.pixelShiftMedian3 = pixelShiftMedian3->getLastActive ();
    pp->raw.bayersensor.pixelShiftMotion = pixelShiftMotion->getIntValue();
    pp->raw.bayersensor.pixelShiftMotionCorrection = (RAWParams::BayerSensor::ePSMotionCorrection)pixelShiftMotionCorrection->get_active_row_number();
    pp->raw.bayersensor.pixelShiftAutomatic = pixelShiftAutomatic->getLastActive ();
    pp->raw.bayersensor.pixelShiftNonGreenHorizontal = pixelShiftNonGreenHorizontal->getLastActive ();
    pp->raw.bayersensor.pixelShiftNonGreenVertical = pixelShiftNonGreenVertical->getLastActive ();
    pp->raw.bayersensor.pixelShiftExp0 = pixelShiftExp0->getLastActive ();
    pp->raw.bayersensor.pixelShiftNonGreenCross2 = pixelShiftNonGreenCross2->getLastActive ();
    pp->raw.bayersensor.pixelShiftNonGreenAmaze = pixelShiftNonGreenAmaze->getLastActive ();
    pp->raw.bayersensor.pixelShiftNreadIso = pixelShiftNreadIso->getValue();
    pp->raw.bayersensor.pixelShiftPrnu = pixelShiftPrnu->getValue();
    pp->raw.bayersensor.pixelShiftRedBlueWeight = pixelShiftRedBlueWeight->getValue();
#endif

    int currentRow = method->get_active_row_number();
    if( currentRow >= 0 && currentRow < std::numeric_limits<int>::max()) {
        pp->raw.bayersensor.method = procparams::RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method(currentRow));
    }

    currentRow = imageNumber->get_active_row_number();
    if (currentRow < 4) {
        pp->raw.bayersensor.imageNum = currentRow;
    }


    if (pedited) {
        pedited->raw.bayersensor.ccSteps = ccSteps->getEditedState ();
        pedited->raw.bayersensor.method = method->get_active_row_number() != std::numeric_limits<int>::max();
        pedited->raw.bayersensor.imageNum = imageNumber->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.dcbIterations = dcbIterations->getEditedState ();
        pedited->raw.bayersensor.dcbEnhance = !dcbEnhance->get_inconsistent();
        //pedited->raw.bayersensor.allEnhance = !allEnhance->get_inconsistent();
        pedited->raw.bayersensor.lmmseIterations = lmmseIterations->getEditedState ();
        pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod = pixelShiftMotionMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getEditedState ();
        pedited->raw.bayersensor.pixelShiftShowMotion = !pixelShiftShowMotion->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftShowMotionMaskOnly = !pixelShiftShowMotionMaskOnly->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftHoleFill = !pixelShiftHoleFill->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftMedian = !pixelShiftMedian->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftGreen = !pixelShiftGreen->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftBlur = !pixelShiftBlur->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftSmooth = pixelShiftSmooth->getEditedState();
        pedited->raw.bayersensor.pixelShiftLmmse = !pixelShiftLmmse->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftOneGreen = !pixelShiftOneGreen->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftEqualBright = !pixelShiftEqualBright->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftEqualBrightChannel = !pixelShiftEqualBrightChannel->get_inconsistent();
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
    method->set_active(std::numeric_limits<int>::max()); // No name
#ifdef PIXELSHIFTDEV
    pixelShiftMotionCorrection->append (M("GENERAL_UNCHANGED"));
    pixelShiftMotionCorrection->set_active_text (M("GENERAL_UNCHANGED"));
#endif
    pixelShiftMotionMethod->append (M("GENERAL_UNCHANGED"));
    pixelShiftMotionMethod->set_active_text (M("GENERAL_UNCHANGED"));
    imageNumber->append (M("GENERAL_UNCHANGED"));
    imageNumber->set_active_text (M("GENERAL_UNCHANGED"));
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
    if (!batchMode) {
        if(pixelShiftMotionCorrection->get_active_row_number() == 5) {
            pixelShiftBlur->set_sensitive(true);
            pixelShiftHoleFill->set_sensitive(true);
            pixelShiftSmooth->set_sensitive(pixelShiftBlur->getValue() != CheckValue::off);
        } else {
            pixelShiftBlur->set_sensitive(false);
            pixelShiftHoleFill->set_sensitive(false);
            pixelShiftSmooth->set_sensitive(false);
        }
    }

    if (listener) {
        listener->panelChanged (EvPixelShiftMotionCorrection, pixelShiftMotionCorrection->get_active_text());
    }
}
#endif

void BayerProcess::methodChanged ()
{
    const int curSelection = method->get_active_row_number();
    const RAWParams::BayerSensor::Method method = RAWParams::BayerSensor::Method(curSelection);

    if (!batchMode) {
        if (method == procparams::RAWParams::BayerSensor::Method::DCB) {
            dcbOptions->show();
        } else {
            dcbOptions->hide();
        }

        if (method == procparams::RAWParams::BayerSensor::Method::LMMSE) {
            lmmseOptions->show();
        } else {
            lmmseOptions->hide();
        }

        if (method == procparams::RAWParams::BayerSensor::Method::PIXELSHIFT) {
            if(pixelShiftMotionMethod->get_active_row_number() == 2) {
                pixelShiftOptions->show();
            } else {
                pixelShiftOptions->hide();
            }
            pixelShiftFrame->show();
        } else {
            pixelShiftFrame->hide();
        }
    }

    Glib::ustring methodName = "";
    bool ppreq = false;

    if (curSelection >= 0 && curSelection < std::numeric_limits<int>::max()) {
        methodName = procparams::RAWParams::BayerSensor::getMethodString(method);

        if (method == procparams::RAWParams::BayerSensor::Method::MONO || RAWParams::BayerSensor::Method(oldMethod) == procparams::RAWParams::BayerSensor::Method::MONO) {
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

void BayerProcess::checkBoxToggled (CheckBox* c, CheckValue newval)
{
    if (c == dcbEnhance) {
        if (listener) {
            listener->panelChanged (EvDemosaicDCBEnhanced, dcbEnhance->getValueAsStr ());
        }
    } else if (c == pixelShiftShowMotion) {
        if (!batchMode) {
            pixelShiftShowMotionMaskOnly->set_sensitive(newval != CheckValue::off);
        }
        if (listener) {
            listener->panelChanged (EvPixelshiftShowMotion, pixelShiftShowMotion->getValueAsStr ());
        }
    } else if (c == pixelShiftShowMotionMaskOnly) {
        if (listener) {
            listener->panelChanged (EvPixelshiftShowMotionMaskOnly, pixelShiftShowMotionMaskOnly->getValueAsStr ());
        }
    } else if (c == pixelShiftHoleFill) {
        if (listener) {
            listener->panelChanged (EvPixelShiftHoleFill, pixelShiftHoleFill->getValueAsStr ());
        }
    } else if (c == pixelShiftMedian) {
#ifdef PIXELSHIFTDEV
        if (!batchMode) {
            pixelShiftMedian3->set_sensitive(newval != CheckValue::off);
        }
#endif
        if (listener) {
            listener->panelChanged (EvPixelShiftMedian, pixelShiftMedian->getValueAsStr ());
        }
    } else if (c == pixelShiftGreen) {
        if (listener) {
            listener->panelChanged (EvPixelShiftGreen, pixelShiftGreen->getValueAsStr ());
        }
    } else if (c == pixelShiftBlur) {
        if (!batchMode) {
            pixelShiftSmooth->set_sensitive(newval != CheckValue::off);
            pixelShiftSigma->set_sensitive(newval != CheckValue::off);
        }
        if (listener) {
            listener->panelChanged (EvPixelShiftBlur, pixelShiftBlur->getValueAsStr ());
        }
    } else if (c == pixelShiftLmmse) {
        if (listener) {
            listener->panelChanged (EvPixelShiftLmmse, pixelShiftLmmse->getValueAsStr ());
        }
    } else if (c == pixelShiftOneGreen) {
        if (listener) {
            listener->panelChanged (EvPixelShiftOneGreen, pixelShiftOneGreen->getValueAsStr ());
        }
    } else if (c == pixelShiftEqualBright) {
        if (!batchMode) {
            pixelShiftEqualBrightChannel->set_sensitive(newval != CheckValue::off);
        }
        if (listener) {
            listener->panelChanged (EvPixelShiftEqualBright, pixelShiftEqualBright->getValueAsStr ());
        }
    } else if (c == pixelShiftEqualBrightChannel) {
        if (listener) {
            listener->panelChanged (EvPixelShiftEqualBrightChannel, pixelShiftEqualBrightChannel->getValueAsStr ());
        }
    } else if (c == pixelShiftNonGreenCross) {
        if (listener) {
            listener->panelChanged (EvPixelShiftNonGreenCross, pixelShiftNonGreenCross->getValueAsStr ());
        }
    }
#ifdef PIXELSHIFTDEV
    else if (c == pixelShiftAutomatic) {
        if (!batchMode) {
            pixelShiftMotion->set_sensitive(newval != CheckValue::off);
            pixelShiftEperIso->set_sensitive(newval != CheckValue::off);
            pixelShiftNreadIso->set_sensitive(newval != CheckValue::off);
            pixelShiftPrnu->set_sensitive(newval != CheckValue::off);
            pixelShiftSigma->set_sensitive(newval != CheckValue::off);
            pixelShiftSum->set_sensitive(newval != CheckValue::off);
            pixelShiftRedBlueWeight->set_sensitive(newval != CheckValue::off);
            pixelShiftStddevFactorGreen->set_sensitive(newval != CheckValue::off);
            pixelShiftStddevFactorRed->set_sensitive(newval != CheckValue::off);
            pixelShiftStddevFactorBlue->set_sensitive(newval != CheckValue::off);
            pixelShiftNonGreenHorizontal->set_sensitive(newval != CheckValue::off);
            pixelShiftNonGreenVertical->set_sensitive(newval != CheckValue::off);
            pixelShiftHoleFill->set_sensitive(newval != CheckValue::off && pixelShiftMotionCorrection->get_active_row_number() == 5);
            pixelShiftMedian3->set_sensitive(newval != CheckValue::off && pixelShiftMedian->getValue () != CheckValue::off);
            pixelShiftGreen->set_sensitive(newval != CheckValue::off);
            pixelShiftBlur->set_sensitive(newval != CheckValue::off && pixelShiftMotionCorrection->get_active_row_number() == 5);
            pixelShiftSmooth->set_sensitive(newval != CheckValue::off && pixelShiftMotionCorrection->get_active_row_number() == 5 && pixelShiftBlur->getValue () != CheckValue::off);
            pixelShiftExp0->set_sensitive(newval != CheckValue::off);
            pixelShiftNonGreenCross->set_sensitive(newval != CheckValue::off);
            pixelShiftNonGreenCross2->set_sensitive(newval != CheckValue::off);
            pixelShiftNonGreenAmaze->set_sensitive(newval != CheckValue::off);
        }

        if (listener) {
            listener->panelChanged (EvPixelShiftAutomatic, pixelShiftAutomatic->getValueAsStr ());
        }
    } else if (c == pixelShiftNonGreenHorizontal) {
        if (listener) {
            listener->panelChanged (EvPixelShiftNonGreenHorizontal, pixelShiftNonGreenHorizontal->getValueAsStr ());
        }
    } else if (c == pixelShiftNonGreenVertical) {
        if (listener) {
            listener->panelChanged (EvPixelShiftNonGreenVertical, pixelShiftNonGreenVertical->getValueAsStr ());
        }
    } else if (c == pixelShiftMedian3) {
        if (listener) {
            listener->panelChanged (EvPixelShiftMedian3, pixelShiftMedian3->getValueAsStr ());
        }
    } else if (c == pixelShiftExp0) {
        if (listener) {
            listener->panelChanged (EvPixelShiftExp0, pixelShiftExp0->getValueAsStr ());
        }
    } else if (c == pixelShiftNonGreenCross2) {
        if (listener) {
            listener->panelChanged (EvPixelShiftGreenAmaze, pixelShiftNonGreenCross2->getValueAsStr ());
        }
    } else if (c == pixelShiftNonGreenAmaze) {
        if (listener) {
            listener->panelChanged (EvPixelShiftNonGreenAmaze, pixelShiftNonGreenAmaze->getValueAsStr ());
        }
    }
#endif
}

void BayerProcess::pixelShiftMotionMethodChanged ()
{
    if (!batchMode) {
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
    }
    if (listener) {
        listener->panelChanged (EvPixelShiftMotionMethod, pixelShiftMotionMethod->get_active_text());
    }
}

void BayerProcess::FrameCountChanged(int n, int frameNum)
{
    struct Data {
        BayerProcess *me;
        int n;
        int frameNum;
    };
    const auto func = [](gpointer data) -> gboolean {
        Data *d = static_cast<Data *>(data);
        BayerProcess *me = d->me;
        me->imageNumber->block (true);
        int n = d->n;
        int frameNum = d->frameNum;

        me->imageNumber->remove_all();
        me->imageNumber->append("1");
        for(int i = 2; i <= std::min(n, 4); ++i) {
            std::ostringstream entry;
            entry << i;
            me->imageNumber->append(entry.str());
        }
        me->imageNumber->set_active(std::min(frameNum, n - 1));
        if(n == 1) {
            me->imageNumberBox->hide();
        } else {
            me->imageNumberBox->show();
        }
        me->imageNumber->block (false);
        delete d;
        return FALSE;
    };

    idle_register.add(func, new Data { this, n, frameNum });
    
    // GThreadLock lock;
    // imageNumber->block (true);

    // imageNumber->remove_all();
    // imageNumber->append("1");
    // for(int i = 2; i <= std::min(n, 4); ++i) {
    //     std::ostringstream entry;
    //     entry << i;
    //     imageNumber->append(entry.str());
    // }
    // imageNumber->set_active(std::min(frameNum, n - 1));
    // if(n == 1) {
    //     imageNumberBox->hide();
    // } else {
    //     imageNumberBox->show();
    // }
    // imageNumber->block (false);
}
