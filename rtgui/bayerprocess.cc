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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "bayerprocess.h"

#include "eventmapper.h"
#include "guiutils.h"
#include "options.h"

#include "rtengine/procparams.h"
#include "rtengine/utils.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring BayerProcess::TOOL_NAME = "bayerprocess";

BayerProcess::BayerProcess () :
    FoldableToolPanel(this, TOOL_NAME, M("TP_RAW_LABEL"), options.prevdemo != PD_Sidecar),
    oldMethod(-1)
{

    auto m = ProcEventMapper::getInstance();
    EvDemosaicBorder = m->newEvent(DARKFRAME, "HISTORY_MSG_RAW_BORDER");
    EvDemosaicContrast = m->newEvent(DEMOSAIC, "HISTORY_MSG_DUALDEMOSAIC_CONTRAST");
    EvDemosaicAutoContrast = m->newEvent(DEMOSAIC, "HISTORY_MSG_DUALDEMOSAIC_AUTO_CONTRAST");
    EvDemosaicPixelshiftDemosaicMethod = m->newEvent(DEMOSAIC, "HISTORY_MSG_PIXELSHIFT_DEMOSAIC");
    EvPixelshiftAverage = m->newEvent(DEMOSAIC, "HISTORY_MSG_PIXELSHIFT_AVERAGE");

    Gtk::Box* hb1 = Gtk::manage (new Gtk::Box ());
    hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_DMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    method = Gtk::manage (new MyComboBoxText ());

    for(const auto method_string : procparams::RAWParams::BayerSensor::getMethodStrings()) {
        method->append(M("TP_RAW_" + Glib::ustring(method_string).uppercase()));
    }

    method->set_active(0);
    hb1->set_tooltip_markup (M("TP_RAW_DMETHOD_TOOLTIP"));

    hb1->pack_end (*method, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *hb1, Gtk::PACK_SHRINK, 4);

    dualDemosaicOptions = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    dualDemosaicContrast = Gtk::manage(new Adjuster(M("TP_RAW_DUALDEMOSAICCONTRAST"), 0, 100, 1, 20));
    dualDemosaicContrast->setAdjusterListener(this);
    dualDemosaicContrast->addAutoButton(M("TP_RAW_DUALDEMOSAICAUTOCONTRAST_TOOLTIP"));
    dualDemosaicContrast->setAutoValue(true);
    dualDemosaicContrast->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    dualDemosaicContrast->show();
    dualDemosaicOptions->pack_start(*dualDemosaicContrast);
    pack_start( *dualDemosaicOptions, Gtk::PACK_SHRINK, 4);

    borderbox = Gtk::manage(new Gtk::Box());
    border = Gtk::manage(new Adjuster(M("TP_RAW_BORDER"), 0, 16, 1, 4));
    border->setAdjusterListener (this);

    border->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    border->show();
    borderbox->pack_start(*border);
    pack_start(*borderbox, Gtk::PACK_SHRINK, 4);

    imageNumberBox = Gtk::manage (new Gtk::Box ());
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

    pack_start( *Gtk::manage( new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)), Gtk::PACK_SHRINK, 0 );
    ccSteps = Gtk::manage (new Adjuster (M("TP_RAW_FALSECOLOR"), 0, 5, 1, 0 ));
    ccSteps->setAdjusterListener (this);

    ccSteps->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    ccSteps->show();
    pack_start( *ccSteps, Gtk::PACK_SHRINK, 4);

    dcbOptions = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    dcbIterations = Gtk::manage (new Adjuster (M("TP_RAW_DCBITERATIONS"), 0, 5, 1, 2));
    dcbIterations->setAdjusterListener (this);

    dcbIterations->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    dcbIterations->show();
    dcbEnhance = Gtk::manage (new CheckBox(M("TP_RAW_DCBENHANCE"), multiImage));
    dcbEnhance->setCheckBoxListener (this);
    dcbOptions->pack_start(*dcbIterations);
    dcbOptions->pack_start(*dcbEnhance);
    pack_start( *dcbOptions, Gtk::PACK_SHRINK, 4);

    lmmseOptions = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    lmmseIterations = Gtk::manage (new Adjuster (M("TP_RAW_LMMSEITERATIONS"), 0, 6, 1, 2));
    lmmseIterations->setAdjusterListener (this);
    lmmseIterations->set_tooltip_markup (M("TP_RAW_LMMSE_TOOLTIP"));

    lmmseIterations->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    lmmseIterations->show();
    lmmseOptions->pack_start(*lmmseIterations);
    pack_start( *lmmseOptions, Gtk::PACK_SHRINK, 4);

    // --------------------  PixelShift  ----------------------


    pixelShiftFrame = Gtk::manage(new Gtk::Frame(M("TP_RAW_PIXELSHIFT")));

    Gtk::Box *pixelShiftMainVBox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    pixelShiftMainVBox->set_border_width(0);

    pixelShiftEqualBright = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTEQUALBRIGHT"), multiImage));
    pixelShiftEqualBright->setCheckBoxListener (this);
    pixelShiftEqualBright->set_tooltip_text (M("TP_RAW_PIXELSHIFTEQUALBRIGHT_TOOLTIP"));
    pixelShiftMainVBox->pack_start(*pixelShiftEqualBright);

    pixelShiftEqualBrightChannel = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTEQUALBRIGHTCHANNEL"), multiImage));
    pixelShiftEqualBrightChannel->setCheckBoxListener (this);
    pixelShiftEqualBrightChannel->set_tooltip_text (M("TP_RAW_PIXELSHIFTEQUALBRIGHTCHANNEL_TOOLTIP"));
    pixelShiftMainVBox->pack_start(*pixelShiftEqualBrightChannel);

    Gtk::Box* hb3 = Gtk::manage (new Gtk::Box ());
    hb3->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_PIXELSHIFTMOTIONMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    pixelShiftMotionMethod = Gtk::manage (new MyComboBoxText ());
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_OFF"));
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_AUTO"));
    pixelShiftMotionMethod->append(M("TP_RAW_PIXELSHIFTMM_CUSTOM"));
    pixelShiftMotionMethod->set_active(toUnderlying(RAWParams::BayerSensor::PSMotionCorrectionMethod::AUTO));
    pixelShiftMotionMethod->show();
    hb3->pack_start(*pixelShiftMotionMethod);
    pixelShiftMainVBox->pack_start(*hb3);

    pixelShiftOptions = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    pixelShiftOptions->set_border_width(0);

    pixelShiftShowMotion = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTSHOWMOTION"), multiImage));
    pixelShiftShowMotion->setCheckBoxListener (this);
    pixelShiftShowMotion->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTION_TOOLTIP"));
    pixelShiftMainVBox->pack_start(*pixelShiftShowMotion);

    pixelShiftShowMotionMaskOnly = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY"), multiImage));
    pixelShiftShowMotionMaskOnly->setCheckBoxListener (this);
    pixelShiftShowMotionMaskOnly->set_tooltip_text (M("TP_RAW_PIXELSHIFTSHOWMOTIONMASKONLY_TOOLTIP"));
    pixelShiftMainVBox->pack_start(*pixelShiftShowMotionMaskOnly);


    Gtk::Box* hb4 = Gtk::manage (new Gtk::Box ());
    hb4->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_PIXELSHIFTDMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    pixelShiftDemosaicMethod = Gtk::manage (new MyComboBoxText ());
    for(const auto method_string : procparams::RAWParams::BayerSensor::getPSDemosaicMethodStrings()) {
        pixelShiftDemosaicMethod->append(M("TP_RAW_" + Glib::ustring(method_string).uppercase()));
    }
    pixelShiftDemosaicMethod->set_active(0);
    hb4->pack_start(*pixelShiftDemosaicMethod);
    pixelShiftOptions->pack_start(*hb4);

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

    pixelShiftSigma->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    pixelShiftSigma->show();
    pixelShiftOptions->pack_start(*pixelShiftSigma);


    pixelShiftSmooth = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTSMOOTH"), 0, 1, 0.05, 0.7));
    pixelShiftSmooth->set_tooltip_text (M("TP_RAW_PIXELSHIFTSMOOTH_TOOLTIP"));
    pixelShiftSmooth->setAdjusterListener (this);

    pixelShiftSmooth->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    pixelShiftSmooth->show();
    pixelShiftOptions->pack_start(*pixelShiftSmooth);

    pixelShiftEperIso = Gtk::manage (new Adjuster (M("TP_RAW_PIXELSHIFTEPERISO"), -5.0, 5.0, 0.05, 0.0));
    pixelShiftEperIso->set_tooltip_text(M("TP_RAW_PIXELSHIFTEPERISO_TOOLTIP"));
    pixelShiftEperIso->setAdjusterListener (this);

    pixelShiftEperIso->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    pixelShiftEperIso->show();
    pixelShiftOptions->pack_start(*pixelShiftEperIso);


    pixelShiftMedian = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTMEDIAN"), multiImage));
    pixelShiftMedian->setCheckBoxListener (this);
    pixelShiftMedian->set_tooltip_text (M("TP_RAW_PIXELSHIFTMEDIAN_TOOLTIP"));

    pixelShiftAverage = Gtk::manage (new CheckBox(M("TP_RAW_PIXELSHIFTAVERAGE"), multiImage));
    pixelShiftAverage->setCheckBoxListener (this);
    pixelShiftAverage->set_tooltip_text (M("TP_RAW_PIXELSHIFTAVERAGE_TOOLTIP"));

    pixelShiftOptions->pack_start(*pixelShiftMedian);
    pixelShiftOptions->pack_start(*pixelShiftAverage);

    pixelShiftMainVBox->pack_start(*pixelShiftOptions);
    pixelShiftFrame->add(*pixelShiftMainVBox);
    pixelShiftOptions->hide();

    pack_start( *pixelShiftFrame, Gtk::PACK_SHRINK, 4);

    method->connect(method->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::methodChanged) ));
    imageNumber->connect(imageNumber->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::imageNumberChanged) ));
    pixelShiftMotionMethod->connect(pixelShiftMotionMethod->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::pixelShiftMotionMethodChanged) ));
    pixelShiftDemosaicMethod->connect(pixelShiftDemosaicMethod->signal_changed().connect( sigc::mem_fun(*this, &BayerProcess::pixelShiftDemosaicMethodChanged) ));

}

BayerProcess::~BayerProcess ()
{
    idle_register.destroy();
}

void BayerProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    method->block (true);
    imageNumber->block (true);
    pixelShiftDemosaicMethod->block(true);
    //allEnhconn.block (true);

    border->setValue(pp->raw.bayersensor.border);
    imageNumber->set_active(pp->raw.bayersensor.imageNum);

    for (size_t i = 0; i < procparams::RAWParams::BayerSensor::getMethodStrings().size(); ++i) {
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodStrings()[i]) {
            method->set_active(i);
            oldMethod = i;
            break;
        }
    }
    for (size_t i = 0; i < procparams::RAWParams::BayerSensor::getPSDemosaicMethodStrings().size(); ++i) {
        if (pp->raw.bayersensor.pixelShiftDemosaicMethod == procparams::RAWParams::BayerSensor::getPSDemosaicMethodStrings()[i]) {
            pixelShiftDemosaicMethod->set_active(i);
            break;
        }
    }


    dcbIterations->setValue (pp->raw.bayersensor.dcb_iterations);
    dcbEnhance->setValue (pp->raw.bayersensor.dcb_enhance);
    pixelShiftShowMotion->setValue (pp->raw.bayersensor.pixelShiftShowMotion);
    pixelShiftShowMotionMaskOnly->setValue (pp->raw.bayersensor.pixelShiftShowMotionMaskOnly);
    if (!batchMode) {
        pixelShiftShowMotionMaskOnly->set_sensitive (pp->raw.bayersensor.pixelShiftShowMotion);
    }
    pixelShiftHoleFill->setValue (pp->raw.bayersensor.pixelShiftHoleFill);
    pixelShiftMedian->setValue (pp->raw.bayersensor.pixelShiftMedian);
    pixelShiftAverage->setValue (pp->raw.bayersensor.pixelShiftAverage);
    pixelShiftGreen->setValue (pp->raw.bayersensor.pixelShiftGreen);
    pixelShiftBlur->setValue (pp->raw.bayersensor.pixelShiftBlur);
    pixelShiftSmooth->setValue (pp->raw.bayersensor.pixelShiftSmoothFactor);
    if (!batchMode) {
        pixelShiftSmooth->set_sensitive (pp->raw.bayersensor.pixelShiftBlur);
    }
    pixelShiftEqualBright->setValue (pp->raw.bayersensor.pixelShiftEqualBright);
    pixelShiftEqualBrightChannel->setValue (pp->raw.bayersensor.pixelShiftEqualBrightChannel);
    if (!batchMode) {
        pixelShiftEqualBrightChannel->set_sensitive (pp->raw.bayersensor.pixelShiftEqualBright);
    }
    pixelShiftNonGreenCross->setValue (pp->raw.bayersensor.pixelShiftNonGreenCross);
    ccSteps->setValue (pp->raw.bayersensor.ccSteps);
    lmmseIterations->setValue (pp->raw.bayersensor.lmmse_iterations);
    dualDemosaicContrast->setAutoValue(pp->raw.bayersensor.dualDemosaicAutoContrast);
    dualDemosaicContrast->setValue (pp->raw.bayersensor.dualDemosaicContrast);

    pixelShiftMotionMethod->set_active ((int)pp->raw.bayersensor.pixelShiftMotionCorrectionMethod);
    pixelShiftEperIso->setValue (pp->raw.bayersensor.pixelShiftEperIso);
    pixelShiftSigma->setValue (pp->raw.bayersensor.pixelShiftSigma);
    if (!batchMode) {
        pixelShiftSigma->set_sensitive (pp->raw.bayersensor.pixelShiftBlur);
    }

    if(pedited) {
        border->setEditedState (pedited->raw.bayersensor.border ? Edited : UnEdited);
        ccSteps->setEditedState (pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
        dcbIterations->setEditedState ( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        dcbEnhance->setEdited (pedited->raw.bayersensor.dcbEnhance);
        pixelShiftShowMotion->setEdited (pedited->raw.bayersensor.pixelShiftShowMotion);
        pixelShiftShowMotionMaskOnly->setEdited (pedited->raw.bayersensor.pixelShiftShowMotionMaskOnly);
        pixelShiftHoleFill->setEdited (pedited->raw.bayersensor.pixelShiftHoleFill);
        pixelShiftMedian->setEdited(pedited->raw.bayersensor.pixelShiftMedian);
        pixelShiftAverage->setEdited(pedited->raw.bayersensor.pixelShiftAverage);
        pixelShiftGreen->setEdited (pedited->raw.bayersensor.pixelShiftGreen);
        pixelShiftBlur->setEdited (pedited->raw.bayersensor.pixelShiftBlur);
        pixelShiftSmooth->setEditedState ( pedited->raw.bayersensor.pixelShiftSmooth ? Edited : UnEdited);
        pixelShiftEqualBright->setEdited (pedited->raw.bayersensor.pixelShiftEqualBright);
        pixelShiftEqualBrightChannel->setEdited (pedited->raw.bayersensor.pixelShiftEqualBrightChannel);
        pixelShiftNonGreenCross->setEdited (pedited->raw.bayersensor.pixelShiftNonGreenCross);
        lmmseIterations->setEditedState ( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        dualDemosaicContrast->setEditedState ( pedited->raw.bayersensor.dualDemosaicContrast ? Edited : UnEdited);
        dualDemosaicContrast->setAutoInconsistent   (multiImage && !pedited->raw.bayersensor.dualDemosaicAutoContrast);
        pixelShiftEperIso->setEditedState ( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftSigma->setEditedState ( pedited->raw.bayersensor.pixelShiftSigma ? Edited : UnEdited);

        if(!pedited->raw.bayersensor.method) {
            method->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if(!pedited->raw.bayersensor.imageNum) {
            imageNumber->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if(!pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod) {
            pixelShiftMotionMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if(!pedited->raw.bayersensor.pixelShiftDemosaicMethod) {
            pixelShiftDemosaicMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

    }
    lastAutoContrast = pp->raw.bayersensor.dualDemosaicAutoContrast;

    if (!batchMode) {
        dcbOptions->set_visible(pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCB) || pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBVNG4));
        lmmseOptions->set_visible(pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::LMMSE));
        dualDemosaicOptions->set_visible(pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEVNG4)
                                         || pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEBILINEAR)
                                         || pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBVNG4)
                                         || pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::DCBBILINEAR)
                                         || pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDBILINEAR)
                                         || pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDVNG4));
        if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::PIXELSHIFT)) {
            pixelShiftOptions->set_visible(pp->raw.bayersensor.pixelShiftMotionCorrectionMethod == RAWParams::BayerSensor::PSMotionCorrectionMethod::CUSTOM);
            pixelShiftFrame->show();
        } else {
            pixelShiftFrame->hide();
        }

        // False color suppression is applied to all demozaicing method, so don't hide anything
        /*if (pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::EAHD) ||
              pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::HPHD) ||
              pp->raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::VNG4))
            ccSteps->show();
        else
            ccSteps->hide();*/
    }

    //lastALLen = pp->raw.bayersensor.all_enhance;

    method->block (false);
    imageNumber->block (false);
    pixelShiftDemosaicMethod->block(false);
    //allEnhconn.block (false);

    enableListener ();
}

void BayerProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.bayersensor.ccSteps = ccSteps->getIntValue();
    pp->raw.bayersensor.dcb_iterations = dcbIterations->getIntValue();
    pp->raw.bayersensor.dcb_enhance = dcbEnhance->getLastActive ();
    pp->raw.bayersensor.border = border->getIntValue();
    pp->raw.bayersensor.lmmse_iterations = lmmseIterations->getIntValue();
    pp->raw.bayersensor.dualDemosaicAutoContrast = dualDemosaicContrast->getAutoValue();
    pp->raw.bayersensor.dualDemosaicContrast = dualDemosaicContrast->getValue();
    pp->raw.bayersensor.pixelShiftMotionCorrectionMethod = (RAWParams::BayerSensor::PSMotionCorrectionMethod)pixelShiftMotionMethod->get_active_row_number();
    pp->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getValue();
    pp->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getValue();
    pp->raw.bayersensor.pixelShiftShowMotion = pixelShiftShowMotion->getLastActive ();
    pp->raw.bayersensor.pixelShiftShowMotionMaskOnly = pixelShiftShowMotionMaskOnly->getLastActive ();
    pp->raw.bayersensor.pixelShiftHoleFill = pixelShiftHoleFill->getLastActive ();
    pp->raw.bayersensor.pixelShiftMedian = pixelShiftMedian->getLastActive ();
    pp->raw.bayersensor.pixelShiftAverage = pixelShiftAverage->getLastActive ();
    pp->raw.bayersensor.pixelShiftGreen = pixelShiftGreen->getLastActive ();
    pp->raw.bayersensor.pixelShiftBlur = pixelShiftBlur->getLastActive ();
    pp->raw.bayersensor.pixelShiftSmoothFactor = pixelShiftSmooth->getValue();
    pp->raw.bayersensor.pixelShiftEqualBright = pixelShiftEqualBright->getLastActive ();
    pp->raw.bayersensor.pixelShiftEqualBrightChannel = pixelShiftEqualBrightChannel->getLastActive ();
    pp->raw.bayersensor.pixelShiftNonGreenCross = pixelShiftNonGreenCross->getLastActive ();

    int currentRow = method->get_active_row_number();
    if( currentRow >= 0 && method->get_active_text() != M("GENERAL_UNCHANGED")) {
        pp->raw.bayersensor.method = procparams::RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method(currentRow));
    }

    currentRow = imageNumber->get_active_row_number();
    if (currentRow >= 0 && imageNumber->get_active_text() != M("GENERAL_UNCHANGED")) {
        pp->raw.bayersensor.imageNum = currentRow;
    }

    currentRow = pixelShiftDemosaicMethod->get_active_row_number();
    if( currentRow >= 0 && pixelShiftDemosaicMethod->get_active_text() != M("GENERAL_UNCHANGED")) {
        pp->raw.bayersensor.pixelShiftDemosaicMethod = procparams::RAWParams::BayerSensor::getPSDemosaicMethodString(RAWParams::BayerSensor::PSDemosaicMethod(currentRow));
    }


    if (pedited) {
        pedited->raw.bayersensor.border = border->getEditedState ();
        pedited->raw.bayersensor.ccSteps = ccSteps->getEditedState ();
        pedited->raw.bayersensor.method = method->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.imageNum = imageNumber->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.dcbIterations = dcbIterations->getEditedState ();
        pedited->raw.bayersensor.dcbEnhance = !dcbEnhance->get_inconsistent();
        //pedited->raw.bayersensor.allEnhance = !allEnhance->get_inconsistent();
        pedited->raw.bayersensor.lmmseIterations = lmmseIterations->getEditedState ();
        pedited->raw.bayersensor.dualDemosaicAutoContrast = !dualDemosaicContrast->getAutoInconsistent ();
        pedited->raw.bayersensor.dualDemosaicContrast = dualDemosaicContrast->getEditedState ();
        pedited->raw.bayersensor.pixelShiftMotionCorrectionMethod = pixelShiftMotionMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.pixelShiftDemosaicMethod = pixelShiftDemosaicMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->raw.bayersensor.pixelShiftEperIso = pixelShiftEperIso->getEditedState ();
        pedited->raw.bayersensor.pixelShiftSigma = pixelShiftSigma->getEditedState ();
        pedited->raw.bayersensor.pixelShiftShowMotion = !pixelShiftShowMotion->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftShowMotionMaskOnly = !pixelShiftShowMotionMaskOnly->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftHoleFill = !pixelShiftHoleFill->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftMedian = !pixelShiftMedian->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftAverage = !pixelShiftAverage->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftGreen = !pixelShiftGreen->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftBlur = !pixelShiftBlur->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftSmooth = pixelShiftSmooth->getEditedState();
        pedited->raw.bayersensor.pixelShiftEqualBright = !pixelShiftEqualBright->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftEqualBrightChannel = !pixelShiftEqualBrightChannel->get_inconsistent();
        pedited->raw.bayersensor.pixelShiftNonGreenCross = !pixelShiftNonGreenCross->get_inconsistent();
    }
}

void BayerProcess::setAdjusterBehavior (bool falsecoloradd, bool iteradd, bool dualdemozecontrastadd, bool pssigmaadd, bool pssmoothadd, bool pseperisoadd)
{
    border->setAddMode(false);
    ccSteps->setAddMode(falsecoloradd);
    dcbIterations->setAddMode(iteradd);
    lmmseIterations->setAddMode(iteradd);
    pixelShiftSmooth->setAddMode(pssmoothadd);
    pixelShiftEperIso->setAddMode(pseperisoadd);
    pixelShiftSigma->setAddMode(pssigmaadd);
    dualDemosaicContrast->setAddMode(dualdemozecontrastadd);
}

void BayerProcess::trimValues (rtengine::procparams::ProcParams* pp)
{
    border->trimValue(pp->raw.bayersensor.border);
    ccSteps->trimValue(pp->raw.bayersensor.ccSteps);
    dcbIterations->trimValue(pp->raw.bayersensor.dcb_iterations);
    lmmseIterations->trimValue(pp->raw.bayersensor.lmmse_iterations);
    pixelShiftSmooth->trimValue(pp->raw.bayersensor.pixelShiftSmoothFactor);
    pixelShiftEperIso->trimValue(pp->raw.bayersensor.pixelShiftEperIso);
    pixelShiftSigma->trimValue(pp->raw.bayersensor.pixelShiftSigma);
    dualDemosaicContrast->trimValue(pp->raw.bayersensor.dualDemosaicContrast);
}

void BayerProcess::setBatchMode(bool batchMode)
{
    method->append (M("GENERAL_UNCHANGED"));
    method->set_active_text(M("GENERAL_UNCHANGED")); // No name
    pixelShiftMotionMethod->append (M("GENERAL_UNCHANGED"));
    pixelShiftMotionMethod->set_active_text (M("GENERAL_UNCHANGED"));
    pixelShiftDemosaicMethod->append (M("GENERAL_UNCHANGED"));
    pixelShiftDemosaicMethod->set_active_text(M("GENERAL_UNCHANGED")); // No name
    imageNumber->append (M("GENERAL_UNCHANGED"));
    imageNumber->set_active_text (M("GENERAL_UNCHANGED"));
    ToolPanel::setBatchMode (batchMode);
    border->showEditedCB ();
    ccSteps->showEditedCB ();
    dcbIterations->showEditedCB ();
    lmmseIterations->showEditedCB ();
    dualDemosaicContrast->showEditedCB ();
    pixelShiftEperIso->showEditedCB ();
    pixelShiftSigma->showEditedCB ();
}

void BayerProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    dcbIterations->setDefault( defParams->raw.bayersensor.dcb_iterations);
    lmmseIterations->setDefault( defParams->raw.bayersensor.lmmse_iterations);
    dualDemosaicContrast->setDefault( defParams->raw.bayersensor.dualDemosaicContrast);
    pixelShiftEperIso->setDefault( defParams->raw.bayersensor.pixelShiftEperIso);
    pixelShiftSigma->setDefault( defParams->raw.bayersensor.pixelShiftSigma);
    border->setDefault (defParams->raw.bayersensor.border);
    ccSteps->setDefault (defParams->raw.bayersensor.ccSteps);

    if (pedited) {
        dcbIterations->setDefaultEditedState( pedited->raw.bayersensor.dcbIterations ? Edited : UnEdited);
        lmmseIterations->setDefaultEditedState( pedited->raw.bayersensor.lmmseIterations ? Edited : UnEdited);
        dualDemosaicContrast->setDefaultEditedState( pedited->raw.bayersensor.dualDemosaicContrast ? Edited : UnEdited);
        pixelShiftEperIso->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftEperIso ? Edited : UnEdited);
        pixelShiftSigma->setDefaultEditedState( pedited->raw.bayersensor.pixelShiftSigma ? Edited : UnEdited);
        border->setDefaultEditedState(pedited->raw.bayersensor.border ? Edited : UnEdited);
        ccSteps->setDefaultEditedState(pedited->raw.bayersensor.ccSteps ? Edited : UnEdited);
    } else {
        dcbIterations->setDefaultEditedState( Irrelevant );
        lmmseIterations->setDefaultEditedState( Irrelevant );
        dualDemosaicContrast->setDefaultEditedState( Irrelevant );
        pixelShiftEperIso->setDefaultEditedState( Irrelevant );
        pixelShiftSigma->setDefaultEditedState( Irrelevant );
        border->setDefaultEditedState(Irrelevant);
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
        } else if (a == dualDemosaicContrast) {
            listener->panelChanged (EvDemosaicContrast, a->getTextValue() );
        } else if (a == pixelShiftEperIso) {
            listener->panelChanged (EvPixelShiftEperIso, a->getTextValue() );
        } else if (a == pixelShiftSigma) {
            listener->panelChanged (EvPixelShiftSigma, a->getTextValue() );
        } else if (a == pixelShiftSmooth) {
            listener->panelChanged (EvPixelShiftSmooth, a->getTextValue() );
        } else if (a == border) {
            listener->panelChanged (EvDemosaicBorder, a->getTextValue() );
        }
    }
}

void BayerProcess::methodChanged ()
{
    const int currentSelection = method->get_active_row_number();
    const RAWParams::BayerSensor::Method currentMethod = RAWParams::BayerSensor::Method(currentSelection);

    if (!batchMode) {
        if (currentMethod == procparams::RAWParams::BayerSensor::Method::DCB || currentMethod == procparams::RAWParams::BayerSensor::Method::DCBVNG4) {
            dcbOptions->show();
        } else {
            dcbOptions->hide();
        }

        if (currentMethod == procparams::RAWParams::BayerSensor::Method::LMMSE) {
            lmmseOptions->show();
        } else {
            lmmseOptions->hide();
        }

        if (currentMethod == procparams::RAWParams::BayerSensor::Method::AMAZEVNG4 ||
            currentMethod == procparams::RAWParams::BayerSensor::Method::DCBVNG4 ||
            currentMethod == procparams::RAWParams::BayerSensor::Method::RCDVNG4 ||
            currentMethod == procparams::RAWParams::BayerSensor::Method::AMAZEBILINEAR ||
            currentMethod == procparams::RAWParams::BayerSensor::Method::DCBBILINEAR ||
            currentMethod == procparams::RAWParams::BayerSensor::Method::RCDBILINEAR) {
            dualDemosaicOptions->show();
        } else {
            dualDemosaicOptions->hide();
        }

        if (currentMethod == procparams::RAWParams::BayerSensor::Method::PIXELSHIFT) {
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

    if (listener && method->get_active_row_number() >= 0) {
        listener->panelChanged (
            currentMethod == procparams::RAWParams::BayerSensor::Method::MONO || RAWParams::BayerSensor::Method(oldMethod) == procparams::RAWParams::BayerSensor::Method::MONO
            ? EvDemosaicMethodPreProc
            : EvDemosaicMethod, method->get_active_text());
    }

    oldMethod = currentSelection;

}

void BayerProcess::pixelShiftDemosaicMethodChanged ()
{

    if (listener && pixelShiftDemosaicMethod->get_active_row_number() >= 0) {
        listener->panelChanged (EvDemosaicPixelshiftDemosaicMethod, pixelShiftDemosaicMethod->get_active_text());
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
        if (pixelShiftMedian->getLastActive()) {
            pixelShiftAverage->setValue(false);
        }
        if (listener) {
            listener->panelChanged (EvPixelShiftMedian, pixelShiftMedian->getValueAsStr ());
        }
    } else if (c == pixelShiftAverage) {
        if (pixelShiftAverage->getLastActive()) {
            pixelShiftMedian->setValue(false);
        }
        if (listener) {
            listener->panelChanged (EvPixelshiftAverage, pixelShiftAverage->getValueAsStr ());
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
}

void BayerProcess::adjusterAutoToggled(Adjuster* a)
{
    if (multiImage) {
        if (dualDemosaicContrast->getAutoInconsistent()) {
            dualDemosaicContrast->setAutoInconsistent (false);
            dualDemosaicContrast->setAutoValue (false);
        } else if (lastAutoContrast) {
            dualDemosaicContrast->setAutoInconsistent (true);
        }

        lastAutoContrast = dualDemosaicContrast->getAutoValue();
    }

    if (listener) {

        if (a == dualDemosaicContrast) {
            if (dualDemosaicContrast->getAutoInconsistent()) {
                listener->panelChanged (EvDemosaicAutoContrast, M ("GENERAL_UNCHANGED"));
            } else if (dualDemosaicContrast->getAutoValue()) {
                listener->panelChanged (EvDemosaicAutoContrast, M ("GENERAL_ENABLED"));
            } else {
                listener->panelChanged (EvDemosaicAutoContrast, M ("GENERAL_DISABLED"));
            }
        }
    }
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
        } else if(pixelShiftMotionMethod->get_active_row_number() > 0) {
            pixelShiftOptions->hide();
            pixelShiftShowMotion->show();
            pixelShiftShowMotionMaskOnly->show();
        }
    }
    if (listener && pixelShiftMotionMethod->get_active_row_number() >= 0) {
        listener->panelChanged (EvPixelShiftMotionMethod, pixelShiftMotionMethod->get_active_text());
    }
}

void BayerProcess::FrameCountChanged(int n, int frameNum)
{
    idle_register.add(
        [this, n, frameNum]() -> bool
        {
            imageNumber->block(true);

            imageNumber->remove_all();
            imageNumber->append("1");
            for (int i = 2; i <= n; ++i) {
                std::ostringstream entry;
                entry << i;
                imageNumber->append(entry.str());
            }
            if (n == 2) {
                imageNumber->append(M("TP_RAW_IMAGENUM_SN"));
            }
            imageNumber->set_active(std::min(frameNum, n == 2 ? n : n - 1));
            if (n == 1) {
                imageNumberBox->hide();
            } else {
                imageNumberBox->show();
            }
            imageNumber->block (false);
            return false;
        }
    );
}

void BayerProcess::autoContrastChanged (double autoContrast)
{
    idle_register.add(
        [this, autoContrast]() -> bool
        {
            disableListener();
            dualDemosaicContrast->setValue(autoContrast);
            enableListener();
            return false;
        }
    );
}
