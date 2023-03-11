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
#include <iomanip>

#include <sigc++/slot.h>

#include "tonecurve.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "eventmapper.h"
#include "options.h"

#include "../rtengine/procparams.h"
#include "../rtengine/utils.h"
#include "editcallbacks.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring ToneCurve::TOOL_NAME = "tonecurve";

ToneCurve::ToneCurve() : FoldableToolPanel(this, TOOL_NAME, M("TP_EXPOSURE_LABEL"))
{
    auto m = ProcEventMapper::getInstance();
    EvHistMatching = m->newEvent(AUTOEXP, "HISTORY_MSG_HISTMATCHING");
    EvHistMatchingBatch = m->newEvent(M_VOID, "HISTORY_MSG_HISTMATCHING");
    EvClampOOG = m->newEvent(DARKFRAME, "HISTORY_MSG_CLAMPOOG");
    EvHLbl = m->newEvent(DEMOSAIC, "HISTORY_MSG_HLBL");
    EvHLth = m->newEvent(DEMOSAIC, "HISTORY_MSG_HLTH");

    CurveListener::setMulti(true);

    std::vector<GradientMilestone> bottomMilestones;
    bottomMilestones.push_back(GradientMilestone(0., 0., 0., 0.));
    bottomMilestones.push_back(GradientMilestone(1., 1., 1., 1.));

//----------- OOG clamping ----------------------------------
    clampOOG = Gtk::manage(new Gtk::CheckButton(M("TP_EXPOSURE_CLAMPOOG")));
    pack_start(*clampOOG);
    pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)));
    clampOOG->signal_toggled().connect(sigc::mem_fun(*this, &ToneCurve::clampOOGChanged));

//----------- Auto Levels ----------------------------------
    abox = Gtk::manage (new Gtk::Box ());
    abox->set_spacing (4);

    autolevels = Gtk::manage(new Gtk::ToggleButton(M("TP_EXPOSURE_AUTOLEVELS")));
    autolevels->set_tooltip_markup(M("TP_EXPOSURE_AUTOLEVELS_TOOLTIP"));
    autoconn = autolevels->signal_toggled().connect(sigc::mem_fun(*this, &ToneCurve::autolevels_toggled));

    lclip = Gtk::manage(new Gtk::Label(M("TP_EXPOSURE_CLIP")));
    lclip->set_tooltip_text(M("TP_EXPOSURE_CLIP_TOOLTIP"));

    sclip = Gtk::manage(new MySpinButton());
    sclip->set_range(0.0, 0.99);
    sclip->set_increments(0.01, 0.10);
    sclip->set_value(0.02);
    sclip->set_digits(2);
    sclip->set_width_chars(4);
    sclip->set_max_width_chars(4);
    sclip->signal_value_changed().connect(sigc::mem_fun(*this, &ToneCurve::clip_changed));

    neutral = Gtk::manage(new Gtk::Button(M("TP_NEUTRAL")));
    neutral->set_tooltip_text(M("TP_NEUTRAL_TOOLTIP"));
    neutralconn = neutral->signal_pressed().connect(sigc::mem_fun(*this, &ToneCurve::neutral_pressed));
    neutral->show();

    abox->pack_start(*autolevels, true, true, 0);
    // pack_end is used for these controls as autolevels is replaceable using pack_start in batchmode
    abox->pack_end(*neutral, true, true, 0);
    abox->pack_end(*sclip, false, false, 0);
    abox->pack_end(*lclip, false, false, 0);
    pack_start(*abox);

//-------------- Highlight Reconstruction -----------------
    pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)));

    hrenabled = Gtk::manage (new Gtk::CheckButton (M("TP_HLREC_LABEL")));
    hrenabled->set_active (false);
    hrenabled->set_tooltip_markup (M("TP_HLREC_ENA_TOOLTIP"));

    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_HLREC_LUMINANCE"));
    method->append (M("TP_HLREC_CIELAB"));
    method->append (M("TP_HLREC_BLEND"));
    method->append (M("TP_HLREC_COLOR"));
    method->append (M("TP_HLREC_COLOROPP"));
    Gtk::Box *hrVBox;
    hrVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    hrVBox->set_spacing(2);

    method->set_active(4);
    Gtk::Frame* const hrFrame = Gtk::manage(new Gtk::Frame());
    hrFrame->set_label_align(0.025, 0.5);
    hrFrame->set_label_widget(*hrenabled);

    hlrbox = Gtk::manage (new Gtk::Box ());
    Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_HLREC_METHOD")));
    hlrbox->pack_start (*lab, Gtk::PACK_SHRINK);
    hlrbox->pack_start (*method);
    hlbl = Gtk::manage(new Adjuster(M("TP_HLREC_HLBLUR"), 0, 4, 1, 0));
    hlth = Gtk::manage(new Adjuster(M("TP_HLREC_HLTH"), 0.25, 1.75, 0.01, 1.));

    hrVBox->pack_start(*hlrbox, Gtk::PACK_SHRINK);
    hrVBox->pack_start(*hlbl);
    hrVBox->pack_start(*hlth);
    hrFrame->add(*hrVBox);
    pack_start(*hrFrame);

    enaconn  = hrenabled->signal_toggled().connect( sigc::mem_fun(*this, &ToneCurve::hrenabledChanged) );
    methconn = method->signal_changed().connect ( sigc::mem_fun(*this, &ToneCurve::methodChanged) );

    //----------- Exposure Compensation ---------------------
    pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)));

    expcomp = Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), -5, 12, 0.05, 0));
    expcomp->setLogScale(2, 0, true);
    pack_start(*expcomp);

    //----------- Highlight recovery & threshold -------------
    hlcompr = Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0));
    pack_start(*hlcompr);
    hlcomprthresh = Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 0));
    pack_start(*hlcomprthresh);

//----------- Black Level & Compression -------------------
    black = Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0));
    black->setLogScale(10, 0, true);
    pack_start(*black);
    shcompr = Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50));
    pack_start(*shcompr);

    pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)));

//---------Brightness / Contrast -------------------------
    brightness = Gtk::manage(new Adjuster(M("TP_EXPOSURE_BRIGHTNESS"), -100, 100, 1, 0));
    pack_start(*brightness);
    contrast = Gtk::manage(new Adjuster(M("TP_EXPOSURE_CONTRAST"), -100, 100, 1, 0));
    pack_start(*contrast);
    saturation = Gtk::manage(new Adjuster(M("TP_EXPOSURE_SATURATION"), -100, 100, 1, 0));
    pack_start(*saturation);

    brightness->setLogScale(2, 0, true);
    contrast->setLogScale(2, 0, true);
    saturation->setLogScale(2, 0, true);

//----------- Curve 1 ------------------------------
    pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)));

    histmatching = Gtk::manage(new Gtk::ToggleButton(M("TP_EXPOSURE_HISTMATCHING")));
    histmatching->set_tooltip_markup(M("TP_EXPOSURE_HISTMATCHING_TOOLTIP"));
    histmatchconn = histmatching->signal_toggled().connect(sigc::mem_fun(*this, &ToneCurve::histmatchingToggled));
    pack_start(*histmatching, true, true, 2);

    toneCurveMode = Gtk::manage(new MyComboBoxText());
    toneCurveMode->append(M("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode->append(M("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode->append(M("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode->append(M("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode->append(M("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode->append(M("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode->set_active(0);
    toneCurveMode->set_tooltip_text(M("TP_EXPOSURE_TCMODE_LABEL1"));

    curveEditorG = new CurveEditorGroup(options.lastToneCurvesDir, M("TP_EXPOSURE_CURVEEDITOR1"));
    curveEditorG->setCurveListener(this);

    shape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "", toneCurveMode));
    shape->setEditID(EUID_ToneCurve1, BT_IMAGEFLOAT);
    shape->setBottomBarBgGradient(bottomMilestones);
    shape->setLeftBarBgGradient(bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG->curveListComplete();

    pack_start(*curveEditorG, Gtk::PACK_SHRINK, 2);

    tcmodeconn = toneCurveMode->signal_changed().connect(sigc::mem_fun(*this, &ToneCurve::curveMode1Changed), true);

//----------- Curve 2 ------------------------------

    toneCurveMode2 = Gtk::manage(new MyComboBoxText());
    toneCurveMode2->append(M("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode2->append(M("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode2->append(M("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode2->append(M("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode2->append(M("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode2->append(M("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode2->set_active(0);
    toneCurveMode2->set_tooltip_text(M("TP_EXPOSURE_TCMODE_LABEL2"));

    curveEditorG2 = new CurveEditorGroup(options.lastToneCurvesDir, M("TP_EXPOSURE_CURVEEDITOR2"));
    curveEditorG2->setCurveListener(this);

    shape2 = static_cast<DiagonalCurveEditor*>(curveEditorG2->addCurve(CT_Diagonal, "", toneCurveMode2));
    shape2->setEditID(EUID_ToneCurve2, BT_IMAGEFLOAT);
    shape2->setBottomBarBgGradient(bottomMilestones);
    shape2->setLeftBarBgGradient(bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG2->curveListComplete();
    curveEditorG2->setTooltip(M("TP_EXPOSURE_CURVEEDITOR2_TOOLTIP"));

    pack_start(*curveEditorG2, Gtk::PACK_SHRINK, 2);

    tcmode2conn = toneCurveMode2->signal_changed().connect(sigc::mem_fun(*this, &ToneCurve::curveMode2Changed), true);

// --------- Set Up Listeners -------------
    expcomp->setAdjusterListener(this);
    brightness->setAdjusterListener(this);
    black->setAdjusterListener(this);
    hlcompr->setAdjusterListener(this);
    hlbl->setAdjusterListener(this);
    hlth->setAdjusterListener(this);
    hlcomprthresh->setAdjusterListener(this);
    shcompr->setAdjusterListener(this);
    contrast->setAdjusterListener(this);
    saturation->setAdjusterListener(this);
}

ToneCurve::~ToneCurve()
{
    idle_register.destroy();

    delete curveEditorG;
    delete curveEditorG2;
}

void ToneCurve::read(const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener();

    tcmodeconn.block(true);
    tcmode2conn.block(true);
    autoconn.block(true);

    autolevels->set_active(pp->toneCurve.autoexp);
    lastAuto = pp->toneCurve.autoexp;
    sclip->set_value(pp->toneCurve.clip);

    expcomp->setValue(pp->toneCurve.expcomp);
    black->setValue(pp->toneCurve.black);
    hlcompr->setValue(pp->toneCurve.hlcompr);
    hlbl->setValue(pp->toneCurve.hlbl);
    hlth->setValue(pp->toneCurve.hlth);
    hlcomprthresh->setValue(pp->toneCurve.hlcomprthresh);
    shcompr->setValue(pp->toneCurve.shcompr);

    if (!black->getAddMode() && !batchMode) {
        shcompr->set_sensitive(!((int)black->getValue() == 0));    //at black=0 shcompr value has no effect
    }

    if (!hlcompr->getAddMode() && !batchMode) {
        hlcomprthresh->set_sensitive(!((int)hlcompr->getValue() == 0));    //at hlcompr=0 hlcomprthresh value has no effect
    }

    brightness->setValue(pp->toneCurve.brightness);
    contrast->setValue(pp->toneCurve.contrast);
    saturation->setValue(pp->toneCurve.saturation);
    shape->setCurve(pp->toneCurve.curve);
    shape2->setCurve(pp->toneCurve.curve2);

    toneCurveMode->set_active(rtengine::toUnderlying(pp->toneCurve.curveMode));
    toneCurveMode2->set_active(rtengine::toUnderlying(pp->toneCurve.curveMode2));

    histmatching->set_active(pp->toneCurve.histmatching);
    fromHistMatching = pp->toneCurve.fromHistMatching;
    clampOOG->set_active(pp->toneCurve.clampOOG);

    if (pedited) {
        expcomp->setEditedState(pedited->toneCurve.expcomp ? Edited : UnEdited);
        black->setEditedState(pedited->toneCurve.black ? Edited : UnEdited);
        hlcompr->setEditedState(pedited->toneCurve.hlcompr ? Edited : UnEdited);
        hlbl->setEditedState(pedited->toneCurve.hlbl ? Edited : UnEdited);
        hlth->setEditedState(pedited->toneCurve.hlth ? Edited : UnEdited);
        hlcomprthresh->setEditedState(pedited->toneCurve.hlcomprthresh ? Edited : UnEdited);
        shcompr->setEditedState(pedited->toneCurve.shcompr ? Edited : UnEdited);
        brightness->setEditedState(pedited->toneCurve.brightness ? Edited : UnEdited);
        contrast->setEditedState(pedited->toneCurve.contrast ? Edited : UnEdited);
        saturation->setEditedState(pedited->toneCurve.saturation ? Edited : UnEdited);
        autolevels->set_inconsistent(!pedited->toneCurve.autoexp);
        clipDirty = pedited->toneCurve.clip;
        shape->setUnChanged(!pedited->toneCurve.curve);
        shape2->setUnChanged(!pedited->toneCurve.curve2);
        hrenabled->set_inconsistent(!pedited->toneCurve.hrenabled);

        if (!pedited->toneCurve.curveMode) {
            toneCurveMode->set_active(6);
        }

        if (!pedited->toneCurve.curveMode2) {
            toneCurveMode2->set_active(6);
        }

        histmatching->set_inconsistent(!pedited->toneCurve.histmatching);
        clampOOG->set_inconsistent(!pedited->toneCurve.clampOOG);
    }

    enaconn.block(true);
    hrenabled->set_active(pp->toneCurve.hrenabled);
    enaconn.block(false);

   if (pedited && !pedited->toneCurve.method) {
        method->set_active(5);
	} else if (pp->toneCurve.method == "Luminance") {
        method->set_active(0);
    } else if (pp->toneCurve.method == "CIELab blending") {
        method->set_active(1);
    } else if (pp->toneCurve.method == "Blend") {
        method->set_active(2);
    } else if (pp->toneCurve.method == "Color") {
        method->set_active(3);
    } else if (pp->toneCurve.method == "Coloropp") {
        method->set_active(4);
    }

    hrenabledChanged();

    lasthrEnabled = pp->toneCurve.hrenabled;

    autoconn.block(false);
    tcmode2conn.block(false);
    tcmodeconn.block(false);

    enableListener();
}

void ToneCurve::autoOpenCurve()
{
    shape->openIfNonlinear();
    shape2->openIfNonlinear();
}

void ToneCurve::setEditProvider(EditDataProvider *provider)
{
    shape->setEditProvider(provider);
    shape2->setEditProvider(provider);
}

void ToneCurve::write(ProcParams* pp, ParamsEdited* pedited)
{

    pp->toneCurve.autoexp = autolevels->get_active();
    pp->toneCurve.clip = sclip->get_value();
    pp->toneCurve.expcomp = expcomp->getValue();
    pp->toneCurve.black = black->getValue();
    pp->toneCurve.hlcompr = hlcompr->getValue();
    pp->toneCurve.hlbl = hlbl->getValue();
    pp->toneCurve.hlth = hlth->getValue();
    pp->toneCurve.hlcomprthresh = hlcomprthresh->getValue();
    pp->toneCurve.shcompr = shcompr->getValue();
    pp->toneCurve.brightness = brightness->getValue();
    pp->toneCurve.contrast = contrast->getValue();
    pp->toneCurve.saturation = saturation->getValue();
    pp->toneCurve.curve = shape->getCurve();
    pp->toneCurve.curve2 = shape2->getCurve();

    int tcMode = toneCurveMode->get_active_row_number();

    if (tcMode == 0) {
        pp->toneCurve.curveMode = ToneCurveMode::STD;
    } else if (tcMode == 1) {
        pp->toneCurve.curveMode = ToneCurveMode::WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->toneCurve.curveMode = ToneCurveMode::FILMLIKE;
    } else if (tcMode == 3) {
        pp->toneCurve.curveMode = ToneCurveMode::SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->toneCurve.curveMode = ToneCurveMode::LUMINANCE;
    } else if (tcMode == 5) {
        pp->toneCurve.curveMode = ToneCurveMode::PERCEPTUAL;
    }

    tcMode = toneCurveMode2->get_active_row_number();

    if (tcMode == 0) {
        pp->toneCurve.curveMode2 = ToneCurveMode::STD;
    } else if (tcMode == 1) {
        pp->toneCurve.curveMode2 = ToneCurveMode::WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->toneCurve.curveMode2 = ToneCurveMode::FILMLIKE;
    } else if (tcMode == 3) {
        pp->toneCurve.curveMode2 = ToneCurveMode::SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->toneCurve.curveMode2 = ToneCurveMode::LUMINANCE;
    } else if (tcMode == 5) {
        pp->toneCurve.curveMode2 = ToneCurveMode::PERCEPTUAL;
    }

    pp->toneCurve.histmatching = histmatching->get_active();
    pp->toneCurve.fromHistMatching = fromHistMatching;
    pp->toneCurve.clampOOG = clampOOG->get_active();

    if (pedited) {
        pedited->toneCurve.expcomp = expcomp->getEditedState();
        pedited->toneCurve.black = black->getEditedState();
        pedited->toneCurve.hlcompr = hlcompr->getEditedState();
        pedited->toneCurve.hlbl = hlbl->getEditedState();
        pedited->toneCurve.hlth = hlth->getEditedState();
        pedited->toneCurve.hlcomprthresh = hlcomprthresh->getEditedState();
        pedited->toneCurve.shcompr = shcompr->getEditedState();
        pedited->toneCurve.brightness = brightness->getEditedState();
        pedited->toneCurve.contrast = contrast->getEditedState();
        pedited->toneCurve.saturation = saturation->getEditedState();
        pedited->toneCurve.autoexp = !autolevels->get_inconsistent();
        pedited->toneCurve.clip = clipDirty;
        pedited->toneCurve.curve = !shape->isUnChanged();
        pedited->toneCurve.curve2 = !shape2->isUnChanged();
        pedited->toneCurve.curveMode = toneCurveMode->get_active_row_number() != 6;
        pedited->toneCurve.curveMode2 = toneCurveMode2->get_active_row_number() != 6;
        pedited->toneCurve.method = method->get_active_row_number() != 5;
        pedited->toneCurve.hrenabled = !hrenabled->get_inconsistent();
        pedited->toneCurve.histmatching = !histmatching->get_inconsistent();
        pedited->toneCurve.fromHistMatching = true;
        pedited->toneCurve.clampOOG = !clampOOG->get_inconsistent();
    }

    pp->toneCurve.hrenabled = hrenabled->get_active();

    if (method->get_active_row_number() == 0) {
        pp->toneCurve.method = "Luminance";
    } else if (method->get_active_row_number() == 1) {
        pp->toneCurve.method = "CIELab blending";
    } else if (method->get_active_row_number() == 2) {
        pp->toneCurve.method = "Blend";
    } else if (method->get_active_row_number() == 3) {
        pp->toneCurve.method = "Color";
    } else if (method->get_active_row_number() == 4) {
        pp->toneCurve.method = "Coloropp";
    }
}

void ToneCurve::hrenabledChanged()
{

    if (multiImage) {
        if (hrenabled->get_inconsistent()) {
            hrenabled->set_inconsistent(false);
            enaconn.block(true);
            hrenabled->set_active(false);
            enaconn.block(false);
        } else if (lasthrEnabled) {
            hrenabled->set_inconsistent(true);
        }

        lasthrEnabled = hrenabled->get_active();
    }

    if (!batchMode) {
        if (hrenabled->get_active()) {
            hlrbox->show();
            hlrbox->set_sensitive(true);
            if (method->get_active_row_number() == 3) {
                hlbl->show();
                hlth->hide();
            } else if (method->get_active_row_number() == 4){
                hlbl->hide();
                hlth->show();
            } else {
                hlbl->hide();
				hlth->hide();	
			}
        } else {
            hlrbox->show();
            hlrbox->set_sensitive(false);
            hlbl->hide();
            hlth->hide();
        }
   }

    if (listener) {
        // Switch off auto exposure if user changes enabled manually
        if (autolevels->get_active()) {
            autoconn.block(true);
            autolevels->set_active(false);
            autoconn.block(false);
            autolevels->set_inconsistent(false);
        }

        setHistmatching(false);

        if (hrenabled->get_active()) {
            listener->panelChanged(EvHREnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvHREnabled, M("GENERAL_DISABLED"));
        }
    }
}
void ToneCurve::methodChanged()
{

    if (method->get_active_row_number() == 3) {
        hlbl->show();
        hlth->hide();	
    } else if (method->get_active_row_number() == 4){
        hlbl->hide();    
        hlth->show();    
    } else {
        hlbl->hide();
		hlth->hide();
	}
    if (listener) {
        setHistmatching(false);
        if (hrenabled->get_active()) {
            listener->panelChanged(EvHRMethod, method->get_active_text());
        }
    }
}

void ToneCurve::clampOOGChanged()
{
    if (listener) {
        listener->panelChanged(EvClampOOG, clampOOG->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void ToneCurve::setRaw(bool raw)
{

    disableListener();
    method->set_sensitive(raw);
    hlbl->set_sensitive(raw);
    hlth->set_sensitive(raw);
    hrenabled->set_sensitive(raw);
    histmatching->set_sensitive(raw);
    enableListener();
}

void ToneCurve::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited)
{

    expcomp->setDefault(defParams->toneCurve.expcomp);
    brightness->setDefault(defParams->toneCurve.brightness);
    black->setDefault(defParams->toneCurve.black);
    hlcompr->setDefault(defParams->toneCurve.hlcompr);
    hlbl->setDefault(defParams->toneCurve.hlbl);
    hlth->setDefault(defParams->toneCurve.hlth);
    hlcomprthresh->setDefault(defParams->toneCurve.hlcomprthresh);
    shcompr->setDefault(defParams->toneCurve.shcompr);
    contrast->setDefault(defParams->toneCurve.contrast);
    saturation->setDefault(defParams->toneCurve.saturation);

    if (pedited) {
        expcomp->setDefaultEditedState(pedited->toneCurve.expcomp ? Edited : UnEdited);
        black->setDefaultEditedState(pedited->toneCurve.black ? Edited : UnEdited);
        hlcompr->setDefaultEditedState(pedited->toneCurve.hlcompr ? Edited : UnEdited);
        hlbl->setDefaultEditedState(pedited->toneCurve.hlbl ? Edited : UnEdited);
        hlth->setDefaultEditedState(pedited->toneCurve.hlth ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState(pedited->toneCurve.hlcomprthresh ? Edited : UnEdited);
        shcompr->setDefaultEditedState(pedited->toneCurve.shcompr ? Edited : UnEdited);
        brightness->setDefaultEditedState(pedited->toneCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState(pedited->toneCurve.contrast ? Edited : UnEdited);
        saturation->setDefaultEditedState(pedited->toneCurve.saturation ? Edited : UnEdited);
    } else {
        expcomp->setDefaultEditedState(Irrelevant);
        black->setDefaultEditedState(Irrelevant);
        hlcompr->setDefaultEditedState(Irrelevant);
        hlbl->setDefaultEditedState(Irrelevant);
        hlth->setDefaultEditedState(Irrelevant);
        hlcomprthresh->setDefaultEditedState(Irrelevant);
        shcompr->setDefaultEditedState(Irrelevant);
        brightness->setDefaultEditedState(Irrelevant);
        contrast->setDefaultEditedState(Irrelevant);
        saturation->setDefaultEditedState(Irrelevant);
    }

}

void ToneCurve::curveChanged(CurveEditor* ce)
{

    if (listener) {
        setHistmatching(false);
        if (ce == shape) {
            listener->panelChanged(EvToneCurve1, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == shape2) {
            listener->panelChanged(EvToneCurve2, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void ToneCurve::curveMode1Changed()
{
    if (listener) {
        setHistmatching(false);
        Glib::signal_idle().connect(sigc::mem_fun(*this, &ToneCurve::curveMode1Changed_));
    }
}

bool ToneCurve::curveMode1Changed_()
{
    if (listener) {
        listener->panelChanged(EvToneCurveMode1, toneCurveMode->get_active_text());
    }

    return false;
}

void ToneCurve::curveMode2Changed()
{
    if (listener) {
        setHistmatching(false);
        Glib::signal_idle().connect(sigc::mem_fun(*this, &ToneCurve::curveMode2Changed_));
    }
}

bool ToneCurve::curveMode2Changed_()
{
    if (listener) {
        listener->panelChanged(EvToneCurveMode2, toneCurveMode2->get_active_text());
    }

    return false;
}

float ToneCurve::blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3)
{
    // assuming that all the channels are used...
    if (ce == shape) {
        if (toneCurveMode->get_active_row_number() == 4) {
            return chan1 * 0.2126729f + chan2 * 0.7151521f + chan3 * 0.0721750f;
        }
    } else if (ce == shape2) {
        if (toneCurveMode2->get_active_row_number() == 4) {
            return chan1 * 0.2126729f + chan2 * 0.7151521f + chan3 * 0.0721750f;
        }
    }

    return CurveListener::blendPipetteValues(ce, chan1, chan2, chan3);
}

void ToneCurve::adjusterChanged(Adjuster* a, double newval)
{
    // Switch off auto exposure if user changes sliders manually
    if (autolevels->get_active() && (a == expcomp || a == brightness || a == contrast || a == black || a == hlcompr || a == hlcomprthresh)) {
        autoconn.block(true);
        autolevels->set_active(false);
        autoconn.block(false);
        autolevels->set_inconsistent(false);
    }

    if (!listener) {
        return;
    }

    if (a != expcomp && a != hlcompr && a != hlcomprthresh) {
        setHistmatching(false);
    }

    Glib::ustring costr;

    if (a == expcomp) {
        costr = Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), a->getValue());
    } else {
        costr = Glib::ustring::format((int)a->getValue());
    }

    if (a == expcomp) {
        listener->panelChanged(EvExpComp, costr);
    } else if (a == brightness) {
        listener->panelChanged(EvBrightness, costr);
    } else if (a == black) {
        listener->panelChanged(EvBlack, costr);

        if (!black->getAddMode() && !batchMode) {
            shcompr->set_sensitive(!((int)black->getValue() == 0));    //at black=0 shcompr value has no effect
        }
    } else if (a == contrast) {
        listener->panelChanged(EvContrast, costr);
    } else if (a == saturation) {
        listener->panelChanged(EvSaturation, costr);
    } else if (a == hlbl) {
        listener->panelChanged(EvHLbl, costr);
    } else if (a == hlth) {
        listener->panelChanged(EvHLth, costr);
    } else if (a == hlcompr) {
        listener->panelChanged(EvHLCompr, costr);

        if (!hlcompr->getAddMode() && !batchMode) {
            hlcomprthresh->set_sensitive(!((int)hlcompr->getValue() == 0));    //at hlcompr=0 hlcomprthresh value has no effect
        }
    } else if (a == hlcomprthresh) {
        listener->panelChanged(EvHLComprThreshold, costr);
    } else if (a == shcompr) {
        listener->panelChanged(EvSHCompr, costr);
    }
}

void ToneCurve::neutral_pressed()
{
// This method deselects auto levels and HL reconstruction auto
// and sets neutral values to params in exposure panel

    setHistmatching(false);

    if (batchMode) {
        autolevels->set_inconsistent(false);
        autoconn.block(true);
        autolevels->set_active(false);
        autoconn.block(false);

        lastAuto = autolevels->get_active();
    } else { //!batchMode
        autolevels->set_active(false);
        autolevels->set_inconsistent(false);
    }

    expcomp->setValue(0);
    hlcompr->setValue(0);
    hlbl->setValue(0);
    hlth->setValue(1.0);
    hlcomprthresh->setValue(0);
    brightness->setValue(0);
    black->setValue(0);
    shcompr->setValue(50);
    enaconn.block(true);
    hrenabled->set_active(false);
    enaconn.block(false);

    if (!batchMode) {
        hlrbox->show();
        hlrbox->set_sensitive(false);
        hlbl->hide();
        hlth->hide();
    }

    if (!black->getAddMode() && !batchMode) {
        shcompr->set_sensitive(!((int)black->getValue() == 0));    //at black=0 shcompr value has no effect
    }

    if (!hlcompr->getAddMode() && !batchMode) {
        hlcomprthresh->set_sensitive(!((int)hlcompr->getValue() == 0));    //at hlcompr=0 hlcomprthresh value has no effect
    }

    contrast->setValue(0);
    //saturation->setValue(0);

    listener->panelChanged(EvNeutralExp, M("GENERAL_ENABLED"));
}

void ToneCurve::autolevels_toggled()
{
    setHistmatching(false);
    if (batchMode) {
        if (autolevels->get_inconsistent()) {
            autolevels->set_inconsistent(false);
            autoconn.block(true);
            autolevels->set_active(false);
            autoconn.block(false);
        } else if (lastAuto) {
            autolevels->set_inconsistent(true);
        }

        lastAuto = autolevels->get_active();

        expcomp->setEditedState(UnEdited);
        brightness->setEditedState(UnEdited);
        contrast->setEditedState(UnEdited);
        black->setEditedState(UnEdited);
        hlcompr->setEditedState(UnEdited);
        hlcomprthresh->setEditedState(UnEdited);

        if (expcomp->getAddMode()) {
            expcomp->setValue(0);
        }

        if (brightness->getAddMode()) {
            brightness->setValue(0);
        }

        if (contrast->getAddMode()) {
            contrast->setValue(0);
        }

        if (black->getAddMode()) {
            black->setValue(0);
        }

        if (hlcompr->getAddMode()) {
            hlcompr->setValue(0);
        }

        if (hlcomprthresh->getAddMode()) {
            hlcomprthresh->setValue(0);
        }

        if (listener) {
            if (!autolevels->get_inconsistent()) {
                if (autolevels->get_active()) {
                    listener->panelChanged(EvAutoExp, M("GENERAL_ENABLED"));
                } else {
                    listener->panelChanged(EvFixedExp, M("GENERAL_DISABLED"));
                }
            }
        }
    } else if (/* !batchMode && */ listener) {
        if (autolevels->get_active()) {
            listener->panelChanged(EvAutoExp, M("GENERAL_ENABLED"));
            waitForAutoExp();

            if (!black->getAddMode()) {
                shcompr->set_sensitive(!((int)black->getValue() == 0));    //at black=0 shcompr value has no effect
            }

            if (!hlcompr->getAddMode() && !batchMode) {
                hlcomprthresh->set_sensitive(!((int)hlcompr->getValue() == 0));    //at hlcompr=0 hlcomprthresh value has no effect
            }

        } else {
            listener->panelChanged(EvFixedExp, M("GENERAL_DISABLED"));
        }
    }
}

void ToneCurve::clip_changed()
{

    clipDirty = true;

    if (autolevels->get_active() && listener) {
        Glib::signal_idle().connect(sigc::mem_fun(*this, &ToneCurve::clip_changed_));
    }
}

bool ToneCurve::clip_changed_()
{

    if (listener) {
        listener->panelChanged(EvClip, Glib::ustring::format(std::setprecision(5), sclip->get_value()));

        if (!batchMode) {
            waitForAutoExp();
        }
    }

    return false;
}

void ToneCurve::waitForAutoExp()
{

    sclip->set_sensitive(false);
    expcomp->setEnabled(false);
    brightness->setEnabled(false);
    contrast->setEnabled(false);
    black->setEnabled(false);
    hlcompr->setEnabled(false);
    hlcomprthresh->setEnabled(false);
    shcompr->setEnabled(false);
    contrast->setEnabled(false);
    saturation->setEnabled(false);
    curveEditorG->set_sensitive(false);
    toneCurveMode->set_sensitive(false);
    curveEditorG2->set_sensitive(false);
    toneCurveMode2->set_sensitive(false);
    hrenabled->set_sensitive(false);
    method->set_sensitive(false);
    hlbl->set_sensitive(false);
    hlth->set_sensitive(false);
    histmatching->set_sensitive(false);
}

void ToneCurve::enableAll()
{

    sclip->set_sensitive(true);
    expcomp->setEnabled(true);
    brightness->setEnabled(true);
    black->setEnabled(true);
    hlcompr->setEnabled(true);
    hlbl->setEnabled(true);
    hlth->setEnabled(true);
    hlcomprthresh->setEnabled(true);
    shcompr->setEnabled(true);
    contrast->setEnabled(true);
    saturation->setEnabled(true);
    curveEditorG->set_sensitive(true);
    toneCurveMode->set_sensitive(true);
    curveEditorG2->set_sensitive(true);
    toneCurveMode2->set_sensitive(true);
    hrenabled->set_sensitive(true);
    method->set_sensitive(true);
    hlbl->set_sensitive(true);
    hlth->set_sensitive(true);
    histmatching->set_sensitive(true);
}

void ToneCurve::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    method->append(M("GENERAL_UNCHANGED"));

    removeIfThere(abox, autolevels, false);
    autolevels = Gtk::manage(new Gtk::CheckButton(M("TP_EXPOSURE_AUTOLEVELS")));
    autolevels->set_tooltip_markup(M("TP_EXPOSURE_AUTOLEVELS_TOOLTIP"));
    autoconn = autolevels->signal_toggled().connect(sigc::mem_fun(*this, &ToneCurve::autolevels_toggled));
    abox->pack_start(*autolevels);

    ToolPanel::setBatchMode(batchMode);
    expcomp->showEditedCB();
    black->showEditedCB();
    hlcompr->showEditedCB();
    hlbl->showEditedCB();
    hlth->showEditedCB();
    hlcomprthresh->showEditedCB();
    shcompr->showEditedCB();
    brightness->showEditedCB();
    contrast->showEditedCB();
    saturation->showEditedCB();

    toneCurveMode->append(M("GENERAL_UNCHANGED"));
    toneCurveMode2->append(M("GENERAL_UNCHANGED"));

    curveEditorG->setBatchMode(batchMode);
    curveEditorG2->setBatchMode(batchMode);
}

void ToneCurve::setAdjusterBehavior(bool expadd, bool hlcompadd, bool hlcompthreshadd, bool bradd, bool blackadd, bool shcompadd, bool contradd, bool satadd)
{

    expcomp->setAddMode(expadd);
    hlcompr->setAddMode(hlcompadd);
    hlcomprthresh->setAddMode(hlcompthreshadd);
    brightness->setAddMode(bradd);
    black->setAddMode(blackadd);
    shcompr->setAddMode(shcompadd);
    contrast->setAddMode(contradd);
    saturation->setAddMode(satadd);
}

void ToneCurve::trimValues(rtengine::procparams::ProcParams* pp)
{

    expcomp->trimValue(pp->toneCurve.expcomp);
    hlcompr->trimValue(pp->toneCurve.hlcompr);
    hlcomprthresh->trimValue(pp->toneCurve.hlcomprthresh);
    brightness->trimValue(pp->toneCurve.brightness);
    black->trimValue(pp->toneCurve.black);
    shcompr->trimValue(pp->toneCurve.shcompr);
    contrast->trimValue(pp->toneCurve.contrast);
    saturation->trimValue(pp->toneCurve.saturation);
}

void ToneCurve::updateCurveBackgroundHistogram(
    const LUTu& histToneCurve,
    const LUTu& histLCurve,
    const LUTu& histCCurve,
    const LUTu& histLCAM,
    const LUTu& histCCAM,
    const LUTu& histRed,
    const LUTu& histGreen,
    const LUTu& histBlue,
    const LUTu& histLuma,
    const LUTu& histLRETI
)
{
    shape->updateBackgroundHistogram(histToneCurve);
}

void ToneCurve::setHistmatching(bool enabled)
{
    fromHistMatching = enabled;
    if (histmatching->get_active()) {
        histmatchconn.block(true);
        histmatching->set_active(enabled);
        histmatchconn.block(false);
        histmatching->set_inconsistent(false);
    }
}

void ToneCurve::histmatchingToggled()
{
    if (listener) {
        if (!batchMode) {
            if (histmatching->get_active()) {
                fromHistMatching = false;
                listener->panelChanged(EvHistMatching, M("GENERAL_ENABLED"));
                waitForAutoExp();
            } else {
                listener->panelChanged(EvHistMatching, M("GENERAL_DISABLED"));
            }
        } else {
            listener->panelChanged(EvHistMatchingBatch, histmatching->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
        }
    }
}

void ToneCurve::autoExpChanged(double expcomp, int bright, int contr, int black, int hlcompr, int hlcomprthresh, bool hlrecons)
{
    nextBlack = black;
    nextExpcomp = expcomp;
    nextBrightness = bright;
    nextContrast = contr;
    nextHlcompr = hlcompr;
    nextHlcomprthresh = hlcomprthresh;
    nextHLRecons = hlrecons;

    idle_register.add(
        [this]() -> bool
        {
            GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected
            // FIXME: We don't need the GThreadLock, don't we?
            disableListener();
            enableAll();
            this->expcomp->setValue(nextExpcomp);
            brightness->setValue(nextBrightness);
            contrast->setValue(nextContrast);
            this->black->setValue(nextBlack);
            this->hlcompr->setValue(nextHlcompr);
            this->hlcomprthresh->setValue(nextHlcomprthresh);
            enaconn.block(true);
            hrenabled->set_active(nextHLRecons);
            enaconn.block(false);

            if (nextHLRecons) {
                hlrbox->show();
                hlrbox->set_sensitive(true);
                if (method->get_active_row_number() == 3) {
                    hlbl->show();
                    hlth->hide();
                } else if (method->get_active_row_number() == 4){
                    hlbl->hide();
                    hlth->show();
                } else {
					hlbl->hide();
					hlth->hide();
				}
            } else if (!batchMode) {
                hlrbox->show();
                hlrbox->set_sensitive(false);
                hlbl->hide();
                hlth->hide();
           }

            if (!this->black->getAddMode() && !batchMode) {
                shcompr->set_sensitive(static_cast<int>(this->black->getValue()));    //at black=0 shcompr value has no effect
            }

            if (!this->hlcompr->getAddMode() && !batchMode) {
                this->hlcomprthresh->set_sensitive(static_cast<int>(this->hlcompr->getValue()));    //at hlcompr=0 hlcomprthresh value has no effect
            }

            enableListener();
            return false;
        }
    );
}

void ToneCurve::autoMatchedToneCurveChanged(rtengine::procparams::ToneCurveMode curveMode, const std::vector<double>& curve)
{
    nextToneCurveMode = curveMode;
    nextToneCurve = curve;

    idle_register.add(
        [this]() -> bool
        {
            GThreadLock lock; // FIXME: Obsolete
            disableListener();
            enableAll();
            brightness->setValue(0);
            contrast->setValue(0);
            black->setValue(0);

            if (!black->getAddMode() && !batchMode) {
                shcompr->set_sensitive(static_cast<int>(black->getValue()));
            }

            if (!hlcompr->getAddMode() && !batchMode) {
                hlcomprthresh->set_sensitive(static_cast<int>(hlcompr->getValue()));    //at hlcompr=0 hlcomprthresh value has no effect
            }

            if (autolevels->get_active()) {
                expcomp->setValue(0);
                autoconn.block(true);
                autolevels->set_active(false);
                autoconn.block(false);
                autolevels->set_inconsistent(false);
            }

            toneCurveMode->set_active(rtengine::toUnderlying(nextToneCurveMode));
            shape->setCurve(nextToneCurve);
            shape2->setCurve({DCT_Linear});
            shape->openIfNonlinear();

            enableListener();
            fromHistMatching = true;

            return false;
        }
    );
}
