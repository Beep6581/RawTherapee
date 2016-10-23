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
#include "tonecurve.h"
#include "adjuster.h"
#include <sigc++/slot.h>
#include <iomanip>
#include "ppversion.h"
#include "edit.h"

using namespace rtengine;
using namespace rtengine::procparams;

ToneCurve::ToneCurve () : FoldableToolPanel(this, "tonecurve", M("TP_EXPOSURE_LABEL"))
{

    CurveListener::setMulti(true);

    std::vector<GradientMilestone> bottomMilestones;
    bottomMilestones.push_back( GradientMilestone(0., 0., 0., 0.) );
    bottomMilestones.push_back( GradientMilestone(1., 1., 1., 1.) );

//----------- Auto Levels ----------------------------------
    abox = Gtk::manage (new Gtk::HBox ());
    abox->set_spacing (10);

    autolevels = Gtk::manage (new Gtk::ToggleButton (M("TP_EXPOSURE_AUTOLEVELS")));
    autolevels->set_tooltip_markup (M("TP_EXPOSURE_AUTOLEVELS_TIP"));
    autoconn = autolevels->signal_toggled().connect( sigc::mem_fun(*this, &ToneCurve::autolevels_toggled) );

    lclip = Gtk::manage (new Gtk::Label (M("TP_EXPOSURE_CLIP")));
    lclip->set_tooltip_text (M("TP_EXPOSURE_CLIP_TIP"));

    sclip = Gtk::manage (new MySpinButton ());
    sclip->set_range (0.0, 0.99);
    sclip->set_increments (0.01, 0.10);
    sclip->set_value (0.02);
    sclip->set_digits (2);
    sclip->set_width_chars(4);
    sclip->set_max_width_chars(4);
    sclip->signal_value_changed().connect( sigc::mem_fun(*this, &ToneCurve::clip_changed) );

    neutral = Gtk::manage (new Gtk::Button (M("TP_NEUTRAL")));
    neutral->set_tooltip_text (M("TP_NEUTRAL_TIP"));
    neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &ToneCurve::neutral_pressed) );
    neutral->show();

    abox->pack_start (*autolevels, true, true, 0);
    // pack_end is used for these controls as autolevels is replaceable using pack_start in batchmode
    abox->pack_end (*neutral, true, true, 0);
    abox->pack_end (*sclip, false, false, 0);
    abox->pack_end (*lclip, false, false, 0);
    pack_start (*abox);

//-------------- Highlight Reconstruction -----------------
    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    hrenabled = Gtk::manage (new Gtk::CheckButton (M("TP_HLREC_LABEL")));
    hrenabled->set_active (false);
    hrenabled->set_tooltip_markup (M("TP_HLREC_ENA_TOOLTIP"));
    pack_start (*hrenabled);

    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_HLREC_LUMINANCE"));
    method->append (M("TP_HLREC_CIELAB"));
    method->append (M("TP_HLREC_COLOR"));
    method->append (M("TP_HLREC_BLEND"));

    method->set_active (0);
    hlrbox = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_HLREC_METHOD")));
    hlrbox->pack_start (*lab, Gtk::PACK_SHRINK, 4);
    hlrbox->pack_start (*method);
    pack_start (*hlrbox);

    enaconn  = hrenabled->signal_toggled().connect( sigc::mem_fun(*this, &ToneCurve::hrenabledChanged) );
    methconn = method->signal_changed().connect ( sigc::mem_fun(*this, &ToneCurve::methodChanged) );

    //----------- Exposure Compensation ---------------------
    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    expcomp   = Gtk::manage (new Adjuster (M("TP_EXPOSURE_EXPCOMP"), -5, 12, 0.05, 0));
    pack_start (*expcomp);

    //----------- Highlight recovery & threshold -------------
    hlcompr = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0));
    pack_start (*hlcompr);
    hlcomprthresh = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 33));
    pack_start (*hlcomprthresh);

//----------- Black Level & Compression -------------------
    black = Gtk::manage (new Adjuster (M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0));
    pack_start (*black);
    shcompr = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50));
    pack_start (*shcompr);

    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

//---------Brightness / Contrast -------------------------
    brightness = Gtk::manage (new Adjuster (M("TP_EXPOSURE_BRIGHTNESS"), -100, 100, 1, 0));
    pack_start (*brightness);
    contrast   = Gtk::manage (new Adjuster (M("TP_EXPOSURE_CONTRAST"), -100, 100, 1, 0));
    pack_start (*contrast);
    saturation = Gtk::manage (new Adjuster (M("TP_EXPOSURE_SATURATION"), -100, 100, 1, 0));
    pack_start (*saturation);

//----------- Curve 1 ------------------------------
    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    toneCurveMode = Gtk::manage (new MyComboBoxText ());
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode->set_active (0);
    toneCurveMode->set_tooltip_text(M("TP_EXPOSURE_TCMODE_LABEL1"));

    curveEditorG = new CurveEditorGroup (options.lastToneCurvesDir, M("TP_EXPOSURE_CURVEEDITOR1"));
    curveEditorG->setCurveListener (this);

    shape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "", toneCurveMode));
    shape->setEditID(EUID_ToneCurve1, BT_IMAGEFLOAT);
    shape->setBottomBarBgGradient(bottomMilestones);
    shape->setLeftBarBgGradient(bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG->curveListComplete();

    pack_start( *curveEditorG, Gtk::PACK_SHRINK, 2);

    tcmodeconn = toneCurveMode->signal_changed().connect( sigc::mem_fun(*this, &ToneCurve::curveMode1Changed), true );

//----------- Curve 2 ------------------------------

    toneCurveMode2 = Gtk::manage (new MyComboBoxText ());
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode2->set_active (0);
    toneCurveMode2->set_tooltip_text(M("TP_EXPOSURE_TCMODE_LABEL2"));

    curveEditorG2 = new CurveEditorGroup (options.lastToneCurvesDir, M("TP_EXPOSURE_CURVEEDITOR2"));
    curveEditorG2->setCurveListener (this);

    shape2 = static_cast<DiagonalCurveEditor*>(curveEditorG2->addCurve(CT_Diagonal, "", toneCurveMode2));
    shape2->setEditID(EUID_ToneCurve2, BT_IMAGEFLOAT);
    shape2->setBottomBarBgGradient(bottomMilestones);
    shape2->setLeftBarBgGradient(bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG2->curveListComplete();
    curveEditorG2->setTooltip(M("TP_EXPOSURE_CURVEEDITOR2_TOOLTIP"));

    pack_start( *curveEditorG2, Gtk::PACK_SHRINK, 2);

    tcmode2conn = toneCurveMode2->signal_changed().connect( sigc::mem_fun(*this, &ToneCurve::curveMode2Changed), true );

// --------- Set Up Listeners -------------
    expcomp->setAdjusterListener (this);
    brightness->setAdjusterListener (this);
    black->setAdjusterListener (this);
    hlcompr->setAdjusterListener (this);
    hlcomprthresh->setAdjusterListener (this);
    shcompr->setAdjusterListener (this);
    contrast->setAdjusterListener (this);
    saturation->setAdjusterListener (this);
}

ToneCurve::~ToneCurve ()
{
    delete curveEditorG;
    delete curveEditorG2;
}

void ToneCurve::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    tcmodeconn.block(true);
    tcmode2conn.block(true);
    autoconn.block (true);

    autolevels->set_active (pp->toneCurve.autoexp);
    lastAuto = pp->toneCurve.autoexp;
    sclip->set_value (pp->toneCurve.clip);

    expcomp->setValue (pp->toneCurve.expcomp);
    black->setValue (pp->toneCurve.black);
    hlcompr->setValue (pp->toneCurve.hlcompr);
    hlcomprthresh->setValue (pp->toneCurve.hlcomprthresh);
    shcompr->setValue (pp->toneCurve.shcompr);

    if (!black->getAddMode()) {
        shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
    }

    brightness->setValue (pp->toneCurve.brightness);
    contrast->setValue (pp->toneCurve.contrast);
    saturation->setValue (pp->toneCurve.saturation);
    shape->setCurve (pp->toneCurve.curve);
    shape2->setCurve (pp->toneCurve.curve2);

    toneCurveMode->set_active(pp->toneCurve.curveMode);
    toneCurveMode2->set_active(pp->toneCurve.curveMode2);

    if (pedited) {
        expcomp->setEditedState (pedited->toneCurve.expcomp ? Edited : UnEdited);
        black->setEditedState (pedited->toneCurve.black ? Edited : UnEdited);
        hlcompr->setEditedState (pedited->toneCurve.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setEditedState (pedited->toneCurve.hlcomprthresh ? Edited : UnEdited);
        shcompr->setEditedState (pedited->toneCurve.shcompr ? Edited : UnEdited);
        brightness->setEditedState (pedited->toneCurve.brightness ? Edited : UnEdited);
        contrast->setEditedState (pedited->toneCurve.contrast ? Edited : UnEdited);
        saturation->setEditedState (pedited->toneCurve.saturation ? Edited : UnEdited);
        autolevels->set_inconsistent (!pedited->toneCurve.autoexp);
        clipDirty = pedited->toneCurve.clip;
        shape->setUnChanged (!pedited->toneCurve.curve);
        shape2->setUnChanged (!pedited->toneCurve.curve2);
        hrenabled->set_inconsistent (!pedited->toneCurve.hrenabled);

        if (!pedited->toneCurve.curveMode) {
            toneCurveMode->set_active(6);
        }

        if (!pedited->toneCurve.curveMode2) {
            toneCurveMode2->set_active(6);
        }
    }

    enaconn.block (true);
    hrenabled->set_active  (pp->toneCurve.hrenabled);
    enaconn.block (false);

    if (pedited && !pedited->toneCurve.method) {
        method->set_active (4);
    } else if (pp->toneCurve.method == "Luminance") {
        method->set_active (0);
    } else if (pp->toneCurve.method == "CIELab blending") {
        method->set_active (1);
    } else if (pp->toneCurve.method == "Color") {
        method->set_active (2);
    } else if (pp->toneCurve.method == "Blend") {
        method->set_active (3);
    }

    if (!batchMode) {
        if (hrenabled->get_active()) {
            hlrbox->show();
        } else {
            hlrbox->hide();
        }
    }

    lasthrEnabled = pp->toneCurve.hrenabled;

    autoconn.block (false);
    tcmode2conn.block(false);
    tcmodeconn.block(false);

    enableListener ();
}

void ToneCurve::autoOpenCurve  ()
{
    shape->openIfNonlinear();
    shape2->openIfNonlinear();
}

void ToneCurve::setEditProvider  (EditDataProvider *provider)
{
    shape->setEditProvider(provider);
    shape2->setEditProvider(provider);
}

void ToneCurve::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->toneCurve.autoexp = autolevels->get_active();
    pp->toneCurve.clip = sclip->get_value ();
    pp->toneCurve.expcomp = expcomp->getValue ();
    pp->toneCurve.black = (int)black->getValue ();
    pp->toneCurve.hlcompr = (int)hlcompr->getValue ();
    pp->toneCurve.hlcomprthresh = (int)hlcomprthresh->getValue ();
    pp->toneCurve.shcompr = (int)shcompr->getValue ();
    pp->toneCurve.brightness = (int)brightness->getValue ();
    pp->toneCurve.contrast = (int)contrast->getValue ();
    pp->toneCurve.saturation = (int)saturation->getValue ();
    pp->toneCurve.curve = shape->getCurve ();
    pp->toneCurve.curve2 = shape2->getCurve ();

    int tcMode = toneCurveMode->get_active_row_number();

    if      (tcMode == 0) {
        pp->toneCurve.curveMode = ToneCurveParams::TC_MODE_STD;
    } else if (tcMode == 1) {
        pp->toneCurve.curveMode = ToneCurveParams::TC_MODE_WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->toneCurve.curveMode = ToneCurveParams::TC_MODE_FILMLIKE;
    } else if (tcMode == 3) {
        pp->toneCurve.curveMode = ToneCurveParams::TC_MODE_SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->toneCurve.curveMode = ToneCurveParams::TC_MODE_LUMINANCE;
    } else if (tcMode == 5) {
        pp->toneCurve.curveMode = ToneCurveParams::TC_MODE_PERCEPTUAL;
    }

    tcMode = toneCurveMode2->get_active_row_number();

    if      (tcMode == 0) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_STD;
    } else if (tcMode == 1) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_FILMLIKE;
    } else if (tcMode == 3) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_LUMINANCE;
    } else if (tcMode == 5) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_PERCEPTUAL;
    }

    if (pedited) {
        pedited->toneCurve.expcomp    = expcomp->getEditedState ();
        pedited->toneCurve.black      = black->getEditedState ();
        pedited->toneCurve.hlcompr    = hlcompr->getEditedState ();
        pedited->toneCurve.hlcomprthresh = hlcomprthresh->getEditedState ();
        pedited->toneCurve.shcompr    = shcompr->getEditedState ();
        pedited->toneCurve.brightness = brightness->getEditedState ();
        pedited->toneCurve.contrast   = contrast->getEditedState ();
        pedited->toneCurve.saturation = saturation->getEditedState ();
        pedited->toneCurve.autoexp    = !autolevels->get_inconsistent();
        pedited->toneCurve.clip       = clipDirty;
        pedited->toneCurve.curve      = !shape->isUnChanged ();
        pedited->toneCurve.curve2     = !shape2->isUnChanged ();
        pedited->toneCurve.curveMode  = toneCurveMode->get_active_row_number() != 6;
        pedited->toneCurve.curveMode2 = toneCurveMode2->get_active_row_number() != 6;
        pedited->toneCurve.method     = method->get_active_row_number() != 4;
        pedited->toneCurve.hrenabled  = !hrenabled->get_inconsistent();
    }

    pp->toneCurve.hrenabled = hrenabled->get_active();

    if (method->get_active_row_number() == 0) {
        pp->toneCurve.method = "Luminance";
    } else if (method->get_active_row_number() == 1) {
        pp->toneCurve.method = "CIELab blending";
    } else if (method->get_active_row_number() == 2) {
        pp->toneCurve.method = "Color";
    } else if (method->get_active_row_number() == 3) {
        pp->toneCurve.method = "Blend";
    }
}

void ToneCurve::hrenabledChanged ()
{

    if (multiImage) {
        if (hrenabled->get_inconsistent()) {
            hrenabled->set_inconsistent (false);
            enaconn.block (true);
            hrenabled->set_active (false);
            enaconn.block (false);
        } else if (lasthrEnabled) {
            hrenabled->set_inconsistent (true);
        }

        lasthrEnabled = hrenabled->get_active ();
    }

    if (!batchMode) {
        if (hrenabled->get_active()) {
            hlrbox->show();
        } else {
            hlrbox->hide();
        }
    }

    if (listener) {
        // Switch off auto exposure if user changes enabled manually
        if (autolevels->get_active() ) {
            autoconn.block(true);
            autolevels->set_active (false);
            autoconn.block(false);
            autolevels->set_inconsistent (false);
        }

        if (hrenabled->get_active ()) {
            listener->panelChanged (EvHREnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvHREnabled, M("GENERAL_DISABLED"));
        }
    }
}
void ToneCurve::methodChanged ()
{

    if (listener) {
        if (hrenabled->get_active ()) {
            listener->panelChanged (EvHRMethod, method->get_active_text ());
        }
    }
}
void ToneCurve::setRaw (bool raw)
{

    disableListener ();
    method->set_sensitive (raw);
    hrenabled->set_sensitive (raw);
    enableListener ();
}


void ToneCurve::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    expcomp->setDefault (defParams->toneCurve.expcomp);
    brightness->setDefault (defParams->toneCurve.brightness);
    black->setDefault (defParams->toneCurve.black);
    hlcompr->setDefault (defParams->toneCurve.hlcompr);
    hlcomprthresh->setDefault (defParams->toneCurve.hlcomprthresh);
    shcompr->setDefault (defParams->toneCurve.shcompr);
    contrast->setDefault (defParams->toneCurve.contrast);
    saturation->setDefault (defParams->toneCurve.saturation);

    if (pedited) {
        expcomp->setDefaultEditedState (pedited->toneCurve.expcomp ? Edited : UnEdited);
        black->setDefaultEditedState (pedited->toneCurve.black ? Edited : UnEdited);
        hlcompr->setDefaultEditedState (pedited->toneCurve.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState (pedited->toneCurve.hlcomprthresh ? Edited : UnEdited);
        shcompr->setDefaultEditedState (pedited->toneCurve.shcompr ? Edited : UnEdited);
        brightness->setDefaultEditedState (pedited->toneCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->toneCurve.contrast ? Edited : UnEdited);
        saturation->setDefaultEditedState (pedited->toneCurve.saturation ? Edited : UnEdited);
    } else {
        expcomp->setDefaultEditedState (Irrelevant);
        black->setDefaultEditedState (Irrelevant);
        hlcompr->setDefaultEditedState (Irrelevant);
        hlcomprthresh->setDefaultEditedState (Irrelevant);
        shcompr->setDefaultEditedState (Irrelevant);
        brightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
        saturation->setDefaultEditedState (Irrelevant);
    }
}

void ToneCurve::curveChanged (CurveEditor* ce)
{

    if (listener) {
        if (ce == shape) {
            listener->panelChanged (EvToneCurve1, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == shape2) {
            listener->panelChanged (EvToneCurve2, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void ToneCurve::curveMode1Changed ()
{
    //if (listener)  listener->panelChanged (EvToneCurveMode, toneCurveMode->get_active_text());
    if (listener) {
        Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::curveMode1Changed_));
    }
}

bool ToneCurve::curveMode1Changed_ ()
{
    if (listener) {
        listener->panelChanged (EvToneCurveMode1, toneCurveMode->get_active_text());
    }

    return false;
}

void ToneCurve::curveMode2Changed ()
{
    //if (listener)  listener->panelChanged (EvToneCurveMode, toneCurveMode->get_active_text());
    if (listener) {
        Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::curveMode2Changed_));
    }
}

bool ToneCurve::curveMode2Changed_ ()
{
    if (listener) {
        listener->panelChanged (EvToneCurveMode2, toneCurveMode2->get_active_text());
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

void ToneCurve::adjusterChanged (Adjuster* a, double newval)
{

    // Switch off auto exposure if user changes sliders manually
    if (autolevels->get_active() && (a == expcomp || a == brightness || a == contrast || a == black || a == hlcompr || a == hlcomprthresh)) {
        autoconn.block(true);
        autolevels->set_active (false);
        autoconn.block(false);
        autolevels->set_inconsistent (false);
    }

    if (!listener) {
        return;
    }

    Glib::ustring costr;

    if (a == expcomp) {
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
    } else {
        costr = Glib::ustring::format ((int)a->getValue());
    }

    if (a == expcomp) {
        listener->panelChanged (EvExpComp, costr);
    } else if (a == brightness) {
        listener->panelChanged (EvBrightness, costr);
    } else if (a == black) {
        listener->panelChanged (EvBlack, costr);

        if (!black->getAddMode()) {
            shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
        }
    } else if (a == contrast) {
        listener->panelChanged (EvContrast, costr);
    } else if (a == saturation) {
        listener->panelChanged (EvSaturation, costr);
    } else if (a == hlcompr) {
        listener->panelChanged (EvHLCompr, costr);
    } else if (a == hlcomprthresh) {
        listener->panelChanged (EvHLComprThreshold, costr);
    } else if (a == shcompr) {
        listener->panelChanged (EvSHCompr, costr);
    }
}

void ToneCurve::neutral_pressed ()
{
// This method deselects auto levels and HL reconstruction auto
// and sets neutral values to params in exposure panel

    if (batchMode) {
        autolevels->set_inconsistent (false);
        autoconn.block (true);
        autolevels->set_active (false);
        autoconn.block (false);

        lastAuto = autolevels->get_active ();
    } else { //!batchMode
        autolevels->set_active (false);
        autolevels->set_inconsistent (false);
    }

    expcomp->setValue(0);
    hlcompr->setValue(0);
    hlcomprthresh->setValue(0);
    brightness->setValue(0);
    black->setValue(0);
    shcompr->setValue(50);
    enaconn.block (true);
    hrenabled->set_active (false);
    enaconn.block (false);

    if (!batchMode) {
        hlrbox->hide();
    }

    if (!black->getAddMode()) {
        shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
    }

    contrast->setValue(0);
    //saturation->setValue(0);

    listener->panelChanged (EvNeutralExp, M("GENERAL_ENABLED"));
}
void ToneCurve::autolevels_toggled ()
{

    if (batchMode) {
        if (autolevels->get_inconsistent()) {
            autolevels->set_inconsistent (false);
            autoconn.block (true);
            autolevels->set_active (false);
            autoconn.block (false);
        } else if (lastAuto) {
            autolevels->set_inconsistent (true);
        }

        lastAuto = autolevels->get_active ();

        expcomp->setEditedState (UnEdited);
        brightness->setEditedState (UnEdited);
        contrast->setEditedState (UnEdited);
        black->setEditedState (UnEdited);
        hlcompr->setEditedState (UnEdited);
        hlcomprthresh->setEditedState (UnEdited);

        if (expcomp->getAddMode()) {
            expcomp->setValue (0);
        }

        if (brightness->getAddMode()) {
            brightness->setValue (0);
        }

        if (contrast->getAddMode()) {
            contrast->setValue (0);
        }

        if (black->getAddMode()) {
            black->setValue (0);
        }

        if (hlcompr->getAddMode()) {
            hlcompr->setValue (0);
        }

        if (hlcomprthresh->getAddMode()) {
            hlcomprthresh->setValue (0);
        }

        if (listener) {
            if (!autolevels->get_inconsistent()) {
                if (autolevels->get_active ()) {
                    listener->panelChanged (EvAutoExp, M("GENERAL_ENABLED"));
                } else {
                    listener->panelChanged (EvFixedExp, M("GENERAL_DISABLED"));
                }
            }
        }
    } else if (/* !batchMode && */ listener) {
        if (autolevels->get_active()) {
            listener->panelChanged (EvAutoExp, M("GENERAL_ENABLED"));
            waitForAutoExp ();

            if (!black->getAddMode()) {
                shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
            }
        } else {
            listener->panelChanged (EvFixedExp, M("GENERAL_DISABLED"));
        }
    }
}

void ToneCurve::clip_changed ()
{

    clipDirty = true;

    if (autolevels->get_active() && listener) {
        Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::clip_changed_));
    }
}

bool ToneCurve::clip_changed_ ()
{

    if (listener) {
        listener->panelChanged (EvClip, Glib::ustring::format (std::setprecision(5), sclip->get_value()));

        if (!batchMode) {
            waitForAutoExp ();
        }
    }

    return false;
}

void ToneCurve::waitForAutoExp ()
{

    sclip->set_sensitive (false);
    expcomp->setEnabled (false);
    brightness->setEnabled (false);
    contrast->setEnabled (false);
    black->setEnabled (false);
    hlcompr->setEnabled (false);
    hlcomprthresh->setEnabled (false);
    shcompr->setEnabled (false);
    contrast->setEnabled (false);
    saturation->setEnabled (false);
    curveEditorG->set_sensitive (false);
    toneCurveMode->set_sensitive (false);
    curveEditorG2->set_sensitive (false);
    toneCurveMode2->set_sensitive (false);
    hrenabled->set_sensitive(false);
    method->set_sensitive(false);
}

int autoExpChangedUI (void* data)
{
    (static_cast<ToneCurve*>(data))->autoExpComputed_ ();
    return 0;
}

void ToneCurve::autoExpChanged (double expcomp, int bright, int contr, int black, int hlcompr, int hlcomprthresh, bool hlrecons)
{

    nextBlack = black;
    nextExpcomp = expcomp;
    nextBrightness = bright;
    nextContrast = contr;
    nextHlcompr = hlcompr;
    nextHlcomprthresh = hlcomprthresh;
    nextHLRecons = hlrecons;
    g_idle_add (autoExpChangedUI, this);
}

void ToneCurve::enableAll ()
{

    sclip->set_sensitive (true);
    expcomp->setEnabled (true);
    brightness->setEnabled (true);
    black->setEnabled (true);
    hlcompr->setEnabled (true);
    hlcomprthresh->setEnabled (true);
    shcompr->setEnabled (true);
    contrast->setEnabled (true);
    saturation->setEnabled (true);
    curveEditorG->set_sensitive (true);
    toneCurveMode->set_sensitive (true);
    curveEditorG2->set_sensitive (true);
    toneCurveMode2->set_sensitive (true);
    hrenabled->set_sensitive(true);
    method->set_sensitive(true);
}

bool ToneCurve::autoExpComputed_ ()
{

    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    disableListener ();
    enableAll ();
    expcomp->setValue (nextExpcomp);
    brightness->setValue (nextBrightness);
    contrast->setValue (nextContrast);
    black->setValue (nextBlack);
    hlcompr->setValue (nextHlcompr);
    hlcomprthresh->setValue (nextHlcomprthresh);
    enaconn.block (true);
    hrenabled->set_active (nextHLRecons);
    enaconn.block (false);

    if (nextHLRecons) {
        hlrbox->show();
    } else if (!batchMode) {
        hlrbox->hide();
    }

    if (!black->getAddMode()) {
        shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
    }

    enableListener ();

    return false;
}

void ToneCurve::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    method->append (M("GENERAL_UNCHANGED"));

    removeIfThere (abox, autolevels, false);
    autolevels = Gtk::manage (new Gtk::CheckButton (M("TP_EXPOSURE_AUTOLEVELS")));
    autolevels->set_tooltip_markup (M("TP_EXPOSURE_AUTOLEVELS_TIP"));
    autoconn = autolevels->signal_toggled().connect( sigc::mem_fun(*this, &ToneCurve::autolevels_toggled) );
    abox->pack_start (*autolevels);

    ToolPanel::setBatchMode (batchMode);
    expcomp->showEditedCB ();
    black->showEditedCB ();
    hlcompr->showEditedCB ();
    hlcomprthresh->showEditedCB ();
    shcompr->showEditedCB ();
    brightness->showEditedCB ();
    contrast->showEditedCB ();
    saturation->showEditedCB ();

    toneCurveMode->append (M("GENERAL_UNCHANGED"));
    toneCurveMode2->append (M("GENERAL_UNCHANGED"));

    curveEditorG->setBatchMode (batchMode);
    curveEditorG2->setBatchMode (batchMode);
}

void ToneCurve::setAdjusterBehavior (bool expadd, bool hlcompadd, bool hlcompthreshadd, bool bradd, bool blackadd, bool shcompadd, bool contradd, bool satadd)
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

void ToneCurve::trimValues (rtengine::procparams::ProcParams* pp)
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

void ToneCurve::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve, /*LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM, LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histLRETI)
{

    shape->updateBackgroundHistogram (histToneCurve);
}
