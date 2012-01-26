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
#include <sigc++/class_slot.h>
#include <iomanip>
#include "ppversion.h"

using namespace rtengine;
using namespace rtengine::procparams;

ToneCurve::ToneCurve () : Gtk::VBox(), FoldableToolPanel(this) {

//----------- Auto Levels ----------------------------------
  abox = Gtk::manage (new Gtk::HBox ());
  abox->set_border_width (2);

  autolevels = Gtk::manage (new Gtk::ToggleButton (M("TP_EXPOSURE_AUTOLEVELS")));
  autolevels->set_tooltip_markup (M("TP_EXPOSURE_AUTOLEVELS_TIP"));
  autoconn = autolevels->signal_toggled().connect( sigc::mem_fun(*this, &ToneCurve::autolevels_toggled) );

  sclip = Gtk::manage (new MySpinButton ());
  sclip->set_range (0.0, 0.9999);
  sclip->set_increments (0.001, 0.01);
  sclip->set_value (0.002);
  sclip->set_digits (4);
  sclip->set_tooltip_text (M("TP_EXPOSURE_CLIP_TIP"));
  sclip->signal_value_changed().connect( sigc::mem_fun(*this, &ToneCurve::clip_changed) );

  neutral = Gtk::manage (new Gtk::Button (M("TP_NEUTRAL")));
  neutral->set_tooltip_text (M("TP_NEUTRAL_TIP"));
  neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &ToneCurve::neutral_pressed) );
  neutral->show();

  abox->pack_start (*autolevels);
  // pack_end is used for these controls as autolevels is replaceable using pack_start in batchmode
  abox->pack_end (*neutral);
  abox->pack_end (*Gtk::manage (new Gtk::Label (" "))); //spacer
  abox->pack_end (*sclip);
  abox->pack_end (*Gtk::manage (new Gtk::Label (M("TP_EXPOSURE_CLIP"))));
  pack_start (*abox);

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

//----------- Exposure Compensation ------------------------
  expcomp   = Gtk::manage (new Adjuster (M("TP_EXPOSURE_EXPCOMP"), -5, 10, 0.05, 0));
  pack_start (*expcomp);
  hlcompr = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 70));
  pack_start (*hlcompr);
  hlcomprthresh = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 0));
  pack_start (*hlcomprthresh);

//----------- Black Level ----------------------------------
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
  
//----------- Curve ------------------------------
  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  curveEditorG = new CurveEditorGroup (M("TP_EXPOSURE_CURVEEDITOR"));
  curveEditorG->setCurveListener (this);

  shape = (DiagonalCurveEditor*)curveEditorG->addCurve(CT_Diagonal, "");

  // This will add the reset button at the end of the curveType buttons
  curveEditorG->curveListComplete();

  pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);

  //curveEditorG->show();

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

void ToneCurve::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

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
    }

    autoconn.block (true);
    autolevels->set_active (pp->toneCurve.autoexp);
    autoconn.block (false);
    lastAuto = pp->toneCurve.autoexp;
    sclip->set_value (pp->toneCurve.clip);

    expcomp->setValue (pp->toneCurve.expcomp);
    black->setValue (pp->toneCurve.black);
    hlcompr->setValue (pp->toneCurve.hlcompr);
    hlcomprthresh->setValue (pp->toneCurve.hlcomprthresh);
    shcompr->setValue (pp->toneCurve.shcompr);
    if (!black->getAddMode()) shcompr->set_sensitive(!((int)black->getValue ()==0)); //at black=0 shcompr value has no effect
    brightness->setValue (pp->toneCurve.brightness);
    contrast->setValue (pp->toneCurve.contrast);
	saturation->setValue (pp->toneCurve.saturation);
	shape->setCurve (pp->toneCurve.curve);

    shape->openIfNonlinear();

    enableListener ();
}

void ToneCurve::write (ProcParams* pp, ParamsEdited* pedited) {

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

    if (pedited) {
        pedited->toneCurve.expcomp = expcomp->getEditedState ();
        pedited->toneCurve.black   = black->getEditedState ();
        pedited->toneCurve.hlcompr = hlcompr->getEditedState ();
        pedited->toneCurve.hlcomprthresh = hlcomprthresh->getEditedState ();
        pedited->toneCurve.shcompr = shcompr->getEditedState ();
        pedited->toneCurve.brightness = brightness->getEditedState ();
        pedited->toneCurve.contrast = contrast->getEditedState ();
		pedited->toneCurve.saturation = saturation->getEditedState ();
		pedited->toneCurve.autoexp  = !autolevels->get_inconsistent();
        pedited->toneCurve.clip     = clipDirty;
        pedited->toneCurve.curve    = !shape->isUnChanged ();
    }
}

void ToneCurve::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

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
    }
    else {
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

void ToneCurve::curveChanged () {

    if (listener) listener->panelChanged (EvToneCurve, M("HISTORY_CUSTOMCURVE"));
}

void ToneCurve::adjusterChanged (Adjuster* a, double newval) {

    // Switch off auto exposure if user changes sliders manually
    if (autolevels->get_active() && (a==expcomp || a==brightness || a==contrast || a==black || a==hlcompr || a==hlcomprthresh)) {
        autolevels->set_active (false);
        autolevels->set_inconsistent (false);
    }

    if (!listener)
        return;

    Glib::ustring costr;
    if (a==expcomp)
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
    else
        costr = Glib::ustring::format ((int)a->getValue());
    
    if (a==expcomp) 
        listener->panelChanged (EvExpComp, costr);
    else if (a==brightness) 
        listener->panelChanged (EvBrightness, costr);
    else if (a==black){
        listener->panelChanged (EvBlack, costr);
        if (!black->getAddMode()) shcompr->set_sensitive(!((int)black->getValue ()==0)); //at black=0 shcompr value has no effect
    }
    else if (a==contrast)
        listener->panelChanged (EvContrast, costr);
	else if (a==saturation)
        listener->panelChanged (EvSaturation, costr);
    else if (a==hlcompr)
        listener->panelChanged (EvHLCompr, costr);
    else if (a==hlcomprthresh)
        listener->panelChanged (EvHLComprThreshold, costr);
    else if (a==shcompr)
        listener->panelChanged (EvSHCompr, costr);
}

void ToneCurve::neutral_pressed () {
// This method deselects auto levels 
// and sets neutral values to params in exposure panel

    if (batchMode) {
        autolevels->set_inconsistent (false);
        autoconn.block (true);
        autolevels->set_active (false);
        autoconn.block (false);

        lastAuto = autolevels->get_active ();
    }
    else { //!batchMode
        autolevels->set_active (false);
        autolevels->set_inconsistent (false);
    }
    
    expcomp->setValue(0);
    hlcompr->setValue(0);
    hlcomprthresh->setValue(0);
    brightness->setValue(0);
    black->setValue(0);
    shcompr->setValue(0);
    if (!black->getAddMode()) shcompr->set_sensitive(!((int)black->getValue ()==0)); //at black=0 shcompr value has no effect
    contrast->setValue(0);
    //saturation->setValue(0);
    
    listener->panelChanged (EvNeutralExp, M("GENERAL_ENABLED"));
}
void ToneCurve::autolevels_toggled () {

    if (batchMode) {
        if (autolevels->get_inconsistent()) {
            autolevels->set_inconsistent (false);
            autoconn.block (true);
            autolevels->set_active (false);
            autoconn.block (false);
        }
        else if (lastAuto)
            autolevels->set_inconsistent (true);

        lastAuto = autolevels->get_active ();
    }

    if (!batchMode && autolevels->get_active() && listener) {
        listener->panelChanged (EvAutoExp, M("GENERAL_ENABLED"));
        waitForAutoExp ();
        if (!black->getAddMode()) shcompr->set_sensitive(!((int)black->getValue ()==0)); //at black=0 shcompr value has no effect
    } 
    
    if (batchMode) {
        expcomp->setEditedState (UnEdited);
		brightness->setEditedState (UnEdited);
		contrast->setEditedState (UnEdited);
        black->setEditedState (UnEdited);
        hlcompr->setEditedState (UnEdited);
        hlcomprthresh->setEditedState (UnEdited);
        if (expcomp->getAddMode())
            expcomp->setValue (0);
		if (brightness->getAddMode())
            brightness->setValue (0);
		if (contrast->getAddMode())
            contrast->setValue (0);
		if (black->getAddMode())
            black->setValue (0);
        if (hlcompr->getAddMode())
        	hlcompr->setValue (0);
        if (hlcomprthresh->getAddMode())
        	hlcomprthresh->setValue (0);
        listener->panelChanged (EvAutoExp, M("GENERAL_ENABLED"));
    }
}

void ToneCurve::clip_changed () {

    clipDirty = true;
    if (autolevels->get_active() && listener)
        Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::clip_changed_));
}

bool ToneCurve::clip_changed_ () {

    if (listener) {
        listener->panelChanged (EvClip, Glib::ustring::format (std::setprecision(5), sclip->get_value()));
        if (!batchMode)
            waitForAutoExp ();
    }
    return false;
}

void ToneCurve::waitForAutoExp () {

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
}

int autoExpChangedUI (void* data) {
    ((ToneCurve*)data)->autoExpComputed_ ();
    return 0;
}

void ToneCurve::autoExpChanged (double expcomp, int bright, int contr, int black, int hlcompr, int hlcomprthresh) {

    nextBlack = black;
    nextExpcomp = expcomp;
	nextBrightness = bright;
	nextContrast = contr;
    nextHlcompr = hlcompr;
    nextHlcomprthresh = hlcomprthresh;
    g_idle_add (autoExpChangedUI, this);

//    Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::autoExpComputed_));
}

void ToneCurve::enableAll () {

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
}

bool ToneCurve::autoExpComputed_ () {

    disableListener ();
    enableAll ();
    expcomp->setValue (nextExpcomp);
	brightness->setValue (nextBrightness);
	contrast->setValue (nextContrast);
    black->setValue (nextBlack);
    hlcompr->setValue (nextHlcompr);
    hlcomprthresh->setValue (nextHlcomprthresh);
    if (!black->getAddMode()) shcompr->set_sensitive(!((int)black->getValue ()==0)); //at black=0 shcompr value has no effect
    enableListener ();

    return false;
}

void ToneCurve::setBatchMode (bool batchMode) {

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

    curveEditorG->setBatchMode (batchMode);
}

void ToneCurve::setAdjusterBehavior (bool expadd, bool hlcompadd, bool hlcompthreshadd, bool bradd, bool blackadd, bool shcompadd, bool contradd, bool satadd) {

	expcomp->setAddMode(expadd);
	hlcompr->setAddMode(hlcompadd);
	hlcomprthresh->setAddMode(hlcompthreshadd);
	brightness->setAddMode(bradd);
	black->setAddMode(blackadd);
	shcompr->setAddMode(shcompadd);
	contrast->setAddMode(contradd);
	saturation->setAddMode(satadd);
}

void ToneCurve::trimValues (rtengine::procparams::ProcParams* pp) {

	expcomp->trimValue(pp->toneCurve.expcomp);
	hlcompr->trimValue(pp->toneCurve.hlcompr);
	hlcomprthresh->trimValue(pp->toneCurve.hlcomprthresh);
	brightness->trimValue(pp->toneCurve.brightness);
	black->trimValue(pp->toneCurve.black);
	shcompr->trimValue(pp->toneCurve.shcompr);
	contrast->trimValue(pp->toneCurve.contrast);
	saturation->trimValue(pp->toneCurve.saturation);
}

void ToneCurve::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma) {

    shape->updateBackgroundHistogram (histToneCurve);
}
