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
#include <tonecurve.h>
#include <adjuster.h>
#include <sigc++/class_slot.h>
#include <iomanip>
#include <guiutils.h>

using namespace rtengine;
using namespace rtengine::procparams;

ToneCurve::ToneCurve () : ToolPanel(), expAdd(false), blackAdd(false), brAdd(false), contrAdd(false) {

//----------- Auto Levels ----------------------------------
  abox = Gtk::manage (new Gtk::HBox ());
  abox->set_border_width (2);

  autolevels = Gtk::manage (new Gtk::ToggleButton (M("TP_EXPOSURE_AUTOLEVELS")));
  autoconn = autolevels->signal_toggled().connect( sigc::mem_fun(*this, &ToneCurve::autolevels_toggled) );

  sclip = Gtk::manage (new Gtk::SpinButton ());
  sclip->set_range (0.0, 0.9999);
  sclip->set_increments (0.0001, 0.01);
  sclip->set_value (0.002);
  sclip->set_digits (4);
  sclip->signal_value_changed().connect( sigc::mem_fun(*this, &ToneCurve::clip_changed) );

  abox->pack_start (*autolevels);
  abox->pack_end (*sclip);
  abox->pack_end (*Gtk::manage (new Gtk::Label (M("TP_EXPOSURE_CLIP"))));
  pack_start (*abox);

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

//----------- Exposure Compensation ------------------------
  expcomp   = new Adjuster (M("TP_EXPOSURE_EXPCOMP"), -5, 5, 0.01, 0);
  pack_start (*expcomp);
  hlcompr = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 150, 1, 0));
  pack_start (*hlcompr);

//----------- Black Level ----------------------------------
  black = Gtk::manage (new Adjuster (M("TP_EXPOSURE_BLACKLEVEL"), 0, 32768, 1, 0));
  pack_start (*black);
  shcompr = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRSHADOWS"), 0, 150, 1, 0));
  pack_start (*shcompr);

  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

//---------Brightness / Contrast -------------------------
  brightness = Gtk::manage (new Adjuster (M("TP_EXPOSURE_BRIGHTNESS"), -100, 100, 1, 0));
  pack_start (*brightness);
  contrast   = Gtk::manage (new Adjuster (M("TP_EXPOSURE_CONTRAST"), -100, 100, 1, 0));
  pack_start (*contrast);

/*  
//----------- Curve ------------------------------
  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  shape = Gtk::manage (new CurveEditor ());
  shape->setCurveListener (this);
  curvexp = Gtk::manage (new Gtk::Expander (M("TP_EXPOSURE_CURVEEDITOR")));
  curvexp->add (*shape);

  pack_start (*curvexp, Gtk::PACK_SHRINK, 4);
*/
// --------- Set Up Listeners -------------
  expcomp->setAdjusterListener (this);
  brightness->setAdjusterListener (this);
  black->setAdjusterListener (this);
  hlcompr->setAdjusterListener (this);
  shcompr->setAdjusterListener (this);
  contrast->setAdjusterListener (this);
}

ToneCurve::~ToneCurve () {

    delete expcomp;
}

void ToneCurve::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        expcomp->setEditedState (pedited->toneCurve.expcomp ? Edited : UnEdited);
        black->setEditedState (pedited->toneCurve.black ? Edited : UnEdited);
        hlcompr->setEditedState (pedited->toneCurve.hlcompr ? Edited : UnEdited);
        shcompr->setEditedState (pedited->toneCurve.shcompr ? Edited : UnEdited);
        brightness->setEditedState (pedited->toneCurve.brightness ? Edited : UnEdited);
        contrast->setEditedState (pedited->toneCurve.contrast ? Edited : UnEdited);
        autolevels->set_inconsistent (!pedited->toneCurve.autoexp);
        clipDirty = pedited->toneCurve.clip;
    }

    autoconn.block (true);
    autolevels->set_active (pp->toneCurve.autoexp);
    autoconn.block (false);
    lastAuto = pp->toneCurve.autoexp;
    sclip->set_value (pp->toneCurve.clip);

    expcomp->setValue (pp->toneCurve.expcomp);
    black->setValue (pp->toneCurve.black);
    hlcompr->setValue (pp->toneCurve.hlcompr);
    shcompr->setValue (pp->toneCurve.shcompr);
    brightness->setValue (pp->toneCurve.brightness);
    contrast->setValue (pp->toneCurve.contrast);
//    shape->setCurve (pp->toneCurve.curve);

    enableListener ();
}

void ToneCurve::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->toneCurve.autoexp = autolevels->get_active();
    pp->toneCurve.clip = sclip->get_value ();
    pp->toneCurve.expcomp = expcomp->getValue ();
    pp->toneCurve.black = (int)black->getValue ();
    pp->toneCurve.hlcompr = (int)hlcompr->getValue ();
    pp->toneCurve.shcompr = (int)shcompr->getValue ();
    pp->toneCurve.brightness = (int)brightness->getValue ();
    pp->toneCurve.contrast = (int)contrast->getValue ();
//    pp->toneCurve.curve = shape->getCurve ();

    if (pedited) {
        pedited->toneCurve.expcomp = expcomp->getEditedState ();
        pedited->toneCurve.black   = black->getEditedState ();
        pedited->toneCurve.hlcompr = hlcompr->getEditedState ();
        pedited->toneCurve.shcompr = shcompr->getEditedState ();
        pedited->toneCurve.brightness = brightness->getEditedState ();
        pedited->toneCurve.contrast = contrast->getEditedState ();
        pedited->toneCurve.autoexp   = !autolevels->get_inconsistent();
        pedited->toneCurve.clip   = clipDirty;
    }
}

void ToneCurve::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    expcomp->setDefault (defParams->toneCurve.expcomp);
    brightness->setDefault (defParams->toneCurve.brightness);
    black->setDefault (defParams->toneCurve.black);
    hlcompr->setDefault (defParams->toneCurve.hlcompr);
    shcompr->setDefault (defParams->toneCurve.shcompr);
    contrast->setDefault (defParams->toneCurve.contrast);

    if (pedited) {
        expcomp->setDefaultEditedState (pedited->toneCurve.expcomp ? Edited : UnEdited);
        black->setDefaultEditedState (pedited->toneCurve.black ? Edited : UnEdited);
        hlcompr->setDefaultEditedState (pedited->toneCurve.hlcompr ? Edited : UnEdited);
        shcompr->setDefaultEditedState (pedited->toneCurve.shcompr ? Edited : UnEdited);
        brightness->setDefaultEditedState (pedited->toneCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->toneCurve.contrast ? Edited : UnEdited);
    }
    else {
        expcomp->setDefaultEditedState (Irrelevant);
        black->setDefaultEditedState (Irrelevant);
        hlcompr->setDefaultEditedState (Irrelevant);
        shcompr->setDefaultEditedState (Irrelevant);
        brightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
    }
}

/*void ToneCurve::curveChanged () {

    if (listener) {
        listener->panelChanged (EvToneCurve, M("HISTORY_CUSTOMCURVE"));
    }
}
*/
void ToneCurve::adjusterChanged (Adjuster* a, double newval) {

    if (autolevels->get_active() && (a==expcomp || a==black || a==hlcompr || a==shcompr)) {
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
    else if (a==black) 
        listener->panelChanged (EvBlack, costr);
    else if (a==contrast)
        listener->panelChanged (EvContrast, costr);
    else if (a==hlcompr)
        listener->panelChanged (EvHLCompr, costr);
    else if (a==shcompr)
        listener->panelChanged (EvSHCompr, costr);
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
    } 
    
    if (batchMode) {
        expcomp->setEditedState (UnEdited);
        black->setEditedState (UnEdited);
        if (expAdd)
            expcomp->setValue (0);
        if (blackAdd)
            black->setValue (0);
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
    black->setEnabled (false);
    hlcompr->setEnabled (false);
    shcompr->setEnabled (false);
    contrast->setEnabled (false);
//    shape->set_sensitive (false);
}

int aexpcomputed (void* data) {

    gdk_threads_enter();
    ((ToneCurve*)data)->autoExpComputed_ ();
    gdk_threads_leave();
    return 0;
}

void ToneCurve::autoExpChanged (double br, int bl) {

    nextBl = bl;
    nextBr = br;
    g_idle_add (aexpcomputed, this);


//    Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::autoExpComputed_));
}

void ToneCurve::enableAll () {

    sclip->set_sensitive (true);
    expcomp->setEnabled (true);
    brightness->setEnabled (true);
    black->setEnabled (true);
    hlcompr->setEnabled (true);
    shcompr->setEnabled (true);
    contrast->setEnabled (true);
//     shape->set_sensitive (true);
}

bool ToneCurve::autoExpComputed_ () {

    disableListener ();
    enableAll ();
    expcomp->setValue (nextBr);
    black->setValue (nextBl);
    enableListener ();

    return false;
}

/*
void ToneCurve::expandCurve (bool isExpanded) {

    curvexp->set_expanded (isExpanded);
}

bool ToneCurve::isCurveExpanded () {

    return curvexp->get_expanded ();
}
*/

void ToneCurve::setBatchMode (bool batchMode) {

    removeIfThere (abox, autolevels, false);
    autolevels = Gtk::manage (new Gtk::CheckButton (M("TP_EXPOSURE_AUTOLEVELS")));
    autoconn = autolevels->signal_toggled().connect( sigc::mem_fun(*this, &ToneCurve::autolevels_toggled) );
    abox->pack_start (*autolevels);

    ToolPanel::setBatchMode (batchMode);
    expcomp->showEditedCB ();
    black->showEditedCB ();
    hlcompr->showEditedCB ();
    shcompr->showEditedCB ();
    brightness->showEditedCB ();
    contrast->showEditedCB ();
}

void ToneCurve::setAdjusterBehavior (bool expadd, bool bradd, bool blackadd, bool contradd) {

    if (!expAdd && expadd || expAdd && !expadd)
        expcomp->setLimits (-5, 5, 0.01, 0);
    if (!blackAdd && blackadd)
        black->setLimits (0, 16384, 1, 0);
    else if (blackAdd && !blackadd)
        black->setLimits (0, 32768, 1, 0);
    if (!brAdd && bradd || brAdd && !bradd)
        brightness->setLimits (-100, 100, 0.01, 0);
    if (!contrAdd && contradd || contrAdd && !contradd)
        contrast->setLimits (-100, 100, 0.01, 0);
    
    expAdd = expadd;
    blackAdd = blackadd;
    brAdd = bradd;
    contrAdd = contradd;
}

