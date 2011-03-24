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
#include <lcurve.h>
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

LCurve::LCurve () : Gtk::VBox(), FoldableToolPanel(this), brAdd(false), contrAdd(false) {

  Gtk::HBox* abox = Gtk::manage (new Gtk::HBox ());
  abox->set_border_width (2);

  brightness = Gtk::manage (new Adjuster (M("TP_LUMACURVE_BRIGHTNESS"), -100, 100, 0.01, 0));
  contrast   = Gtk::manage (new Adjuster (M("TP_LUMACURVE_CONTRAST"), -100, 100, 1, 0));

  pack_start (*brightness);
  brightness->show ();

  pack_start (*contrast);
  contrast->show ();

  Gtk::HSeparator *hsep3 = Gtk::manage (new  Gtk::HSeparator());
  hsep3->show ();
  pack_start (*hsep3);

  shape = Gtk::manage (new DiagonalCurveEditor ());
  shape->show ();
  shape->setCurveListener (this);

  pack_start (*shape, Gtk::PACK_SHRINK, 4);

  brightness->setAdjusterListener (this);
  contrast->setAdjusterListener (this);
}

void LCurve::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        brightness->setEditedState (pedited->lumaCurve.brightness ? Edited : UnEdited);
        contrast->setEditedState (pedited->lumaCurve.contrast ? Edited : UnEdited);
        shape->setUnChanged (!pedited->lumaCurve.curve);
    }

    brightness->setValue    (pp->lumaCurve.brightness);
    contrast->setValue      (pp->lumaCurve.contrast);
    shape->setCurve         (pp->lumaCurve.curve);

    enableListener ();
}

void LCurve::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->lumaCurve.brightness    = brightness->getValue ();
    pp->lumaCurve.contrast      = (int)contrast->getValue ();
    pp->lumaCurve.curve         = shape->getCurve ();

    if (pedited) {
        pedited->lumaCurve.brightness = brightness->getEditedState ();
        pedited->lumaCurve.contrast = contrast->getEditedState ();
        pedited->lumaCurve.curve    = !shape->isUnChanged ();
    }
}

void LCurve::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    brightness->setDefault (defParams->lumaCurve.brightness);
    contrast->setDefault (defParams->lumaCurve.contrast);

    if (pedited) {
        brightness->setDefaultEditedState (pedited->lumaCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->lumaCurve.contrast ? Edited : UnEdited);
    }
    else {
        brightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
    }
}

void LCurve::curveChanged () {

    if (listener) 
        listener->panelChanged (EvLCurve, M("HISTORY_CUSTOMCURVE"));
}

void LCurve::adjusterChanged (Adjuster* a, double newval) {

    if (!listener)
        return;

    Glib::ustring costr;
    if (a==brightness)
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
    else
        costr = Glib::ustring::format ((int)a->getValue());

    if (a==brightness)
        listener->panelChanged (EvLBrightness, costr);
    else if (a==contrast)
        listener->panelChanged (EvLContrast, costr);
}

void LCurve::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    brightness->showEditedCB ();
    contrast->showEditedCB ();
       
    shape->setBatchMode (batchMode);
}

void LCurve::setAdjusterBehavior (bool bradd, bool contradd) {

    if ((!brAdd && bradd) || (brAdd && !bradd))
        brightness->setLimits (-100, 100, 1, 0);
    if ((!contrAdd && contradd) || (contrAdd && !contradd))
        contrast->setLimits (-100, 100, 1, 0);
    
    brAdd = bradd;
    contrAdd = contradd;
}

void LCurve::updateCurveBackgroundHistogram (unsigned* hist) {
    
    shape->updateBackgroundHistogram (hist);
}
