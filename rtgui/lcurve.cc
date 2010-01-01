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

LCurve::LCurve () : ToolPanel() {

  Gtk::HBox* abox = Gtk::manage (new Gtk::HBox ());
  abox->set_border_width (2);

  brightness = Gtk::manage (new Adjuster (M("TP_LUMACURVE_BRIGHTNESS"), -2, 2, 0.01, 0));
  hlcompr    = Gtk::manage (new Adjuster (M("TP_LUMACURVE_COMPRHIGHLIGHTS"), 0, 100, 1, 0));
  black      = Gtk::manage (new Adjuster (M("TP_LUMACURVE_BLACKLEVEL"), 0, 32768, 1, 0));
  shcompr    = Gtk::manage (new Adjuster (M("TP_LUMACURVE_COMPRSHADOWS"), 0, 100, 1, 0));
  contrast   = Gtk::manage (new Adjuster (M("TP_LUMACURVE_CONTRAST"), -50, 50, 1, 0));

  pack_start (*brightness);
  brightness->show ();

  pack_start (*hlcompr);
  hlcompr->show ();

  pack_start (*black);
  black->show ();

  pack_start (*shcompr);
  shcompr->show ();

  pack_start (*contrast);
  contrast->show ();

  Gtk::HSeparator *hsep3 = Gtk::manage (new  Gtk::HSeparator());
  hsep3->show ();
  pack_start (*hsep3);

  shape = Gtk::manage (new CurveEditor ());
  shape->show ();
  shape->setCurveListener (this);

  curvexp = Gtk::manage (new Gtk::Expander (M("TP_LUMACURVE_CURVEEDITOR")));
  curvexp->show ();
  curvexp->add (*shape);

  pack_start (*curvexp, Gtk::PACK_SHRINK, 4);

  brightness->setAdjusterListener (this);
  hlcompr->setAdjusterListener (this);
  black->setAdjusterListener (this);
  shcompr->setAdjusterListener (this);
  contrast->setAdjusterListener (this);
}

void LCurve::read (const ProcParams* pp) {

    disableListener ();

    brightness->setValue    (pp->lumaCurve.brightness);
    black->setValue         (pp->lumaCurve.black);
    hlcompr->setValue       (pp->lumaCurve.hlcompr);
    shcompr->setValue       (pp->lumaCurve.shcompr);

    contrast->setValue      (pp->lumaCurve.contrast);
    shape->setCurve         (pp->lumaCurve.curve);

    enableListener ();
}

void LCurve::write (ProcParams* pp) {

    pp->lumaCurve.brightness    = brightness->getValue ();
    pp->lumaCurve.black         = (int)black->getValue ();
    pp->lumaCurve.hlcompr       = (int)hlcompr->getValue ();
    pp->lumaCurve.shcompr       = (int)shcompr->getValue ();
    pp->lumaCurve.contrast      = (int)contrast->getValue ();
    pp->lumaCurve.curve         = shape->getCurve ();
}

void LCurve::setDefaults (const ProcParams* defParams) {

    brightness->setDefault (defParams->lumaCurve.brightness);
    black->setDefault (defParams->lumaCurve.black);
    hlcompr->setDefault (defParams->lumaCurve.hlcompr);
    shcompr->setDefault (defParams->lumaCurve.shcompr);
    contrast->setDefault (defParams->lumaCurve.contrast);
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
    else if (a==black)
        listener->panelChanged (EvLBlack, costr);
    else if (a==contrast)
        listener->panelChanged (EvLContrast, costr);
    else if (a==hlcompr)
        listener->panelChanged (EvLHLCompr, costr);
    else if (a==shcompr)
        listener->panelChanged (EvLSHCompr, costr);
}
void LCurve::expandCurve (bool isExpanded) {

    curvexp->set_expanded (isExpanded);
}

bool LCurve::isCurveExpanded () {

    return curvexp->get_expanded ();
}

