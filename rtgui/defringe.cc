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
#include "defringe.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

Defringe::Defringe () : FoldableToolPanel(this, "defringe", M("TP_DEFRINGE_LABEL"), true, true)
{

    std::vector<GradientMilestone> bottomMilestones;
    float R, G, B;

    for (int i = 0; i < 7; i++) {
        float x = float(i) * (1.0f / 6.0);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        bottomMilestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
    }

    curveEditorPF = new CurveEditorGroup (options.lastPFCurvesDir);
    curveEditorPF->setCurveListener (this);
    chshape = static_cast<FlatCurveEditor*>(curveEditorPF->addCurve(CT_Flat, M("TP_PFCURVE_CURVEEDITOR_CH")));
    chshape->setTooltip(M("TP_PFCURVE_CURVEEDITOR_CH_TOOLTIP"));
    chshape->setIdentityValue(0.);
    chshape->setBottomBarBgGradient(bottomMilestones);
    chshape->setCurveColorProvider(this, 1);

    //edgConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Defringe::edgeChanged) );

    radius  = Gtk::manage (new Adjuster (M("TP_DEFRINGE_RADIUS"), 0.5, 5.0, 0.1, 2.0));
    threshold    = Gtk::manage (new Adjuster (M("TP_DEFRINGE_THRESHOLD"), 0, 100, 1, 13));
    radius->setAdjusterListener (this);
    threshold->setAdjusterListener (this);
// radius->show();
// threshold->show();

    pack_start (*radius);
    pack_start (*threshold);
    curveEditorPF->curveListComplete();

    pack_start (*curveEditorPF, Gtk::PACK_SHRINK, 4);

}

Defringe::~Defringe ()
{
    delete curveEditorPF;
}

void Defringe::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R = 0.f, G = 0.f, B = 0.f;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) { // ch
        Color::hsv2rgb01(float(valX), float(valY), 0.5f, R, G, B);
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}

void Defringe::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        radius->setEditedState    ( pedited->defringe.radius ? Edited : UnEdited);
        threshold->setEditedState ( pedited->defringe.threshold ? Edited : UnEdited);
        set_inconsistent          (multiImage && !pedited->defringe.enabled);
        chshape->setUnChanged     (!pedited->defringe.huecurve);
    }

    setEnabled(pp->defringe.enabled);

    radius->setValue    (pp->defringe.radius);
    threshold->setValue (pp->defringe.threshold);
    chshape->setCurve   (pp->defringe.huecurve);

    enableListener ();
}


void Defringe::autoOpenCurve ()
{
    // WARNING: The following line won't work, since linear is a flat curve at 0.
    // chshape->openIfNonlinear();
}

void Defringe::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->defringe.radius    = radius->getValue ();
    pp->defringe.threshold = (int)threshold->getValue ();
    pp->defringe.enabled   = getEnabled();
    pp->defringe.huecurve  = chshape->getCurve ();

    if (pedited) {
        pedited->defringe.radius    = radius->getEditedState ();
        pedited->defringe.threshold = threshold->getEditedState ();
        pedited->defringe.enabled   = !get_inconsistent();
        pedited->defringe.huecurve  = !chshape->isUnChanged ();
    }
}

void Defringe::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    radius->setDefault (defParams->defringe.radius);
    threshold->setDefault (defParams->defringe.threshold);

    if (pedited) {
        radius->setDefaultEditedState (pedited->defringe.radius ? Edited : UnEdited);
        threshold->setDefaultEditedState (pedited->defringe.threshold ? Edited : UnEdited);
    } else {
        radius->setDefaultEditedState (Irrelevant);
        threshold->setDefaultEditedState (Irrelevant);
    }
}
void Defringe::curveChanged ()
{

    if (listener && getEnabled()) {
        listener->panelChanged (EvPFCurve, M("HISTORY_CUSTOMCURVE"));
    }
}

void Defringe::adjusterChanged (Adjuster* a, double newval)
{

    if (listener && getEnabled()) {

        if (a == radius) {
            listener->panelChanged (EvDefringeRadius, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
        } else if (a == threshold) {
            listener->panelChanged (EvDefringeThreshold, Glib::ustring::format ((int)a->getValue()));
        }
    }
}

void Defringe::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvDefringeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvDefringeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDefringeEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Defringe::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    radius->showEditedCB ();
    threshold->showEditedCB ();
}
