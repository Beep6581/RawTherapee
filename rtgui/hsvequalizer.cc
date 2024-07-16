/*
 *  This file is part of RawTherapee.
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
 *
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 */
#include "hsvequalizer.h"

#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "options.h"

#include "../rtengine/color.h"
#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring HSVEqualizer::TOOL_NAME = "hsvequalizer";

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HSVEqualizer::HSVEqualizer () : FoldableToolPanel(this, TOOL_NAME, M("TP_HSVEQUALIZER_LABEL"), false, true)
{

    std::vector<GradientMilestone> bottomMilestones;
    float R, G, B;

    // -0.1 rad < Hue < 1.6 rad
    for (int i = 0; i < 7; i++) {
        float x = i / 6.0;
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        bottomMilestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
    }

    curveEditorG = new CurveEditorGroup (options.lastHsvCurvesDir, M("TP_HSVEQUALIZER_CHANNEL"));
    curveEditorG->setCurveListener (this);

    hshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_HUE")));
    hshape->setEditID(EUID_HSV_H, BT_SINGLEPLANE_FLOAT);
    hshape->setBottomBarBgGradient(bottomMilestones);
    //hshape->setLeftBarColorProvider(this);  Not working yet
    hshape->setCurveColorProvider(this, 1);

    sshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_SAT")));
    sshape->setEditID(EUID_HSV_S, BT_SINGLEPLANE_FLOAT);
    sshape->setBottomBarBgGradient(bottomMilestones);
    //sshape->setLeftBarColorProvider(this);  Not working yet
    sshape->setCurveColorProvider(this, 2);

    vshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_VAL")));
    vshape->setEditID(EUID_HSV_V, BT_SINGLEPLANE_FLOAT);
    vshape->setBottomBarBgGradient(bottomMilestones);
    //vshape->setLeftBarColorProvider(this);  Not working yet
    vshape->setCurveColorProvider(this, 3);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG->curveListComplete();

    pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);

    //curveEditorG->show();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HSVEqualizer::~HSVEqualizer ()
{
    delete curveEditorG;
}


void HSVEqualizer::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        hshape->setUnChanged (!pedited->hsvequalizer.hcurve);
        sshape->setUnChanged (!pedited->hsvequalizer.scurve);
        vshape->setUnChanged (!pedited->hsvequalizer.vcurve);
        set_inconsistent(multiImage && !pedited->hsvequalizer.enabled);
    }

    hshape->setCurve         (pp->hsvequalizer.hcurve);
    sshape->setCurve         (pp->hsvequalizer.scurve);
    vshape->setCurve         (pp->hsvequalizer.vcurve);
    setEnabled(pp->hsvequalizer.enabled);

    enableListener ();
}

void HSVEqualizer::setEditProvider (EditDataProvider *provider)
{
    hshape->setEditProvider(provider);
    sshape->setEditProvider(provider);
    vshape->setEditProvider(provider);
}

void HSVEqualizer::autoOpenCurve ()
{
    // Open up the first curve if selected
    bool active = hshape->openIfNonlinear();

    if (!active) {
        active = sshape->openIfNonlinear();
    }

    if (!active) {
        vshape->openIfNonlinear();
    }
}

void HSVEqualizer::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->hsvequalizer.enabled = getEnabled();
    pp->hsvequalizer.hcurve = hshape->getCurve ();
    pp->hsvequalizer.scurve = sshape->getCurve ();
    pp->hsvequalizer.vcurve = vshape->getCurve ();

    if (pedited) {

        pedited->hsvequalizer.hcurve = !hshape->isUnChanged ();
        pedited->hsvequalizer.scurve = !sshape->isUnChanged ();
        pedited->hsvequalizer.vcurve = !vshape->isUnChanged ();
        pedited->hsvequalizer.enabled = !get_inconsistent();
    }
}

/*
 * Curve listener
 *
 * If more than one curve has been added, the curve listener is automatically
 * set to 'multi=true', and send a pointer of the modified curve in a parameter
 */
void HSVEqualizer::curveChanged (CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == hshape) {
            listener->panelChanged (EvHSVEqualizerH, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == sshape) {
            listener->panelChanged (EvHSVEqualizerS, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == vshape) {
            listener->panelChanged (EvHSVEqualizerV, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void HSVEqualizer::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller)
{

    float r, g, b;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {        // Hue = f(Hue)

        float h = float((valY - 0.5) * 2. + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01(h, 0.5f, 0.5f, r, g, b);
        caller->ccRed = double(r);
        caller->ccGreen = double(g);
        caller->ccBlue = double(b);
    } else if (callerId == 2) { // Saturation = f(Hue)
        Color::hsv2rgb01(float(valX), float(valY), 0.5f, r, g, b);
        caller->ccRed = double(r);
        caller->ccGreen = double(g);
        caller->ccBlue = double(b);
    } else if (callerId == 3) { // Value = f(Hue)
        Color::hsv2rgb01(float(valX), 0.5f, float(valY), r, g, b);
        caller->ccRed = double(r);
        caller->ccGreen = double(g);
        caller->ccBlue = double(b);
    } else {
        printf("Error: no curve displayed!\n");
    }

}

void HSVEqualizer::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);

    curveEditorG->setBatchMode (batchMode);
}


void HSVEqualizer::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvHSVEqEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvHSVEqEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvHSVEqEnabled, M("GENERAL_DISABLED"));
        }
    }
}
