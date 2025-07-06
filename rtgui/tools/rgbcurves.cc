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
#include "rgbcurves.h"

#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "options.h"

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring RGBCurves::TOOL_NAME = "rgbcurves";

RGBCurves::RGBCurves () : FoldableToolPanel(this, TOOL_NAME, M("TP_RGBCURVES_LABEL"), false, true), lastLumamode(false)
{

    lumamode = Gtk::manage (new Gtk::CheckButton (M("TP_RGBCURVES_LUMAMODE")));
    lumamode->set_tooltip_markup (M("TP_RGBCURVES_LUMAMODE_TOOLTIP"));
    lumamode->set_active (false);
    lumamode->show ();
    pack_start (*lumamode);

    Gtk::Separator *hsep1 = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    hsep1->show ();
    pack_start (*hsep1);

    lumamodeConn = lumamode->signal_toggled().connect( sigc::mem_fun(*this, &RGBCurves::lumamodeChanged) );

    std::vector<GradientMilestone> milestones;

    curveEditorG = new CurveEditorGroup (options.lastRgbCurvesDir, M("TP_RGBCURVES_CHANNEL"));
    curveEditorG->setCurveListener (this);

    Rshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_RGBCURVES_RED")));
    Rshape->setEditID(EUID_RGB_R, BT_SINGLEPLANE_FLOAT);
    milestones.push_back( GradientMilestone(0.0, 0.0, 0.0, 0.0) );
    milestones.push_back( GradientMilestone(1.0, 1.0, 0.0, 0.0) );
    Rshape->setBottomBarBgGradient(milestones);
    Rshape->setLeftBarBgGradient(milestones);

    milestones[1].r = 0.0;
    milestones[1].g = 1.0;
    Gshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_RGBCURVES_GREEN")));
    Gshape->setEditID(EUID_RGB_G, BT_SINGLEPLANE_FLOAT);
    Gshape->setBottomBarBgGradient(milestones);
    Gshape->setLeftBarBgGradient(milestones);

    milestones[1].g = 0.0;
    milestones[1].b = 1.0;
    Bshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_RGBCURVES_BLUE")));
    Bshape->setEditID(EUID_RGB_B, BT_SINGLEPLANE_FLOAT);
    Bshape->setBottomBarBgGradient(milestones);
    Bshape->setLeftBarBgGradient(milestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG->curveListComplete();

    pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);

}

RGBCurves::~RGBCurves ()
{
    delete curveEditorG;
}

void RGBCurves::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        Rshape->setUnChanged (!pedited->rgbCurves.rcurve);
        Gshape->setUnChanged (!pedited->rgbCurves.gcurve);
        Bshape->setUnChanged (!pedited->rgbCurves.bcurve);
        lumamode->set_inconsistent (!pedited->rgbCurves.lumamode);
        set_inconsistent(multiImage && !pedited->rgbCurves.enabled);
    }

    lumamodeConn.block (true);
    lumamode->set_active (pp->rgbCurves.lumamode);
    lumamodeConn.block (false);

    lastLumamode = pp->rgbCurves.lumamode;

    Rshape->setCurve         (pp->rgbCurves.rcurve);
    Gshape->setCurve         (pp->rgbCurves.gcurve);
    Bshape->setCurve         (pp->rgbCurves.bcurve);

    setEnabled(pp->rgbCurves.enabled);

    enableListener ();
}

void RGBCurves::setEditProvider (EditDataProvider *provider)
{
    Rshape->setEditProvider(provider);
    Gshape->setEditProvider(provider);
    Bshape->setEditProvider(provider);
}

void RGBCurves::autoOpenCurve  ()
{
    // Open up the first curve if selected
    bool active = Rshape->openIfNonlinear();

    if (!active) {
        active = Gshape->openIfNonlinear();
    }

    if (!active) {
        Bshape->openIfNonlinear();
    }
}

void RGBCurves::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->rgbCurves.enabled = getEnabled();
    pp->rgbCurves.rcurve         = Rshape->getCurve ();
    pp->rgbCurves.gcurve         = Gshape->getCurve ();
    pp->rgbCurves.bcurve         = Bshape->getCurve ();
    pp->rgbCurves.lumamode       = lumamode->get_active();

    if (pedited) {
        pedited->rgbCurves.enabled = !get_inconsistent();
        pedited->rgbCurves.rcurve    = !Rshape->isUnChanged ();
        pedited->rgbCurves.gcurve    = !Gshape->isUnChanged ();
        pedited->rgbCurves.bcurve    = !Bshape->isUnChanged ();
        pedited->rgbCurves.lumamode  = !lumamode->get_inconsistent();
    }
}


/*
 * Curve listener
 *
 * If more than one curve has been added, the curve listener is automatically
 * set to 'multi=true', and send a pointer of the modified curve in a parameter
 */
void RGBCurves::curveChanged (CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == Rshape) {
            listener->panelChanged (EvRGBrCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == Gshape) {
            listener->panelChanged (EvRGBgCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == Bshape) {
            listener->panelChanged (EvRGBbCurve, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void RGBCurves::lumamodeChanged ()
{

    if (batchMode) {
        if (lumamode->get_inconsistent()) {
            lumamode->set_inconsistent (false);
            lumamodeConn.block (true);
            lumamode->set_active (false);
            lumamodeConn.block (false);
        } else if (lastLumamode) {
            lumamode->set_inconsistent (true);
        }

        lastLumamode = lumamode->get_active ();
    }

    if (listener && getEnabled()) {
        if (lumamode->get_active ()) {
            listener->panelChanged (EvRGBrCurveLumamode, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvRGBrCurveLumamode, M("GENERAL_DISABLED"));
        }
    }
}

void RGBCurves::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    curveEditorG->setBatchMode (batchMode);
}


void RGBCurves::updateCurveBackgroundHistogram(
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
}


void RGBCurves::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvRGBEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvRGBEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvRGBEnabled, M("GENERAL_DISABLED"));
        }
    }
}
