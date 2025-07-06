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
#pragma once

#include <gtkmm.h>

#include "adjuster.h"
#include "curvelistener.h"
#include "thresholdadjuster.h"
#include "toolpanel.h"

class CurveEditorGroup;
class DiagonalCurveEditor;

class Vibrance final :
    public ToolParamBlock,
    public AdjusterListener,
    public ThresholdCurveProvider,
    public ThresholdAdjusterListener,
    public FoldableToolPanel,
    public CurveListener
{

protected:
    CurveEditorGroup* curveEditorGG;

    Adjuster* pastels;
    Adjuster* saturated;
    ThresholdAdjuster* psThreshold;
    Gtk::CheckButton* protectSkins;
    Gtk::CheckButton* avoidColorShift;
    Gtk::CheckButton* pastSatTog;
    DiagonalCurveEditor* skinTonesCurve;

    bool lastProtectSkins;
    bool lastAvoidColorShift;
    bool lastPastSatTog;

    sigc::connection pskinsconn;
    sigc::connection ashiftconn;
    sigc::connection pastsattogconn;

public:
    static const Glib::ustring TOOL_NAME;

    Vibrance                 ();
    ~Vibrance                () override;

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode        (bool batchMode) override;
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void setAdjusterBehavior (bool pastelsadd, bool saturatedadd);
    void adjusterChanged     (Adjuster* a, double newval) override;
    void curveChanged        () override;
    void autoOpenCurve       () override;

    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override;
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;

    void enabledChanged          () override;
    void protectskins_toggled    ();
    void avoidcolorshift_toggled ();
    void pastsattog_toggled      ();
    std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const override;
};
