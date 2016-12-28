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
#ifndef _VIBRANCE_
#define _VIBRANCE_

#include <gtkmm.h>
#include "adjuster.h"
#include "thresholdadjuster.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "toolpanel.h"

class Vibrance : public ToolParamBlock, public AdjusterListener, public ThresholdCurveProvider, public ThresholdAdjusterListener,
    public FoldableToolPanel, public CurveListener
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

    Vibrance                 ();
    ~Vibrance                ();

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void setBatchMode        (bool batchMode);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void setAdjusterBehavior (bool pastelsadd, bool saturatedadd);
    void adjusterChanged     (Adjuster* a, double newval);
    void adjusterChanged     (ThresholdAdjuster* a, int newBottom, int newTop);
    void curveChanged        ();
    void autoOpenCurve       ();

    void enabledChanged          ();
    void protectskins_toggled    ();
    void avoidcolorshift_toggled ();
    void pastsattog_toggled      ();
    std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const;
};


#endif
