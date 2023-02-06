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
#include "colorprovider.h"
#include "curvelistener.h"
#include "toolpanel.h"

class CurveEditorGroup;
class DiagonalCurveEditor;
class EditDataProvider;
class FlatCurveEditor;

class LCurve final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public CurveListener,
    public ColorProvider
{

protected:
    CurveEditorGroup* curveEditorG;
//    CurveEditorGroup* curveEditorGD;
    Adjuster* brightness;
    Adjuster* contrast;
    Adjuster* chromaticity;
    DiagonalCurveEditor* lshape;
    DiagonalCurveEditor* ashape;
    DiagonalCurveEditor* bshape;
    DiagonalCurveEditor* ccshape;
    DiagonalCurveEditor* lcshape;
    FlatCurveEditor*   chshape;
    FlatCurveEditor*   lhshape;
    FlatCurveEditor*   hhshape;
    Gtk::Label* labmdh;
    Gtk::Box* dhbox;

    DiagonalCurveEditor* clshape;
    DiagonalCurveEditor* cdshape;

    //%%%%%%%%%%%%%%%%
    Gtk::CheckButton* lcredsk;

    MyComboBoxText* gamutmunselmethod;
    sigc::connection   gamutmunselmethodconn;
    rtengine::ProcEvent Evgamutmunsell;

    Adjuster* rstprotection;
    sigc::connection  bwtconn, lcconn;
    bool lastACVal, lastLCVal;

    //%%%%%%%%%%%%%%%%

public:
    static const Glib::ustring TOOL_NAME;

    LCurve ();
    ~LCurve () override;

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;
    void autoOpenCurve  () override;
    void setEditProvider     (EditDataProvider *provider) override;
    void setAdjusterBehavior (bool bradd, bool contradd, bool satadd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;

    void curveChanged (CurveEditor* ce) override;
    void adjusterChanged (Adjuster* a, double newval) override;
    void lcredsk_toggled();
    void gamutmunselChanged();

    void updateCurveBackgroundHistogram(
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
    );

    void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller) override;

    void enabledChanged() override;
};
