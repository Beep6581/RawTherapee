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

#include "colorprovider.h"
#include "curvelistener.h"
#include "toolpanel.h"

class CurveEditorGroup;
class DiagonalCurveEditor;

class RGBCurves final :
    public ToolParamBlock,
    public FoldableToolPanel,
    public CurveListener,
    public ColorProvider
{

protected:
    CurveEditorGroup* curveEditorG;
    DiagonalCurveEditor* Rshape;
    DiagonalCurveEditor* Gshape;
    DiagonalCurveEditor* Bshape;

    Gtk::CheckButton* lumamode;
    bool lastLumamode;
    sigc::connection lumamodeConn;

public:
    static const Glib::ustring TOOL_NAME;

    RGBCurves ();
    ~RGBCurves () override;

    void read            (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write           (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setBatchMode    (bool batchMode) override;
    void setEditProvider (EditDataProvider *provider) override;
    void autoOpenCurve   () override;

    void curveChanged (CurveEditor* ce) override;
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
    void lumamodeChanged  ();
    void enabledChanged() override;
};
