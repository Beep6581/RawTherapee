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
#ifndef _RGBCURVES_H_
#define _RGBCURVES_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "colorprovider.h"

class RGBCurves : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public CurveListener, public ColorProvider {

  protected:
	CurveEditorGroup* curveEditorG;
	DiagonalCurveEditor* Rshape;
	DiagonalCurveEditor* Gshape;
	DiagonalCurveEditor* Bshape;

    Gtk::CheckButton* lumamode;
    bool lastLumamode;
    sigc::connection lumamodeConn;

  public:

    RGBCurves ();
    ~RGBCurves ();

    void read            (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
    void write           (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setBatchMode    (bool batchMode);
    void setEditProvider (EditDataProvider *provider);
    void autoOpenCurve   ();

    void curveChanged (CurveEditor* ce);
    void updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve, /*LUTu & histCLurve,  LUTu & histLLCurve,*/ LUTu & histLCAM,  LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma);
    void lumamodeChanged  ();
};

#endif
