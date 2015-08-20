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
#ifndef _LABCURVE_H_
#define _LABCURVE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "colorprovider.h"

class LCurve : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public CurveListener, public ColorProvider
{

protected:
    CurveEditorGroup* curveEditorG;
    Adjuster* brightness;
    Adjuster* contrast;
    Adjuster* chromaticity;
    Adjuster* str;
    Adjuster* scal;
    Adjuster* neigh;
    Adjuster* gain;
    Adjuster* offs;
    Adjuster* vart;
   
    DiagonalCurveEditor* lshape;
    DiagonalCurveEditor* ashape;
    DiagonalCurveEditor* bshape;
    DiagonalCurveEditor* ccshape;
    DiagonalCurveEditor* lcshape;
    FlatCurveEditor*   chshape;
    FlatCurveEditor*   lhshape;
    FlatCurveEditor*   hhshape;
    Gtk::Label* labmdh;
    Gtk::HBox* dhbox;
    MyComboBoxText*   dehazmet;

    DiagonalCurveEditor* clshape;

    //%%%%%%%%%%%%%%%%
    Gtk::CheckButton* avoidcolorshift;
    Gtk::CheckButton* lcredsk;

    Adjuster* rstprotection;
    sigc::connection  bwtconn, acconn, lcconn, dehazmetConn;
    bool lastACVal, lastLCVal;

    //%%%%%%%%%%%%%%%%

public:

    LCurve ();
    ~LCurve ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = NULL);
    void setBatchMode   (bool batchMode);
    void autoOpenCurve  ();
    void setEditProvider     (EditDataProvider *provider);
    void setAdjusterBehavior (bool bradd, bool contradd, bool satadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);

    void curveChanged (CurveEditor* ce);
    void adjusterChanged (Adjuster* a, double newval);
    void avoidcolorshift_toggled ();
    void lcredsk_toggled();

    void updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve,/* LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM,  LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma);

    virtual void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);
    
 private:	
    void dehazmetChanged();


   
};

#endif
