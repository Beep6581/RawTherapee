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
#ifndef _BAYERPROCESS_H_
#define _BAYERPROCESS_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "guiutils.h"
#include "toolpanel.h"


class BayerProcess : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel
{

protected:

    MyComboBoxText* method;
    Gtk::HBox *imageNumberBox;
    MyComboBoxText* imageNumber;
    Adjuster* ccSteps;
    Gtk::VBox *dcbOptions;
    Adjuster* dcbIterations;
    MyCheckButton* dcbEnhance;
    Gtk::VBox *lmmseOptions;
    Adjuster* lmmseIterations;
    Gtk::VBox *pixelShiftFrame;
    Gtk::VBox *pixelShiftOptions;
    MyComboBoxText* pixelShiftMotionMethod;
    MyCheckButton* pixelShiftShowMotion;
    MyCheckButton* pixelShiftShowMotionMaskOnly;
    MyCheckButton* pixelShiftNonGreenCross;
    MyCheckButton* pixelShiftGreen;
    MyCheckButton* pixelShiftBlur;
    MyCheckButton* pixelShiftHoleFill;
    MyCheckButton* pixelShiftMedian;
    MyCheckButton* pixelShiftLmmse;
    MyCheckButton* pixelShiftEqualBright;
    Adjuster* pixelShiftSmooth;
    Adjuster* pixelShiftEperIso;
    Adjuster* pixelShiftSigma;
#ifdef PIXELSHIFTDEV
    Adjuster* pixelShiftSum;
    Adjuster* pixelShiftMotion;
    MyComboBoxText* pixelShiftMotionCorrection;
    MyCheckButton* pixelShiftAutomatic;
    MyCheckButton* pixelShiftNonGreenHorizontal;
    MyCheckButton* pixelShiftNonGreenVertical;
    MyCheckButton* pixelShiftNonGreenCross2;
    MyCheckButton* pixelShiftNonGreenAmaze;
    MyCheckButton* pixelShiftExp0;
    MyCheckButton* pixelShiftMedian3;
    Adjuster* pixelShiftStddevFactorGreen;
    Adjuster* pixelShiftStddevFactorRed;
    Adjuster* pixelShiftStddevFactorBlue;
    Adjuster* pixelShiftNreadIso;
    Adjuster* pixelShiftPrnu;
    Adjuster* pixelShiftRedBlueWeight;
#endif
    int oldMethod;
public:

    BayerProcess ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setBatchMode   (bool batchMode);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void methodChanged ();
    void imageNumberChanged ();
    void adjusterChanged     (Adjuster* a, double newval);
    void dcbEnhanceChanged();
    void pixelShiftShowMotionChanged();
    void pixelShiftShowMotionMaskOnlyChanged();
    void pixelShiftHoleFillChanged();
    void pixelShiftMedianChanged();
    void pixelShiftMedian3Changed();
    void pixelShiftGreenChanged();
    void pixelShiftBlurChanged();
    void pixelShiftLmmseChanged();
    void pixelShiftEqualBrightChanged();
    void pixelShiftNonGreenCrossChanged();
    void pixelShiftMotionMethodChanged();
#ifdef PIXELSHIFTDEV
    void psMotionCorrectionChanged ();
    void pixelShiftAutomaticChanged();
    void pixelShiftNonGreenHorizontalChanged();
    void pixelShiftNonGreenVerticalChanged();
    void pixelShiftExp0Changed();
    void pixelShiftNonGreenCross2Changed();
    void pixelShiftNonGreenAmazeChanged();
#endif
};

#endif
