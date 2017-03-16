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
#include "checkbox.h"
#include "guiutils.h"
#include "toolpanel.h"


class BayerProcess : public ToolParamBlock, public AdjusterListener, public CheckBoxListener, public FoldableToolPanel
{

protected:

    MyComboBoxText* method;
    Gtk::HBox *imageNumberBox;
    MyComboBoxText* imageNumber;
    Adjuster* ccSteps;
    Gtk::VBox *dcbOptions;
    Adjuster* dcbIterations;
    CheckBox* dcbEnhance;
    Gtk::VBox *lmmseOptions;
    Adjuster* lmmseIterations;
    Gtk::VBox *pixelShiftFrame;
    Gtk::VBox *pixelShiftOptions;
    MyComboBoxText* pixelShiftMotionMethod;
    CheckBox* pixelShiftShowMotion;
    CheckBox* pixelShiftShowMotionMaskOnly;
    CheckBox* pixelShiftNonGreenCross;
    CheckBox* pixelShiftGreen;
    CheckBox* pixelShiftBlur;
    CheckBox* pixelShiftHoleFill;
    CheckBox* pixelShiftMedian;
    CheckBox* pixelShiftLmmse;
    CheckBox* pixelShiftEqualBright;
    Adjuster* pixelShiftSmooth;
    Adjuster* pixelShiftEperIso;
    Adjuster* pixelShiftSigma;
#ifdef PIXELSHIFTDEV
    Adjuster* pixelShiftSum;
    Adjuster* pixelShiftMotion;
    MyComboBoxText* pixelShiftMotionCorrection;
    CheckBox* pixelShiftAutomatic;
    CheckBox* pixelShiftNonGreenHorizontal;
    CheckBox* pixelShiftNonGreenVertical;
    CheckBox* pixelShiftNonGreenCross2;
    CheckBox* pixelShiftNonGreenAmaze;
    CheckBox* pixelShiftExp0;
    CheckBox* pixelShiftMedian3;
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
    void adjusterChanged (Adjuster* a, double newval);
    void checkBoxToggled (CheckBox* c, CheckValue newval);
    void pixelShiftMotionMethodChanged();
#ifdef PIXELSHIFTDEV
    void psMotionCorrectionChanged ();
#endif
};

#endif
