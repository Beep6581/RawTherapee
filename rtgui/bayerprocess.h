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

class BayerProcess : public ToolParamBlock, public AdjusterListener, public CheckBoxListener, public FoldableToolPanel, public rtengine::FrameCountListener
{

protected:

    MyComboBoxText* method;
    Gtk::HBox* borderbox;
    Gtk::HBox *imageNumberBox;
    Adjuster* border;
    MyComboBoxText* imageNumber;
    Adjuster* ccSteps;
    Gtk::VBox *dcbOptions;
    Adjuster* dcbIterations;
    CheckBox* dcbEnhance;
    Gtk::VBox *lmmseOptions;
    Adjuster* lmmseIterations;
    Gtk::Frame *pixelShiftFrame;
    Gtk::VBox *pixelShiftOptions;
    MyComboBoxText* pixelShiftMotionMethod;
    MyComboBoxText* pixelShiftDemosaicMethod;
    CheckBox* pixelShiftShowMotion;
    CheckBox* pixelShiftShowMotionMaskOnly;
    CheckBox* pixelShiftNonGreenCross;
    CheckBox* pixelShiftGreen;
    CheckBox* pixelShiftBlur;
    CheckBox* pixelShiftHoleFill;
    CheckBox* pixelShiftMedian;
    CheckBox* pixelShiftEqualBright;
    CheckBox* pixelShiftEqualBrightChannel;
    Adjuster* pixelShiftSmooth;
    Adjuster* pixelShiftEperIso;
    Adjuster* pixelShiftSigma;
    Gtk::VBox *dualDemosaicOptions;
    Adjuster* dualDemosaicContrast;
    int oldMethod;

    IdleRegister idle_register;

    rtengine::ProcEvent EvDemosaicBorder;
    rtengine::ProcEvent EvDemosaicContrast;
    rtengine::ProcEvent EvDemosaicPixelshiftDemosaicMethod;
public:

    BayerProcess ();

    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setAdjusterBehavior(bool falsecoloradd, bool iteradd, bool dualdemozecontrastadd, bool pssigmaadd, bool pssmoothadd, bool pseperisoadd);
    void trimValues(rtengine::procparams::ProcParams* pp);
    void setBatchMode(bool batchMode);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void methodChanged();
    void imageNumberChanged();
    void adjusterChanged(Adjuster* a, double newval);
    void adjusterAutoToggled (Adjuster* a, bool newval);
    void checkBoxToggled(CheckBox* c, CheckValue newval);
    void pixelShiftMotionMethodChanged();
    void pixelShiftDemosaicMethodChanged();
    void FrameCountChanged(int n, int frameNum);
};

#endif
