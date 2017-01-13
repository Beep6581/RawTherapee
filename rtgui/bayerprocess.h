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
    Gtk::CheckButton* dcbEnhance;
    //Gtk::VBox *allOptions;
    //Gtk::CheckButton* allEnhance;
    Gtk::VBox *lmmseOptions;
    Adjuster* lmmseIterations;
    Gtk::VBox *pixelShiftOptions;
    Adjuster* pixelShiftMotion;
    MyComboBoxText* pixelShiftMotionCorrection;
    Gtk::CheckButton* pixelShiftShowMotion;
    Gtk::CheckButton* pixelShiftShowMotionMaskOnly;
    Gtk::CheckButton* pixelShiftAutomatic;
    Gtk::CheckButton* pixelShiftNonGreenHorizontal;
    Gtk::CheckButton* pixelShiftNonGreenVertical;
    Gtk::CheckButton* pixelShiftNonGreenCross;
    Gtk::CheckButton* pixelShiftNonGreenCross2;
    Gtk::CheckButton* pixelShiftNonGreenAmaze;
    Gtk::CheckButton* pixelShiftGreen;
    Gtk::CheckButton* pixelShiftBlur;
    Gtk::CheckButton* pixelShiftExp0;
    Gtk::CheckButton* pixelShiftHoleFill;
    Gtk::CheckButton* pixelShiftMedian;
    Gtk::CheckButton* pixelShiftMedian3;
    Adjuster* pixelShiftStddevFactorGreen;
    Adjuster* pixelShiftStddevFactorRed;
    Adjuster* pixelShiftStddevFactorBlue;
    Adjuster* pixelShiftEperIso;
    Adjuster* pixelShiftNreadIso;
    Adjuster* pixelShiftPrnu;
    Adjuster* pixelShiftSigma;
    Adjuster* pixelShiftSum;
    Adjuster* pixelShiftRedBlueWeight;
    bool lastDCBen;
    int oldMethod;
    //bool lastALLen;
    sigc::connection methodconn, imagenumberconn, psmcconn, dcbEnhconn,
                     pixelShiftShowMotionconn, pixelShiftShowMotionMaskOnlyconn, pixelShiftAutomaticconn,
                     pixelShiftNonGreenHorizontalconn, pixelShiftNonGreenVerticalconn, pixelShiftHoleFillconn, pixelShiftMedianconn, pixelShiftMedian3conn, pixelShiftNonGreenCrossconn,
                     pixelShiftNonGreenCross2conn, pixelShiftNonGreenAmazeconn, pixelShiftGreenconn, pixelShiftBlurconn, pixelShiftExp0conn;
public:

    BayerProcess ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setBatchMode   (bool batchMode);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void methodChanged ();
    void psMotionCorrectionChanged ();
    void imageNumberChanged ();
    void adjusterChanged     (Adjuster* a, double newval);
    void dcbEnhanceChanged();
    void pixelShiftShowMotionChanged();
    void pixelShiftShowMotionMaskOnlyChanged();
    void pixelShiftAutomaticChanged();
    void pixelShiftNonGreenHorizontalChanged();
    void pixelShiftNonGreenVerticalChanged();
    void pixelShiftHoleFillChanged();
    void pixelShiftMedianChanged();
    void pixelShiftMedian3Changed();
    void pixelShiftGreenChanged();
    void pixelShiftBlurChanged();
    void pixelShiftExp0Changed();
    void pixelShiftNonGreenCrossChanged();
    void pixelShiftNonGreenCross2Changed();
    void pixelShiftNonGreenAmazeChanged();
};

#endif
