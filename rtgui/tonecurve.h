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
#include "guiutils.h"
#include "toolpanel.h"

class CurveEditor;
class CurveEditorGroup;
class DiagonalCurveEditor;
class EditDataProvider;

class ToneCurve final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::AutoExpListener,
    public CurveListener
{
private:
    IdleRegister idle_register;

protected:
    // from HLRecovery
    Gtk::CheckButton*   hrenabled;
    MyComboBoxText*     method;
    sigc::connection    methconn;
    sigc::connection    enaconn;
    bool                lasthrEnabled;
    Adjuster* hlbl;
    Adjuster* hlth;

    Gtk::Box* abox;
    Gtk::Box* hlrbox;

    Gtk::ToggleButton* autolevels;
    Gtk::Label* lclip;
    MySpinButton* sclip;
    Gtk::Button* neutral;
    Adjuster* expcomp;
    Adjuster* brightness;
    Adjuster* black;
    Adjuster* hlcompr;
    Adjuster* hlcomprthresh;
    Adjuster* shcompr;
    Adjuster* contrast;
    Adjuster* saturation;
    MyComboBoxText* toneCurveMode;
    MyComboBoxText* toneCurveMode2;
    Gtk::ToggleButton *histmatching;
    bool fromHistMatching;
    Gtk::CheckButton *clampOOG;

    bool clipDirty, lastAuto;
    sigc::connection autoconn, neutralconn, tcmodeconn, tcmode2conn;
    sigc::connection histmatchconn;
    CurveEditorGroup* curveEditorG;
    CurveEditorGroup* curveEditorG2;
    DiagonalCurveEditor* shape;
    DiagonalCurveEditor* shape2;

    rtengine::ProcEvent EvHistMatching;
    rtengine::ProcEvent EvHistMatchingBatch;
    rtengine::ProcEvent EvClampOOG;
    rtengine::ProcEvent EvHLbl;
    rtengine::ProcEvent EvHLth;

    // used temporarily in eventing
    double nextExpcomp;
    int nextBrightness;
    int nextContrast;
    int nextBlack;
    int nextHlcompr;
    int nextHlcomprthresh;
    bool nextHLRecons;
    rtengine::procparams::ToneCurveMode nextToneCurveMode;
    std::vector<double> nextToneCurve;

    void setHistmatching(bool enabled);

public:
    static const Glib::ustring TOOL_NAME;

    ToneCurve ();
    ~ToneCurve () override;

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode        (bool batchMode) override;
    void setAdjusterBehavior (bool expadd, bool hlcompadd, bool hlcompthreshadd, bool bradd, bool blackadd, bool shcompadd, bool contradd, bool satadd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void autoOpenCurve       () override;
    void setEditProvider     (EditDataProvider *provider) override;

    float blendPipetteValues (CurveEditor *ce, float chan1, float chan2, float chan3) override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void neutral_pressed ();
    void autolevels_toggled ();
    void clip_changed ();
    bool clip_changed_ ();
    void waitForAutoExp ();
    void enableAll ();
    void curveChanged (CurveEditor* ce) override;
    void curveMode1Changed ();
    bool curveMode1Changed_ ();
    void curveMode2Changed ();
    bool curveMode2Changed_ ();
    void expandCurve (bool isExpanded);
    bool isCurveExpanded ();
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

    void histmatchingToggled();

    void autoExpChanged(double expcomp, int bright, int contr, int black, int hlcompr, int hlcomprthresh, bool hlrecons) override;
    void autoMatchedToneCurveChanged(rtengine::procparams::ToneCurveMode curveMode, const std::vector<double>& curve) override;

    void setRaw (bool raw);

    void hrenabledChanged ();
    void methodChanged ();
    void clampOOGChanged();
};
