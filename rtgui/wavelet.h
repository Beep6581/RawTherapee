/*
 *  This file is part of RawTherapee.
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
 *
 *  2014 Jacques Desmis <jdesmis@gmail.com>
 */

#pragma once

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "thresholdadjuster.h"
#include "colorprovider.h"
#include "guiutils.h"
#include "options.h"

class Wavelet :
    public ToolParamBlock,
    public ThresholdAdjusterListener,
    public AdjusterListener,
    public CurveListener,
    public ColorProvider,
    public rtengine::WaveletListener,
    public FoldableToolPanel
{
public:
    Wavelet ();
    virtual ~Wavelet ();

    bool wavComputed_ ();
    void adjusterChanged (ThresholdAdjuster* a, double newBottom, double newTop);
    void adjusterChanged (Adjuster* a, double newval);
    void adjusterChanged2 (ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR);
    void autoOpenCurve ();
    void curveChanged (CurveEditor* ce);
    void read (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool threshold2add, bool thresadd, bool chroadd, bool chromaadd, bool contrastadd, bool skinadd, bool reschroadd, bool tmrsadd, bool resconadd, bool resconHadd, bool thradd, bool thrHadd, bool skyadd, bool edgradadd, bool edgvaladd, bool strengthadd, bool gammaadd, bool edgedetectadd, bool edgedetectthradd, bool edgedetectthr2add);
    void setBatchMode (bool batchMode);
    void setDefaults  (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void setEditProvider (EditDataProvider *provider);
    void updateToolState (std::vector<int> &tpOpen);
    void write (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void writeOptions (std::vector<int> &tpOpen);

private:
    void foldAllButMe (GdkEventButton* event, MyExpander *expander);

    virtual void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);
    void BAmethodChanged ();
    void NPmethodChanged ();
    void BackmethodChanged ();
    void CHSLmethodChanged ();
    void CHmethodChanged ();
    void CLmethodChanged ();
    void DirmethodChanged ();
    void EDmethodChanged ();
    void HSmethodChanged ();
    void LmethodChanged ();
    void MedgreinfChanged ();
    void TMmethodChanged ();
    void TilesmethodChanged ();
    void avoidToggled ();
    void cbenabToggled ();
    void contrastMinusPressed ();
    void contrastPlusPressed ();
    void daubcoeffmethodChanged ();
    void enabledChanged ();
    void linkedgToggled ();
    void lipstToggled ();
    void medianToggled ();
    void medianlevToggled ();
    void neutralPressed ();
    void neutral_pressed ();
    void neutralchPressed ();
    void tmrToggled ();
    void updatewavLabel ();
    void wavChanged (double nlevel);

    void HSmethodUpdateUI();
    void CHmethodUpdateUI();
//  void CHSLmethodChangedUI();
    void EDmethodUpdateUI();
    void NPmethodUpdateUI();
    void BAmethodUpdateUI();
    void TMmethodUpdateUI();
//  void BackmethodUpdateUI();
    void CLmethodUpdateUI();
//  void TilesmethodUpdateUI();
//  void daubcoeffmethodUpdateUI();
//  void MedgreinfUpdateUI();
//  void DirmethodUpdateUI();
//  void LmethodUpdateUI();
    void adjusterUpdateUI (Adjuster* a);
    void enabledUpdateUI ();
    void medianlevUpdateUI ();
    void cbenabUpdateUI ();
    void lipstUpdateUI ();

    void enableToggled(MyExpander *expander);

    CurveEditorGroup* const curveEditorG;

    CurveEditorGroup* const CCWcurveEditorG;
    CurveEditorGroup* const curveEditorRES;
    CurveEditorGroup* const curveEditorGAM;
    Gtk::HSeparator* const separatorNeutral;
    Gtk::HSeparator* const separatoredge;

    CurveEditorGroup* const opaCurveEditorG;
    FlatCurveEditor* opacityShapeRG;
    CurveEditorGroup* const opacityCurveEditorG;
    FlatCurveEditor* opacityShapeBY;
    CurveEditorGroup* const opacityCurveEditorW;
    CurveEditorGroup* const opacityCurveEditorWL;
    FlatCurveEditor* opacityShape;
    FlatCurveEditor* opacityShapeWL;
    FlatCurveEditor* hhshape;
    FlatCurveEditor* Chshape;
    DiagonalCurveEditor* clshape;

    FlatCurveEditor* ccshape;
    Gtk::CheckButton* const median;
    Gtk::CheckButton* const medianlev;
    Gtk::CheckButton* const linkedg;
    Gtk::CheckButton* const cbenab;
    Gtk::CheckButton* const lipst;
    Gtk::CheckButton* const avoid;
    Gtk::CheckButton* const tmr;

    Gtk::Button* const neutralchButton;
    Adjuster* correction[9];
    Adjuster* correctionch[9];
    Adjuster* const rescon;
    Adjuster* const resconH;
    Adjuster* const reschro;
    Adjuster* const tmrs;
    Adjuster* const gamma;
    Adjuster* const sup;
    Adjuster* const sky;
    Adjuster* const thres;
    Adjuster* const chroma;
    Adjuster* const chro;
    Adjuster* const contrast;
    Adjuster* const thr;
    Adjuster* const thrH;
    Adjuster* const skinprotect;
    Adjuster* const edgrad;
    Adjuster* const edgval;
    Adjuster* const edgthresh;
    Adjuster* const strength;
    Adjuster* const balance;
    Adjuster* const iter;
    Adjuster* greenlow;
    Adjuster* bluelow;
    Adjuster* greenmed;
    Adjuster* bluemed;
    Adjuster* greenhigh;
    Adjuster* bluehigh;

    ThresholdAdjuster* const hueskin;
    ThresholdAdjuster* const hueskin2;
    ThresholdAdjuster* const hllev;
    ThresholdAdjuster* const bllev;
    ThresholdAdjuster* const pastlev;
    ThresholdAdjuster* const satlev;
    ThresholdAdjuster* const edgcont;
    ThresholdAdjuster* const level0noise;
    ThresholdAdjuster* const level1noise;
    ThresholdAdjuster* const level2noise;
    ThresholdAdjuster* const level3noise;

    Adjuster* const threshold;
    Adjuster* const threshold2;
    Adjuster* const edgedetect;
    Adjuster* const edgedetectthr;
    Adjuster* const edgedetectthr2;
    Adjuster* const edgesensi;
    Adjuster* const edgeampli;
    MyComboBoxText* const Lmethod;
    sigc::connection  Lmethodconn;
    MyComboBoxText* const CHmethod;
    sigc::connection  CHmethodconn;
    MyComboBoxText* const CHSLmethod;
    sigc::connection  CHSLmethodconn;
    MyComboBoxText* const EDmethod;
    sigc::connection  EDmethodconn;
    MyComboBoxText* const BAmethod;
    sigc::connection  BAmethodconn;
    MyComboBoxText* const NPmethod;
    sigc::connection  NPmethodconn;
    MyComboBoxText* const TMmethod;
    sigc::connection  TMmethodconn;
    MyComboBoxText* const HSmethod;
    sigc::connection  HSmethodconn;
    MyComboBoxText* const CLmethod;
    sigc::connection  CLmethodconn;
    MyComboBoxText* const Backmethod;
    sigc::connection  Backmethodconn;
    MyComboBoxText* const Tilesmethod;
    sigc::connection  Tilesmethodconn;
    MyComboBoxText* const daubcoeffmethod;
    sigc::connection  daubcoeffmethodconn;
    MyComboBoxText* const Dirmethod;
    sigc::connection  Dirmethodconn;
    MyComboBoxText* const Medgreinf;
    sigc::connection  MedgreinfConn;
    Gtk::Frame* const chanMixerHLFrame;
    Gtk::Frame* const chanMixerMidFrame;
    Gtk::Frame* const chanMixerShadowsFrame;

    Gtk::Label* const wavLabels;
    Gtk::Label* const labmC;
    Gtk::Label* const labmNP;
    MyExpander* const expchroma;
    MyExpander* const expcontrast;
    MyExpander* const expedge;
    MyExpander* const expfinal;
    MyExpander* const expgamut;
    MyExpander* const expnoise;
    MyExpander* const expresid;
    MyExpander* const expsettings;
    MyExpander* const exptoning;

    Gtk::HBox* const neutrHBox;

    sigc::connection enableChromaConn, enableContrastConn, enableEdgeConn, enableFinalConn;
    sigc::connection enableNoiseConn, enableResidConn, enableToningConn;
    sigc::connection medianConn, avoidConn, tmrConn, medianlevConn, linkedgConn, lipstConn, cbenabConn, neutralconn;
    sigc::connection neutralPressedConn;
    sigc::connection contrastPlusPressedConn;
    sigc::connection contrastMinusPressedConn;
    sigc::connection neutralchPressedConn;

    bool lastmedian, lastmedianlev, lastlinkedg, lastavoid, lastlipst, lasttmr, lastcbenab;
    int nextnlevel;
};
