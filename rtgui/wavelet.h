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

#ifndef WAVELET_H_INCLUDED
#define WAVELET_H_INCLUDED

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "thresholdadjuster.h"
#include "colorprovider.h"
#include "guiutils.h"
#include "options.h"

class Wavelet : public ToolParamBlock, public ThresholdAdjusterListener, public AdjusterListener, public CurveListener,  public ColorProvider, public rtengine::WaveletListener, public FoldableToolPanel
{
protected:
    Glib::RefPtr<Gtk::Tooltip> bgTTips;
    Glib::RefPtr<Gtk::Tooltip> srTTips;
    Glib::RefPtr<Gdk::Pixbuf> bgPixbuf;
    Glib::RefPtr<Gdk::Pixbuf> srPixbuf;
    CurveEditorGroup* curveEditorG;

    CurveEditorGroup* CCWcurveEditorG;
    CurveEditorGroup* curveEditorRES;
    CurveEditorGroup* curveEditorGAM;
    Gtk::HSeparator* colorSep;
    Gtk::HSeparator* separator3;
    Gtk::HSeparator* separatorCB;
    Gtk::HSeparator* separatorNeutral;
    Gtk::HSeparator* separatoredge;

    CurveEditorGroup* opaCurveEditorG;
    FlatCurveEditor* opacityShapeRG;
    CurveEditorGroup* opacityCurveEditorG;
    FlatCurveEditor* opacityShapeBY;
    CurveEditorGroup* opacityCurveEditorW;
    CurveEditorGroup* opacityCurveEditorWL;
    FlatCurveEditor* opacityShape;
    FlatCurveEditor* opacityShapeWL;
    FlatCurveEditor*   hhshape;
    FlatCurveEditor*   Chshape;
    DiagonalCurveEditor* clshape;
    Gtk::VBox* chanMixerBox;

    FlatCurveEditor* ccshape;
    Gtk::CheckButton * display;
    Gtk::CheckButton * displaylevel;
    Gtk::CheckButton * displaychro;
    Gtk::CheckButton * displaygam;
    Gtk::CheckButton * displayres;
    Gtk::CheckButton * median;
    Gtk::CheckButton * medianlev;
    Gtk::CheckButton * linkedg;
    Gtk::CheckButton * cbenab;
    Gtk::CheckButton * lipst;
    Gtk::CheckButton * avoid;
    Gtk::CheckButton * tmr;

    Gtk::Button * neutralchButton;
    Adjuster* correction[9];
    Adjuster* correctionch[9];
    Adjuster* rescon;
    Adjuster* resconH;
    Adjuster* reschro;
    Adjuster* tmrs;
    Adjuster* gamma;
    Adjuster* sup;
    Adjuster* sky;
    Adjuster* thres;
    Adjuster* chroma;
    Adjuster* chro;
    Adjuster* contrast;
    Adjuster* thr;
    Adjuster* thrH;
    Adjuster* skinprotect;
    Adjuster* edgrad;
    Adjuster* edgval;
    Adjuster* edgthresh;
    Adjuster* strength;
    Adjuster* balance;
    Adjuster* iter;
    Adjuster* greenlow;
    Adjuster* bluelow;
    Adjuster* greenmed;
    Adjuster* bluemed;
    Adjuster* greenhigh;
    Adjuster* bluehigh;

    ThresholdAdjuster* hueskin;
    ThresholdAdjuster* hueskin2;
    ThresholdAdjuster* hllev;
    ThresholdAdjuster* bllev;
    ThresholdAdjuster* pastlev;
    ThresholdAdjuster* satlev;
    ThresholdAdjuster* edgcont;
    ThresholdAdjuster* level0noise;
    ThresholdAdjuster* level1noise;
    ThresholdAdjuster* level2noise;
    ThresholdAdjuster* level3noise;

    Adjuster* threshold;
    Adjuster* threshold2;
    Adjuster* edgedetect;
    Adjuster* edgedetectthr;
    Adjuster* edgedetectthr2;
    Adjuster* edgesensi;
    Adjuster* edgeampli;
    MyComboBoxText*   Lmethod;
    sigc::connection  Lmethodconn;
    MyComboBoxText*   CHmethod;
    sigc::connection  CHmethodconn;
    MyComboBoxText*   CHSLmethod;
    sigc::connection  CHSLmethodconn;
    MyComboBoxText*   EDmethod;
    sigc::connection  EDmethodconn;
    MyComboBoxText*   BAmethod;
    sigc::connection  BAmethodconn;
    MyComboBoxText*   NPmethod;
    sigc::connection  NPmethodconn;
    MyComboBoxText*   TMmethod;
    sigc::connection  TMmethodconn;
    MyComboBoxText*   HSmethod;
    sigc::connection  HSmethodconn;
    MyComboBoxText*   CLmethod;
    sigc::connection  CLmethodconn;
    MyComboBoxText*   Backmethod;
    sigc::connection  Backmethodconn;
    MyComboBoxText*   Tilesmethod;
    sigc::connection  Tilesmethodconn;
    MyComboBoxText*   daubcoeffmethod;
    sigc::connection  daubcoeffmethodconn;
    MyComboBoxText*   Dirmethod;
    sigc::connection  Dirmethodconn;
    MyComboBoxText*   Medgreinf;
    sigc::connection  MedgreinfConn;
    Gtk::Frame* settingsFrame;
    Gtk::Frame* toningFrame;
    Gtk::Frame* residualFrame;
    Gtk::Frame* dispFrame;
    Gtk::Frame* levelFrame;
    Gtk::Frame* chromaFrame;
    Gtk::Frame* controlFrame;
    Gtk::Frame* edgeFrame;
    Gtk::Frame* noiseFrame;
    Gtk::Frame* contrastSHFrame;
    Gtk::Frame* finalFrame;
    Gtk::Frame *chanMixerHLFrame;
    Gtk::Frame *chanMixerMidFrame;
    Gtk::Frame *chanMixerShadowsFrame;
    Gtk::Frame *dFrame;

    Gtk::Label* colLabel;
    Gtk::Label* interLabel;
    Gtk::Label* wavLabels;
    Gtk::Label* hsmethodLabel;
    Gtk::Label* daubcoeffLabel;
    Gtk::Label* ColorBalanceLabel;
    Gtk::Label* labmC;
    Gtk::Label* labmch;
    Gtk::Label* labmED;
    Gtk::Label* labmTM;
    Gtk::Label* labmBA;
    Gtk::Label* labmNP;
    Gtk::Label* labmedgr;
    Gtk::Label* labmednois;
    MyExpander* expchroma;
    MyExpander* expcontrast;
    MyExpander* expedge;
    MyExpander* expfinal;
    MyExpander* expgamut;
    MyExpander* expnoise;
    MyExpander* expresid;
    MyExpander* expsettings;
    MyExpander* exptoning;
    Gtk::HBox* ctboxCB;
    Gtk::HBox* ctboxCH;
    Gtk::HBox* ctboxED;
    Gtk::HBox* ctboxTM;
    Gtk::HBox* hbresid;
    Gtk::HBox* backgroundHBox;
    Gtk::HBox* daubcoeffHBox;
    Gtk::HBox* hsmethodHBox;
    Gtk::HBox* levdirMainHBox;
    Gtk::HBox* levdirSubHBox;
    Gtk::HBox* tilesizeHBox;

    Gtk::HBox* ctboxFI;
    Gtk::HBox* ctboxNP;
    Gtk::HBox* ctboxch;
    Gtk::HBox* edbox;
    Gtk::HBox* ednoisbox;
    Gtk::HBox* eddebox;
    Gtk::VBox* settingsVBox;
    Gtk::VBox* contrastSHVBox;
    Gtk::Label* tilesizeLabel;
    Gtk::Label* levdirMainLabel;
    Gtk::Label* backgroundLabel;
    Gtk::Button* neutral;
    Gtk::HBox* neutrHBox;

    sigc::connection enableChromaConn, enableContrastConn, enableEdgeConn, enableFinalConn;
    sigc::connection enableNoiseConn, enableResidConn, enableToningConn;
    sigc::connection expConn,  medianConn, avoidConn, tmrConn, medianlevConn, linkedgConn, lipstConn, cbenabConn, neutralconn;
    sigc::connection neutralPressedConn;
    sigc::connection contrastPlusPressedConn;
    sigc::connection contrastMinusPressedConn;
    sigc::connection neutralchPressedConn;

    bool lastdisplay, lastdisplaygam, lastdisplayres, lastdisplaychro, lastdisplaylevel, lastmedian, lastmedianlev, lastlinkedg, lastavoid, lastlipst, lasttmr, lastcbenab;
    int nextnlevel;
    double tr;
    double br;
    double tl;
    double bl;

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
};

#endif
