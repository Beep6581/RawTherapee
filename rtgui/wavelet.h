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

    CurveEditorGroup* CCWcurveEditorG;
    CurveEditorGroup* curveEditorRES;
    CurveEditorGroup* curveEditorGAM;
    Gtk::HSeparator* colorSep;
    Gtk::HSeparator* separator3;
    CurveEditorGroup* opaCurveEditorG;
    FlatCurveEditor* opacityShapeRG;
    CurveEditorGroup* opacityCurveEditorG;
    FlatCurveEditor* opacityShapeBY;
    FlatCurveEditor*   hhshape;
    FlatCurveEditor*   Chshape;
	
    FlatCurveEditor* ccshape;
    Gtk::CheckButton * display;
    Gtk::CheckButton * displaylevel;
    Gtk::CheckButton * displaychro;
    Gtk::CheckButton * displaygam;
    Gtk::CheckButton * displayres;
    Gtk::CheckButton * median;
    Gtk::CheckButton * medianlev;
    Gtk::CheckButton * linkedg;
//    Gtk::CheckButton * edgreinf;
    Gtk::CheckButton * lipst;
    Gtk::CheckButton * avoid;
    Gtk::ToggleButton * tbresid;
    Gtk::ToggleButton * tbcontrast;
    Gtk::ToggleButton * tbgamut;
    Gtk::ToggleButton * tbchroma;
    Gtk::ToggleButton * tbtoning;
    Gtk::ToggleButton * tbnoise;
    Gtk::ToggleButton * tbdisplay;
    Gtk::ToggleButton * tbedge;
    Gtk::CheckButton * cbresid;
	Gtk::Image* igRes;
	Gtk::Button * neutralchButton;
	
    Adjuster* correction[9]; 
    Adjuster* correctionch[9]; 
    Adjuster* rescon; 
    Adjuster* resconH; 
    Adjuster* reschro; 
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
	
    Adjuster* threshold;
    Adjuster* threshold2;
    Adjuster* edgedetect;
    Adjuster* edgedetectthr;
    Adjuster* edgedetectthr2;
    MyComboBoxText*   Lmethod;
    sigc::connection  Lmethodconn;
    MyComboBoxText*   CHmethod;
    sigc::connection  CHmethodconn;
    MyComboBoxText*   CHSLmethod;
    sigc::connection  CHSLmethodconn;
    MyComboBoxText*   EDmethod;
    sigc::connection  EDmethodconn;
    MyComboBoxText*   HSmethod;
    sigc::connection  HSmethodconn;
    MyComboBoxText*   CLmethod;
    sigc::connection  CLmethodconn;
    MyComboBoxText*   Backmethod;
    sigc::connection  Backmethodconn;
    MyComboBoxText*   Tilesmethod;
    sigc::connection  Tilesmethodconn;
    MyComboBoxText*   choicemethod;
    sigc::connection  choicemethodconn;
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
    Gtk::Label* colLabel;
    Gtk::Label* interLabel;
	Gtk::Label* wavLabels;
	Gtk::Label* wavLabelsch;
    Gtk::Label* hsmethodLabel;
	Gtk::Label* choiceLabel;
	Gtk::Label* labmC;
	Gtk::Label* labmch;
	Gtk::Label* labmED;
	Gtk::Label* labmedgr;
	Gtk::Label* labmednois;	
	Gtk::Expander* expcontrast;	
	Gtk::Expander* expresid;
	Gtk::Expander* expgamut;	
	Gtk::Expander* expchroma;	
	Gtk::Expander* exptoning;	
	Gtk::Expander* expdisplay;	
	Gtk::Expander* expnoise;	
	Gtk::Expander* expedge;	
	Gtk::HBox* hbresid;
    Gtk::HBox* tilesizeHBox;
    Gtk::HBox* previewLevelsHBox;
    Gtk::HBox* previewBackHBox;
    Gtk::HBox* previewLDirHBox;
    Gtk::HBox* hsmethodHBox;
    Gtk::HBox* choiceHBox;
	Gtk::HBox* ctboxCH;
	Gtk::HBox* ctboxED;
	Gtk::HBox* ctboxch;
	Gtk::HBox* edbox;
	Gtk::HBox* ednoisbox;
	Gtk::HBox* eddebox;
    Gtk::VBox* settingsVBox;
    Gtk::VBox* contrastSHVBox;
    Gtk::Label* tilesizeLabel;
    Gtk::Label* previewLevelsLabel;
    Gtk::Label* previewBackLabel;
	
    sigc::connection expConn,  medianConn, avoidConn, medianlevConn, linkedgConn, lipstConn;
    sigc::connection neutralPressedConn;
    sigc::connection contrastPlusPressedConn;
    sigc::connection contrastMinusPressedConn;
    sigc::connection neutralchPressedConn;

    bool lastdisplay, lastdisplaygam,lastdisplayres,lastdisplaychro, lastdisplaylevel,lastmedian, lastmedianlev, lastlinkedg, lastavoid, lastlipst;
	int nextnlevel;

public:

    Wavelet ();
    virtual ~Wavelet ();
    void curveChanged 	(CurveEditor* ce);
    void setEditProvider     (EditDataProvider *provider);
    void autoOpenCurve  ();

 //   virtual Gtk::Expander* getExpander     () { return NULL; }
 //   virtual void           setExpanded     (bool expanded) {}
 //   virtual bool           getExpanded     () { return false; }



     /*   Gtk::Expander * getExpander() { return expresid; }
        void setExpanded (bool expanded) { if (expresid) expresid->set_expanded( expanded ); }
        bool getExpanded () { if (expresid) return expresid->get_expanded(); return false; }
	*/
    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);
    void adjusterChanged2 (ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR);
	void setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool threshold2add, bool thresadd, bool chroadd,bool chromaadd, bool contrastadd, bool skinadd, bool reschroadd, bool resconadd, bool resconHadd, bool thradd, bool thrHadd, bool skyadd, bool edgradadd, bool edgvaladd, bool strengthadd, bool edgedetectadd, bool edgedetectthradd, bool edgedetectthr2add);
   
    void adjusterChanged (Adjuster* a, double newval);
	void adjusterChanged       (ThresholdAdjuster* a, double newBottom, double newTop);
	
    void enabledChanged ();
    void medianToggled ();
    void medianlevToggled ();
    void linkedgToggled ();
    void lipstToggled ();
    void expcontrastTog ();
    void expresidTog ();
    void expdisplayTog ();
    void expgamutTog ();
    void expchromaTog ();
    void exptoningTog ();
    void avoidToggled ();
    void neutralPressed ();
    void neutralchPressed ();
    void contrastPlusPressed ();
    void contrastMinusPressed ();
    void LmethodChanged      ();
    void choicemethodChanged      ();
    void CHmethodChanged      ();
    void MedgreinfChanged      ();
    void CHSLmethodChanged      ();
    void EDmethodChanged      ();
    void HSmethodChanged      ();
    void CLmethodChanged      ();
    void BackmethodChanged      ();
    void TilesmethodChanged      ();
    void DirmethodChanged      ();
    virtual void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);
    void wavChanged (double nlevel);
    bool wavComputed_ ();
	void updatewavLabel      ();
	
};

#endif
