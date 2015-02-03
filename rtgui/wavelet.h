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

    CurveEditorGroup* CLVcurveEditorG;
    Gtk::HSeparator* colorSep;
    CurveEditorGroup* opaCurveEditorG;
    FlatCurveEditor* opacityShapeRG;
    CurveEditorGroup* opacityCurveEditorG;
    FlatCurveEditor* opacityShapeBY;
	
    FlatCurveEditor* ccshape;
    Gtk::CheckButton * enabled;
    Gtk::CheckButton * display;
    Gtk::CheckButton * displaylevel;
    Gtk::CheckButton * displaychro;
    Gtk::CheckButton * displaygam;
    Gtk::CheckButton * displayres;
    Gtk::CheckButton * median;
    Gtk::CheckButton * avoid;
    Gtk::ToggleButton * tbresid;
    Gtk::ToggleButton * tbcontrast;
    Gtk::ToggleButton * tbgamut;
    Gtk::ToggleButton * tbchroma;
    Gtk::ToggleButton * tbtoning;
    Gtk::ToggleButton * tbutil;
    Gtk::ToggleButton * tbdisplay;
    Gtk::CheckButton * cbresid;
	Gtk::Image* igRes;
	
    Adjuster* correction[9]; 
    Adjuster* tiles; 
    Adjuster* rescon; 
    Adjuster* resconH; 
    Adjuster* reschro; 
    Adjuster* sup; 
    Adjuster* sky; 
    Adjuster* thres; 
    Adjuster* chroma; 
    Adjuster* chro; 
    Adjuster* unif; 
    Adjuster* thr; 
    Adjuster* thrH; 
    Adjuster* skinprotect;
    ThresholdAdjuster* hueskin;
    ThresholdAdjuster* hueskin2;
    ThresholdAdjuster* hllev;
    ThresholdAdjuster* bllev;
    ThresholdAdjuster* pastlev;
    ThresholdAdjuster* satlev;
    Adjuster* threshold;
    Adjuster* threshold2;
    MyComboBoxText*   Lmethod;
    sigc::connection  Lmethodconn;
    MyComboBoxText*   CHmethod;
    sigc::connection  CHmethodconn;
    MyComboBoxText*   HSmethod;
    sigc::connection  HSmethodconn;
    MyComboBoxText*   CLmethod;
    sigc::connection  CLmethodconn;
    MyComboBoxText*   Tilesmethod;
    sigc::connection  Tilesmethodconn;
    MyComboBoxText*   Dirmethod;
    sigc::connection  Dirmethodconn;
	Gtk::Frame* toningFrame;
	Gtk::Frame* residualFrame;
	Gtk::Frame* utilFrame;
	Gtk::Frame* dispFrame;
	Gtk::Frame* levelFrame;
	Gtk::Frame* chromaFrame;
	Gtk::Frame* controlFrame;
    Gtk::Label* colLabel;
    Gtk::Label* interLabel;
	Gtk::Label*   wavLabels;
	Gtk::Expander* expcontrast;	
	Gtk::Expander* expresid;
	Gtk::Expander* expgamut;	
	Gtk::Expander* expchroma;	
	Gtk::Expander* exptoning;	
	Gtk::Expander* exputil;	
	Gtk::Expander* expdisplay;	
	Gtk::HBox* hbresid;
	
    sigc::connection enaConn, expConn,  medianConn, avoidConn;
    sigc::connection neutralPressedConn;
    sigc::connection contrastPlusPressedConn;
    sigc::connection contrastMinusPressedConn;
    
    bool lastEnabled, lastdisplay, lastdisplaygam,lastdisplayres,lastdisplaychro, lastdisplaylevel,lastmedian, lastavoid;
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
	void setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool threshold2add, bool thresadd, bool chroadd,bool chromaadd, bool unifadd, bool skinadd, bool reschroadd, bool resconadd, bool resconHadd, bool thradd, bool thrHadd, bool skyadd);
   
    void adjusterChanged (Adjuster* a, double newval);
    void enabledToggled ();
    void medianToggled ();
    void expcontrastTog ();
    void expresidTog ();
    void expdisplayTog ();
    void expgamutTog ();
    void exputilTog ();
    void expchromaTog ();
    void exptoningTog ();
    void avoidToggled ();
    void neutralPressed ();
    void contrastPlusPressed ();
    void contrastMinusPressed ();
    void LmethodChanged      ();
    void CHmethodChanged      ();
    void HSmethodChanged      ();
    void CLmethodChanged      ();
    void TilesmethodChanged      ();
    void DirmethodChanged      ();
    virtual void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);
    void wavChanged (double nlevel);
    bool wavComputed_ ();
	void updatewavLabel      ();
	
};

#endif
