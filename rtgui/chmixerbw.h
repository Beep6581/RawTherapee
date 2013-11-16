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
#ifndef _CHMIXERBW_H_
#define _CHMIXERBW_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "mycurve.h"
#include "colorprovider.h"

class ChMixerbw : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel, public rtengine::AutoBWListener, public CurveListener, public ColorProvider{

  protected:
	FlatCurveEditor*   vshape;
 	CurveEditorGroup*  curveEditorG;
 	CurveEditorGroup*  curveEditorGBW;
    DiagonalCurveEditor* shape;
 	CurveEditorGroup*  curveEditorGBW2;
    DiagonalCurveEditor* shape2;
    Gtk::ToggleButton* autoch;
    Gtk::HBox* abox;
    Gtk::Button* neutral;

    Adjuster *bwred;
    Adjuster *bwgreen;
    Adjuster *bwblue;
    Adjuster *bwredgam;
    Adjuster *bwgreengam;
    Adjuster *bwbluegam;
    Adjuster *bworan;
    Adjuster *bwyell;
    Adjuster *bwcyan;
    Adjuster *bwmag;
    Adjuster *bwpur;
    MyComboBoxText*   met;
    sigc::connection  metconn;
    MyComboBoxText*   fil;
    sigc::connection  filconn;
    MyComboBoxText*   set;
    sigc::connection  setconn;
	Gtk::Label* rlabel;	
	Gtk::Label* glabel;	
	Gtk::Label* blabel;	
	Gtk::Label* rglabel;	
	Gtk::Label* gglabel;	
	Gtk::Label* bglabel;	
	Gtk::Label* Gamlabel;	
	Gtk::Label* orlabel;	
	Gtk::Label* ylabel;	
	Gtk::Label* clabel;	
	Gtk::Label* mlabel;	
	Gtk::Label* plabel;	
	Gtk::Label* setLabel;	
	Gtk::Label* filLabel;	
	
    Gtk::Image *imgIcon[11];
	Gtk::CheckButton* enabled;
	bool lastEnabled;
	sigc::connection enaconn;
	
	Gtk::CheckButton* enabledcc;
	bool lastEnabledcc, lastAuto;
	sigc::connection enaccconn,tcmodeconn,tcmodeconn2, autoconn, neutralconn;
    MyComboBoxText* toneCurveBW;
    MyComboBoxText* toneCurveBW2;
	
	double nextredbw;
	double nextgreenbw;
	double nextbluebw;
	
	

  public:

    ChMixerbw ();
    ~ChMixerbw ();

    void read            (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write           (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults     (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode    (bool batchMode);
    void autoOpenCurve  ();
	
    void autoch_toggled ();
    void neutral_pressed ();
	
    void adjusterChanged (Adjuster* a, double newval);
    void setAdjusterBehavior (bool bwadd, bool bwgadd, bool bwfadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
	void enabledcc_toggled     ();
	void enabled_toggled     ();
    void metChanged         ();	
    void filChanged         ();
    void setChanged         ();
    virtual void colorForValue (double valX, double valY, int callerId, ColorCaller* caller);
    void BWChanged (double redbw, double greenbw, double bluebw);
	bool BWComputed_ ();
	void curveChanged   (CurveEditor* ce);
    void curveMode1Changed ();
    bool curveMode1Changed_ ();
    void curveMode1Changed2 ();
    bool curveMode1Changed2_ ();
};

#endif
