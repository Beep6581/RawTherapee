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
#ifndef _COLORAPPEARANCE_H_
#define _COLORAPPEARANCE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "mycurve.h"
#include "guiutils.h"
#include "colorprovider.h"

class ColorAppearance : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel, public CurveListener,  public ColorProvider {

  protected:
    Glib::RefPtr<Gtk::Tooltip> bgTTips;
    Glib::RefPtr<Gtk::Tooltip> srTTips;
    Glib::RefPtr<Gdk::Pixbuf> bgPixbuf;
    Glib::RefPtr<Gdk::Pixbuf> srPixbuf;

    Adjuster* degree;
    Adjuster* adapscen;
    Adjuster* adaplum;
    Adjuster* jlight;
    Adjuster* qbright;
    Adjuster* chroma;
    Adjuster* schroma;
    Adjuster* mchroma;
    Adjuster* rstprotection;
    Adjuster* contrast;
    Adjuster* qcontrast;
    Adjuster* colorh;
    MyComboBoxText* toneCurveMode;
    MyComboBoxText* toneCurveMode2;
    MyComboBoxText* toneCurveMode3;

    //Adjuster* edge;
    Gtk::CheckButton* surrsource;
    Gtk::CheckButton* gamut;

    Gtk::CheckButton* enabled;
    MyComboBoxText*   surround;
    sigc::connection  surroundconn;
    MyComboBoxText*   wbmodel;
    sigc::connection  wbmodelconn;
    MyComboBoxText*   algo;
    sigc::connection  algoconn;
    sigc::connection  surrconn;
    sigc::connection  gamutconn;
    sigc::connection  tcmodeconn, tcmode2conn, tcmode3conn;
    CurveEditorGroup* curveEditorG;
    CurveEditorGroup* curveEditorG2;
    CurveEditorGroup* curveEditorG3;

    DiagonalCurveEditor* shape;
    DiagonalCurveEditor* shape2;
    DiagonalCurveEditor* shape3;

    bool lastEnabled;
    bool lastAutoDegree;
    sigc::connection enaConn;
    bool lastsurr;
    bool lastgamut;
    bool bgTTipQuery(int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip);
    bool srTTipQuery(int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip);

  public:

    ColorAppearance ();
    ~ColorAppearance ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged     (Adjuster* a, double newval);
    void adjusterAutoToggled (Adjuster* a, bool newval);
    void enabledChanged      ();
    void surroundChanged     ();
    void wbmodelChanged      ();
    void algoChanged         ();
    void surrsource_toggled  ();
    void gamut_toggled       ();
    void curveChanged        (CurveEditor* ce);
    void curveMode1Changed   ();
    bool curveMode1Changed_  ();
    void curveMode2Changed   ();
    bool curveMode2Changed_  ();
    void curveMode3Changed   ();
    bool curveMode3Changed_  ();

    void expandCurve         (bool isExpanded);
    bool isCurveExpanded     ();
    void autoOpenCurve       ();

    void setAdjusterBehavior (bool degreeadd, bool adapscenadd, bool adaplumadd, bool jlightadd, bool chromaadd, bool contrastadd, bool rstprotectionadd, bool qbrightadd, bool qcontrastadd, bool schromaadd, bool mchromaadd, bool colorhadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);

    virtual void colorForValue (double valX, double valY, int callerId, ColorCaller *caller);
};

#endif
