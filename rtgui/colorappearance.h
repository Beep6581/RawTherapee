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
#include "colorprovider.h"
#include "curvelistener.h"
#include "guiutils.h"
#include "toolpanel.h"

class DiagonalCurveEditor;
class CurveEditorGroup;
class CurveEditor;

class ColorAppearance final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::AutoCamListener,
    public CurveListener,
    public ColorProvider
{
public:
    ColorAppearance ();
    ~ColorAppearance () override;

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;
    void adjusterChanged     (Adjuster* a, double newval) override;
    void adjusterAutoToggled (Adjuster* a, bool newval) override;
//    void adjusterAdapToggled (Adjuster* a, bool newval);
    void enabledChanged      () override;
    void surroundChanged     ();
    void surrsrcChanged     ();
    void wbmodelChanged      ();
    void algoChanged         ();
    void surrsource_toggled  ();
    void gamut_toggled       ();
//   void badpix_toggled       ();
    void datacie_toggled     ();
    void tonecie_toggled     ();
//    void sharpcie_toggled     ();
    void autoCamChanged (double ccam, double ccamout) override;
    bool autoCamComputed_ ();
    void adapCamChanged (double cadap) override;
    bool adapCamComputed_ ();
    void ybCamChanged (int yb) override;
    bool ybCamComputed_ ();

    void curveChanged        (CurveEditor* ce) override;
    void curveMode1Changed   ();
    bool curveMode1Changed_  ();
    void curveMode2Changed   ();
    bool curveMode2Changed_  ();
    void curveMode3Changed   ();
    bool curveMode3Changed_  ();
    void neutral_pressed       ();

    void expandCurve         (bool isExpanded);
    bool isCurveExpanded     ();
    void autoOpenCurve       () override;

    void setAdjusterBehavior (bool degreeadd, bool adapscenadd, bool adaplumadd, bool badpixsladd, bool jlightadd, bool chromaadd, bool contrastadd, bool rstprotectionadd, bool qbrightadd, bool qcontrastadd, bool schromaadd, bool mchromaadd, bool colorhadd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
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
    void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller) override;
    void updateToolState (std::vector<int> &tpOpen);
    void writeOptions (std::vector<int> &tpOpen);

private:
    bool bgTTipQuery (int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip);
    bool srTTipQuery (int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip);
    void foldAllButMe (GdkEventButton* event, MyExpander *expander);

    Glib::RefPtr<Gtk::Tooltip> bgTTips;
    Glib::RefPtr<Gtk::Tooltip> srTTips;
    Glib::RefPtr<Gdk::Pixbuf> bgPixbuf;
    Glib::RefPtr<Gdk::Pixbuf> srPixbuf;

    Adjuster* degree;
    Adjuster* adapscen;
    Adjuster* ybscen;
    Adjuster* adaplum;
    Adjuster* degreeout;
    Adjuster* badpixsl;
    Adjuster* jlight;
    Adjuster* qbright;
    Adjuster* chroma;
    Adjuster* schroma;
    Adjuster* mchroma;
    Adjuster* rstprotection;
    Adjuster* contrast;
    Adjuster* qcontrast;
    Adjuster* colorh;
    Adjuster* tempout;
    Adjuster* greenout;
    Adjuster* ybout;
    Adjuster* tempsc;
    Adjuster* greensc;

    MyExpander* expadjust;

    MyComboBoxText* toneCurveMode;
    MyComboBoxText* toneCurveMode2;
    MyComboBoxText* toneCurveMode3;

    //Adjuster* edge;
    Gtk::CheckButton* surrsource;
    Gtk::CheckButton* gamut;
//   Gtk::CheckButton* badpix;
    Gtk::CheckButton* datacie;
    Gtk::CheckButton* tonecie;
    //  Gtk::CheckButton* sharpcie;
    Gtk::Button* neutral;
    MyComboBoxText* surrsrc;
    sigc::connection  surrsrcconn;

    MyComboBoxText*   surround;
    sigc::connection  surroundconn;
    MyComboBoxText*   wbmodel;
    sigc::connection  wbmodelconn;
    MyComboBoxText*   algo;
    sigc::connection  algoconn;
    sigc::connection  surrconn;
    sigc::connection  gamutconn, datacieconn, tonecieconn /*,badpixconn , sharpcieconn*/;
    sigc::connection  tcmodeconn, tcmode2conn, tcmode3conn, neutralconn;
    CurveEditorGroup* curveEditorG;
    CurveEditorGroup* curveEditorG2;
    CurveEditorGroup* curveEditorG3;

    DiagonalCurveEditor* shape;
    DiagonalCurveEditor* shape2;
    DiagonalCurveEditor* shape3;
    double nextCcam, nextCcamout, nextCadap;
    int nextYbscn;
    bool lastAutoDegree;
    bool lastAutoAdapscen;
    bool lastAutoDegreeout;
    bool lastAutoybscen;
    bool lastsurr;
    bool lastgamut;
    bool lastdatacie;
    bool lasttonecie;

    IdleRegister idle_register;
};
