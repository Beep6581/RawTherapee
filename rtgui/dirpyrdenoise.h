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
#ifndef _DIRPYRDENOISE_H_
#define _DIRPYRDENOISE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "colorprovider.h"
#include "guiutils.h"
#include "options.h"

class DirPyrDenoise final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::AutoChromaListener,
    public CurveListener,
    public ColorProvider
{
public:
    DirPyrDenoise ();
    ~DirPyrDenoise ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void setBatchMode   (bool batchMode);
    void curveChanged   (CurveEditor* ce);
    void setEditProvider     (EditDataProvider *provider);
    void autoOpenCurve  ();

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged  ();
    void enhanceChanged  ();
    void medianChanged  ();
    void autochromaChanged  ();
    void chromaChanged (double autchroma, double autred, double autblue);
    bool chromaComputed_ ();
    void noiseChanged (double nresid, double highresid);
    bool noiseComputed_ ();
    void noiseTilePrev (int tileX, int tileY, int prevX, int prevY, int sizeT, int sizeP);
    bool TilePrevComputed_ ();

//    void perform_toggled  ();
    void updateNoiseLabel      ();
    void LmethodChanged      ();
    void CmethodChanged      ();
    void C2methodChanged      ();
    void updateTileLabel      ();
    void updatePrevLabel      ();

    void dmethodChanged      ();
    void medmethodChanged      ();
    void methodmedChanged      ();
    void rgbmethodChanged      ();
    void smethodChanged      ();
    virtual void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);

    void setAdjusterBehavior (bool lumaadd, bool lumdetadd, bool chromaadd, bool chromaredadd, bool chromablueadd, bool gammaadd, bool passesadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    Glib::ustring getSettingString ();

private:
    CurveEditorGroup* NoiscurveEditorG;
    CurveEditorGroup* CCcurveEditorG;
    Adjuster* luma;
    Adjuster* Ldetail;
    Adjuster* chroma;
    Adjuster* redchro;
    Adjuster* bluechro;
    Adjuster* gamma;
    Adjuster* passes;
    FlatCurveEditor* lshape;
    FlatCurveEditor* ccshape;

    Gtk::CheckButton* enhance;
    bool lastenhance;
    sigc::connection enhanConn, medianConn, autochromaConn;
    Gtk::CheckButton* median;
    bool lastmedian;
    Gtk::CheckButton* autochroma;
    bool lastautochroma;
    Gtk::Label*    NoiseLabels;
    Gtk::Label*    TileLabels;
    Gtk::Label*    PrevLabels;

//    Gtk::CheckButton* perform;
//    bool lastperform;
//    sigc::connection perfconn;
    MyComboBoxText*   dmethod;
    sigc::connection  dmethodconn;
    MyComboBoxText*   Lmethod;
    sigc::connection  Lmethodconn;
    MyComboBoxText*   Cmethod;
    sigc::connection  Cmethodconn;
    MyComboBoxText*   C2method;
    sigc::connection  C2methodconn;
    MyComboBoxText*   smethod;
    sigc::connection  smethodconn;
    MyComboBoxText*   medmethod;
    sigc::connection  medmethodconn;
    Gtk::HBox* ctbox;
    MyComboBoxText*   methodmed;
    sigc::connection  methodmedconn;
    Gtk::HBox* ctboxm;
    MyComboBoxText*   rgbmethod;
    sigc::connection  rgbmethodconn;
    Gtk::HBox* ctboxrgb;
    double nextchroma;
    double nextred;
    double nextblue;
    double nextnresid;
    double nexthighresid;
    Gtk::HBox* ctboxL;
    Gtk::HBox* ctboxC;
    Gtk::HBox* ctboxC2;
    int nexttileX;
    int nexttileY;
    int nextprevX;
    int nextprevY;
    int nextsizeT;
    int nextsizeP;

    IdleRegister idle_register;
};

#endif
