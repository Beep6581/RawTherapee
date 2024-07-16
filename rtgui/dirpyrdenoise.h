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
#include "checkbox.h"
#include "colorprovider.h"
#include "curvelistener.h"
#include "guiutils.h"
#include "toolpanel.h"

class CurveEditor;
class CurveEditorGroup;
class FlatCurveEditor;
class EditDataProvider;

class DirPyrDenoise final :
    public ToolParamBlock,
    public AdjusterListener,
    public CheckBoxListener,
    public FoldableToolPanel,
    public rtengine::AutoChromaListener,
    public CurveListener,
    public ColorProvider
{
public:
    static const Glib::ustring TOOL_NAME;

    DirPyrDenoise ();
    ~DirPyrDenoise () override;

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;
    void curveChanged   (CurveEditor* ce) override;
    void setEditProvider     (EditDataProvider *provider) override;
    void autoOpenCurve  () override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void checkBoxToggled(CheckBox* c, CheckValue newval) override;
    void enabledChanged  () override;
    void medianChanged  ();
    void chromaChanged (double autchroma, double autred, double autblue) override;
    bool chromaComputed_ ();
    void noiseChanged (double nresid, double highresid) override;
    bool noiseComputed_ ();
    void noiseTilePrev (int tileX, int tileY, int prevX, int prevY, int sizeT, int sizeP) override;
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
    void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller) override;

    void setAdjusterBehavior (bool lumaadd, bool lumdetadd, bool chromaadd, bool chromaredadd, bool chromablueadd, bool gammaadd, bool passesadd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    Glib::ustring getSettingString ();

private:
    rtengine::ProcEvent EvDPDNGain;
    CurveEditorGroup* NoiscurveEditorG;
    CurveEditorGroup* CCcurveEditorG;
    Adjuster* luma;
    Adjuster* Ldetail;
    Adjuster* chroma;
    Adjuster* redchro;
    Adjuster* bluechro;
    CheckBox* autoGain;
    Adjuster* gamma;
    Adjuster* passes;
    FlatCurveEditor* lshape;
    FlatCurveEditor* ccshape;

    sigc::connection medianConn;
    Gtk::CheckButton* median;
    bool lastmedian;
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
    Gtk::Box* ctbox;
    MyComboBoxText*   methodmed;
    sigc::connection  methodmedconn;
    Gtk::Box* ctboxm;
    MyComboBoxText*   rgbmethod;
    sigc::connection  rgbmethodconn;
    Gtk::Box* ctboxrgb;
    double nextchroma;
    double nextred;
    double nextblue;
    double nextnresid;
    double nexthighresid;
    Gtk::Box* ctboxL;
    Gtk::Box* ctboxC;
    Gtk::Box* ctboxC2;
    int nexttileX;
    int nexttileY;
    int nextprevX;
    int nextprevY;
    int nextsizeT;
    int nextsizeP;

    IdleRegister idle_register;
};
