/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>frame
 *
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
 *  2019-2020 Pierre Cabrera <pierre.cab@gmail.com>
 */
#ifndef _LOCALLABTOOLS_H_
#define _LOCALLABTOOLS_H_

#include "curveeditorgroup.h"
#include "curveeditor.h"
#include "labgrid.h"
#include "thresholdadjuster.h"
#include "toolpanel.h"
#include "adjuster.h"

/* ==== LocallabToolListener ==== */
class LocallabTool;
class LocallabToolListener
{
public:
    LocallabToolListener() {};
    virtual ~LocallabToolListener() {};

    virtual void resetOtherMaskView(LocallabTool* current) = 0;
    virtual void toolRemoved(LocallabTool* current) = 0;
};


/* ==== LocallabTool ==== */
class LocallabTool:
    public ToolPanel,
    public CurveListener,
    public ColorProvider,
    public AdjusterListener
{
protected:
    // LocallabTool mode enumeration
    enum modeType {
        Expert = 0,
        Normal = 1,
        Simple = 2
    };

    // LocallabTool parameters
    bool needMode;
    bool isLocActivated;
    const Glib::ustring *spotNameSource;
    LocallabToolListener* locToolListener;

    // LocallabTool generic widgets
    MyExpander* exp;
    MyComboBoxText* const complexity;

    sigc::connection enaExpConn, complexityConn;

    IdleRegister idle_register;

    Glib::ustring getSpotName() const;

public:
    // Locallab tool constructor/destructor
    LocallabTool(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11 = false, bool needMode = true);
    virtual ~LocallabTool();

    // Getter for Locallab tool expander
    MyExpander* getExpander() override
    {
        return exp;
    }

    // Getter/setter for Locallab tool expanded status
    void setExpanded(bool expanded) override
    {
        exp->set_expanded(expanded);
    }

    bool getExpanded() override
    {
        return exp->get_expanded();
    }

    // Setter for Locallab activation indicator
    void isLocallabActivated(bool cond)
    {
        isLocActivated = cond;
    }

    // Setter for spot name source
    void setSpotNameSource(const Glib::ustring *source)
    {
        spotNameSource = source;
    }

    // Setter for Locallab tool listener
    void setLocallabToolListener(LocallabToolListener* ltl)
    {
        locToolListener = ltl;
    }

    // Management functions to add/remove Locallab tool
    void addLocallabTool(bool raiseEvent);
    void removeLocallabTool(bool raiseEvent);
    bool isLocallabToolAdded();

    // Mask background management function
    void refChanged(const double huer, const double lumar, const double chromar, const float fab);

    // Mask preview functions
    virtual bool isMaskViewActive()
    {
        return false;
    };
    virtual void resetMaskView() {};
    virtual void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) {};

    // Advice tooltips management function
    virtual void updateAdviceTooltips(const bool showTooltips) {};

    /* Notes:
     - callerId #1: Mask CC shape (bottom bar) + Color CC/LC shape (left bar)
     - callerId #2: Mask HH shape (main curve and bottom bar)
     - callerId #3: Color LH/HH shape (main curve)
     - callerId #4: Color CC/LC shape (bottom bar)
    */
    void colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller) override;

    // To be implemented
    virtual void setDefaultExpanderVisibility() {};
    virtual void disableListener();
    virtual void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override {};
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override {};
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override {};
    void adjusterChanged(Adjuster* a, double newval) override {};
    void curveChanged(CurveEditor* ce) override {};

protected:
    // To be implemented
    virtual void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) {}; // Only necessary when using mask

private:
    // Remove button event function
    bool on_remove_change(GdkEventButton* event);

    // Tool expander event function
    void foldThemAll(GdkEventButton* event);

    // Complexity mode event function
    void complexityModeChanged();

    // To be implemented
    virtual void enabledChanged() {};
    virtual void convertParamToNormal() {}; // From Expert mode to Normal mode; Only necessary when using mode
    virtual void convertParamToSimple() {}; // From Normal mode to Simple mode; Only necessary when using mode
    virtual void updateGUIToMode(const modeType new_type) {}; // Only necessary when using mode
};

/* ==== LocallabColor ==== */
class LocallabColor:
    public Gtk::Box,
    public LocallabTool,
    public ThresholdAdjusterListener
{
private:
    // Color & Light specific widgets
    Gtk::Frame* const lumFrame;
    Adjuster* const reparcol;
    Adjuster* const gamc;
    Adjuster* const lightness;
    Adjuster* const contrast;
    Adjuster* const chroma;
    Gtk::CheckButton* const curvactiv;
    Gtk::Frame* const gridFrame;
    LabGrid* const labgrid;
    MyComboBoxText* const gridMethod;
    Adjuster* const strengthgrid;
    Adjuster* const sensi;
    Adjuster* const structcol;
    Adjuster* const blurcolde;
    Adjuster* const softradiuscol;
    MyExpander* const exprecov;
    Gtk::Label* const maskusablec;
    Gtk::Label* const maskunusablec;
    Adjuster* const recothresc;
    Adjuster* const lowthresc;
    Adjuster* const higthresc;
    Adjuster* const decayc;
    Gtk::CheckButton* const invers;
    MyExpander* const expgradcol;
    Adjuster* const strcol;
    Adjuster* const strcolab;
    Adjuster* const strcolh;
    Adjuster* const angcol;
    MyExpander* const expcurvcol;
    Gtk::Label* const labqualcurv;
    MyComboBoxText* const qualitycurveMethod;
    CurveEditorGroup* const llCurveEditorG;
    DiagonalCurveEditor* const llshape;
    DiagonalCurveEditor* const ccshape;
    CurveEditorGroup* const clCurveEditorG;
    DiagonalCurveEditor* const clshape;
    DiagonalCurveEditor* const lcshape;
    CurveEditorGroup* const HCurveEditorG;
    FlatCurveEditor* const LHshape;
    CurveEditorGroup* const H3CurveEditorG;
    FlatCurveEditor* const CHshape;
    CurveEditorGroup* const H2CurveEditorG;
    FlatCurveEditor* const HHshape;
    CurveEditorGroup* const rgbCurveEditorG;
    MyComboBoxText* const toneMethod;
    DiagonalCurveEditor* const rgbshape;
    Gtk::CheckButton* const special;
    MyExpander* const expmaskcol1;
    MyComboBoxText* const merMethod;
    ToolParamBlock* const mask7;
    MyComboBoxText* const mergecolMethod;
    Adjuster* const mercol;
    Adjuster* const opacol;
    Adjuster* const conthrcol;
    Gtk::Frame* const gridmerFrame;
    LabGrid* const labgridmerg;
    Adjuster* const merlucol;
    MyExpander* const expmaskcol;
    Gtk::Frame* const mergecolFrame ;
    MyComboBoxText* const showmaskcolMethod;
    MyComboBoxText* const showmaskcolMethodinv;
    Gtk::CheckButton* const enaColorMask;
    CurveEditorGroup* const maskCurveEditorG;
    FlatCurveEditor* const CCmaskshape;
    FlatCurveEditor* const LLmaskshape;
    FlatCurveEditor* const HHmaskshape;
    Gtk::Frame* const struFrame;
    Adjuster* const strumaskcol;
    Gtk::CheckButton* const toolcol;
    Gtk::Frame* const blurFrame;
    Gtk::CheckButton* const fftColorMask;
    Adjuster* const contcol;
    Adjuster* const blurcol;
    Adjuster* const blendmaskcol;
    Gtk::Frame* const toolcolFrame;
    Gtk::Frame* const toolcolFrame2;
    Adjuster* const radmaskcol;
    Adjuster* const lapmaskcol;
    Adjuster* const chromaskcol;
    Adjuster* const gammaskcol;
    Adjuster* const slomaskcol;
    Adjuster* const shadmaskcol;
    CurveEditorGroup* const maskHCurveEditorG;
    FlatCurveEditor* const HHhmaskshape;
    CurveEditorGroup* const mask2CurveEditorG;
    DiagonalCurveEditor* const Lmaskshape;
    CurveEditorGroup* const mask2CurveEditorGwav;
    FlatCurveEditor* const LLmaskcolshapewav;
    ThresholdAdjuster* const csThresholdcol;

    sigc::connection curvactivConn, gridMethodConn, inversConn, qualitycurveMethodConn, toneMethodConn, specialConn, merMethodConn, mergecolMethodConn, showmaskcolMethodConn, showmaskcolMethodConninv, enaColorMaskConn, toolcolConn, fftColorMaskConn;

public:
    LocallabColor();
    ~LocallabColor();

    void setListener(ToolPanelListener* tpl) override;

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
//    void adjusterChanged3(ThresholdAdjuster* a, double newBottom, double newTop) override {};
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void curvactivChanged();
    void gridMethodChanged();
    void inversChanged();
    void qualitycurveMethodChanged();
    void toneMethodChanged();
    void specialChanged();
    void merMethodChanged();
    void mergecolMethodChanged();
    void showmaskcolMethodChanged();
    void showmaskcolMethodChangedinv();
    void enaColorMaskChanged();
    void toolcolChanged();
    void fftColorMaskChanged();

    void updateColorGUI1();
    void updateColorGUI2();
    void updateColorGUI3();
};

/* ==== LocallabExposure ==== */
class LocallabExposure:
    public Gtk::Box,
    public LocallabTool
{
private:
    // Exposure specific widgets
    MyComboBoxText* const expMethod;
//    Gtk::Frame* const pdeFrame;
    MyExpander* const exppde;
    Adjuster* const laplacexp;
    Adjuster* const reparexp;
    Adjuster* const linear;
    Adjuster* const balanexp;
    Adjuster* const gamm;
    Gtk::Label* const labelexpmethod;
    MyComboBoxText* const exnoiseMethod;
//    Gtk::Frame* const fatFrame;
    MyExpander* const expfat;
    Adjuster* const fatamount;
    Adjuster* const fatdetail;
    Gtk::CheckButton* const norm;
    Adjuster* const fatlevel;
    Adjuster* const fatanchor;
    Adjuster* const gamex;
    Adjuster* const sensiex;
    Adjuster* const structexp;
    Adjuster* const blurexpde;
    MyExpander* const exptoolexp;
    Adjuster* const expcomp;
    Adjuster* const black;
    Adjuster* const hlcompr;
    Adjuster* const hlcomprthresh;
    Adjuster* const shadex;
    Adjuster* const shcompr;
    Adjuster* const expchroma;
    CurveEditorGroup* const curveEditorG;
    DiagonalCurveEditor* shapeexpos;
    MyExpander* const exprecove;
    Gtk::Label* const maskusablee;
    Gtk::Label* const maskunusablee;
    Adjuster* const recothrese;
    Adjuster* const lowthrese;
    Adjuster* const higthrese;
    Adjuster* const decaye;
    
    MyExpander* const expgradexp;
    Adjuster* const strexp;
    Adjuster* const angexp;
    Adjuster* const softradiusexp;
    Gtk::CheckButton* const inversex;
    MyExpander* const expmaskexp;
    MyComboBoxText* const showmaskexpMethod;
    MyComboBoxText* const showmaskexpMethodinv;
    Gtk::CheckButton* const enaExpMask;
    Gtk::CheckButton* const enaExpMaskaft;
    CurveEditorGroup* const maskexpCurveEditorG;
    FlatCurveEditor* const CCmaskexpshape;
    FlatCurveEditor* const LLmaskexpshape;
    FlatCurveEditor* const HHmaskexpshape;
    Adjuster* const blendmaskexp;
    Adjuster* const radmaskexp;
    Adjuster* const lapmaskexp;
    Adjuster* const chromaskexp;
    Adjuster* const gammaskexp;
    Adjuster* const slomaskexp;
    Gtk::Frame* const gradFramemask;
    Adjuster* const strmaskexp;
    Adjuster* const angmaskexp;
    CurveEditorGroup* const mask2expCurveEditorG;
    DiagonalCurveEditor* const Lmaskexpshape;

    sigc::connection expMethodConn, exnoiseMethodConn, inversexConn, normConn, showmaskexpMethodConn, showmaskexpMethodConninv, enaExpMaskConn, enaExpMaskaftConn;

public:
    LocallabExposure();
    ~LocallabExposure();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void expMethodChanged();
    void exnoiseMethodChanged();
    void inversexChanged();
    void normChanged();
    void showmaskexpMethodChanged();
    void showmaskexpMethodChangedinv();
    void enaExpMaskChanged();
    void enaExpMaskaftChanged();

    void updateExposureGUI1();
    void updateExposureGUI2();
    void updateExposureGUI3();
};


/* ==== LocallabjShadow ==== */
class LocallabShadow:
    public Gtk::Box,
    public LocallabTool
{
private:
    // Shadow highlight specific widgets
    MyComboBoxText* const shMethod;
    Adjuster* const reparsh;
    const std::array<Adjuster*, 5> multipliersh;
    Adjuster* const detailSH;
    Adjuster* const tePivot;
    Adjuster* const highlights;
    Adjuster* const h_tonalwidth;
    Adjuster* const shadows;
    Adjuster* const s_tonalwidth;
    Adjuster* const sh_radius;
    Adjuster* const sensihs;
    Adjuster* const blurSHde;
    MyExpander* const exprecovs;
    Gtk::Label* const maskusables;
    Gtk::Label* const maskunusables;
    Adjuster* const recothress;
    Adjuster* const lowthress;
    Adjuster* const higthress;
    Adjuster* const decays;
    Gtk::Frame* const gamFrame;
    Adjuster* const gamSH;
    Adjuster* const sloSH;
    MyExpander* const expgradsh;
    Adjuster* const strSH;
    Adjuster* const angSH;
    Gtk::CheckButton* const inverssh;
    MyExpander* const expmasksh;
    MyComboBoxText* const showmaskSHMethod;
    MyComboBoxText* const showmaskSHMethodinv;
    Gtk::CheckButton* const enaSHMask;
    CurveEditorGroup* const maskSHCurveEditorG;
    FlatCurveEditor* const CCmaskSHshape;
    FlatCurveEditor* const LLmaskSHshape;
    FlatCurveEditor* const HHmaskSHshape;
    Adjuster* const blendmaskSH;
    Adjuster* const radmaskSH;
    Adjuster* const lapmaskSH;
    Adjuster* const chromaskSH;
    Adjuster* const gammaskSH;
    Adjuster* const slomaskSH;
    CurveEditorGroup* const mask2SHCurveEditorG;
    DiagonalCurveEditor* const LmaskSHshape;
    Gtk::Frame* const fatSHFrame;
    Adjuster* const fatamountSH;
    Adjuster* const fatanchorSH;

    rtengine::ProcEvent EvlocallabTePivot;

    sigc::connection shMethodConn, inversshConn, showmaskSHMethodConn, showmaskSHMethodConninv, enaSHMaskConn;

public:
    LocallabShadow();
    ~LocallabShadow();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void shMethodChanged();
    void inversshChanged();
    void showmaskSHMethodChanged();
    void showmaskSHMethodChangedinv();
    void enaSHMaskChanged();

    void updateShadowGUI1();
    void updateShadowGUI2();
};

/* ==== LocallabVibrance ==== */
class LocallabVibrance:
    public Gtk::Box,
    public LocallabTool,
    public ThresholdAdjusterListener,
    public ThresholdCurveProvider
{
private:
    // Vibrance specific widgets
    Adjuster* const saturated;
    Adjuster* const pastels;
    Adjuster* const vibgam;
    Adjuster* const warm;
    ThresholdAdjuster* const psThreshold;
    Gtk::CheckButton* const protectSkins;
    Gtk::CheckButton* const avoidColorShift;
    Gtk::CheckButton* const pastSatTog;
    Adjuster* const sensiv;
    CurveEditorGroup* const curveEditorGG;
    DiagonalCurveEditor* const skinTonesCurve;
    MyExpander* const exprecovv;
    Gtk::Label* const maskusablev;
    Gtk::Label* const maskunusablev;
    Adjuster* const recothresv;
    Adjuster* const lowthresv;
    Adjuster* const higthresv;
    Adjuster* const decayv;
    MyExpander* const expgradvib;
    Adjuster* const strvib;
    Adjuster* const strvibab;
    Adjuster* const strvibh;
    Adjuster* const angvib;
    MyExpander* const expmaskvib;
    MyComboBoxText* const showmaskvibMethod;
    Gtk::CheckButton* const enavibMask;
    CurveEditorGroup* const maskvibCurveEditorG;
    FlatCurveEditor* const CCmaskvibshape;
    FlatCurveEditor* const LLmaskvibshape;
    FlatCurveEditor* const HHmaskvibshape;
    Adjuster* const blendmaskvib;
    Adjuster* const radmaskvib;
    Adjuster* const lapmaskvib;
    Adjuster* const chromaskvib;
    Adjuster* const gammaskvib;
    Adjuster* const slomaskvib;
    CurveEditorGroup* const mask2vibCurveEditorG;
    DiagonalCurveEditor* const Lmaskvibshape;

    sigc::connection pskinsConn, ashiftConn, pastsattogConn, showmaskvibMethodConn, enavibMaskConn;

public:
    LocallabVibrance();
    ~LocallabVibrance();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
//    void adjusterChanged3(ThresholdAdjuster* a, double newBottom, double newTop) override {};
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override {}; // Not used
    std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void protectskins_toggled();
    void avoidcolorshift_toggled();
    void pastsattog_toggled();
    void showmaskvibMethodChanged();
    void enavibMaskChanged();

    void updateVibranceGUI();
};

/* ==== LocallabSoft ==== */
class LocallabSoft:
    public Gtk::Box,
    public LocallabTool
{
private:
    // Soft light specific widgets
    MyComboBoxText* const softMethod;
    Gtk::Box* const ctboxsoftmethod;
    MyComboBoxText* const showmasksoftMethod;
    Adjuster* const streng;
    Adjuster* const laplace;
    Adjuster* const sensisf;

    sigc::connection softMethodConn, showmasksoftMethodConn;

public:
    LocallabSoft();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;

private:
    void complexityModeChanged();

    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void softMethodChanged();
    void showmasksoftMethodChanged();

    void updateSoftGUI();
};

/* ==== LocallabBlur ==== */
class LocallabBlur:
    public Gtk::Box,
    public LocallabTool,
    public ThresholdAdjusterListener
//    public ThresholdCurveProvider
    
{
private:
    // Blur & Noise specific widgets
    MyExpander* const expblnoise;
    MyComboBoxText* const blMethod;
    Gtk::CheckButton* const fftwbl;
    Adjuster* const radius;
    Adjuster* const strength;
    Gtk::Frame* const grainFrame;
    Gtk::Frame* const grainFrame2;
    Adjuster* const isogr;
    Adjuster* const strengr;
    Adjuster* const scalegr;
    Adjuster* const divgr;
    MyComboBoxText* const medMethod;
    Adjuster* const itera;
    Adjuster* const guidbl;
    Adjuster* const strbl;
    Adjuster* const epsbl;
    MyExpander* const expdenoise2;
    Adjuster* const recothres;
    Adjuster* const lowthres;
    Adjuster* const higthres;
    Adjuster* const sensibn;
    MyComboBoxText* const blurMethod;
    Gtk::CheckButton* const invbl;
    MyComboBoxText* const chroMethod;
    Gtk::CheckButton* const activlum;
    MyExpander* const expdenoise;
    MyComboBoxText* const quamethod;
    MyExpander* const expdenoisenl;
    MyExpander* const expdenoiselum;
    MyExpander* const expdenoisech;
    CurveEditorGroup* const LocalcurveEditorwavden;
    FlatCurveEditor* const wavshapeden;
    Gtk::Label* const lCLabels;
    Gtk::Label* const lumLabels;
    Gtk::Label* const lum46Labels;
    Gtk::Label* const chroLabels;
    Gtk::Label* const chro46Labels;
    
    MyExpander* const expdenoise1;
    Gtk::Label* const maskusable;
    Gtk::Label* const maskunusable;
    Gtk::Label* const maskusable2;
    Gtk::Label* const maskunusable2;
    Gtk::Label* const maskusable3;
    Gtk::Label* const maskunusable3;

    Gtk::CheckButton* const usemask;
    Adjuster* const lnoiselow;
    Adjuster* const levelthr;
    Adjuster* const levelthrlow;
    Adjuster* const noiselumf0;
    Adjuster* const noiselumf;
    Adjuster* const noiselumf2;
    Adjuster* const noiselumc;
    Adjuster* const noiselumdetail;
    Adjuster* const noiselequal;
    Adjuster* const noisegam;
    CurveEditorGroup* const LocalcurveEditorwavhue;
    FlatCurveEditor* wavhue;
    Adjuster* const noisechrof;
    Adjuster* const noisechroc;
    Adjuster* const noisechrodetail;
    Gtk::Frame* const detailFrame;
    Adjuster* const detailthr;
    Adjuster* const adjblur;
    MyExpander* const expdenoise3;
    Adjuster* const recothresd;
    Adjuster* const lowthresd;
    Adjuster* const midthresd;
    Adjuster* const midthresdch;
    Adjuster* const higthresd;
    Adjuster* const decayd;
    
    Gtk::CheckButton* const invmaskd;
    Gtk::CheckButton* const invmask;
    Gtk::Frame* const prevFrame;
    Adjuster* const nlstr;
    Adjuster* const nldet;
    Adjuster* const nlpat;
    Adjuster* const nlrad;
    Adjuster* const nlgam;
    Adjuster* const bilateral;
    Adjuster* const sensiden;
    Adjuster* const reparden;
    Gtk::Button* neutral;
    MyExpander* const expmaskbl;
    MyComboBoxText* const showmaskblMethod;
    MyComboBoxText* const showmaskblMethodtyp;
    Gtk::CheckButton* const enablMask;
    CurveEditorGroup* const maskblCurveEditorG;
    FlatCurveEditor* const CCmaskblshape;
    FlatCurveEditor* const LLmaskblshape;
    FlatCurveEditor* const HHmaskblshape;
    Adjuster* const strumaskbl;
    Gtk::CheckButton* const toolbl;
    Gtk::Frame* const toolblFrame;
    Gtk::Frame* const toolblFrame2;
    Adjuster* const blendmaskbl;
    Adjuster* const radmaskbl;
    Adjuster* const lapmaskbl;
    Adjuster* const chromaskbl;
    Adjuster* const gammaskbl;
    Adjuster* const slomaskbl;
    Adjuster* const shadmaskbl;
    Adjuster* const shadmaskblsha;
    CurveEditorGroup* const mask2blCurveEditorG;
    DiagonalCurveEditor* const Lmaskblshape;
    CurveEditorGroup* const mask2blCurveEditorGwav;
    FlatCurveEditor* const LLmaskblshapewav;
    Gtk::Box* const quaHBox;
    ThresholdAdjuster* const csThresholdblur;

    sigc::connection blMethodConn, fftwblConn, invblConn, medMethodConn, blurMethodConn, chroMethodConn, activlumConn, showmaskblMethodConn, showmaskblMethodtypConn, enablMaskConn, toolblConn;
    sigc::connection  quamethodconn, usemaskConn, invmaskdConn, invmaskConn, neutralconn;
public:
    LocallabBlur();
    ~LocallabBlur();
    void updatedenlc(const double highres, const double nres, const double highres46, const double nres46, const double Lhighres, const double Lnres, const double Lhighres46, const double Lnres46);


    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;
    void neutral_pressed();

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
//    void adjusterChanged3(ThresholdAdjuster* a, double newBotto, double newTo) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void blMethodChanged();
    void fftwblChanged();
    void usemaskChanged();
    void invmaskdChanged();
    void invmaskChanged();
    void invblChanged();
    void medMethodChanged();
    void blurMethodChanged();
    void chroMethodChanged();
    void activlumChanged();
    void showmaskblMethodChanged();
    void showmaskblMethodtypChanged();
    void enablMaskChanged();
    void toolblChanged();
    void quamethodChanged();

    void updateBlurGUI();
};

/* ==== LocallabTone ==== */
class LocallabTone:
    public Gtk::Box,
    public LocallabTool
{
private:
    // Tone Mapping specific widgets
    Adjuster* const repartm;
    Adjuster* const amount;
    Adjuster* const stren;
    Gtk::CheckButton* const equiltm;
    Adjuster* const gamma;
    Adjuster* const satur;
    Adjuster* const estop;
    Adjuster* const scaltm;
    Adjuster* const rewei;
    Adjuster* const softradiustm;
    Adjuster* const sensitm;
    MyExpander* const exprecovt;
    Gtk::Label* const maskusablet;
    Gtk::Label* const maskunusablet;
    Adjuster* const recothrest;
    Adjuster* const lowthrest;
    Adjuster* const higthrest;
    Adjuster* const decayt;
    MyExpander* const expmasktm;
    MyComboBoxText* const showmasktmMethod;
    Gtk::CheckButton* const enatmMask;
    Gtk::CheckButton* const enatmMaskaft;
    CurveEditorGroup* const masktmCurveEditorG;
    FlatCurveEditor* const CCmasktmshape;
    FlatCurveEditor* const LLmasktmshape;
    FlatCurveEditor* const HHmasktmshape;
    Adjuster* const blendmasktm;
    Adjuster* const lapmasktm;
    Adjuster* const radmasktm;
    Adjuster* const chromasktm;
    Adjuster* const gammasktm;
    Adjuster* const slomasktm;
    CurveEditorGroup* const mask2tmCurveEditorG;
    DiagonalCurveEditor* const Lmasktmshape;

    sigc::connection equiltmConn, showmasktmMethodConn, enatmMaskConn, enatmMaskaftConn;

public:
    LocallabTone();
    ~LocallabTone();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void equiltmChanged();
    void showmasktmMethodChanged();
    void enatmMaskChanged();
    void enatmMaskaftChanged();
};

/* ==== LocallabRetinex ==== */
class LocallabRetinex:
    public Gtk::Box,
    public LocallabTool
{
private:
    // Retinex specific widgets
    Gtk::Frame* const dehaFrame;
    Adjuster* const dehaz;
    Adjuster* const depth;
    Adjuster* const dehazeSaturation;
    Gtk::Frame* const retiFrame;
    Adjuster* const str;
    Gtk::CheckButton* const loglin;
    Adjuster* const sensih;
    Gtk::Frame* const retitoolFrame;
    MyComboBoxText* const retinexMethod;
    Gtk::CheckButton* const fftwreti;
    Gtk::CheckButton* const equilret;
    Adjuster* const neigh;
    Adjuster* const vart;
    Adjuster* const scalereti;
    Adjuster* const limd;
    Adjuster* const offs;
    MyExpander* const expretitools;
    Adjuster* const chrrt;
    Adjuster* const darkness;
    Adjuster* const lightnessreti;
    Adjuster* const cliptm;
    Adjuster* const softradiusret;
    CurveEditorGroup* const LocalcurveEditortransT;
    FlatCurveEditor* const cTtransshape;
    Gtk::Label* const mMLabels;
    Gtk::Label* const transLabels;
    Gtk::Label* const transLabels2;
    CurveEditorGroup* const LocalcurveEditorgainT;
    FlatCurveEditor* const cTgainshape;
    MyExpander* const exprecovr;
    Gtk::Label* const maskusabler;
    Gtk::Label* const maskunusabler;
    Adjuster* const recothresr;
    Adjuster* const lowthresr;
    Adjuster* const higthresr;
    Adjuster* const decayr;
    MyExpander* const expmaskreti;
    MyComboBoxText* const showmaskretiMethod;
    Gtk::CheckButton* const enaretiMask;
    Gtk::CheckButton* const enaretiMasktmap;
    CurveEditorGroup* const maskretiCurveEditorG;
    FlatCurveEditor* const CCmaskretishape;
    FlatCurveEditor* const LLmaskretishape;
    FlatCurveEditor* const HHmaskretishape;
    Adjuster* const blendmaskreti;
    Adjuster* const radmaskreti;
    Adjuster* const lapmaskreti;
    Adjuster* const chromaskreti;
    Adjuster* const gammaskreti;
    Adjuster* const slomaskreti;
    CurveEditorGroup* const mask2retiCurveEditorG;
    DiagonalCurveEditor* const Lmaskretishape;
    Gtk::CheckButton* const inversret;

    sigc::connection loglinConn, retinexMethodConn, fftwretiConn, equilretConn, showmaskretiMethodConn, enaretiMaskConn, enaretiMasktmapConn, inversretConn;

public:
    LocallabRetinex();
    ~LocallabRetinex();

    void updateMinMax(const double cdma, const double cdmin, const double mini, const double maxi, const double Tmean, const double Tsigma, const double Tmin, const double Tmax);

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void loglinChanged();
    void retinexMethodChanged();
    void fftwretiChanged();
    void equilretChanged();
    void showmaskretiMethodChanged();
    void enaretiMaskChanged();
    void enaretiMasktmapChanged();
    void inversretChanged();

    void updateRetinexGUI1();
    void updateRetinexGUI2();
    void updateRetinexGUI3();
};

/* ==== LocallabSharp ==== */
class LocallabSharp:
    public Gtk::Box,
    public LocallabTool
{
private:
    Adjuster* const sharcontrast;
    Adjuster* const sharblur;
    Adjuster* const shargam;
    Adjuster* const sharamount;
    Adjuster* const shardamping;
    Adjuster* const shariter;
    Adjuster* const sharradius;
    Adjuster* const sensisha;
    Gtk::CheckButton* const inverssha;
    Gtk::Frame* const sharFrame;
    MyComboBoxText* const showmasksharMethod;

    sigc::connection inversshaConn, showmasksharMethodConn;

public:
    LocallabSharp();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void inversshaChanged();
    void showmasksharMethodChanged();
};

/* ==== LocallabContrast ==== */
class LocallabContrast:
    public Gtk::Box,
    public LocallabTool,
    public ThresholdAdjusterListener

{
private:
    MyComboBoxText* const localcontMethod;
    Adjuster* const lcradius;
    Adjuster* const lcamount;
    Adjuster* const lcdarkness;
    Adjuster* const lclightness;
    Gtk::Frame* const contFrame;
    Adjuster* const sigmalc;
    CurveEditorGroup* const LocalcurveEditorwav;
    FlatCurveEditor* const wavshape;
    ThresholdAdjuster* const csThreshold;
    Adjuster* const levelwav;
    MyExpander* const expresidpyr;
    Adjuster* const residcont;
    Adjuster* const residchro;
    Adjuster* const residsha;
    Adjuster* const residshathr;
    Adjuster* const residhi;
    Adjuster* const residhithr;
    Adjuster* const gamlc;
    Adjuster* const residgam;
    Adjuster* const residslop;
    Adjuster* const sensilc;
    Adjuster* const reparw;
    Gtk::Frame* const clariFrame;
    Adjuster* const clarilres;
    Adjuster* const claricres;
    Adjuster* const clarisoft;
    Gtk::CheckButton* const origlc;
    MyExpander* const expcontrastpyr;
    Gtk::CheckButton* const wavgradl;
    Adjuster* const sigmalc2;
    Adjuster* const strwav;
    Adjuster* const angwav;
    Gtk::CheckButton* const wavedg;
    Adjuster* const strengthw;
    Adjuster* const sigmaed;
    CurveEditorGroup* const LocalcurveEditorwavedg;
    FlatCurveEditor* const wavshapeedg;
    Adjuster* const gradw;
    Gtk::CheckButton* const waveshow;
    ToolParamBlock* const edgsBoxshow;
    Adjuster* const radiusw;
    Adjuster* const detailw;
    MyComboBoxText* const localedgMethod;
    Adjuster* const tloww;
    Adjuster* const thigw;
    Adjuster* const edgw;
    Adjuster* const basew;
    MyComboBoxText* const localneiMethod;
    Gtk::CheckButton* const wavblur;
    Adjuster* const levelblur;
    Adjuster* const sigmabl;
    Adjuster* const chromablu;
    CurveEditorGroup* const LocalcurveEditorwavlev;
    FlatCurveEditor* const wavshapelev;
    Adjuster* const residblur;
    Gtk::CheckButton* const blurlc;
    MyExpander* const expcontrastpyr2;
    Gtk::CheckButton* const wavcont;
    Adjuster* const sigma;
    Adjuster* const offset;
    Adjuster* const chromalev;
    CurveEditorGroup* const LocalcurveEditorwavcon;
    FlatCurveEditor* const wavshapecon;
    Gtk::CheckButton* const wavcompre;
    CurveEditorGroup* const LocalcurveEditorwavcompre;
    FlatCurveEditor* const wavshapecompre;
    Adjuster* const sigmadr;
    Adjuster* const threswav;
    Adjuster* const residcomp;
    Gtk::CheckButton* const wavcomp;
    Adjuster* const sigmadc;
    Adjuster* const deltad;
    CurveEditorGroup* const LocalcurveEditorwavcomp;
    FlatCurveEditor* const wavshapecomp;
    //Adjuster* const fatres;
    Gtk::CheckButton* const fftwlc;
    MyExpander* const exprecovw;
    Gtk::Label* const maskusablew;
    Gtk::Label* const maskunusablew;
    Adjuster* const recothresw;
    Adjuster* const lowthresw;
    Adjuster* const higthresw;
    Adjuster* const decayw;
    MyExpander* const expmasklc;
    MyComboBoxText* const showmasklcMethod;
    Gtk::CheckButton* const enalcMask;
    CurveEditorGroup* const masklcCurveEditorG;
    FlatCurveEditor* const CCmasklcshape;
    FlatCurveEditor* const LLmasklcshape;
    FlatCurveEditor* const HHmasklcshape;
    Adjuster* const blendmasklc;
    Adjuster* const radmasklc;
    Adjuster* const chromasklc;
    CurveEditorGroup* const mask2lcCurveEditorG;
    DiagonalCurveEditor* const Lmasklcshape;

    sigc::connection localcontMethodConn, origlcConn, wavgradlConn, wavedgConn, localedgMethodConn, waveshowConn, localneiMethodConn, wavblurConn, blurlcConn, wavcontConn, wavcompreConn, wavcompConn, fftwlcConn, showmasklcMethodConn, enalcMaskConn;

public:
    LocallabContrast();
    ~LocallabContrast();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
//    void adjusterChanged3(ThresholdAdjuster* a, double newBottom, double newTop) override {};
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void localcontMethodChanged();
    void origlcChanged();
    void wavgradlChanged();
    void wavedgChanged();
    void localedgMethodChanged();
    void waveshowChanged();
    void localneiMethodChanged();
    void wavblurChanged();
    void blurlcChanged();
    void wavcontChanged();
    void wavcompreChanged();
    void wavcompChanged();
    void fftwlcChanged();
    void showmasklcMethodChanged();
    void enalcMaskChanged();

    void updateContrastGUI1();
    void updateContrastGUI2();
    void updateContrastGUI3();
};

/* ==== LocallabCBDL ==== */
class LocallabCBDL:
    public Gtk::Box,
    public LocallabTool
{
private:
    Gtk::Frame* const levFrame;
    const std::array<Adjuster*, 6> multiplier;
    Adjuster* const chromacbdl;
    Adjuster* const threshold;
    Adjuster* const clarityml;
    Adjuster* const contresid;
    Adjuster* const softradiuscb;
    Adjuster* const sensicb;
    MyExpander* const exprecovcb;
    Gtk::Label* const maskusablecb;
    Gtk::Label* const maskunusablecb;
    Adjuster* const recothrescb;
    Adjuster* const lowthrescb;
    Adjuster* const higthrescb;
    Adjuster* const decaycb;
    MyExpander* const expmaskcb;
    MyComboBoxText* const showmaskcbMethod;
    Gtk::CheckButton* const enacbMask;
    CurveEditorGroup* const maskcbCurveEditorG;
    FlatCurveEditor* const CCmaskcbshape;
    FlatCurveEditor* const LLmaskcbshape;
    FlatCurveEditor* const HHmaskcbshape;
    Adjuster* const blendmaskcb;
    Adjuster* const radmaskcb;
    Adjuster* const lapmaskcb;
    Adjuster* const chromaskcb;
    Adjuster* const gammaskcb;
    Adjuster* const slomaskcb;
    CurveEditorGroup* const mask2cbCurveEditorG;
    DiagonalCurveEditor* const Lmaskcbshape;

    sigc::connection showmaskcbMethodConn, enacbMaskConn;

    Gtk::Button* const lumacontrastMinusButton;
    Gtk::Button* const lumaneutralButton;
    Gtk::Button* const lumacontrastPlusButton;

    sigc::connection lumacontrastMinusPressedConn, lumaneutralPressedConn, lumacontrastPlusPressedConn;

public:
    LocallabCBDL();
    ~LocallabCBDL();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void showmaskcbMethodChanged();
    void enacbMaskChanged();

    void lumacontrastMinusPressed();
    void lumaneutralPressed();
    void lumacontrastPlusPressed();
};

/* ==== LocallabLog ==== */
class LocallabLog:
    public Gtk::Box,
    public LocallabTool
{
private:
    Adjuster* const repar;
    Gtk::CheckButton* const ciecam;
    Gtk::ToggleButton* const autocompute;
    Gtk::Frame* const logPFrame;
    Adjuster* const blackEv;
    Adjuster* const whiteEv;
    Gtk::CheckButton* const fullimage;
    Gtk::Frame* const logFrame;
    Gtk::CheckButton* const Autogray;
    Adjuster* const sourceGray;
    Adjuster* const sourceabs;
    MyComboBoxText*  const sursour;
    Gtk::Box* const surHBox;
    Gtk::Frame* const log1Frame;
    Gtk::Frame* const log2Frame;
    Adjuster* const targetGray;
    Adjuster* const detail;
    Adjuster* const catad;
    Adjuster* const lightl;
    Adjuster* const lightq;
    Adjuster* const contl;
    Adjuster* const contq;
    Adjuster* const contthres;
    Adjuster* const colorfl;
    Adjuster* const saturl;
    Adjuster* const chroml;
    MyExpander* const expL;
    //CurveEditorGroup* const CurveEditorL;
    //DiagonalCurveEditor* const LshapeL;
    Adjuster* const targabs;
    MyComboBoxText*  const surround;
    Gtk::Box* const surrHBox;
    
    Adjuster* const baselog;
    MyExpander* const exprecovl;
    Gtk::Label* const maskusablel;
    Gtk::Label* const maskunusablel;
    Adjuster* const recothresl;
    Adjuster* const lowthresl;
    Adjuster* const higthresl;
    Adjuster* const decayl;
    
    Adjuster* const sensilog;
    Gtk::Frame* const gradlogFrame;
    Adjuster* const strlog;
    Adjuster* const anglog;
    MyExpander* const expmaskL;
    MyComboBoxText* const showmaskLMethod;
    Gtk::CheckButton* const enaLMask;
    CurveEditorGroup* const maskCurveEditorL;
    FlatCurveEditor* const CCmaskshapeL;
    FlatCurveEditor* const LLmaskshapeL;
    FlatCurveEditor* const HHmaskshapeL;
    Adjuster* const blendmaskL;
    Adjuster* const radmaskL;
    Adjuster* const chromaskL;
    CurveEditorGroup* const mask2CurveEditorL;
    DiagonalCurveEditor* const LmaskshapeL;

    sigc::connection autoconn, ciecamconn, fullimageConn, AutograyConn;
    sigc::connection  surroundconn, sursourconn;
    sigc::connection showmaskLMethodConn, enaLMaskConn;
public:
    LocallabLog();
    ~LocallabLog();
    
    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;
    void surroundChanged();
    void sursourChanged();
    void setDefaultExpanderVisibility() override;

    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void curveChanged(CurveEditor* ce) override;

    void updateAutocompute(const float blackev, const float whiteev, const float sourceg, const float sourceab, const float targetg, const float jz1);

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;
    void complexityModeChanged();

    void autocomputeToggled();
    void fullimageChanged();
    void AutograyChanged();
    void ciecamChanged();
    void showmaskLMethodChanged();
    void enaLMaskChanged();
    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void updateLogGUI();
    void updateLogGUI2();
};


/* ==== LocallabMask ==== */
class LocallabMask:
    public Gtk::Box,
    public LocallabTool,
    public ThresholdAdjusterListener
{
private:
    Adjuster* const sensimask;
    Adjuster* const blendmask;
    Adjuster* const blendmaskab;
    Adjuster* const softradiusmask;
    MyComboBoxText* const showmask_Method;
    Gtk::CheckButton* const enamask;
    CurveEditorGroup* const mask_CurveEditorG;
    FlatCurveEditor* const CCmask_shape;
    FlatCurveEditor* const LLmask_shape;
    FlatCurveEditor* const HHmask_shape;
    Gtk::Frame* const struFrame;
    Adjuster* const strumaskmask;
    Gtk::CheckButton* const toolmask;
    Gtk::Frame* const blurFrame;
    Gtk::CheckButton* const fftmask;
    Adjuster* const contmask;
    Adjuster* const blurmask;
    Gtk::Frame* const toolmaskFrame;
    Adjuster* const radmask;
    Adjuster* const lapmask;
    Adjuster* const chromask;
    Adjuster* const gammask;
    Adjuster* const slopmask;
    Adjuster* const shadmask;
    CurveEditorGroup* const mask_HCurveEditorG;
    FlatCurveEditor* const HHhmask_shape;
    CurveEditorGroup* const mask2CurveEditorG;
    DiagonalCurveEditor* const Lmask_shape;
    CurveEditorGroup* const mask2CurveEditorGwav;
    FlatCurveEditor* const LLmask_shapewav;
    ThresholdAdjuster* const csThresholdmask;
    Gtk::Frame* const gradFramemask;
    Adjuster* const str_mask;
    Adjuster* const ang_mask;

    sigc::connection showmask_MethodConn, enamaskConn, toolmaskConn, fftmaskConn;

public:
    LocallabMask();
    ~LocallabMask();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
//    void adjusterChanged3(ThresholdAdjuster* a, double newBottom, double newTop) override {};
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void complexityModeChanged();

    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;

    void showmask_MethodChanged();
    void enamaskChanged();
    void toolmaskChanged();
    void fftmaskChanged();

    void updateMaskGUI();
};


/* ==== Locallabcie ==== */
class Locallabcie:
    public Gtk::Box,
    public ThresholdAdjusterListener,
    public LocallabTool
{
private:
    Adjuster* const sensicie;
    Adjuster* const reparcie;
    Gtk::CheckButton* const jabcie;
    MyComboBoxText*  const modecam;
    MyComboBoxText*  const modecie;
    Gtk::Frame* const jzFrame;
    Gtk::Box* const modeHBoxcam;
    Gtk::Box* const modeHBoxcie;
    Gtk::Frame* const cieFrame;
    Gtk::CheckButton* const Autograycie;
    Adjuster* const sourceGraycie;
    Adjuster* const sourceabscie;
    MyComboBoxText*  const sursourcie;
    Gtk::Box* const surHBoxcie;
    Gtk::Frame* const cie1Frame;
    Gtk::Frame* const cie1lightFrame;
    Gtk::Frame* const cie1contFrame;
    Gtk::Frame* const cie1colorFrame;
    Gtk::Frame* const czlightFrame;
//    Gtk::Frame* const czcontFrame;
    Gtk::Frame* const czcolorFrame;
    Gtk::Frame* const PQFrame;
    Gtk::CheckButton* const qtoj;
    Adjuster* const lightlcie;
    Adjuster* const lightjzcie;
    Adjuster* const contjzcie;
    Adjuster* const adapjzcie;
    Adjuster* const jz100;
    Adjuster* const pqremap;
    Adjuster* const pqremapcam16;
    Gtk::CheckButton* const forcejz;
    MyExpander* const expjz;
    Gtk::Frame* const jzshFrame;
    Adjuster* const hljzcie;
    Adjuster* const hlthjzcie;
    Adjuster* const shjzcie;
    Adjuster* const shthjzcie;
    Adjuster* const radjzcie;
    
    MyExpander* const expwavjz;
    
    Gtk::Frame* const contFramejz;
    Adjuster* const sigmalcjz;
    CurveEditorGroup* const LocalcurveEditorwavjz;
    FlatCurveEditor* const wavshapejz;
    ThresholdAdjuster* const csThresholdjz;
    Gtk::Frame* const clariFramejz;
    Adjuster* const clarilresjz;
    Adjuster* const claricresjz;
    Adjuster* const clarisoftjz;
    
    MyExpander* const expcam16;
    
    Adjuster* const lightqcie;
    Adjuster* const contlcie;
    Adjuster* const contqcie;
    Adjuster* const contthrescie;
    Gtk::Frame* const logjzFrame;
    Gtk::CheckButton* const logjz;
    Adjuster* const blackEvjz;
    Adjuster* const whiteEvjz;
    Adjuster* const targetjz;
    Gtk::Frame* const bevwevFrame;
    Gtk::CheckButton* const forcebw;
    
    Gtk::Frame* const sigmoidFrame;
    Gtk::CheckButton* const sigq;
    Adjuster* const sigmoidldacie;
    Adjuster* const sigmoidthcie;
    Adjuster* const sigmoidblcie;
    Gtk::CheckButton* const sigmoidqjcie;
    Gtk::CheckButton* const logcie;
    Gtk::Frame* const sigmoidjzFrame;
    Gtk::CheckButton* const sigjz;
    Adjuster* const sigmoidldajzcie;
    Adjuster* const sigmoidthjzcie;
    Adjuster* const sigmoidbljzcie;
    
    Adjuster* const colorflcie;
    Adjuster* const saturlcie;
    Adjuster* const rstprotectcie;
    Adjuster* const chromlcie;
    Adjuster* const huecie;
    CurveEditorGroup* const cieCurveEditorG;
    MyComboBoxText* const toneMethodcie;
    DiagonalCurveEditor* const shapecie;
    CurveEditorGroup* const cieCurveEditorG2;
    MyComboBoxText* const toneMethodcie2;
    DiagonalCurveEditor* const shapecie2;
    
    Adjuster* const chromjzcie;
    Adjuster* const saturjzcie;
    Adjuster* const huejzcie;
    CurveEditorGroup* const jz1CurveEditorG;
    DiagonalCurveEditor* const shapejz;
    DiagonalCurveEditor* const shapecz;

    
    Gtk::Frame* const HFramejz;
    Gtk::Frame* const JzHFramejz;
    CurveEditorGroup* const jz2CurveEditorG;
    CurveEditorGroup* const jz3CurveEditorG;
    DiagonalCurveEditor* const shapeczjz;
    FlatCurveEditor* const HHshapejz;
    FlatCurveEditor* const CHshapejz;
    FlatCurveEditor* const LHshapejz;
    Adjuster* const softjzcie;
    Adjuster* const thrhjzcie;
    Gtk::CheckButton* const chjzcie;
    Adjuster* const strsoftjzcie;
   
/*
    Gtk::Frame* const ciezFrame;
    Adjuster* const lightlzcam;
    Adjuster* const lightqzcam;
    Adjuster* const contlzcam;
    Adjuster* const contqzcam;
    Adjuster* const contthreszcam;
    Adjuster* const colorflzcam;
    Adjuster* const saturzcam;
    Adjuster* const chromzcam;
*/
    MyExpander* const expLcie;
    Gtk::Frame* const cie2Frame;
    Adjuster* const targetGraycie;
    Adjuster* const targabscie;
    Adjuster* const detailcie;
    Adjuster* const catadcie;
    MyComboBoxText*  const surroundcie;
    Gtk::Box* const surrHBoxcie;

    MyExpander* const exprecovcie;
    Gtk::Label* const maskusablecie;
    Gtk::Label* const maskunusablecie;
    Adjuster* const recothrescie;
    Adjuster* const lowthrescie;
    Adjuster* const higthrescie;
    Adjuster* const decaycie;

    MyExpander* const expmaskcie;
    MyComboBoxText* const showmaskcieMethod;
    Gtk::CheckButton* const enacieMask;
    CurveEditorGroup* const maskcieCurveEditorG;
    FlatCurveEditor* const CCmaskcieshape;
    FlatCurveEditor* const LLmaskcieshape;
    FlatCurveEditor* const HHmaskcieshape;
    Adjuster* const blendmaskcie;
    Adjuster* const radmaskcie;
    Adjuster* const lapmaskcie;
    Adjuster* const chromaskcie;
    Adjuster* const gammaskcie;
    Adjuster* const slomaskcie;
    
    CurveEditorGroup* const mask2cieCurveEditorG;
    DiagonalCurveEditor* const Lmaskcieshape;
    
    sigc::connection AutograycieConn,  forcejzConn, forcebwConn, qtojConn, showmaskcieMethodConn, enacieMaskConn, jabcieConn, sursourcieconn, surroundcieconn, modecieconn, modecamconn, sigmoidqjcieconn, logcieconn, logjzconn, sigjzconn, sigqconn, chjzcieconn, toneMethodcieConn, toneMethodcieConn2;
public:
    Locallabcie();
    ~Locallabcie();
    
    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;
    void setDefaultExpanderVisibility() override;

    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
//    void adjusterChanged3(ThresholdAdjuster* a, double newBottom, double newTop) override {};
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void sursourcieChanged();
    void surroundcieChanged();
    void modecieChanged();
    void modecamChanged();
    void curveChanged(CurveEditor* ce) override;
    void toneMethodcieChanged();
    void toneMethodcie2Changed();
    void updateAutocompute(const float blackev, const float whiteev, const float sourceg, const float sourceab, const float targetg, const float jz1);

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void convertParamToSimple() override;
    void updateGUIToMode(const modeType new_type) override;
    void complexityModeChanged();
    void AutograycieChanged();
    void forcejzChanged();
    void forcebwChanged();
    void qtojChanged();
    void jabcieChanged();
    void sigmoidqjcieChanged();
    void logcieChanged();
    void logjzChanged();
    void sigjzChanged();
    void sigqChanged();
    void chjzcieChanged();
    void updatecieGUI();
    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz) override;
    void showmaskcieMethodChanged();
    void enacieMaskChanged();

};

#endif
