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
 *  2019 Pierre Cabrera <pierre.cab@gmail.com>
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
        Normal = 1
    };

    // LocallabTool parameters
    bool needMode;
    bool isLocActivated;
    Glib::ustring spotName;
    LocallabToolListener* locToolListener;

    // LocallabTool generic widgets
    MyExpander* exp;
    MyComboBoxText* const complexity;

    sigc::connection enaExpConn, complexityConn;

    IdleRegister idle_register;

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

    // Setter for spot name
    void setSpotName(const Glib::ustring &spotname)
    {
        spotName = spotname;
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
    void refChanged(const double huer, const double lumar, const double chromar);

    // Mask preview functions
    virtual bool isMaskViewActive()
    {
        return false;
    };
    virtual void resetMaskView() {};
    virtual void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) {};

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
    virtual void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) {}; // Only necessary when using mask

private:
    // Remove button event function
    bool on_remove_change(GdkEventButton* event);

    // Tool expander event function
    void foldThemAll(GdkEventButton* event);

    // Complexity mode event function
    void complexityModeChanged();

    // To be implemented
    virtual void enabledChanged() {};
    virtual void convertParamToNormal() {}; // Only necessary when using mode
    virtual void updateGUIToMode(const modeType new_type) {}; // Only necessary when using mode
};

/* ==== LocallabColor ==== */
class LocallabColor:
    public Gtk::VBox,
    public LocallabTool,
    public ThresholdAdjusterListener
{
private:
    // Color & Light specific widgets
    Gtk::CheckButton* const curvactiv;
    Adjuster* const lightness;
    Adjuster* const contrast;
    Adjuster* const chroma;
    Gtk::Frame* const gridFrame;
    LabGrid* const labgrid;
    MyComboBoxText* const gridMethod;
    Adjuster* const strengthgrid;
    Adjuster* const sensi;
    Adjuster* const structcol;
    Adjuster* const blurcolde;
    Adjuster* const softradiuscol;
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
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

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
    public Gtk::VBox,
    public LocallabTool
{
private:
    // Exposure specific widgets
    MyComboBoxText* const expMethod;
    Gtk::Frame* const pdeFrame;
    Adjuster* const laplacexp;
    Adjuster* const linear;
    Adjuster* const balanexp;
    Adjuster* const gamm;
    MyComboBoxText* const exnoiseMethod;
    Gtk::Frame* const fatFrame;
    Adjuster* const fatamount;
    Adjuster* const fatdetail;
    Adjuster* const fatlevel;
    Adjuster* const fatanchor;
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

    sigc::connection expMethodConn, exnoiseMethodConn, inversexConn, showmaskexpMethodConn, showmaskexpMethodConninv, enaExpMaskConn, enaExpMaskaftConn;

public:
    LocallabExposure();
    ~LocallabExposure();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

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
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

    void expMethodChanged();
    void exnoiseMethodChanged();
    void inversexChanged();
    void showmaskexpMethodChanged();
    void showmaskexpMethodChangedinv();
    void enaExpMaskChanged();
    void enaExpMaskaftChanged();

    void updateExposureGUI1();
    void updateExposureGUI2();
    void updateExposureGUI3();
};

/* ==== LocallabShadow ==== */
class LocallabShadow:
    public Gtk::VBox,
    public LocallabTool
{
private:
    // Shadow highlight specific widgets
    MyComboBoxText* const shMethod;
    const std::array<Adjuster*, 5> multipliersh;
    Adjuster* const detailSH;
    Adjuster* const highlights;
    Adjuster* const h_tonalwidth;
    Adjuster* const shadows;
    Adjuster* const s_tonalwidth;
    Adjuster* const sh_radius;
    Adjuster* const sensihs;
    Adjuster* const blurSHde;
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

    sigc::connection shMethodConn, inversshConn, showmaskSHMethodConn, showmaskSHMethodConninv, enaSHMaskConn;

public:
    LocallabShadow();
    ~LocallabShadow();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

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
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

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
    public Gtk::VBox,
    public LocallabTool,
    public ThresholdAdjusterListener,
    public ThresholdCurveProvider
{
private:
    // Vibrance specific widgets
    Adjuster* const saturated;
    Adjuster* const pastels;
    Adjuster* const warm;
    ThresholdAdjuster* const psThreshold;
    Gtk::CheckButton* const protectSkins;
    Gtk::CheckButton* const avoidColorShift;
    Gtk::CheckButton* const pastSatTog;
    Adjuster* const sensiv;
    CurveEditorGroup* const curveEditorGG;
    DiagonalCurveEditor* const skinTonesCurve;
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
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override {}; // Not used
    std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

    void protectskins_toggled();
    void avoidcolorshift_toggled();
    void pastsattog_toggled();
    void showmaskvibMethodChanged();
    void enavibMaskChanged();

    void updateVibranceGUI();
};

/* ==== LocallabSoft ==== */
class LocallabSoft:
    public Gtk::VBox,
    public LocallabTool
{
private:
    // Soft light specific widgets
    MyComboBoxText* const softMethod;
    Gtk::HBox* const ctboxsoftmethod;
    MyComboBoxText* const showmasksoftMethod;
    Adjuster* const streng;
    Adjuster* const laplace;
    Adjuster* const sensisf;

    sigc::connection softMethodConn, showmasksoftMethodConn;

public:
    LocallabSoft();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

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
    void updateGUIToMode(const modeType new_type) override;

    void softMethodChanged();
    void showmasksoftMethodChanged();

    void updateSoftGUI();
};

/* ==== LocallabBlur ==== */
class LocallabBlur:
    public Gtk::VBox,
    public LocallabTool,
    public ThresholdAdjusterListener
{
private:
    // Blur & Noise specific widgets
    MyExpander* const expblnoise;
    MyComboBoxText* const blMethod;
    Gtk::CheckButton* const fftwbl;
    Adjuster* const radius;
    Adjuster* const strength;
    Gtk::Frame* const grainFrame;
    Adjuster* const isogr;
    Adjuster* const strengr;
    Adjuster* const scalegr;
    MyComboBoxText* const medMethod;
    Adjuster* const itera;
    Adjuster* const guidbl;
    Adjuster* const strbl;
    Adjuster* const epsbl;
    Adjuster* const sensibn;
    MyComboBoxText* const blurMethod;
    MyComboBoxText* const chroMethod;
    Gtk::CheckButton* const activlum;
    MyExpander* const expdenoise;
    CurveEditorGroup* const LocalcurveEditorwavden;
    FlatCurveEditor* const wavshapeden;
    Adjuster* const noiselumf0;
    Adjuster* const noiselumf;
    Adjuster* const noiselumf2;
    Adjuster* const noiselumc;
    Adjuster* const noiselumdetail;
    Adjuster* const noiselequal;
    Adjuster* const noisechrof;
    Adjuster* const noisechroc;
    Adjuster* const noisechrodetail;
    Adjuster* const detailthr;
    Adjuster* const adjblur;
    Adjuster* const bilateral;
    Adjuster* const sensiden;
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
    ThresholdAdjuster* const csThresholdblur;

    sigc::connection blMethodConn, fftwblConn, medMethodConn, blurMethodConn, chroMethodConn, activlumConn, showmaskblMethodConn, showmaskblMethodtypConn, enablMaskConn, toolblConn;

public:
    LocallabBlur();
    ~LocallabBlur();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

    void blMethodChanged();
    void fftwblChanged();
    void medMethodChanged();
    void blurMethodChanged();
    void chroMethodChanged();
    void activlumChanged();
    void showmaskblMethodChanged();
    void showmaskblMethodtypChanged();
    void enablMaskChanged();
    void toolblChanged();

    void updateBlurGUI();
};

/* ==== LocallabTone ==== */
class LocallabTone:
    public Gtk::VBox,
    public LocallabTool
{
private:
    // Tone Mapping specific widgets
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
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

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
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

    void equiltmChanged();
    void showmasktmMethodChanged();
    void enatmMaskChanged();
    void enatmMaskaftChanged();
};

/* ==== LocallabRetinex ==== */
class LocallabRetinex:
    public Gtk::VBox,
    public LocallabTool
{
private:
    // Retinex specific widgets
    Adjuster* const dehaz;
    Adjuster* const depth;
    Gtk::CheckButton* const lumonly;
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

    sigc::connection lumonlyConn, loglinConn, retinexMethodConn, fftwretiConn, equilretConn, showmaskretiMethodConn, enaretiMaskConn, enaretiMasktmapConn, inversretConn;

public:
    LocallabRetinex();
    ~LocallabRetinex();

    void updateMinMax(const double cdma, const double cdmin, const double mini, const double maxi, const double Tmean, const double Tsigma, const double Tmin, const double Tmax);

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

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
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

    void lumonlyChanged();
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
    public Gtk::VBox,
    public LocallabTool
{
private:
    Adjuster* const sharcontrast;
    Adjuster* const sharblur;
    Adjuster* const sharamount;
    Adjuster* const shardamping;
    Adjuster* const shariter;
    Adjuster* const sharradius;
    Adjuster* const sensisha;
    Gtk::CheckButton* const inverssha;
    MyComboBoxText* const showmasksharMethod;

    sigc::connection inversshaConn, showmasksharMethodConn;

public:
    LocallabSharp();

    bool isMaskViewActive() override;
    void resetMaskView() override;
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

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
    void updateGUIToMode(const modeType new_type) override;

    void inversshaChanged();
    void showmasksharMethodChanged();
};

/* ==== LocallabContrast ==== */
class LocallabContrast:
    public Gtk::VBox,
    public LocallabTool,
    public ThresholdAdjusterListener
{
private:
    MyComboBoxText* const localcontMethod;
    Adjuster* const lcradius;
    Adjuster* const lcamount;
    Adjuster* const lcdarkness;
    Adjuster* const lclightness;
    Adjuster* const sigmalc;
    CurveEditorGroup* const LocalcurveEditorwav;
    FlatCurveEditor* const wavshape;
    Adjuster* const levelwav;
    ThresholdAdjuster* const csThreshold;
    MyExpander* const expresidpyr;
    Adjuster* const residcont;
    Adjuster* const residchro;
    Adjuster* const residsha;
    Adjuster* const residshathr;
    Adjuster* const residhi;
    Adjuster* const residhithr;
    Adjuster* const sensilc;
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
    Adjuster* const fatres;
    Gtk::CheckButton* const fftwlc;
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
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

    void updateAdviceTooltips(const bool showTooltips) override;

    void setDefaultExpanderVisibility() override;
    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}; // Not used
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}; // Not used
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
    void curveChanged(CurveEditor* ce) override;

private:
    void enabledChanged() override;
    void convertParamToNormal() override;
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

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
    public Gtk::VBox,
    public LocallabTool
{
private:
    const std::array<Adjuster*, 6> multiplier;
    Adjuster* const chromacbdl;
    Adjuster* const threshold;
    Adjuster* const blurcbdl;
    Adjuster* const clarityml;
    Adjuster* const contresid;
    Adjuster* const softradiuscb;
    Adjuster* const sensicb;
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
    void getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask) override;

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
    void updateGUIToMode(const modeType new_type) override;

    void updateMaskBackground(const double normChromar, const double normLumar, const double normHuer) override;

    void showmaskcbMethodChanged();
    void enacbMaskChanged();

    void lumacontrastMinusPressed();
    void lumaneutralPressed();
    void lumacontrastPlusPressed();
};

/* ==== LocallabLog ==== */
class LocallabLog:
    public Gtk::VBox,
    public LocallabTool
{
private:
    Gtk::ToggleButton* const autocompute;
    Gtk::Frame* const logPFrame;
    Adjuster* const blackEv;
    Adjuster* const whiteEv;
    Gtk::CheckButton* const fullimage;
    Gtk::CheckButton* const Autogray;
    Adjuster* const sourceGray;
    Adjuster* const targetGray;
    Adjuster* const detail;
    Adjuster* const baselog;
    Adjuster* const sensilog;
    Adjuster* const strlog;
    Adjuster* const anglog;

    sigc::connection autoconn, fullimageConn, AutograyConn;

public:
    LocallabLog();
    void updateAdviceTooltips(const bool showTooltips) override;

    void disableListener() override;
    void enableListener() override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster* a, double newval) override;

    void updateAutocompute(const float blackev, const float whiteev, const float sourceg, const float targetg);

private:
    void enabledChanged() override;

    void autocomputeToggled();
    void fullimageChanged();
    void AutograyChanged();

    void updateLogGUI();
};

#endif
