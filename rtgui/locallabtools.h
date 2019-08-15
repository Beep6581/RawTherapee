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

#include "toolpanel.h"
#include "curveeditorgroup.h"
#include "curveeditor.h"
#include "../rtengine/rtengine.h"
#include "labgrid.h"
#include "thresholdadjuster.h"

/* ==== LocallabToolListener ==== */
class LocallabTool;
class LocallabToolListener
{
public:
    LocallabToolListener() {};
    virtual ~LocallabToolListener() {};

    virtual void resetOtherMaskView(LocallabTool* current) {};
};


/* ==== LocallabTool ==== */
class LocallabTool:
    public ToolPanel,
    public CurveListener,
    public ColorProvider,
    public AdjusterListener,
    public rtengine::LocallabListener
{
protected:
    // Enumeration to manage mask type
    enum maskType {
        MaskNone = 0,
        MaskNormal = 1,
        MaskWithTrMap = 2
    };

    // LocallabTool parameters
    const maskType useMask;
    bool isLocActivated;
    Glib::ustring spotName;
    LocallabToolListener* locToolListener;

    // LocallabTool generic widgets
    MyExpander* exp;

    // Mask widgets
    MyExpander* maskExp;
    Gtk::CheckButton* enaMask;
    Gtk::CheckButton* enaMaskTrMap;
    MyComboBoxText* showMaskMethod;
    CurveEditorGroup* maskCurveEditorG;
    FlatCurveEditor* CCMaskShape;
    FlatCurveEditor* LLMaskShape;
    FlatCurveEditor* HHMaskShape;
    Adjuster* blendMask;
    Adjuster* radMask;
    Adjuster* chroMask;
    Adjuster* gamMask;
    Adjuster* sloMask;
    sigc::connection enaExpConn, enaMaskConn, enaMaskTrMapConn, showMaskMethodConn;

private:
    IdleRegister idle_register;

public:
    // Locallab tool constructor/destructor
    LocallabTool(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11 = false, maskType usemask = MaskNone);
    virtual ~LocallabTool();

    // Getter for Locallab tool expander
    MyExpander* getExpander()
    {
        return exp;
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
    void addLocallabTool(bool cond);
    bool isLocallabToolAdded();

    // Mask background management function
    void refChanged(double huer, double lumar, double chromar);

    // Mask preview reset function
    void resetMaskView();

    /* To be completed
    Notes:
     - parent class functions shall be called when using mask
     - callerId #1 and #2 are reserved for mask curve color management
    */
    virtual void colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);

    // To be implemented
    virtual void disableListener();
    virtual void enableListener();
    virtual void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) {};
    virtual void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) {};
    virtual void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) {};
    virtual void adjusterChanged(Adjuster* a, double newval) {};  // At least when using mask
    virtual void curveChanged(CurveEditor* ce) {}; // At least when using mask

protected:
    // To be implemented
    virtual void enaMaskChanged() {}; // Only necessary when using mask
    virtual void enaMaskTrMapChanged() {}; // Only necessary when using mask with transmission map
    virtual void showMaskMethodChanged() {}; // Only necessary when using mask

private:
    // Remove button event function
    bool on_remove_change(GdkEventButton* event);

    // To be implemented
    virtual void enabledChanged() {};
};

/* ==== LocallabColor ==== */
class LocallabColor:
    public Gtk::VBox,
    public LocallabTool
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
    Gtk::Label* const labqualcurv;
    MyComboBoxText* const qualitycurveMethod;
    CurveEditorGroup* const llCurveEditorG;
    DiagonalCurveEditor* llshape;
    DiagonalCurveEditor* ccshape;
    CurveEditorGroup* const HCurveEditorG;
    FlatCurveEditor* LHshape;
    FlatCurveEditor* HHshape;
    Gtk::CheckButton* const invers;
    sigc::connection curvactivConn, gridMethodConn, qualitycurveMethodConn, inversConn;

public:
    LocallabColor();
    ~LocallabColor();

    void setListener(ToolPanelListener* tpl);

    void colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void enaMaskChanged();
    void showMaskMethodChanged();

    void curvactivChanged();
    void inversChanged();
    void qualitycurveMethodChanged();
    void gridMethodChanged();

    void updateColorGUI();
};

/* ==== LocallabExposure ==== */
class LocallabExposure:
    public Gtk::VBox,
    public LocallabTool
{
private:
    MyComboBoxText* const expMethod;
    Gtk::Frame* const pdeFrame;
    Adjuster* const laplacexp;
    Adjuster* const linear;
    Adjuster* const balanexp;
    Adjuster* const expcomp;
    Adjuster* const hlcompr;
    Adjuster* const hlcomprthresh;
    Adjuster* const black;
    Adjuster* const shadex;
    Adjuster* const shcompr;
    Adjuster* const expchroma;
    Adjuster* const warm;
    Adjuster* const sensiex;
    Adjuster* const structexp;
    Adjuster* const blurexpde;
    Adjuster* const softradiusexp;
    CurveEditorGroup* const curveEditorG;
    DiagonalCurveEditor* shapeexpos;
    Gtk::CheckButton* const inversex;
    sigc::connection expMethodConn, inversexConn;

public:
    LocallabExposure();
    ~LocallabExposure();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void enaMaskChanged();
    void showMaskMethodChanged();

    void expMethodChanged();
    void inversexChanged();

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
    Adjuster* const highlights;
    Adjuster* const h_tonalwidth;
    Adjuster* const shadows;
    Adjuster* const s_tonalwidth;
    Adjuster* const sh_radius;
    Adjuster* const sensihs;
    Adjuster* const blurSHde;
    Gtk::CheckButton* const inverssh;
    sigc::connection inversshConn;

public:
    LocallabShadow();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void enaMaskChanged();
    void showMaskMethodChanged();

    void inversshChanged();

    void updateShadowGUI();
};

/* ==== LocallabVibrance ==== */
class LocallabVibrance:
    public Gtk::VBox,
    public LocallabTool,
    public ThresholdAdjusterListener,
    public ThresholdCurveProvider
{
private:
    Adjuster* const saturated;
    Adjuster* const pastels;
    ThresholdAdjuster* const psThreshold;
    Gtk::CheckButton* const protectSkins;
    Gtk::CheckButton* const avoidColorShift;
    Gtk::CheckButton* const pastSatTog;
    Adjuster* const sensiv;
    CurveEditorGroup* const curveEditorGG;
    DiagonalCurveEditor* skinTonesCurve;
    sigc::connection pskinsConn, ashiftConn, pastsattogConn;

public:
    LocallabVibrance();
    ~LocallabVibrance();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) {};
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) {};
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop);
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) {};
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) {};
    std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const;
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void protectskins_toggled();
    void avoidcolorshift_toggled();
    void pastsattog_toggled();

    void updateVibranceGUI();
};

/* ==== LocallabSoft ==== */
class LocallabSoft:
    public Gtk::VBox,
    public LocallabTool
{
private:
    MyComboBoxText* const softMethod;
    Gtk::HBox* const ctboxsoftmethod;
    MyComboBoxText* const showmasksoftMethod;
    Adjuster* const streng;
    Adjuster* const laplace;
    Adjuster* const sensisf;
    sigc::connection softMethodConn, showmasksoftMethodConn;

public:
    LocallabSoft();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);

private:
    void enabledChanged();

    void softMethodChanged();
    void showmasksoftMethodChanged();

    void updateSoftGUI();
};

/* ==== LocallabBlur ==== */
class LocallabBlur:
    public Gtk::VBox,
    public LocallabTool
{
private:
    Adjuster* const radius;
    Adjuster* const strength;
    Adjuster* const sensibn;
    MyComboBoxText* const blurMethod;
    Gtk::CheckButton* const activlum;
    sigc::connection blurMethodConn, activlumConn;

public:
    LocallabBlur();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);

private:
    void enabledChanged();

    void blurMethodChanged();
    void activlumChanged();
};

/* ==== LocallabTone ==== */
class LocallabTone:
    public Gtk::VBox,
    public LocallabTool
{
private:
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
    sigc::connection equiltmConn;

public:
    LocallabTone();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void enaMaskChanged();
    void showMaskMethodChanged();

    void equiltmChanged();
};

/* ==== LocallabRetinex ==== */
class LocallabRetinex:
    public Gtk::VBox,
    public LocallabTool
{
private:
    MyComboBoxText* const retinexMethod;
    Gtk::CheckButton* const fftwreti;
    Gtk::CheckButton* const equilret;
    Adjuster* const str;
    Adjuster* const chrrt;
    Adjuster* const neigh;
    Adjuster* const vart;
    Adjuster* const scalereti;
    Adjuster* const limd;
    Adjuster* const darkness;
    Adjuster* const lightnessreti;
    Adjuster* const dehaz;
    Adjuster* const softradiusret;
    Adjuster* const sensih;
    CurveEditorGroup* const LocalcurveEditorgainT;
    FlatCurveEditor* cTgainshape;
    Gtk::CheckButton* const inversret;
    sigc::connection retinexMethodConn, fftwretiConn, equilretConn, inversretConn;

public:
    LocallabRetinex();
    ~LocallabRetinex();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void enaMaskChanged();
    void enaMaskTrMapChanged();
    void showMaskMethodChanged();

    void retinexMethodChanged();
    void fftwretiChanged();
    void equilretChanged();
    void inversretChanged();

    void updateRetinexGUI();
    void updateRetinexGUI2();
};

/* ==== LocallabSharp ==== */
class LocallabSharp:
    public Gtk::VBox,
    public LocallabTool
{
private:
    Adjuster* const sharcontrast;
    Adjuster* const sharradius;
    Adjuster* const sharamount;
    Adjuster* const shardamping;
    Adjuster* const shariter;
    Adjuster* const sharblur;
    Adjuster* const sensisha;
    Gtk::CheckButton* const inverssha;
    sigc::connection inversshaConn;

public:
    LocallabSharp();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);

private:
    void enabledChanged();

    void inversshaChanged();
};

/* ==== LocallabContrast ==== */
class LocallabContrast:
    public Gtk::VBox,
    public LocallabTool
{
private:
    MyComboBoxText* const localcontMethod;
    Adjuster* const lcradius;
    Adjuster* const lcamount;
    Adjuster* const lcdarkness;
    Adjuster* const lclightness;
    CurveEditorGroup* const LocalcurveEditorwav;
    FlatCurveEditor* wavshape;
    Adjuster* const levelwav;
    Adjuster* const residcont;
    Adjuster* const sensilc;
    Gtk::CheckButton* const fftwlc;
    sigc::connection localcontMethodConn, fftwlcConn;

public:
    LocallabContrast();
    ~LocallabContrast();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void localcontMethodChanged();
    void fftwlcChanged();

    void updateContrastGUI();
};

/* ==== LocallabCBDL ==== */
class LocallabCBDL:
    public Gtk::VBox,
    public LocallabTool
{
private:
    Adjuster* multiplier[6];
    Adjuster* const chromacbdl;
    Adjuster* const threshold;
    Adjuster* const blurcbdl;
    Adjuster* const clarityml;
    Adjuster* const contresid;
    Adjuster* const softradiuscb;
    Adjuster* const sensicb;

    Gtk::Button* const lumacontrastMinusButton;
    Gtk::Button* const lumaneutralButton;
    Gtk::Button* const lumacontrastPlusButton;
    sigc::connection lumacontrastMinusPressedConn, lumaneutralPressedConn, lumacontrastPlusPressedConn;

public:
    LocallabCBDL();

    void disableListener();
    void enableListener();
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);
    void curveChanged(CurveEditor* ce);

private:
    void enabledChanged();

    void enaMaskChanged();
    void showMaskMethodChanged();

    void lumacontrastMinusPressed();
    void lumaneutralPressed();
    void lumacontrastPlusPressed();
};

/* ==== LocallabDenoise ==== */
class LocallabDenoise:
    public Gtk::VBox,
    public LocallabTool
{
private:
    Adjuster* const noiselumf0;
    Adjuster* const noiselumf;
    Adjuster* const noiselumf2;
    Adjuster* const noiselumc;
    Adjuster* const noiselumdetail;
    Adjuster* const noiselequal;
    Adjuster* const noisechrof;
    Adjuster* const noisechroc;
    Adjuster* const noisechrodetail;
    Adjuster* const adjblur;
    Adjuster* const bilateral;
    Adjuster* const sensiden;

public:
    LocallabDenoise();

    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged(Adjuster* a, double newval);

private:
    void enabledChanged();
};

#endif
