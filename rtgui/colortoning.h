/*
 *  This file is part of RawTherapee.
 */
#ifndef _COLORTONING_H_
#define _COLORTONING_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "thresholdadjuster.h"
#include "colorprovider.h"
#include "labgrid.h"

class ColorToning final :
    public ToolParamBlock,
    public FoldableToolPanel,
    public rtengine::AutoColorTonListener,
    public CurveListener,
    public ColorProvider,
    public ThresholdAdjusterListener,
    public AdjusterListener
{
public:
    ColorToning ();
    ~ColorToning();
    void read                  (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write                 (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setBatchMode          (bool batchMode);
    void setDefaults           (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void trimValues            (rtengine::procparams::ProcParams* pp);
    void adjusterChanged       (Adjuster* a, double newval);
    void adjusterAutoToggled   (Adjuster* a, bool newval);
    void setAdjusterBehavior   (bool splitAdd, bool satThresholdAdd, bool satOpacityAdd, bool strprotectAdd, bool balanceAdd);
    void neutral_pressed       ();
    //void neutralCurves_pressed ();
    void autoColorTonChanged   (int bwct, int satthres, int satprot);
    bool CTComp_               ();

    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop);
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight);
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop);
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight);
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR);

    void enabledChanged        ();
    void curveChanged          (CurveEditor* ce);
    void autosatChanged        ();
    void autoOpenCurve         ();
    void methodChanged         ();
    void twocolorChanged       (bool changedbymethod);
    void twoColorChangedByGui  ();
    void lumamodeChanged       ();

    void colorForValue         (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);

    void setListener(ToolPanelListener *tpl);

    void setEditProvider(EditDataProvider *provider);
    float blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3);

private:
    void onLabRegionSelectionChanged();
    void labRegionAddPressed();
    void labRegionRemovePressed();
    void labRegionUpPressed();
    void labRegionDownPressed();
    void labRegionShowMaskChanged();
    void labRegionPopulateList();
    void labRegionShow(int idx, bool list_only=false);
    void labRegionGet(int idx);

    //Gtk::HSeparator* satLimiterSep;
    Gtk::HSeparator* colorSep;
    CurveEditorGroup* colorCurveEditorG;
    CurveEditorGroup* opacityCurveEditorG;
    CurveEditorGroup* clCurveEditorG;
    CurveEditorGroup* cl2CurveEditorG;
    FlatCurveEditor* opacityShape;
    FlatCurveEditor* colorShape;
    DiagonalCurveEditor* clshape;
    DiagonalCurveEditor* cl2shape;
    Gtk::HBox* ctbox;
    Gtk::Frame *p1Frame;

    Gtk::VBox* chanMixerBox;
    MyComboBoxText* method;
    sigc::connection methodconn;
    MyComboBoxText* twocolor;
    Adjuster* redlow;
    Adjuster* greenlow;
    Adjuster* bluelow;
    Adjuster* redmed;
    Adjuster* greenmed;
    Adjuster* bluemed;
    Adjuster* redhigh;
    Adjuster* greenhigh;
    Adjuster* bluehigh;
    Adjuster* balance;
    Gtk::CheckButton* autosat;
    ThresholdAdjuster* shadowsColSat;
    ThresholdAdjuster* hlColSat;
    Adjuster* satProtectionThreshold;
    Adjuster* saturatedOpacity;
    Adjuster* strength;
    Gtk::Image* iby;
    Gtk::Image* irg;

    Gtk::Button* neutral;
    Gtk::HBox* neutrHBox;
    int nextbw;
    int nextsatth;
    int nextsatpr;
    Glib::ustring nextbalcolor;
    Glib::ustring balcolor;
    sigc::connection neutralconn, twocconn; //, neutralcurvesconn;
    bool lastautosat;
    sigc::connection autosatConn;

    Gtk::CheckButton* lumamode;
    bool lastLumamode;
    sigc::connection lumamodeConn;

    rtengine::ProcEvent EvColorToningLabGridValue;
    LabGrid *labgrid;

    rtengine::ProcEvent EvLabRegionList;
    rtengine::ProcEvent EvLabRegionAB;
    rtengine::ProcEvent EvLabRegionSaturation;
    rtengine::ProcEvent EvLabRegionLightness;
    rtengine::ProcEvent EvLabRegionHueMask;
    rtengine::ProcEvent EvLabRegionChromaticityMask;
    rtengine::ProcEvent EvLabRegionLightnessMask;
    rtengine::ProcEvent EvLabRegionShowMask;

    Gtk::VBox *labRegionBox;
    Gtk::ListViewText *labRegionList;
    Gtk::Button *labRegionAdd;
    Gtk::Button *labRegionRemove;
    Gtk::Button *labRegionUp;
    Gtk::Button *labRegionDown;
    LabGrid *labRegionAB;
    Adjuster *labRegionSaturation;
    Adjuster *labRegionLightness;
    FlatCurveEditor *labRegionHueMask;
    FlatCurveEditor *labRegionChromaticityMask;
    FlatCurveEditor *labRegionLightnessMask;
    Gtk::CheckButton *labRegionShowMask;
    std::vector<rtengine::ColorToningParams::LabCorrectionRegion> labRegionData;
    int labRegionSelected;
    sigc::connection labRegionSelectionConn;

    IdleRegister idle_register;
};

#endif
