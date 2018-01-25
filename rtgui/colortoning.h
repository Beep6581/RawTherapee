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
    void adjusterChanged       (ThresholdAdjuster* a, double newBottom, double newTop);
    void setAdjusterBehavior   (bool splitAdd, bool satThresholdAdd, bool satOpacityAdd, bool strprotectAdd, bool balanceAdd);
    void neutral_pressed       ();
    //void neutralCurves_pressed ();
    void autoColorTonChanged   (int bwct, int satthres, int satprot);
    bool CTComp_               ();

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
    
private:
    bool resetPressed(GdkEventButton* event);

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
    Gtk::Button *labgridReset;
    LabGrid *labgrid;
    Gtk::HBox *labgridBox;

    IdleRegister idle_register;
};

#endif
