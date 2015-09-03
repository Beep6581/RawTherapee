/*
 *  This file is part of RawTherapee.
 */
#ifndef _DEHAZ_H_
#define _DEHAZ_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "thresholdadjuster.h"
#include "colorprovider.h"

class Dehaz : public ToolParamBlock, public FoldableToolPanel,  public CurveListener,
    public AdjusterListener
{

protected:
    CurveEditorGroup* curveEditorGD;
    Adjuster* str;
    Adjuster* scal;
    Adjuster* neigh;
    Adjuster* gain;
    Adjuster* offs;
    Adjuster* vart;
    Adjuster* limd;

    Gtk::Label* labmdh;
    Gtk::HBox* dhbox;
    MyComboBoxText*   dehazmet;
    Gtk::CheckButton* retinex;
    Gtk::Frame* dehazFrame;
    Gtk::CheckButton* medianmap;

    DiagonalCurveEditor* cdshape;
    CurveEditorGroup* transmissionCurveEditorG;
    sigc::connection dehazmetConn;
    FlatCurveEditor* transmissionShape;
    bool lastretinex, lastmedianmap;
    sigc::connection retinexConn, medianmapConn;

public:
    Dehaz   ();
    ~Dehaz  ();

    void read                  (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write                 (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);
    void setBatchMode          (bool batchMode);
    void setDefaults           (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = NULL);
    void trimValues            (rtengine::procparams::ProcParams* pp);
    void adjusterChanged       (Adjuster* a, double newval);
//    void setAdjusterBehavior   (bool splitAdd, bool satThresholdAdd, bool satOpacityAdd, bool strprotectAdd, bool balanceAdd);
    void autoOpenCurve         ();
    void retinexChanged        ();
    void medianmapChanged        ();

    void enabledChanged        ();
    void curveChanged          (CurveEditor* ce);
    void dehazmetChanged();
    void retinexUpdateUI();

};

#endif
