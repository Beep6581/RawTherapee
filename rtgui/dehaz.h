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

class Dehaz : public ToolParamBlock, public FoldableToolPanel,  public rtengine::DehazListener, public CurveListener,
    public AdjusterListener

{

protected:
    CurveEditorGroup* curveEditorGD;
    CurveEditorGroup* curveEditorGDH;
    Adjuster* str;
    Adjuster* scal;
    Adjuster* neigh;
    Adjuster* gain;
    Adjuster* offs;
    Adjuster* vart;
    Adjuster* limd;
    MyExpander* expretinex;

    Gtk::Label* labmdh;
    Gtk::HBox* dhbox;
    MyComboBoxText*   dehazmet;
    MyComboBoxText*   dehazcolorspace;
    Gtk::CheckButton* retinex;
    Gtk::Frame* dehazFrame;
    Gtk::CheckButton* medianmap;
    double nextmin;
    double nextmax;
    double nextminiT;
    double nextmaxiT;
    double nextmeanT;
    double nextminT;
    double nextmaxT;
    double nextsigma;

    Gtk::Label* mMLabels;
    Gtk::Label* transLabels;

    DiagonalCurveEditor* cdshape;
    DiagonalCurveEditor* cdshapeH;
    CurveEditorGroup* transmissionCurveEditorG;
    sigc::connection dehazmetConn;
    sigc::connection dehazColorSpaceConn;
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
    void autoOpenCurve         ();
    void retinexChanged        ();
    void medianmapChanged        ();
    void minmaxChanged (double cdma, double cdmin, double mini, double maxi, double Tmean, double Tsigma, double Tmin, double Tmax);
    bool minmaxComputed_ ();
    void updateLabel      ();
    void updateTrans      ();

    void enabledChanged        ();
    void curveChanged          (CurveEditor* ce);
    void dehazmetChanged();
    void dehazColorSpaceChanged();
    void retinexUpdateUI();
    void ColorSpaceUpdateUI();
    void writeOptions (std::vector<int> &tpOpen);
    void updateToolState (std::vector<int> &tpOpen);
    void setAdjusterBehavior (bool strAdd, bool neighAdd, bool scalAdd, bool limdAdd, bool gainAdd, bool offsAdd, bool vartAdd);
 
private:
    void foldAllButMe (GdkEventButton* event, MyExpander *expander);

};

#endif
