/*
 *  This file is part of RawTherapee.
 */
#ifndef _GAMMA_H_
#define _GAMMA_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "thresholdadjuster.h"
#include "colorprovider.h"

class Gamma : public ToolParamBlock, public FoldableToolPanel, public CurveListener,
    public AdjusterListener

{

protected:
    Adjuster* gamm;
    Adjuster* slop;

    Gtk::Label* labmga;
    Gtk::HBox* gabox;

    MyComboBoxText*   gammaMethod;
    bool lastoutp;
    Gtk::CheckButton*  outp;

    sigc::connection gammaMethodConn, outpconn;

public:
    Gamma   ();
    ~Gamma  ();

    void read                  (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write                 (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);
    void setBatchMode          (bool batchMode);
    void setDefaults           (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = NULL);
    void trimValues            (rtengine::procparams::ProcParams* pp);
    void adjusterChanged       (Adjuster* a, double newval);

    void enabledChanged        ();
    void gammaMethodChanged();
    void setAdjusterBehavior (bool gammAdd, bool slopeAdd);
    void outpChanged ();

private:

};

#endif
