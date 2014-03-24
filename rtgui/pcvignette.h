/*
 *  This file is part of RawTherapee.
 */
#ifndef _PCVIGNETTE_H_
#define _PCVIGNETTE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class PCVignette : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel {

  protected:
    Gtk::CheckButton* enabled;
    Adjuster* strength;
    Adjuster* feather;
    Adjuster* roundness;
    bool lastEnabled;
    sigc::connection enaConn;

  public:

    PCVignette ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool strengthadd, bool featheradd, bool roundnessadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
};

#endif
