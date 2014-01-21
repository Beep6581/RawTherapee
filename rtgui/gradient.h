/*
 *  This file is part of RawTherapee.
 */
#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class Gradient : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel {

  protected:
    Gtk::CheckButton* enabled;
    Adjuster* degree;
    Adjuster* feather;
    Adjuster* strength;
    Adjuster* centerX;
    Adjuster* centerY;
    bool lastEnabled;
    sigc::connection enaConn;

  public:

    Gradient ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool degreeadd, bool featheradd, bool strengthadd, bool centeradd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
};

#endif
