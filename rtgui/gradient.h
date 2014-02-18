/*
 *  This file is part of RawTherapee.
 */
#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "edit.h"

class Gradient : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel, public EditSubscriber {

  private:
    int lastObject;

  protected:
    Gtk::CheckButton* enabled;
    Gtk::ToggleButton* edit;
    Adjuster* degree;
    Adjuster* feather;
    Adjuster* strength;
    Adjuster* centerX;
    Adjuster* centerY;
    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    Coord draggedCenter;
    bool lastEnabled;
    sigc::connection enaConn, editConn;

    void editToggled ();

  public:

    Gradient ();
    ~Gradient ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void updateGeometry (int centerX_, int centerY_, double strength_, double degree_);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool degreeadd, bool featheradd, bool strengthadd, bool centeradd);
    void trimValues          (rtengine::procparams::ProcParams* pp);

    void setEditProvider (EditDataProvider* provider);

    // EditSubscriber interface
    CursorShape getCursor(int objectID);
    bool mouseOver(int modifierKey);
    bool button1Pressed(int modifierKey);
    bool button1Released();
    bool drag(int modifierKey);
    void switchOffEditMode ();
};

#endif
