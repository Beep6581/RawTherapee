/*
 *  This file is part of RawTherapee.
 */
#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "edit.h"
#include "guiutils.h"

class Gradient : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public EditSubscriber
{

private:
    int lastObject;

protected:
    Gtk::HBox *editHBox;
    Gtk::ToggleButton* edit;
    Adjuster* degree;
    Adjuster* feather;
    Adjuster* strength;
    Adjuster* centerX;
    Adjuster* centerY;
    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    rtengine::Coord draggedCenter;
    sigc::connection editConn;

    void editToggled ();

public:

    Gradient ();
    ~Gradient ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged  ();
    void setAdjusterBehavior (bool degreeadd, bool featheradd, bool strengthadd, bool centeradd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void updateGeometry  (const int centerX, const int centerY, const double feather, const double degree, const int fullWidth=-1, const int fullHeight=-1);

    void setEditProvider (EditDataProvider* provider);

    // EditSubscriber interface
    CursorShape getCursor(const int objectID);
    bool mouseOver(const int modifierKey);
    bool button1Pressed(const int modifierKey);
    bool button1Released();
    bool drag1(const int modifierKey);
    void switchOffEditMode ();
};

#endif
