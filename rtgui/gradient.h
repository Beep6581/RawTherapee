/*
 *  This file is part of RawTherapee.
 */
#pragma once

#include <gtkmm.h>

#include "adjuster.h"
#include "editcallbacks.h"
#include "guiutils.h"
#include "toolpanel.h"

class Gradient final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public EditSubscriber
{

private:
    int lastObject;

protected:
    Gtk::Box *editHBox;
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
    void releaseEdit();

public:
    static const Glib::ustring TOOL_NAME;

    Gradient ();
    ~Gradient () override;

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void enabledChanged  () override;
    void setAdjusterBehavior (bool degreeadd, bool featheradd, bool strengthadd, bool centeradd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void updateGeometry  (const int centerX, const int centerY, const double feather, const double degree, const int fullWidth=-1, const int fullHeight=-1);

    void setEditProvider (EditDataProvider* provider) override;

    // EditSubscriber interface
    CursorShape getCursor(int objectID, int xPos, int yPos) const override;
    bool mouseOver(int modifierKey) override;
    bool button1Pressed(int modifierKey) override;
    bool button1Released() override;
    bool drag1(int modifierKey) override;
    void switchOffEditMode () override;
};
