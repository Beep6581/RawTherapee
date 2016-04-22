/*
 *  This file is part of RawTherapee.
 */
#ifndef _SPOT_H_
#define _SPOT_H_

#include <gtkmm.h>
#include "toolpanel.h"
#include "edit.h"

/**
 * @brief Let the user create/edit/delete points for Spot Removal tool
 *
 * User Interface:
 *
 * For the rest of this documentation, T represent a "target" point (where the image is edited) and
 * S represent the "source" location (where the edition takes its source data).
 *
 * When the edit button is active, all T points are shown by a small "dot". When the user
 * move the cursor over one of them, a circle is displayed to show the radius of the brush, as well
 * as a second circle representing the source data (S point). The user can then use the left mouse button
 * over the icon to drag the T point. The left mouse button can be used over the S circle or the right
 * mouse button can be used over the T point to move the S point.
 *
 * Using the left mouse button over the circle of the T point will let the user adjust its radius.
 *
 * Using the left mouse button over the feather circle will let the user adjust its radius by setting
 * a coefficient (0.0 = same radius than the inner circle ; 1.0 = 2 times the inner radius).
 *
 * To create a new point, just move over a free area, and press the left mouse button while holding
 * the CTRL key. This will create a new S and T pair of points. The CTRL key can be released, but keep
 * the left mouse button pressed and move away to position the S point.
 *
 * To delete a point, move your mouse over any of its geometry press the middle or right mouse button
 * (the point will be deleted on button release).
 */
class Spot : public ToolParamBlock, public FoldableToolPanel, public EditSubscriber
{

private:
    int lastObject;                // current object that is hovered
    int activeSpot;                // currently active spot, being edited
    std::vector<rtengine::SpotEntry> spots; // list of edited spots
    OPIcon sourceIcon;             // to show the source location
    Circle sourceCircle;           // to show and change the Source radius
    Circle sourceMODisc;           // to change the Source position
    Circle targetCircle;           // to show and change the Target radius
    Circle targetMODisc;           // to change the Target position
    Circle sourceFeatherCircle;    // to show the Feather radius at the Source position
    Circle targetFeatherCircle;    // to show the Feather radius at the Target position
    Line link;                     // to show the link between the Source and Target position

    OPIcon *getActiveSpotIcon ();
    void updateGeometry ();
    void createGeometry ();
    void addNewEntry ();
    void deleteSelectedEntry ();
    void resetPressed ();

protected:
    Gtk::HBox* labelBox;
    Gtk::CheckButton* editedCheckBox;
    Gtk::Label* countLabel;
    Gtk::ToggleButton* edit;
    Gtk::Button* reset;
    sigc::connection editConn, editedConn;

    void editToggled ();
    void editedToggled ();
    Geometry* getVisibleGeometryFromMO (int MOID);

public:

    Spot ();
    ~Spot ();

    void read (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);

    void enabledChanged ();

    void setEditProvider (EditDataProvider* provider);

    void setBatchMode (bool batchMode);

    // EditSubscriber interface
    CursorShape getCursor (const int objectID);
    bool mouseOver (const int modifierKey);
    bool button1Pressed (const int modifierKey);
    bool button1Released ();
    bool button2Pressed (const int modifierKey);
    bool button3Pressed (const int modifierKey);
    bool button3Released ();
    bool drag1 (const int modifierKey);
    bool drag3 (const int modifierKey);
    bool pick2 (const bool picked);
    bool pick3 (const bool picked);
    void switchOffEditMode ();
};

#endif
