/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>

#include "adjuster.h"
#include "editcallbacks.h"
#include "editwidgets.h"
#include "lensgeomlistener.h"
#include "toolpanel.h"
#include "../rtengine/perspectivecorrection.h"

struct ControlLine
{
    static constexpr int OBJ_COUNT = 4;
    std::unique_ptr<Line> line;
    std::shared_ptr<OPIcon> icon;
    std::shared_ptr<OPIcon> icon_h, icon_v;
    std::unique_ptr<Circle> begin, end;
    rtengine::ControlLine::Type type;
};

class ControlLineManager: EditSubscriber
{

protected:
    /** Determine how horizontal and vertical lines are displayed. */
    bool active_h, active_v;
    /** Hidden object for capturing mouse events. */
    std::unique_ptr<Rectangle> canvas_area;
    rtengine::Coord drag_delta;
    std::vector<std::unique_ptr<ControlLine>> control_lines;
    CursorShape cursor;
    bool draw_mode;
    Cairo::RefPtr<RTSurface> line_icon_h, line_icon_v;
    Cairo::RefPtr<RTSurface> line_icon_h_prelight, line_icon_v_prelight;
    int prev_obj;
    int selected_object;

    void addLine (rtengine::Coord begin, rtengine::Coord end);
    Geometry::State calcLineState(const ControlLine& line) const;
    void removeLine (size_t line_id);

public:
    class Callbacks
    {
    public:
        virtual ~Callbacks() {};
        /** Called when the EditSubscriber's switchOffEditMode is called. */
        virtual void switchOffEditMode (void) {};
    };

    /** Callbacks to invoke. */
    std::shared_ptr<Callbacks> callbacks;
    /** Type of line for newly drawn lines. */
    rtengine::ControlLine::Type draw_line_type;

    ControlLineManager();

    void removeAll (void);
    /** Sets whether or not the lines are visible and interact-able. */
    void setActive (bool active);
    /** Set whether or not lines can be drawn and deleted. */
    void setDrawMode (bool draw);
    void setEditProvider (EditDataProvider* provider);
    /** Determines how each line type is displayed. */
    void setLinesState (bool horiz_active, bool vert_active);
    /** Returns the number of lines. */
    size_t size (void) const;
    /**
     * Allocates a new array and populates it with copies of the control lines.
     */
    void toControlLines (std::vector<rtengine::ControlLine>& converted) const;

    // EditSubscriber overrides
    bool button1Pressed (int modifierKey) override;
    bool button1Released (void) override;
    bool button3Pressed (int modifierKey) override;
    bool pick1 (bool picked) override;
    bool pick3 (bool picked) override;
    bool drag1 (int modifierKey) override;
    CursorShape getCursor (int objectID) const override;
    bool mouseOver (int modifierKey) override;
    void switchOffEditMode (void) override;
};

class PerspCorrection final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:
    bool render = true;
    MyComboBoxText* method;
    Gtk::VBox* simple;
    Adjuster* horiz;
    Adjuster* vert;
    Gtk::Button* auto_pitch;
    Gtk::Button* auto_yaw;
    Gtk::Button* auto_pitch_yaw;
    Gtk::VBox* camera_based;
    Adjuster* camera_crop_factor;
    Adjuster* camera_focal_length;
    Adjuster* camera_pitch;
    Adjuster* camera_roll;
    Adjuster* camera_shift_horiz;
    Adjuster* camera_shift_vert;
    Adjuster* camera_yaw;
    std::unique_ptr<Gtk::Image> img_ctrl_lines_edit;
    std::unique_ptr<Gtk::Image> img_ctrl_lines_apply;
    std::unique_ptr<ControlLineManager> lines;
    Gtk::ToggleButton* lines_button_edit;
    Gtk::Button* lines_button_erase;
    Gtk::ToggleButton* lines_button_h;
    Gtk::ToggleButton* lines_button_v;
    Adjuster* projection_pitch;
    Adjuster* projection_rotate;
    Adjuster* projection_shift_horiz;
    Adjuster* projection_shift_vert;
    Adjuster* projection_yaw;
    rtengine::ProcEvent EvPerspCamFocalLength;
    rtengine::ProcEvent EvPerspCamShift;
    rtengine::ProcEvent EvPerspCamAngle;
    rtengine::ProcEvent EvPerspMethod;
    rtengine::ProcEvent EvPerspProjShift;
    rtengine::ProcEvent EvPerspProjRotate;
    rtengine::ProcEvent EvPerspProjAngle;
    rtengine::ProcEvent EvPerspRender;
    rtengine::ProcEvent EvPerspCamFocalLengthVoid;
    rtengine::ProcEvent EvPerspCamShiftVoid;
    rtengine::ProcEvent EvPerspCamAngleVoid;
    rtengine::ProcEvent EvPerspProjShiftVoid;
    rtengine::ProcEvent EvPerspProjRotateVoid;
    rtengine::ProcEvent EvPerspProjAngleVoid;
    rtengine::ProcEvent* event_persp_cam_focal_length;
    rtengine::ProcEvent* event_persp_cam_shift;
    rtengine::ProcEvent* event_persp_cam_angle;
    rtengine::ProcEvent* event_persp_proj_shift;
    rtengine::ProcEvent* event_persp_proj_rotate;
    rtengine::ProcEvent* event_persp_proj_angle;
    LensGeomListener* lens_geom_listener;
    const rtengine::FramesMetaData* metadata;

    void applyControlLines (void);
    void setCamBasedEventsActive (bool active = true);
    void setFocalLengthValue (const rtengine::procparams::ProcParams* pparams, const rtengine::FramesMetaData* metadata);

public:

    PerspCorrection ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void autoCorrectionPressed (Gtk::Button* b);
    void linesButtonPressed (Gtk::ToggleButton* button);
    void linesEditButtonPressed (void);
    void linesEraseButtonPressed (void);
    void methodChanged (void);
    void setAdjusterBehavior (bool badd, bool camera_focal_length_add, bool camera_shift_add, bool camera_angle_add, bool projection_angle_add, bool projection_shift_add, bool projection_rotate_add);
    void setEditProvider (EditDataProvider* provider) override;
    void setLensGeomListener (LensGeomListener* listener)
    {
        lens_geom_listener = listener;
    }
    void setMetadata (const rtengine::FramesMetaData* metadata);
    void switchOffEditMode (void);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
};

class LinesCallbacks: public ControlLineManager::Callbacks
{
protected:
    PerspCorrection* tool;

public:
    LinesCallbacks(PerspCorrection* tool);
    void switchOffEditMode (void) override;
};
