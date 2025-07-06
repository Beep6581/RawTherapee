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
#include "controllines.h"
#include "lensgeomlistener.h"
#include "toolpanel.h"
#include "rtengine/tweakoperator.h"

class PerspCorrectionPanelListener
{
public:
    virtual ~PerspCorrectionPanelListener() = default;

    virtual void controlLineEditModeChanged(bool active) = 0;
};

class PerspCorrection final :
    public rtengine::TweakOperator,
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:
    bool render = true;
    MyComboBoxText* method;
    Gtk::Box* simple;
    Adjuster* horiz;
    Adjuster* vert;
    Gtk::Button* auto_pitch;
    Gtk::Button* auto_yaw;
    Gtk::Button* auto_pitch_yaw;
    Gtk::Box* camera_based;
    Adjuster* camera_crop_factor;
    Adjuster* camera_focal_length;
    Adjuster* camera_pitch;
    Adjuster* camera_roll;
    Adjuster* camera_shift_horiz;
    Adjuster* camera_shift_vert;
    Adjuster* camera_yaw;
    std::unique_ptr<ControlLineManager> lines;
    Gtk::Button* lines_button_apply;
    Gtk::ToggleButton* lines_button_edit;
    Gtk::Button* lines_button_erase;
    Adjuster* projection_pitch;
    Adjuster* projection_rotate;
    Adjuster* projection_shift_horiz;
    Adjuster* projection_shift_vert;
    Adjuster* projection_yaw;
    rtengine::ProcEvent EvPerspCamFocalLength;
    rtengine::ProcEvent EvPerspCamShift;
    rtengine::ProcEvent EvPerspCamAngle;
    rtengine::ProcEvent EvPerspControlLines;
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
    PerspCorrectionPanelListener* panel_listener;
    const rtengine::FramesMetaData* metadata;

    void applyControlLines (void);
    void tweakParams(rtengine::procparams::ProcParams &pparams) override;
    void setCamBasedEventsActive (bool active = true);
    void setFocalLengthValue (const rtengine::procparams::ProcParams* pparams, const rtengine::FramesMetaData* metadata);
    void updateApplyDeleteButtons();

public:

    /** Minimum number of horizontal lines for horizontal/full correction. */
    static constexpr std::size_t MIN_HORIZ_LINES = 2;
    /** Minimum number of vertical lines for vertical/full correction. */
    static constexpr std::size_t MIN_VERT_LINES = 2;
    static const Glib::ustring TOOL_NAME;

    PerspCorrection ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void autoCorrectionPressed (Gtk::Button* b);
    void lineChanged (void);
    void linesApplyButtonPressed (void);
    void linesEditButtonPressed (void);
    void linesEraseButtonPressed (void);
    void methodChanged (void);
    void requestApplyControlLines(void);
    void setAdjusterBehavior (bool badd,
        bool camera_focal_length_add,
        bool camera_shift_add,
        bool camera_angle_add,
        bool projection_angle_add,
        bool projection_shift_add,
        bool projection_rotate_add);
    void setControlLineEditMode(bool active);
    void setEditProvider (EditDataProvider* provider) override;
    void setLensGeomListener (LensGeomListener* listener)
    {
        lens_geom_listener = listener;
    }
    void setPerspCorrectionPanelListener(PerspCorrectionPanelListener* listener)
    {
        panel_listener = listener;
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
    explicit LinesCallbacks(PerspCorrection* tool);
    void lineChanged (void) override;
    void switchOffEditMode (void) override;
};
