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
#include "lensgeomlistener.h"
#include "toolpanel.h"

class PerspCorrection final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:
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
    LensGeomListener* lens_geom_listener;
    const rtengine::FramesMetaData* metadata;

    void setFocalLengthValue (const rtengine::procparams::ProcParams* pparams, const rtengine::FramesMetaData* metadata);

public:

    PerspCorrection ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void autoCorrectionPressed (Gtk::Button* b);
    void methodChanged (void);
    void setAdjusterBehavior (bool badd, bool camera_focal_length_add, bool camera_shift_add, bool camera_angle_add, bool projection_angle_add, bool projection_shift_add, bool projection_rotate_add);
    void setLensGeomListener (LensGeomListener* listener)
    {
        lens_geom_listener = listener;
    }
    void setMetadata (const rtengine::FramesMetaData* metadata);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
};
