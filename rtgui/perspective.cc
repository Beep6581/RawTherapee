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
#include "perspective.h"

#include "rtimage.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

PerspCorrection::PerspCorrection () : FoldableToolPanel(this, "perspective", M("TP_PERSPECTIVE_LABEL"))
{

    Gtk::Image* ipersHL =   Gtk::manage (new RTImage ("perspective-horizontal-left-small.png"));
    Gtk::Image* ipersHR =   Gtk::manage (new RTImage ("perspective-horizontal-right-small.png"));
    Gtk::Image* ipersVL =   Gtk::manage (new RTImage ("perspective-vertical-bottom-small.png"));
    Gtk::Image* ipersVR =   Gtk::manage (new RTImage ("perspective-vertical-top-small.png"));

    Gtk::Image* ipers_proj_yaw_left = Gtk::manage (new RTImage ("perspective-horizontal-left-small.png"));
    Gtk::Image* ipers_proj_yaw_right = Gtk::manage (new RTImage ("perspective-horizontal-right-small.png"));
    Gtk::Image* ipers_proj_pitch_left = Gtk::manage (new RTImage ("perspective-vertical-bottom-small.png"));
    Gtk::Image* ipers_proj_pitch_right = Gtk::manage (new RTImage ("perspective-vertical-top-small.png"));
    Gtk::Image* ipers_rotate_left = Gtk::manage(new RTImage("rotate-right-small.png"));
    Gtk::Image* ipers_rotate_right = Gtk::manage(new RTImage("rotate-left-small.png"));

    vert = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_VERTICAL"), -85, 85, 0.1, 0, ipersVL, ipersVR));
    vert->setAdjusterListener (this);

    horiz = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_HORIZONTAL"), -85, 85, 0.1, 0, ipersHL, ipersHR));
    horiz->setAdjusterListener (this);

    camera_focal_length = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_CAMERA_FOCAL_LENGTH"), 0.5, 2000, 0.01, 24));
    camera_focal_length->setAdjusterListener (this);

    camera_crop_factor = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_CAMERA_CROP_FACTOR"), 0.1, 30, 0.01, 1));
    camera_crop_factor->setAdjusterListener (this);

    camera_shift_horiz = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_CAMERA_SHIFT_HORIZONTAL"), -100, 100, 0.01, 0));
    camera_shift_horiz->setAdjusterListener (this);

    camera_shift_vert = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_CAMERA_SHIFT_VERTICAL"), -100, 100, 0.01, 0));
    camera_shift_vert->setAdjusterListener (this);

    projection_shift_horiz = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_PROJECTION_SHIFT_HORIZONTAL"), -100, 100, 0.01, 0));
    projection_shift_horiz->setAdjusterListener (this);

    projection_shift_vert = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_PROJECTION_SHIFT_VERTICAL"), -100, 100, 0.01, 0));
    projection_shift_vert->setAdjusterListener (this);

    projection_rotate = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_PROJECTION_ROTATE"), -180, 180, 0.01, 0, ipers_rotate_left, ipers_rotate_right));
    projection_rotate->setAdjusterListener (this);

    projection_scale = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_PROJECTION_SCALE"), 0.5, 2, 0.01, 1));
    projection_scale->setAdjusterListener (this);

    projection_pitch = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_PROJECTION_PITCH"), -85, 85, 0.1, 0, ipers_proj_pitch_left, ipers_proj_pitch_right));
    projection_pitch->setAdjusterListener (this);

    projection_yaw = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_PROJECTION_YAW"), -85, 85, 0.1, 0, ipers_proj_yaw_left, ipers_proj_yaw_right));
    projection_yaw->setAdjusterListener (this);

    pack_start (*vert);
    pack_start (*horiz);
    pack_start (*camera_focal_length);
    pack_start (*camera_crop_factor);
    pack_start (*camera_shift_horiz);
    pack_start (*camera_shift_vert);
    pack_start (*projection_shift_horiz);
    pack_start (*projection_shift_vert);
    pack_start (*projection_rotate);
    pack_start (*projection_scale);
    pack_start (*projection_yaw);
    pack_start (*projection_pitch);

    horiz->setLogScale(2, 0);
    vert->setLogScale(2, 0);
    camera_focal_length->setLogScale(4000, 0.5);
    camera_crop_factor->setLogScale(300, 0.1);
    projection_scale->setLogScale(4, 0.5);

    show_all();
}

void PerspCorrection::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        horiz->setEditedState (pedited->perspective.horizontal ? Edited : UnEdited);
        vert->setEditedState (pedited->perspective.vertical ? Edited : UnEdited);
        camera_crop_factor->setEditedState (pedited->perspective.camera_crop_factor ? Edited : UnEdited);
        camera_focal_length->setEditedState (pedited->perspective.camera_focal_length ? Edited : UnEdited);
        camera_shift_horiz->setEditedState (pedited->perspective.camera_shift_horiz ? Edited : UnEdited);
        camera_shift_vert->setEditedState (pedited->perspective.camera_shift_vert ? Edited : UnEdited);
        projection_pitch->setEditedState (pedited->perspective.projection_pitch ? Edited : UnEdited);
        projection_rotate->setEditedState (pedited->perspective.projection_rotate ? Edited : UnEdited);
        projection_scale->setEditedState (pedited->perspective.projection_scale ? Edited : UnEdited);
        projection_shift_horiz->setEditedState (pedited->perspective.projection_shift_horiz ? Edited : UnEdited);
        projection_shift_vert->setEditedState (pedited->perspective.projection_shift_vert ? Edited : UnEdited);
        projection_yaw->setEditedState (pedited->perspective.projection_yaw ? Edited : UnEdited);
    }

    horiz->setValue (pp->perspective.horizontal);
    vert->setValue (pp->perspective.vertical);
    camera_crop_factor->setValue (pp->perspective.camera_crop_factor);
    camera_focal_length->setValue (pp->perspective.camera_focal_length);
    camera_shift_horiz->setValue (pp->perspective.camera_shift_horiz);
    camera_shift_vert->setValue (pp->perspective.camera_shift_vert);
    projection_pitch->setValue (pp->perspective.projection_pitch);
    projection_rotate->setValue (pp->perspective.projection_rotate);
    projection_scale->setValue (pp->perspective.projection_scale);
    projection_shift_horiz->setValue (pp->perspective.projection_shift_horiz);
    projection_shift_vert->setValue (pp->perspective.projection_shift_vert);
    projection_yaw->setValue (pp->perspective.projection_yaw);

    enableListener ();
}

void PerspCorrection::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->perspective.horizontal  = horiz->getValue ();
    pp->perspective.vertical = vert->getValue ();
    pp->perspective.camera_crop_factor= camera_crop_factor->getValue ();
    pp->perspective.camera_focal_length = camera_focal_length->getValue ();
    pp->perspective.camera_shift_horiz = camera_shift_horiz->getValue ();
    pp->perspective.camera_shift_vert = camera_shift_vert->getValue ();
    pp->perspective.projection_pitch = projection_pitch->getValue ();
    pp->perspective.projection_rotate = projection_rotate->getValue ();
    pp->perspective.projection_scale = projection_scale->getValue ();
    pp->perspective.projection_shift_horiz = projection_shift_horiz->getValue ();
    pp->perspective.projection_shift_vert = projection_shift_vert->getValue ();
    pp->perspective.projection_yaw = projection_yaw->getValue ();

    if (pedited) {
        pedited->perspective.horizontal = horiz->getEditedState ();
        pedited->perspective.vertical = vert->getEditedState ();
        pedited->perspective.camera_crop_factor= camera_crop_factor->getEditedState ();
        pedited->perspective.camera_focal_length = camera_focal_length->getEditedState ();
        pedited->perspective.camera_shift_horiz = camera_shift_horiz->getEditedState();
        pedited->perspective.camera_shift_vert = camera_shift_vert->getEditedState();
        pedited->perspective.projection_pitch = projection_pitch->getEditedState();
        pedited->perspective.projection_rotate = projection_rotate->getEditedState();
        pedited->perspective.projection_scale = projection_scale->getEditedState();
        pedited->perspective.projection_shift_horiz = projection_shift_horiz->getEditedState();
        pedited->perspective.projection_shift_vert = projection_shift_vert->getEditedState();
        pedited->perspective.projection_yaw = projection_yaw->getEditedState();
    }
}

void PerspCorrection::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    horiz->setDefault (defParams->perspective.horizontal);
    vert->setDefault (defParams->perspective.vertical);
    camera_crop_factor->setDefault (defParams->perspective.camera_crop_factor);
    camera_focal_length->setDefault (defParams->perspective.camera_focal_length);
    camera_shift_horiz->setDefault (defParams->perspective.camera_shift_horiz);
    camera_shift_vert->setDefault (defParams->perspective.camera_shift_vert);
    projection_pitch->setDefault (defParams->perspective.projection_pitch);
    projection_rotate->setDefault (defParams->perspective.projection_rotate);
    projection_scale->setDefault (defParams->perspective.projection_scale);
    projection_shift_horiz->setDefault (defParams->perspective.projection_shift_horiz);
    projection_shift_vert->setDefault (defParams->perspective.projection_shift_vert);
    projection_yaw->setDefault (defParams->perspective.projection_yaw);

    if (pedited) {
        horiz->setDefaultEditedState (pedited->perspective.horizontal ? Edited : UnEdited);
        vert->setDefaultEditedState (pedited->perspective.vertical ? Edited : UnEdited);
        camera_crop_factor->setDefaultEditedState (pedited->perspective.camera_crop_factor ? Edited : UnEdited);
        camera_focal_length->setDefaultEditedState (pedited->perspective.camera_focal_length ? Edited : UnEdited);
        camera_shift_horiz->setDefaultEditedState (pedited->perspective.camera_shift_horiz ? Edited : UnEdited);
        camera_shift_vert->setDefaultEditedState (pedited->perspective.camera_shift_vert ? Edited : UnEdited);
        projection_pitch->setDefaultEditedState (pedited->perspective.projection_pitch ? Edited : UnEdited);
        projection_rotate->setDefaultEditedState (pedited->perspective.projection_rotate ? Edited : UnEdited);
        projection_scale->setDefaultEditedState (pedited->perspective.projection_scale ? Edited : UnEdited);
        projection_shift_horiz->setDefaultEditedState (pedited->perspective.projection_shift_horiz ? Edited : UnEdited);
        projection_shift_vert->setDefaultEditedState (pedited->perspective.projection_shift_vert ? Edited : UnEdited);
        projection_yaw->setDefaultEditedState (pedited->perspective.projection_yaw ? Edited : UnEdited);
    } else {
        horiz->setDefaultEditedState (Irrelevant);
        vert->setDefaultEditedState (Irrelevant);
        camera_crop_factor->setDefaultEditedState (Irrelevant);
        camera_focal_length->setDefaultEditedState (Irrelevant);
        camera_shift_horiz->setDefaultEditedState (Irrelevant);
        camera_shift_vert->setDefaultEditedState (Irrelevant);
        projection_pitch->setDefaultEditedState (Irrelevant);
        projection_rotate->setDefaultEditedState (Irrelevant);
        projection_scale->setDefaultEditedState (Irrelevant);
        projection_shift_horiz->setDefaultEditedState (Irrelevant);
        projection_shift_vert->setDefaultEditedState (Irrelevant);
        projection_yaw->setDefaultEditedState (Irrelevant);
    }
}

void PerspCorrection::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        if (a == horiz || a == vert) {
            listener->panelChanged (EvPerspCorr,
                    Glib::ustring::compose("%1=%2\n%3=%4",
                        M("TP_PERSPECTIVE_HORIZONTAL"),
                        horiz->getValue(),
                        M("TP_PERSPECTIVE_VERTICAL"),
                        vert->getValue()));
        } else if (a == camera_focal_length || a == camera_crop_factor) {
            listener->panelChanged (EvPerspCamFocalLength,
                    Glib::ustring::compose("%1=%2\n%3=%4",
                        M("TP_PERSPECTIVE_CAMERA_FOCAL_LENGTH"),
                        camera_focal_length->getValue(),
                        M("TP_PERSPECTIVE_CAMERA_CROP_FACTOR"),
                        camera_crop_factor->getValue()));
        } else if (a == camera_shift_horiz || a == camera_shift_vert) {
            listener->panelChanged (EvPerspCamShift,
                    Glib::ustring::compose("%1=%2\n%3=%4",
                        M("TP_PERSPECTIVE_CAMERA_SHIFT_HORIZONTAL"),
                        camera_shift_horiz->getValue(),
                        M("TP_PERSPECTIVE_CAMERA_SHIFT_VERTICAL"),
                        camera_shift_vert->getValue()));
        } else if (a == projection_shift_horiz || a == projection_shift_vert) {
            listener->panelChanged (EvPerspProjShift,
                    Glib::ustring::compose("%1=%2\n%3=%4",
                        M("TP_PERSPECTIVE_PROJECTION_SHIFT_HORIZONTAL"),
                        projection_shift_horiz->getValue(),
                        M("TP_PERSPECTIVE_PROJECTION_SHIFT_VERTICAL"),
                        projection_shift_vert->getValue()));
        } else if (a == projection_rotate) {
            listener->panelChanged (EvPerspProjRotate,
                    Glib::ustring::format(projection_rotate->getValue()));
        } else if (a == projection_scale) {
            listener->panelChanged (EvPerspProjScale,
                    Glib::ustring::format(projection_scale->getValue()));
        } else if (a == projection_pitch || a == projection_yaw) {
            listener->panelChanged (EvPerspProjAngle,
                    Glib::ustring::compose("%1=%2\n%3=%4",
                        M("TP_PERSPECTIVE_PROJECTION_PITCH"),
                        projection_pitch->getValue(),
                        M("TP_PERSPECTIVE_PROJECTION_YAW"),
                        projection_yaw->getValue()));
        }
    }
}

void PerspCorrection::setAdjusterBehavior (bool badd, bool camera_focal_length_add, bool camera_shift_add, bool projection_angle_add, bool projection_shift_add, bool projection_rotate_add, bool projection_scale_add)
{

    horiz->setAddMode(badd);
    vert->setAddMode(badd);
    camera_crop_factor->setAddMode(camera_focal_length_add);
    camera_focal_length->setAddMode(camera_focal_length_add);
    camera_shift_horiz->setAddMode(camera_shift_add);
    camera_shift_vert->setAddMode(camera_shift_add);
    projection_pitch->setAddMode(projection_angle_add);
    projection_rotate->setAddMode(projection_rotate_add);
    projection_scale->setAddMode(projection_scale_add);
    projection_shift_horiz->setAddMode(projection_shift_add);
    projection_shift_vert->setAddMode(projection_shift_add);
    projection_yaw->setAddMode(projection_angle_add);
}

void PerspCorrection::trimValues (rtengine::procparams::ProcParams* pp)
{

    horiz->trimValue(pp->perspective.horizontal);
    vert->trimValue(pp->perspective.vertical);
    camera_crop_factor->trimValue(pp->perspective.camera_crop_factor);
    camera_focal_length->trimValue(pp->perspective.camera_focal_length);
    camera_shift_horiz->trimValue(pp->perspective.camera_shift_horiz);
    camera_shift_vert->trimValue(pp->perspective.camera_shift_vert);
    projection_pitch->trimValue(pp->perspective.projection_pitch);
    projection_rotate->trimValue(pp->perspective.projection_rotate);
    projection_scale->trimValue(pp->perspective.projection_scale);
    projection_shift_horiz->trimValue(pp->perspective.projection_shift_horiz);
    projection_shift_vert->trimValue(pp->perspective.projection_shift_vert);
    projection_yaw->trimValue(pp->perspective.projection_yaw);
}

void PerspCorrection::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    horiz->showEditedCB ();
    vert->showEditedCB ();
    camera_crop_factor->showEditedCB ();
    camera_focal_length->showEditedCB ();
    camera_shift_horiz->showEditedCB ();
    camera_shift_vert->showEditedCB ();
    projection_pitch->showEditedCB ();
    projection_rotate->showEditedCB ();
    projection_scale->showEditedCB ();
    projection_shift_horiz->showEditedCB ();
    projection_shift_vert->showEditedCB ();
    projection_yaw->showEditedCB ();
}
