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

class Distortion final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:
    Gtk::Button*   autoDistor;
    Adjuster* distor;
    Adjuster* focal_length;
    sigc::connection    idConn;
    LensGeomListener * rlistener;
    Gtk::CheckButton* defish;

public:
    rtengine::ProcEvent EvDistortionDefish;
    rtengine::ProcEvent EvDistortionDefishVoid;
    rtengine::ProcEvent* event_distortion_defish;

    rtengine::ProcEvent EvDistortionDefishFocalLength;
    rtengine::ProcEvent EvDistortionDefishFocalLengthVoid;
    rtengine::ProcEvent* event_distortion_defish_focal_length;

    static const Glib::ustring TOOL_NAME;

    Distortion ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged     (Adjuster* a, double newval) override;
    void setAdjusterBehavior (bool vadd, bool focal_length_add);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void idPressed           ();
    void setLensGeomListener (LensGeomListener* l)
    {
        rlistener = l;
    }
    void defishChanged (void);

    void setCamBasedEventsActive(bool active = true);
};
