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

#include "lensgeomlistener.h"
#include "toolpanel.h"
#include "adjuster.h"

class LensGeometry final :
    public ToolParamBlock,
    public FoldableToolPanel,
    public AdjusterListener
{

protected:
    MyComboBoxText*     method;
    Gtk::Button*        autoCrop;
    LensGeomListener*   rlistener;
    Adjuster*           scale;
    Gtk::CheckButton*   fill;
    bool                lastFill;
    sigc::connection    fillConn;

    rtengine::ProcEvent EvTransMethod;
    rtengine::ProcEvent EvTransScale;
public:
    static const Glib::ustring TOOL_NAME;

    LensGeometry ();
    ~LensGeometry () override;

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void methodChanged();
    void fillPressed            ();
    void autoCropPressed        ();
    void setLensGeomListener    (LensGeomListener* l)
    {
        rlistener = l;
    }

    void adjusterChanged (Adjuster* a, double newval) override;

private:
    IdleRegister idle_register;
};
