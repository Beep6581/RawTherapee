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
#include "guiutils.h"
#include "toolpanel.h"

class PreProcess final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:
    Gtk::CheckButton* hotPixel;
    Gtk::CheckButton* deadPixel;
    bool lastHot, lastDead;
    sigc::connection hpixelconn;
    sigc::connection dpixelconn;
    Adjuster* hdThreshold;
public:
    static const Glib::ustring TOOL_NAME;

    PreProcess ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    //void setBatchMode   (bool batchMode);
    //void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);

    void hotPixelChanged();
    void deadPixelChanged();
    void adjusterChanged(Adjuster* a, double newval) override;


    //void adjusterChanged     (Adjuster* a, double newval);
    //void setAdjusterBehavior (bool linedenoiseadd, bool greenequiladd);
    //void trimValues          (rtengine::procparams::ProcParams* pp);
};
