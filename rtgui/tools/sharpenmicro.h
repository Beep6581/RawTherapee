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
 *
 *
 *  Manuel Llorens' algorithm of micro-contrast sharpening
 *
 *
 */
#pragma once

#include <gtkmm.h>

#include "adjuster.h"
#include "toolpanel.h"

class SharpenMicro final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:

    Gtk::CheckButton* matrix;
    Adjuster* amount;
    Adjuster* uniformity;
    Adjuster* contrast;

    rtengine::ProcEvent EvSharpenMicroContrast;

    sigc::connection matrixconn;
    bool lastmatrix;

public:
    static const Glib::ustring TOOL_NAME;

    SharpenMicro           ();

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode        (bool batchMode) override;
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void setAdjusterBehavior (bool amountadd, bool contrastadd, bool uniformityadd);
    void adjusterChanged     (Adjuster* a, double newval) override;

    void enabledChanged      () override;
    void matrix_toggled      ();


};
