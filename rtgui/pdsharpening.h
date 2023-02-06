/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Ingo Weyrich (heckflosse67@gmx.de)
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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include "adjuster.h"
#include "checkbox.h"
#include "toolpanel.h"

class PdSharpening final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::AutoContrastListener,
    public rtengine::AutoRadiusListener,
    public CheckBoxListener
{

protected:
    Adjuster* contrast;
    Adjuster* dradius;
    Adjuster* dradiusOffset;
    Adjuster* diter;
    CheckBox* itercheck;

    bool lastAutoContrast;
    bool lastAutoRadius;
    rtengine::ProcEvent EvPdShrContrast;
    rtengine::ProcEvent EvPdShrCheckIter;
    rtengine::ProcEvent EvPdShrDRadius;
    rtengine::ProcEvent EvPdShrDRadiusOffset;
    rtengine::ProcEvent EvPdShrDIterations;
    rtengine::ProcEvent EvPdShrAutoContrast;
    rtengine::ProcEvent EvPdShrAutoRadius;
    IdleRegister idle_register;

public:
    static const Glib::ustring TOOL_NAME;

    PdSharpening ();
    ~PdSharpening () override;

    void read (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode (bool batchMode) override;

    void adjusterAutoToggled (Adjuster* a) override;
    void adjusterChanged (Adjuster* a, double newval) override;
    void enabledChanged () override;

    void autoContrastChanged (double autoContrast) override;
    void autoRadiusChanged (double autoRadius) override;

    void setAdjusterBehavior (bool contrastadd, bool radiusadd, bool iteradd);
    void trimValues (rtengine::procparams::ProcParams* pp) override;
    void checkBoxToggled(CheckBox* c, CheckValue newval) override;
};
