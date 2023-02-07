/* -*- C++ -*-
 *
 *  Adapted from ART.
 *
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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>
#include "adjuster.h"
#include "checkbox.h"
#include "toolpanel.h"

class ToneEqualizer: public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public CheckBoxListener {
public:
    static const Glib::ustring TOOL_NAME;

    ToneEqualizer();

    void read(const rtengine::procparams::ProcParams *pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams *pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams *defParams, const ParamsEdited* pedited = nullptr) override;
    void adjusterChanged(Adjuster *a, double newval) override;
    void adjusterAutoToggled(Adjuster *a) override;
    void enabledChanged() override;
    void setBatchMode(bool batchMode) override;
    void setAdjusterBehavior(bool bands_add, bool regularization_add, bool pivot_add);
    void checkBoxToggled(CheckBox* c, CheckValue newval) override;

    void trimValues(rtengine::procparams::ProcParams *pp) override;

private:
    void colormapToggled();

    std::array<Adjuster *, 5> bands;
    Adjuster *regularization;
    Adjuster *pivot;
    CheckBox *show_colormap;

    rtengine::ProcEvent EvEnabled;
    rtengine::ProcEvent EvBands;
    rtengine::ProcEvent EvRegularization;
    rtengine::ProcEvent EvColormap;
    rtengine::ProcEvent EvPivot;

    rtengine::procparams::ToneEqualizerParams inital_params;
};
