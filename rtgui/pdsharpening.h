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

#include <gtkmm.h>
#include "adjuster.h"
#include "thresholdadjuster.h"
#include "toolpanel.h"

class PdSharpening : public ToolParamBlock, public ThresholdAdjusterListener, public AdjusterListener, public FoldableToolPanel
{

protected:
    Adjuster* contrast;
    Adjuster* gamma;
    Adjuster* dradius;
    Adjuster* diter;
    Gtk::VBox* rld;

    rtengine::ProcEvent EvPdShrContrast;
    rtengine::ProcEvent EvPdShrDRadius;
    rtengine::ProcEvent EvPdSharpenGamma;
    rtengine::ProcEvent EvPdShrDIterations;
public:

    PdSharpening ();
    ~PdSharpening () override;

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void enabledChanged  () override;

    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override;
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;

    void setAdjusterBehavior (bool contrastadd, bool gammaadd, bool radiusadd, bool amountadd, bool dampingadd, bool iteradd, bool edgetoladd, bool haloctrladd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
};
