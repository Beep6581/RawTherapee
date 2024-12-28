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
#include "thresholdadjuster.h"
#include "toolpanel.h"

class Sharpening final:
    public ToolParamBlock,
    public ThresholdAdjusterListener,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:
    Adjuster* contrast;
    Adjuster* blur;
    MyComboBoxText* method;
    Adjuster* dradius;
    Adjuster* damount;
    Adjuster* ddamping;
    Adjuster* diter;
    Gtk::Box* usm;
    Gtk::Box* rld;

    Adjuster* radius;
    Adjuster* amount;
    Adjuster* eradius;
    Adjuster* etolerance;
    Adjuster* hcamount;
    Gtk::Box* edgebin;
    Gtk::Box* hcbin;
    Gtk::Box* edgebox;
    Gtk::Box* hcbox;
    ThresholdAdjuster* threshold;
    Gtk::CheckButton* edgesonly;
    bool lastEdgesOnly;
    sigc::connection eonlyConn;
    Gtk::CheckButton* halocontrol;
    bool lastHaloControl;
    sigc::connection hcConn;

    rtengine::ProcEvent EvSharpenContrast;
    rtengine::ProcEvent EvSharpenBlur;
public:
    static const Glib::ustring TOOL_NAME;

    Sharpening ();
    ~Sharpening () override;

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged(Adjuster* a, double newval) override;
    void enabledChanged  () override;
    void edgesonly_toggled ();
    void halocontrol_toggled ();
    void method_changed ();

    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override;
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;

    void setAdjusterBehavior (bool contrastadd, bool radiusadd, bool amountadd, bool dampingadd, bool iteradd, bool edgetoladd, bool haloctrladd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
};
