/*
 *  This file is part of RawTherapee.
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
 *  (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
 */
#pragma once

#include <gtkmm.h>

#include "adjuster.h"
#include "colorprovider.h"
#include "thresholdadjuster.h"
#include "toolpanel.h"

class DirPyrEqualizer final :
    public ToolParamBlock,
    public ThresholdAdjusterListener,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:

    Gtk::CheckButton * gamutlab;
    Adjuster* multiplier[6];
    Adjuster* threshold;
    Adjuster* skinprotect;
    ThresholdAdjuster* hueskin;
    //  MyComboBoxText*   algo;
    //  sigc::connection  algoconn;
    //  Gtk::Label*       alLabel;
    //  Gtk::Box*         algoHBox;

    sigc::connection  gamutlabConn;
    sigc::connection lumaneutralPressedConn;
    sigc::connection lumacontrastPlusPressedConn;
    sigc::connection lumacontrastMinusPressedConn;
    sigc::connection  cbdlMethodConn;
    Gtk::Label* labmcd;
    Gtk::Box* cdbox;
    MyComboBoxText*   cbdlMethod;

    bool lastgamutlab;

public:
    static const Glib::ustring TOOL_NAME;

    DirPyrEqualizer ();
    ~DirPyrEqualizer () override;

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode        (bool batchMode) override;
    void setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool skinadd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void cbdlMethodChanged();
    void adjusterChanged (Adjuster* a, double newval) override;
    void enabledChanged() override;
    void gamutlabToggled ();
    void lumaneutralPressed ();
    void lumacontrastPlusPressed ();
    void lumacontrastMinusPressed ();

    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override;
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override;
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override;
};
