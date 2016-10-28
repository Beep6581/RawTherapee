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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
 */

#ifndef DIRPYREQUALIZER_H_INCLUDED
#define DIRPYREQUALIZER_H_INCLUDED

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "thresholdadjuster.h"
#include "colorprovider.h"

class DirPyrEqualizer : public ToolParamBlock, public ThresholdAdjusterListener, public AdjusterListener, public FoldableToolPanel
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
    //  Gtk::HBox*        algoHBox;

    sigc::connection  gamutlabConn;
    sigc::connection lumaneutralPressedConn;
    sigc::connection lumacontrastPlusPressedConn;
    sigc::connection lumacontrastMinusPressedConn;
    sigc::connection  cbdlMethodConn;
    Gtk::Label* labmcd;
    Gtk::HBox* cdbox;
    MyComboBoxText*   cbdlMethod;

    bool lastgamutlab;

public:

    DirPyrEqualizer ();
    virtual ~DirPyrEqualizer ();

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void setBatchMode        (bool batchMode);
    void setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool skinadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void adjusterChanged (ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight);
//    void algoChanged         ();
    void cbdlMethodChanged();
    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged();
    void gamutlabToggled ();
    void lumaneutralPressed ();
    void lumacontrastPlusPressed ();
    void lumacontrastMinusPressed ();
};

#endif
