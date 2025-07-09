/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Jacques DESMIS <jdesmis@gmail.com>
 *  Copyright (c) 2018 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

#include <lcms2.h>

#include <gtkmm.h>
#include "adjuster.h"
#include <vector>

class RTWindow;

class ICCProfileCreator final : public Gtk::Dialog, public AdjusterListener
{

private:

    enum class ColorTemp {
        D50 = 5003,  // for Widegamut, Prophoto Best, Beta -> D50
        D60 = 6005,  // for ACESc                          -> D60
        D65 = 6504   // for sRGB, AdobeRGB, Bruce Rec2020  -> D65
    };

    cmsFloat64Number ga[7]; // 7 parameters for smoother curves

    //------------------------ Params -----------------------
    Glib::ustring primariesPreset;
    double redPrimaryX;
    double redPrimaryY;
    double greenPrimaryX;
    double greenPrimaryY;
    double bluePrimaryX;
    double bluePrimaryY;
    Glib::ustring gammaPreset;
    double gamma;
    double slope;
    bool appendParamsToDesc;
    bool v2except;
    Glib::ustring profileVersion;
    Glib::ustring illuminant;
    Glib::ustring description;
    Glib::ustring copyright;
    //-------------------------------------------------------

    RTWindow *parent;

    Adjuster* aGamma;
    Adjuster* aSlope;
    Adjuster* aPrimariesRedX;
    Adjuster* aPrimariesRedY;
    Adjuster* aPrimariesGreenX;
    Adjuster* aPrimariesGreenY;
    Adjuster* aPrimariesBlueX;
    Adjuster* aPrimariesBlueY;

    Gtk::Grid* primariesGrid;
    MyComboBoxText* iccVersion;
    MyComboBoxText* trcPresets;
    sigc::connection trcpresetsconn;
    MyComboBoxText* primaries;
    sigc::connection primariesconn;
    MyComboBoxText* cIlluminant;
    sigc::connection illconn;
    Gtk::Entry* eDescription;
    Gtk::Entry* eCopyright;
    Gtk::Button* resetCopyright;
    Gtk::CheckButton *cAppendParamsToDesc;

    //Glib::ustring lastPath;

    void initWithDefaults ();
    void storeDefaults ();
    void storeValues();

    void updateICCVersion();
    void primariesChanged();
    void illuminantChanged();
    void trcPresetsChanged();
    void adjusterChanged(Adjuster* a, double newval) override;
    static std::vector<Glib::ustring> getGamma();
    Glib::ustring getPrimariesPresetName(const Glib::ustring &preset);
    void getPrimaries(const Glib::ustring &preset, double *p, ColorTemp &temp);
    Glib::ustring getGammaPresetName(const Glib::ustring &preset);
    void getGamma(const Glib::ustring &preset, double &gamma, double &slope);
    void savePressed();
    void closePressed();
    void onResetCopyright();

public:
    explicit ICCProfileCreator (RTWindow *rtwindow);
};
