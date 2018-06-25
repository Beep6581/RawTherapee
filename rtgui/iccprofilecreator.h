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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>
#include "adjuster.h"
#include "options.h"
#include <vector>
#include "rtwindow.h"

class ICCProfileCreator : public Gtk::Dialog, public AdjusterListener
{

private:

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

    Glib::ustring profileVersion;
    Glib::ustring illuminant;
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

    //Glib::ustring lastPath;

    void initWithDefaults ();
    void storeDefaults ();
    void storeValues();

    void updateICCVersion();
    void primariesChanged();
    void illuminantChanged();
    void trcPresetsChanged();
    static std::vector<Glib::ustring> getGamma();
    void getGammaArray();
    void savePressed();
    void closePressed();

public:
    explicit ICCProfileCreator (RTWindow *rtwindow);
};
