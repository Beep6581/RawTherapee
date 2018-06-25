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
#include <sigc++/slot.h>
#include "iccprofilecreator.h"
#include "multilangmgr.h"
#include "cachemanager.h"
#include "addsetids.h"
#include "../rtengine/icons.h"
#include "../rtengine/color.h"
#include "rtimage.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern Options options;

namespace rtengine
{

extern const Settings* settings;

}

const char* sTRCPreset[] = {"BT709_g2.2_s4.5", "sRGB_g2.4_s12.92", "linear_g1.0", "standard_g2.2", "standard_g1.8", "High_g1.3_s3.35", "Low_g2.6_s6.9", "Lab_g3.0s9.03296"}; //gamma free

ICCProfileCreator::ICCProfileCreator(RTWindow *rtwindow)
    : Gtk::Dialog (M ("MAIN_BUTTON_ICCPROFCREATOR"), *rtwindow, true)
    , primariesPreset(options.ICCPC_primariesPreset)
    , redPrimaryX(options.ICCPC_redPrimaryX)
    , redPrimaryY(options.ICCPC_redPrimaryY)
    , greenPrimaryX(options.ICCPC_greenPrimaryX)
    , greenPrimaryY(options.ICCPC_greenPrimaryY)
    , bluePrimaryX(options.ICCPC_bluePrimaryX)
    , bluePrimaryY(options.ICCPC_bluePrimaryY)
    , gammaPreset(options.ICCPC_gammaPreset)
    , gamma(options.ICCPC_gamma)
    , slope(options.ICCPC_slope)
    , profileVersion(options.ICCPC_profileVersion)
    , illuminant(options.ICCPC_illuminant)
    , parent(rtwindow)
{

    set_default_size(600, -1);

    Gtk::Grid* mainGrid = Gtk::manage(new Gtk::Grid());

    //--------------------------------- primaries

    Gtk::Label* prilab = Gtk::manage(new Gtk::Label(M("ICCPROFCREATOR_PRIMARIES")));
    setExpandAlignProperties(prilab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    mainGrid->attach(*prilab, 0, 0, 1, 1);

    primaries = Gtk::manage(new MyComboBoxText());
    setExpandAlignProperties(primaries, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    primaries->append(M("ICCPROFCREATOR_PRIM_ACESP0"));
    primaries->append(M("ICCPROFCREATOR_PRIM_ACESP1"));
    primaries->append(M("ICCPROFCREATOR_PRIM_ADOBE"));
    primaries->append(M("ICCPROFCREATOR_PRIM_PROPH"));
    primaries->append(M("ICCPROFCREATOR_PRIM_REC2020"));
    primaries->append(M("ICCPROFCREATOR_PRIM_SRGB"));
    primaries->append(M("ICCPROFCREATOR_PRIM_WIDEG"));
    primaries->append(M("ICCPROFCREATOR_PRIM_BEST"));
    primaries->append(M("ICCPROFCREATOR_PRIM_BETA"));
    primaries->append(M("ICCPROFCREATOR_PRIM_BRUCE"));
    primaries->append(M("ICCPROFCREATOR_CUSTOM"));
    primaries->set_tooltip_text(M("ICCPROFCREATOR_PRIM_TOOLTIP"));
    mainGrid->attach(*primaries, 1, 0, 1, 1);

    primariesGrid = Gtk::manage(new Gtk::Grid());
    setExpandAlignProperties(primariesGrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    primariesGrid->set_column_spacing(5);

    /*
    Gtk::Image* gamuts0 =  Gtk::manage(new RTImage("rt-logo-tiny.png"));
    Gtk::Image* gamutl0 =  Gtk::manage(new RTImage("rt-logo-small.png"));
    Gtk::Image* gamuts1 =  Gtk::manage(new RTImage("rt-logo-tiny.png"));
    Gtk::Image* gamutl1 =  Gtk::manage(new RTImage("rt-logo-small.png"));
    Gtk::Image* gamuts2 =  Gtk::manage(new RTImage("rt-logo-tiny.png"));
    Gtk::Image* gamutl2 =  Gtk::manage(new RTImage("rt-logo-small.png"));
    Gtk::Image* gamuts3 =  Gtk::manage(new RTImage("rt-logo-tiny.png"));
    Gtk::Image* gamutl3 =  Gtk::manage(new RTImage("rt-logo-small.png"));
    Gtk::Image* gamuts4 =  Gtk::manage(new RTImage("rt-logo-tiny.png"));
    Gtk::Image* gamutl4 =  Gtk::manage(new RTImage("rt-logo-small.png"));
    Gtk::Image* gamuts5 =  Gtk::manage(new RTImage("rt-logo-tiny.png"));
    Gtk::Image* gamutl5 =  Gtk::manage(new RTImage("rt-logo-small.png"));
    */

    aPrimariesRedX = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_PRIM_REDX"), 0.6300, 0.7350, 0.0001, 0.6400/*, gamuts0, gamutl0*/));
    setExpandAlignProperties(aPrimariesRedX, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    aPrimariesRedY = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_PRIM_REDY"), 0.2650, 0.3350, 0.0001, 0.3300/*, gamutl1, gamuts1*/));
    setExpandAlignProperties(aPrimariesRedY, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    aPrimariesGreenX = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_PRIM_GREX"), 0.0000, 0.3100, 0.0001, 0.3000/*, gamutl2, gamuts2*/));
    setExpandAlignProperties(aPrimariesGreenX, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    aPrimariesGreenY = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_PRIM_GREY"), 0.5900, 1.0000, 0.0001, 0.6000/*, gamuts3, gamutl3*/));
    setExpandAlignProperties(aPrimariesGreenY, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    aPrimariesBlueX = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_PRIM_BLUX"), 0.0001, 0.1600, 0.0001, 0.1500/*, gamutl4, gamuts4*/));
    setExpandAlignProperties(aPrimariesBlueX, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    aPrimariesBlueY = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_PRIM_BLUY"), -0.0700, 0.0700, 0.0001, 0.060/*, gamutl5, gamuts5*/));
    setExpandAlignProperties(aPrimariesBlueY, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    primariesGrid->attach(*aPrimariesRedX,   0, 0, 1, 1);
    primariesGrid->attach(*aPrimariesRedY,   1, 0, 1, 1);

    primariesGrid->attach(*aPrimariesGreenX, 0, 1, 1, 1);
    primariesGrid->attach(*aPrimariesGreenY, 1, 1, 1, 1);

    primariesGrid->attach(*aPrimariesBlueX,  0, 2, 1, 1);
    primariesGrid->attach(*aPrimariesBlueY,  1, 2, 1, 1);

    mainGrid->attach(*primariesGrid, 1, 1, 1, 1);

    //--------------------------------- output gamma

    Gtk::Label* galab = Gtk::manage(new Gtk::Label(M("ICCPROFCREATOR_TRC_PRESET")));
    setExpandAlignProperties(galab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    mainGrid->attach(*galab, 0, 2, 1, 1);

    trcPresets = Gtk::manage(new MyComboBoxText());
    setExpandAlignProperties(trcPresets, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    std::vector<Glib::ustring> outputTRCPresets;
    outputTRCPresets.push_back(M("ICCPROFCREATOR_CUSTOM"));
    for (unsigned int i = 0; i < sizeof(sTRCPreset) / sizeof(sTRCPreset[0]); i++) {
        outputTRCPresets.push_back(sTRCPreset[i]);
    }
    for (size_t i = 0; i < outputTRCPresets.size(); i++) {
        trcPresets->append(outputTRCPresets[i]);
    }
    mainGrid->attach(*trcPresets, 1, 2, 1, 1);

    //--------------------------------- sliders gampos and slpos

    aGamma = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_GAMMA"), 1, 3.5, 0.01, 2.4));
    setExpandAlignProperties(aGamma, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);

    if (aGamma->delay < options.adjusterMaxDelay) {
        aGamma->delay = options.adjusterMaxDelay;
    }
    aGamma->show();
    mainGrid->attach(*aGamma, 1, 3, 1, 1); //gamma

    aSlope = Gtk::manage(new Adjuster(M("ICCPROFCREATOR_SLOPE"), 0, 15, 0.00001, 12.92310));
    setExpandAlignProperties(aSlope, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);

    if (aSlope->delay < options.adjusterMaxDelay) {
        aSlope->delay = options.adjusterMaxDelay;
    }
    aSlope->show();
    mainGrid->attach(*aSlope, 1, 4, 1, 1); //slope

    //--------------------------------- temperature

    Gtk::Label* illlab = Gtk::manage(new Gtk::Label(M("ICCPROFCREATOR_ILL")));
    setExpandAlignProperties(illlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    mainGrid->attach(*illlab, 0, 5, 1, 1); //slope
    cIlluminant = Gtk::manage(new MyComboBoxText());
    setExpandAlignProperties(cIlluminant, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    cIlluminant->append(M("ICCPROFCREATOR_ILL_DEF"));
    cIlluminant->append(M("ICCPROFCREATOR_ILL_41"));
    cIlluminant->append(M("ICCPROFCREATOR_ILL_50"));
    cIlluminant->append(M("ICCPROFCREATOR_ILL_55"));
    cIlluminant->append(M("ICCPROFCREATOR_ILL_60"));
    cIlluminant->append(M("ICCPROFCREATOR_ILL_65"));
    cIlluminant->append(M("ICCPROFCREATOR_ILL_80"));
    cIlluminant->append(M("ICCPROFCREATOR_ILL_INC"));
    cIlluminant->set_tooltip_text(M("ICCPROFCREATOR_ILL_TOOLTIP"));
    mainGrid->attach(*cIlluminant, 1, 5, 1, 1);

    //--------------------------------- V2  or V4 profiles

    Gtk::Label* proflab = Gtk::manage(new Gtk::Label(M("ICCPROFCREATOR_ICCVERSION")));
    setExpandAlignProperties(proflab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    mainGrid->attach(*proflab, 0, 6, 1, 1);
    iccVersion = Gtk::manage(new MyComboBoxText());
    setExpandAlignProperties(iccVersion, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    iccVersion->append(M("ICCPROFCREATOR_PROF_V4"));
    iccVersion->append(M("ICCPROFCREATOR_PROF_V2"));
    mainGrid->attach(*iccVersion, 1, 6, 1, 1);

    //--------------------------------- Adding the mainGrid

    get_content_area()->add(*mainGrid);

    //--------------------------------- Setting default values for Adjusters

    aGamma->setDefault(options.ICCPC_gamma);
    aSlope->setDefault(options.ICCPC_slope);
    aPrimariesRedX->setDefault(options.ICCPC_redPrimaryX);
    aPrimariesRedY->setDefault(options.ICCPC_redPrimaryY);
    aPrimariesGreenX->setDefault(options.ICCPC_greenPrimaryX);
    aPrimariesGreenY->setDefault(options.ICCPC_greenPrimaryY);
    aPrimariesBlueX->setDefault(options.ICCPC_bluePrimaryX);
    aPrimariesBlueY->setDefault(options.ICCPC_bluePrimaryY);

    //--------------- Updating widgets with actual values (from options)

    if (primariesPreset == "ACES-AP0") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_ACESP0"));
    } else if (primariesPreset == "ACES-AP1") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_ACESP1"));
    } else if (primariesPreset == "Adobe") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_ADOBE"));
    } else if (primariesPreset == "ProPhoto") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_PROPH"));
    } else if (primariesPreset == "Rec2020") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_REC2020"));
    } else if (primariesPreset == "sRGB") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_SRGB"));
    } else if (primariesPreset == "Widegamut") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_WIDEG"));
    } else if (primariesPreset == "BestRGB") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_BEST"));
    } else if (primariesPreset == "BetaRGB") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_BETA"));
    } else if (primariesPreset == "BruceRGB") {
        primaries->set_active_text(M("ICCPROFCREATOR_PRIM_BRUCE"));
    } else if (primariesPreset == "custom") {
        primaries->set_active_text(M("ICCPROFCREATOR_CUSTOM"));
    }

    trcPresets->set_active(0);
    if (gammaPreset != "Custom") {
        trcPresets->set_active_text(gammaPreset);
    }

    aGamma->setValue(options.ICCPC_gamma);
    aSlope->setValue(options.ICCPC_slope);
    aPrimariesRedX->setValue(redPrimaryX);
    aPrimariesRedY->setValue(redPrimaryY);
    aPrimariesGreenX->setValue(greenPrimaryX);
    aPrimariesGreenY->setValue(greenPrimaryY);
    aPrimariesBlueX->setValue(bluePrimaryX);
    aPrimariesBlueY->setValue(bluePrimaryY);

    if (illuminant == "DEF") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_DEF"));
    } else if (illuminant == "D41") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_41"));
    } else if (illuminant == "D50") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_50"));
    } else if (illuminant == "D55") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_55"));
    } else if (illuminant == "D60") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_60"));
    } else if (illuminant == "D65") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_65"));
    } else if (illuminant == "D80") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_80"));
    } else if (illuminant == "stdA") {
        cIlluminant->set_active_text(M("ICCPROFCREATOR_ILL_INC"));
    }

    iccVersion->set_active(0);
    if (profileVersion == "v2") {
        iccVersion->set_active(1);
    }

    trcPresetsChanged();
    illuminantChanged();
    primariesChanged();

    //--------------- Action area button

    Gtk::Button* save = Gtk::manage (new Gtk::Button (M ("GENERAL_SAVE_AS")));
    save->signal_clicked().connect ( sigc::mem_fun (*this, &ICCProfileCreator::savePressed) );
    get_action_area()->pack_start (*save);

    Gtk::Button* close = Gtk::manage (new Gtk::Button (M ("GENERAL_CLOSE")));
    close->signal_clicked().connect ( sigc::mem_fun (*this, &ICCProfileCreator::closePressed) );
    get_action_area()->pack_start (*close);

    //--------------- Show childrens

    show_all_children ();

    //--------------- Connecting the signals

    aPrimariesRedX->setAdjusterListener(this);
    aPrimariesRedY->setAdjusterListener(this);
    aPrimariesGreenX->setAdjusterListener(this);
    aPrimariesGreenY->setAdjusterListener(this);
    aPrimariesBlueX->setAdjusterListener(this);
    aPrimariesBlueY->setAdjusterListener(this);
    aGamma->setAdjusterListener(this);
    aSlope->setAdjusterListener(this);
    primariesconn = primaries->signal_changed().connect(sigc::mem_fun(*this, &ICCProfileCreator::primariesChanged));
    trcpresetsconn = trcPresets->signal_changed().connect(sigc::mem_fun(*this, &ICCProfileCreator::trcPresetsChanged));
    illconn = cIlluminant->signal_changed().connect(sigc::mem_fun(*this, &ICCProfileCreator::illuminantChanged));
}

void ICCProfileCreator::closePressed()
{
    storeValues();
    hide();
}

void ICCProfileCreator::updateICCVersion()
{
    if (cIlluminant->get_active_text() != M("ICCPROFCREATOR_ILL_DEF") || primaries->get_active_text() == M("ICCPROFCREATOR_CUSTOM")) {
        iccVersion->set_active_text(M("ICCPROFCREATOR_PROF_V4"));
        iccVersion->set_sensitive(false);
    } else {
        iccVersion->set_sensitive(true);
    }
}

void ICCProfileCreator::primariesChanged()
{
    if (primaries->get_active_text() == M("ICCPROFCREATOR_CUSTOM")) {
        primariesGrid->set_sensitive(true);
    } else {
        primariesGrid->set_sensitive(false);
    }
    updateICCVersion();
}

void ICCProfileCreator::illuminantChanged()
{
    updateICCVersion();
}

void ICCProfileCreator::trcPresetsChanged()
{
    aGamma->block(true);
    aSlope->block(true);

    bool sensitive = trcPresets->get_active_row_number() == 0;
    aGamma->set_sensitive(sensitive);
    aSlope->set_sensitive(sensitive);

    aGamma->block(false);
    aSlope->block(false);
}

void ICCProfileCreator::storeValues()
{
    if (iccVersion->get_active_text() == M("ICCPROFCREATOR_PROF_V4")) {
        options.ICCPC_profileVersion = profileVersion = "v4";
    } else if (iccVersion->get_active_text() == M("ICCPROFCREATOR_PROF_V2")) {
        options.ICCPC_profileVersion = profileVersion = "v2";
    }

    options.ICCPC_gammaPreset = gammaPreset = trcPresets->get_active_text();
    if (gammaPreset == M("ICCPROFCREATOR_CUSTOM")) {
        options.ICCPC_gammaPreset = gammaPreset = "Custom";
    }

    if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_DEF")) {
        options.ICCPC_illuminant = illuminant = "DEF";
    } else if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_41")) {
        options.ICCPC_illuminant = illuminant = "D41";
    } else if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_50")) {
        options.ICCPC_illuminant = illuminant = "D50";
    } else if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_55")) {
        options.ICCPC_illuminant = illuminant = "D55";
    } else if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_60")) {
        options.ICCPC_illuminant = illuminant = "D60";
    } else if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_65")) {
        options.ICCPC_illuminant = illuminant = "D65";
    } else if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_80")) {
        options.ICCPC_illuminant = illuminant = "D80";
    } else if (cIlluminant->get_active_text() == M("ICCPROFCREATOR_ILL_INC")) {
        options.ICCPC_illuminant = illuminant = "stdA";
    }

    if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_ACESP0")) {
        options.ICCPC_primariesPreset = primariesPreset = "ACES-AP0";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_ACESP1")) {
        options.ICCPC_primariesPreset = primariesPreset = "ACES-AP1";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_ADOBE")) {
        options.ICCPC_primariesPreset = primariesPreset = "Adobe";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_PROPH")) {
        options.ICCPC_primariesPreset = primariesPreset = "ProPhoto";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_REC2020")) {
        options.ICCPC_primariesPreset = primariesPreset = "Rec2020";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_SRGB")) {
        options.ICCPC_primariesPreset = primariesPreset = "sRGB";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_WIDEG")) {
        options.ICCPC_primariesPreset = primariesPreset = "Widegamut";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_BEST")) {
        options.ICCPC_primariesPreset = primariesPreset = "BestRGB";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_BETA")) {
        options.ICCPC_primariesPreset = primariesPreset = "BetaRGB";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_BRUCE")) {
        options.ICCPC_primariesPreset = primariesPreset = "BruceRGB";
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_CUSTOM")) {
        options.ICCPC_primariesPreset = primariesPreset = "custom";
    }

    options.ICCPC_gamma = gamma = aGamma->getValue();
    options.ICCPC_slope = slope = aSlope->getValue();
    options.ICCPC_redPrimaryX = redPrimaryX = aPrimariesRedX->getValue();
    options.ICCPC_redPrimaryY = redPrimaryY = aPrimariesRedY->getValue();
    options.ICCPC_greenPrimaryX = greenPrimaryX = aPrimariesGreenX->getValue();
    options.ICCPC_greenPrimaryY = greenPrimaryY = aPrimariesGreenY->getValue();
    options.ICCPC_bluePrimaryX = bluePrimaryX = aPrimariesBlueX->getValue();
    options.ICCPC_bluePrimaryY = bluePrimaryY = aPrimariesBlueY->getValue();

}

// Copyright (c) 2018 Jacques DESMIS <jdesmis@gmail.com>
// WARNING: the caller must lock lcmsMutex
void ICCProfileCreator::savePressed()
{
    bool pro = false;
    cmsHPROFILE newProfile = nullptr;
    Glib::ustring sNewProfile;
    Glib::ustring sPrimariesPreset;
    Glib::ustring sGammaPreset;

    storeValues();

    // -------------------------------------------- Compute de default file name

    if (gammaPreset == "linear_g1.0" || (gammaPreset == "High_g1.3_s3.35")) {
        pro = true;    //pro=0  RT_sRGB || Prophoto
    }

    //necessary for V2 profile
    if (primariesPreset == "ACES-AP0"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.ACESp0)) {
        sNewProfile = options.rtSettings.ACESp0;
        sPrimariesPreset = "ACES-AP0_";
    } else if (primariesPreset == "ACES-AP1"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.ACESp1)) {
        sNewProfile = options.rtSettings.ACESp1;
        sPrimariesPreset = "ACES-AP1_";
    } else if (primariesPreset == "Adobe"      && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.adobe)) {
        sNewProfile = options.rtSettings.adobe;
        sPrimariesPreset = "Medium_";
    } else if (primariesPreset == "ProPhoto"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.prophoto)   && !pro) {
        sNewProfile = options.rtSettings.prophoto;
        sPrimariesPreset = "Large_";
    } else if (primariesPreset == "ProPhoto"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.prophoto10) &&  pro) {
        sNewProfile = options.rtSettings.prophoto10;
        sPrimariesPreset = "Large_";
    } else if (primariesPreset == "Rec2020"    && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.rec2020)) {
        sNewProfile = options.rtSettings.rec2020;
        sPrimariesPreset = "Rec2020_";
    } else if (primariesPreset == "sRGB"       && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.srgb)       && !pro) {
        sNewProfile = options.rtSettings.srgb;
        sPrimariesPreset = "sRGB_";
    } else if (primariesPreset == "sRGB"       && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.srgb10)     &&  pro) {
        sNewProfile = options.rtSettings.srgb10;
        sPrimariesPreset = "sRGB_";
    } else if (primariesPreset == "Widegamut"  && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.widegamut)) {
        sNewProfile = options.rtSettings.widegamut;
        sPrimariesPreset = "Wide_";
    } else if (primariesPreset == "BestRGB"    && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.best)) {
        sNewProfile = options.rtSettings.best;
        sPrimariesPreset = "Best_";
    } else if (primariesPreset == "BetaRGB"    && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.beta)) {
        sNewProfile = options.rtSettings.beta;
        sPrimariesPreset = "Beta_";
    } else if (primariesPreset == "BruceRGB"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.bruce)) {
        sNewProfile = options.rtSettings.bruce;
        sPrimariesPreset = "Bruce_";
    } else if (primariesPreset == "custom") {
        sNewProfile = options.rtSettings.srgb;
        sPrimariesPreset = "Custom_";
    } else {
        // Should not occurs
        if (rtengine::settings->verbose) {
            printf("\"%s\": unknown working profile! - use LCMS2 substitution\n", primariesPreset.c_str());
        }
        return;
    }

    //begin adaptation rTRC gTRC bTRC
    //"newProfile" profile has the same characteristics than RGB values, but TRC are adapted... for applying profile
    if (rtengine::settings->verbose) {
        printf("Output Gamma - profile Primaries as RT profile: \"%s\"\n", sNewProfile.c_str());
    }

    newProfile = rtengine::ICCStore::getInstance()->getProfile(sNewProfile); //get output profile

    if (newProfile == nullptr) {

        if (rtengine::settings->verbose) {
            printf("\"%s\" ICC output profile not found!\n", sNewProfile.c_str());
        }
        return;
    }

    //change desc Tag , to "free gamma", or "BT709", etc.
    Glib::ustring fName;
    Glib::ustring sPrimariesAndIlluminant;
    double presetGamma = 2.4;
    double presetSlope = 12.92310;
    const double eps = 0.000000001; // not divide by zero
    if (gammaPreset == "High_g1.3_s3.35") {
        sGammaPreset = "_High_g=1.3_s=3.35";
        presetGamma = 1.3;
        presetSlope = 3.35;
        ga[0] = 1.3 ;    //for high dynamic images
        ga[1] = 0.998279;
        ga[2] = 0.001721;
        ga[3] = 0.298507;
        ga[4] = 0.005746;
    } else if (gammaPreset == "Low_g2.6_s6.9") {
        sGammaPreset = "_Low_g=2.6_s=6.9";
        presetGamma = 2.6;
        presetSlope = 6.9;
        ga[0] = 2.6 ;    //gamma 2.6 variable : for low contrast images
        ga[1] = 0.891161;
        ga[2] = 0.108839;
        ga[3] = 0.144928;
        ga[4] = 0.076332;
    } else if (gammaPreset == "sRGB_g2.4_s12.92") {
        sGammaPreset = "_sRGB_g=2.4_s=12.92310";
        presetGamma = 2.4;
        presetSlope = 12.92310;
        ga[0] = 2.40;    //sRGB 2.4 12.92  - RT default as Lightroom
        ga[1] = 0.947858;
        ga[2] = 0.052142;
        ga[3] = 0.077399;
        ga[4] = 0.039293;
    } else if (gammaPreset == "BT709_g2.2_s4.5") {
        sGammaPreset = "_BT709_g=2.2_s=4.5";
        presetGamma = 2.22;
        presetSlope = 4.5;
        ga[0] = 2.22;    //BT709  2.2  4.5  - my preferred as D.Coffin
        ga[1] = 0.909995;
        ga[2] = 0.090005;
        ga[3] = 0.222222;
        ga[4] = 0.081071;
    } else if (gammaPreset == "linear_g1.0") {
        sGammaPreset = "_Linear_g=1.0";
        presetGamma = 1.;
        presetSlope = 0.;
        ga[0] = 1.0;    //gamma=1 linear : for high dynamic images (cf D.Coffin...)
        ga[1] = 1.;
        ga[2] = 0.;
        ga[3] = 1. / eps;
        ga[4] = 0.;
    } else if (gammaPreset == "standard_g2.2") {
        sGammaPreset = "_g=2.2";
        presetGamma = 2.2;
        presetSlope = 0.;
        ga[0] = 2.2;    //gamma=2.2(as gamma of Adobe, Widegamut...)
        ga[1] = 1.;
        ga[2] = 0.;
        ga[3] = 1. / eps;
        ga[4] = 0.;
    } else if (gammaPreset == "standard_g1.8") {
        sGammaPreset = "_g=1.8";
        presetGamma = 1.8;
        presetSlope = 0.;
        ga[0] = 1.8;    //gamma=1.8(as gamma of Prophoto)
        ga[1] = 1.;
        ga[2] = 0.;
        ga[3] = 1. / eps;
        ga[4] = 0.;
    } else if (gammaPreset == "Lab_g3.0s9.03296") {
        sGammaPreset = "_LAB_g3.0_s9.03296";
        presetGamma = 3.0;
        presetSlope = 9.03296;
        ga[0] = 3.0;    //Lab gamma =3 slope=9.03296
        ga[1] = 0.8621;
        ga[2] = 0.1379;
        ga[3] = 0.1107;
        ga[4] = 0.08;
    } else if (gammaPreset == "Custom") {
        rtengine::GammaValues g_a; //gamma parameters
        double pwr = 1.0 / gamma;
        double ts = slope;
        double slope2 = slope == 0 ? eps : slope;

        int mode = 0;
        rtengine::Color::calcGamma(pwr, ts, mode, g_a); // call to calcGamma with selected gamma and slope : return parameters for LCMS2
        ga[4] = g_a[3] * ts;
        //printf("g_a.gamma0=%f g_a.gamma1=%f g_a.gamma2=%f g_a.gamma3=%f g_a.gamma4=%f\n", g_a.gamma0,g_a.gamma1,g_a.gamma2,g_a.gamma3,g_a.gamma4);
        ga[0] = gamma;
        ga[1] = 1. /(1.0 + g_a[4]);
        ga[2] = g_a[4] /(1.0 + g_a[4]);
        ga[3] = 1. / slope2;
        //printf("ga[0]=%f ga[1]=%f ga[2]=%f ga[3]=%f ga[4]=%f\n", ga[0],ga[1],ga[2],ga[3],ga[4]);

        sGammaPreset = Glib::ustring::compose("_g%1_s%2",
                       Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), gamma),
                       Glib::ustring::format (std::setw(6), std::fixed, std::setprecision(5), slope));
        presetGamma = gamma;
        presetSlope = slope;
    }
    ga[5] = 0.0;
    ga[6] = 0.0;


    sPrimariesAndIlluminant = sPrimariesPreset;

    if (profileVersion == "v4" && illuminant != "DEF") {
         sPrimariesPreset += illuminant;
        //  printf("outpr=%s \n",outPr.c_str());
    }

    // create description with gamma + slope + primaries
    std::wostringstream gammaWs;
    std::wstring gammaStrICC;

    Glib::ustring gammaGS;//to save gamma and slope in a tag

    if (gammaPreset == "Custom") {
        Glib::ustring sGamma(Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), gamma));
        Glib::ustring sSlope(Glib::ustring::format (std::setw(6), std::fixed, std::setprecision(5), slope));
        fName = (profileVersion == "v4" ? "RTv4_" : "RTv2_") + sPrimariesAndIlluminant + sGamma + " " + sSlope + ".icc";
        gammaWs << sPrimariesPreset << " g=" << Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), gamma) << " s=" << Glib::ustring::format (std::setw(6), std::fixed, std::setprecision(5), slope);
        gammaGS = "g" + sGamma + "s" + sSlope + "!";


    } else {
        Glib::ustring sGamma(Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), presetGamma));
        Glib::ustring sSlope(Glib::ustring::format (std::setw(6), std::fixed, std::setprecision(5), presetSlope));
        fName = (profileVersion == "v4" ? "RTv4_" : "RTv2_") + sPrimariesAndIlluminant + sGammaPreset + ".icc";
        gammaWs << sPrimariesPreset << sGammaPreset;
        gammaGS = "g" + sGamma + "s" + sSlope + "!";
    }

    // -------------------------------------------- Asking the file name

    Gtk::FileChooserDialog dialog(getToplevelWindow(this), M("ICCPROFCREATOR_SAVEDIALOG_TITLE"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    bindCurrentFolder(dialog, options.lastICCProfCreatorDir);
    dialog.set_current_name(fName);
    //dialog.set_current_folder(lastPath);

    dialog.add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(M("GENERAL_SAVE"), Gtk::RESPONSE_OK);

    Glib::RefPtr<Gtk::FileFilter> filter_icc = Gtk::FileFilter::create();
    filter_icc->set_name(M("FILECHOOSER_FILTER_COLPROF"));
    filter_icc->add_pattern("*.icc");
    dialog.add_filter(filter_icc);

    /*
    Glib::RefPtr<Gtk::FileFilter> filter_any = Gtk::FileFilter::create();
    filter_any->set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any->add_pattern("*");
    dialog.add_filter(filter_any);
    */

    dialog.show_all_children();
    //dialog.set_do_overwrite_confirmation (true);

    Glib::ustring absoluteFName;

    do {
        int result = dialog.run();

        if (result != Gtk::RESPONSE_OK) {
            return;
        } else {
            absoluteFName = dialog.get_filename();
            Glib::ustring ext = getExtension(absoluteFName);

            if (ext != "icc") {
                absoluteFName += ".icc";
            }

            if (confirmOverwrite(dialog, absoluteFName)) {
                //lastPath = Glib::path_get_dirname(absoluteFName);
                break;
            }
        }
    } while (1);

    // -----------------------------------------------------------------

    //write in tag 'dmdd' values of current gamma and slope to retrive after in Output profile
    wchar_t *wGammaGS = (wchar_t*)g_utf8_to_utf16 (gammaGS.c_str(), -1, NULL, NULL, NULL);
    if (!wGammaGS) {
        printf("Error: lab2rgbOut  /  g_utf8_to_utf16 failed!\n");
    }

    cmsMLU *description = cmsMLUalloc(NULL, 1);
    // Language code (3 letters code) : https://www.iso.org/obp/ui/
    // Country code (3 letters code)  : http://www.loc.gov/standards/iso639-2/php/code_list.php
    if (cmsMLUsetWide(description, "eng", "USA", wGammaGS)) {
        cmsWriteTag(newProfile, cmsSigDeviceModelDescTag, description); //save 'dmdd' in description
    } else {
        printf("Error: lab2rgbOut  /  cmsMLUsetWide failed for \"%s\" !\n", gammaGS.c_str());
    }
    cmsMLU *mlu;
    cmsContext ContextID = cmsGetProfileContextID(newProfile); // create context to modify some TAGs
    mlu = cmsMLUalloc(ContextID, 1);
    cmsMLUsetWide(mlu,  "en", "US", gammaWs.str().c_str());

    cmsMLUfree(description);

    // instruction with //ICC are used to generate ICC profile
    if (mlu == nullptr) {
        printf("Description error\n");
    } else {

        if (profileVersion == "v4") {
            cmsSetProfileVersion(newProfile, 4.3);
        } else {
            cmsSetProfileVersion(newProfile, 2.0);
        }

//change
        float p[6]; //primaries
        ga[6] = 0.0;

        enum class ColorTemp {
            D50 = 5003,  // for Widegamut, Prophoto Best, Beta -> D50
            D65 = 6504,   // for sRGB, AdobeRGB, Bruce Rec2020  -> D65
            D60 = 6005        //for ACESc->D60
        };
        ColorTemp temp = ColorTemp::D50;

        if (primariesPreset == "Widegamut") {
            p[0] = 0.7350;    //Widegamut primaries
            p[1] = 0.2650;
            p[2] = 0.1150;
            p[3] = 0.8260;
            p[4] = 0.1570;
            p[5] = 0.0180;

        } else if (primariesPreset == "Adobe") {
            p[0] = 0.6400;    //Adobe primaries
            p[1] = 0.3300;
            p[2] = 0.2100;
            p[3] = 0.7100;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (primariesPreset == "sRGB") {
            p[0] = 0.6400;    // sRGB primaries
            p[1] = 0.3300;
            p[2] = 0.3000;
            p[3] = 0.6000;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (primariesPreset == "BruceRGB") {
            p[0] = 0.6400;    // Bruce primaries
            p[1] = 0.3300;
            p[2] = 0.2800;
            p[3] = 0.6500;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (primariesPreset == "BetaRGB") {
            p[0] = 0.6888;    // Beta primaries
            p[1] = 0.3112;
            p[2] = 0.1986;
            p[3] = 0.7551;
            p[4] = 0.1265;
            p[5] = 0.0352;
        } else if (primariesPreset == "BestRGB") {
            p[0] = 0.7347;    // Best primaries
            p[1] = 0.2653;
            p[2] = 0.2150;
            p[3] = 0.7750;
            p[4] = 0.1300;
            p[5] = 0.0350;
        } else if (primariesPreset == "Rec2020") {
            p[0] = 0.7080;    // Rec2020 primaries
            p[1] = 0.2920;
            p[2] = 0.1700;
            p[3] = 0.7970;
            p[4] = 0.1310;
            p[5] = 0.0460;
            temp = ColorTemp::D65;
        } else if (primariesPreset == "ACES-AP0") {
            p[0] = 0.7347;    // ACES P0 primaries
            p[1] = 0.2653;
            p[2] = 0.0000;
            p[3] = 1.0;
            p[4] = 0.0001;
            p[5] = -0.0770;
            temp = ColorTemp::D60;
        } else if (primariesPreset == "ACES-AP1") {
            p[0] = 0.713;    // ACES P1 primaries
            p[1] = 0.293;
            p[2] = 0.165;
            p[3] = 0.830;
            p[4] = 0.128;
            p[5] = 0.044;
            temp = ColorTemp::D60;
        } else if (primariesPreset == "ProPhoto") {
            p[0] = 0.7347;    // ProPhoto and default primaries
            p[1] = 0.2653;
            p[2] = 0.1596;
            p[3] = 0.8404;
            p[4] = 0.0366;
            p[5] = 0.0001;
        } else if (primariesPreset == "custom") {
            p[0] = redPrimaryX;
            p[1] = redPrimaryY;
            p[2] = greenPrimaryX;
            p[3] = greenPrimaryY;
            p[4] = bluePrimaryX;
            p[5] = bluePrimaryY;

        } else {
            p[0] = 0.7347;    //default primaries
            p[1] = 0.2653;
            p[2] = 0.1596;
            p[3] = 0.8404;
            p[4] = 0.0366;
            p[5] = 0.0001;
        }

        cmsCIExyY xyD;
        cmsCIExyYTRIPLE Primaries = {
            {p[0], p[1], 1.0}, // red
            {p[2], p[3], 1.0}, // green
            {p[4], p[5], 1.0}  // blue
        };
        double tempv4 = 5000.;

        if (profileVersion == "v4" && illuminant != "DEF") {
            if (illuminant == "D41") {
                tempv4 = 4100.;
            } else if (illuminant == "D50") {
                tempv4 = 5003.;
            } else if (illuminant == "D55") {
                tempv4 = 5500.;
            } else if (illuminant == "D60") {
                tempv4 = 6004.;
            } else if (illuminant == "D65") {
                tempv4 = 6504.;
            } else if (illuminant == "D80") {
                tempv4 = 8000.;
            } else if (illuminant == "stdA") {
                tempv4 = 5003.;
            }

            //printf("tempv4=%f \n", tempv4);

        }

        if (profileVersion == "v4" && illuminant != "DEF") {
            cmsWhitePointFromTemp(&xyD, tempv4);
        } else {
            cmsWhitePointFromTemp(&xyD, (double)temp);
        }

        if (illuminant == "stdA") xyD = {0.447573, 0.407440, 1.0};

        // Calculate output profile's rTRC gTRC bTRC
        cmsToneCurve* GammaTRC[3];
        GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(nullptr, 5, ga);

        if (profileVersion == "v4") {
            newProfile = cmsCreateRGBProfile(&xyD, &Primaries, GammaTRC);
        }

        cmsWriteTag(newProfile, cmsSigRedTRCTag, GammaTRC[0]);
        cmsWriteTag(newProfile, cmsSigGreenTRCTag, GammaTRC[1]);
        cmsWriteTag(newProfile, cmsSigBlueTRCTag, GammaTRC[2]);
        cmsWriteTag(newProfile, cmsSigProfileDescriptionTag, mlu);//desc changed

        /*  //to read XYZ values
            cmsCIEXYZ *redT = static_cast<cmsCIEXYZ*>(cmsReadTag(newProfile, cmsSigRedMatrixColumnTag));
            cmsCIEXYZ *greenT  = static_cast<cmsCIEXYZ*>(cmsReadTag(newProfile, cmsSigGreenMatrixColumnTag));
            cmsCIEXYZ *blueT  = static_cast<cmsCIEXYZ*>(cmsReadTag(newProfile, cmsSigBlueMatrixColumnTag));
            printf("rx=%f gx=%f bx=%f ry=%f gy=%f by=%f rz=%f gz=%f bz=%f\n", redT->X, greenT->X, blueT->X, redT->Y, greenT->Y, blueT->Y, redT->Z, greenT->Z, blueT->Z);
        */

        cmsMLUfree(mlu);
        cmsMLU *copyright = cmsMLUalloc(NULL, 1);
        cmsMLUsetASCII(copyright, "eng", "USA", "Copyright RawTherapee 2018, CC0");
        cmsWriteTag(newProfile, cmsSigCopyrightTag, copyright);
        cmsMLUfree(copyright);
        //cmsWriteTag(newProfile, cmsSigProfileDescriptionTag,  mlu);//desc changed
        cmsMLU *MfgDesc;
        MfgDesc   = cmsMLUalloc(NULL, 1);
        cmsMLUsetASCII(MfgDesc, "eng", "USA", "RawTherapee");
        cmsWriteTag(newProfile, cmsSigDeviceMfgDescTag, MfgDesc);
        cmsMLUfree(MfgDesc);

        /*
        Glib::ustring realoutPro;
        realoutPro = options.cacheBaseDir + "/" + fName;//ICC profile in cache
        */

        if (profileVersion == "v2" || profileVersion == "v4") {
            cmsSaveProfileToFile(newProfile,  absoluteFName.c_str());

        }

        //if (GammaTRC) {
        cmsFreeToneCurve(GammaTRC[0]);
        //}
    }
}
