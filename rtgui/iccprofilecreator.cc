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
    , appendParamsToDesc(options.ICCPC_appendParamsToDesc)
    , profileVersion(options.ICCPC_profileVersion)
    , illuminant(options.ICCPC_illuminant)
    , description(options.ICCPC_description)
    , copyright(options.ICCPC_copyright)
    , parent(rtwindow)
{

    set_default_size(600, -1);

    Gtk::Grid* mainGrid = Gtk::manage(new Gtk::Grid());
    mainGrid->set_column_spacing(3);
    mainGrid->set_row_spacing(3);

    //--------------------------------- primaries

    Gtk::Label* prilab = Gtk::manage(new Gtk::Label(M("ICCPROFCREATOR_PRIMARIES")));
    setExpandAlignProperties(prilab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    mainGrid->attach(*prilab, 0, 0, 1, 1);

    primaries = Gtk::manage(new MyComboBoxText());
    setExpandAlignProperties(primaries, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    primaries->append(M("ICCPROFCREATOR_CUSTOM"));
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

    //--------------------------------- Description

    Gtk::Label* desclab = Gtk::manage(new Gtk::Label(M("ICCPROFCREATOR_DESCRIPTION")));
    setExpandAlignProperties(desclab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_START);
    mainGrid->attach(*desclab, 0, 7, 1, 2);
    eDescription = Gtk::manage(new Gtk::Entry());
    setExpandAlignProperties(eDescription, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    eDescription->set_tooltip_text(M("ICCPROFCREATOR_DESCRIPTION_TOOLTIP"));
    mainGrid->attach(*eDescription, 1, 7, 1, 1);
    cAppendParamsToDesc = Gtk::manage(new Gtk::CheckButton(M("ICCPROFCREATOR_DESCRIPTION_ADDPARAM")));
    setExpandAlignProperties(cAppendParamsToDesc, true, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    mainGrid->attach(*cAppendParamsToDesc, 1, 8, 1, 1);

    //--------------------------------- Copyright

    Gtk::Label* copylab = Gtk::manage(new Gtk::Label(M("ICCPROFCREATOR_COPYRIGHT")));
    setExpandAlignProperties(copylab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    mainGrid->attach(*copylab, 0, 9, 1, 1);
    Gtk::Grid* copygrid = Gtk::manage(new Gtk::Grid());
    setExpandAlignProperties(copygrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    eCopyright = Gtk::manage(new Gtk::Entry());
    setExpandAlignProperties(eCopyright, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    copygrid->attach(*eCopyright, 0, 0, 1, 1);
    resetCopyright = Gtk::manage(new Gtk::Button());
    resetCopyright->add (*Gtk::manage (new RTImage ("undo-small.png", "redo-small.png")));
    setExpandAlignProperties(resetCopyright, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    resetCopyright->set_relief (Gtk::RELIEF_NONE);
    resetCopyright->set_tooltip_markup (M("ICCPROFCREATOR_COPYRIGHT_RESET_TOOLTIP"));
    resetCopyright->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
    resetCopyright->set_can_focus(false);
    copygrid->attach(*resetCopyright, 1, 0, 1, 1);
    mainGrid->attach(*copygrid, 1, 9, 1, 1);

    //--------------------------------- Adding the mainGrid

    get_content_area()->add(*mainGrid);

    //--------------------------------- Setting default values for Adjusters

    aGamma->setDefault(gamma);
    aSlope->setDefault(slope);
    aPrimariesRedX->setDefault(redPrimaryX);
    aPrimariesRedY->setDefault(redPrimaryY);
    aPrimariesGreenX->setDefault(greenPrimaryX);
    aPrimariesGreenY->setDefault(greenPrimaryY);
    aPrimariesBlueX->setDefault(bluePrimaryX);
    aPrimariesBlueY->setDefault(bluePrimaryY);

    //--------------- Updating widgets with actual values (from options)

    if (primariesPreset == "custom") {
        primaries->set_active_text(M("ICCPROFCREATOR_CUSTOM"));
    } else if (primariesPreset == "ACES-AP0") {
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
    }

    trcPresets->set_active(0);
    if (gammaPreset != "Custom") {
        trcPresets->set_active_text(gammaPreset);
    }

    aGamma->setValue(gamma);
    aSlope->setValue(slope);
    aPrimariesRedX->setValue(redPrimaryX);
    aPrimariesRedY->setValue(redPrimaryY);
    aPrimariesGreenX->setValue(greenPrimaryX);
    aPrimariesGreenY->setValue(greenPrimaryY);
    aPrimariesBlueX->setValue(bluePrimaryX);
    aPrimariesBlueY->setValue(bluePrimaryY);

    eDescription->set_text(description);
    eCopyright->set_text(copyright);
    cAppendParamsToDesc->set_active(appendParamsToDesc);

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
    resetCopyright->signal_clicked().connect(sigc::mem_fun(*this, &ICCProfileCreator::onResetCopyright));
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

void ICCProfileCreator::adjusterChanged (Adjuster* a, double newval)
{
    if (a == aPrimariesRedX   || a == aPrimariesRedY   ||
        a == aPrimariesGreenX || a == aPrimariesGreenY ||
        a == aPrimariesBlueX  || a == aPrimariesBlueY)
    {
        if (primaries->get_active_row_number() > 0) {
            ConnectionBlocker blocker(primariesconn);
            primaries->set_active(0);
            updateICCVersion();
        }
    } else if (a == aGamma || a == aSlope) {
        if (trcPresets->get_active_row_number() > 0) {
            ConnectionBlocker blocker(trcpresetsconn);
            trcPresets->set_active(0);
        }
    }
}

void ICCProfileCreator::primariesChanged()
{
    if (primaries->get_active_row_number() > 0) {
        float p[6];
        ColorTemp temp;
        Glib::ustring activeValue = primaries->get_active_text();
        Glib::ustring primPresetName = getPrimariesPresetName(activeValue);
        getPrimaries(primPresetName, p, temp);
        aPrimariesRedX->setValue(p[0]);
        aPrimariesRedY->setValue(p[1]);
        aPrimariesGreenX->setValue(p[2]);
        aPrimariesGreenY->setValue(p[3]);
        aPrimariesBlueX->setValue(p[4]);
        aPrimariesBlueY->setValue(p[5]);
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

    double gamma;
    double slope;
    getGamma(getGammaPresetName(trcPresets->get_active_text()), gamma, slope);
    aGamma->setValue(gamma);
    aSlope->setValue(slope);

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

    options.ICCPC_primariesPreset = primariesPreset = getPrimariesPresetName(primaries->get_active_text());
    options.ICCPC_gammaPreset = gammaPreset = getGammaPresetName(trcPresets->get_active_text());
    options.ICCPC_gamma = gamma = aGamma->getValue();
    options.ICCPC_slope = slope = aSlope->getValue();
    options.ICCPC_redPrimaryX = redPrimaryX = aPrimariesRedX->getValue();
    options.ICCPC_redPrimaryY = redPrimaryY = aPrimariesRedY->getValue();
    options.ICCPC_greenPrimaryX = greenPrimaryX = aPrimariesGreenX->getValue();
    options.ICCPC_greenPrimaryY = greenPrimaryY = aPrimariesGreenY->getValue();
    options.ICCPC_bluePrimaryX = bluePrimaryX = aPrimariesBlueX->getValue();
    options.ICCPC_bluePrimaryY = bluePrimaryY = aPrimariesBlueY->getValue();
    options.ICCPC_description = description = eDescription->get_text();
    options.ICCPC_copyright = copyright = eCopyright->get_text();
    options.ICCPC_appendParamsToDesc = appendParamsToDesc = cAppendParamsToDesc->get_active();
}

Glib::ustring ICCProfileCreator::getPrimariesPresetName(const Glib::ustring &preset)
{
    if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_ACESP0")) {
        return Glib::ustring("ACES-AP0");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_ACESP1")) {
        return Glib::ustring("ACES-AP1");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_ADOBE")) {
        return Glib::ustring("Adobe");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_PROPH")) {
        return Glib::ustring("ProPhoto");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_REC2020")) {
        return Glib::ustring("Rec2020");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_SRGB")) {
        return Glib::ustring("sRGB");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_WIDEG")) {
        return Glib::ustring("Widegamut");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_BEST")) {
        return Glib::ustring("BestRGB");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_BETA")) {
        return Glib::ustring("BetaRGB");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_PRIM_BRUCE")) {
        return Glib::ustring("BruceRGB");
    } else if (primaries->get_active_text() == M("ICCPROFCREATOR_CUSTOM")) {
        return Glib::ustring("custom");
    } else {
        return Glib::ustring();
    }
}

void ICCProfileCreator::getPrimaries(const Glib::ustring &preset, float *p, ColorTemp &temp)
{
    temp = ColorTemp::D50;
    if (preset == "Widegamut") {
        p[0] = 0.7350;    //Widegamut primaries
        p[1] = 0.2650;
        p[2] = 0.1150;
        p[3] = 0.8260;
        p[4] = 0.1570;
        p[5] = 0.0180;

    } else if (preset == "Adobe") {
        p[0] = 0.6400;    //Adobe primaries
        p[1] = 0.3300;
        p[2] = 0.2100;
        p[3] = 0.7100;
        p[4] = 0.1500;
        p[5] = 0.0600;
        temp = ColorTemp::D65;
    } else if (preset == "sRGB") {
        p[0] = 0.6400;    // sRGB primaries
        p[1] = 0.3300;
        p[2] = 0.3000;
        p[3] = 0.6000;
        p[4] = 0.1500;
        p[5] = 0.0600;
        temp = ColorTemp::D65;
    } else if (preset == "BruceRGB") {
        p[0] = 0.6400;    // Bruce primaries
        p[1] = 0.3300;
        p[2] = 0.2800;
        p[3] = 0.6500;
        p[4] = 0.1500;
        p[5] = 0.0600;
        temp = ColorTemp::D65;
    } else if (preset == "BetaRGB") {
        p[0] = 0.6888;    // Beta primaries
        p[1] = 0.3112;
        p[2] = 0.1986;
        p[3] = 0.7551;
        p[4] = 0.1265;
        p[5] = 0.0352;
    } else if (preset == "BestRGB") {
        p[0] = 0.7347;    // Best primaries
        p[1] = 0.2653;
        p[2] = 0.2150;
        p[3] = 0.7750;
        p[4] = 0.1300;
        p[5] = 0.0350;
    } else if (preset == "Rec2020") {
        p[0] = 0.7080;    // Rec2020 primaries
        p[1] = 0.2920;
        p[2] = 0.1700;
        p[3] = 0.7970;
        p[4] = 0.1310;
        p[5] = 0.0460;
        temp = ColorTemp::D65;
    } else if (preset == "ACES-AP0") {
        p[0] = 0.7347;    // ACES P0 primaries
        p[1] = 0.2653;
        p[2] = 0.0000;
        p[3] = 1.0;
        p[4] = 0.0001;
        p[5] = -0.0770;
        temp = ColorTemp::D60;
    } else if (preset == "ACES-AP1") {
        p[0] = 0.713;    // ACES P1 primaries
        p[1] = 0.293;
        p[2] = 0.165;
        p[3] = 0.830;
        p[4] = 0.128;
        p[5] = 0.044;
        temp = ColorTemp::D60;
    } else if (preset == "ProPhoto") {
        p[0] = 0.7347;    // ProPhoto and default primaries
        p[1] = 0.2653;
        p[2] = 0.1596;
        p[3] = 0.8404;
        p[4] = 0.0366;
        p[5] = 0.0001;
    } else if (preset == "custom") {
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
}

Glib::ustring ICCProfileCreator::getGammaPresetName(const Glib::ustring &preset)
{
    Glib::ustring name(trcPresets->get_active_text());
    if (name == M("ICCPROFCREATOR_CUSTOM")) {
        name = "Custom";
    }
    return name;
}

void ICCProfileCreator::getGamma(const Glib::ustring &preset, double &presetGamma, double &presetSlope)
{
    if (preset == "High_g1.3_s3.35") {
        presetGamma = 1.3;
        presetSlope = 3.35;
    } else if (preset == "Low_g2.6_s6.9") {
        presetGamma = 2.6;
        presetSlope = 6.9;
    } else if (preset == "sRGB_g2.4_s12.92") {
        presetGamma = 2.4;
        presetSlope = 12.92310;
    } else if (preset == "BT709_g2.2_s4.5") {
        presetGamma = 2.22;
        presetSlope = 4.5;
    } else if (preset == "linear_g1.0") {
        presetGamma = 1.;
        presetSlope = 0.;
    } else if (preset == "standard_g2.2") {
        presetGamma = 2.2;
        presetSlope = 0.;
    } else if (preset == "standard_g1.8") {
        presetGamma = 1.8;
        presetSlope = 0.;
    } else if (preset == "Lab_g3.0s9.03296") {
        presetGamma = 3.0;
        presetSlope = 9.03296;
    } else if (preset == "Custom") {
        presetGamma = gamma;
        presetSlope = slope;
    } else {
        presetGamma = 2.4;
        presetSlope = 12.92310;
    }
}

void ICCProfileCreator::onResetCopyright()
{
    eCopyright->set_text(Options::getICCProfileCopyright());
}

// Copyright (c) 2018 Jacques DESMIS <jdesmis@gmail.com>
// WARNING: the caller must lock lcmsMutex
void ICCProfileCreator::savePressed()
{
    cmsHPROFILE newProfile = nullptr;
    Glib::ustring sNewProfile;
    Glib::ustring sPrimariesPreset;
    Glib::ustring sGammaPreset;

    storeValues();

    // -------------------------------------------- Compute the default file name

    //necessary for V2 profile
    if (primariesPreset == "ACES-AP0"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.ACESp0)) {
        sNewProfile = options.rtSettings.ACESp0;
        sPrimariesPreset = "ACES-AP0";
    } else if (primariesPreset == "ACES-AP1"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.ACESp1)) {
        sNewProfile = options.rtSettings.ACESp1;
        sPrimariesPreset = "ACES-AP1";
    } else if (primariesPreset == "Adobe"      && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.adobe)) {
        sNewProfile = options.rtSettings.adobe;
        sPrimariesPreset = "Medium";
    } else if (primariesPreset == "ProPhoto"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.prophoto)) {
        sNewProfile = options.rtSettings.prophoto;
        sPrimariesPreset = "Large";
    } else if (primariesPreset == "Rec2020"    && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.rec2020)) {
        sNewProfile = options.rtSettings.rec2020;
        sPrimariesPreset = "Rec2020";
    } else if (primariesPreset == "sRGB"       && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.srgb)) {
        sNewProfile = options.rtSettings.srgb;
        sPrimariesPreset = "sRGB";
    } else if (primariesPreset == "Widegamut"  && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.widegamut)) {
        sNewProfile = options.rtSettings.widegamut;
        sPrimariesPreset = "Wide";
    } else if (primariesPreset == "BestRGB"    && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.best)) {
        sNewProfile = options.rtSettings.best;
        sPrimariesPreset = "Best";
    } else if (primariesPreset == "BetaRGB"    && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.beta)) {
        sNewProfile = options.rtSettings.beta;
        sPrimariesPreset = "Beta";
    } else if (primariesPreset == "BruceRGB"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.bruce)) {
        sNewProfile = options.rtSettings.bruce;
        sPrimariesPreset = "Bruce";
    } else if (primariesPreset == "custom") {
        sNewProfile = options.rtSettings.srgb;
        sPrimariesPreset = "Custom";
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
    getGamma(gammaPreset, presetGamma, presetSlope);
    if (gammaPreset == "High_g1.3_s3.35") {
        sGammaPreset = "High_g=1.3_s=3.35";
        ga[0] = 1.3 ;    //for high dynamic images
        ga[1] = 0.998279;
        ga[2] = 0.001721;
        ga[3] = 0.298507;
        ga[4] = 0.005746;
    } else if (gammaPreset == "Low_g2.6_s6.9") {
        sGammaPreset = "Low_g=2.6_s=6.9";
        ga[0] = 2.6 ;    //gamma 2.6 variable : for low contrast images
        ga[1] = 0.891161;
        ga[2] = 0.108839;
        ga[3] = 0.144928;
        ga[4] = 0.076332;
    } else if (gammaPreset == "sRGB_g2.4_s12.92") {
        sGammaPreset = "sRGB_g=2.4_s=12.92310";
        ga[0] = 2.40;    //sRGB 2.4 12.92  - RT default as Lightroom
        ga[1] = 0.947858;
        ga[2] = 0.052142;
        ga[3] = 0.077399;
        ga[4] = 0.039293;
    } else if (gammaPreset == "BT709_g2.2_s4.5") {
        sGammaPreset = "BT709_g=2.2_s=4.5";
        ga[0] = 2.22;    //BT709  2.2  4.5  - my preferred as D.Coffin
        ga[1] = 0.909995;
        ga[2] = 0.090005;
        ga[3] = 0.222222;
        ga[4] = 0.081071;
    } else if (gammaPreset == "linear_g1.0") {
        sGammaPreset = "Linear_g=1.0";
        ga[0] = 1.0;    //gamma=1 linear : for high dynamic images (cf D.Coffin...)
        ga[1] = 1.;
        ga[2] = 0.;
        ga[3] = 1. / eps;
        ga[4] = 0.;
    } else if (gammaPreset == "standard_g2.2") {
        sGammaPreset = "g=2.2";
        ga[0] = 2.2;    //gamma=2.2(as gamma of Adobe, Widegamut...)
        ga[1] = 1.;
        ga[2] = 0.;
        ga[3] = 1. / eps;
        ga[4] = 0.;
    } else if (gammaPreset == "standard_g1.8") {
        sGammaPreset = "g=1.8";
        ga[0] = 1.8;    //gamma=1.8(as gamma of Prophoto)
        ga[1] = 1.;
        ga[2] = 0.;
        ga[3] = 1. / eps;
        ga[4] = 0.;
    } else if (gammaPreset == "Lab_g3.0s9.03296") {
        sGammaPreset = "LAB_g3.0_s9.03296";
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

        sGammaPreset = Glib::ustring::compose("g%1_s%2",
                       Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), gamma),
                       Glib::ustring::format (std::setw(6), std::fixed, std::setprecision(5), slope));
        presetGamma = gamma;
        presetSlope = slope;
    }
    ga[5] = 0.0;
    ga[6] = 0.0;


    sPrimariesAndIlluminant = sPrimariesPreset;

    if (profileVersion == "v4" && illuminant != "DEF") {
         sPrimariesPreset += "-" + illuminant;
        //  printf("outpr=%s \n",outPr.c_str());
    }

    Glib::ustring profileDesc;
    Glib::ustring sGammaSlopeParam;//to save gamma and slope in a dmdd
    Glib::ustring sGammaSlopeDesc; //to save gamma and slope in a desc
    Glib::ustring sGamma;
    Glib::ustring sSlope;

    if (gammaPreset == "Custom") {
        sGamma = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), gamma);
        sSlope = Glib::ustring::format (std::setw(6), std::fixed, std::setprecision(5), slope);
        fName = Glib::ustring::compose("RT%1_%2_g%3_s%4.icc", profileVersion, sPrimariesAndIlluminant, sGamma, sSlope);
        profileDesc = sPrimariesPreset;
    } else {
        sGamma = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), presetGamma);
        sSlope = Glib::ustring::format (std::setw(6), std::fixed, std::setprecision(5), presetSlope);
        fName = Glib::ustring::compose("RT%1_%2_%3.icc", profileVersion, sPrimariesAndIlluminant, sGammaPreset);
        profileDesc == sPrimariesPreset + sGammaPreset;
    }
    sGammaSlopeParam = Glib::ustring::compose("g%1s%2!", sGamma, sSlope);
    sGammaSlopeDesc = Glib::ustring::compose("g=%1 s=%2", sGamma, sSlope);

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

    // --------------- main tags ------------------

    if (profileVersion == "v4") {
        cmsSetProfileVersion(newProfile, 4.3);
    } else {
        cmsSetProfileVersion(newProfile, 2.0);
    }

//change
    float p[6]; //primaries
    ga[6] = 0.0;

    ColorTemp temp;
    getPrimaries(primariesPreset, p, temp);

    cmsCIExyY xyD;
    cmsCIExyYTRIPLE Primaries = {
        {p[0], p[1], 1.0}, // red
        {p[2], p[3], 1.0}, // green
        {p[4], p[5], 1.0}  // blue
    };

    if (profileVersion == "v4" && illuminant != "DEF") {
        double tempv4 = 5000.;
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
        cmsWhitePointFromTemp(&xyD, tempv4);
    } else {
        cmsWhitePointFromTemp(&xyD, (double)temp);
    }

    if (illuminant == "stdA") {
        xyD = {0.447573, 0.407440, 1.0};
    }

    // Calculate output profile's rTRC gTRC bTRC
    cmsToneCurve* GammaTRC[3];
    GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(nullptr, 5, ga);

    if (profileVersion == "v4") {
        newProfile = cmsCreateRGBProfile(&xyD, &Primaries, GammaTRC);
    }

    cmsWriteTag(newProfile, cmsSigRedTRCTag, GammaTRC[0]);
    cmsWriteTag(newProfile, cmsSigGreenTRCTag, GammaTRC[1]);
    cmsWriteTag(newProfile, cmsSigBlueTRCTag, GammaTRC[2]);

    // --------------- set dmnd tag ------------------

    cmsMLU *dmnd;
    dmnd = cmsMLUalloc(nullptr, 1);
    cmsMLUsetASCII(dmnd, "en", "US", "RawTherapee");
    cmsWriteTag(newProfile, cmsSigDeviceMfgDescTag, dmnd);
    cmsMLUfree(dmnd);

    // --------------- set dmdd tag ------------------

    if (profileVersion == "v2") {
        //write in tag 'dmdd' values of current gamma and slope to retrieve after in Output profile
        std::wostringstream wGammaSlopeParam;
        wGammaSlopeParam << sGammaSlopeParam;

        cmsMLU *dmdd = cmsMLUalloc(nullptr, 1);
        // Language code (2 letters code) : https://www.iso.org/obp/ui/
        // Country code (2 letters code)  : http://www.loc.gov/standards/iso639-2/php/code_list.php
        if (sGammaSlopeParam.is_ascii()) {
            if (cmsMLUsetASCII(dmdd, "en", "US", sGammaSlopeParam.c_str())) {
                if (!cmsWriteTag(newProfile, cmsSigProfileDescriptionTag, dmdd)) {
                    printf("Error: Can't write cmsSigProfileDescriptionTag!\n");
                }
            }
        } else if (cmsMLUsetWide(dmdd, "en", "US", wGammaSlopeParam.str().c_str())) {
            if (!cmsWriteTag(newProfile, cmsSigDeviceModelDescTag, dmdd)) {
                printf("Error: Can't write cmsSigDeviceModelDescTag!\n");
            }
        } else {
            printf("Error: cmsMLUsetWide failed for dmdd \"%s\" !\n", sGammaSlopeParam.c_str());
        }
        cmsMLUfree(dmdd);
    }

    // --------------- set desc tag ------------------

    Glib::ustring sDescription;
    if (!description.empty()) {
        if (cAppendParamsToDesc->get_active()) {
            sDescription = description + " / " + sGammaSlopeDesc;
        } else {
            sDescription = description;
        }
    } else {
        if (cAppendParamsToDesc->get_active()) {
            sDescription = profileDesc + " / " + sGammaSlopeDesc;
        } else {
            sDescription = profileDesc;
        }
    }

    //write in tag 'dmdd' values of current gamma and slope to retrieve after in Output profile
    std::wostringstream wDescription;
    wDescription << sDescription;

    cmsMLU *descMLU = cmsMLUalloc(nullptr, 1);
    // Language code (2 letters code) : https://www.iso.org/obp/ui/
    // Country code (2 letters code)  : http://www.loc.gov/standards/iso639-2/php/code_list.php
    if (sDescription.is_ascii()) {
        if (cmsMLUsetASCII(descMLU, "en", "US", sDescription.c_str())) {
            if (!cmsWriteTag(newProfile, cmsSigProfileDescriptionTag, descMLU)) {
                printf("Error: Can't write cmsSigProfileDescriptionTag!\n");
            }
        }
    } else if (cmsMLUsetWide(descMLU, "en", "US", wDescription.str().c_str())) {
        if (!cmsWriteTag(newProfile, cmsSigProfileDescriptionTag, descMLU)) {
            printf("Error: Can't write cmsSigProfileDescriptionTag!\n");
        }
    } else {
        printf("Error: cmsMLUsetWide failed for desc \"%s\" !\n", sDescription.c_str());
    }
    cmsMLUfree(descMLU);

    // --------------- set cprt tag ------------------

    std::wostringstream wCopyright;
    wCopyright << copyright;

    cmsMLU *copyMLU = cmsMLUalloc(nullptr, 1);
    if (cmsMLUsetWide(copyMLU, "en", "US", wCopyright.str().c_str())) {
        if (!cmsWriteTag(newProfile, cmsSigCopyrightTag, copyMLU)) {
            printf("Error: Can't write cmsSigCopyrightTag!\n");
        }
    } else {
        printf("Error: cmsMLUsetWide failed for cprt \"%s\" !\n", copyright.c_str());
    }
    cmsMLUfree(copyMLU);


    /*  //to read XYZ values
        cmsCIEXYZ *redT = static_cast<cmsCIEXYZ*>(cmsReadTag(newProfile, cmsSigRedMatrixColumnTag));
        cmsCIEXYZ *greenT  = static_cast<cmsCIEXYZ*>(cmsReadTag(newProfile, cmsSigGreenMatrixColumnTag));
        cmsCIEXYZ *blueT  = static_cast<cmsCIEXYZ*>(cmsReadTag(newProfile, cmsSigBlueMatrixColumnTag));
        printf("rx=%f gx=%f bx=%f ry=%f gy=%f by=%f rz=%f gz=%f bz=%f\n", redT->X, greenT->X, blueT->X, redT->Y, greenT->Y, blueT->Y, redT->Z, greenT->Z, blueT->Z);
    */

    cmsSaveProfileToFile(newProfile,  absoluteFName.c_str());

    cmsFreeToneCurve(GammaTRC[0]);
}
