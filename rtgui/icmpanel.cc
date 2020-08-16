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
#include <iomanip>

#include "icmpanel.h"

#include "eventmapper.h"
#include "guiutils.h"
#include "options.h"
#include "pathutils.h"
#include "rtimage.h"

#include "../rtengine/dcp.h"
#include "../rtengine/iccstore.h"
#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

ICMPanel::ICMPanel() : FoldableToolPanel(this, "icm", M("TP_ICM_LABEL")), iunchanged(nullptr), icmplistener(nullptr)
{
    auto m = ProcEventMapper::getInstance();
    EvICMprimariMethod = m->newEvent(GAMMA, "HISTORY_MSG_ICM_OUTPUT_PRIMARIES");
    EvICMprofileMethod = m->newEvent(GAMMA, "HISTORY_MSG_ICM_OUTPUT_TYPE");
    EvICMtempMethod = m->newEvent(GAMMA, "HISTORY_MSG_ICM_OUTPUT_TEMP");
    EvICMpredx = m->newEvent(GAMMA, "HISTORY_MSG_ICMPREDX");
    EvICMpredy = m->newEvent(GAMMA, "HISTORY_MSG_ICMPREDY");
    EvICMpgrex = m->newEvent(GAMMA, "HISTORY_MSG_ICMPGREX");
    EvICMpgrey = m->newEvent(GAMMA, "HISTORY_MSG_ICMPGREY");
    EvICMpblux = m->newEvent(GAMMA, "HISTORY_MSG_ICMPBLUX");
    EvICMpbluy = m->newEvent(GAMMA, "HISTORY_MSG_ICMPBLUY");
    EvICMgamm = m->newEvent(AUTOEXP, "HISTORY_MSG_ICM_WORKING_GAMMA");
    EvICMslop = m->newEvent(AUTOEXP, "HISTORY_MSG_ICM_WORKING_SLOPE");
    EvICMtrcinMethod = m->newEvent(AUTOEXP, "HISTORY_MSG_ICM_WORKING_TRC_METHOD");

    isBatchMode = lastToneCurve = lastApplyLookTable = lastApplyBaselineExposureOffset = lastApplyHueSatMap = false;

    ipDialog = Gtk::manage(new MyFileChooserButton(M("TP_ICM_INPUTDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    ipDialog->set_tooltip_text(M("TP_ICM_INPUTCUSTOM_TOOLTIP"));
    bindCurrentFolder(*ipDialog, options.lastIccDir);


    // ------------------------------- Input profile


    Gtk::Frame *iFrame = Gtk::manage(new Gtk::Frame(M("TP_ICM_INPUTPROFILE")));
    iFrame->set_label_align(0.025, 0.5);

    iVBox = Gtk::manage(new Gtk::VBox());

    inone = Gtk::manage(new Gtk::RadioButton(M("TP_ICM_INPUTNONE")));
    inone->set_tooltip_text(M("TP_ICM_INPUTNONE_TOOLTIP"));
    iVBox->pack_start(*inone, Gtk::PACK_SHRINK);

    iembedded = Gtk::manage(new Gtk::RadioButton(M("TP_ICM_INPUTEMBEDDED")));
    iembedded->set_tooltip_text(M("TP_ICM_INPUTEMBEDDED_TOOLTIP"));
    iVBox->pack_start(*iembedded, Gtk::PACK_SHRINK);

    icamera = Gtk::manage(new Gtk::RadioButton(M("TP_ICM_INPUTCAMERA")));
    icamera->set_tooltip_text(M("TP_ICM_INPUTCAMERA_TOOLTIP"));
    iVBox->pack_start(*icamera, Gtk::PACK_SHRINK);

    icameraICC = Gtk::manage(new Gtk::RadioButton(M("TP_ICM_INPUTCAMERAICC")));
    icameraICC->set_tooltip_text(M("TP_ICM_INPUTCAMERAICC_TOOLTIP"));
    iVBox->pack_start(*icameraICC, Gtk::PACK_SHRINK);

    ifromfile = Gtk::manage(new Gtk::RadioButton(M("TP_ICM_INPUTCUSTOM") + ":"));
    Gtk::HBox* ffbox = Gtk::manage(new Gtk::HBox());
    ifromfile->set_tooltip_text(M("TP_ICM_INPUTCUSTOM_TOOLTIP"));
    ffbox->pack_start(*ifromfile, Gtk::PACK_SHRINK);
    ffbox->pack_start(*ipDialog);

    iVBox->pack_start(*ffbox, Gtk::PACK_SHRINK);

    opts = icamera->get_group();
    icameraICC->set_group(opts);
    iembedded->set_group(opts);
    ifromfile->set_group(opts);
    inone->set_group(opts);

    dcpFrame = Gtk::manage(new Gtk::Frame("DCP"));

    Gtk::Grid* dcpGrid = Gtk::manage(new Gtk::Grid());
    dcpGrid->set_column_homogeneous(false);
    dcpGrid->set_row_homogeneous(false);
    dcpGrid->set_column_spacing(2);
    dcpGrid->set_row_spacing(2);

    Gtk::Grid* dcpIllGrid = Gtk::manage(new Gtk::Grid());
    dcpIllGrid->set_column_homogeneous(false);
    dcpIllGrid->set_row_homogeneous(false);
    dcpIllGrid->set_column_spacing(2);
    dcpIllGrid->set_row_spacing(2);

    dcpIllLabel = Gtk::manage(new Gtk::Label(M("TP_ICM_DCPILLUMINANT") + ":"));
    setExpandAlignProperties(dcpIllLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    dcpIllLabel->set_tooltip_text(M("TP_ICM_DCPILLUMINANT_TOOLTIP"));
    dcpIllLabel->show();
    dcpIll = Gtk::manage(new MyComboBoxText());
    setExpandAlignProperties(dcpIll, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    dcpIll->set_tooltip_text(M("TP_ICM_DCPILLUMINANT_TOOLTIP"));
    dcpIll->append(M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
    dcpIll->append(M("TP_ICM_DCPILLUMINANT") + " 1");
    dcpIll->append(M("TP_ICM_DCPILLUMINANT") + " 2");
    dcpIll->show();
    dcpTemperatures[0] = 0;
    dcpTemperatures[1] = 0;
    dcpIllGrid->attach_next_to(*dcpIllLabel, Gtk::POS_LEFT, 1, 1);
    dcpIllGrid->attach_next_to(*dcpIll, *dcpIllLabel, Gtk::POS_RIGHT, 1, 1);

    ckbToneCurve = Gtk::manage(new Gtk::CheckButton(M("TP_ICM_TONECURVE")));
    ckbToneCurve->set_sensitive(false);
    ckbToneCurve->set_tooltip_text(M("TP_ICM_TONECURVE_TOOLTIP"));
    setExpandAlignProperties(ckbToneCurve, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    ckbApplyLookTable = Gtk::manage(new Gtk::CheckButton(M("TP_ICM_APPLYLOOKTABLE")));
    ckbApplyLookTable->set_sensitive(false);
    ckbApplyLookTable->set_tooltip_text(M("TP_ICM_APPLYLOOKTABLE_TOOLTIP"));
    setExpandAlignProperties(ckbApplyLookTable, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    ckbApplyHueSatMap = Gtk::manage(new Gtk::CheckButton(M("TP_ICM_APPLYHUESATMAP")));
    ckbApplyHueSatMap->set_sensitive(false);
    ckbApplyHueSatMap->set_tooltip_text(M("TP_ICM_APPLYHUESATMAP_TOOLTIP"));
    setExpandAlignProperties(ckbApplyHueSatMap, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    ckbApplyBaselineExposureOffset = Gtk::manage(new Gtk::CheckButton(M("TP_ICM_APPLYBASELINEEXPOSUREOFFSET")));
    ckbApplyBaselineExposureOffset->set_sensitive(false);
    ckbApplyBaselineExposureOffset->set_tooltip_text(M("TP_ICM_APPLYBASELINEEXPOSUREOFFSET_TOOLTIP"));
    setExpandAlignProperties(ckbApplyBaselineExposureOffset, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    dcpGrid->attach_next_to(*dcpIllGrid, Gtk::POS_BOTTOM, 1, 1);
    dcpGrid->attach_next_to(*ckbToneCurve, Gtk::POS_BOTTOM, 1, 1);
    dcpGrid->attach_next_to(*ckbApplyHueSatMap, Gtk::POS_BOTTOM, 1, 1);
    dcpGrid->attach_next_to(*ckbApplyLookTable, Gtk::POS_BOTTOM, 1, 1);
    dcpGrid->attach_next_to(*ckbApplyBaselineExposureOffset, Gtk::POS_BOTTOM, 1, 1);

    dcpFrame->add(*dcpGrid);
    dcpFrame->set_sensitive(false);
    iVBox->pack_start(*dcpFrame);

    saveRef = Gtk::manage(new Gtk::Button(M("TP_ICM_SAVEREFERENCE")));
    saveRef->set_image(*Gtk::manage(new RTImage("save-small.png")));
    saveRef->set_alignment(0.5f, 0.5f);
    saveRef->set_tooltip_markup(M("TP_ICM_SAVEREFERENCE_TOOLTIP"));
    iVBox->pack_start(*saveRef, Gtk::PACK_SHRINK);

    iFrame->add(*iVBox);
    pack_start(*iFrame, Gtk::PACK_EXPAND_WIDGET);


    // ---------------------------- Working profile


    Gtk::Frame *wFrame = Gtk::manage(new Gtk::Frame(M("TP_ICM_WORKINGPROFILE")));
    wFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *wProfVBox = Gtk::manage(new Gtk::VBox());

    wProfNames = Gtk::manage(new MyComboBoxText());
    wProfVBox->pack_start(*wProfNames, Gtk::PACK_SHRINK);

    std::vector<Glib::ustring> wpnames = rtengine::ICCStore::getInstance()->getWorkingProfiles();

    for (size_t i = 0; i < wpnames.size(); i++) {
        wProfNames->append(wpnames[i]);
    }

    wProfNames->set_active(0);

    // wFrame->add(*wVBox);

    //-----------------gamma TRC working

    wTRCHBox = Gtk::manage(new Gtk::HBox());

    Gtk::Label* wtrclab = Gtk::manage(new Gtk::Label(M("TP_ICM_WORKING_TRC")));

    wTRCHBox->pack_start(*wtrclab, Gtk::PACK_SHRINK);
    wTRC = Gtk::manage(new MyComboBoxText());
    wTRCHBox->pack_start(*wTRC, Gtk::PACK_EXPAND_WIDGET);
    wProfVBox->pack_start(*wTRCHBox, Gtk::PACK_EXPAND_WIDGET);
    wTRC->append(M("TP_ICM_WORKING_TRC_NONE"));
    wTRC->append(M("TP_ICM_WORKING_TRC_CUSTOM"));

    wTRC->set_active(0);
    wTRC->set_tooltip_text(M("TP_ICM_WORKING_TRC_TOOLTIP"));


    wGamma = Gtk::manage(new Adjuster(M("TP_ICM_WORKING_TRC_GAMMA"), 0.40, 15.0, 0.001, 2.4));
    wSlope = Gtk::manage(new Adjuster(M("TP_ICM_WORKING_TRC_SLOPE"), 0., 150., 0.01, 12.92310));
    wProfVBox->pack_start(*wGamma, Gtk::PACK_SHRINK);
    wGamma->show();

    wProfVBox->pack_start(*wSlope, Gtk::PACK_SHRINK);
    wSlope->show();


    wGamma->setAdjusterListener(this);
    wSlope->setAdjusterListener(this);

    wGamma->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    wSlope->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    wFrame->add(*wProfVBox);

    pack_start(*wFrame, Gtk::PACK_EXPAND_WIDGET);


    // ---------------------------- Output profile


    Gtk::Frame *oFrame = Gtk::manage(new Gtk::Frame(M("TP_ICM_OUTPUTPROFILE")));
    oFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *oProfVBox = Gtk::manage(new Gtk::VBox());

    oProfNames = Gtk::manage(new MyComboBoxText());
    oProfVBox->pack_start(*oProfNames, Gtk::PACK_SHRINK);

    oProfNames->append(M("TP_ICM_NOICM"));
    oProfNames->set_active(0);

    std::vector<Glib::ustring> opnames = ICCStore::getInstance()->getProfiles(rtengine::ICCStore::ProfileType::OUTPUT);

    for (size_t i = 0; i < opnames.size(); i++) {
        oProfNames->append(opnames[i]);
    }

    oProfNames->set_active(0);

    // Rendering intent
    Gtk::HBox *riHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* outputIntentLbl = Gtk::manage(new Gtk::Label(M("TP_ICM_PROFILEINTENT")));
    riHBox->pack_start(*outputIntentLbl, Gtk::PACK_SHRINK);
    oRendIntent.reset(new PopUpButton());
    oRendIntent->addEntry("intent-perceptual.png", M("PREFERENCES_INTENT_PERCEPTUAL"));
    oRendIntent->addEntry("intent-relative.png", M("PREFERENCES_INTENT_RELATIVE"));
    oRendIntent->addEntry("intent-saturation.png", M("PREFERENCES_INTENT_SATURATION"));
    oRendIntent->addEntry("intent-absolute.png", M("PREFERENCES_INTENT_ABSOLUTE"));
    oRendIntent->setSelected(1);
    oRendIntent->show();
    riHBox->pack_start(*oRendIntent->buttonGroup, Gtk::PACK_EXPAND_PADDING);
    oProfVBox->pack_start(*riHBox, Gtk::PACK_SHRINK);

    // Black Point Compensation
    obpc = Gtk::manage(new Gtk::CheckButton((M("TP_ICM_BPC"))));
    obpc->set_active(true);
    oProfVBox->pack_start(*obpc, Gtk::PACK_SHRINK);

    oFrame->add(*oProfVBox);

    pack_start(*oFrame, Gtk::PACK_EXPAND_WIDGET);

    // ---------------------------- Output gamma list entries

    Glib::RefPtr<Gtk::FileFilter> filter_icc = Gtk::FileFilter::create();
    filter_icc->set_name(M("FILECHOOSER_FILTER_COLPROF"));
    filter_icc->add_pattern("*.dcp");
    filter_icc->add_pattern("*.DCP");
    filter_icc->add_pattern("*.icc");
    filter_icc->add_pattern("*.icm");
    filter_icc->add_pattern("*.ICC");
    filter_icc->add_pattern("*.ICM");
    Glib::RefPtr<Gtk::FileFilter> filter_iccdng = Gtk::FileFilter::create();
    filter_iccdng->set_name(M("FILECHOOSER_FILTER_COLPROF") + " + DNG");
    filter_iccdng->add_pattern("*.dcp");
    filter_iccdng->add_pattern("*.DCP");
    filter_iccdng->add_pattern("*.dng");
    filter_iccdng->add_pattern("*.DNG");
    filter_iccdng->add_pattern("*.icc");
    filter_iccdng->add_pattern("*.icm");
    filter_iccdng->add_pattern("*.ICC");
    filter_iccdng->add_pattern("*.ICM");
    Glib::RefPtr<Gtk::FileFilter> filter_any = Gtk::FileFilter::create();
    filter_any->set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any->add_pattern("*");

    ipDialog->add_filter(filter_icc);
    ipDialog->add_filter(filter_iccdng);
    ipDialog->add_filter(filter_any);
#ifdef WIN32
    ipDialog->set_show_hidden(true);  // ProgramData is hidden on Windows
#endif

    wprofnamesconn = wProfNames->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::wpChanged));
    oprofnamesconn = oProfNames->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::opChanged));
    orendintentconn = oRendIntent->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::oiChanged));
    dcpillconn = dcpIll->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::dcpIlluminantChanged));
    wtrcconn = wTRC->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::wtrcinChanged));

    obpcconn = obpc->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::oBPCChanged));
    tcurveconn = ckbToneCurve->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::toneCurveChanged));
    ltableconn = ckbApplyLookTable->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::applyLookTableChanged));
    beoconn = ckbApplyBaselineExposureOffset->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::applyBaselineExposureOffsetChanged));
    hsmconn = ckbApplyHueSatMap->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::applyHueSatMapChanged));

    icamera->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::ipChanged));
    icameraICC->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::ipChanged));
    iembedded->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::ipChanged));
    ifromfile->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::ipChanged));

    ipc = ipDialog->signal_selection_changed().connect(sigc::mem_fun(*this, &ICMPanel::ipSelectionChanged));
    saveRef->signal_pressed().connect(sigc::mem_fun(*this, &ICMPanel::saveReferencePressed));

    show_all();
}

void ICMPanel::updateRenderingIntent(const Glib::ustring &profile)
{
    const uint8_t supportedIntents = rtengine::ICCStore::getInstance()->getOutputIntents(profile);
    const bool supportsPerceptual = supportedIntents & 1 << INTENT_PERCEPTUAL;
    const bool supportsRelative   = supportedIntents & 1 << INTENT_RELATIVE_COLORIMETRIC;
    const bool supportsSaturation = supportedIntents & 1 << INTENT_SATURATION;
    const bool supportsAbsolute   = supportedIntents & 1 << INTENT_ABSOLUTE_COLORIMETRIC;

    //printf("Intents: %d / Perceptual: %d  Relative: %d  Saturation: %d  Absolute: %d\n", supportedIntents, supportsPerceptual, supportsRelative, supportsSaturation, supportsAbsolute);

    if (!profile.empty() && (supportsPerceptual || supportsRelative || supportsSaturation || supportsAbsolute)) {
        oRendIntent->set_sensitive(true);
        oRendIntent->setItemSensitivity(0, supportsPerceptual);
        oRendIntent->setItemSensitivity(1, supportsRelative);
        oRendIntent->setItemSensitivity(2, supportsSaturation);
        oRendIntent->setItemSensitivity(3, supportsAbsolute);
    } else {
        oRendIntent->setItemSensitivity(0, true);
        oRendIntent->setItemSensitivity(1, true);
        oRendIntent->setItemSensitivity(2, true);
        oRendIntent->setItemSensitivity(3, true);
        oRendIntent->set_sensitive(false);
        oRendIntent->setSelected(1);
    }
}

void ICMPanel::updateDCP(int dcpIlluminant, Glib::ustring dcp_name)
{
    ConnectionBlocker dcpillconn_(dcpillconn);

    if (isBatchMode) {
        dcpFrame->set_sensitive(true);
        ckbToneCurve->set_sensitive(true);
        ckbApplyLookTable->set_sensitive(true);
        ckbApplyBaselineExposureOffset->set_sensitive(true);
        ckbApplyHueSatMap->set_sensitive(true);
        dcpIllLabel->set_sensitive(true);
        dcpIll->set_sensitive(true);

        if (dcpTemperatures[0] != 0 || dcpTemperatures[1] != 0) {
            int curr_active = dcpIll->get_active_row_number();
            dcpIll->remove_all();
            dcpIll->append(M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
            dcpIll->append(M("TP_ICM_DCPILLUMINANT") + " 1");
            dcpIll->append(M("TP_ICM_DCPILLUMINANT") + " 2");
            dcpIll->append(M("GENERAL_UNCHANGED"));
            dcpTemperatures[0] = 0;
            dcpTemperatures[1] = 0;
            dcpIll->set_active(curr_active);
        }

        if (dcpIll->get_active_row_number() == -1 && dcpIlluminant == -1) {
            dcpIll->set_active(0);
        } else if (dcpIlluminant >= 0 && dcpIlluminant != dcpIll->get_active_row_number()) {
            dcpIll->set_active(dcpIlluminant);
        }

        dcpIll->set_sensitive(true);
        dcpIllLabel->set_sensitive(true);
        return;
    }

    ckbToneCurve->set_sensitive(false);
    ckbApplyLookTable->set_sensitive(false);
    ckbApplyBaselineExposureOffset->set_sensitive(false);
    ckbApplyHueSatMap->set_sensitive(false);
    dcpIllLabel->set_sensitive(false);
    dcpIll->set_sensitive(false);
    dcpFrame->set_sensitive(false);

    DCPProfile* dcp = nullptr;

    if (dcp_name == "(cameraICC)") {
        dcp = DCPStore::getInstance()->getStdProfile(camName);
    } else if (ifromfile->get_active() && DCPStore::getInstance()->isValidDCPFileName(dcp_name)) {
        dcp = DCPStore::getInstance()->getProfile(dcp_name);
    }

    if (dcp) {
        dcpFrame->set_sensitive(true);

        if (dcp->getHasToneCurve()) {
            ckbToneCurve->set_sensitive(true);
        }

        if (dcp->getHasLookTable()) {
            ckbApplyLookTable->set_sensitive(true);
        }

        if (dcp->getHasBaselineExposureOffset()) {
            ckbApplyBaselineExposureOffset->set_sensitive(true);
        }

        if (dcp->getHasHueSatMap()) {
            ckbApplyHueSatMap->set_sensitive(true);
        }

        const DCPProfile::Illuminants illuminants = dcp->getIlluminants();

        if (illuminants.will_interpolate) {
            if (dcpTemperatures[0] != illuminants.temperature_1 || dcpTemperatures[1] != illuminants.temperature_2) {
                char tempstr1[64], tempstr2[64];
                sprintf(tempstr1, "%.0fK", illuminants.temperature_1);
                sprintf(tempstr2, "%.0fK", illuminants.temperature_2);
                int curr_active = dcpIll->get_active_row_number();
                dcpIll->remove_all();
                dcpIll->append(M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
                dcpIll->append(tempstr1);
                dcpIll->append(tempstr2);
                dcpTemperatures[0] = illuminants.temperature_1;
                dcpTemperatures[1] = illuminants.temperature_2;
                dcpIll->set_active(curr_active);
            }

            if (dcpIlluminant > 2) {
                dcpIlluminant = 0;
            }

            if (dcpIll->get_active_row_number() == -1 && dcpIlluminant == -1) {
                dcpIll->set_active(0);
            } else if (dcpIlluminant >= 0 && dcpIlluminant != dcpIll->get_active_row_number()) {
                dcpIll->set_active(dcpIlluminant);
            }

            dcpIll->set_sensitive(true);
            dcpIllLabel->set_sensitive(true);
        } else {
            if (dcpIll->get_active_row_number() != -1) {
                dcpIll->set_active(-1);
            }
        }
    }

    if (!dcpIllLabel->get_sensitive() && dcpIll->get_active_row_number() != 0) {
        if (dcpTemperatures[0] != 0 || dcpTemperatures[1] != 0) {
            int curr_active = dcpIll->get_active_row_number();
            dcpIll->remove_all();
            dcpIll->append(M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
            dcpIll->append(M("TP_ICM_DCPILLUMINANT") + " 1");
            dcpIll->append(M("TP_ICM_DCPILLUMINANT") + " 2");

            if (isBatchMode) {
                dcpIll->append(M("GENERAL_UNCHANGED"));
            }

            dcpTemperatures[0] = 0;
            dcpTemperatures[1] = 0;
            dcpIll->set_active(curr_active);
        }
    }
}

void ICMPanel::read(const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener();

    ConnectionBlocker obpcconn_(obpcconn);
    ConnectionBlocker ipc_(ipc);
    ConnectionBlocker tcurveconn_(tcurveconn);
    ConnectionBlocker ltableconn_(ltableconn);
    ConnectionBlocker beoconn_(beoconn);
    ConnectionBlocker hsmconn_(hsmconn);
    ConnectionBlocker wprofnamesconn_(wprofnamesconn);
    ConnectionBlocker oprofnamesconn_(oprofnamesconn);
    ConnectionBlocker orendintentconn_(orendintentconn);
    ConnectionBlocker dcpillconn_(dcpillconn);

    if (pp->icm.inputProfile.substr(0, 5) != "file:" && !ipDialog->get_filename().empty()) {
        ipDialog->set_filename(pp->icm.inputProfile);
    }

    if (pp->icm.inputProfile == "(none)") {
        inone->set_active(true);
        updateDCP(pp->icm.dcpIlluminant, "");
    } else if (pp->icm.inputProfile == "(embedded)" || ((pp->icm.inputProfile == "(camera)" || pp->icm.inputProfile.empty()) && icamera->get_state() == Gtk::STATE_INSENSITIVE)) {
        iembedded->set_active(true);
        updateDCP(pp->icm.dcpIlluminant, "");
    } else if ((pp->icm.inputProfile == "(cameraICC)") && icameraICC->get_state() != Gtk::STATE_INSENSITIVE) {
        icameraICC->set_active(true);
        updateDCP(pp->icm.dcpIlluminant, "(cameraICC)");
    } else if ((pp->icm.inputProfile == "(cameraICC)") && icamera->get_state() != Gtk::STATE_INSENSITIVE && icameraICC->get_state() == Gtk::STATE_INSENSITIVE) {
        // this is the case when (cameraICC) is instructed by packaged profiles, but ICC file is not found
        // therefore falling back UI to explicitly reflect the (camera) option
        icamera->set_active(true);
        updateDCP(pp->icm.dcpIlluminant, "");
    } else if ((pp->icm.inputProfile == "(cameraICC)") && icamera->get_state() == Gtk::STATE_INSENSITIVE && icameraICC->get_state() == Gtk::STATE_INSENSITIVE) {
        // If neither (camera) nor (cameraICC) are available, as is the case when loading a non-raw, activate (embedded).
        iembedded->set_active(true);
        updateDCP(pp->icm.dcpIlluminant, "(cameraICC)");
    } else if ((pp->icm.inputProfile == "(camera)" || pp->icm.inputProfile.empty()) && icamera->get_state() != Gtk::STATE_INSENSITIVE) {
        icamera->set_active(true);
        updateDCP(pp->icm.dcpIlluminant, "");
    } else {
        ifromfile->set_active(true);
        oldip = pp->icm.inputProfile.substr(5);  // cut of "file:"
        ipDialog->set_filename(pp->icm.inputProfile.substr(5));
        updateDCP(pp->icm.dcpIlluminant, pp->icm.inputProfile.substr(5));
    }

    wProfNames->set_active_text(pp->icm.workingProfile);

    if (pp->icm.workingTRC == "none") {
        wTRC->set_active(0);
    } else if (pp->icm.workingTRC == "Custom") {
        wTRC->set_active(1);
    }

    wtrcinChanged();

    if (pp->icm.outputProfile == ColorManagementParams::NoICMString) {
        oProfNames->set_active_text(M("TP_ICM_NOICM"));
    } else {
        oProfNames->set_active_text(pp->icm.outputProfile);
    }

    if (oProfNames->get_active_row_number() == -1) {
        oProfNames->set_active_text(M("TP_ICM_NOICM"));
    }

    oRendIntent->setSelected(pp->icm.outputIntent);

    obpc->set_active(pp->icm.outputBPC);
    ckbToneCurve->set_active(pp->icm.toneCurve);
    lastToneCurve = pp->icm.toneCurve;
    ckbApplyLookTable->set_active(pp->icm.applyLookTable);
    lastApplyLookTable = pp->icm.applyLookTable;
    ckbApplyBaselineExposureOffset->set_active(pp->icm.applyBaselineExposureOffset);
    lastApplyBaselineExposureOffset = pp->icm.applyBaselineExposureOffset;
    ckbApplyHueSatMap->set_active(pp->icm.applyHueSatMap);
    lastApplyHueSatMap = pp->icm.applyHueSatMap;

    wGamma->setValue(pp->icm.workingTRCGamma);
    wSlope->setValue(pp->icm.workingTRCSlope);

    if (pedited) {
        iunchanged->set_active(!pedited->icm.inputProfile);
        obpc->set_inconsistent(!pedited->icm.outputBPC);
        ckbToneCurve->set_inconsistent(!pedited->icm.toneCurve);
        ckbApplyLookTable->set_inconsistent(!pedited->icm.applyLookTable);
        ckbApplyBaselineExposureOffset->set_inconsistent(!pedited->icm.applyBaselineExposureOffset);
        ckbApplyHueSatMap->set_inconsistent(!pedited->icm.applyHueSatMap);

        if (!pedited->icm.workingProfile) {
            wProfNames->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.outputProfile) {
            oProfNames->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.outputIntent) {
            oRendIntent->setSelected(4);
        }

        if (!pedited->icm.dcpIlluminant) {
            dcpIll->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.workingTRC) {
            wTRC->set_active_text(M("GENERAL_UNCHANGED"));
        }

        wGamma->setEditedState(pedited->icm.workingTRCGamma ? Edited : UnEdited);
        wSlope->setEditedState(pedited->icm.workingTRCSlope  ? Edited : UnEdited);

    }

    enableListener();
}

void ICMPanel::write(ProcParams* pp, ParamsEdited* pedited)
{

    if (inone->get_active()) {
        pp->icm.inputProfile = "(none)";
    } else if (iembedded->get_active()) {
        pp->icm.inputProfile = "(embedded)";
    } else if (icamera->get_active()) {
        pp->icm.inputProfile = "(camera)";
    } else if (icameraICC->get_active()) {
        pp->icm.inputProfile = "(cameraICC)";
    } else {
        if (Glib::file_test(ipDialog->get_filename(), Glib::FILE_TEST_EXISTS) && !Glib::file_test(ipDialog->get_filename(), Glib::FILE_TEST_IS_DIR)) {
            pp->icm.inputProfile = "file:" + ipDialog->get_filename();
        } else {
            pp->icm.inputProfile = "";    // just a directory
        }
    }

    pp->icm.workingProfile = wProfNames->get_active_text();
    pp->icm.dcpIlluminant = rtengine::max<int>(dcpIll->get_active_row_number(), 0);
    pp->icm.workingTRC = wTRC->get_active_text();

    if (oProfNames->get_active_text() == M("TP_ICM_NOICM")) {
        pp->icm.outputProfile  = ColorManagementParams::NoICMString;
    } else {
        pp->icm.outputProfile  = oProfNames->get_active_text();
    }

    int ointentVal = oRendIntent->getSelected();

    if (ointentVal >= 0 && ointentVal < RI__COUNT) {
        pp->icm.outputIntent  = static_cast<RenderingIntent>(ointentVal);
    } else {
        pp->icm.outputIntent  = rtengine::RI_RELATIVE;
    }

    if (wTRC->get_active_row_number() == 0) {
        pp->icm.workingTRC = "none";
    } else if (wTRC->get_active_row_number() == 1) {
        pp->icm.workingTRC = "Custom";
    }

    pp->icm.toneCurve = ckbToneCurve->get_active();
    pp->icm.applyLookTable = ckbApplyLookTable->get_active();
    pp->icm.applyBaselineExposureOffset = ckbApplyBaselineExposureOffset->get_active();
    pp->icm.applyHueSatMap = ckbApplyHueSatMap->get_active();
    pp->icm.outputBPC = obpc->get_active();
    pp->icm.workingTRCGamma = (double) wGamma->getValue();
    pp->icm.workingTRCSlope = (double) wSlope->getValue();
    pp->toneCurve.fromHistMatching = false;

    if (pedited) {
        pedited->icm.inputProfile = !iunchanged->get_active();
        pedited->icm.workingProfile = wProfNames->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.outputProfile = oProfNames->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.outputIntent = oRendIntent->getSelected() < 4;
        pedited->icm.outputBPC = !obpc->get_inconsistent();
        pedited->icm.dcpIlluminant = dcpIll->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.toneCurve = !ckbToneCurve->get_inconsistent();
        pedited->icm.applyLookTable = !ckbApplyLookTable->get_inconsistent();
        pedited->icm.applyBaselineExposureOffset = !ckbApplyBaselineExposureOffset->get_inconsistent();
        pedited->icm.applyHueSatMap = !ckbApplyHueSatMap->get_inconsistent();
        pedited->icm.workingTRCGamma = wGamma->getEditedState();
        pedited->icm.workingTRCSlope = wSlope->getEditedState();
        pedited->icm.workingTRC = wTRC->get_active_text() != M("GENERAL_UNCHANGED");
    }
}

void ICMPanel::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited)
{
    wGamma->setDefault(defParams->icm.workingTRCGamma);
    wSlope->setDefault(defParams->icm.workingTRCSlope);

    if (pedited) {
        wGamma->setDefaultEditedState(pedited->icm.workingTRCGamma ? Edited : UnEdited);
        wSlope->setDefaultEditedState(pedited->icm.workingTRCSlope ? Edited : UnEdited);

    } else {
        wGamma->setDefaultEditedState(Irrelevant);
        wSlope->setDefaultEditedState(Irrelevant);

    }
}

void ICMPanel::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        Glib::ustring costr2 = Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), newval);

        if (a == wGamma) {
            listener->panelChanged(EvICMgamm, costr2);
        } else if (a == wSlope) {
            listener->panelChanged(EvICMslop, costr2);
        }

    }
}

void ICMPanel::wpChanged()
{
    if (listener) {
        listener->panelChanged(EvWProfile, wProfNames->get_active_text());
    }
}

void ICMPanel::wtrcinChanged()
{
    if (wTRC->get_active_row_number() == 0) {
        wGamma->set_sensitive(false);
        wSlope->set_sensitive(false);

    } else {
        wGamma->set_sensitive(true);
        wSlope->set_sensitive(true);
    }

    if (listener) {
        listener->panelChanged(EvICMtrcinMethod, wTRC->get_active_text());
    }
}

void ICMPanel::dcpIlluminantChanged()
{
    if (listener) {
        listener->panelChanged(EvDCPIlluminant, dcpIll->get_active_text());
    }
}

void ICMPanel::toneCurveChanged()
{
    if (multiImage) {
        if (ckbToneCurve->get_inconsistent()) {
            ckbToneCurve->set_inconsistent(false);
            tcurveconn.block(true);
            ckbToneCurve->set_active(false);
            tcurveconn.block(false);
        } else if (lastToneCurve) {
            ckbToneCurve->set_inconsistent(true);
        }

        lastToneCurve = ckbToneCurve->get_active();
    }

    if (listener) {
        if (ckbToneCurve->get_inconsistent()) {
            listener->panelChanged(EvDCPToneCurve, M("GENERAL_UNCHANGED"));
        } else if (ckbToneCurve->get_active()) {
            listener->panelChanged(EvDCPToneCurve, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDCPToneCurve, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::applyLookTableChanged()
{
    if (multiImage) {
        if (ckbApplyLookTable->get_inconsistent()) {
            ckbApplyLookTable->set_inconsistent(false);
            ltableconn.block(true);
            ckbApplyLookTable->set_active(false);
            ltableconn.block(false);
        } else if (lastApplyLookTable) {
            ckbApplyLookTable->set_inconsistent(true);
        }

        lastApplyLookTable = ckbApplyLookTable->get_active();
    }

    if (listener) {
        if (ckbApplyLookTable->get_inconsistent()) {
            listener->panelChanged(EvDCPApplyLookTable, M("GENERAL_UNCHANGED"));
        } else if (ckbApplyLookTable->get_active()) {
            listener->panelChanged(EvDCPApplyLookTable, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDCPApplyLookTable, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::applyBaselineExposureOffsetChanged()
{
    if (multiImage) {
        if (ckbApplyBaselineExposureOffset->get_inconsistent()) {
            ckbApplyBaselineExposureOffset->set_inconsistent(false);
            beoconn.block(true);
            ckbApplyBaselineExposureOffset->set_active(false);
            beoconn.block(false);
        } else if (lastApplyBaselineExposureOffset) {
            ckbApplyBaselineExposureOffset->set_inconsistent(true);
        }

        lastApplyBaselineExposureOffset = ckbApplyBaselineExposureOffset->get_active();
    }

    if (listener) {
        if (ckbApplyBaselineExposureOffset->get_inconsistent()) {
            listener->panelChanged(EvDCPApplyBaselineExposureOffset, M("GENERAL_UNCHANGED"));
        } else if (ckbApplyBaselineExposureOffset->get_active()) {
            listener->panelChanged(EvDCPApplyBaselineExposureOffset, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDCPApplyBaselineExposureOffset, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::applyHueSatMapChanged()
{
    if (multiImage) {
        if (ckbApplyHueSatMap->get_inconsistent()) {
            ckbApplyHueSatMap->set_inconsistent(false);
            hsmconn.block(true);
            ckbApplyHueSatMap->set_active(false);
            hsmconn.block(false);
        } else if (lastApplyHueSatMap) {
            ckbApplyHueSatMap->set_inconsistent(true);
        }

        lastApplyHueSatMap = ckbApplyHueSatMap->get_active();
    }

    if (listener) {
        if (ckbApplyHueSatMap->get_inconsistent()) {
            listener->panelChanged(EvDCPApplyHueSatMap, M("GENERAL_UNCHANGED"));
        } else if (ckbApplyHueSatMap->get_active()) {
            listener->panelChanged(EvDCPApplyHueSatMap, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDCPApplyHueSatMap, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::ipChanged()
{

    Glib::ustring profname;

    if (inone->get_active()) {
        profname = "(none)";
    } else if (iembedded->get_active()) {
        profname = "(embedded)";
    } else if (icamera->get_active()) {
        profname = "(camera)";
    } else if (icameraICC->get_active()) {
        profname = "(cameraICC)";
    } else {
        profname = ipDialog->get_filename();
    }

    updateDCP(-1, profname);

    if (listener && profname != oldip) {
        listener->panelChanged(EvIProfile, profname);
    }

    oldip = profname;
}

void ICMPanel::opChanged()
{
    if (!batchMode) {
        updateRenderingIntent(oProfNames->get_active_text());
    }

    if (listener) {
        listener->panelChanged(EvOProfile, oProfNames->get_active_text());
    }
}

void ICMPanel::oiChanged(int n)
{

    if (listener) {
        Glib::ustring str;

        switch (n) {
            case 0:
                str = M("PREFERENCES_INTENT_PERCEPTUAL");
                break;

            case 1:
                str = M("PREFERENCES_INTENT_RELATIVE");
                break;

            case 2:
                str = M("PREFERENCES_INTENT_SATURATION");
                break;

            case 3:
                str = M("PREFERENCES_INTENT_ABSOLUTE");
                break;

            case 4:
            default:
                str = M("GENERAL_UNCHANGED");
                break;
        }

        listener->panelChanged(EvOIntent, str);
    }
}

void ICMPanel::oBPCChanged()
{
    if (multiImage) {
        if (obpc->get_inconsistent()) {
            obpc->set_inconsistent(false);
            obpcconn.block(true);
            obpc->set_active(false);
            obpcconn.block(false);
        } else if (lastobpc) {
            obpc->set_inconsistent(true);
        }

        lastobpc = obpc->get_active();
    }

    if (listener) {
        if (obpc->get_inconsistent()) {
            listener->panelChanged(EvOBPCompens, M("GENERAL_UNCHANGED"));
        } else if (obpc->get_active()) {
            listener->panelChanged(EvOBPCompens, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvOBPCompens, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::setRawMeta(bool raw, const rtengine::FramesData* pMeta)
{

    disableListener();

    icamera->set_active(raw);
    iembedded->set_active(!raw);
    icamera->set_sensitive(raw);
    camName = pMeta->getCamera();
    icameraICC->set_sensitive(raw && (ICCStore::getInstance()->getStdProfile(pMeta->getCamera()) != nullptr || DCPStore::getInstance()->getStdProfile(pMeta->getCamera()) != nullptr));
    iembedded->set_sensitive(!raw);

    enableListener();
}

void ICMPanel::ipSelectionChanged()
{

    if (ipDialog->get_filename().empty()) {
        return;
    }

    ipChanged();
}

void ICMPanel::saveReferencePressed()
{

    if (!icmplistener) {
        return;
    }

    Gtk::FileChooserDialog dialog(getToplevelWindow(this), M("TP_ICM_SAVEREFERENCE"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    bindCurrentFolder(dialog, options.lastProfilingReferenceDir);
    dialog.set_current_name(lastRefFilename);

    dialog.add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(M("GENERAL_SAVE"), Gtk::RESPONSE_OK);

    Gtk::CheckButton applyWB(M("TP_ICM_SAVEREFERENCE_APPLYWB"));
    applyWB.set_tooltip_text(M("TP_ICM_SAVEREFERENCE_APPLYWB_TOOLTIP"));
    Gtk::HBox* hbox = Gtk::manage(new Gtk::HBox());
    hbox->pack_end(applyWB, Gtk::PACK_SHRINK, 2);
    Gtk::Box *box = dialog.get_content_area();
    box->pack_end(*hbox, Gtk::PACK_SHRINK, 2);

    Glib::RefPtr<Gtk::FileFilter> filter_tif = Gtk::FileFilter::create();
    filter_tif->set_name(M("FILECHOOSER_FILTER_TIFF"));
    filter_tif->add_pattern("*.tif");
    filter_tif->add_pattern("*.tiff");
    dialog.add_filter(filter_tif);

    Glib::RefPtr<Gtk::FileFilter> filter_any = Gtk::FileFilter::create();
    filter_any->set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any->add_pattern("*");
    dialog.add_filter(filter_any);

    dialog.show_all_children();
    //dialog.set_do_overwrite_confirmation (true);

    bool done = false;

    do {
        int result = dialog.run();

        if (result != Gtk::RESPONSE_OK) {
            done = true;
        } else {
            std::string fname = dialog.get_filename();
            Glib::ustring ext = getExtension(fname);

            if (ext != "tif" && ext != "tiff") {
                fname += ".tif";
            }

            if (confirmOverwrite(dialog, fname)) {
                icmplistener->saveInputICCReference(fname, applyWB.get_active());
                lastRefFilename = Glib::path_get_basename(fname);
                done = true;
            }
        }
    } while (!done);

    return;
}

void ICMPanel::setBatchMode(bool batchMode)
{

    isBatchMode = true;
    ToolPanel::setBatchMode(batchMode);
    iunchanged = Gtk::manage(new Gtk::RadioButton(M("GENERAL_UNCHANGED")));
    iunchanged->set_group(opts);
    iVBox->pack_start(*iunchanged, Gtk::PACK_SHRINK, 4);
    iVBox->reorder_child(*iunchanged, 5);
    removeIfThere(this, saveRef);
    oProfNames->append(M("GENERAL_UNCHANGED"));
    oRendIntent->addEntry("template-24.png", M("GENERAL_UNCHANGED"));
    oRendIntent->show();
    wProfNames->append(M("GENERAL_UNCHANGED"));
    wTRC->append(M("GENERAL_UNCHANGED"));
    dcpIll->append(M("GENERAL_UNCHANGED"));
    wGamma->showEditedCB();
    wSlope->showEditedCB();
}

