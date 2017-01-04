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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iomanip>
#include "icmpanel.h"
#include "options.h"
#include "guiutils.h"
#include "../rtengine/iccstore.h"
#include "../rtengine/dcp.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

extern Options options;

ICMPanel::ICMPanel () : FoldableToolPanel(this, "icm", M("TP_ICM_LABEL")), iunchanged(nullptr), icmplistener(nullptr), lastRefFilename(""), camName("")
{

    isBatchMode = lastToneCurve = lastApplyLookTable = lastApplyBaselineExposureOffset = lastApplyHueSatMap = lastBlendCMSMatrix = lastgamfree = false;

    ipDialog = Gtk::manage (new MyFileChooserButton (M("TP_ICM_INPUTDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    ipDialog->set_tooltip_text (M("TP_ICM_INPUTCUSTOM_TOOLTIP"));
    bindCurrentFolder (*ipDialog, options.lastIccDir);


    // ------------------------------- Input profile


    Gtk::Frame *iFrame = Gtk::manage (new Gtk::Frame(M("TP_ICM_INPUTPROFILE")) );
    iFrame->set_border_width(0);
    iFrame->set_label_align(0.025, 0.5);

    iVBox = Gtk::manage ( new Gtk::VBox());
    iVBox->set_border_width(4);
    iVBox->set_spacing(2);

    inone = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTNONE")));
    inone->set_tooltip_text (M("TP_ICM_INPUTNONE_TOOLTIP"));
    iVBox->pack_start (*inone, Gtk::PACK_SHRINK, 2);

    iembedded = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTEMBEDDED")));
    iembedded->set_tooltip_text (M("TP_ICM_INPUTEMBEDDED_TOOLTIP"));
    iVBox->pack_start (*iembedded, Gtk::PACK_SHRINK, 2);

    icamera = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCAMERA")));
    icamera->set_tooltip_text (M("TP_ICM_INPUTCAMERA_TOOLTIP"));
    iVBox->pack_start (*icamera, Gtk::PACK_SHRINK, 2);

    icameraICC = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCAMERAICC")));
    icameraICC->set_tooltip_text (M("TP_ICM_INPUTCAMERAICC_TOOLTIP"));
    iVBox->pack_start (*icameraICC, Gtk::PACK_SHRINK, 2);

    ifromfile = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCUSTOM") + ":"));
    Gtk::HBox* ffbox = Gtk::manage (new Gtk::HBox ());
    ifromfile->set_tooltip_text (M("TP_ICM_INPUTCUSTOM_TOOLTIP"));
    ffbox->pack_start (*ifromfile, Gtk::PACK_SHRINK);
    ffbox->pack_start (*ipDialog);

    iVBox->pack_start (*ffbox, Gtk::PACK_SHRINK, 2);

    opts = icamera->get_group();
    icameraICC->set_group (opts);
    iembedded->set_group (opts);
    ifromfile->set_group (opts);
    inone->set_group (opts);

    dcpFrame = Gtk::manage (new Gtk::Frame ("DCP"));
    Gtk::VBox* dcpFrameVBox = Gtk::manage (new Gtk::VBox ());
    dcpFrameVBox->set_border_width(0);
    dcpFrameVBox->set_spacing(2);

    Gtk::HBox* dcpIllHBox = Gtk::manage (new Gtk::HBox ());
    dcpIllLabel = Gtk::manage (new Gtk::Label (M("TP_ICM_DCPILLUMINANT") + ":"));
    dcpIllLabel->set_tooltip_text (M("TP_ICM_DCPILLUMINANT_TOOLTIP"));
    dcpIllLabel->show ();
    dcpIll = Gtk::manage (new MyComboBoxText ());
    dcpIll->set_tooltip_text (M("TP_ICM_DCPILLUMINANT_TOOLTIP"));
    dcpIll->append_text (M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
    dcpIll->append_text (M("TP_ICM_DCPILLUMINANT") + " 1");
    dcpIll->append_text (M("TP_ICM_DCPILLUMINANT") + " 2");
    dcpIll->show ();
    dcpTemperatures[0] = 0;
    dcpTemperatures[1] = 0;
    ignoreDcpSignal = true;
    dcpIllHBox->pack_start(*dcpIllLabel, Gtk::PACK_SHRINK, 4);
    dcpIllHBox->pack_start(*dcpIll);

    ckbToneCurve = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_TONECURVE")));
    ckbToneCurve->set_sensitive (false);
    ckbToneCurve->set_tooltip_text (M("TP_ICM_TONECURVE_TOOLTIP"));
    ckbApplyHueSatMap = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_APPLYHUESATMAP")));
    ckbApplyHueSatMap->set_sensitive (false);
    ckbApplyHueSatMap->set_tooltip_text (M("TP_ICM_APPLYHUESATMAP_TOOLTIP"));

    ckbApplyLookTable = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_APPLYLOOKTABLE")));
    ckbApplyLookTable->set_sensitive (false);
    ckbApplyLookTable->set_tooltip_text (M("TP_ICM_APPLYLOOKTABLE_TOOLTIP"));
    ckbApplyBaselineExposureOffset = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_APPLYBASELINEEXPOSUREOFFSET")));
    ckbApplyBaselineExposureOffset->set_sensitive (false);
    ckbApplyBaselineExposureOffset->set_tooltip_text (M("TP_ICM_APPLYBASELINEEXPOSUREOFFSET_TOOLTIP"));

    dcpFrameVBox->pack_start(*dcpIllHBox, Gtk::PACK_SHRINK, 0);
    dcpFrameVBox->pack_start(*ckbToneCurve, Gtk::PACK_SHRINK,0);
    dcpFrameVBox->pack_start(*ckbApplyHueSatMap, Gtk::PACK_SHRINK,0);
    dcpFrameVBox->pack_start(*ckbApplyLookTable, Gtk::PACK_SHRINK,0);
    dcpFrameVBox->pack_start(*ckbApplyBaselineExposureOffset, Gtk::PACK_SHRINK,0);
    dcpFrame->add(*dcpFrameVBox);
    dcpFrame->set_sensitive(false);
    iVBox->pack_start (*dcpFrame);

    ckbBlendCMSMatrix = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_BLENDCMSMATRIX")));
    ckbBlendCMSMatrix->set_sensitive (false);
    ckbBlendCMSMatrix->set_tooltip_text (M("TP_ICM_BLENDCMSMATRIX_TOOLTIP"));
    // blend cms matrix no longer used
    //iVBox->pack_start (*ckbBlendCMSMatrix, Gtk::PACK_SHRINK, 2);

    saveRef = Gtk::manage (new Gtk::Button (M("TP_ICM_SAVEREFERENCE")));
    saveRef->set_alignment (0.5f, 0.5f);
    saveRef->set_tooltip_markup (M("TP_ICM_SAVEREFERENCE_TOOLTIP"));
    iVBox->pack_start (*saveRef, Gtk::PACK_SHRINK, 2);

    iFrame->add(*iVBox);
    pack_start (*iFrame, Gtk::PACK_EXPAND_WIDGET, 4);


    // ---------------------------- Working profile


    Gtk::Frame *wFrame = Gtk::manage (new Gtk::Frame(M("TP_ICM_WORKINGPROFILE")) );
    wFrame->set_border_width(0);
    wFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *wVBox = Gtk::manage ( new Gtk::VBox());
    wVBox->set_border_width(4);

    wnames = Gtk::manage (new MyComboBoxText ());
    wVBox->pack_start (*wnames, Gtk::PACK_SHRINK);

    std::vector<Glib::ustring> wpnames = rtengine::getWorkingProfiles ();

    for (size_t i = 0; i < wpnames.size(); i++) {
        wnames->append_text (wpnames[i]);
    }

    wnames->set_active (0);

    wFrame->add(*wVBox);
    pack_start (*wFrame, Gtk::PACK_EXPAND_WIDGET, 4);


    // ---------------------------- Output profile


    Gtk::Frame *oFrame = Gtk::manage (new Gtk::Frame(M("TP_ICM_OUTPUTPROFILE")) );
    oFrame->set_border_width(0);
    oFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *oVBox = Gtk::manage ( new Gtk::VBox());
    oVBox->set_border_width(4);
    oVBox->set_spacing(2);

    onames = Gtk::manage (new MyComboBoxText ());
    oVBox->pack_start (*onames, Gtk::PACK_SHRINK);

    onames->append_text (M("TP_ICM_NOICM"));
    onames->set_active (0);

    std::vector<Glib::ustring> opnames = iccStore->getProfiles (rtengine::ICCStore::ProfileType::OUTPUT);

    for (size_t i = 0; i < opnames.size(); i++) {
        onames->append_text (opnames[i]);
    }

    onames->set_active (0);

    // Rendering intent
    Gtk::HBox *riHBox = Gtk::manage ( new Gtk::HBox());
    Gtk::Label* outputIntentLbl = Gtk::manage (new Gtk::Label(M("TP_ICM_PROFILEINTENT")+":"));
    riHBox->pack_start (*outputIntentLbl, Gtk::PACK_SHRINK);
    ointent = Gtk::manage (new PopUpButton ());
    ointent->addEntry("intent-perceptual.png", M("PREFERENCES_INTENT_PERCEPTUAL"));
    ointent->addEntry("intent-relative.png", M("PREFERENCES_INTENT_RELATIVE"));
    ointent->addEntry("intent-saturation.png", M("PREFERENCES_INTENT_SATURATION"));
    ointent->addEntry("intent-absolute.png", M("PREFERENCES_INTENT_ABSOLUTE"));
    ointent->setSelected (1);
    ointent->show();
    riHBox->pack_start (*ointent->buttonGroup, Gtk::PACK_EXPAND_PADDING);
    oVBox->pack_start(*riHBox, Gtk::PACK_SHRINK);

    // Black Point Compensation
    obpc = Gtk::manage(new Gtk::CheckButton((M("TP_ICM_BPC"))));
    obpc->set_active (true);
    oVBox->pack_start(*obpc, Gtk::PACK_SHRINK);

    // Output gamma

    Gtk::HBox* gaHBox = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* galab = Gtk::manage (new Gtk::Label (M("TP_GAMMA_OUTPUT") + ":"));
    //galab->set_alignment (0.0, 0.5);

    gaHBox->pack_start (*galab, Gtk::PACK_SHRINK, 4);
    wgamma = Gtk::manage (new MyComboBoxText ());
    gaHBox->pack_start (*wgamma, Gtk::PACK_EXPAND_WIDGET);

    oVBox->pack_start(*gaHBox, Gtk::PACK_EXPAND_WIDGET, 2);

    std::vector<Glib::ustring> wpgamma = rtengine::getGamma ();

    for (size_t i = 0; i < wpgamma.size(); i++) {
        wgamma->append_text (wpgamma[i]);
    }

    wgamma->set_active (0);

    Gtk::Frame* fgFrame = Gtk::manage (new Gtk::Frame ());

    Gtk::VBox *fgVBox = Gtk::manage ( new Gtk::VBox());
    fgVBox->set_spacing(0);
    fgVBox->set_border_width(4);

    freegamma = Gtk::manage(new Gtk::CheckButton((M("TP_GAMMA_FREE"))));
    freegamma->set_active (false);
    fgFrame->set_label_widget(*freegamma);

    gampos = Gtk::manage(new Adjuster (M("TP_GAMMA_CURV"), 1, 3.5, 0.01, 2.22));
    gampos->setAdjusterListener (this);

    if (gampos->delay < options.adjusterMaxDelay) {
        gampos->delay = options.adjusterMaxDelay;
    }

    gampos->show();
    slpos = Gtk::manage(new Adjuster (M("TP_GAMMA_SLOP"), 0, 15, 0.01, 4.5));
    slpos->setAdjusterListener (this);

    if (slpos->delay < options.adjusterMaxDelay) {
        slpos->delay = options.adjusterMaxDelay;
    }

    slpos->show();
    fgVBox->pack_start( *gampos, Gtk::PACK_SHRINK);//gamma
    fgVBox->pack_start( *slpos, Gtk::PACK_SHRINK);//slope

    fgFrame->add(*fgVBox);
    oVBox->pack_start(*fgFrame, Gtk::PACK_EXPAND_WIDGET, 2);

    oFrame->add(*oVBox);
    pack_start (*oFrame, Gtk::PACK_EXPAND_WIDGET, 4);


    // ---------------------------- Output gamma list entries


    Gtk::FileFilter filter_icc;
    filter_icc.set_name(M("FILECHOOSER_FILTER_COLPROF"));
    filter_icc.add_pattern("*.dcp");
    filter_icc.add_pattern("*.DCP");
    filter_icc.add_pattern("*.icc");
    filter_icc.add_pattern("*.icm");
    filter_icc.add_pattern("*.ICC");
    filter_icc.add_pattern("*.ICM");
    Gtk::FileFilter filter_iccdng;
    filter_iccdng.set_name(M("FILECHOOSER_FILTER_COLPROF") + " + DNG");
    filter_iccdng.add_pattern("*.dcp");
    filter_iccdng.add_pattern("*.DCP");
    filter_iccdng.add_pattern("*.dng");
    filter_iccdng.add_pattern("*.DNG");
    filter_iccdng.add_pattern("*.icc");
    filter_iccdng.add_pattern("*.icm");
    filter_iccdng.add_pattern("*.ICC");
    filter_iccdng.add_pattern("*.ICM");
    Gtk::FileFilter filter_any;
    filter_any.set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any.add_pattern("*");

    ipDialog->add_filter (filter_icc);
    ipDialog->add_filter (filter_iccdng);
    ipDialog->add_filter (filter_any);

    oldip = "";

    wnames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::wpChanged) );
    onames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::opChanged) );
    ointent->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::oiChanged) );
    wgamma->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::gpChanged) );
    dcpIll->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::dcpIlluminantChanged) );

    obpcconn = obpc->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::oBPCChanged) );
    gamcsconn = freegamma->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::GamChanged));
    tcurveconn = ckbToneCurve->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::toneCurveChanged));
    ltableconn = ckbApplyLookTable->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::applyLookTableChanged));
    beoconn = ckbApplyBaselineExposureOffset->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::applyBaselineExposureOffsetChanged));
    hsmconn = ckbApplyHueSatMap->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::applyHueSatMapChanged));
    blendcmsconn = ckbBlendCMSMatrix->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::blendCMSMatrixChanged));

    icamera->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    icameraICC->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    iembedded->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    ifromfile->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );

    ipc = ipDialog->signal_selection_changed().connect( sigc::mem_fun(*this, &ICMPanel::ipSelectionChanged) );
    saveRef->signal_pressed().connect( sigc::mem_fun(*this, &ICMPanel::saveReferencePressed) );

    show_all ();
}

void ICMPanel::updateRenderingIntent (const Glib::ustring &profile) {
    const uint8_t supportedIntents = rtengine::iccStore->getOutputIntents (profile);
    const bool supportsPerceptual = supportedIntents & 1 << INTENT_PERCEPTUAL;
    const bool supportsRelative   = supportedIntents & 1 << INTENT_RELATIVE_COLORIMETRIC;
    const bool supportsSaturation = supportedIntents & 1 << INTENT_SATURATION;
    const bool supportsAbsolute   = supportedIntents & 1 << INTENT_ABSOLUTE_COLORIMETRIC;

    //printf("Intents: %d / Perceptual: %d  Relative: %d  Saturation: %d  Absolute: %d\n", supportedIntents, supportsPerceptual, supportsRelative, supportsSaturation, supportsAbsolute);

    if (!profile.empty() && (supportsPerceptual || supportsRelative || supportsSaturation || supportsAbsolute)) {
        ointent->set_sensitive (true);
        ointent->setItemSensitivity(0, supportsPerceptual);
        ointent->setItemSensitivity(1, supportsRelative);
        ointent->setItemSensitivity(2, supportsSaturation);
        ointent->setItemSensitivity(3, supportsAbsolute);
    } else {
        ointent->setItemSensitivity(0, true);
        ointent->setItemSensitivity(1, true);
        ointent->setItemSensitivity(2, true);
        ointent->setItemSensitivity(3, true);
        ointent->set_sensitive (false);
        ointent->setSelected (1);
    }
}

void ICMPanel::updateDCP (int dcpIlluminant, Glib::ustring dcp_name)
{

    if (isBatchMode) {
        dcpFrame->set_sensitive(true);
        ckbToneCurve->set_sensitive (true);
        ckbApplyLookTable->set_sensitive (true);
        ckbApplyBaselineExposureOffset->set_sensitive (true);
        ckbApplyHueSatMap->set_sensitive (true);
        dcpIllLabel->set_sensitive (true);
        dcpIll->set_sensitive (true);

        if (dcpTemperatures[0] != 0 || dcpTemperatures[1] != 0) {
            int curr_active = dcpIll->get_active_row_number();
            ignoreDcpSignal = true;
            dcpIll->clear_items ();
            dcpIll->append_text (M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
            dcpIll->append_text (M("TP_ICM_DCPILLUMINANT") + " 1");
            dcpIll->append_text (M("TP_ICM_DCPILLUMINANT") + " 2");
            dcpIll->append_text (M("GENERAL_UNCHANGED"));
            dcpTemperatures[0] = 0;
            dcpTemperatures[1] = 0;
            dcpIll->set_active (curr_active);
            ignoreDcpSignal = false;
        }

        if (dcpIll->get_active_row_number() == -1 && dcpIlluminant == -1) {
            dcpIll->set_active(0);
        } else if (dcpIlluminant >= 0 && dcpIlluminant != dcpIll->get_active_row_number()) {
            dcpIll->set_active(dcpIlluminant);
        }

        dcpIll->set_sensitive (true);
        dcpIllLabel->set_sensitive (true);
        return;
    }

    ckbToneCurve->set_sensitive (false);
    ckbApplyLookTable->set_sensitive (false);
    ckbApplyBaselineExposureOffset->set_sensitive (false);
    ckbApplyHueSatMap->set_sensitive (false);
    dcpIllLabel->set_sensitive (false);
    dcpIll->set_sensitive (false);
    dcpFrame->set_sensitive(false);

    DCPProfile* dcp = nullptr;

    if(dcp_name == "(cameraICC)") {
        dcp = DCPStore::getInstance()->getStdProfile(camName);
    } else if (ifromfile->get_active() && DCPStore::getInstance()->isValidDCPFileName(dcp_name)) {
        dcp = DCPStore::getInstance()->getProfile(dcp_name);
    }

    if (dcp) {
        dcpFrame->set_sensitive(true);

        if (dcp->getHasToneCurve()) {
            ckbToneCurve->set_sensitive (true);
        }

        if (dcp->getHasLookTable()) {
            ckbApplyLookTable->set_sensitive (true);
        }

        if (dcp->getHasBaselineExposureOffset()) {
            ckbApplyBaselineExposureOffset->set_sensitive (true);
        }

        if (dcp->getHasHueSatMap()) {
            ckbApplyHueSatMap->set_sensitive (true);
        }

        const DCPProfile::Illuminants illuminants = dcp->getIlluminants();

        if (illuminants.will_interpolate) {
            if (dcpTemperatures[0] != illuminants.temperature_1 || dcpTemperatures[1] != illuminants.temperature_2) {
                char tempstr1[64], tempstr2[64];
                sprintf(tempstr1, "%.0fK", illuminants.temperature_1);
                sprintf(tempstr2, "%.0fK", illuminants.temperature_2);
                int curr_active = dcpIll->get_active_row_number();
                ignoreDcpSignal = true;
                dcpIll->clear_items ();
                dcpIll->append_text (M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
                dcpIll->append_text (tempstr1);
                dcpIll->append_text (tempstr2);
                dcpTemperatures[0] = illuminants.temperature_1;
                dcpTemperatures[1] = illuminants.temperature_2;
                dcpIll->set_active (curr_active);
                ignoreDcpSignal = false;
            }

            if (dcpIlluminant > 2) {
                dcpIlluminant = 0;
            }

            if (dcpIll->get_active_row_number() == -1 && dcpIlluminant == -1) {
                ignoreDcpSignal = true;
                dcpIll->set_active(0);
                ignoreDcpSignal = false;
            } else if (dcpIlluminant >= 0 && dcpIlluminant != dcpIll->get_active_row_number()) {
                ignoreDcpSignal = true;
                dcpIll->set_active(dcpIlluminant);
                ignoreDcpSignal = false;
            }

            dcpIll->set_sensitive (true);
            dcpIllLabel->set_sensitive (true);
        } else {
            if (dcpIll->get_active_row_number() != -1) {
                ignoreDcpSignal = true;
                dcpIll->set_active(-1);
                ignoreDcpSignal = false;
            }
        }
    }

    if (!dcpIllLabel->get_sensitive() && dcpIll->get_active_row_number() != 0) {
        if (dcpTemperatures[0] != 0 || dcpTemperatures[1] != 0) {
            int curr_active = dcpIll->get_active_row_number();
            ignoreDcpSignal = true;
            dcpIll->clear_items ();
            dcpIll->append_text (M("TP_ICM_DCPILLUMINANT_INTERPOLATED"));
            dcpIll->append_text (M("TP_ICM_DCPILLUMINANT") + " 1");
            dcpIll->append_text (M("TP_ICM_DCPILLUMINANT") + " 2");

            if (isBatchMode) {
                dcpIll->append_text (M("GENERAL_UNCHANGED"));
            }

            dcpTemperatures[0] = 0;
            dcpTemperatures[1] = 0;
            dcpIll->set_active (curr_active);
            ignoreDcpSignal = false;
        }
    }
}

void ICMPanel::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    obpcconn.block (true);
    ipc.block (true);
    gamcsconn.block (true);
    tcurveconn.block(true);
    ltableconn.block(true);
    beoconn.block(true);
    hsmconn.block(true);
    blendcmsconn.block(true);

    if(pp->icm.input.substr(0, 5) != "file:" && !ipDialog->get_filename().empty()) {
        ipDialog->set_filename(pp->icm.input);
    }

    if (pp->icm.input == "(none)") {
        inone->set_active (true);
        if (!batchMode) {
            ckbBlendCMSMatrix->set_sensitive (false);
        }
        updateDCP(pp->icm.dcpIlluminant, "");
    } else if (pp->icm.input == "(embedded)" || ((pp->icm.input == "(camera)" || pp->icm.input == "") && icamera->get_state() == Gtk::STATE_INSENSITIVE)) {
        iembedded->set_active (true);
        if (!batchMode) {
            ckbBlendCMSMatrix->set_sensitive (false);
        }
        updateDCP(pp->icm.dcpIlluminant, "");
    } else if ((pp->icm.input == "(cameraICC)") && icameraICC->get_state() != Gtk::STATE_INSENSITIVE) {
        icameraICC->set_active (true);
        if (!batchMode) {
            ckbBlendCMSMatrix->set_sensitive (true);
        }
        updateDCP(pp->icm.dcpIlluminant, "(cameraICC)");
    } else if ((pp->icm.input == "(cameraICC)") && icamera->get_state() != Gtk::STATE_INSENSITIVE && icameraICC->get_state() == Gtk::STATE_INSENSITIVE) {
        // this is the case when (cameraICC) is instructed by packaged profiles, but ICC file is not found
        // therefore falling back UI to explicitly reflect the (camera) option
        icamera->set_active (true);
        if (!batchMode) {
            ckbBlendCMSMatrix->set_sensitive (false);
        }
        updateDCP(pp->icm.dcpIlluminant, "");
    } else if ((pp->icm.input == "(cameraICC)") && icamera->get_state() == Gtk::STATE_INSENSITIVE && icameraICC->get_state() == Gtk::STATE_INSENSITIVE) {
        // If neither (camera) nor (cameraICC) are available, as is the case when loading a non-raw, activate (embedded).
        iembedded->set_active (true);
        if (!batchMode) {
            ckbBlendCMSMatrix->set_sensitive (false);
        }
        updateDCP(pp->icm.dcpIlluminant, "(cameraICC)");
    } else if ((pp->icm.input == "(camera)" || pp->icm.input == "") && icamera->get_state() != Gtk::STATE_INSENSITIVE) {
        icamera->set_active (true);
        if (!batchMode) {
            ckbBlendCMSMatrix->set_sensitive (false);
        }
        updateDCP(pp->icm.dcpIlluminant, "");
    } else {
        ifromfile->set_active (true);
        oldip = pp->icm.input.substr(5);  // cut of "file:"
        ipDialog->set_filename (pp->icm.input.substr(5));
        if (!batchMode) {
            ckbBlendCMSMatrix->set_sensitive (true);
        }
        updateDCP(pp->icm.dcpIlluminant, pp->icm.input.substr(5));
    }

    wnames->set_active_text (pp->icm.working);
    wgamma->set_active_text (pp->icm.gamma);

    if (pp->icm.output == ColorManagementParams::NoICMString) {
        onames->set_active_text (M("TP_ICM_NOICM"));
    } else {
        onames->set_active_text (pp->icm.output);
    }

    if (onames->get_active_row_number() == -1) {
        onames->set_active_text (M("TP_ICM_NOICM"));
    }
    ointent->setSelected (pp->icm.outputIntent);

    obpc->set_active (pp->icm.outputBPC);
    ckbToneCurve->set_active (pp->icm.toneCurve);
    lastToneCurve = pp->icm.toneCurve;
    ckbApplyLookTable->set_active (pp->icm.applyLookTable);
    lastApplyLookTable = pp->icm.applyLookTable;
    ckbApplyBaselineExposureOffset->set_active (pp->icm.applyBaselineExposureOffset);
    lastApplyBaselineExposureOffset = pp->icm.applyBaselineExposureOffset;
    ckbApplyHueSatMap->set_active (pp->icm.applyHueSatMap);
    lastApplyHueSatMap = pp->icm.applyHueSatMap;

    ckbBlendCMSMatrix->set_active (pp->icm.blendCMSMatrix);
    lastBlendCMSMatrix = pp->icm.blendCMSMatrix;

    freegamma->set_active (pp->icm.freegamma);
    lastgamfree = pp->icm.freegamma;

    if (!batchMode) {
        onames->set_sensitive(wgamma->get_active_row_number() == 0 && !pp->icm.freegamma); //"default"
        wgamma->set_sensitive(!pp->icm.freegamma);
        gampos->set_sensitive(pp->icm.freegamma);
        slpos->set_sensitive(pp->icm.freegamma);
        updateRenderingIntent(pp->icm.output);
    }

    gampos->setValue (pp->icm.gampos);
    slpos->setValue (pp->icm.slpos);

    if (pedited) {
        iunchanged->set_active (!pedited->icm.input);
        obpc->set_inconsistent(!pedited->icm.outputBPC);
        ckbToneCurve->set_inconsistent(!pedited->icm.toneCurve);
        ckbApplyLookTable->set_inconsistent(!pedited->icm.applyLookTable);
        ckbApplyBaselineExposureOffset->set_inconsistent(!pedited->icm.applyBaselineExposureOffset);
        ckbApplyHueSatMap->set_inconsistent(!pedited->icm.applyHueSatMap);
        ckbBlendCMSMatrix->set_inconsistent(!pedited->icm.blendCMSMatrix);
        freegamma->set_inconsistent (!pedited->icm.freegamma);

        if (!pedited->icm.working) {
            wnames->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.output) {
            onames->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.outputIntent) {
            ointent->setSelected (4);
        }

        if (!pedited->icm.dcpIlluminant) {
            dcpIll->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.gamma) {
            wgamma->set_active_text(M("GENERAL_UNCHANGED"));
            wgamma->set_active_text(M("GENERAL_UNCHANGED"));
        }

        gampos->setEditedState (pedited->icm.gampos ? Edited : UnEdited);
        slpos->setEditedState  (pedited->icm.slpos  ? Edited : UnEdited);

    }

    blendcmsconn.block(false);
    tcurveconn.block(false);
    ltableconn.block(false);
    beoconn.block(false);
    hsmconn.block(false);
    gamcsconn.block (false);
    ipc.block (false);
    obpcconn.block (false);

    enableListener ();
}

void ICMPanel::write (ProcParams* pp, ParamsEdited* pedited)
{

    if (inone->get_active()) {
        pp->icm.input = "(none)";
    } else if (iembedded->get_active ()) {
        pp->icm.input = "(embedded)";
    } else if (icamera->get_active ()) {
        pp->icm.input = "(camera)";
    } else if (icameraICC->get_active ()) {
        pp->icm.input = "(cameraICC)";
    } else {
        if (Glib::file_test (ipDialog->get_filename (), Glib::FILE_TEST_EXISTS) && !Glib::file_test (ipDialog->get_filename (), Glib::FILE_TEST_IS_DIR)) {
            pp->icm.input = "file:" + ipDialog->get_filename ();
        } else {
            pp->icm.input = "";    // just a directory
        }

        Glib::ustring p = Glib::path_get_dirname(ipDialog->get_filename ());
    }

    pp->icm.working = wnames->get_active_text ();
    pp->icm.gamma = wgamma->get_active_text ();
    pp->icm.dcpIlluminant = dcpIll->get_active_row_number();

    if (pp->icm.dcpIlluminant < 0) {
        pp->icm.dcpIlluminant = 0;
    }

    if (onames->get_active_text() == M("TP_ICM_NOICM")) {
        pp->icm.output  = ColorManagementParams::NoICMString;
    } else {
        pp->icm.output  = onames->get_active_text();
    }

    int ointentVal = ointent->getSelected ();
    if (ointentVal >= 0 && ointentVal < RI__COUNT) {
        pp->icm.outputIntent  = static_cast<RenderingIntent>(ointentVal);
    } else {
        pp->icm.outputIntent  = rtengine::RI_RELATIVE;
    }

    pp->icm.freegamma = freegamma->get_active();

    DCPProfile* dcp = nullptr;

    if (ifromfile->get_active() && pp->icm.input.substr(0, 5) == "file:" && DCPStore::getInstance()->isValidDCPFileName(pp->icm.input.substr(5))) {
        dcp = DCPStore::getInstance()->getProfile(pp->icm.input.substr(5));
    } else if(icameraICC->get_active()) {
        dcp = DCPStore::getInstance()->getStdProfile(camName);
    }

    if (dcp) {
        if (dcp->getHasToneCurve()) {
            pp->icm.toneCurve = ckbToneCurve->get_active ();
        }

        if (dcp->getHasLookTable()) {
            pp->icm.applyLookTable = ckbApplyLookTable->get_active ();
        }

        if (dcp->getHasBaselineExposureOffset()) {
            pp->icm.applyBaselineExposureOffset = ckbApplyBaselineExposureOffset->get_active ();
        }

        if (dcp->getHasHueSatMap()) {
            pp->icm.applyHueSatMap = ckbApplyHueSatMap->get_active ();
        }
    }

    pp->icm.blendCMSMatrix = ckbBlendCMSMatrix->get_active ();
    pp->icm.gampos = (double) gampos->getValue();
    pp->icm.slpos = (double) slpos->getValue();
    pp->icm.outputBPC = obpc->get_active ();

    if (pedited) {
        pedited->icm.input = !iunchanged->get_active ();
        pedited->icm.working = wnames->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.output = onames->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.outputIntent = ointent->getSelected () < 4;
        pedited->icm.outputBPC = !obpc->get_inconsistent ();
        pedited->icm.dcpIlluminant = dcpIll->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.toneCurve = !ckbToneCurve->get_inconsistent ();
        pedited->icm.applyLookTable = !ckbApplyLookTable->get_inconsistent ();
        pedited->icm.applyBaselineExposureOffset = !ckbApplyBaselineExposureOffset->get_inconsistent ();
        pedited->icm.applyHueSatMap = !ckbApplyHueSatMap->get_inconsistent ();
        pedited->icm.blendCMSMatrix = !ckbBlendCMSMatrix->get_inconsistent ();
        pedited->icm.gamma = wgamma->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.freegamma = !freegamma->get_inconsistent();
        pedited->icm.gampos = gampos->getEditedState ();
        pedited->icm.slpos = slpos->getEditedState ();
    }
}

void ICMPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
    gampos->setDefault (defParams->icm.gampos);
    slpos->setDefault (defParams->icm.slpos);

    if (pedited) {
        gampos->setDefaultEditedState (pedited->icm.gampos ? Edited : UnEdited);
        slpos->setDefaultEditedState (pedited->icm.slpos ? Edited : UnEdited);
    } else {
        gampos->setDefaultEditedState (Irrelevant);
        slpos->setDefaultEditedState (Irrelevant);
    }
}

void ICMPanel::setAdjusterBehavior (bool gammaadd, bool slopeadd)
{
    gampos->setAddMode (gammaadd);
    slpos->setAddMode (slopeadd);
}

void ICMPanel::adjusterChanged (Adjuster* a, double newval)
{

    if (listener && (freegamma->get_active() || batchMode)) {

        Glib::ustring costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), newval);

        if (a == gampos) {
            listener->panelChanged (EvGAMPOS, costr);
        } else if (a == slpos) {
            listener->panelChanged (EvSLPOS, costr);
        }
    }
}

void ICMPanel::wpChanged ()
{

    if (listener) {
        listener->panelChanged (EvWProfile, wnames->get_active_text ());
    }
}

void ICMPanel::gpChanged ()
{

    if (listener) {
        listener->panelChanged (EvGAMMA, wgamma->get_active_text ());
        onames->set_sensitive(wgamma->get_active_row_number() == 0); //"default"
    }
}

void ICMPanel::dcpIlluminantChanged()
{
    if (listener && !ignoreDcpSignal) {
        listener->panelChanged (EvDCPIlluminant, dcpIll->get_active_text ());
    }
}

void ICMPanel::toneCurveChanged()
{
    if (multiImage) {
        if (ckbToneCurve->get_inconsistent()) {
            ckbToneCurve->set_inconsistent (false);
            tcurveconn.block (true);
            ckbToneCurve->set_active (false);
            tcurveconn.block (false);
        } else if (lastToneCurve) {
            ckbToneCurve->set_inconsistent (true);
        }

        lastToneCurve = ckbToneCurve->get_active ();
    }

    if (listener) {
        if (ckbToneCurve->get_inconsistent()) {
            listener->panelChanged (EvDCPToneCurve, M("GENERAL_UNCHANGED"));
        } else if (ckbToneCurve->get_active()) {
            listener->panelChanged (EvDCPToneCurve, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDCPToneCurve, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::applyLookTableChanged()
{
    if (multiImage) {
        if (ckbApplyLookTable->get_inconsistent()) {
            ckbApplyLookTable->set_inconsistent (false);
            ltableconn.block (true);
            ckbApplyLookTable->set_active (false);
            ltableconn.block (false);
        } else if (lastApplyLookTable) {
            ckbApplyLookTable->set_inconsistent (true);
        }

        lastApplyLookTable = ckbApplyLookTable->get_active ();
    }

    if (listener) {
        if (ckbApplyLookTable->get_inconsistent()) {
            listener->panelChanged (EvDCPApplyLookTable, M("GENERAL_UNCHANGED"));
        } else if (ckbApplyLookTable->get_active()) {
            listener->panelChanged (EvDCPApplyLookTable, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDCPApplyLookTable, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::applyBaselineExposureOffsetChanged()
{
    if (multiImage) {
        if (ckbApplyBaselineExposureOffset->get_inconsistent()) {
            ckbApplyBaselineExposureOffset->set_inconsistent (false);
            beoconn.block (true);
            ckbApplyBaselineExposureOffset->set_active (false);
            beoconn.block (false);
        } else if (lastApplyBaselineExposureOffset) {
            ckbApplyBaselineExposureOffset->set_inconsistent (true);
        }

        lastApplyBaselineExposureOffset = ckbApplyBaselineExposureOffset->get_active ();
    }

    if (listener) {
        if (ckbApplyBaselineExposureOffset->get_inconsistent()) {
            listener->panelChanged (EvDCPApplyBaselineExposureOffset, M("GENERAL_UNCHANGED"));
        } else if (ckbApplyBaselineExposureOffset->get_active()) {
            listener->panelChanged (EvDCPApplyBaselineExposureOffset, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDCPApplyBaselineExposureOffset, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::applyHueSatMapChanged()
{
    if (multiImage) {
        if (ckbApplyHueSatMap->get_inconsistent()) {
            ckbApplyHueSatMap->set_inconsistent (false);
            hsmconn.block (true);
            ckbApplyHueSatMap->set_active (false);
            hsmconn.block (false);
        } else if (lastApplyHueSatMap) {
            ckbApplyHueSatMap->set_inconsistent (true);
        }

        lastApplyHueSatMap = ckbApplyHueSatMap->get_active ();
    }

    if (listener) {
        if (ckbApplyHueSatMap->get_inconsistent()) {
            listener->panelChanged (EvDCPApplyHueSatMap, M("GENERAL_UNCHANGED"));
        } else if (ckbApplyHueSatMap->get_active()) {
            listener->panelChanged (EvDCPApplyHueSatMap, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDCPApplyHueSatMap, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::ipChanged ()
{

    Glib::ustring profname;

    if (inone->get_active()) {
        profname = "(none)";
        ckbBlendCMSMatrix->set_sensitive(false);
    } else if (iembedded->get_active ()) {
        profname = "(embedded)";
        ckbBlendCMSMatrix->set_sensitive(false);
    } else if (icamera->get_active ()) {
        profname = "(camera)";
        ckbBlendCMSMatrix->set_sensitive(false);
    } else if (icameraICC->get_active ()) {
        profname = "(cameraICC)";
        ckbBlendCMSMatrix->set_sensitive(true);
    } else {
        profname = ipDialog->get_filename ();
        ckbBlendCMSMatrix->set_sensitive(true);
    }

    updateDCP(-1, profname);

    if (listener && profname != oldip) {
        listener->panelChanged (EvIProfile, profname);
    }

    oldip = profname;
}

void ICMPanel::blendCMSMatrixChanged()
{
    if (multiImage) {
        if (ckbBlendCMSMatrix->get_inconsistent()) {
            ckbBlendCMSMatrix->set_inconsistent (false);
            blendcmsconn.block (true);
            ckbBlendCMSMatrix->set_active (false);
            blendcmsconn.block (false);
        } else if (lastBlendCMSMatrix) {
            ckbBlendCMSMatrix->set_inconsistent (true);
        }

        lastBlendCMSMatrix = ckbBlendCMSMatrix->get_active ();
    }

    if (listener) {
        if (ckbBlendCMSMatrix->get_inconsistent()) {
            listener->panelChanged (EvBlendCMSMatrix, M("GENERAL_UNCHANGED"));
        } else if (ckbBlendCMSMatrix->get_active()) {
            listener->panelChanged (EvBlendCMSMatrix, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvBlendCMSMatrix, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::GamChanged()
{
    if (multiImage) {
        if (freegamma->get_inconsistent()) {
            freegamma->set_inconsistent (false);
            gamcsconn.block (true);
            freegamma->set_active (false);
            gamcsconn.block (false);
        } else if (lastgamfree) {
            freegamma->set_inconsistent (true);
        }

        lastgamfree = freegamma->get_active ();
    }

    if (listener) {
        if (freegamma->get_inconsistent()) {
            listener->panelChanged (EvGAMFREE, M("GENERAL_UNCHANGED"));
        }
        else if (freegamma->get_active()) {
            listener->panelChanged (EvGAMFREE, M("GENERAL_ENABLED"));
            if (!batchMode) {
                onames->set_sensitive(false);//disabled choice
                wgamma->set_sensitive(false);
                gampos->set_sensitive(true);
                slpos->set_sensitive(true);
            }
        } else {
            listener->panelChanged (EvGAMFREE, M("GENERAL_DISABLED"));
            if (!batchMode) {
                onames->set_sensitive(wgamma->get_active_row_number() == 0);
                wgamma->set_sensitive(true);
                gampos->set_sensitive(false);
                slpos->set_sensitive(false);
            }
        }
    }
}

void ICMPanel::opChanged ()
{
    if (!batchMode) {
        updateRenderingIntent(onames->get_active_text());
    }

    if (listener) {
        listener->panelChanged (EvOProfile, onames->get_active_text());
    }
}

void ICMPanel::oiChanged (int n)
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
        listener->panelChanged (EvOIntent, str);
    }
}

void ICMPanel::oBPCChanged ()
{
    if (multiImage) {
        if (obpc->get_inconsistent()) {
            obpc->set_inconsistent (false);
            obpcconn.block (true);
            obpc->set_active (false);
            obpcconn.block (false);
        } else if (lastobpc) {
            obpc->set_inconsistent (true);
        }

        lastobpc = obpc->get_active ();
    }

    if (listener) {
        if (obpc->get_inconsistent()) {
            listener->panelChanged (EvOBPCompens, M("GENERAL_UNCHANGED"));
        } else if (obpc->get_active()) {
            listener->panelChanged (EvOBPCompens, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvOBPCompens, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::setRawMeta (bool raw, const rtengine::ImageData* pMeta)
{

    disableListener ();

    icamera->set_active (raw);
    iembedded->set_active (!raw);
    icamera->set_sensitive (raw);
    camName = pMeta->getCamera();
    icameraICC->set_sensitive (raw && (iccStore->getStdProfile(pMeta->getCamera()) != nullptr || DCPStore::getInstance()->getStdProfile(pMeta->getCamera()) != nullptr));
    iembedded->set_sensitive (!raw);

    enableListener ();
}

void ICMPanel::ipSelectionChanged()
{

    if (ipDialog->get_filename() == "") {
        return;
    }

    ipChanged();
}

void ICMPanel::saveReferencePressed ()
{

    if (!icmplistener) {
        return;
    }

    Gtk::FileChooserDialog dialog (getToplevelWindow (this), M("TP_ICM_SAVEREFERENCE"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    bindCurrentFolder (dialog, options.lastProfilingReferenceDir);
    dialog.set_current_name (lastRefFilename);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    Gtk::CheckButton applyWB(M("TP_ICM_SAVEREFERENCE_APPLYWB"));
    applyWB.set_tooltip_text (M("TP_ICM_SAVEREFERENCE_APPLYWB_TOOLTIP"));
    Gtk::HBox* hbox = Gtk::manage( new Gtk::HBox() );
    hbox->pack_end(applyWB, Gtk::PACK_SHRINK, 2);
    Gtk::VBox *vbox = dialog.get_vbox();
    vbox->pack_end(*hbox, Gtk::PACK_SHRINK, 2);

    Gtk::FileFilter filter_tif;
    filter_tif.set_name(M("FILECHOOSER_FILTER_TIFF"));
    filter_tif.add_pattern("*.tif");
    filter_tif.add_pattern("*.tiff");
    dialog.add_filter(filter_tif);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any.add_pattern("*");
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
                icmplistener->saveInputICCReference (fname, applyWB.get_active());
                lastRefFilename = Glib::path_get_basename (fname);
                done = true;
            }
        }
    } while (!done);

    return;
}

void ICMPanel::setBatchMode (bool batchMode)
{

    isBatchMode = true;
    ignoreDcpSignal = false;
    ToolPanel::setBatchMode (batchMode);
    iunchanged = Gtk::manage (new Gtk::RadioButton (M("GENERAL_UNCHANGED")));
    iunchanged->set_group (opts);
    iVBox->pack_start (*iunchanged, Gtk::PACK_SHRINK, 4);
    iVBox->reorder_child (*iunchanged, 5);
    removeIfThere (this, saveRef);
    onames->append_text (M("GENERAL_UNCHANGED"));
    ointent->addEntry("unchanged-22.png", M("GENERAL_UNCHANGED"));
    ointent->show();
    wnames->append_text (M("GENERAL_UNCHANGED"));
    wgamma->append_text (M("GENERAL_UNCHANGED"));
    dcpIll->append_text (M("GENERAL_UNCHANGED"));
    gampos->showEditedCB ();
    slpos->showEditedCB ();
}

