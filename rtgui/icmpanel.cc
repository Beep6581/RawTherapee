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
#include "labgrid.h"
#include "options.h"
#include "pathutils.h"
#include "rtimage.h"

#include "../rtengine/dcp.h"
#include "../rtengine/iccstore.h"
#include "../rtengine/procparams.h"
#include "../rtengine/utils.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring ICMPanel::TOOL_NAME = "icm";

ICMPanel::ICMPanel() : FoldableToolPanel(this, TOOL_NAME, M("TP_ICM_LABEL")), iunchanged(nullptr), icmplistener(nullptr)
{
    auto m = ProcEventMapper::getInstance();
    EvICMprimariMethod = m->newEvent(GAMMA, "HISTORY_MSG_ICM_OUTPUT_PRIMARIES");
    EvICMprofileMethod = m->newEvent(GAMMA, "HISTORY_MSG_ICM_OUTPUT_TYPE");
    EvICMtempMethod = m->newEvent(GAMMA, "HISTORY_MSG_ICM_OUTPUT_TEMP");
    //EvICMpredx = m->newEvent(GAMMA, "HISTORY_MSG_ICMPREDX");
    //EvICMpredy = m->newEvent(GAMMA, "HISTORY_MSG_ICMPREDY");
    //EvICMpgrex = m->newEvent(GAMMA, "HISTORY_MSG_ICMPGREX");
    //EvICMpgrey = m->newEvent(GAMMA, "HISTORY_MSG_ICMPGREY");
    //EvICMpblux = m->newEvent(GAMMA, "HISTORY_MSG_ICMPBLUX");
    //EvICMpbluy = m->newEvent(GAMMA, "HISTORY_MSG_ICMPBLUY");
    EvICMgamm = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_WORKING_GAMMA");
    EvICMslop = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_WORKING_SLOPE");
    EvICMtrcinMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_WORKING_TRC_METHOD");
    EvICMwillMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_WORKING_ILLUM_METHOD");
    EvICMwprimMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_WORKING_PRIM_METHOD");
    EvICMredx = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_REDX");
    EvICMredy = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_REDY");
    EvICMgrex = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_GREX");
    EvICMgrey = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_GREY");
    EvICMblux = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_BLUX");
    EvICMbluy = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_BLUY");
    EvaIntent = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_AINTENT");
    EvICMpreser = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_PRESER");
    EvICMLabGridciexy = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICL_LABGRIDCIEXY");
    EvICMfbw = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_FBW");
    EvICMgamut = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_ICM_GAMUT");
    isBatchMode = lastToneCurve = lastApplyLookTable = lastApplyBaselineExposureOffset = lastApplyHueSatMap = false;

    ipDialog = Gtk::manage(new MyFileChooserButton(M("TP_ICM_INPUTDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    ipDialog->set_tooltip_text(M("TP_ICM_INPUTCUSTOM_TOOLTIP"));
    bindCurrentFolder(*ipDialog, options.lastIccDir);
    labgridcie = Gtk::manage(new LabGrid(EvICMLabGridciexy, M("TP_ICM_LABGRID_CIEXY"), true, true));


    // ------------------------------- Input profile


    Gtk::Frame *iFrame = Gtk::manage(new Gtk::Frame(M("TP_ICM_INPUTPROFILE")));
    iFrame->set_label_align(0.025, 0.5);

    iVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

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
    Gtk::Box* ffbox = Gtk::manage(new Gtk::Box());
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
    dcpFrame->set_label_align(0.025, 0.5);

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

    Gtk::Box* wProfVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    wProfNames = Gtk::manage(new MyComboBoxText());
    wProfVBox->pack_start(*wProfNames, Gtk::PACK_SHRINK);

    std::vector<Glib::ustring> wpnames = rtengine::ICCStore::getInstance()->getWorkingProfiles();

    for (size_t i = 0; i < wpnames.size(); i++) {
        wProfNames->append(wpnames[i]);
    }

    wProfNames->set_active(0);

    wFrame->add(*wProfVBox);

    //-----------------gamma TRC working
//    Gtk::Frame *trcFrame = Gtk::manage(new Gtk::Frame(M("TP_ICM_TRCFRAME")));
    trcExp = Gtk::manage(new MyExpander(false, M("TP_ICM_TRCFRAME")));
    setExpandAlignProperties(trcExp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
//    trcFrame->set_label_align(0.025, 0.5);
    Gtk::Box *trcProfVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    trcExp->set_tooltip_text(M("TP_ICM_TRCFRAME_TOOLTIP"));
    trcExp->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &ICMPanel::foldAllButMe), trcExp) );

    wTRCBox = Gtk::manage(new Gtk::Box());

    wTRC = Gtk::manage(new MyComboBoxText());
    wTRCBox->pack_start(*wTRC, Gtk::PACK_EXPAND_WIDGET);
    trcProfVBox->pack_start(*wTRCBox, Gtk::PACK_EXPAND_WIDGET);
    wTRC->append(M("TP_ICM_WORKING_TRC_NONE"));
    wTRC->append(M("TP_ICM_WORKING_TRC_CUSTOM"));
    wTRC->append(M("TP_ICM_WORKING_TRC_BT709"));
    wTRC->append(M("TP_ICM_WORKING_TRC_SRGB"));
    wTRC->append(M("TP_ICM_WORKING_TRC_22"));
    wTRC->append(M("TP_ICM_WORKING_TRC_18"));
    wTRC->append(M("TP_ICM_WORKING_TRC_LIN"));

    wTRC->set_active(0);
    wTRC->set_tooltip_text(M("TP_ICM_TRC_TOOLTIP"));

    wFrame->set_tooltip_text(M("TP_ICM_WORKING_TRC_TOOLTIP"));


    wGamma = Gtk::manage(new Adjuster(M("TP_ICM_WORKING_TRC_GAMMA"), 0.40, 15.0, 0.001, 2.222));
    wSlope = Gtk::manage(new Adjuster(M("TP_ICM_WORKING_TRC_SLOPE"), 0., 300., 0.01, 4.5));
    trcProfVBox->pack_start(*wGamma, Gtk::PACK_SHRINK);
    wGamma->show();

    trcProfVBox->pack_start(*wSlope, Gtk::PACK_SHRINK);
    wSlope->show();

    willuBox = Gtk::manage(new Gtk::Box());
    willulab = Gtk::manage(new Gtk::Label(M("TP_ICM_WORKING_ILLU") + ":"));

    willuBox->pack_start(*willulab, Gtk::PACK_SHRINK);
    will = Gtk::manage(new MyComboBoxText());
    willuBox->pack_start(*will, Gtk::PACK_EXPAND_WIDGET);
    trcProfVBox->pack_start(*willuBox, Gtk::PACK_EXPAND_WIDGET);
    will->append(M("TP_ICM_WORKING_ILLU_NONE"));
    will->append(M("TP_ICM_WORKING_ILLU_D41"));
    will->append(M("TP_ICM_WORKING_ILLU_D50"));
    will->append(M("TP_ICM_WORKING_ILLU_D55"));
    will->append(M("TP_ICM_WORKING_ILLU_D60"));
    will->append(M("TP_ICM_WORKING_ILLU_D65"));
    will->append(M("TP_ICM_WORKING_ILLU_D80"));
    will->append(M("TP_ICM_WORKING_ILLU_D120"));
    will->append(M("TP_ICM_WORKING_ILLU_STDA"));
    will->append(M("TP_ICM_WORKING_ILLU_2000"));
    will->append(M("TP_ICM_WORKING_ILLU_1500"));
    will->set_active(0);
    will->set_tooltip_text(M("TP_ICM_ILLUMPRIM_TOOLTIP"));


    wprimBox = Gtk::manage(new Gtk::Box());
    wprimlab = Gtk::manage(new Gtk::Label(M("TP_ICM_WORKING_PRIM") + ":"));

    wprimBox->pack_start(*wprimlab, Gtk::PACK_SHRINK);
    wprim = Gtk::manage(new MyComboBoxText());
    wprimBox->pack_start(*wprim, Gtk::PACK_EXPAND_WIDGET);
    fbw = Gtk::manage(new Gtk::CheckButton((M("TP_ICM_FBW"))));
    fbw->set_active(true);
    gamut = Gtk::manage(new Gtk::CheckButton((M("TP_ICM_GAMUT"))));
    gamut->set_active(false);

    trcProfVBox->pack_start(*wprimBox, Gtk::PACK_EXPAND_WIDGET);
    trcProfVBox->pack_start(*fbw, Gtk::PACK_EXPAND_WIDGET);
    trcProfVBox->pack_start(*gamut, Gtk::PACK_EXPAND_WIDGET);

    neutral = Gtk::manage (new Gtk::Button (M ("TP_ICM_NEUTRAL")));
    setExpandAlignProperties (neutral, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    RTImage *resetImg = Gtk::manage (new RTImage ("undo-small.png", "redo-small.png"));
    setExpandAlignProperties (resetImg, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    neutral->set_image (*resetImg);
    neutralconn = neutral->signal_pressed().connect ( sigc::mem_fun (*this, &ICMPanel::neutral_pressed) );
    neutral->show();

    trcProfVBox->pack_start (*neutral);


    wprim->append(M("TP_ICM_WORKING_PRIM_NONE"));
    wprim->append(M("TP_ICM_WORKING_PRIM_SRGB"));
    wprim->append(M("TP_ICM_WORKING_PRIM_ADOB"));
    wprim->append(M("TP_ICM_WORKING_PRIM_PROP"));
    wprim->append(M("TP_ICM_WORKING_PRIM_REC"));
    wprim->append(M("TP_ICM_WORKING_PRIM_ACE"));
    wprim->append(M("TP_ICM_WORKING_PRIM_WID"));
    wprim->append(M("TP_ICM_WORKING_PRIM_AC0"));
    wprim->append(M("TP_ICM_WORKING_PRIM_JDCMAX"));
    wprim->append(M("TP_ICM_WORKING_PRIM_BRU"));
    wprim->append(M("TP_ICM_WORKING_PRIM_BET"));
    wprim->append(M("TP_ICM_WORKING_PRIM_BST"));
    wprim->append(M("TP_ICM_WORKING_PRIM_CUS"));
    wprim->append(M("TP_ICM_WORKING_PRIM_CUSGR"));
    wprim->set_active(0);

    wprim->set_tooltip_text(M("TP_ICM_PRIMILLUM_TOOLTIP"));


    redx = Gtk::manage(new Adjuster(M("TC_PRIM_REDX"), 0.41, 1.0, 0.0001, 0.7347));
    setExpandAlignProperties(redx, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    redy = Gtk::manage(new Adjuster(M("TC_PRIM_REDY"), 0.0, 0.70, 0.0001, 0.2653));
    setExpandAlignProperties(redy, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    grex = Gtk::manage(new Adjuster(M("TC_PRIM_GREX"), -0.1, 0.4, 0.0001, 0.1596));
    setExpandAlignProperties(grex, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    grey = Gtk::manage(new Adjuster(M("TC_PRIM_GREY"), 0.50, 1.0, 0.0001, 0.8404));
    setExpandAlignProperties(grey, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    blux = Gtk::manage(new Adjuster(M("TC_PRIM_BLUX"), -0.1, 0.4, 0.0001, 0.0366));
    setExpandAlignProperties(blux, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    bluy = Gtk::manage(new Adjuster(M("TC_PRIM_BLUY"), -0.1, 0.49, 0.0001, 0.0001));
    setExpandAlignProperties(bluy, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    redx->set_tooltip_text(M("TP_ICM_PRIMRED_TOOLTIP"));
    grex->set_tooltip_text(M("TP_ICM_PRIMGRE_TOOLTIP"));
    blux->set_tooltip_text(M("TP_ICM_PRIMBLU_TOOLTIP"));

    redFrame = Gtk::manage(new Gtk::Frame(M("TP_ICM_REDFRAME")));
    redFrame->set_label_align(0.025, 0.5);
    redFrame->set_tooltip_text(M("TP_ICM_WORKING_PRIMFRAME_TOOLTIP"));

    Gtk::Box *redVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    primCoordGrid = Gtk::manage(new Gtk::Grid());
    primCoordGrid->set_column_homogeneous(true);
    primCoordGrid->attach(*redx, 0, 0, 1, 1);
    primCoordGrid->attach_next_to(*redy, *redx, Gtk::PositionType::POS_RIGHT, 1, 1);
    primCoordGrid->attach_next_to(*grex, *redx, Gtk::PositionType::POS_BOTTOM, 1, 1);
    primCoordGrid->attach_next_to(*grey, *grex, Gtk::PositionType::POS_RIGHT, 1, 1);
    primCoordGrid->attach_next_to(*blux, *grex, Gtk::PositionType::POS_BOTTOM, 1, 1);
    primCoordGrid->attach_next_to(*bluy, *blux, Gtk::PositionType::POS_RIGHT, 1, 1);
    redVBox->pack_start(*primCoordGrid, Gtk::PACK_EXPAND_WIDGET);

    Gtk::Separator* const separator1 = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_VERTICAL));
    Gtk::Separator* const separator2 = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_VERTICAL));

    preser = Gtk::manage(new Adjuster(M("TP_ICM_WORKING_PRESER"), 0., 100., 0.5, 0.));
    preser->setAdjusterListener(this);
    
    preBox = Gtk::manage(new Gtk::Box());
    preBox->pack_start(*preser, Gtk::PACK_EXPAND_WIDGET);
    redVBox->pack_start(*separator1, Gtk::PACK_SHRINK);
    redVBox->pack_start(*preBox, Gtk::PACK_EXPAND_WIDGET);
    redVBox->pack_start(*separator2, Gtk::PACK_SHRINK);

    cielab = Gtk::manage(new Gtk::Label(M("TP_ICM_WORKING_CIEDIAG") + ":"));

    redVBox->pack_start(*cielab, Gtk::PACK_SHRINK);

    redVBox->pack_start(*labgridcie, Gtk::PACK_EXPAND_WIDGET, 4);
    
    redFrame->add(*redVBox);

    wGamma->setAdjusterListener(this);
    wSlope->setLogScale(16, 0);
    wSlope->setAdjusterListener(this);
    redx->setAdjusterListener(this);
    redy->setAdjusterListener(this);
    grex->setAdjusterListener(this);
    grey->setAdjusterListener(this);
    blux->setAdjusterListener(this);
    bluy->setAdjusterListener(this);

    wGamma->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    wSlope->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    // Rendering intent
    riaHBox = Gtk::manage(new Gtk::Box());
    Gtk::Label* abIntentLbl = Gtk::manage(new Gtk::Label(M("TP_ICM_PROFILEINTENT")));
    riaHBox->pack_start(*abIntentLbl, Gtk::PACK_SHRINK);
    aRendIntent.reset(new PopUpButton());
    aRendIntent->addEntry("intent-perceptual.png", M("PREFERENCES_INTENT_PERCEPTUAL"));
    aRendIntent->addEntry("intent-relative.png", M("PREFERENCES_INTENT_RELATIVE"));
    aRendIntent->addEntry("intent-saturation.png", M("PREFERENCES_INTENT_SATURATION"));
    aRendIntent->addEntry("intent-absolute.png", M("PREFERENCES_INTENT_ABSOLUTE"));
    aRendIntent->setSelected(1);
    aRendIntent->show();
    riaHBox->pack_start(*aRendIntent->buttonGroup, Gtk::PACK_EXPAND_PADDING);


    pack_start(*wFrame, Gtk::PACK_EXPAND_WIDGET);
    trcProfVBox->pack_start(*redFrame, Gtk::PACK_EXPAND_WIDGET);
    trcExp->add(*trcProfVBox, false);
    trcExp->setLevel (2);
    pack_start(*trcExp, Gtk::PACK_EXPAND_WIDGET);
    trcExp->set_expanded(false);


    // ---------------------------- Output profile


    Gtk::Frame *oFrame = Gtk::manage(new Gtk::Frame(M("TP_ICM_OUTPUTPROFILE")));
    oFrame->set_label_align(0.025, 0.5);
    oFrame->set_tooltip_text(M("TP_ICM_OUTPUTPROFILE_TOOLTIP"));

    Gtk::Box* oProfVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

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
    Gtk::Box *riHBox = Gtk::manage(new Gtk::Box());
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
#ifdef _WIN32
    ipDialog->set_show_hidden(true);  // ProgramData is hidden on Windows
#endif

    wprofnamesconn = wProfNames->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::wpChanged));
    oprofnamesconn = oProfNames->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::opChanged));
    orendintentconn = oRendIntent->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::oiChanged));
    arendintentconn = aRendIntent->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::aiChanged));
    dcpillconn = dcpIll->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::dcpIlluminantChanged));
    wtrcconn = wTRC->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::wtrcinChanged));
    willconn = will->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::willChanged));
    wprimconn = wprim->signal_changed().connect(sigc::mem_fun(*this, &ICMPanel::wprimChanged));

    fbwconn = fbw->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::fbwChanged));
    gamutconn = gamut->signal_toggled().connect(sigc::mem_fun(*this, &ICMPanel::gamutChanged));
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

void ICMPanel::foldAllButMe (GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        trcExp->set_expanded (trcExp == expander);
    }
}


void ICMPanel::neutral_pressed ()
{   //find working profile and set the same destination proile
    if (wProfNames->get_active_text() == "Rec2020") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::REC2020));
    } else if (wProfNames->get_active_text() == "sRGB") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::SRGB));
    } else if (wProfNames->get_active_text() == "Adobe RGB") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::ADOBE_RGB));
    } else if (wProfNames->get_active_text() == "ProPhoto") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::PRO_PHOTO));
    } else if (wProfNames->get_active_text() == "ACESp1") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::ACES_P1));
    } else if (wProfNames->get_active_text() == "WideGamut") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::WIDE_GAMUT));
    } else if (wProfNames->get_active_text() == "ACESp0") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::ACES_P0));
    } else if (wProfNames->get_active_text() == "JDCmax") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::JDC_MAX));
    } else if (wProfNames->get_active_text() == "BruceRGB") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::BRUCE_RGB));
    } else if (wProfNames->get_active_text() == "Beta RGB") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::BETA_RGB));
    } else if (wProfNames->get_active_text() == "BestRGB") {
        wprim->set_active(toUnderlying(ColorManagementParams::Primaries::BEST_RGB));
    }
    const ColorManagementParams defPar;
    wGamma->setValue(defPar.workingTRCGamma);//2.4
    wSlope->setValue(defPar.workingTRCSlope);//12.92
    preser->setValue(defPar.preser);
    fbw->set_active(defPar.fbw);
    gamut->set_active(defPar.gamut);	
    wTRC->set_active(toUnderlying(ColorManagementParams::WorkingTrc::NONE));//reset to none
    will->set_active(toUnderlying(ColorManagementParams::Illuminant::DEFAULT));//reset to default - after wprim
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

        aRendIntent->set_sensitive(true);
        aRendIntent->setItemSensitivity(0, supportsPerceptual);
        aRendIntent->setItemSensitivity(1, supportsRelative);
        aRendIntent->setItemSensitivity(2, supportsSaturation);
        aRendIntent->setItemSensitivity(3, supportsAbsolute);

    } else {
        oRendIntent->setItemSensitivity(0, true);
        oRendIntent->setItemSensitivity(1, true);
        oRendIntent->setItemSensitivity(2, true);
        oRendIntent->setItemSensitivity(3, true);
        oRendIntent->set_sensitive(false);
        oRendIntent->setSelected(1);

        aRendIntent->setItemSensitivity(0, true);
        aRendIntent->setItemSensitivity(1, true);
        aRendIntent->setItemSensitivity(2, true);
        aRendIntent->setItemSensitivity(3, true);
        aRendIntent->set_sensitive(false);
        aRendIntent->setSelected(1);

    }
}

ICMPanel::~ICMPanel()
{
    idle_register.destroy();
}

void ICMPanel::primChanged (float rx, float ry, float bx, float by, float gx, float gy)
{ //update sliders R G B Ciexy
    nextrx = rx;
    nextry = ry;
    nextbx = bx;
    nextby = by;
    nextgx = gx;
    nextgy = gy;

    idle_register.add(
        [this]() -> bool
        {
            disableListener();
            redx->setValue(nextrx);
            redy->setValue(nextry);
            blux->setValue(nextbx);
            bluy->setValue(nextby);
            grex->setValue(nextgx);
            grey->setValue(nextgy);

            enableListener();
            return false;
        }
    );
}

void ICMPanel::iprimChanged (float r_x, float r_y, float b_x, float b_y, float g_x, float g_y, float w_x, float w_y)
{//update CIE xy graph
    nextrx = r_x;
    nextry = r_y;
    nextbx = b_x;
    nextby = b_y;
    nextgx = g_x;
    nextgy = g_y;
    nextwx = w_x;
    nextwy = w_y;
    //convert xy datas in datas for labgrid areas
    nextrx = 1.81818f * (nextrx + 0.1f) - 1.f;
    nextry = 1.81818f * (nextry + 0.1f) - 1.f;
    nextbx = 1.81818f * (nextbx + 0.1f) - 1.f;
    nextby = 1.81818f * (nextby + 0.1f) - 1.f;
    nextgx = 1.81818f * (nextgx + 0.1f) - 1.f;
    nextgy = 1.81818f * (nextgy + 0.1f) - 1.f;
    nextwx = 1.81818f * (nextwx + 0.1f) - 1.f;
    nextwy = 1.81818f * (nextwy + 0.1f) - 1.f;

    idle_register.add(
        [this]() -> bool
        {
            disableListener();
            labgridcie->setParams(nextrx, nextry, nextbx, nextby, nextgx, nextgy, nextwx, nextwy, false);
            enableListener();
            return false;
        }
    );
}


void ICMPanel::setEditProvider(EditDataProvider *provider)
{
    //in case of
}

void ICMPanel::setListener(ToolPanelListener *tpl)
{//enable listener Toolpanel and Labgridcie
        ToolPanel::setListener(tpl);
        labgridcie->setListener(tpl);
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
                snprintf(tempstr1, sizeof(tempstr1), "%.0fK", illuminants.temperature_1);
                snprintf(tempstr2, sizeof(tempstr2), "%.0fK", illuminants.temperature_2);
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
    ConnectionBlocker fbwconn_(fbwconn);
    ConnectionBlocker gamutconn_(gamutconn);
    ConnectionBlocker ipc_(ipc);
    ConnectionBlocker tcurveconn_(tcurveconn);
    ConnectionBlocker ltableconn_(ltableconn);
    ConnectionBlocker beoconn_(beoconn);
    ConnectionBlocker hsmconn_(hsmconn);
    ConnectionBlocker wprofnamesconn_(wprofnamesconn);
    ConnectionBlocker oprofnamesconn_(oprofnamesconn);
    ConnectionBlocker orendintentconn_(orendintentconn);
    ConnectionBlocker arendintentconn_(arendintentconn);
    ConnectionBlocker dcpillconn_(dcpillconn);
    ConnectionBlocker wtrcconn_(wtrcconn);
    ConnectionBlocker willconn_(willconn);
    ConnectionBlocker wprimconn_(wprimconn);
    trcExp->set_expanded(false);

    if (pp->icm.inputProfile.substr(0, 5) != "file:") {
        ipDialog->set_filename(" ");
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

    wTRC->set_active(rtengine::toUnderlying(pp->icm.workingTRC));

    will->set_active(rtengine::toUnderlying(pp->icm.will));

    wprim->set_active(rtengine::toUnderlying(pp->icm.wprim));

    wtrcinChanged();
    willChanged();
    wprimChanged();

    if (pp->icm.outputProfile == ColorManagementParams::NoICMString) {
        oProfNames->set_active_text(M("TP_ICM_NOICM"));
    } else {
        oProfNames->set_active_text(pp->icm.outputProfile);
    }

    if (oProfNames->get_active_row_number() == -1) {
        oProfNames->set_active_text(M("TP_ICM_NOICM"));
    }

    oRendIntent->setSelected(pp->icm.outputIntent);
    aRendIntent->setSelected(pp->icm.aRendIntent);

    obpc->set_active(pp->icm.outputBPC);
    fbw->set_active(pp->icm.fbw);
    gamut->set_active(pp->icm.gamut);
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
    redx->setValue(pp->icm.redx);
    redy->setValue(pp->icm.redy);
    grex->setValue(pp->icm.grex);
    grey->setValue(pp->icm.grey);
    blux->setValue(pp->icm.blux);
    bluy->setValue(pp->icm.bluy);
    preser->setValue(pp->icm.preser);
    labgridcie->setParams(pp->icm.labgridcieALow, pp->icm.labgridcieBLow, pp->icm.labgridcieAHigh, pp->icm.labgridcieBHigh, pp->icm.labgridcieGx, pp->icm.labgridcieGy, pp->icm.labgridcieWx, pp->icm.labgridcieWy, false);

    if (pedited) {
        iunchanged->set_active(!pedited->icm.inputProfile);
        obpc->set_inconsistent(!pedited->icm.outputBPC);
        fbw->set_inconsistent(!pedited->icm.fbw);
        gamut->set_inconsistent(!pedited->icm.gamut);
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

        if (!pedited->icm.aRendIntent) {
            aRendIntent->setSelected(4);
        }

        if (!pedited->icm.dcpIlluminant) {
            dcpIll->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.workingTRC) {
            wTRC->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.will) {
            will->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->icm.wprim) {
            wprim->set_active_text(M("GENERAL_UNCHANGED"));
        }
        labgridcie->setEdited(pedited->icm.labgridcieALow || pedited->icm.labgridcieBLow || pedited->icm.labgridcieAHigh || pedited->icm.labgridcieBHigh  || pedited->icm.labgridcieGx  || pedited->icm.labgridcieGy || pedited->icm.labgridcieWx  || pedited->icm.labgridcieWy);

        wGamma->setEditedState(pedited->icm.workingTRCGamma ? Edited : UnEdited);
        wSlope->setEditedState(pedited->icm.workingTRCSlope  ? Edited : UnEdited);
        redx->setEditedState(pedited->icm.redx  ? Edited : UnEdited);
        redy->setEditedState(pedited->icm.redy  ? Edited : UnEdited);
        grex->setEditedState(pedited->icm.grex  ? Edited : UnEdited);
        grey->setEditedState(pedited->icm.grey  ? Edited : UnEdited);
        blux->setEditedState(pedited->icm.blux  ? Edited : UnEdited);
        bluy->setEditedState(pedited->icm.bluy  ? Edited : UnEdited);
        preser->setEditedState(pedited->icm.preser  ? Edited : UnEdited);

    }

    switch (ColorManagementParams::WorkingTrc(wTRC->get_active_row_number())) {
        case ColorManagementParams::WorkingTrc::NONE: {
            wSlope->set_sensitive(false);
            wGamma->set_sensitive(false);
            will->set_sensitive(false);
            willulab->set_sensitive(false);
            wprim->set_sensitive(false);
            fbw->set_sensitive(false);
            gamut->set_sensitive(false);
            wprimlab->set_sensitive(false);
            riaHBox->set_sensitive(false);
            redFrame->hide();
            break;
        }

        case ColorManagementParams::WorkingTrc::CUSTOM: {
            will->set_sensitive(false);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
                if (
                    ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM
                    && ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM_GRID
                ) {
                  will->set_sensitive(false);
                  primCoordGrid->set_sensitive(false);
                  labgridcie->set_sensitive(false);

                } else {
                  will->set_sensitive(false);
                  if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::CUSTOM) {
                    will->set_sensitive(true);
                  }
                  primCoordGrid->set_sensitive(true);
                  labgridcie->set_sensitive(true);
                }

            }
            riaHBox->set_sensitive(true);

            if (pp->icm.workingTRCGamma <= 1.) {
                wGamma->set_sensitive(true);
                wSlope->set_sensitive(false);
            } else {
                wGamma->set_sensitive(true);
                wSlope->set_sensitive(true);
            }
            break;
        }

        case ColorManagementParams::WorkingTrc::BT709:
            wGamma->setValue(2.222);
            wSlope->setValue(4.5);
            will->set_sensitive(true);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
            }
            riaHBox->set_sensitive(true);
            break;
        case ColorManagementParams::WorkingTrc::SRGB:
            wGamma->setValue(2.4);
            wSlope->setValue(12.92);
            will->set_sensitive(true);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
            }
            break;
        case ColorManagementParams::WorkingTrc::GAMMA_2_2:
            wGamma->setValue(2.2);
            wSlope->setValue(0.);
            will->set_sensitive(true);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            redFrame->show();
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
            }
            break;
        case ColorManagementParams::WorkingTrc::GAMMA_1_8:
            wGamma->setValue(1.8);
            wSlope->setValue(0.);
            will->set_sensitive(true);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
            }
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            break;
        case ColorManagementParams::WorkingTrc::LINEAR:
            wGamma->setValue(1.);
            wSlope->setValue(1.);
            will->set_sensitive(true);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
            }
            break;
    }

    switch (ColorManagementParams::Primaries(wprim->get_active_row_number())) {
        case ColorManagementParams::Primaries::DEFAULT:
        case ColorManagementParams::Primaries::SRGB:
        case ColorManagementParams::Primaries::ADOBE_RGB:
        case ColorManagementParams::Primaries::PRO_PHOTO:
        case ColorManagementParams::Primaries::REC2020:
        case ColorManagementParams::Primaries::ACES_P1:
        case ColorManagementParams::Primaries::WIDE_GAMUT:
        case ColorManagementParams::Primaries::ACES_P0:
        case ColorManagementParams::Primaries::JDC_MAX:
        case ColorManagementParams::Primaries::BRUCE_RGB:
        case ColorManagementParams::Primaries::BETA_RGB:
        case ColorManagementParams::Primaries::BEST_RGB: {
            labgridcie->set_sensitive(false);
            will->set_sensitive(false);
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM: {
            will->set_sensitive(true);
            labgridcie->set_sensitive(false);
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM_GRID: {
            labgridcie->set_sensitive(true);
            primCoordGrid->set_sensitive(false);
            will->set_sensitive(false);
            break;
        }
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
    labgridcie->getParams(pp->icm.labgridcieALow, pp->icm.labgridcieBLow, pp->icm.labgridcieAHigh, pp->icm.labgridcieBHigh, pp->icm.labgridcieGx, pp->icm.labgridcieGy, pp->icm.labgridcieWx, pp->icm.labgridcieWy);

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

    int aintentVal = aRendIntent->getSelected();

    if (aintentVal >= 0 && aintentVal < RI__COUNT) {
        pp->icm.aRendIntent  = static_cast<RenderingIntent>(aintentVal);
    } else {
        pp->icm.aRendIntent  = rtengine::RI_RELATIVE;
    }

    pp->icm.workingTRC = ColorManagementParams::WorkingTrc(wTRC->get_active_row_number());
    pp->icm.will = ColorManagementParams::Illuminant(will->get_active_row_number());
    pp->icm.wprim = ColorManagementParams::Primaries(wprim->get_active_row_number());

    pp->icm.toneCurve = ckbToneCurve->get_active();
    pp->icm.applyLookTable = ckbApplyLookTable->get_active();
    pp->icm.applyBaselineExposureOffset = ckbApplyBaselineExposureOffset->get_active();
    pp->icm.applyHueSatMap = ckbApplyHueSatMap->get_active();
    pp->icm.outputBPC = obpc->get_active();
    pp->icm.fbw = fbw->get_active();
    pp->icm.gamut = gamut->get_active();
    pp->icm.workingTRCGamma =  wGamma->getValue();
    pp->icm.workingTRCSlope =  wSlope->getValue();
    pp->icm.redx =  redx->getValue();
    pp->icm.redy =  redy->getValue();
    pp->icm.grex =  grex->getValue();
    pp->icm.grey =  grey->getValue();
    pp->icm.blux =  blux->getValue();
    pp->icm.bluy =  bluy->getValue();
    pp->toneCurve.fromHistMatching = false;
    pp->icm.preser =  preser->getValue();

    if (pedited) {
        pedited->icm.inputProfile = !iunchanged->get_active();
        pedited->icm.workingProfile = wProfNames->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.outputProfile = oProfNames->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.outputIntent = oRendIntent->getSelected() < 4;
        pedited->icm.aRendIntent = aRendIntent->getSelected() < 4;
        pedited->icm.outputBPC = !obpc->get_inconsistent();
        pedited->icm.fbw = !fbw->get_inconsistent();
        pedited->icm.gamut = !gamut->get_inconsistent();
        pedited->icm.dcpIlluminant = dcpIll->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.toneCurve = !ckbToneCurve->get_inconsistent();
        pedited->icm.applyLookTable = !ckbApplyLookTable->get_inconsistent();
        pedited->icm.applyBaselineExposureOffset = !ckbApplyBaselineExposureOffset->get_inconsistent();
        pedited->icm.applyHueSatMap = !ckbApplyHueSatMap->get_inconsistent();
        pedited->icm.workingTRCGamma = wGamma->getEditedState();
        pedited->icm.workingTRCSlope = wSlope->getEditedState();
        pedited->icm.workingTRC = wTRC->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.will = will->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.wprim = wprim->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->icm.redx = redx->getEditedState();
        pedited->icm.redy = redy->getEditedState();
        pedited->icm.labgridcieALow = pedited->icm.labgridcieBLow = pedited->icm.labgridcieAHigh = pedited->icm.labgridcieBHigh = pedited->icm.labgridcieGx = pedited->icm.labgridcieGy = pedited->icm.labgridcieWx = pedited->icm.labgridcieWy = labgridcie->getEdited();
   }
}

void ICMPanel::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited)
{
    wGamma->setDefault(defParams->icm.workingTRCGamma);
    wSlope->setDefault(defParams->icm.workingTRCSlope);
    redx->setDefault(defParams->icm.redx);
    redy->setDefault(defParams->icm.redy);
    grex->setDefault(defParams->icm.grex);
    grey->setDefault(defParams->icm.grey);
    blux->setDefault(defParams->icm.blux);
    bluy->setDefault(defParams->icm.bluy);
    preser->setDefault(defParams->icm.preser);
    labgridcie->setDefault(defParams->icm.labgridcieALow, defParams->icm.labgridcieBLow , defParams->icm.labgridcieAHigh, defParams->icm.labgridcieBHigh, defParams->icm.labgridcieGx, defParams->icm.labgridcieGy, defParams->icm.labgridcieWx, defParams->icm.labgridcieWy);

    if (pedited) {
        wGamma->setDefaultEditedState(pedited->icm.workingTRCGamma ? Edited : UnEdited);
        wSlope->setDefaultEditedState(pedited->icm.workingTRCSlope ? Edited : UnEdited);
        redx->setDefaultEditedState(pedited->icm.redx ? Edited : UnEdited);
        redy->setDefaultEditedState(pedited->icm.redy ? Edited : UnEdited);
        grex->setDefaultEditedState(pedited->icm.grex ? Edited : UnEdited);
        grey->setDefaultEditedState(pedited->icm.grey ? Edited : UnEdited);
        blux->setDefaultEditedState(pedited->icm.blux ? Edited : UnEdited);
        bluy->setDefaultEditedState(pedited->icm.bluy ? Edited : UnEdited);
        labgridcie->setEdited((pedited->icm.labgridcieALow || pedited->icm.labgridcieBLow || pedited->icm.labgridcieAHigh || pedited->icm.labgridcieBHigh || pedited->icm.labgridcieGx || pedited->icm.labgridcieGy || pedited->icm.labgridcieWx || pedited->icm.labgridcieWy) ? Edited : UnEdited);
        preser->setDefaultEditedState(pedited->icm.preser ? Edited : UnEdited);

    } else {
        wGamma->setDefaultEditedState(Irrelevant);
        wSlope->setDefaultEditedState(Irrelevant);
        redx->setDefaultEditedState(Irrelevant);
        redy->setDefaultEditedState(Irrelevant);
        grex->setDefaultEditedState(Irrelevant);
        grey->setDefaultEditedState(Irrelevant);
        blux->setDefaultEditedState(Irrelevant);
        bluy->setDefaultEditedState(Irrelevant);
        preser->setDefaultEditedState(Irrelevant);
        labgridcie->setEdited(Edited);

    }
}

void ICMPanel::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        Glib::ustring costr2 = Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), newval);

        if (a == wGamma) {
            if (wGamma->getValue() <= 1.) {
                wSlope->set_sensitive(false);
            } else {
                wSlope->set_sensitive(true);
            }
            listener->panelChanged(EvICMgamm, costr2);
        } else if (a == wSlope) {
            listener->panelChanged(EvICMslop, costr2);
        } else if (a == redx) {
            listener->panelChanged(EvICMredx, costr2);
        } else if (a == redy) {
            listener->panelChanged(EvICMredy, costr2);
        } else if (a == grex) {
            listener->panelChanged(EvICMgrex, costr2);
        } else if (a == grey) {
            listener->panelChanged(EvICMgrey, costr2);
        } else if (a == blux) {
            listener->panelChanged(EvICMblux, costr2);
        } else if (a == bluy) {
            listener->panelChanged(EvICMbluy, costr2);
        } else if (a == preser) {
            listener->panelChanged(EvICMpreser, costr2);
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
    switch (ColorManagementParams::WorkingTrc(wTRC->get_active_row_number())) {
        case ColorManagementParams::WorkingTrc::NONE: {
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            will->set_sensitive(false);
            willulab->set_sensitive(false);
            wprim->set_sensitive(false);
            fbw->set_sensitive(false);
            gamut->set_sensitive(false);
            wprimlab->set_sensitive(false);
            redFrame->hide();
            riaHBox->set_sensitive(false);
            break;
        }

        case ColorManagementParams::WorkingTrc::CUSTOM: {
            will->set_sensitive(false);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            willulab->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
                if (
                    ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM
                    && ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM_GRID
                ) {
                   primCoordGrid->set_sensitive(false);
                } else {
                   primCoordGrid->set_sensitive(true);
                }
            }
            riaHBox->set_sensitive(true);
            if (wGamma->getValue() <= 1.) {
                wGamma->set_sensitive(true);
                wSlope->set_sensitive(false);
            } else {
                wGamma->set_sensitive(true);
                wSlope->set_sensitive(true);
            }
            break;
        }

        case ColorManagementParams::WorkingTrc::BT709: {
            wGamma->setValue(2.222);
            wSlope->setValue(4.5);
            will->set_sensitive(false);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
                if (
                    ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM
                    && ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM_GRID
                ) {
                    primCoordGrid->set_sensitive(false);
                }
            }
            riaHBox->set_sensitive(true);
            break;
        }

        case ColorManagementParams::WorkingTrc::SRGB: {
            wGamma->setValue(2.4);
            wSlope->setValue(12.92);
            will->set_sensitive(false);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
                if (
                    ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM
                    && ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM_GRID
                ) {
                    primCoordGrid->set_sensitive(false);
                } else {
                    primCoordGrid->set_sensitive(true);
                }
            }
            break;
        }

        case ColorManagementParams::WorkingTrc::GAMMA_2_2: {
            wGamma->setValue(2.2);
            wSlope->setValue(0.);
            will->set_sensitive(false);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
                if (
                    ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM
                    && ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM_GRID
                ) {
                    primCoordGrid->set_sensitive(false);
                } else {
                    primCoordGrid->set_sensitive(true);
                }
            }
            break;
        }

        case ColorManagementParams::WorkingTrc::GAMMA_1_8: {
            wGamma->setValue(1.8);
            wSlope->setValue(0.);
            will->set_sensitive(false);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
                if (
                    ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM
                    && ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM_GRID
                ) {
                    primCoordGrid->set_sensitive(false);
                } else {
                    primCoordGrid->set_sensitive(true);
                }
            }
            break;
        }

        case ColorManagementParams::WorkingTrc::LINEAR: {
            wGamma->setValue(1.0);
            wSlope->setValue(1.);
            will->set_sensitive(false);
            willulab->set_sensitive(true);
            wprim->set_sensitive(true);
            fbw->set_sensitive(true);
            gamut->set_sensitive(true);
            wprimlab->set_sensitive(true);
            wGamma->set_sensitive(false);
            wSlope->set_sensitive(false);
            riaHBox->set_sensitive(true);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
                redFrame->hide();
            } else {
                redFrame->show();
                if (
                    ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM
                    && ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM_GRID
                ) {
                    primCoordGrid->set_sensitive(false);
                } else {
                    primCoordGrid->set_sensitive(true);
                }
            }
            break;
        }
    }
    wprimChanged();

    switch (ColorManagementParams::Primaries(wprim->get_active_row_number())) {
        case ColorManagementParams::Primaries::DEFAULT:
        case ColorManagementParams::Primaries::SRGB:
        case ColorManagementParams::Primaries::ADOBE_RGB:
        case ColorManagementParams::Primaries::PRO_PHOTO:
        case ColorManagementParams::Primaries::REC2020:
        case ColorManagementParams::Primaries::ACES_P1:
        case ColorManagementParams::Primaries::WIDE_GAMUT:
        case ColorManagementParams::Primaries::ACES_P0:
        case ColorManagementParams::Primaries::JDC_MAX:
        case ColorManagementParams::Primaries::BRUCE_RGB:
        case ColorManagementParams::Primaries::BETA_RGB:
        case ColorManagementParams::Primaries::BEST_RGB: {
            labgridcie->set_sensitive(false);
            will->set_sensitive(false);
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM: {
            will->set_sensitive(true);
            labgridcie->set_sensitive(false);
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM_GRID: {
            labgridcie->set_sensitive(true);
            will->set_sensitive(false);
            break;
        }
    }

    if (ColorManagementParams::WorkingTrc(wTRC->get_active_row_number()) == ColorManagementParams::WorkingTrc::NONE) {
        redFrame->hide();
    }

    if (listener) {
        listener->panelChanged(EvICMtrcinMethod, wTRC->get_active_text());
    }
}

void ICMPanel::willChanged()
{
    switch (ColorManagementParams::Primaries(wprim->get_active_row_number() )) {
        case ColorManagementParams::Primaries::DEFAULT:
        case ColorManagementParams::Primaries::SRGB:
        case ColorManagementParams::Primaries::ADOBE_RGB:
        case ColorManagementParams::Primaries::PRO_PHOTO:
        case ColorManagementParams::Primaries::REC2020:
        case ColorManagementParams::Primaries::ACES_P1:
        case ColorManagementParams::Primaries::WIDE_GAMUT:
        case ColorManagementParams::Primaries::ACES_P0:
        case ColorManagementParams::Primaries::JDC_MAX:
        case ColorManagementParams::Primaries::BRUCE_RGB:
        case ColorManagementParams::Primaries::BETA_RGB:
        case ColorManagementParams::Primaries::BEST_RGB: {
            labgridcie->set_sensitive(false);
            will->set_sensitive(false);
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM: {
            will->set_sensitive(true);
            labgridcie->set_sensitive(false);
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM_GRID: {
            labgridcie->set_sensitive(true);
            will->set_sensitive(false);
            break;
        }
    }

    if (listener) {
        listener->panelChanged(EvICMwillMethod, will->get_active_text());
    }
}



void ICMPanel::wprimChanged()
{
    switch (ColorManagementParams::Primaries(wprim->get_active_row_number())) {
        case ColorManagementParams::Primaries::DEFAULT:
        case ColorManagementParams::Primaries::CUSTOM:
        case ColorManagementParams::Primaries::CUSTOM_GRID: {
            break;
        }

        case ColorManagementParams::Primaries::SRGB: {
            redx->setValue(0.64);
            redy->setValue(0.33);
            grex->setValue(0.30);
            grey->setValue(0.60);
            blux->setValue(0.15);
            bluy->setValue(0.06);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65));
            break;
        }

        case ColorManagementParams::Primaries::ADOBE_RGB: {
            redx->setValue(0.64);
            redy->setValue(0.33);
            grex->setValue(0.21);
            grey->setValue(0.71);
            blux->setValue(0.15);
            bluy->setValue(0.06);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65));
            break;
        }

        case ColorManagementParams::Primaries::PRO_PHOTO: {
            redx->setValue(0.7347);
            redy->setValue(0.2653);
            grex->setValue(0.1596);
            grey->setValue(0.8404);
            blux->setValue(0.0366);
            bluy->setValue(0.0001);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
            break;
        }

        case ColorManagementParams::Primaries::REC2020: {
            redx->setValue(0.708);
            redy->setValue(0.292);
            grex->setValue(0.17);
            grey->setValue(0.797);
            blux->setValue(0.131);
            bluy->setValue(0.046);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65));
            break;
        }

        case ColorManagementParams::Primaries::ACES_P1: {
            redx->setValue(0.713);
            redy->setValue(0.293);
            grex->setValue(0.165);
            grey->setValue(0.830);
            blux->setValue(0.128);
            bluy->setValue(0.044);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D60));
            break;
        }

        case ColorManagementParams::Primaries::WIDE_GAMUT: {
            redx->setValue(0.735);
            redy->setValue(0.265);
            grex->setValue(0.115);
            grey->setValue(0.826);
            blux->setValue(0.1570);
            bluy->setValue(0.018);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
            break;
        }

        case ColorManagementParams::Primaries::ACES_P0: {
            redx->setValue(0.7347);
            redy->setValue(0.2653);
            grex->setValue(0.);
            grey->setValue(1.0);
            blux->setValue(0.0001);
            bluy->setValue(-0.077);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D60));
            break;
        }

        case ColorManagementParams::Primaries::JDC_MAX: {
            redx->setValue(0.734702);
            redy->setValue(0.265302);
            grex->setValue(0.021908);
            grey->setValue(0.930288);
            blux->setValue(0.120593);
            bluy->setValue(0.001583);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
            break;
        }

        case ColorManagementParams::Primaries::BRUCE_RGB: {
            redx->setValue(0.64);
            redy->setValue(0.33);
            grex->setValue(0.28);
            grey->setValue(0.65);
            blux->setValue(0.15);
            bluy->setValue(0.06);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65));
            break;
        }

        case ColorManagementParams::Primaries::BETA_RGB: {
            redx->setValue(0.6888);
            redy->setValue(0.3112);
            grex->setValue(0.1986);
            grey->setValue(0.7551);
            blux->setValue(0.1265);
            bluy->setValue(0.0352);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
            break;
        }

        case ColorManagementParams::Primaries::BEST_RGB: {
            redx->setValue(0.7347);
            redy->setValue(0.2653);
            grex->setValue(0.2150);
            grey->setValue(0.7750);
            blux->setValue(0.130);
            bluy->setValue(0.035);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
            break;
        }
    }
    
    
    if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::DEFAULT) {
        if (wProfNames->get_active_text() == "Rec2020") {
            redx->setValue(0.708);
            redy->setValue(0.292);
            grex->setValue(0.17);
            grey->setValue(0.797);
            blux->setValue(0.131);
            bluy->setValue(0.046);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65));
        } else if (wProfNames->get_active_text() == "sRGB") {
            redx->setValue(0.64);
            redy->setValue(0.33);
            grex->setValue(0.30);
            grey->setValue(0.60);
            blux->setValue(0.15);
            bluy->setValue(0.06);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65));
        } else if (wProfNames->get_active_text() == "Adobe RGB") {
            redx->setValue(0.64);
            redy->setValue(0.33);
            grex->setValue(0.21);
            grey->setValue(0.71);
            blux->setValue(0.15);
            bluy->setValue(0.06);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65));
        } else if (wProfNames->get_active_text() == "ProPhoto") {
            redx->setValue(0.7347);
            redy->setValue(0.2653);
            grex->setValue(0.1596);
            grey->setValue(0.8404);
            blux->setValue(0.0366);
            bluy->setValue(0.0001);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
        } else if (wProfNames->get_active_text() == "ACESp1") {
            redx->setValue(0.713);
            redy->setValue(0.293);
            grex->setValue(0.165);
            grey->setValue(0.830);
            blux->setValue(0.128);
            bluy->setValue(0.044);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D60));
        } else if (wProfNames->get_active_text() == "WideGamut") {
            redx->setValue(0.735);
            redy->setValue(0.265);
            grex->setValue(0.115);
            grey->setValue(0.826);
            blux->setValue(0.1570);
            bluy->setValue(0.018);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
        } else if (wProfNames->get_active_text() == "ACESp0") {
            redx->setValue(0.7347);
            redy->setValue(0.2653);
            grex->setValue(0.);
            grey->setValue(1.0);
            blux->setValue(0.0001);
            bluy->setValue(-0.077);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D60));
        } else if (wProfNames->get_active_text() == "JDCmax") {
            redx->setValue(0.734702);
            redy->setValue(0.265302);
            grex->setValue(0.021908);
            grey->setValue(0.930288);
            blux->setValue(0.120593);
            bluy->setValue(0.001583);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
        } else if (wProfNames->get_active_text() == "BruceRGB") {
            redx->setValue(0.64);
            redy->setValue(0.33);
            grex->setValue(0.28);
            grey->setValue(0.65);
            blux->setValue(0.15);
            bluy->setValue(0.06);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D65)); 
        } else if (wProfNames->get_active_text() == "Beta RGB") {
            redx->setValue(0.6888);
            redy->setValue(0.3112);
            grex->setValue(0.1986);
            grey->setValue(0.7551);
            blux->setValue(0.1265);
            bluy->setValue(0.0352);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
        } else if (wProfNames->get_active_text() == "BestRGB") {
            redx->setValue(0.7347);
            redy->setValue(0.2653);
            grex->setValue(0.2150);
            grey->setValue(0.7750);
            blux->setValue(0.130);
            bluy->setValue(0.035);
            will->set_active(toUnderlying(ColorManagementParams::Illuminant::D50));
        }

        redFrame->hide();
    } else {
        redFrame->show();

        if (ColorManagementParams::Primaries(wprim->get_active_row_number()) != ColorManagementParams::Primaries::CUSTOM) {
            primCoordGrid->set_sensitive(false);
            labgridcie->set_sensitive(false);
            will->set_sensitive(false);
            if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::CUSTOM_GRID) {
                labgridcie->set_sensitive(true);
            }
        } else {
            primCoordGrid->set_sensitive(true);
            labgridcie->set_sensitive(false);
            will->set_sensitive(true);
        }
        
    }
    willChanged ();

    if (ColorManagementParams::Primaries(wprim->get_active_row_number()) == ColorManagementParams::Primaries::CUSTOM_GRID) {
        labgridcie->set_sensitive(true);
    } else {
        labgridcie->set_sensitive(false);
    }

    if (listener) {
        listener->panelChanged(EvICMwprimMethod, wprim->get_active_text());
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
    Glib::ustring localized_profname;

    if (inone->get_active()) {
        profname = "(none)";
        localized_profname = inone->get_label();
    } else if (iembedded->get_active()) {
        profname = "(embedded)";
        localized_profname = iembedded->get_label();
    } else if (icamera->get_active()) {
        profname = "(camera)";
        localized_profname = icamera->get_label();
    } else if (icameraICC->get_active()) {
        profname = "(cameraICC)";
        localized_profname = icameraICC->get_label();
    } else {
        profname = ipDialog->get_filename();
        localized_profname = profname;
    }

    updateDCP(-1, profname);

    if (listener && profname != oldip) {
        listener->panelChanged(EvIProfile, localized_profname);
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

void ICMPanel::aiChanged(int n)
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

        listener->panelChanged(EvaIntent, str);
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

void ICMPanel::fbwChanged()
{
    if (multiImage) {
        if (fbw->get_inconsistent()) {
            fbw->set_inconsistent(false);
            fbwconn.block(true);
            fbw->set_active(false);
            fbwconn.block(false);
        } else if (lastfbw) {
            fbw->set_inconsistent(true);
        }

        lastfbw = fbw->get_active();
    }

    if (listener) {
        if (fbw->get_inconsistent()) {
            listener->panelChanged(EvICMfbw, M("GENERAL_UNCHANGED"));
        } else if (fbw->get_active()) {
            listener->panelChanged(EvICMfbw, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvICMfbw, M("GENERAL_DISABLED"));
        }
    }
}

void ICMPanel::gamutChanged()
{
    if (multiImage) {
        if (gamut->get_inconsistent()) {
            gamut->set_inconsistent(false);
            gamutconn.block(true);
            gamut->set_active(false);
            gamutconn.block(false);
        } else if (lastgamut) {
            gamut->set_inconsistent(true);
        }

        lastgamut = gamut->get_active();
    }

    if (listener) {
        if (gamut->get_inconsistent()) {
            listener->panelChanged(EvICMgamut, M("GENERAL_UNCHANGED"));
        } else if (fbw->get_active()) {
            listener->panelChanged(EvICMgamut, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvICMgamut, M("GENERAL_DISABLED"));
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
    Gtk::Box* hbox = Gtk::manage(new Gtk::Box());
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
    aRendIntent->addEntry("template-24.png", M("GENERAL_UNCHANGED"));
    aRendIntent->show();
    wProfNames->append(M("GENERAL_UNCHANGED"));
    wTRC->append(M("GENERAL_UNCHANGED"));
    will->append(M("GENERAL_UNCHANGED"));
    wprim->append(M("GENERAL_UNCHANGED"));
    dcpIll->append(M("GENERAL_UNCHANGED"));
    wGamma->showEditedCB();
    wSlope->showEditedCB();
    redx->showEditedCB();
    redy->showEditedCB();
    grex->showEditedCB();
    grey->showEditedCB();
    blux->showEditedCB();
    bluy->showEditedCB();
    preser->showEditedCB();
}

