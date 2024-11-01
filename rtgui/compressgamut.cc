/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2024 RawTherapee team
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
 /*
 * tweaked from the original from https://github.com/jedypod/gamut-compress
 * https://docs.acescentral.com/specifications/rgc/
*/ 
#include "compressgamut.h"

#include "eventmapper.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Compressgamut::TOOL_NAME = "compressgamut";

Compressgamut::Compressgamut () : FoldableToolPanel(this, TOOL_NAME, M("TP_COMPRESSGAMUT_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvcgColorspace = m->newEvent(COMPR, "HISTORY_MSG_CG_COLORSPACE");//FIRST to reinit histogram
    Evcgthc = m->newEvent(COMPR, "HISTORY_MSG_CG_CYANTH");
    Evcgthm = m->newEvent(COMPR, "HISTORY_MSG_CG_MAGENTATH");
    Evcgthy = m->newEvent(COMPR, "HISTORY_MSG_CG_YELLOWTH");
    Evcgdc = m->newEvent(COMPR, "HISTORY_MSG_CG_CYANLIM");
    Evcgdm = m->newEvent(COMPR, "HISTORY_MSG_CG_MAGENTALIM");
    Evcgdy = m->newEvent(COMPR, "HISTORY_MSG_CG_YELLOWLIM");
    Evcgroll = m->newEvent(COMPR, "HISTORY_MSG_CG_ROLLOFF");
    Evcgpwr = m->newEvent(COMPR, "HISTORY_MSG_CG_VALUE");
    Evcgenabled = m->newEvent(COMPR, "HISTORY_MSG_CG_ENABLED");
    Evcgkeepset = m->newEvent(COMPR, "HISTORY_MSG_CG_KEEPSET");



    Gtk::Frame *iFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_MAIN_COLORSPACE")));

    iFrame->set_label_align(0.025f, 0.5);
    iFrame->set_tooltip_markup (M("TP_COMPRESSGAMUT_COLORSPACE_TOOLTIP"));

    Gtk::Box *iVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
   
    colorspace = Gtk::manage(new MyComboBoxText());
    colorspace->append(M("TP_COMPRESSGAMUT_REC2020"));
    colorspace->append(M("TP_COMPRESSGAMUT_PROPHOTO"));
    colorspace->append(M("TP_COMPRESSGAMUT_ADOBE"));
    colorspace->append(M("TP_COMPRESSGAMUT_SRGB"));
  //  colorspace->append(M("TP_COMPRESSGAMUT_DCIP3"));
    colorspace->append(M("TP_COMPRESSGAMUT_ACESP1"));
    colorspace->set_active(3);

    //keep settings when changing target gamut workspace
    keepset = Gtk::manage(new Gtk::CheckButton(M("TP_COMPRESSGAMUT_KEEPSET")));
    keepsetconn = keepset->signal_toggled().connect (sigc::mem_fun (*this, &Compressgamut::keepset_change));
    keepset->set_active(false); 

    iVBox->pack_start(*colorspace);
    iVBox->pack_start(*keepset);
    iFrame->add(*iVBox);
    pack_start(*iFrame);
    colorspaceconn = colorspace->signal_changed().connect(sigc::mem_fun(*this, &Compressgamut::colorspaceChanged));

    // Percentage of the core gamut to protect Limits
    // Values calculated to protect all the colors of the ColorChecker Classic 24 as given by
    // ISO 17321-1 and Ohta (1997)

    //others values calculated (estimated) with my color chart Jacques Desmis (468 colors) a "super" Colorchecker and CIExy distance between white point and borders 
    //and others difficult images : flowers, submarine

    //max th is important to avoid artifacts and segmentation fault - with some images max th = 1 = segmentation fault or artifacts
    th_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANTH"), 0., 1., 0.001, 0.815));//0.25 sRGB 0.999 to avoid: 1 - th = 0
    th_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTATH"), 0., 1., 0.001, 0.803));//0.925 sRGB
    th_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWTH"), 0., 1., 0.001, 0.880));//0.934 sRGB
    //see others values for Target workspace in Procparams.cc and in updategamutGUI
    Gtk::Frame *thFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_THRESHOLD")));
    thFrame->set_label_align(0.025f, 0.5);
    Gtk::Box *thVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    thFrame->set_tooltip_markup (M("TP_COMPRESSGAMUT_THRESHOLD_TOOLTIP"));

    thVBox->pack_start (*th_c);
    thVBox->pack_start (*th_m);
    thVBox->pack_start (*th_y);
    thFrame->add(*thVBox);
    pack_start(*thFrame, Gtk::PACK_SHRINK);
    // from ACES https://docs.acescentral.com/specifications/rgc/#appendix-c-illustrations
    // https://docs.acescentral.com/specifications/rgc/#appendix-d-ctl-reference-implementation
    // Distance from achromatic which will be compressed to the gamut boundary
    // Values calculated to encompass the encoding gamuts of common digital cinema cameras
    d_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANLIM"), 1.001, 2.0, 0.001, 1.147));//1.05 sRGB
    d_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTALIM"), 1.001, 2.0, 0.001, 1.264));//1.08 sRGB
    d_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWLIM"), 1.001, 2.0, 0.001, 1.312));//1.10 sRGB

    //see others values for Target workspace in Procparams.cc and in updategamutGUI
    //made by estimation using my color chart (468 colors) a "super" Colorchecker
    //and others difficult images : flowers, submarine
    // and the CIExy diagram - position of the white point relative to the 3 edges of the triangle cyan, magenta, yellow
    // of course to refine by testing 

    Gtk::Frame *limFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_LIMIT")));
    limFrame->set_label_align(0.025f, 0.5);
    limFrame->set_tooltip_markup (M("TP_COMPRESSGAMUT_LIMIT_TOOLTIP"));

    Gtk::Box *limVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    limVBox->pack_start (*d_c);
    limVBox->pack_start (*d_m);
    limVBox->pack_start (*d_y);
    limFrame->add(*limVBox);
    pack_start(*limFrame, Gtk::PACK_SHRINK);

    // Aggressiveness of the compression curve Pwr
    rolloff = Gtk::manage(new Gtk::CheckButton(M("TP_COMPRESSGAMUT_ROLLOFF")));
    pwr = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_PWR"), 0.5, 2.0, 0.01, 1.2));
    rolloffconn = rolloff->signal_toggled().connect (sigc::mem_fun (*this, &Compressgamut::rolloff_change));

    Gtk::Frame *rollFrame = Gtk::manage(new Gtk::Frame());
    rollFrame->set_label_align(0.025f, 0.5);
    rollFrame->set_label_widget(*rolloff);
    rolloff->set_active(true); 
    Gtk::Box *rollVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    rollFrame->set_tooltip_markup (M("TP_COMPRESSGAMUT_POWER_TOOLTIP"));

    rollVBox->pack_start (*pwr);
    rollFrame->add(*rollVBox);
    pack_start(*rollFrame, Gtk::PACK_SHRINK);
    th_c->setAdjusterListener (this);
    th_m->setAdjusterListener (this);
    th_y->setAdjusterListener (this);
    d_c->setAdjusterListener (this);
    d_m->setAdjusterListener (this);
    d_y->setAdjusterListener (this);
    pwr->setAdjusterListener (this);
    
    show_all_children ();
}

void Compressgamut::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
   
    if (pedited) {
        th_c->setEditedState       (pedited->cg.th_c ? Edited : UnEdited);
        th_m->setEditedState       (pedited->cg.th_m ? Edited : UnEdited);
        th_y->setEditedState       (pedited->cg.th_y ? Edited : UnEdited);
        d_c->setEditedState        (pedited->cg.d_c ? Edited : UnEdited);
        d_m->setEditedState        (pedited->cg.d_m ? Edited : UnEdited);
        d_y->setEditedState        (pedited->cg.d_y ? Edited : UnEdited);
        pwr->setEditedState        (pedited->cg.pwr ? Edited : UnEdited);
        set_inconsistent           (multiImage && !pedited->cg.enabled);
        rolloff->set_inconsistent  (!pedited->cg.rolloff);
        keepset->set_inconsistent  (!pedited->cg.keepset);
    }

    setEnabled (pp->cg.enabled);

    keepsetconn.block (true);
    keepset->set_active (pp->cg.keepset);
    keepsetconn.block (false);

    rolloffconn.block (true);
    rolloff->set_active (pp->cg.rolloff);
    rolloffconn.block (false);
    updategamutGUI();    
    
    th_c->setValue(pp->cg.th_c);
    th_m->setValue(pp->cg.th_m);
    th_y->setValue(pp->cg.th_y);
    d_c->setValue(pp->cg.d_c);
    d_m->setValue(pp->cg.d_m);
    d_y->setValue(pp->cg.d_y);
    pwr->setValue(pp->cg.pwr);

    colorspaceconn.block (true);

    if (pp->cg.colorspace == "rec2020") {
        colorspace->set_active(0);
    } else if (pp->cg.colorspace == "prophoto") {
        colorspace->set_active(1);
    } else if (pp->cg.colorspace == "adobe") {
        colorspace->set_active(2);
     } else if (pp->cg.colorspace == "srgb") {
        colorspace->set_active(3);
//    } else if (pp->cg.colorspace == "dcip3") {
//        colorspace->set_active(4);
    } else if (pp->cg.colorspace == "acesp1") {
        colorspace->set_active(4);
    }
    colorspaceconn.block (false);

    keepsetconn.block (true);
    keepset->set_active (pp->cg.keepset);
    keepsetconn.block (false);

    lastkeepset = pp->cg.keepset;
    keepset_change();


    rolloffconn.block (true);
    rolloff->set_active (pp->cg.rolloff);
    rolloffconn.block (false);

    lastrolloff = pp->cg.rolloff;
    rolloff_change();
    colorspaceChanged();
    enabledChanged ();

    enableListener ();
}

void Compressgamut::updategamutGUI()
{
    // Update default slider value GUI according to colorspace
    //new last factor vDef1 vDef2 try to take into account the colorspace gamut 
    // th_xx->setLimits(0., 0.999, 0.001, vDef1); I have modified (a little) adjuster.cc (I hope no border effects)
    // d_xx->setLimits(1.001, 2.0, 0.001, vDef2);
    //made by estimation using my color chart (468 colors) a "super" Colorchecker
    //and others difficult images : flowers, submarine
    // and the CIExy diagram - position of the white point relative to the 3 edges of the triangle cyan, magenta, yellow
    // of course to refine by testing 
    //save values in case of 
    const double temp_tc = th_c->getValue();
    const double temp_tm = th_m->getValue();
    const double temp_ty = th_y->getValue();
    const double temp_dc = d_c->getValue();
    const double temp_dm = d_m->getValue();
    const double temp_dy = d_y->getValue();
    double maxth = 1.;// 1.0 (or 0.9999) leeds in some cases to artifacts or segmentation fault - I change in color.cc aces_reference_gamut_compression by limiting theshold to 0.999
    
     if (colorspace->get_active_row_number() == 0) {//rec2020
        th_c->setLimits(0., maxth, 0.001, 0.71);
        th_m->setLimits(0., maxth, 0.001, 0.803);
        th_y->setLimits(0., maxth, 0.001, 0.870);
        d_c->setLimits(1.001, 2., 0.001, 1.12);
        d_m->setLimits(1.001, 2., 0.001, 1.26);
        d_y->setLimits(1.001, 2., 0.001, 1.31);       
    } else if (colorspace->get_active_row_number() == 1){//prophoto
        th_c->setLimits(0., maxth, 0.001, 0.85);
        th_m->setLimits(0., maxth, 0.001, 0.79);
        th_y->setLimits(0., maxth, 0.001, 0.90);
        d_c->setLimits(1.001, 2., 0.001, 1.17);
        d_m->setLimits(1.001, 2., 0.001, 1.27);
        d_y->setLimits(1.001, 2., 0.001, 1.35);
    } else if (colorspace->get_active_row_number() == 2){//Adobe
        th_c->setLimits(0., maxth, 0.001, 0.45);
        th_m->setLimits(0., maxth, 0.001, 0.86);
        th_y->setLimits(0., maxth, 0.001, 0.92);
        d_c->setLimits(1.001, 2., 0.001, 1.09);
        d_m->setLimits(1.001, 2., 0.001, 1.20);
        d_y->setLimits(1.001, 2., 0.001, 1.09);
    } else if (colorspace->get_active_row_number() == 3){//srgb
        th_c->setLimits(0., maxth, 0.001, 0.25);
        th_m->setLimits(0., maxth, 0.001, 0.90);
        th_y->setLimits(0., maxth, 0.001, 0.934);
        d_c->setLimits(1.001, 2., 0.001, 1.05);
        d_m->setLimits(1.001, 2., 0.001, 1.08);
        d_y->setLimits(1.001, 2., 0.001, 1.10);
//    } else if (colorspace->get_active_row_number() == 4){//dci-p3
//        th_c->setLimits(0., maxth, 0.001, 0.40);
//        th_m->setLimits(0., maxth, 0.001, 0.91);
//        th_y->setLimits(0., maxth, 0.001, 0.916);
//        d_c->setLimits(1.001, 2., 0.001, 1.08);
//        d_m->setLimits(1.001, 2., 0.001, 1.18);
//        d_y->setLimits(1.001, 2., 0.001, 1.26);
    } else if (colorspace->get_active_row_number() == 4){//acesp1
        th_c->setLimits(0., maxth, 0.001, 0.815);
        th_m->setLimits(0., maxth, 0.001, 0.803);
        th_y->setLimits(0., maxth, 0.001, 0.880);
        d_c->setLimits(1.001, 2., 0.001, 1.147);
        d_m->setLimits(1.001, 2., 0.001, 1.264);
        d_y->setLimits(1.001, 2., 0.001, 1.312);
    }
    //restore values
    if(keepset->get_active()){
        th_c->setValue(temp_tc);
        th_m->setValue(temp_tm);
        th_y->setValue(temp_ty);
        d_c->setValue(temp_dc);
        d_m->setValue(temp_dm);
        d_y->setValue(temp_dy);
    }

}


void Compressgamut::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->cg.th_c = th_c->getValue ();
    pp->cg.th_m = th_m->getValue ();
    pp->cg.th_y = th_y->getValue ();
    pp->cg.d_c = d_c->getValue ();
    pp->cg.d_m = d_m->getValue ();
    pp->cg.d_y = d_y->getValue ();
    pp->cg.pwr = pwr->getValue ();
    pp->cg.rolloff = rolloff->get_active();
    pp->cg.keepset = keepset->get_active();
    pp->cg.enabled = getEnabled();

    if (colorspace->get_active_row_number() == 0) {
        pp->cg.colorspace = "rec2020";
    } else if (colorspace->get_active_row_number() == 1){
        pp->cg.colorspace = "prophoto";
    } else if (colorspace->get_active_row_number() == 2){
        pp->cg.colorspace = "adobe";
    } else if (colorspace->get_active_row_number() == 3){
        pp->cg.colorspace = "srgb";
 //   } else if (colorspace->get_active_row_number() == 4){
 //       pp->cg.colorspace = "dcip3";
    } else if (colorspace->get_active_row_number() == 4){
        pp->cg.colorspace = "acesp1";
    }


    if (pedited) {
        pedited->cg.th_c          = th_c->getEditedState ();
        pedited->cg.th_m          = th_m->getEditedState ();
        pedited->cg.th_y          = th_y->getEditedState ();
        pedited->cg.d_c           = d_c->getEditedState ();
        pedited->cg.d_m           = d_m->getEditedState ();
        pedited->cg.d_y           = d_y->getEditedState ();
        pedited->cg.pwr           = pwr->getEditedState ();
        pedited->cg.enabled       = !get_inconsistent();
        pedited->cg.colorspace = colorspace->get_active_row_number() != 6;
        pedited->cg.rolloff       = !rolloff->get_inconsistent();
        pedited->cg.keepset       = !keepset->get_inconsistent();
    }

}

void Compressgamut::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    th_c->setDefault (defParams->cg.th_c);
    th_m->setDefault (defParams->cg.th_m);
    th_y->setDefault (defParams->cg.th_y);
    d_c->setDefault (defParams->cg.d_c);
    d_m->setDefault (defParams->cg.d_m);
    d_y->setDefault (defParams->cg.d_y);
    pwr->setDefault (defParams->cg.pwr);

    if (pedited) {
        th_c->setDefaultEditedState       (pedited->cg.th_c ? Edited : UnEdited);
        th_m->setDefaultEditedState       (pedited->cg.th_m ? Edited : UnEdited);
        th_y->setDefaultEditedState       (pedited->cg.th_y ? Edited : UnEdited);
        d_c->setDefaultEditedState        (pedited->cg.d_c ? Edited : UnEdited);
        d_m->setDefaultEditedState        (pedited->cg.d_m ? Edited : UnEdited);
        d_y->setDefaultEditedState        (pedited->cg.d_y ? Edited : UnEdited);
        pwr->setDefaultEditedState        (pedited->cg.pwr ? Edited : UnEdited);
        
    } else {
        th_c->setDefaultEditedState       (Irrelevant);
        th_m->setDefaultEditedState       (Irrelevant);
        th_y->setDefaultEditedState       (Irrelevant);
        d_c->setDefaultEditedState       (Irrelevant);
        d_m->setDefaultEditedState       (Irrelevant);
        d_y->setDefaultEditedState       (Irrelevant);
        pwr->setDefaultEditedState       (Irrelevant);
    }

}

void Compressgamut::adjusterChanged (Adjuster* a, double newval)
{
    
    if (listener && getEnabled()) {
        const Glib::ustring costr = Glib::ustring::format (a->getValue());

        if (a == th_c) {
            listener->panelChanged (Evcgthc, costr);
        } else if (a == th_m) {
            listener->panelChanged (Evcgthm, costr);
        } else if (a == th_y) {
            listener->panelChanged (Evcgthy, costr);
        } else if (a == d_c) {
            listener->panelChanged (Evcgdc, costr);
        } else if (a == d_m) {
            listener->panelChanged (Evcgdm, costr);
        } else if (a == d_y) {
            listener->panelChanged (Evcgdy, costr);
        } else if (a == pwr) {
            listener->panelChanged (Evcgpwr, costr);
        }
    }
   
}

void Compressgamut::keepset_change()
{

    if (batchMode) {
        if (keepset->get_inconsistent()) {
            keepset->set_inconsistent (false);
            keepsetconn.block (true);
            keepset->set_active (false);
            keepsetconn.block (false);
        } else if (lastkeepset) {
            keepset->set_inconsistent (true);
        }

        lastkeepset = keepset->get_active ();
    }

    if (listener) {
        if (keepset->get_active()) {
            listener->panelChanged(Evcgkeepset, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evcgkeepset, M("GENERAL_DISABLED"));
        }
    }
}


void Compressgamut::rolloff_change()
{
    if (rolloff->get_active()) {
        pwr->set_sensitive(true);
    } else {
        pwr->set_sensitive(false);
    }

    if (batchMode) {
        if (rolloff->get_inconsistent()) {
            rolloff->set_inconsistent (false);
            rolloffconn.block (true);
            rolloff->set_active (false);
            rolloffconn.block (false);
        } else if (lastrolloff) {
            rolloff->set_inconsistent (true);
        }

        lastrolloff = rolloff->get_active ();
    }

    if (listener) {
        if (rolloff->get_active()) {
            listener->panelChanged(Evcgroll, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evcgroll, M("GENERAL_DISABLED"));
        }
    }
}

void Compressgamut::colorspaceChanged()
{
    updategamutGUI();
    if (listener && getEnabled()) {
        listener->panelChanged(EvcgColorspace, colorspace->get_active_text());
    }
}

void Compressgamut::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (Evcgenabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (Evcgenabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evcgenabled, M("GENERAL_DISABLED"));
        }
    }
}
void Compressgamut::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    colorspace->append (M("GENERAL_UNCHANGED"));
    th_c->showEditedCB ();
    th_m->showEditedCB ();
    th_y->showEditedCB ();
    d_c->showEditedCB ();
    d_m->showEditedCB ();
    d_y->showEditedCB ();
    pwr->showEditedCB ();
}


void Compressgamut::trimValues (rtengine::procparams::ProcParams* pp)
{
    th_c->trimValue(pp->cg.th_c);
    th_m->trimValue(pp->cg.th_m);
    th_y->trimValue(pp->cg.th_y);
    d_c->trimValue(pp->cg.d_c);
    d_m->trimValue(pp->cg.d_m);
    d_y->trimValue(pp->cg.d_y);
    pwr->trimValue(pp->cg.pwr);
}
