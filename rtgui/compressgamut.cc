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
#include "compressgamut.h"

#include "eventmapper.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Compressgamut::TOOL_NAME = "compressgamut";

Compressgamut::Compressgamut () : FoldableToolPanel(this, TOOL_NAME, M("TP_COMPRESSGAMUT_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvcgColorspace = m->newEvent(WB, "HISTORY_MSG_CG_COLORSPACE");
    Evcgthc = m->newEvent(WB, "HISTORY_MSG_CG_CYANTH");
    Evcgthm = m->newEvent(WB, "HISTORY_MSG_CG_MAGENTATH");
    Evcgthy = m->newEvent(WB, "HISTORY_MSG_CG_YELLOWTH");
    Evcgdc = m->newEvent(WB, "HISTORY_MSG_CG_CYANLIM");
    Evcgdm = m->newEvent(WB, "HISTORY_MSG_CG_MAGENTALIM");
    Evcgdy = m->newEvent(WB, "HISTORY_MSG_CG_YELLOWLIM");
    Evcgroll = m->newEvent(WB, "HISTORY_MSG_CG_ROLLOFF");
    Evcgpwr = m->newEvent(WB, "HISTORY_MSG_CG_VALUE");
    Evcgenabled = m->newEvent(WB, "HISTORY_MSG_CG_ENABLED");


    Gtk::Frame *iFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_MAIN_COLORSPACE")));

    iFrame->set_label_align(0.025, 0.5);
    iFrame->set_tooltip_markup (M("PREFERENCES_WBACORR_TOOLTIP"));

    ToolParamBlock* const iBox = Gtk::manage(new ToolParamBlock());
    
    colorspace = Gtk::manage(new MyComboBoxText());
    colorspace->append(M("TP_COMPRESSGAMUT_REC2020"));
    colorspace->append(M("TP_COMPRESSGAMUT_PROPHOTO"));
    colorspace->append(M("TP_COMPRESSGAMUT_SRGB"));
    colorspace->append(M("TP_COMPRESSGAMUT_DCIP3"));
    colorspace->append(M("TP_COMPRESSGAMUT_ACESP1"));
    iBox->pack_start(*colorspace);
    iFrame->add(*iBox);
    pack_start(*iFrame);
    colorspaceconn = colorspace->signal_changed().connect(sigc::mem_fun(*this, &Compressgamut::colorspaceChanged));


    th_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANTH"), 0., 0.999, 0.001, 0.815));
    th_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTATH"), 0., 0.999, 0.001, 0.803));
    th_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWTH"), 0., 0.999, 0.001, 0.880));

    Gtk::Frame *thFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_THRESHOLD")));
    thFrame->set_label_align(0.025, 0.5);
    Gtk::Box *thVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    thVBox->pack_start (*th_c, Gtk::PACK_SHRINK);
    thVBox->pack_start (*th_m, Gtk::PACK_SHRINK);
    thVBox->pack_start (*th_y, Gtk::PACK_SHRINK);
    thFrame->add(*thVBox);
    pack_start(*thFrame, Gtk::PACK_SHRINK);//, Gtk::PACK_EXPAND_WIDGET);

    d_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANLIM"), 1.001, 2.0, 0.001, 1.147));
    d_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTALIM"), 1.001, 2.0, 0.001, 1.264));
    d_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWLIM"), 1.001, 2.0, 0.001, 1.312));

    Gtk::Frame *limFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_LIMIT")));
    limFrame->set_label_align(0.025, 0.5);
    Gtk::Box *limVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    limVBox->pack_start (*d_c, Gtk::PACK_SHRINK);
    limVBox->pack_start (*d_m, Gtk::PACK_SHRINK);
    limVBox->pack_start (*d_y, Gtk::PACK_SHRINK);
    limFrame->add(*limVBox);
    pack_start(*limFrame, Gtk::PACK_SHRINK);//, Gtk::PACK_EXPAND_WIDGET);


    rolloff = Gtk::manage(new Gtk::CheckButton(M("TP_COMPRESSGAMUT_ROLLOFF")));
    pwr = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_PWR"), 0.5, 2.0, 0.01, 1.2));
    rolloffconn = rolloff->signal_pressed().connect ( sigc::mem_fun (*this, &Compressgamut::rolloff_change) );

    Gtk::Frame *rollFrame = Gtk::manage(new Gtk::Frame());
    rollFrame->set_label_align(0.025, 0.5);
    rollFrame->set_label_widget(*rolloff);
    Gtk::Box *rollVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    
    rollVBox->pack_start (*pwr, Gtk::PACK_SHRINK);
    rollFrame->add(*rollVBox);
    pack_start(*rollFrame, Gtk::PACK_SHRINK);//, Gtk::PACK_EXPAND_WIDGET);
    th_c->setAdjusterListener (this);
    th_m->setAdjusterListener (this);
    th_y->setAdjusterListener (this);
    d_c->setAdjusterListener (this);
    d_m->setAdjusterListener (this);
    d_y->setAdjusterListener (this);
    pwr->setAdjusterListener (this);
    

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

    }

    setEnabled (pp->cg.enabled);
    
    
    rolloffconn.block (true);
    rolloff->set_active (pp->cg.rolloff);
    rolloffconn.block (false);
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
     } else if (pp->cg.colorspace == "srgb") {
        colorspace->set_active(2);
    } else if (pp->cg.colorspace == "dcip3") {
        colorspace->set_active(3);
    } else if (pp->cg.colorspace == "acesp1") {
        colorspace->set_active(4);
    }
    colorspaceconn.block (false);
    rolloff_change();
    colorspaceChanged();
    enabledChanged ();
    enableListener ();
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
    pp->cg.enabled = getEnabled();

    if (colorspace->get_active_row_number() == 0) {
        pp->cg.colorspace = "rec2020";
    } else if (colorspace->get_active_row_number() == 1){
        pp->cg.colorspace = "prophoto";
    } else if (colorspace->get_active_row_number() == 2){
        pp->cg.colorspace = "srgb";
    } else if (colorspace->get_active_row_number() == 3){
        pp->cg.colorspace = "dcip3";
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
        pedited->cg.colorspace = colorspace->get_active_row_number() != 5;
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

void Compressgamut::rolloff_change()
{
    if (rolloff->get_active()) {
        pwr->set_sensitive(true);
    } else {
        pwr->set_sensitive(false);
    }

    if (multiImage) {
        if (rolloff->get_inconsistent()) {
            rolloff->set_inconsistent(false);
            rolloffconn.block(true);
            rolloff->set_active(false);
            rolloffconn.block(false);
        } else if (lastrolloff) {
            rolloff->set_inconsistent(true);
        }

        lastrolloff = rolloff->get_active();
    }
    
    
    if (listener) {
        if (rolloff->get_inconsistent()) {
            listener->panelChanged(Evcgroll, M("GENERAL_UNCHANGED"));
        } else if (rolloff->get_active()) {
            listener->panelChanged(Evcgroll, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evcgroll, M("GENERAL_DISABLED"));
        }
    }
    
}

void Compressgamut::colorspaceChanged()
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvcgColorspace, M("GENERAL_ENABLED"));
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
