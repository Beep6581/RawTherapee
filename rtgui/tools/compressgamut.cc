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

 //Jacques Desmis - December 2025
*/ 
#include "compressgamut.h"

#include "eventmapper.h"
#include <iomanip>
#include "rtengine/procparams.h"


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
    Evcgdcautoon = m->newEvent(COMPR, "HISTORY_MSG_CG_CYANDC_AUTO");
    Evcgdmautoon = m->newEvent(COMPR, "HISTORY_MSG_CG_CYANDM_AUTO");
    Evcgdyautoon = m->newEvent(COMPR, "HISTORY_MSG_CG_CYANDY_AUTO");
    Gtk::Frame *iFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_MAIN_COLORSPACE")));

    iFrame->set_label_align(0.025f, 0.5);
    iFrame->set_tooltip_markup (M("TP_COMPRESSGAMUT_COLORSPACE_TOOLTIP"));

    Gtk::Box *iVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
   
    colorspace = Gtk::manage(new MyComboBoxText());
    colorspace->append(M("TP_COMPRESSGAMUT_REC2020"));
    colorspace->append(M("TP_COMPRESSGAMUT_PROPHOTO"));
    colorspace->append(M("TP_COMPRESSGAMUT_ADOBE"));
    colorspace->append(M("TP_COMPRESSGAMUT_SRGB"));
    colorspace->append(M("TP_COMPRESSGAMUT_DCIP3"));
    colorspace->append(M("TP_COMPRESSGAMUT_ACESP1"));
    colorspace->append(M("TP_COMPRESSGAMUT_BETA"));
    colorspace->set_active(3);
    iVBox->pack_start(*colorspace);
    iFrame->add(*iVBox);
    pack_start(*iFrame);
    colorspaceconn = colorspace->signal_changed().connect(sigc::mem_fun(*this, &Compressgamut::colorspaceChanged));
 
    acLabel = Gtk::manage (new Gtk::Label ("---"));
    setExpandAlignProperties (acLabel, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
    acLabel->set_tooltip_markup (M ("TP_COMPRESSGAMUT_MACLABEL_TOOLTIP"));
    acLabel->show ();

    acLabelrgb = Gtk::manage (new Gtk::Label ("---"));
    setExpandAlignProperties (acLabelrgb, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
    acLabelrgb->set_tooltip_markup (M ("TP_COMPRESSGAMUT_MACLABELRGB_TOOLTIP"));
    acLabelrgb->show (); 

    acLabelcmy = Gtk::manage (new Gtk::Label ("---"));
    setExpandAlignProperties (acLabelcmy, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
    acLabelcmy->set_tooltip_markup (M ("TP_COMPRESSGAMUT_MACLABELRGB_TOOLTIP"));
    acLabelcmy->show ();

    // Percentage of the core gamut to protect Limits
    // Values calculated to protect all the colors of the ColorChecker Classic 24 as given by
    // ISO 17321-1 and Ohta (1997)

    th_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANTH"), 0., 1., 0.001, 0.815));//0.25 sRGB 0.999 to avoid: 1 - th = 0
    th_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTATH"), 0., 1., 0.001, 0.803));//0.925 sRGB
    th_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWTH"), 0., 1., 0.001, 0.880));//0.934 sRGB

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
    d_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANLIM"), 1.001, 10., 0.001, 1.147));//1.05 sRGB
    d_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTALIM"), 1.001, 10., 0.001, 1.264));//1.08 sRGB
    d_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWLIM"), 1.001, 10., 0.001, 1.312));//1.10 sRGB

    Gtk::Frame *limFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_LIMIT")));
    limFrame->set_label_align(0.025f, 0.5);
    limFrame->set_tooltip_markup (M("TP_COMPRESSGAMUT_LIMIT_TOOLTIP"));

    Gtk::Box *limVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    limVBox->pack_start (*d_c);
    limVBox->pack_start (*d_m);
    limVBox->pack_start (*d_y);
    d_c->addAutoButton();
    d_c->setAutoValue(true);
    d_m->addAutoButton();
    d_m->setAutoValue(true);
    d_y->addAutoButton();
    d_y->setAutoValue(true);
    limVBox->pack_start (*acLabel);
    limVBox->pack_start (*acLabelrgb);
    limVBox->pack_start (*acLabelcmy);
    limFrame->add(*limVBox);
    pack_start(*limFrame, Gtk::PACK_SHRINK);

    // Aggressiveness of the compression curve Pwr
    rolloff = Gtk::manage(new Gtk::CheckButton(M("TP_COMPRESSGAMUT_ROLLOFF")));
    pwr = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_PWR"), 0.2, 2.0, 0.01, 1.2));
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
    d_c->setLogScale(100, 1);
    d_m->setLogScale(100, 1);
    d_y->setLogScale(100, 1);
    show_all_children ();
}

Compressgamut::~Compressgamut()
{
    idle_register.destroy();
}

void Compressgamut::achromaticChanged (double acmax, double acmax0, double acmax1, double acmax2, bool auto_dc, bool auto_dm, bool auto_dy)
{


    idle_register.add(
         [this, acmax, acmax0, acmax1, acmax2, auto_dc, auto_dm, auto_dy]() -> bool

        {
            GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

            disableListener();
            acLabel->set_text(
                Glib::ustring::compose(M("TP_COMPRESSGAMUT_MACLABEL"),

                                        Glib::ustring::format(std::fixed, std::setprecision(2), acmax))//achromatic maximum value
            );
            acLabelrgb->set_text (
                Glib::ustring::compose (M ("TP_COMPRESSGAMUT_MACLABELRGB"),
                                        Glib::ustring::format (std::fixed, std::setprecision (2), acmax0),//red
                                        Glib::ustring::format (std::fixed, std::setprecision (2), acmax1),//green
                                        Glib::ustring::format (std::fixed, std::setprecision (2), acmax2))//blue
            );
            acLabelcmy->set_text (
                Glib::ustring::compose (M ("TP_COMPRESSGAMUT_MACLABELCMY"),
                                        Glib::ustring::format (std::fixed, std::setprecision (1), (acmax1 + acmax2) * 0.43),//estimated value for Cyan - about 2 * (acmax1 + acmax2) * Sin (PI/3)
                                        Glib::ustring::format (std::fixed, std::setprecision (1), (acmax0 + acmax2) * 0.43),//estimated value for Magenta
                                        Glib::ustring::format (std::fixed, std::setprecision (1), (acmax0 + acmax1) * 0.43))//estimated value for Yellow
            );
            if(auto_dc) {
                double valdc = max(1.02, (acmax1 + acmax2) * 0.4);//about 90% max
                d_c->setValue(valdc);
            }
            if(auto_dm) {
                double valdm = max(1.02, (acmax0 + acmax2) * 0.4);//about 90% max
                d_m->setValue(valdm);
            }
            if(auto_dy) {
                double valdy = max(1.02, (acmax0 + acmax1) * 0.4);//about 90% max
                d_y->setValue(valdy);
            }                       
            enableListener();
            return false;
        }
    );

}


void Compressgamut::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
   
    if (pedited) {
        th_c->setEditedState       (pedited->cg.th_c ? Edited : UnEdited);
        th_m->setEditedState       (pedited->cg.th_m ? Edited : UnEdited);
        th_y->setEditedState       (pedited->cg.th_y ? Edited : UnEdited);
        d_c->setEditedState        (pedited->cg.d_c ? Edited : UnEdited);
        d_c->setAutoInconsistent   (multiImage && !pedited->cg.autodc);
        d_m->setEditedState        (pedited->cg.d_m ? Edited : UnEdited);
        d_m->setAutoInconsistent   (multiImage && !pedited->cg.autodm);
        d_y->setEditedState        (pedited->cg.d_y ? Edited : UnEdited);
        d_y->setAutoInconsistent   (multiImage && !pedited->cg.autody);
        pwr->setEditedState        (pedited->cg.pwr ? Edited : UnEdited);
        set_inconsistent           (multiImage && !pedited->cg.enabled);
        rolloff->set_inconsistent  (!pedited->cg.rolloff);
    }

    setEnabled (pp->cg.enabled);

    rolloffconn.block (true);
    rolloff->set_active (pp->cg.rolloff);
    rolloffconn.block (false);

    colorspaceconn.block (true);

    if (pp->cg.colorspace == "rec2020") {
        colorspace->set_active(0);
    } else if (pp->cg.colorspace == "prophoto") {
        colorspace->set_active(1);
    } else if (pp->cg.colorspace == "adobe") {
        colorspace->set_active(2);
     } else if (pp->cg.colorspace == "srgb") {
        colorspace->set_active(3);
    } else if (pp->cg.colorspace == "dcip3") {
        colorspace->set_active(4);
    } else if (pp->cg.colorspace == "acesp1") {
        colorspace->set_active(5);
    } else if (pp->cg.colorspace == "beta") {
        colorspace->set_active(6);    
    }
    colorspaceconn.block (false);


    updategamutGUI(); 


    th_c->setValue(pp->cg.th_c);
    th_m->setValue(pp->cg.th_m);
    th_y->setValue(pp->cg.th_y);
    d_c->setValue(pp->cg.d_c);
    d_c->setAutoValue(pp->cg.autodc);
    d_m->setValue(pp->cg.d_m);
    d_m->setAutoValue(pp->cg.autodm);
    d_y->setValue(pp->cg.d_y);
    d_y->setAutoValue(pp->cg.autody);
    pwr->setValue(pp->cg.pwr);
    
    rolloffconn.block (true);
    rolloff->set_active (pp->cg.rolloff);
    rolloffconn.block (false);

 

    lastrolloff = pp->cg.rolloff;
    lastAutodc = pp->cg.autodc;
    lastAutodm = pp->cg.autodm;
    lastAutody = pp->cg.autody;  
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
    pp->cg.autodc = d_c->getAutoValue();
    pp->cg.d_m = d_m->getValue ();
    pp->cg.autodm = d_m->getAutoValue();
    pp->cg.d_y = d_y->getValue ();
    pp->cg.autody = d_y->getAutoValue();
    pp->cg.pwr = pwr->getValue ();
    pp->cg.rolloff = rolloff->get_active();
    pp->cg.enabled = getEnabled();

    if (colorspace->get_active_row_number() == 0) {
        pp->cg.colorspace = "rec2020";
    } else if (colorspace->get_active_row_number() == 1){
        pp->cg.colorspace = "prophoto";
    } else if (colorspace->get_active_row_number() == 2){
        pp->cg.colorspace = "adobe";
    } else if (colorspace->get_active_row_number() == 3){
        pp->cg.colorspace = "srgb";
    } else if (colorspace->get_active_row_number() == 4){
        pp->cg.colorspace = "dcip3";
    } else if (colorspace->get_active_row_number() == 5){
        pp->cg.colorspace = "acesp1";
    } else if (colorspace->get_active_row_number() == 6){
        pp->cg.colorspace = "beta";
    }


    if (pedited) {
        pedited->cg.th_c          = th_c->getEditedState ();
      
        pedited->cg.th_m          = th_m->getEditedState ();
        pedited->cg.th_y          = th_y->getEditedState ();
        pedited->cg.d_c           = d_c->getEditedState ();
        pedited->cg.autodc        = !d_c->getAutoInconsistent();
        pedited->cg.d_m           = d_m->getEditedState ();
        pedited->cg.autodm        = !d_m->getAutoInconsistent();
        pedited->cg.d_y           = d_y->getEditedState ();
        pedited->cg.autody        = !d_y->getAutoInconsistent();
        pedited->cg.pwr           = pwr->getEditedState ();
        pedited->cg.enabled       = !get_inconsistent();
        pedited->cg.colorspace = colorspace->get_active_row_number() != 6;
        pedited->cg.rolloff       = !rolloff->get_inconsistent();
    }

}

void Compressgamut::updategamutGUI()
{
    // Update default slider value GUI according to colorspace
    // Only for Working profile = Rec2020
/* Calculating the 'Threshold' values ​​is anything but straightforward. Of course, one might say, "Just match them with the colors in ColorChecker24." But there are many unknown parameters:
* What is the actual illuminant and its temperature and green settings?
* The image colors may have been compressed using "Maximum limits," but there's no guarantee they're within the new, reduced gamut.
* Differences in white points, such as the atypical white point of DCI-P3, complicate matters.

* Furthermore, what about the necessary color adaptation, both at this stage of compression and at the stage implied by the output temperature (see CIECAM)?

Therefore, these parameters correspond roughly to the ratio of the distances between the white point of the Working Profile, 
* The white point of the "Target compression gamut" corrected by chromatic adaptation and the Yellow, Magenta, and Cyan values ​​of the primary color triangle.
* These calculated values ​​are for 'normal' cases where the illuminants are not exhaustive, and the colors to be recovered are within reasonable limits.
* In other cases, the 3 Threshold sliders need to be adjusted (often lowered). Pay attention to the artifacts.
*/



   
        if (colorspace->get_active_row_number() == 2){//Adobe D65
            th_c->setLimits(0., 1., 0.001, 0.79);
            th_m->setLimits(0., 1., 0.001, 0.82);
            th_y->setLimits(0., 1., 0.001, 0.85);
 
        } else if (colorspace->get_active_row_number() == 3){//srgb D65
            th_c->setLimits(0., 1., 0.001, 0.51);
            th_m->setLimits(0., 1., 0.001, 0.82);
            th_y->setLimits(0., 1., 0.001, 0.82);
        } else if (colorspace->get_active_row_number() == 4){//dci-p3 - WP 0.314 0.351 
            th_c->setLimits(0., 1., 0.001, 0.68);
            th_m->setLimits(0., 1., 0.001, 0.89);
            th_y->setLimits(0., 1., 0.001, 0.9);
        } else if (colorspace->get_active_row_number() == 6){//beta RGB D50
            th_c->setLimits(0., 1., 0.001, 0.90);
            th_m->setLimits(0., 1., 0.001, 0.95);
            th_y->setLimits(0., 1., 0.001, 0.95);
        } else {
            th_c->setLimits(0., 1., 0.001, 0.815);
            th_m->setLimits(0., 1., 0.001, 0.803);
            th_y->setLimits(0., 1., 0.001, 0.880);
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

void Compressgamut::adjusterAutoToggled(Adjuster* a, bool newval)
{
    if (multiImage) {
       if (d_c->getAutoInconsistent()) {
            d_c->setAutoInconsistent (false);
            d_c->setAutoValue (false);
        } else if (lastAutodc) {
            d_c->setAutoInconsistent (true);
        }

        lastAutodc = d_c->getAutoValue();
    }

    if (multiImage) {
       if (d_m->getAutoInconsistent()) {
            d_m->setAutoInconsistent (false);
            d_m->setAutoValue (false);
        } else if (lastAutodm) {
            d_m->setAutoInconsistent (true);
        }

        lastAutodm = d_m->getAutoValue(); 
    }

    if (multiImage) {
       if (d_y->getAutoInconsistent()) {
            d_y->setAutoInconsistent (false);
            d_y->setAutoValue (false);
        } else if (lastAutody) {
            d_y->setAutoInconsistent (true);
        }

        lastAutody = d_y->getAutoValue(); 
    }
    if (listener && (multiImage || getEnabled()) ) {

        if (a == d_c) {
            if (d_c->getAutoInconsistent()) {
                listener->panelChanged (Evcgdcautoon, M ("GENERAL_UNCHANGED"));
            } else if (d_c->getAutoValue()) {
                listener->panelChanged (Evcgdcautoon , M ("GENERAL_ENABLED"));
            } else {
                listener->panelChanged (Evcgdcautoon, M ("GENERAL_DISABLED"));
            }
        }

        if (a == d_m) {
            if (d_m->getAutoInconsistent()) {
                listener->panelChanged (Evcgdmautoon, M ("GENERAL_UNCHANGED"));
            } else if (d_m->getAutoValue()) {
                listener->panelChanged (Evcgdmautoon , M ("GENERAL_ENABLED"));
            } else {
                listener->panelChanged (Evcgdmautoon, M ("GENERAL_DISABLED"));
            }
        }

        if (a == d_y) {
            if (d_y->getAutoInconsistent()) {
                listener->panelChanged (Evcgdyautoon, M ("GENERAL_UNCHANGED"));
            } else if (d_y->getAutoValue()) {
                listener->panelChanged (Evcgdyautoon , M ("GENERAL_ENABLED"));
            } else {
                listener->panelChanged (Evcgdyautoon, M ("GENERAL_DISABLED"));
            }
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
