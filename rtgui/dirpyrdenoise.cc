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
#include "dirpyrdenoise.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

DirPyrDenoise::DirPyrDenoise () : FoldableToolPanel(this), lastenhance(false) {
	
	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_tooltip_text (M("TP_DIRPYRDENOISE_ENABLED_TOOLTIP"));
	
	enabled->set_active (false);
	enabled->show ();
	pack_start (*enabled);
	
	Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
	hsep1->show ();
	pack_start (*hsep1);
	
	enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::enabledChanged) );
	

	luma  = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_LUMA"), 0, 100, 0.01, 0));
	Ldetail  = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_LDETAIL"), 0, 100, 0.01, 50));
	chroma    = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_CHROMA"), 0, 100, 0.01, 15));
	redchro    = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_RED"), -100, 100, 0.1, 0));
	bluechro    = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_BLUE"), -100, 100, 0.1, 0));
	
	gamma	= Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_GAMMA"), 1.0, 3.0, 0.01, 1.7));
	gamma->set_tooltip_text (M("TP_DIRPYRDENOISE_GAMMA_TOOLTIP"));

	passes	= Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_PASSE"), 1.0, 3.0, 1., 1.));
	passes->set_tooltip_text (M("TP_DIRPYRDENOISE_PASSES_TOOLTIP"));

	
	Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
	hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_DIRPYRDENOISE_METHOD") +": ")),Gtk::PACK_SHRINK, 4);
	hb1->set_tooltip_markup (M("TP_DIRPYRDENOISE_METHOD_TOOLTIP"));
   
	dmethod = Gtk::manage (new MyComboBoxText ());
	dmethod->append_text (M("TP_DIRPYRDENOISE_RGB"));
	dmethod->append_text (M("TP_DIRPYRDENOISE_LAB"));
	dmethod->set_active(0);
	hb1->pack_end (*dmethod, Gtk::PACK_EXPAND_WIDGET, 4);
	pack_start( *hb1, Gtk::PACK_SHRINK, 4);
	
	
//	perform = Gtk::manage (new Gtk::CheckButton (M("TP_DIRPYRDENOISE_PERF")));
//	perform->set_tooltip_text (M("TP_DIRPYRDENOISE_PERF_TOOLTIP"));

//	perfconn = perform->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::perform_toggled) );
	dmethodconn = dmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::dmethodChanged) );
	
	luma->setAdjusterListener (this);
	Ldetail->setAdjusterListener (this);
	chroma->setAdjusterListener (this); 
	redchro->setAdjusterListener (this); 
	bluechro->setAdjusterListener (this); 
	
    gamma->setAdjusterListener (this); 
    passes->setAdjusterListener (this); 
	
    luma->show();
    Ldetail->show();
    chroma->show();
    redchro->show();
    bluechro->show();
//	perform->show();
	gamma->show();
//	perform->set_active (true);
	passes->show();

	enhance = Gtk::manage (new Gtk::CheckButton (M("TP_DIRPYRDENOISE_ENH")));
	enhance->set_active (false);
	enhance->set_tooltip_text (M("TP_DIRPYRDENOISE_ENH_TOOLTIP"));

	median = Gtk::manage (new Gtk::CheckButton (M("TP_DIRPYRDENOISE_MED")));
	median->set_active (true);
	median->set_tooltip_text (M("TP_DIRPYRDENOISE_MED_TOOLTIP"));

	
	Gtk::HSeparator *hsep2 = Gtk::manage (new  Gtk::HSeparator());
	hsep2->show ();
	
	methodmed = Gtk::manage (new MyComboBoxText ());
	methodmed->append_text (M("TP_DIRPYRDENOISE_NONE"));
	methodmed->append_text (M("TP_DIRPYRDENOISE_LM"));
	methodmed->append_text (M("TP_DIRPYRDENOISE_LABM"));
	methodmed->append_text (M("TP_DIRPYRDENOISE_RGBM"));
	methodmed->set_active (0);
	methodmed->set_tooltip_text (M("TP_DIRPYRDENOISE_METM_TOOLTIP"));
	methodmedconn = methodmed->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::methodmedChanged) );
	
	rgbmethod = Gtk::manage (new MyComboBoxText ());
	rgbmethod->append_text (M("TP_DIRPYRDENOISE_SOFT"));
	rgbmethod->append_text (M("TP_DIRPYRDENOISE_33"));
	rgbmethod->append_text (M("TP_DIRPYRDENOISE_55SOFT"));
	rgbmethod->set_active (0);
	rgbmethod->set_tooltip_text (M("TP_DIRPYRDENOISE_MET_TOOLTIP"));
	rgbmethodconn = rgbmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::rgbmethodChanged) );
	
	
	medmethod = Gtk::manage (new MyComboBoxText ());
	medmethod->append_text (M("TP_DIRPYRDENOISE_SOFT"));
	medmethod->append_text (M("TP_DIRPYRDENOISE_33"));
	medmethod->append_text (M("TP_DIRPYRDENOISE_55SOFT"));
	medmethod->append_text (M("TP_DIRPYRDENOISE_55"));
	medmethod->append_text (M("TP_DIRPYRDENOISE_77"));
	medmethod->set_active (0);
	medmethod->set_tooltip_text (M("TP_DIRPYRDENOISE_MET_TOOLTIP"));
	medmethodconn = medmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::medmethodChanged) );

	ctboxm = Gtk::manage (new Gtk::HBox ());
	Gtk::Label* labmm = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_METHODMED")));
	ctboxm->pack_start (*labmm, Gtk::PACK_SHRINK, 4);
	
	ctbox = Gtk::manage (new Gtk::HBox ());
	Gtk::Label* labm = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_MEDMETHOD")));
	ctbox->pack_start (*labm, Gtk::PACK_SHRINK, 4);

	ctboxrgb = Gtk::manage (new Gtk::HBox ());
	Gtk::Label* labrgb = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_MEDMETHOD")));
	ctboxrgb->pack_start (*labrgb, Gtk::PACK_SHRINK, 4);
	
	pack_start (*luma);
	pack_start (*Ldetail);
	pack_start (*chroma);
	pack_start (*redchro);
	pack_start (*bluechro);
	
	pack_start (*gamma);
	pack_start (*enhance);
	
	pack_start (*hsep2);
	
	ctboxm->pack_start (*methodmed);
	ctbox->pack_start (*medmethod);
	ctboxrgb->pack_start (*rgbmethod);
	pack_start (*ctboxm);
	pack_start (*ctbox);
	pack_start (*ctboxrgb);
	pack_start (*passes);
	
//	pack_start (*median);
	
//	pack_start (*perform);
	enhanConn = enhance->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::enhanceChanged) );
	medianConn = median->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::medianChanged) );
	ctboxrgb->hide();
}

void DirPyrDenoise::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
    dmethodconn.block(true);
    enaConn.block (true);
	medmethodconn.block(true);
	rgbmethodconn.block(true);
	methodmedconn.block(true);

    dmethod->set_active (0);
    if (pp->dirpyrDenoise.dmethod=="RGB")
        dmethod->set_active (0);
    else if (pp->dirpyrDenoise.dmethod=="Lab")
        dmethod->set_active (1);

    methodmed->set_active (0);
    if (pp->dirpyrDenoise.methodmed=="none")
        methodmed->set_active (0);
    else if (pp->dirpyrDenoise.methodmed=="Lonly")
        methodmed->set_active (1);
    else if (pp->dirpyrDenoise.methodmed=="Lab")
        methodmed->set_active (2);
    else if (pp->dirpyrDenoise.methodmed=="RGB")
        methodmed->set_active (3);
	methodmedChanged();
		
    medmethod->set_active (0);
//    if (pp->dirpyrDenoise.medmethod=="none")
//        medmethod->set_active (0);
    if (pp->dirpyrDenoise.medmethod=="soft")
        medmethod->set_active (0);
    else if (pp->dirpyrDenoise.medmethod=="33")
        medmethod->set_active (1);
    else if (pp->dirpyrDenoise.medmethod=="55soft")
        medmethod->set_active (2);
    else if (pp->dirpyrDenoise.medmethod=="55")
        medmethod->set_active (3);
    else if (pp->dirpyrDenoise.medmethod=="77")
        medmethod->set_active (4);
	medmethodChanged();
		
    rgbmethod->set_active (0);
//    if (pp->dirpyrDenoise.medmethod=="none")
//        medmethod->set_active (0);
    if (pp->dirpyrDenoise.rgbmethod=="soft")
        rgbmethod->set_active (0);
    else if (pp->dirpyrDenoise.rgbmethod=="33")
        rgbmethod->set_active (1);
    else if (pp->dirpyrDenoise.rgbmethod=="55soft")
        rgbmethod->set_active (2);
	rgbmethodChanged();
		
		
    if (pedited) {
        if (!pedited->dirpyrDenoise.dmethod)
            dmethod->set_active (2);
        if (!pedited->dirpyrDenoise.rgbmethod)
            rgbmethod->set_active (2);
         if (!pedited->dirpyrDenoise.medmethod)
            medmethod->set_active (4);
		if (!pedited->dirpyrDenoise.methodmed)
            methodmed->set_active (3);
			
        luma->setEditedState      (pedited->dirpyrDenoise.luma ? Edited : UnEdited);
        Ldetail->setEditedState   (pedited->dirpyrDenoise.Ldetail ? Edited : UnEdited);
        chroma->setEditedState    (pedited->dirpyrDenoise.chroma ? Edited : UnEdited);
        redchro->setEditedState    (pedited->dirpyrDenoise.redchro ? Edited : UnEdited);
        bluechro->setEditedState    (pedited->dirpyrDenoise.bluechro ? Edited : UnEdited);

        gamma->setEditedState     (pedited->dirpyrDenoise.gamma ? Edited : UnEdited);
        passes->setEditedState     (pedited->dirpyrDenoise.passes ? Edited : UnEdited);
        enabled->set_inconsistent (!pedited->dirpyrDenoise.enabled);
        enhance->set_inconsistent (!pedited->dirpyrDenoise.enhance);
        median->set_inconsistent (!pedited->dirpyrDenoise.median);
  //      perform->set_inconsistent (!pedited->dirpyrDenoise.perform);
    }
//	perfconn.block (true);
    enabled->set_active (pp->dirpyrDenoise.enabled);
    enhance->set_active (pp->dirpyrDenoise.enhance);
 //   perform->set_active (pp->dirpyrDenoise.perform);
    median->set_active (pp->dirpyrDenoise.median);

 //   perfconn.block (false);
    lastEnabled = pp->dirpyrDenoise.enabled;
    lastmedian = pp->dirpyrDenoise.median;
    lastenhance = pp->dirpyrDenoise.enhance;
//	lastperform = pp->dirpyrDenoise.perform;	
    luma->setValue    (pp->dirpyrDenoise.luma);
    Ldetail->setValue (pp->dirpyrDenoise.Ldetail);
    chroma->setValue  (pp->dirpyrDenoise.chroma);
    redchro->setValue  (pp->dirpyrDenoise.redchro);
    bluechro->setValue  (pp->dirpyrDenoise.bluechro);

    gamma->setValue   (pp->dirpyrDenoise.gamma);
    passes->setValue   (pp->dirpyrDenoise.passes);

    enaConn.block (false);
    dmethodconn.block(false);
    medmethodconn.block(false);
    rgbmethodconn.block(false);
    methodmedconn.block(false);
    enableListener ();

}

void DirPyrDenoise::write (ProcParams* pp, ParamsEdited* pedited) {
	
	pp->dirpyrDenoise.luma      = luma->getValue ();
	pp->dirpyrDenoise.Ldetail   = Ldetail->getValue ();
	pp->dirpyrDenoise.chroma	= chroma->getValue ();
	pp->dirpyrDenoise.redchro	= redchro->getValue ();
	pp->dirpyrDenoise.bluechro	= bluechro->getValue ();
	pp->dirpyrDenoise.gamma		= gamma->getValue ();
	pp->dirpyrDenoise.passes	= passes->getValue ();
	pp->dirpyrDenoise.enabled   = enabled->get_active();
	pp->dirpyrDenoise.enhance   = enhance->get_active();
//	pp->dirpyrDenoise.perform   = perform->get_active();
	pp->dirpyrDenoise.median   = median->get_active();
	
    if (pedited) {
        pedited->dirpyrDenoise.dmethod  = dmethod->get_active_row_number() != 2;
        pedited->dirpyrDenoise.medmethod  = medmethod->get_active_row_number() != 4;
        pedited->dirpyrDenoise.rgbmethod  = rgbmethod->get_active_row_number() != 2;
        pedited->dirpyrDenoise.methodmed  = methodmed->get_active_row_number() != 3;
        pedited->dirpyrDenoise.luma     = luma->getEditedState ();
        pedited->dirpyrDenoise.Ldetail  = Ldetail->getEditedState ();
        pedited->dirpyrDenoise.chroma   = chroma->getEditedState ();
        pedited->dirpyrDenoise.redchro  = redchro->getEditedState ();
        pedited->dirpyrDenoise.bluechro = bluechro->getEditedState ();
        pedited->dirpyrDenoise.gamma    = gamma->getEditedState ();
        pedited->dirpyrDenoise.passes    = passes->getEditedState ();
        pedited->dirpyrDenoise.enabled  = !enabled->get_inconsistent();
        pedited->dirpyrDenoise.enhance  = !enhance->get_inconsistent();
        pedited->dirpyrDenoise.median  = !median->get_inconsistent();
    //    pedited->dirpyrDenoise.perform  = !perform->get_inconsistent();
    }
    if (dmethod->get_active_row_number()==0)
        pp->dirpyrDenoise.dmethod = "RGB";
    else if (dmethod->get_active_row_number()==1)
        pp->dirpyrDenoise.dmethod = "Lab";

		if (methodmed->get_active_row_number()==0)
        pp->dirpyrDenoise.methodmed = "none";	
    else if (methodmed->get_active_row_number()==1)
        pp->dirpyrDenoise.methodmed = "Lonly";
    else if (methodmed->get_active_row_number()==2)
        pp->dirpyrDenoise.methodmed = "Lab";
    else if (methodmed->get_active_row_number()==3)
        pp->dirpyrDenoise.methodmed = "RGB";


		
 //   if (medmethod->get_active_row_number()==0)
 //       pp->dirpyrDenoise.medmethod = "none";	
    if (medmethod->get_active_row_number()==0)
        pp->dirpyrDenoise.medmethod = "soft";
    else if (medmethod->get_active_row_number()==1)
        pp->dirpyrDenoise.medmethod = "33";
    else if (medmethod->get_active_row_number()==2)
        pp->dirpyrDenoise.medmethod = "55soft";
    else if (medmethod->get_active_row_number()==3)
        pp->dirpyrDenoise.medmethod = "55";
    else if (medmethod->get_active_row_number()==4)
        pp->dirpyrDenoise.medmethod = "77";

    if (rgbmethod->get_active_row_number()==0)
        pp->dirpyrDenoise.rgbmethod = "soft";
    else if (rgbmethod->get_active_row_number()==1)
        pp->dirpyrDenoise.rgbmethod = "33";
    else if (rgbmethod->get_active_row_number()==2)
        pp->dirpyrDenoise.rgbmethod = "55soft";
		
}

void DirPyrDenoise::dmethodChanged () {

	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvDPDNmet, dmethod->get_active_text ());
	}
}

void DirPyrDenoise::medmethodChanged () {
	
	
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvDPDNmedmet, medmethod->get_active_text ());
	}
}

void DirPyrDenoise::rgbmethodChanged () {
	ctboxrgb->hide();
	if(methodmed->get_active_row_number()==3) ctboxrgb->show();
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvDPDNrgbmet, rgbmethod->get_active_text ());
	}
}



void DirPyrDenoise::methodmedChanged () {
	if(methodmed->get_active_row_number()==3) {ctboxrgb->show();ctbox->hide();}
	else {ctboxrgb->hide();ctbox->show();}

	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvDPDNmetmed, methodmed->get_active_text ());
	}
}

void DirPyrDenoise::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {
	
	luma->setDefault    (defParams->dirpyrDenoise.luma);
	Ldetail->setDefault (defParams->dirpyrDenoise.Ldetail);
	chroma->setDefault  (defParams->dirpyrDenoise.chroma);
	redchro->setDefault  (defParams->dirpyrDenoise.redchro);
	bluechro->setDefault  (defParams->dirpyrDenoise.bluechro);
	gamma->setDefault   (defParams->dirpyrDenoise.gamma);
	passes->setDefault   (defParams->dirpyrDenoise.passes);

    if (pedited) {
    	luma->setDefaultEditedState     (pedited->dirpyrDenoise.luma ? Edited : UnEdited);
    	Ldetail->setDefaultEditedState	(pedited->dirpyrDenoise.Ldetail ? Edited : UnEdited);
    	chroma->setDefaultEditedState   (pedited->dirpyrDenoise.chroma ? Edited : UnEdited);
    	redchro->setDefaultEditedState   (pedited->dirpyrDenoise.redchro ? Edited : UnEdited);
    	bluechro->setDefaultEditedState   (pedited->dirpyrDenoise.bluechro ? Edited : UnEdited);
		gamma->setDefaultEditedState    (pedited->dirpyrDenoise.gamma ? Edited : UnEdited);
		passes->setDefaultEditedState    (pedited->dirpyrDenoise.passes ? Edited : UnEdited);
   }
    else {
    	luma->setDefaultEditedState     (Irrelevant);
    	Ldetail->setDefaultEditedState  (Irrelevant);
    	chroma->setDefaultEditedState   (Irrelevant);
    	redchro->setDefaultEditedState   (Irrelevant);
    	bluechro->setDefaultEditedState   (Irrelevant);
        gamma->setDefaultEditedState    (Irrelevant);
        passes->setDefaultEditedState    (Irrelevant);
    }
}

void DirPyrDenoise::adjusterChanged (Adjuster* a, double newval) {
	
	Glib::ustring costr;
    costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
        
    if (listener && enabled->get_active()) {
        if (a==Ldetail)
            listener->panelChanged (EvDPDNLdetail, costr);
        else if (a==luma)
            listener->panelChanged (EvDPDNLuma, costr);
        else if (a==chroma)
            listener->panelChanged (EvDPDNChroma, costr);
       else if (a==redchro)
            listener->panelChanged (EvDPDNredchro, costr);
        else if (a==bluechro)
            listener->panelChanged (EvDPDNbluechro, costr);		
        else if (a==gamma)
			listener->panelChanged (EvDPDNGamma, costr);
        else if (a==passes)
			listener->panelChanged (EvDPDNpasses, costr);
	}
}

void DirPyrDenoise::enabledChanged () {
	
    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            enaConn.block (true);
            enabled->set_active (false);
            enaConn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);
		
        lastEnabled = enabled->get_active ();
    }
	
    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvDPDNEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvDPDNEnabled, M("GENERAL_DISABLED"));
    }  
}

void DirPyrDenoise::enhanceChanged () {
	
    if (batchMode) {
        if (enhance->get_inconsistent()) {
            enhance->set_inconsistent (false);
            enhanConn.block (true);
            enhance->set_active (false);
            enhanConn.block (false);
        }
        else if (lastenhance)
            enhance->set_inconsistent (true);
		
        lastenhance = enhance->get_active ();
    }
	
    if (listener) {
        if (enhance->get_active ())
            listener->panelChanged (EvDPDNenhance, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvDPDNenhance, M("GENERAL_DISABLED"));
    }  
}

void DirPyrDenoise::medianChanged () {
	
    if (batchMode) {
        if (median->get_inconsistent()) {
            median->set_inconsistent (false);
            medianConn.block (true);
            median->set_active (false);
            medianConn.block (false);
        }
        else if (lastmedian)
            median->set_inconsistent (true);
		
        lastmedian = median->get_active ();
    }
	
    if (listener) {
        if (median->get_active ())
            listener->panelChanged (EvDPDNmedian, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvDPDNmedian, M("GENERAL_DISABLED"));
    }  
}

/*
void DirPyrDenoise::perform_toggled () {

    if (batchMode) {
        if (perform->get_inconsistent()) {
            perform->set_inconsistent (false);
            perfconn.block (true);
            perform->set_active (false);
            perfconn.block (false);
        }
        else if (lastperform)
            perform->set_inconsistent (true);

        lastperform = perform->get_active ();
    }

    if (listener) {
        if (perform->get_active ())
            listener->panelChanged (EvDPDNperform, M("GENERAL_ENABLED"));
        else            
            listener->panelChanged (EvDPDNperform, M("GENERAL_DISABLED"));
    }
}
*/

void DirPyrDenoise::setBatchMode (bool batchMode) {
	
    ToolPanel::setBatchMode (batchMode);
    luma->showEditedCB ();
    Ldetail->showEditedCB ();
    chroma->showEditedCB ();
    redchro->showEditedCB ();
    bluechro->showEditedCB ();
    gamma->showEditedCB ();
    passes->showEditedCB ();
    dmethod->append_text (M("GENERAL_UNCHANGED"));
    medmethod->append_text (M("GENERAL_UNCHANGED"));
    methodmed->append_text (M("GENERAL_UNCHANGED"));
    rgbmethod->append_text (M("GENERAL_UNCHANGED"));
	
}

void DirPyrDenoise::setAdjusterBehavior (bool lumaadd, bool lumdetadd, bool chromaadd, bool chromaredadd, bool chromablueadd, bool gammaadd, bool passesadd) {

	luma->setAddMode(lumaadd);
	Ldetail->setAddMode(lumdetadd);
	chroma->setAddMode(chromaadd);
	redchro->setAddMode(chromaredadd);
	bluechro->setAddMode(chromablueadd);
	gamma->setAddMode(gammaadd);
	passes->setAddMode(passesadd);
	
}

void DirPyrDenoise::trimValues (rtengine::procparams::ProcParams* pp) {

	luma->trimValue(pp->dirpyrDenoise.luma);
	Ldetail->trimValue(pp->dirpyrDenoise.Ldetail);
	chroma->trimValue(pp->dirpyrDenoise.chroma);
	redchro->trimValue(pp->dirpyrDenoise.redchro);
	bluechro->trimValue(pp->dirpyrDenoise.bluechro);
	gamma->trimValue(pp->dirpyrDenoise.gamma);
	passes->trimValue(pp->dirpyrDenoise.passes);
}
