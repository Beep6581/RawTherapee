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
#include <clarity.h>
#include <guiutils.h>
#include <safegtk.h>
#include <sstream>
#include <iomanip>
#include <math.h>



using namespace rtengine;
using namespace rtengine::procparams;



Clarity::Clarity () : Gtk::VBox(), FoldableToolPanel(this)
{

   Gtk::HSeparator *hsep66aa = Gtk::manage (new  Gtk::HSeparator());
   pack_start(*hsep66aa, Gtk::PACK_SHRINK, 2);
   hsep66aa->show ();
    Gtk::Label* clsh = Gtk::manage (new Gtk::Label ());
    clsh->set_alignment (0.0, 0.5);
    clsh->set_markup (Glib::ustring("<b>") + M("TP_CLARITY_SHARPEN") + "</b>");
	clsh->show ();
	pack_start (*clsh);
   enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
   enabled->set_active (true);
   pack_start(*enabled);
   enabled->show ();
	
    Claritypasses = Gtk::manage(new Adjuster (M("TP_CLARITY_PASSES"),1,4,1,1));
	Claritypasses->setAdjusterListener (this);
	if (Claritypasses->delay < 1000) Claritypasses->delay = 1000;
	Claritypasses->show();
    Claritystrength = Gtk::manage(new Adjuster (M("TP_CLARITY_STRENGTH"),0,100,1,40));
	Claritystrength->setAdjusterListener (this);
	if (Claritystrength->delay < 1000) Claritystrength->delay = 1000;
	Claritystrength->show();
	
	Claritythreechannels = Gtk::manage(new Gtk::CheckButton((M("TP_CLARITY_THREE"))));// L + a + b	
	Claritythreechannels->set_active (false);
	pack_start( *Claritypasses, Gtk::PACK_SHRINK, 4);//passes
	pack_start( *Claritystrength, Gtk::PACK_SHRINK, 4);//strength
	pack_start( *Claritythreechannels, Gtk::PACK_SHRINK, 4);//one or 3 channels Lab
	Claritythreechannels->show();
	Gtk::HSeparator *hsep67aa = Gtk::manage (new  Gtk::HSeparator());
	pack_start(*hsep67aa, Gtk::PACK_SHRINK, 2);
	hsep67aa->show ();
	Gtk::Label* clmc = Gtk::manage (new Gtk::Label ());
    clmc->set_alignment (0.0, 0.5);
    clmc->set_markup (Glib::ustring("<b>") + M("TP_CLARITY_MICRO") + "</b>");
	clmc->show ();
	pack_start (*clmc);
	enabledtwo = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabledtwo->set_active (true);
	pack_start(*enabledtwo);
	enabledtwo->show ();
	
	
	MLmicrostrength= Gtk::manage(new Adjuster (M("TP_MLMICRO_STRENGTH"),0,100,1,25));
	MLmicrostrength->setAdjusterListener (this);
	if (MLmicrostrength->delay < 1000) MLmicrostrength->delay = 1000;
	MLmicrostrength->show();
	uniformity= Gtk::manage(new Adjuster (M("TP_MLMICRO_UNIFORMITY"),0,100,10,50));
	
	uniformity->setAdjusterListener (this);
	if (uniformity->delay < 1000) uniformity->delay = 1000;
	uniformity->show();
	MLmicromatrix = Gtk::manage (new Gtk::CheckButton (M("TP_CLARITY_MATRIX")));
	MLmicromatrix->set_active (true);
	pack_start(*MLmicromatrix);
	MLmicromatrix->show ();

	pack_start( *MLmicrostrength, Gtk::PACK_SHRINK, 4);//microcontraste strength	
	pack_start( *uniformity, Gtk::PACK_SHRINK, 4);//uniformity
 
    enaconn 	= enabled->signal_toggled().connect( sigc::mem_fun(*this, &Clarity::enabled_toggled) );
    chanthreeconn 	= Claritythreechannels->signal_toggled().connect( sigc::mem_fun(*this, &Clarity::chanthree_toggled) );
    enatwoconn 	= enabledtwo->signal_toggled().connect( sigc::mem_fun(*this, &Clarity::enabledtwo_toggled) );
    matrixconn 	= MLmicromatrix->signal_toggled().connect( sigc::mem_fun(*this, &Clarity::MLmicromatrix_toggled) );

}
Clarity::~Clarity () {

}

void Clarity::read(const ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();

	if(pedited ){
		Claritypasses->setEditedState 		(pedited->clarity.clpasses ? Edited : UnEdited);
		Claritystrength->setEditedState 		(pedited->clarity.clstrength ? Edited : UnEdited);
		MLmicrostrength->setEditedState 	(pedited->clarity.mlstrength ? Edited : UnEdited);
		uniformity->setEditedState 	(pedited->clarity.uniformity ? Edited : UnEdited);
		enabled->set_inconsistent 		(!pedited->clarity.enabled);
		Claritythreechannels->set_inconsistent 	(!pedited->clarity.clthreechannels);
		enabledtwo->set_inconsistent 		(!pedited->clarity.enabledtwo);
		MLmicromatrix->set_inconsistent 		(!pedited->clarity.MLmicromatrix);

	}
    enaconn.block (true);
    enabled->set_active (pp->clarity.enabled);
    enaconn.block (false);   
    lastEnabled = pp->clarity.enabled;
	
    enatwoconn.block (true);
    enabledtwo->set_active (pp->clarity.enabledtwo);
    enatwoconn.block (false);   
    lastEnabledtwo = pp->clarity.enabledtwo;
	
    chanthreeconn.block (true);
    Claritythreechannels->set_active (pp->clarity.clthreechannels);
    chanthreeconn.block (false);   
    lastchanthree = pp->clarity.clthreechannels;

    matrixconn.block (true);
    MLmicromatrix->set_active (pp->clarity.MLmicromatrix);
    matrixconn.block (false);   
    lastmatrix = pp->clarity.MLmicromatrix;
	
    Claritypasses->setValue        	    	(pp->clarity.clpasses);
    Claritystrength->setValue       	 (pp->clarity.clstrength);
    MLmicrostrength->setValue     		(pp->clarity.mlstrength);
    uniformity->setValue       (pp->clarity.uniformity);

	enableListener ();
}

void Clarity::write( ProcParams* pp, ParamsEdited* pedited)
{
	pp->clarity.clpasses          = (int)Claritypasses->getValue();
    pp->clarity.enabled        	  = enabled->get_active ();
    pp->clarity.clstrength        = Claritystrength->getValue ();
    pp->clarity.mlstrength        = MLmicrostrength->getValue ();
    pp->clarity.uniformity        = uniformity->getValue ();
    pp->clarity.clthreechannels   = Claritythreechannels->get_active ();
    pp->clarity.enabledtwo        = enabledtwo->get_active ();
    pp->clarity.MLmicromatrix    = MLmicromatrix->get_active ();


	if (pedited) {
		pedited->clarity.clpasses			= Claritypasses->getEditedState ();
		pedited->clarity.enabled			= !enabled->get_inconsistent();
		pedited->clarity.clthreechannels		= !Claritythreechannels->get_inconsistent();
		pedited->clarity.clstrength			= Claritystrength->getEditedState ();
		pedited->clarity.mlstrength			= MLmicrostrength->getEditedState ();
		pedited->clarity.uniformity			= uniformity->getEditedState ();
		pedited->clarity.enabledtwo			= !enabledtwo->get_inconsistent();
		pedited->clarity.MLmicromatrix	    = !MLmicromatrix->get_inconsistent();
		
	}

}

void Clarity::enabled_toggled () {

    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            enaconn.block (true);
            enabled->set_active (false);
            enaconn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);

        lastEnabled = enabled->get_active ();
    }

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvClarityEnabled, M("GENERAL_ENABLED"));
        //    listener->panelChanged (EvMLunifor, M("GENERAL_ENABLED"));
			
        else
            listener->panelChanged (EvClarityEnabled, M("GENERAL_DISABLED"));
    }
}
void Clarity::enabledtwo_toggled () {

    if (batchMode) {
        if (enabledtwo->get_inconsistent()) {
            enabledtwo->set_inconsistent (false);
            enatwoconn.block (true);
            enabledtwo->set_active (false);
            enatwoconn.block (false);
        }
        else if (lastEnabledtwo)
            enabledtwo->set_inconsistent (true);

        lastEnabledtwo = enabledtwo->get_active ();
    }

    if (listener) {
        if (enabledtwo->get_active ())
            listener->panelChanged (EvClarityEnabledtwo, M("GENERAL_ENABLED"));
			
        else
            listener->panelChanged (EvClarityEnabledtwo, M("GENERAL_DISABLED"));
			
    }
}

void Clarity::chanthree_toggled () {

    if (batchMode) {
        if (Claritythreechannels->get_inconsistent()) {
            Claritythreechannels->set_inconsistent (false);
            chanthreeconn.block (true);
            Claritythreechannels->set_active (false);
            chanthreeconn.block (false);
        }
        else if (lastchanthree)
            Claritythreechannels->set_inconsistent (true);

        lastchanthree = Claritythreechannels->get_active ();
    }

    if (listener && enabled->get_active ()) {
        if (Claritythreechannels->get_active ())
            listener->panelChanged (EvClaritythreechannels, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvClaritythreechannels, M("GENERAL_DISABLED"));
    }
}

void Clarity::MLmicromatrix_toggled () {
    if (batchMode) {
        if (MLmicromatrix->get_inconsistent()) {
            MLmicromatrix->set_inconsistent (false);
            matrixconn.block (true);
            MLmicromatrix->set_active (false);
            matrixconn.block (false);
        }
        else if (lastmatrix)
            MLmicromatrix->set_inconsistent (true);

        lastmatrix = MLmicromatrix->get_active ();
    }

    if (listener && enabledtwo->get_active ()) {
        if (MLmicromatrix->get_active ())
            listener->panelChanged (EvClaritymatrix, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvClaritymatrix, M("GENERAL_DISABLED"));
    }

}


void Clarity::adjusterChanged (Adjuster* a, double newval)
{
	if (listener && enabled->get_active()) {
		Glib::ustring value = a->getTextValue();
		{ 
		
		if (a == Claritypasses ) 
			listener->panelChanged (EvClaritypasses,  value );
		else if (a == Claritystrength)
			listener->panelChanged (EvClaritystrength,  value );
			}
			
	}
	if (listener && enabledtwo->get_active()) {
		Glib::ustring value = a->getTextValue();
		{ 
			if (a == MLmicrostrength)	
			listener->panelChanged (EvMLmicrostrength,  value );
			else if (a == uniformity)	
			listener->panelChanged (EvMLuniformity,  value );
		
			}
			
	}



}

void Clarity::setBatchMode(bool batchMode)
{
    Claritypasses->showEditedCB ();
    Claritystrength->showEditedCB ();
    MLmicrostrength->showEditedCB ();
    uniformity->showEditedCB ();
	
}

void Clarity::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited)
{
    Claritypasses->setDefault (defParams->clarity.clpasses);
    Claritystrength->setDefault (defParams->clarity.clstrength);
	MLmicrostrength->setDefault (defParams->clarity.mlstrength);
	uniformity->setDefault (defParams->clarity.uniformity);	

	if (pedited) {
       Claritypasses->setDefaultEditedState 		(pedited->clarity.clpasses ? Edited : UnEdited);
       Claritystrength->setDefaultEditedState 		(pedited->clarity.clstrength ? Edited : UnEdited);
       MLmicrostrength->setDefaultEditedState 		(pedited->clarity.mlstrength ? Edited : UnEdited);
       uniformity->setDefaultEditedState 		(pedited->clarity.uniformity ? Edited : UnEdited);
		
	} else {
       Claritypasses->setDefaultEditedState 		(Irrelevant);
       Claritystrength->setDefaultEditedState 		(Irrelevant);
       MLmicrostrength->setDefaultEditedState 		(Irrelevant);
       uniformity->setDefaultEditedState 		(Irrelevant);
		
	}
}

void Clarity::setAdjusterBehavior (bool strengthadd, bool mlstrentghadd, bool passadd, bool uniformityadd ) {
	Claritystrength->setAddMode(strengthadd);
	MLmicrostrength->setAddMode(mlstrentghadd);
	Claritypasses->setAddMode(passadd);
	uniformity->setAddMode(uniformityadd);
}

void Clarity::trimValues (ProcParams* pp) {
	Claritystrength->trimValue(pp->clarity.clstrength);
	MLmicrostrength->trimValue(pp->clarity.mlstrength);
	Claritypasses->trimValue(pp->clarity.clpasses);
	uniformity->trimValue(pp->clarity.uniformity);

}
