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
#include "rawprocess.h"
#include "options.h"
#include "guiutils.h"
using namespace rtengine;
using namespace rtengine::procparams;

RawProcess::RawProcess () : FoldableToolPanel(this)
{
   Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
   hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_DMETHOD") +": ")),Gtk::PACK_SHRINK, 4);
   dmethod = Gtk::manage (new MyComboBoxText ());
   for( size_t i=0; i<procparams::RAWParams::numMethods;i++)
	   dmethod->append_text(procparams::RAWParams::methodstring[i]);

   dmethod->set_active(0);
   hb1->set_tooltip_markup (M("TP_RAW_DMETHOD_TOOLTIP"));

   hb1->pack_end (*dmethod, Gtk::PACK_EXPAND_WIDGET, 4);
   pack_start( *hb1, Gtk::PACK_SHRINK, 4);

   dcbOptions = Gtk::manage (new Gtk::VBox ());
   dcbOptions->set_border_width(4);

   dcbIterations = Gtk::manage (new Adjuster (M("TP_RAW_DCBITERATIONS"),0,5,1,2));
   dcbIterations->setAdjusterListener (this);
   if (dcbIterations->delay < 1000) dcbIterations->delay = 1000;
   dcbIterations->show();
   dcbEnhance = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_DCBENHANCE")));
   dcbOptions->pack_start(*dcbIterations);
   dcbOptions->pack_start(*dcbEnhance);
   pack_start( *dcbOptions, Gtk::PACK_SHRINK, 4);

   lmmseOptions = Gtk::manage (new Gtk::VBox ());
   lmmseOptions->set_border_width(4);

   lmmseIterations = Gtk::manage (new Adjuster (M("TP_RAW_LMMSEITERATIONS"),0,6,1,2));
   lmmseIterations->setAdjusterListener (this);
   lmmseIterations->set_tooltip_markup (M("TP_RAW_LMMSE_TOOLTIP"));

   if (lmmseIterations->delay < 1000) lmmseIterations->delay = 1000;
   lmmseIterations->show();
   lmmseOptions->pack_start(*lmmseIterations);
   pack_start( *lmmseOptions, Gtk::PACK_SHRINK, 4);

   pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
   ccSteps = Gtk::manage (new Adjuster (M("TP_RAW_FALSECOLOR"),0,5,1,0 ));
   ccSteps->setAdjusterListener (this);
   if (ccSteps->delay < 1000) ccSteps->delay = 1000;
   ccSteps->show();
   pack_start( *ccSteps, Gtk::PACK_SHRINK, 4);

   //pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
   //allOptions = Gtk::manage (new Gtk::VBox ());
   //allOptions->set_border_width(2);
   //allEnhance = Gtk::manage (new Gtk::CheckButton(M("TP_RAW_ALLENHANCE")));
   //allOptions->pack_start(*allEnhance);
   //pack_start( *allOptions, Gtk::PACK_SHRINK, 4);

   methodconn = dmethod->signal_changed().connect( sigc::mem_fun(*this, &RawProcess::methodChanged) );
   dcbEnhconn = dcbEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &RawProcess::dcbEnhanceChanged), true);
   //allEnhconn = allEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &RawProcess::allEnhanceChanged), true);
}


void RawProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
   disableListener ();
   methodconn.block (true);
   dcbEnhconn.block (true);
   //allEnhconn.block (true);

   dmethod->set_active(procparams::RAWParams::numMethods);
   for( size_t i=0; i< procparams::RAWParams::numMethods;i++)
       if( pp->raw.dmethod == procparams::RAWParams::methodstring[i]){
           dmethod->set_active(i);
	   oldSelection = i;
           break;
       }

   if(pedited ){
       ccSteps->setEditedState (pedited->raw.ccSteps ? Edited : UnEdited);
       dcbIterations->setEditedState ( pedited->raw.dcbIterations ? Edited : UnEdited);
       dcbEnhance->set_inconsistent(!pedited->raw.dcbEnhance);
       //allEnhance->set_inconsistent(!pedited->raw.allEnhance);
       lmmseIterations->setEditedState ( pedited->raw.lmmseIterations ? Edited : UnEdited);

       if( !pedited->raw.dmethod )
           dmethod->set_active(procparams::RAWParams::numMethods); // No name
   }

   //allEnhance->set_active(pp->raw.all_enhance);

   dcbIterations->setValue (pp->raw.dcb_iterations);
   dcbEnhance->set_active(pp->raw.dcb_enhance);
   ccSteps->setValue (pp->raw.ccSteps);
   if (pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::dcb] ||
       dmethod->get_active_row_number() == procparams::RAWParams::numMethods)
       dcbOptions->show();
   else
       dcbOptions->hide();

   lmmseIterations->setValue (pp->raw.lmmse_iterations);
   if (pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::lmmse] ||
       dmethod->get_active_row_number() == procparams::RAWParams::numMethods)
       lmmseOptions->show();
   else
       lmmseOptions->hide();

   // Flase color suppression is applied to all demozaicing method, so don't hide anything
   /*if (pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::eahd] ||
       pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::hphd] ||
       pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::vng4])
       ccSteps->show();
   else
       ccSteps->hide();*/

   lastDCBen = pp->raw.dcb_enhance;
   //lastALLen = pp->raw.all_enhance;

   methodconn.block (false);
   dcbEnhconn.block (false);
   //allEnhconn.block (false);
   
   enableListener ();
}

void RawProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.ccSteps = ccSteps->getIntValue();
	pp->raw.dcb_iterations = dcbIterations->getIntValue();
	pp->raw.dcb_enhance = dcbEnhance->get_active();
	//pp->raw.all_enhance = allEnhance->get_active();
	pp->raw.lmmse_iterations = lmmseIterations->getIntValue();

	int currentRow = dmethod->get_active_row_number();
	if( currentRow>=0 && currentRow < procparams::RAWParams::numMethods)
		pp->raw.dmethod = procparams::RAWParams::methodstring[currentRow];

	if (pedited) {
		pedited->raw.ccSteps = ccSteps->getEditedState ();
		pedited->raw.dmethod = dmethod->get_active_row_number() != procparams::RAWParams::numMethods;
		pedited->raw.dcbIterations = dcbIterations->getEditedState ();
		pedited->raw.dcbEnhance = !dcbEnhance->get_inconsistent();
		//pedited->raw.allEnhance = !allEnhance->get_inconsistent();
		pedited->raw.lmmseIterations = lmmseIterations->getEditedState ();
		
	}
}

void RawProcess::setBatchMode(bool batchMode)
{
   dmethod->append_text (M("GENERAL_UNCHANGED"));
   dmethod->set_active(procparams::RAWParams::numMethods); // No name
   dcbOptions->hide();
   lmmseOptions->hide();
   ToolPanel::setBatchMode (batchMode);
   ccSteps->showEditedCB ();
   dcbIterations->showEditedCB ();
   lmmseIterations->showEditedCB ();
   
}

void RawProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	dcbIterations->setDefault( defParams->raw.dcb_iterations);
	lmmseIterations->setDefault( defParams->raw.lmmse_iterations);
	ccSteps->setDefault (defParams->raw.ccSteps);
	if (pedited) {
		dcbIterations->setDefaultEditedState( pedited->raw.dcbIterations ? Edited : UnEdited);
		lmmseIterations->setDefaultEditedState( pedited->raw.lmmseIterations ? Edited : UnEdited);
		ccSteps->setDefaultEditedState(pedited->raw.ccSteps ? Edited : UnEdited);
	}else{
		dcbIterations->setDefaultEditedState( Irrelevant );
		lmmseIterations->setDefaultEditedState( Irrelevant );
		ccSteps->setDefaultEditedState(Irrelevant );
	}
}

void RawProcess::adjusterChanged (Adjuster* a, double newval)
{
    if (listener) {
    	if (a == dcbIterations)
    		listener->panelChanged (EvDemosaicDCBIter, a->getTextValue() );
    	else if (a == ccSteps)
    		listener->panelChanged (EvDemosaicFalseColorIter, a->getTextValue() );
    	else if (a == lmmseIterations)
    		listener->panelChanged (EvDemosaicLMMSEIter, a->getTextValue() );

			}
}

void RawProcess::methodChanged ()
{
	int  curSelection = dmethod->get_active_row_number();
	if ( curSelection == procparams::RAWParams::dcb){
		dcbOptions->show();
	}else{
		dcbOptions->hide();
	}
	if ( curSelection == procparams::RAWParams::lmmse){
		lmmseOptions->show();
	}else{
		lmmseOptions->hide();
	}
	
	Glib::ustring methodName="";
	bool ppreq = false;
	if( curSelection>=0 && curSelection < procparams::RAWParams::numMethods) {
	    methodName = procparams::RAWParams::methodstring[curSelection];
	    if (curSelection == procparams::RAWParams::mono || oldSelection == procparams::RAWParams::mono) {
		    ppreq = true;
	    }
	}
	oldSelection = curSelection;

    if (listener)
	    listener->panelChanged (ppreq ? EvDemosaicMethodPreProc : EvDemosaicMethod, methodName);
}

void RawProcess::dcbEnhanceChanged ()
{
    if (batchMode) {
        if (dcbEnhance->get_inconsistent()) {
        	dcbEnhance->set_inconsistent (false);
        	dcbEnhconn.block (true);
        	dcbEnhance->set_active (false);
        	dcbEnhconn.block (false);
        }
        else if (lastDCBen)
        	dcbEnhance->set_inconsistent (true);

        lastDCBen = dcbEnhance->get_active ();
    }
    if (listener)
        listener->panelChanged (EvDemosaicDCBEnhanced, dcbEnhance->get_active()?M("GENERAL_ENABLED"):M("GENERAL_DISABLED"));
}

/*void RawProcess::allEnhanceChanged ()
{
    if (batchMode) {
        if (allEnhance->get_inconsistent()) {
        	allEnhance->set_inconsistent (false);
        	allEnhconn.block (true);
        	allEnhance->set_active (false);
        	allEnhconn.block (false);
        }
        else if (lastALLen)
        	allEnhance->set_inconsistent (true);

        lastALLen = allEnhance->get_active ();
    }
    if (listener)
        listener->panelChanged (EvDemosaicALLEnhanced, allEnhance->get_active()?M("GENERAL_ENABLED"):M("GENERAL_DISABLED"));
}*/
