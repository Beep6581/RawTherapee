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
#include <rawprocess.h>
#include <options.h>
#include <guiutils.h>
using namespace rtengine;
using namespace rtengine::procparams;

RawProcess::RawProcess () : Gtk::VBox(), FoldableToolPanel(this)
{
   Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
   hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("PREFERENCES_DMETHOD") +": ")));
   dmethod = Gtk::manage (new Gtk::ComboBoxText ());
   for( size_t i=0; i< procparams::RAWParams::numMethods;i++)
	   dmethod->append_text(procparams::RAWParams::methodstring[i]);

   dmethod->set_active(0);
   hb1->pack_end (*dmethod);
   pack_start( *hb1, Gtk::PACK_SHRINK, 4);

   dcbOptions = Gtk::manage (new Gtk::VBox ());
   dcbOptions->set_border_width(4);

   dcbIterations = Gtk::manage (new Adjuster (M("PREFERENCES_DCBITERATIONS"),0,5,1,2));
   dcbIterations ->setAdjusterListener (this);
   dcbIterations ->show();
   dcbEnhance = Gtk::manage (new Gtk::CheckButton(M("PREFERENCES_DCBENHANCE")));
   dcbOptions->pack_start(*dcbIterations);
   dcbOptions->pack_start(*dcbEnhance);
   pack_start( *dcbOptions, Gtk::PACK_SHRINK, 4);
   pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 4 );

   ccOptions = Gtk::manage (new Gtk::VBox ());
   ccOptions->set_border_width(4);
   ccSteps = Gtk::manage (new Adjuster (M("PREFERENCES_FALSECOLOR"),0,5,1,2 ));
   ccSteps->setAdjusterListener (this);
   ccSteps->show();
   pack_start( *ccSteps, Gtk::PACK_SHRINK, 4);

   methodconn = dmethod->signal_changed().connect( sigc::mem_fun(*this, &RawProcess::methodChanged) );
   dcbEnhconn = dcbEnhance->signal_toggled().connect ( sigc::mem_fun(*this, &RawProcess::dcbEnhanceChanged), true);
}


void RawProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
   disableListener ();
   methodconn.block (true);
   dcbEnhconn.block (true);

   dmethod->set_active(procparams::RAWParams::numMethods);
   for( size_t i=0; i< procparams::RAWParams::numMethods;i++)
	   if( pp->raw.dmethod == procparams::RAWParams::methodstring[i]){
		   dmethod->set_active(i);
		   break;
	   }

   dcbIterations->setValue (pp->raw.dcb_iterations);
   dcbEnhance->set_active(pp->raw.dcb_enhance);
   ccSteps->setValue (pp->raw.ccSteps);
   if (pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::dcb])
	   dcbOptions->show();
   else
	   dcbOptions->hide();

   if (pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::eahd] ||
	   pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::hphd] ||
	   pp->raw.dmethod == procparams::RAWParams::methodstring[procparams::RAWParams::vng4])
	   ccOptions->show();
   else
	   ccOptions->hide();

   lastDCBen = pp->raw.dcb_enhance;

   if(pedited ){
	   ccSteps->setEditedState (pedited->raw.ccSteps ? Edited : UnEdited);
	   dcbIterations->setEditedState ( pedited->raw.dcbIterations ? Edited : UnEdited);
	   dcbEnhance->set_inconsistent(!pedited->raw.dcbEnhance);
	   if( !pedited->raw.dmethod )
		   dmethod->set_active(procparams::RAWParams::numMethods); // No name
   }

   methodconn.block (false);
   dcbEnhconn.block (false);
   enableListener ();
}

void RawProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.ccSteps = (int)ccSteps->getValue();
	pp->raw.dcb_iterations = (int)dcbIterations->getValue();
	pp->raw.dcb_enhance = dcbEnhance->get_active();

	int currentRow = dmethod->get_active_row_number();
	if( currentRow>=0 && currentRow < procparams::RAWParams::numMethods)
		pp->raw.dmethod = procparams::RAWParams::methodstring[currentRow];

	if (pedited) {
		pedited->raw.ccSteps = ccSteps->getEditedState ();
		pedited->raw.dmethod = dmethod->get_active_row_number() != procparams::RAWParams::numMethods;
		pedited->raw.dcbIterations = dcbIterations->getEditedState ();
		pedited->raw.dcbEnhance = !dcbEnhance->get_inconsistent();
	}
}

void RawProcess::setBatchMode(bool batchMode)
{
   dmethod->set_active(procparams::RAWParams::numMethods); // No name
   dcbOptions->hide();
   ToolPanel::setBatchMode (batchMode);
   ccSteps->showEditedCB ();
   dcbIterations->showEditedCB ();
}

void RawProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	dcbIterations->setDefault( defParams->raw.dcb_iterations);
	ccSteps->setDefault (defParams->raw.ccSteps);
	if (pedited) {
		dcbIterations->setDefaultEditedState( pedited->raw.dcbIterations ? Edited : UnEdited);
		ccSteps->setDefaultEditedState(pedited->raw.ccSteps ? Edited : UnEdited);
	}else{
		dcbIterations->setDefaultEditedState( Irrelevant );
		ccSteps->setDefaultEditedState(Irrelevant );
	}
}

void RawProcess::adjusterChanged (Adjuster* a, double newval)
{
    if (listener)
        listener->panelChanged (EvDemosaic, Glib::ustring("params") );
}

void RawProcess::methodChanged ()
{
	int  curSelection = dmethod->get_active_row_number();
	if ( curSelection == procparams::RAWParams::dcb){
		dcbOptions->show();
	}else{
		dcbOptions->hide();
	}
	Glib::ustring s="";
	if( curSelection>=0 && curSelection < procparams::RAWParams::numMethods)
	    s = procparams::RAWParams::methodstring[curSelection];

    if (listener)
        listener->panelChanged (EvDemosaic, Glib::ustring(M("PREFERENCES_DMETHOD"))+ "="+ s);
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
        listener->panelChanged (EvDemosaic, Glib::ustring("params") );
}
