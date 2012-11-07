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
#include "rawexposure.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

RAWExposure::RAWExposure () : Gtk::VBox(), FoldableToolPanel(this)
{
	set_border_width(4);

	PexPos = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_LINEAR"),0.1,16.0,0.01,1));
	PexPos->setAdjusterListener (this);
	if (PexPos->delay < 1000) PexPos->delay = 1000;
	PexPos->show();
	PexPreser = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_PRESER"),0,2.5,0.1,0));
	PexPreser->setAdjusterListener (this);
	if (PexPreser->delay < 1000) PexPreser->delay = 1000;
	PexPreser->show();
	PexBlackone = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACKONE"),-2048,2048,0.1,0));//black level
	PexBlackone->setAdjusterListener (this);
	if (PexBlackone->delay < 1000) PexBlackone->delay = 1000;
	PexBlackone->show();
	PexBlacktwo = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACKTWO"),-2048,2048,0.1,0));//black level
	PexBlacktwo->setAdjusterListener (this);
	if (PexBlacktwo->delay < 1000) PexBlacktwo->delay = 1000;
	PexBlacktwo->show();
	PexBlackthree = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACKTHREE"),-2048,2048,0.1,0));//black level
	PexBlackthree->setAdjusterListener (this);
	if (PexBlackthree->delay < 1000) PexBlackthree->delay = 1000;
	PexBlackthree->show();
	PexBlackzero = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACKZERO"),-2048,2048,0.1,0));//black level
	PexBlackzero->setAdjusterListener (this);
	if (PexBlackzero->delay < 1000) PexBlackzero->delay = 1000;
	PexBlackzero->show();
	PextwoGreen = Gtk::manage(new Gtk::CheckButton((M("TP_RAWEXPOS_TWOGREEN"))));// two green	
	PextwoGreen->set_active (true);
	greenconn = PextwoGreen->signal_toggled().connect ( sigc::mem_fun(*this, &RAWExposure::GreenChanged));

	

	pack_start( *PexPos, Gtk::PACK_SHRINK, 4);//exposi
	pack_start( *PexPreser, Gtk::PACK_SHRINK, 4);
	pack_start( *PexBlackone, Gtk::PACK_SHRINK, 4);//black
	pack_start( *PexBlackzero, Gtk::PACK_SHRINK, 4);//black	
	pack_start( *PexBlacktwo, Gtk::PACK_SHRINK, 4);//black
	pack_start( *PexBlackthree, Gtk::PACK_SHRINK, 4);//black
	pack_start( *PextwoGreen, Gtk::PACK_SHRINK, 4);//black 2 green
	
}

void RAWExposure::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();

	if(pedited ){
		PexPos->setEditedState( pedited->raw.exPos ? Edited : UnEdited );
		PexPreser->setEditedState( pedited->raw.exPreser ? Edited : UnEdited );
		PexBlackzero->setEditedState( pedited->raw.exBlackzero ? Edited : UnEdited );
		PexBlackone->setEditedState( pedited->raw.exBlackone ? Edited : UnEdited );
		PexBlacktwo->setEditedState( pedited->raw.exBlacktwo ? Edited : UnEdited );
		PexBlackthree->setEditedState( pedited->raw.exBlackthree ? Edited : UnEdited );
		
	}
	greenconn.block (true);
    PextwoGreen->set_active (pp->raw.twogreen);
    greenconn.block (false);
	lastPextwoGreen = pp->raw.twogreen;
	
	
	PexPos->setValue (pp->raw.expos);
	PexPreser->setValue (pp->raw.preser);//exposi
	PexBlackzero->setValue (pp->raw.blackzero);//black
	PexBlackone->setValue (pp->raw.blackone);//black
	PexBlacktwo->setValue (pp->raw.blacktwo);//black
	
	if(!PextwoGreen->get_active())PexBlackthree->setValue (pp->raw.blackthree);else PexBlackthree->setValue (PexBlackzero->getValue());

	enableListener ();
}

void RAWExposure::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.expos = PexPos->getValue();
	pp->raw.preser = PexPreser->getValue();//exposi
	pp->raw.blackzero = PexBlackzero->getValue();// black
	pp->raw.blackone = PexBlackone->getValue();// black
	pp->raw.blacktwo = PexBlacktwo->getValue();// black 
	pp->raw.twogreen=PextwoGreen->get_active();
	if(PextwoGreen->get_active()){pp->raw.blackthree=pp->raw.blackzero;} else {pp->raw.blackthree = PexBlackthree->getValue();}// active or desactive 2 green together
	

	if (pedited) {
		pedited->raw.exPos = PexPos->getEditedState ();
		pedited->raw.exPreser = PexPreser->getEditedState ();//exposi
		pedited->raw.exBlackzero = PexBlackzero->getEditedState ();//black
		pedited->raw.exBlackone = PexBlackone->getEditedState ();//black
		pedited->raw.exBlacktwo = PexBlacktwo->getEditedState ();//black
		pedited->raw.exBlackthree = PexBlackthree->getEditedState ();//black
		pedited->raw.exTwoGreen =!PextwoGreen->get_inconsistent();

		
	}

}

void RAWExposure::adjusterChanged (Adjuster* a, double newval)
{
	if (listener) {
		Glib::ustring value = a->getTextValue();
		{ 
		
		if (a == PexPos ) 
			listener->panelChanged (EvPreProcessExpCorrLinear,  value );
		else if (a == PexPreser && ABS(PexPos->getValue()-1.0)>0.0001)  // update takes long, only do it if it would have an effect
			listener->panelChanged (EvPreProcessExpCorrPH,  value );
		else if (a == PexBlackzero)	{if(!PextwoGreen->get_active())
			listener->panelChanged (EvPreProcessExpBlackzero,  value ); else {listener->panelChanged (EvPreProcessExpBlackzero,  value );PexBlackthree->setValue (PexBlackzero->getValue());}}
		else if (a == PexBlackone)	
			listener->panelChanged (EvPreProcessExpBlackone,  value );
		else if (a == PexBlacktwo)	
			listener->panelChanged (EvPreProcessExpBlacktwo,  value );
		else if (a == PexBlackthree)	{if(!PextwoGreen->get_active())
			listener->panelChanged (EvPreProcessExpBlackthree,  value ); else {listener->panelChanged (EvPreProcessExpBlackthree,  value );PexBlackzero->setValue (PexBlackthree->getValue());}}
			}
			
	}
}
void RAWExposure::GreenChanged() {
   if (batchMode) {
        if (PextwoGreen->get_inconsistent()) {
            PextwoGreen->set_inconsistent (false);
            greenconn.block (true);
            PextwoGreen->set_active (false);
            greenconn.block (false);
        }
        else if (lastPextwoGreen)
            PextwoGreen->set_inconsistent (true);
        lastPextwoGreen = PextwoGreen->get_active ();
    }
    
    if (listener) {
        if (PextwoGreen->get_active())
           { listener->panelChanged (EvPreProcessExptwoGreen, M("GENERAL_ENABLED"));
				PexBlackthree->setValue (PexBlackzero->getValue());//two green together
			}
		
        else
           { listener->panelChanged (EvPreProcessExptwoGreen, M("GENERAL_DISABLED"));
			} 
		
    }

}

void RAWExposure::setBatchMode(bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	PexPos->showEditedCB ();
	PexPreser->showEditedCB ();//exposure
	PexBlackzero->showEditedCB ();//black
	PexBlackone->showEditedCB ();//black
	PexBlacktwo->showEditedCB ();//black
	PexBlackthree->showEditedCB ();//black
	
}

void RAWExposure::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	PexPos->setDefault( defParams->raw.expos);
	PexPreser->setDefault( defParams->raw.preser);
	PexBlackzero->setDefault( defParams->raw.blackzero);
	PexBlackone->setDefault( defParams->raw.blackone);
	PexBlacktwo->setDefault( defParams->raw.blacktwo);
	PexBlackthree->setDefault( defParams->raw.blackthree);

	if (pedited) {
		PexPos->setDefaultEditedState( pedited->raw.exPos ? Edited : UnEdited);
		PexPreser->setDefaultEditedState( pedited->raw.exPreser ? Edited : UnEdited);
		PexBlackzero->setDefaultEditedState( pedited->raw.exBlackzero ? Edited : UnEdited);
		PexBlackone->setDefaultEditedState( pedited->raw.exBlackone ? Edited : UnEdited);
		PexBlacktwo->setDefaultEditedState( pedited->raw.exBlacktwo ? Edited : UnEdited);
		PexBlackthree->setDefaultEditedState( pedited->raw.exBlackthree ? Edited : UnEdited);
		
	} else {
		PexPos->setDefaultEditedState( Irrelevant );
		PexPreser->setDefaultEditedState( Irrelevant );
		PexBlackzero->setDefaultEditedState( Irrelevant );
		PexBlackone->setDefaultEditedState( Irrelevant );
		PexBlacktwo->setDefaultEditedState( Irrelevant );
		PexBlackthree->setDefaultEditedState( Irrelevant );
		
	}
}

void RAWExposure::setAdjusterBehavior (bool pexposadd, bool pexpreseradd, bool pexblackadd) {

	PexPos->setAddMode(pexposadd);
	PexPreser->setAddMode(pexpreseradd);
	PexBlackzero->setAddMode(pexblackadd);
	PexBlackone->setAddMode(pexblackadd);
	PexBlacktwo->setAddMode(pexblackadd);
	PexBlackthree->setAddMode(pexblackadd);
}

void RAWExposure::trimValues (rtengine::procparams::ProcParams* pp) {

	PexPos->trimValue(pp->raw.expos);
	PexPreser->trimValue(pp->raw.preser);
	PexBlackzero->trimValue(pp->raw.blackzero);
	PexBlackone->trimValue(pp->raw.blackone);
	PexBlacktwo->trimValue(pp->raw.blacktwo);
	PexBlackthree->trimValue(pp->raw.blackthree);
}
