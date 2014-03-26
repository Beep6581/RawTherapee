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
#include "epd.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

EdgePreservingDecompositionUI::EdgePreservingDecompositionUI () : FoldableToolPanel(this){

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	enabled->set_tooltip_markup (M("TP_EPD_TOOLTIP"));
	enabled->show ();
	pack_start (*enabled);

	Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
	hsep1->show ();
	pack_start (*hsep1);

	enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &EdgePreservingDecompositionUI::enabledChanged) );

	strength = Gtk::manage(new Adjuster (M("TP_EPD_STRENGTH"), -2.0, 2.0, 0.01, 0.25));
	edgeStopping = Gtk::manage(new Adjuster (M("TP_EPD_EDGESTOPPING"), 0.1, 4.0, 0.01, 1.4));
	scale = Gtk::manage(new Adjuster (M("TP_EPD_SCALE"), 0.1, 10.0, 0.01, 1.0));
	reweightingIterates	= Gtk::manage(new Adjuster (M("TP_EPD_REWEIGHTINGITERATES"), 0, 9, 1, 0));

	strength->setAdjusterListener(this);
	edgeStopping->setAdjusterListener(this);
	scale->setAdjusterListener(this);
	reweightingIterates->setAdjusterListener(this);

	strength->show();
	edgeStopping->show();
	scale->show();
	reweightingIterates->show();

	pack_start(*strength);
	pack_start(*edgeStopping);
	pack_start(*scale);
	pack_start(*reweightingIterates);
}

void EdgePreservingDecompositionUI::read(const ProcParams *pp, const ParamsEdited *pedited){
	disableListener();

	if(pedited){
		strength->setEditedState(pedited->epd.strength ? Edited : UnEdited);
		edgeStopping->setEditedState(pedited->epd.edgeStopping ? Edited : UnEdited);
		scale->setEditedState(pedited->epd.scale ? Edited : UnEdited);
		reweightingIterates->setEditedState(pedited->epd.reweightingIterates ? Edited : UnEdited);

		enabled->set_inconsistent(!pedited->epd.enabled);
	}

	enaConn.block(true);
	enabled->set_active(pp->epd.enabled);
	enaConn.block (false);

	lastEnabled = pp->epd.enabled;

	strength->setValue(pp->epd.strength);
	edgeStopping->setValue(pp->epd.edgeStopping);
	scale->setValue(pp->epd.scale);
	reweightingIterates->setValue(pp->epd.reweightingIterates);

	enableListener();
}

void EdgePreservingDecompositionUI::write(ProcParams *pp, ParamsEdited *pedited){
	pp->epd.strength = strength->getValue();
	pp->epd.edgeStopping = edgeStopping->getValue();
	pp->epd.scale = scale->getValue();
	pp->epd.reweightingIterates = reweightingIterates->getValue();
	pp->epd.enabled = enabled->get_active();
	
	if(pedited){
		pedited->epd.strength = strength->getEditedState();
		pedited->epd.edgeStopping = edgeStopping->getEditedState();
		pedited->epd.scale = scale->getEditedState();
		pedited->epd.reweightingIterates = reweightingIterates->getEditedState();
		pedited->epd.enabled = !enabled->get_inconsistent();
	}
}

void EdgePreservingDecompositionUI::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited){
	strength->setDefault(defParams->epd.strength);
	edgeStopping->setDefault(defParams->epd.edgeStopping);
	scale->setDefault(defParams->epd.scale);
	reweightingIterates->setDefault(defParams->epd.reweightingIterates);

	if(pedited){
		strength->setDefaultEditedState(pedited->epd.strength ? Edited : UnEdited);
		edgeStopping->setDefaultEditedState(pedited->epd.edgeStopping ? Edited : UnEdited);
		scale->setDefaultEditedState(pedited->epd.scale ? Edited : UnEdited);
		reweightingIterates->setDefaultEditedState(pedited->epd.reweightingIterates ? Edited : UnEdited);
	}else{
		strength->setDefaultEditedState(Irrelevant);
		edgeStopping->setDefaultEditedState(Irrelevant);
		scale->setDefaultEditedState(Irrelevant);
		reweightingIterates->setDefaultEditedState(Irrelevant);
	}
}

void EdgePreservingDecompositionUI::adjusterChanged(Adjuster* a, double newval){
	if(listener && enabled->get_active()){
		if(a == strength)
			listener->panelChanged(EvEPDStrength, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
		else if(a == edgeStopping)
			listener->panelChanged(EvEPDEdgeStopping, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
		else if(a == scale)
			listener->panelChanged(EvEPDScale, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
		else if(a == reweightingIterates)
			listener->panelChanged(EvEPDReweightingIterates, Glib::ustring::format((int)a->getValue()));
	}
}

void EdgePreservingDecompositionUI::enabledChanged(){
	if(batchMode){
		if(enabled->get_inconsistent()){
			enabled->set_inconsistent (false);
			enaConn.block (true);
			enabled->set_active (false);
			enaConn.block (false);
		}
		else if(lastEnabled)
			enabled->set_inconsistent (true);

		lastEnabled = enabled->get_active ();
	}
	
	if(listener){
		if(enabled->get_active ())
			listener->panelChanged (EvEPDEnabled, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvEPDEnabled, M("GENERAL_DISABLED"));
	}
}

void EdgePreservingDecompositionUI::setBatchMode(bool batchMode){
	ToolPanel::setBatchMode(batchMode);

	strength->showEditedCB();
	edgeStopping->showEditedCB();
	scale->showEditedCB();
	reweightingIterates->showEditedCB();
}

