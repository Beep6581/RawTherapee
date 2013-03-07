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

EdgePreservingDecompositionUI::EdgePreservingDecompositionUI () : Gtk::VBox(), FoldableToolPanel(this){

	set_border_width(4);

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	enabled->set_tooltip_markup (M("TP_EPD_TOOLTIP"));
	enabled->show ();
	pack_start (*enabled);

	Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
	hsep1->show ();
	pack_start (*hsep1);

	enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &EdgePreservingDecompositionUI::enabledChanged) );

	Strength					= Gtk::manage(new Adjuster (M("TP_EPD_STRENGTH"), -2.0, 2.0, 0.01, 0.25));
	EdgeStopping			= Gtk::manage(new Adjuster (M("TP_EPD_EDGESTOPPING"), 0.1, 4.0, 0.01, 1.4));
	Scale						= Gtk::manage(new Adjuster (M("TP_EPD_SCALE"), 0.1, 10.0, 0.01, 1.0));
	ReweightingIterates	= Gtk::manage(new Adjuster (M("TP_EPD_REWEIGHTINGITERATES"), 0, 9, 1, 0));

	Strength->setAdjusterListener(this);
	EdgeStopping->setAdjusterListener(this);
	Scale->setAdjusterListener(this);
	ReweightingIterates->setAdjusterListener(this);

	Strength->show();
	EdgeStopping->show();
	Scale->show();
	ReweightingIterates->show();

	pack_start(*Strength);
	pack_start(*EdgeStopping);
	pack_start(*Scale);
	pack_start(*ReweightingIterates);
}

void EdgePreservingDecompositionUI::read(const ProcParams *pp, const ParamsEdited *pedited){
	disableListener();

	if(pedited){
		Strength->setEditedState(pedited->edgePreservingDecompositionUI.Strength ? Edited : UnEdited);
		EdgeStopping->setEditedState(pedited->edgePreservingDecompositionUI.EdgeStopping ? Edited : UnEdited);
		Scale->setEditedState(pedited->edgePreservingDecompositionUI.Scale ? Edited : UnEdited);
		ReweightingIterates->setEditedState(pedited->edgePreservingDecompositionUI.ReweightingIterates ? Edited : UnEdited);

		enabled->set_inconsistent(!pedited->edgePreservingDecompositionUI.enabled);
	}

	enaConn.block(true);
	enabled->set_active(pp->edgePreservingDecompositionUI.enabled);
	enaConn.block (false);

	lastEnabled = pp->edgePreservingDecompositionUI.enabled;

	Strength->setValue(pp->edgePreservingDecompositionUI.Strength);
	EdgeStopping->setValue(pp->edgePreservingDecompositionUI.EdgeStopping);
	Scale->setValue(pp->edgePreservingDecompositionUI.Scale);
	ReweightingIterates->setValue(pp->edgePreservingDecompositionUI.ReweightingIterates);

	enableListener();
}

void EdgePreservingDecompositionUI::write(ProcParams *pp, ParamsEdited *pedited){
	pp->edgePreservingDecompositionUI.Strength = Strength->getValue();
	pp->edgePreservingDecompositionUI.EdgeStopping = EdgeStopping->getValue();
	pp->edgePreservingDecompositionUI.Scale = Scale->getValue();
	pp->edgePreservingDecompositionUI.ReweightingIterates = ReweightingIterates->getValue();
	pp->edgePreservingDecompositionUI.enabled = enabled->get_active();
	
	if(pedited){
		pedited->edgePreservingDecompositionUI.Strength = Strength->getEditedState();
		pedited->edgePreservingDecompositionUI.EdgeStopping = EdgeStopping->getEditedState();
		pedited->edgePreservingDecompositionUI.Scale = Scale->getEditedState();
		pedited->edgePreservingDecompositionUI.ReweightingIterates = ReweightingIterates->getEditedState();
		pedited->edgePreservingDecompositionUI.enabled = !enabled->get_inconsistent();
	}
}

void EdgePreservingDecompositionUI::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited){
	Strength->setDefault(defParams->edgePreservingDecompositionUI.Strength);
	EdgeStopping->setDefault(defParams->edgePreservingDecompositionUI.EdgeStopping);
	Scale->setDefault(defParams->edgePreservingDecompositionUI.Scale);
	ReweightingIterates->setDefault(defParams->edgePreservingDecompositionUI.ReweightingIterates);

	if(pedited){
		Strength->setDefaultEditedState(pedited->edgePreservingDecompositionUI.Strength ? Edited : UnEdited);
		EdgeStopping->setDefaultEditedState(pedited->edgePreservingDecompositionUI.EdgeStopping ? Edited : UnEdited);
		Scale->setDefaultEditedState(pedited->edgePreservingDecompositionUI.Scale ? Edited : UnEdited);
		ReweightingIterates->setDefaultEditedState(pedited->edgePreservingDecompositionUI.ReweightingIterates ? Edited : UnEdited);
	}else{
		Strength->setDefaultEditedState(Irrelevant);
		EdgeStopping->setDefaultEditedState(Irrelevant);
		Scale->setDefaultEditedState(Irrelevant);
		ReweightingIterates->setDefaultEditedState(Irrelevant);
	}
}

void EdgePreservingDecompositionUI::adjusterChanged(Adjuster* a, double newval){
	if(listener && enabled->get_active()){
		if(a == Strength)
			listener->panelChanged(EvEPDStrength, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
		else if(a == EdgeStopping)
			listener->panelChanged(EvEPDEdgeStopping, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
		else if(a == Scale) 
			listener->panelChanged(EvEPDScale, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
		else if(a == ReweightingIterates)
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

	Strength->showEditedCB();
	EdgeStopping->showEditedCB();
	Scale->showEditedCB();
	ReweightingIterates->showEditedCB();
}

