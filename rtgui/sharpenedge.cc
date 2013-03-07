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
#include "sharpenedge.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;


SharpenEdge::SharpenEdge () : Gtk::VBox(), FoldableToolPanel(this) {

	set_border_width(4);

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (true);
	pack_start(*enabled, Gtk::PACK_SHRINK, 0);

	Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
	hsep1->show ();
	pack_start (*hsep1);

	passes = Gtk::manage(new Adjuster (M("TP_SHARPENEDGE_PASSES"),1,4,1,2));
	passes->setAdjusterListener (this);
	if (passes->delay < 1000) passes->delay = 1000;
	amount = Gtk::manage(new Adjuster (M("TP_SHARPENEDGE_AMOUNT"),0,100,1,50));
	amount->setAdjusterListener (this);
	if (amount->delay < 1000) amount->delay = 1000;

	threechannels = Gtk::manage(new Gtk::CheckButton((M("TP_SHARPENEDGE_THREE"))));// L + a + b
	threechannels->set_active (false);
	pack_start( *passes, Gtk::PACK_SHRINK, 0);//passes
	pack_start( *amount, Gtk::PACK_SHRINK, 0);//amount
	pack_start( *threechannels, Gtk::PACK_SHRINK, 0);//one or 3 channels Lab

	show ();

	enaconn       = enabled->signal_toggled().connect( sigc::mem_fun(*this, &SharpenEdge::enabled_toggled) );
	chanthreeconn = threechannels->signal_toggled().connect( sigc::mem_fun(*this, &SharpenEdge::chanthree_toggled) );
}

void SharpenEdge::read(const ProcParams* pp, const ParamsEdited* pedited) {
	disableListener ();

	if(pedited ){
		passes->setEditedState          (pedited->sharpenEdge.passes ? Edited : UnEdited);
		amount->setEditedState          (pedited->sharpenEdge.amount ? Edited : UnEdited);
		enabled->set_inconsistent       (!pedited->sharpenEdge.enabled);
		threechannels->set_inconsistent (!pedited->sharpenEdge.threechannels);
	}
	enaconn.block (true);
	enabled->set_active (pp->sharpenEdge.enabled);
	enaconn.block (false);
	lastEnabled = pp->sharpenEdge.enabled;

	chanthreeconn.block (true);
	threechannels->set_active (pp->sharpenEdge.threechannels);
	chanthreeconn.block (false);
	lastchanthree = pp->sharpenEdge.threechannels;

	passes->setValue (pp->sharpenEdge.passes);
	amount->setValue (pp->sharpenEdge.amount);
	
	enableListener ();
}

void SharpenEdge::write( ProcParams* pp, ParamsEdited* pedited) {
	pp->sharpenEdge.enabled       = enabled->get_active ();
	pp->sharpenEdge.passes        = (int)passes->getValue();
	pp->sharpenEdge.amount        = amount->getValue ();
	pp->sharpenEdge.threechannels = threechannels->get_active ();

	if (pedited) {
		pedited->sharpenEdge.enabled       = !enabled->get_inconsistent();
		pedited->sharpenEdge.passes        = passes->getEditedState ();
		pedited->sharpenEdge.amount        = amount->getEditedState ();
		pedited->sharpenEdge.threechannels = !threechannels->get_inconsistent();
	}

}

void SharpenEdge::enabled_toggled () {

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
			listener->panelChanged (EvSharpenEdgeEnabled, M("GENERAL_ENABLED"));
			//listener->panelChanged (EvMLunifor, M("GENERAL_ENABLED"));
			
		else
			listener->panelChanged (EvSharpenEdgeEnabled, M("GENERAL_DISABLED"));
	}
}

void SharpenEdge::chanthree_toggled () {

	if (batchMode) {
		if (threechannels->get_inconsistent()) {
			threechannels->set_inconsistent (false);
			chanthreeconn.block (true);
			threechannels->set_active (false);
			chanthreeconn.block (false);
		}
		else if (lastchanthree)
			threechannels->set_inconsistent (true);

		lastchanthree = threechannels->get_active ();
	}

	if (listener && enabled->get_active ()) {
		if (threechannels->get_active ())
			listener->panelChanged (EvSharpenEdgeThreechannels, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvSharpenEdgeThreechannels, M("GENERAL_DISABLED"));
	}
}

void SharpenEdge::adjusterChanged (Adjuster* a, double newval) {
	if (listener && enabled->get_active()) {
		Glib::ustring value = a->getTextValue();
		
		if (a == passes )
			listener->panelChanged (EvSharpenEdgePasses,   value );
		else if (a == amount)
			listener->panelChanged (EvSharpenEdgeAmount, value );
	}
}

void SharpenEdge::setBatchMode(bool batchMode) {
	passes->showEditedCB   ();
	amount->showEditedCB ();
}

void SharpenEdge::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited) {
	passes->setDefault (defParams->sharpenEdge.passes);
	amount->setDefault (defParams->sharpenEdge.amount);

	if (pedited) {
		passes->setDefaultEditedState   (pedited->sharpenEdge.passes ? Edited : UnEdited);
		amount->setDefaultEditedState (pedited->sharpenEdge.amount ? Edited : UnEdited);
		
	} else {
		passes->setDefaultEditedState   (Irrelevant);
		amount->setDefaultEditedState (Irrelevant);
		
	}
}

void SharpenEdge::setAdjusterBehavior (bool amountadd, bool passadd) {
	amount->setAddMode (amountadd);
	passes->setAddMode (passadd);
}

void SharpenEdge::trimValues (ProcParams* pp) {
	amount->trimValue (pp->sharpenEdge.amount);
	passes->trimValue (pp->sharpenEdge.passes);

}
