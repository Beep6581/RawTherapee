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
#include "xtransprocess.h"
#include "options.h"
#include "guiutils.h"
using namespace rtengine;
using namespace rtengine::procparams;

XTransProcess::XTransProcess () : FoldableToolPanel(this)
{
	Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
	hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_DMETHOD") +": ")),Gtk::PACK_SHRINK, 4);
	method = Gtk::manage (new MyComboBoxText ());
	for( size_t i=0; i<procparams::RAWParams::XTransSensor::numMethods;i++)
		method->append_text(procparams::RAWParams::XTransSensor::methodstring[i]);

	method->set_active(0);
	hb1->set_tooltip_markup (M("TP_RAW_SENSOR_XTRANS_DMETHOD_TOOLTIP"));

	hb1->pack_end (*method, Gtk::PACK_EXPAND_WIDGET, 4);
	pack_start( *hb1, Gtk::PACK_SHRINK, 4);

	pack_start( *Gtk::manage( new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0 );
	ccSteps = Gtk::manage (new Adjuster (M("TP_RAW_FALSECOLOR"),0,5,1,0 ));
	ccSteps->setAdjusterListener (this);
	if (ccSteps->delay < 1000) ccSteps->delay = 1000;
	ccSteps->show();
	pack_start( *ccSteps, Gtk::PACK_SHRINK, 4);

	methodconn = method->signal_changed().connect( sigc::mem_fun(*this, &XTransProcess::methodChanged) );
}


void XTransProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();
	methodconn.block (true);

	method->set_active(procparams::RAWParams::XTransSensor::numMethods);
	for( size_t i=0; i< procparams::RAWParams::XTransSensor::numMethods; i++)
		if( pp->raw.xtranssensor.method == procparams::RAWParams::XTransSensor::methodstring[i]) {
			method->set_active(i);
			oldSelection = i;
			break;
		}

	if(pedited ){
		ccSteps->setEditedState (pedited->raw.xtranssensor.ccSteps ? Edited : UnEdited);
		if( !pedited->raw.xtranssensor.method )
			method->set_active(procparams::RAWParams::XTransSensor::numMethods); // No name
	}

	ccSteps->setValue (pp->raw.xtranssensor.ccSteps);

	methodconn.block (false);

	enableListener ();
}

void XTransProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.xtranssensor.ccSteps = ccSteps->getIntValue();

	int currentRow = method->get_active_row_number();
	if( currentRow>=0 && currentRow < procparams::RAWParams::XTransSensor::numMethods)
		pp->raw.xtranssensor.method = procparams::RAWParams::XTransSensor::methodstring[currentRow];

	if (pedited) {
		pedited->raw.xtranssensor.method = method->get_active_row_number() != procparams::RAWParams::XTransSensor::numMethods;
		pedited->raw.xtranssensor.ccSteps = ccSteps->getEditedState ();
	}
}

void XTransProcess::setBatchMode(bool batchMode)
{
	method->append_text (M("GENERAL_UNCHANGED"));
	method->set_active(procparams::RAWParams::XTransSensor::numMethods); // No name
	ToolPanel::setBatchMode (batchMode);
	ccSteps->showEditedCB ();
}

void XTransProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	ccSteps->setDefault (defParams->raw.xtranssensor.ccSteps);
	if (pedited) {
		ccSteps->setDefaultEditedState(pedited->raw.xtranssensor.ccSteps ? Edited : UnEdited);
	}else{
		ccSteps->setDefaultEditedState(Irrelevant );
	}
}

void XTransProcess::adjusterChanged (Adjuster* a, double newval)
{
	if (listener) {
		if (a == ccSteps)
			listener->panelChanged (EvDemosaicFalseColorIter, a->getTextValue() );
	}
}

void XTransProcess::methodChanged ()
{
	int  curSelection = method->get_active_row_number();

	Glib::ustring methodName="";
	bool ppreq = false;
	if( curSelection>=0 && curSelection < procparams::RAWParams::XTransSensor::numMethods) {
		methodName = procparams::RAWParams::XTransSensor::methodstring[curSelection];
		if (curSelection == procparams::RAWParams::XTransSensor::mono || oldSelection == procparams::RAWParams::XTransSensor::mono) {
			ppreq = true;
		}
	}
	oldSelection = curSelection;

	if (listener)
		listener->panelChanged (ppreq ? EvDemosaicMethodPreProc : EvDemosaicMethod, methodName);
}
