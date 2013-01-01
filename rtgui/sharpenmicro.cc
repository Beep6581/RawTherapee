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
#include "sharpenmicro.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;


SharpenMicro::SharpenMicro () : Gtk::VBox(), FoldableToolPanel(this) {

	set_border_width(4);

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (true);
	enabled->set_tooltip_markup (M("TP_SHARPENING_TOOLTIP"));
	
	pack_start(*enabled, Gtk::PACK_SHRINK, 0);
	enabled->show ();

	Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
	hsep1->show ();
	pack_start (*hsep1);

	amount= Gtk::manage(new Adjuster (M("TP_SHARPENMICRO_AMOUNT"),0,100,1,20));
	amount->setAdjusterListener (this);
	if (amount->delay < 1000) amount->delay = 1000;
	amount->show();
	uniformity= Gtk::manage(new Adjuster (M("TP_SHARPENMICRO_UNIFORMITY"),0,100,10,50));

	uniformity->setAdjusterListener (this);
	if (uniformity->delay < 1000) uniformity->delay = 1000;
	uniformity->show();
	matrix = Gtk::manage (new Gtk::CheckButton (M("TP_SHARPENMICRO_MATRIX")));
	matrix->set_active (true);
	pack_start(*matrix, Gtk::PACK_SHRINK, 0);
	matrix->show ();

	pack_start( *amount, Gtk::PACK_SHRINK, 0);
	pack_start( *uniformity, Gtk::PACK_SHRINK, 0);

	show ();

	enaconn    = enabled->signal_toggled().connect( sigc::mem_fun(*this, &SharpenMicro::enabled_toggled) );
	matrixconn = matrix->signal_toggled().connect( sigc::mem_fun(*this, &SharpenMicro::matrix_toggled) );
}

void SharpenMicro::read(const ProcParams* pp, const ParamsEdited* pedited) {
	disableListener ();

	if(pedited ){
		enabled->set_inconsistent  (!pedited->sharpenMicro.enabled);
		matrix->set_inconsistent   (!pedited->sharpenMicro.matrix);
		amount->setEditedState     (pedited->sharpenMicro.amount ? Edited : UnEdited);
		uniformity->setEditedState (pedited->sharpenMicro.uniformity ? Edited : UnEdited);
	}
	enaconn.block (true);
	enabled->set_active (pp->sharpenMicro.enabled);
	enaconn.block (false);
	lastEnabled = pp->sharpenMicro.enabled;

	matrixconn.block (true);
	matrix->set_active (pp->sharpenMicro.matrix);
	matrixconn.block (false);
	lastmatrix = pp->sharpenMicro.matrix;

	amount->setValue     (pp->sharpenMicro.amount);
	uniformity->setValue (pp->sharpenMicro.uniformity);

	enableListener ();
}

void SharpenMicro::write( ProcParams* pp, ParamsEdited* pedited) {
	pp->sharpenMicro.enabled    = enabled->get_active ();
	pp->sharpenMicro.matrix     = matrix->get_active ();
	pp->sharpenMicro.amount     = amount->getValue ();
	pp->sharpenMicro.uniformity = uniformity->getValue ();

	if (pedited) {
		pedited->sharpenMicro.enabled    = !enabled->get_inconsistent();
		pedited->sharpenMicro.matrix     = !matrix->get_inconsistent();
		pedited->sharpenMicro.amount     = amount->getEditedState ();
		pedited->sharpenMicro.uniformity = uniformity->getEditedState ();
	}
}

void SharpenMicro::enabled_toggled () {
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
			listener->panelChanged (EvSharpenMicroEnabled, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvSharpenMicroEnabled, M("GENERAL_DISABLED"));
			
	}
}

void SharpenMicro::matrix_toggled () {
	if (batchMode) {
		if (matrix->get_inconsistent()) {
			matrix->set_inconsistent (false);
			matrixconn.block (true);
			matrix->set_active (false);
			matrixconn.block (false);
		}
		else if (lastmatrix)
			matrix->set_inconsistent (true);

		lastmatrix = matrix->get_active ();
	}

	if (listener && enabled->get_active ()) {
		if (matrix->get_active ())
			listener->panelChanged (EvSharpenMicroMatrix, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvSharpenMicroMatrix, M("GENERAL_DISABLED"));
	}
}

void SharpenMicro::adjusterChanged (Adjuster* a, double newval) {
	if (listener && enabled->get_active()) {
		Glib::ustring value = a->getTextValue();
		if (a == amount)
			listener->panelChanged (EvSharpenMicroAmount,     value );
		else if (a == uniformity)
			listener->panelChanged (EvSharpenMicroUniformity, value );
	}
}

void SharpenMicro::setBatchMode(bool batchMode) {
	amount->showEditedCB     ();
	uniformity->showEditedCB ();
}

void SharpenMicro::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited) {
	amount->setDefault     (defParams->sharpenMicro.amount);
	uniformity->setDefault (defParams->sharpenMicro.uniformity);

	if (pedited) {
		amount->setDefaultEditedState     (pedited->sharpenMicro.amount ? Edited : UnEdited);
		uniformity->setDefaultEditedState (pedited->sharpenMicro.uniformity ? Edited : UnEdited);
	} else {
		amount->setDefaultEditedState     (Irrelevant);
		uniformity->setDefaultEditedState (Irrelevant);
	}
}

void SharpenMicro::setAdjusterBehavior (bool amountadd, bool uniformityadd ) {
	amount->setAddMode     (amountadd);
	uniformity->setAddMode (uniformityadd);
}

void SharpenMicro::trimValues (ProcParams* pp) {
	amount->trimValue     (pp->sharpenMicro.amount);
	uniformity->trimValue (pp->sharpenMicro.uniformity);
}
