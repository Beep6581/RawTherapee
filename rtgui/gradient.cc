/*
 *  This file is part of RawTherapee.
 */
#include "gradient.h"
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

Gradient::Gradient () : Gtk::VBox(), FoldableToolPanel(this)
{
	set_border_width(4);

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	enaConn  = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Gradient::enabledChanged) );

	strength = Gtk::manage (new Adjuster (M("TP_GRADIENT_STRENGTH"), -5, 5, 0.01, 0));
	strength->set_tooltip_text (M("TP_GRADIENT_STRENGTH_TOOLTIP"));
	strength->setAdjusterListener (this);

	degree = Gtk::manage (new Adjuster (M("TP_GRADIENT_DEGREE"), -180, 180, 1, 0));
	degree->set_tooltip_text (M("TP_GRADIENT_DEGREE_TOOLTIP"));
	degree->setAdjusterListener (this);

	feather = Gtk::manage (new Adjuster (M("TP_GRADIENT_FEATHER"), 0, 100, 1, 25));
	feather->set_tooltip_text (M("TP_GRADIENT_FEATHER_TOOLTIP"));
	feather->setAdjusterListener (this);

	centerX = Gtk::manage (new Adjuster (M("TP_GRADIENT_CENTER_X"), -100, 100, 1, 0));
	centerX->set_tooltip_text (M("TP_GRADIENT_CENTER_X_TOOLTIP"));
	centerX->setAdjusterListener (this);

	centerY = Gtk::manage (new Adjuster (M("TP_GRADIENT_CENTER_Y"), -100, 100, 1, 0));
	centerY->set_tooltip_text (M("TP_GRADIENT_CENTER_Y_TOOLTIP"));
	centerY->setAdjusterListener (this);

	pack_start(*enabled);
	pack_start(*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_EXPAND_WIDGET, 4);
	pack_start (*strength);
	pack_start (*degree);
	pack_start (*feather);
	pack_start (*centerX);
	pack_start (*centerY);

	show_all();
}

void Gradient::read (const ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();

	if (pedited) {
		degree->setEditedState (pedited->gradient.degree ? Edited : UnEdited);
		feather->setEditedState (pedited->gradient.feather ? Edited : UnEdited);
		strength->setEditedState (pedited->gradient.strength ? Edited : UnEdited);
		centerX->setEditedState (pedited->gradient.centerX ? Edited : UnEdited);
		centerY->setEditedState (pedited->gradient.centerY ? Edited : UnEdited);
		enabled->set_inconsistent (!pedited->gradient.enabled);
	}

	enaConn.block (true);
	enabled->set_active (pp->gradient.enabled);
	enaConn.block (false);
	degree->setValue (pp->gradient.degree);
	feather->setValue (pp->gradient.feather);
	strength->setValue (pp->gradient.strength);
	centerX->setValue (pp->gradient.centerX);
	centerY->setValue (pp->gradient.centerY);

	lastEnabled = pp->gradient.enabled;

	enableListener ();
}

void Gradient::write (ProcParams* pp, ParamsEdited* pedited)
{
	pp->gradient.degree = degree->getValue ();
	pp->gradient.feather = feather->getIntValue ();
	pp->gradient.strength = strength->getValue ();
	pp->gradient.centerX = centerX->getIntValue ();
	pp->gradient.centerY = centerY->getIntValue ();
	pp->gradient.enabled = enabled->get_active();

	if (pedited) {
		pedited->gradient.degree = degree->getEditedState ();
		pedited->gradient.feather = feather->getEditedState ();
		pedited->gradient.strength = strength->getEditedState ();
		pedited->gradient.centerX = centerX->getEditedState ();
		pedited->gradient.centerY = centerY->getEditedState ();
		pedited->gradient.enabled = !enabled->get_inconsistent();
	}
}

void Gradient::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
	degree->setDefault (defParams->gradient.degree);
	feather->setDefault (defParams->gradient.feather);
	strength->setDefault (defParams->gradient.strength);
	centerX->setDefault (defParams->gradient.centerX);
	centerY->setDefault (defParams->gradient.centerY);

	if (pedited) {
		degree->setDefaultEditedState (pedited->gradient.degree ? Edited : UnEdited);
		feather->setDefaultEditedState (pedited->gradient.feather ? Edited : UnEdited);
		strength->setDefaultEditedState (pedited->gradient.strength ? Edited : UnEdited);
		centerX->setDefaultEditedState (pedited->gradient.centerX ? Edited : UnEdited);
		centerY->setDefaultEditedState (pedited->gradient.centerY ? Edited : UnEdited);
	} else {
		degree->setDefaultEditedState (Irrelevant);
		feather->setDefaultEditedState (Irrelevant);
		strength->setDefaultEditedState (Irrelevant);
		centerX->setDefaultEditedState (Irrelevant);
		centerY->setDefaultEditedState (Irrelevant);
	}
}

void Gradient::adjusterChanged (Adjuster* a, double newval) {

	if (listener && enabled->get_active()) {
		if (a == degree)
			listener->panelChanged (EvGradientDegree, degree->getTextValue());
		else if (a == feather)
			listener->panelChanged (EvGradientFeather, feather->getTextValue());
		else if (a == strength)
			listener->panelChanged (EvGradientStrength, strength->getTextValue());
		else if (a == centerX || a == centerY)
			listener->panelChanged (EvGradientCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
	}
}

void Gradient::enabledChanged () {

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
		if (enabled->get_active())
			listener->panelChanged (EvGradientEnabled, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvGradientEnabled, M("GENERAL_DISABLED"));
	}
}

void Gradient::setAdjusterBehavior (bool degreeadd, bool featheradd, bool strengthadd, bool centeradd)
{
	degree->setAddMode(degreeadd);
	feather->setAddMode(featheradd);
	strength->setAddMode(strengthadd);
	centerX->setAddMode(centeradd);
	centerY->setAddMode(centeradd);
}

void Gradient::trimValues (rtengine::procparams::ProcParams* pp)
{
	degree->trimValue(pp->gradient.degree);
	feather->trimValue(pp->gradient.feather);
	strength->trimValue(pp->gradient.strength);
	centerX->trimValue(pp->gradient.centerX);
	centerY->trimValue(pp->gradient.centerY);
}

void Gradient::setBatchMode (bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	degree->showEditedCB ();
	feather->showEditedCB ();
	strength->showEditedCB ();
	centerX->showEditedCB ();
	centerY->showEditedCB ();
}
