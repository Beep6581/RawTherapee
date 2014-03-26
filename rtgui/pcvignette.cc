/*
 *  This file is part of RawTherapee.
 */
#include "pcvignette.h"
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

PCVignette::PCVignette () : FoldableToolPanel(this)
{
	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	enaConn  = enabled->signal_toggled().connect( sigc::mem_fun(*this, &PCVignette::enabledChanged) );

	strength = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_STRENGTH"), -6, 6, 0.01, 0));
	strength->set_tooltip_text (M("TP_PCVIGNETTE_STRENGTH_TOOLTIP"));
	strength->setAdjusterListener (this);

	feather = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_FEATHER"), 0, 100, 1, 50));
	feather->set_tooltip_text (M("TP_PCVIGNETTE_FEATHER_TOOLTIP"));
	feather->setAdjusterListener (this);

	roundness = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_ROUNDNESS"), 0, 100, 1, 50));
	roundness->set_tooltip_text (M("TP_PCVIGNETTE_ROUNDNESS_TOOLTIP"));
	roundness->setAdjusterListener (this);

	pack_start(*enabled);
	pack_start(*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_EXPAND_WIDGET, 4);
	pack_start (*strength);
	pack_start (*feather);
	pack_start (*roundness);

	show_all();
}

void PCVignette::read (const ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();

	if (pedited) {
		strength->setEditedState (pedited->pcvignette.strength ? Edited : UnEdited);
		feather->setEditedState (pedited->pcvignette.feather ? Edited : UnEdited);
		roundness->setEditedState (pedited->pcvignette.roundness ? Edited : UnEdited);
		enabled->set_inconsistent (!pedited->pcvignette.enabled);
	}

	enaConn.block (true);
	enabled->set_active (pp->pcvignette.enabled);
	enaConn.block (false);
	strength->setValue (pp->pcvignette.strength);
	feather->setValue (pp->pcvignette.feather);
	roundness->setValue (pp->pcvignette.roundness);

	lastEnabled = pp->pcvignette.enabled;

	enableListener ();
}

void PCVignette::write (ProcParams* pp, ParamsEdited* pedited)
{
	pp->pcvignette.strength = strength->getValue ();
	pp->pcvignette.feather = feather->getIntValue ();
	pp->pcvignette.roundness = roundness->getIntValue ();
	pp->pcvignette.enabled = enabled->get_active();

	if (pedited) {
		pedited->pcvignette.strength = strength->getEditedState ();
		pedited->pcvignette.feather = feather->getEditedState ();
		pedited->pcvignette.roundness = roundness->getEditedState ();
		pedited->pcvignette.enabled = !enabled->get_inconsistent();
	}
}

void PCVignette::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
	strength->setDefault (defParams->pcvignette.strength);
	feather->setDefault (defParams->pcvignette.feather);
	roundness->setDefault (defParams->pcvignette.roundness);

	if (pedited) {
		strength->setDefaultEditedState (pedited->pcvignette.strength ? Edited : UnEdited);
		feather->setDefaultEditedState (pedited->pcvignette.feather ? Edited : UnEdited);
		roundness->setDefaultEditedState (pedited->pcvignette.roundness ? Edited : UnEdited);
	} else {
		strength->setDefaultEditedState (Irrelevant);
		feather->setDefaultEditedState (Irrelevant);
		roundness->setDefaultEditedState (Irrelevant);
	}
}

void PCVignette::adjusterChanged (Adjuster* a, double newval) {

	if (listener && enabled->get_active()) {
		if (a == strength)
			listener->panelChanged (EvPCVignetteStrength, strength->getTextValue());
		else if (a == feather)
			listener->panelChanged (EvPCVignetteFeather, feather->getTextValue());
		else if (a == roundness)
			listener->panelChanged (EvPCVignetteRoundness, roundness->getTextValue());
	}
}

void PCVignette::enabledChanged () {

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
            listener->panelChanged (EvPCVignetteEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvPCVignetteEnabled, M("GENERAL_DISABLED"));
    }
}

void PCVignette::setAdjusterBehavior (bool strengthadd, bool featheradd, bool roundnessadd)
{
	strength->setAddMode(strengthadd);
	feather->setAddMode(featheradd);
	roundness->setAddMode(roundnessadd);
}

void PCVignette::trimValues (rtengine::procparams::ProcParams* pp)
{
	strength->trimValue(pp->pcvignette.strength);
	feather->trimValue(pp->pcvignette.feather);
	roundness->trimValue(pp->pcvignette.roundness);
}

void PCVignette::setBatchMode (bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	strength->showEditedCB ();
	feather->showEditedCB ();
	roundness->showEditedCB ();
}
