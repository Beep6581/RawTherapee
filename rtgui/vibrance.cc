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

#include "vibrance.h"

using namespace rtengine;
using namespace rtengine::procparams;

Vibrance::Vibrance () : Gtk::VBox(), FoldableToolPanel(this) {

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	pack_start(*enabled, Gtk::PACK_SHRINK, 0);

	pastels = Gtk::manage(new Adjuster (M("TP_VIBRANCE_PASTELS"),-100,100,5,0));
	pastels->setAdjusterListener (this);
	//if (pastels->delay < 1000) pastels->delay = 1000;
	pack_start( *pastels, Gtk::PACK_SHRINK, 0);

	saturated = Gtk::manage(new Adjuster (M("TP_VIBRANCE_SATURATED"),-100,100,5,0));
	saturated->setAdjusterListener (this);
	saturated->set_sensitive(false);
	//if (saturated->delay < 1000) saturated->delay = 1000;
	pack_start( *saturated, Gtk::PACK_SHRINK, 0);

	psThreshold = Gtk::manage (new ThresholdAdjuster (M("TP_VIBRANCE_PSTHRESHOLD"), 0., 100., 75., 75., 0, false));
	psThreshold->setAdjusterListener (this);
	psThreshold->set_sensitive(false);
	//if (psThreshold->delay < 1000) psThreshold->delay = 1000;
	pack_start( *psThreshold, Gtk::PACK_SHRINK, 0);

	protectSkins = Gtk::manage (new Gtk::CheckButton (M("TP_VIBRANCE_PROTECTSKINS")));
	protectSkins->set_active (true);
	pack_start(*protectSkins, Gtk::PACK_SHRINK, 0);

	avoidColorShift = Gtk::manage (new Gtk::CheckButton (M("TP_VIBRANCE_AVOIDCOLORSHIFT")));
	avoidColorShift->set_active (true);
	pack_start(*avoidColorShift, Gtk::PACK_SHRINK, 0);

	pastSatTog = Gtk::manage (new Gtk::CheckButton (M("TP_VIBRANCE_PASTSATTOG")));
	pastSatTog->set_active (true);
	pack_start(*pastSatTog, Gtk::PACK_SHRINK, 0);

	show ();

	enaconn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Vibrance::enabled_toggled) );
	pskinsconn = protectSkins->signal_toggled().connect( sigc::mem_fun(*this, &Vibrance::protectskins_toggled) );
	ashiftconn = avoidColorShift->signal_toggled().connect( sigc::mem_fun(*this, &Vibrance::avoidcolorshift_toggled) );
	pastsattogconn = pastSatTog->signal_toggled().connect( sigc::mem_fun(*this, &Vibrance::pastsattog_toggled) );
}

void Vibrance::read(const ProcParams* pp, const ParamsEdited* pedited) {
	disableListener ();

	if(pedited ){
		enabled->set_inconsistent         (!pedited->vibrance.enabled);
		pastels->setEditedState           (pedited->vibrance.pastels ? Edited : UnEdited);
		saturated->setEditedState         (pedited->vibrance.saturated ? Edited : UnEdited);
		psThreshold->setEditedState       (pedited->vibrance.psthreshold ? Edited : UnEdited);
		protectSkins->set_inconsistent    (!pedited->vibrance.protectskins);
		avoidColorShift->set_inconsistent (!pedited->vibrance.avoidcolorshift);
		pastSatTog->set_inconsistent   		 (!pedited->vibrance.pastsattog);
	}

	enaconn.block (true);
	enabled->set_active (pp->vibrance.enabled);
	enaconn.block (false);
	lastEnabled = pp->vibrance.enabled;

	pskinsconn.block (true);
	protectSkins->set_active (pp->vibrance.protectskins);
	pskinsconn.block (false);
	lastProtectSkins = pp->vibrance.protectskins;

	ashiftconn.block (true);
	avoidColorShift->set_active (pp->vibrance.avoidcolorshift);
	ashiftconn.block (false);
	lastAvoidColorShift = pp->vibrance.avoidcolorshift;

	pastsattogconn.block (true);
	pastSatTog->set_active (pp->vibrance.pastsattog);
	pastsattogconn.block (false);
	lastPastSatTog = pp->vibrance.pastsattog;

	pastels->setValue (pp->vibrance.pastels);
	psThreshold->setValue<int> (pp->vibrance.psthreshold);

	if (lastPastSatTog) {
		// Link both slider, so we set saturated and psThresholds unsensitive
		psThreshold->set_sensitive(false);
		saturated->set_sensitive(false);
		saturated->setValue (pp->vibrance.pastels);    // Pastels and Saturated are linked
	}
	else {
		// Separate sliders, so we set saturated and psThresholds sensitive again
		psThreshold->set_sensitive(true);
		saturated->set_sensitive(true);
		saturated->setValue (pp->vibrance.saturated);  // Pastels and Saturated are separate
	}

	enableListener ();
}

void Vibrance::write( ProcParams* pp, ParamsEdited* pedited) {
	pp->vibrance.enabled         = enabled->get_active ();
	pp->vibrance.pastels         = pastels->getIntValue();
	pp->vibrance.saturated       = pastSatTog->get_active() ? pp->vibrance.pastels : saturated->getIntValue ();
	pp->vibrance.psthreshold     = psThreshold->getValue<int> ();
	pp->vibrance.protectskins    = protectSkins->get_active ();
	pp->vibrance.avoidcolorshift = avoidColorShift->get_active ();
	pp->vibrance.pastsattog      = pastSatTog->get_active ();

	if (pedited) {
		pedited->vibrance.enabled         = !enabled->get_inconsistent();
		pedited->vibrance.pastels         = pastels->getEditedState ();
		pedited->vibrance.saturated       = saturated->getEditedState ();
		pedited->vibrance.psthreshold     = psThreshold->getEditedState ();
		pedited->vibrance.protectskins    = !protectSkins->get_inconsistent();
		pedited->vibrance.avoidcolorshift = !avoidColorShift->get_inconsistent();
		pedited->vibrance.pastsattog      = !pastSatTog->get_inconsistent();
	}

}

void Vibrance::enabled_toggled () {
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
			listener->panelChanged (EvVibranceEnabled, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvVibranceEnabled, M("GENERAL_DISABLED"));
	}
}

void Vibrance::protectskins_toggled () {
	if (batchMode) {
		if (protectSkins->get_inconsistent()) {
			protectSkins->set_inconsistent (false);
			pskinsconn.block (true);
			protectSkins->set_active (false);
			pskinsconn.block (false);
		}
		else if (lastProtectSkins)
			protectSkins->set_inconsistent (true);

		lastProtectSkins = protectSkins->get_active ();
	}

	if (listener && enabled->get_active()) {
		if (protectSkins->get_active ())
			listener->panelChanged (EvVibranceProtectSkins, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvVibranceProtectSkins, M("GENERAL_DISABLED"));
	}
}

void Vibrance::avoidcolorshift_toggled () {
	if (batchMode) {
		if (avoidColorShift->get_inconsistent()) {
			avoidColorShift->set_inconsistent (false);
			ashiftconn.block (true);
			avoidColorShift->set_active (false);
			ashiftconn.block (false);
		}
		else if (lastAvoidColorShift)
			avoidColorShift->set_inconsistent (true);

		lastAvoidColorShift = avoidColorShift->get_active ();
	}

	if (listener && enabled->get_active()) {
		if (avoidColorShift->get_active ())
			listener->panelChanged (EvVibranceAvoidColorShift, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvVibranceAvoidColorShift, M("GENERAL_DISABLED"));
	}
}

void Vibrance::pastsattog_toggled () {
	if (batchMode) {
		if (pastSatTog->get_inconsistent()) {
			pastSatTog->set_inconsistent (false);
			pastsattogconn.block (true);
			pastSatTog->set_active (false);
			pastsattogconn.block (false);
		}
		else if (lastPastSatTog)
			pastSatTog->set_inconsistent (true);

		lastPastSatTog = pastSatTog->get_active ();
	}

	if (pastSatTog->get_active()) {
		// Link both slider, so we set saturated and psThresholds unsensitive
		psThreshold->set_sensitive(false);
		saturated->set_sensitive(false);
		saturated->setValue (pastels->getValue());     // Pastels and Saturated are linked
	}
	else {
		// Separate sliders, so we set saturated and psThresholds sensitive again
		psThreshold->set_sensitive(true);
		saturated->set_sensitive(true);
	}

	if (listener && enabled->get_active()) {
		if (pastSatTog->get_active ()) {
			listener->panelChanged (EvVibrancePastSatTog, M("GENERAL_ENABLED"));
		}
		else
			listener->panelChanged (EvVibrancePastSatTog, M("GENERAL_DISABLED"));
	}
}

void Vibrance::adjusterChanged (Adjuster* a, double newval) {
	if (a == pastels && pastSatTog->get_active())
		saturated->setValue (newval);

	if (listener && enabled->get_active()) {
		Glib::ustring value = a->getTextValue();

		if (a == pastels )
			listener->panelChanged (EvVibrancePastels,  value );
		else if (a == saturated && !pastSatTog->get_active())
			listener->panelChanged (EvVibranceSaturated, value );
	}
}

void Vibrance::adjusterChanged (ThresholdAdjuster* a, int newBottom, int newTop) {
    if (listener && enabled->get_active()) {
        listener->panelChanged (EvVibrancePastSatThreshold, psThreshold->getHistoryString());
    }
}


void Vibrance::setBatchMode(bool batchMode) {

    ToolPanel::setBatchMode (batchMode);

	pastels->showEditedCB   ();
	saturated->showEditedCB ();
	psThreshold->showEditedCB ();
}

void Vibrance::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited) {
	pastels->setDefault (defParams->vibrance.pastels);
	saturated->setDefault (defParams->vibrance.saturated);
	psThreshold->setDefault<int> (defParams->vibrance.psthreshold);

	if (pedited) {
		pastels->setDefaultEditedState     (pedited->vibrance.pastels ? Edited : UnEdited);
		saturated->setDefaultEditedState   (pedited->vibrance.saturated ? Edited : UnEdited);
		psThreshold->setDefaultEditedState (pedited->vibrance.psthreshold ? Edited : UnEdited);

	} else {
		pastels->setDefaultEditedState     (Irrelevant);
		saturated->setDefaultEditedState   (Irrelevant);
		psThreshold->setDefaultEditedState (Irrelevant);
	}
}

void Vibrance::setAdjusterBehavior (bool pastelsadd, bool saturatedadd, bool psthreshdadd) {
	pastels->setAddMode (pastelsadd);
	saturated->setAddMode (saturatedadd);
}

void Vibrance::trimValues (ProcParams* pp) {
	pastels->trimValue (pp->vibrance.pastels);
	saturated->trimValue (pp->vibrance.saturated);
}
