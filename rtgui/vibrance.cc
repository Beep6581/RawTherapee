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
#include "../rtengine/color.h"
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

Vibrance::Vibrance () : FoldableToolPanel(this, "vibrance", M("TP_VIBRANCE_LABEL"), false, true) {

	std::vector<GradientMilestone> milestones;
	float R, G, B;
	// -0.1 rad < Hue < 1.6 rad
	Color::hsv2rgb01(0.92f, 0.45f, 0.6f, R, G, B);
	milestones.push_back( GradientMilestone(0.0, double(R), double(G), double(B)) );
	Color::hsv2rgb01(0.14056f, 0.45f, 0.6f, R, G, B);
	milestones.push_back( GradientMilestone(1.0, double(R), double(G), double(B)) );

	saturated = Gtk::manage(new Adjuster (M("TP_VIBRANCE_SATURATED"),-100.,100.,1.,0.));
	saturated->setAdjusterListener (this);
	saturated->set_sensitive(false);
	//if (saturated->delay < 1000) saturated->delay = 1000;
	pack_start( *saturated, Gtk::PACK_SHRINK, 0);

	pastels = Gtk::manage(new Adjuster (M("TP_VIBRANCE_PASTELS"),-100.,100.,1.,0.));
	pastels->setAdjusterListener (this);
	//if (pastels->delay < 1000) pastels->delay = 1000;
	pack_start( *pastels, Gtk::PACK_SHRINK, 0);

	psThreshold = Gtk::manage (new ThresholdAdjuster (M("TP_VIBRANCE_PSTHRESHOLD"), -100., 100., 0., M("TP_VIBRANCE_PSTHRESHOLD_WEIGTHING"), 0, 0., 100., 75., M("TP_VIBRANCE_PSTHRESHOLD_SATTHRESH"), 0, this, false));
	psThreshold->setAdjusterListener (this);
	psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
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

	curveEditorGG = new CurveEditorGroup (options.lastVibranceCurvesDir, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_LABEL"));
	curveEditorGG->setCurveListener (this);

	skinTonesCurve = static_cast<DiagonalCurveEditor*>(curveEditorGG->addCurve(CT_Diagonal, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES")));
	skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
	skinTonesCurve->setBottomBarBgGradient(milestones);
	skinTonesCurve->setLeftBarBgGradient(milestones);
	skinTonesCurve->setRangeLabels(
			M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE1"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE2"),
			M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE3"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE4")
	);
	skinTonesCurve->setRangeDefaultMilestones(0.1, 0.4, 0.85);
	curveEditorGG->curveListComplete();

	pack_start (*curveEditorGG, Gtk::PACK_SHRINK, 4);

	show ();

	pskinsconn = protectSkins->signal_toggled().connect( sigc::mem_fun(*this, &Vibrance::protectskins_toggled) );
	ashiftconn = avoidColorShift->signal_toggled().connect( sigc::mem_fun(*this, &Vibrance::avoidcolorshift_toggled) );
	pastsattogconn = pastSatTog->signal_toggled().connect( sigc::mem_fun(*this, &Vibrance::pastsattog_toggled) );
}

Vibrance::~Vibrance () {
	delete curveEditorGG;
}

void Vibrance::read(const ProcParams* pp, const ParamsEdited* pedited) {
	disableListener ();

	if(pedited ){
		set_inconsistent                  (multiImage && !pedited->vibrance.enabled);
		pastels->setEditedState           (pedited->vibrance.pastels ? Edited : UnEdited);
		saturated->setEditedState         (pedited->vibrance.saturated ? Edited : UnEdited);
		psThreshold->setEditedState       (pedited->vibrance.psthreshold ? Edited : UnEdited);
		protectSkins->set_inconsistent    (!pedited->vibrance.protectskins);
		avoidColorShift->set_inconsistent (!pedited->vibrance.avoidcolorshift);
		pastSatTog->set_inconsistent      (!pedited->vibrance.pastsattog);
		skinTonesCurve->setUnChanged      (!pedited->vibrance.skintonescurve);
	}

	setEnabled (pp->vibrance.enabled);

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
	skinTonesCurve->setCurve (pp->vibrance.skintonescurve);

	enableListener ();
}

void Vibrance::autoOpenCurve  () {
	skinTonesCurve->openIfNonlinear();
}

void Vibrance::write( ProcParams* pp, ParamsEdited* pedited) {
	pp->vibrance.enabled         = getEnabled ();
	pp->vibrance.pastels         = pastels->getIntValue();
	pp->vibrance.saturated       = pastSatTog->get_active() ? pp->vibrance.pastels : saturated->getIntValue ();
	pp->vibrance.psthreshold     = psThreshold->getValue<int> ();
	pp->vibrance.protectskins    = protectSkins->get_active ();
	pp->vibrance.avoidcolorshift = avoidColorShift->get_active ();
	pp->vibrance.pastsattog      = pastSatTog->get_active ();
	pp->vibrance.skintonescurve  = skinTonesCurve->getCurve ();

	if (pedited) {
		pedited->vibrance.enabled         = !get_inconsistent();
		pedited->vibrance.pastels         = pastels->getEditedState ();
		pedited->vibrance.saturated       = saturated->getEditedState ();
		pedited->vibrance.psthreshold     = psThreshold->getEditedState ();
		pedited->vibrance.protectskins    = !protectSkins->get_inconsistent();
		pedited->vibrance.avoidcolorshift = !avoidColorShift->get_inconsistent();
		pedited->vibrance.pastsattog      = !pastSatTog->get_inconsistent();
		pedited->vibrance.skintonescurve  = !skinTonesCurve->isUnChanged ();
	}

}
void Vibrance::curveChanged () {

	if (listener && getEnabled ()) listener->panelChanged (EvVibranceSkinTonesCurve, M("HISTORY_CUSTOMCURVE"));
}

void Vibrance::enabledChanged () {
	if (listener) {
		if (get_inconsistent())
			listener->panelChanged (EvVibranceEnabled, M("GENERAL_UNCHANGED"));
		if (getEnabled())
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

	if (listener && getEnabled()) {
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

	if (listener && getEnabled()) {
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

	if (listener && getEnabled()) {
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

	if (listener && getEnabled()) {
		Glib::ustring value = a->getTextValue();

		if (a == pastels )
			listener->panelChanged (EvVibrancePastels,  value );
		else if (a == saturated && !pastSatTog->get_active())
			listener->panelChanged (EvVibranceSaturated, value );
	}
}

void Vibrance::adjusterChanged (ThresholdAdjuster* a, int newBottom, int newTop) {
	if (listener && getEnabled()) {
		listener->panelChanged (EvVibrancePastSatThreshold, psThreshold->getHistoryString());
	}
}


void Vibrance::setBatchMode(bool batchMode) {

    ToolPanel::setBatchMode (batchMode);

	pastels->showEditedCB   ();
	saturated->showEditedCB ();
	psThreshold->showEditedCB ();

	curveEditorGG->setBatchMode (batchMode);
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

void Vibrance::setAdjusterBehavior (bool pastelsadd, bool saturatedadd) {
	pastels->setAddMode (pastelsadd);
	saturated->setAddMode (saturatedadd);
}

void Vibrance::trimValues (ProcParams* pp) {
	pastels->trimValue (pp->vibrance.pastels);
	saturated->trimValue (pp->vibrance.saturated);
}

std::vector<double> Vibrance::getCurvePoints(ThresholdSelector* tAdjuster) const {
	std::vector<double> points;
	double threshold, transitionWeighting;
	tAdjuster->getPositions<double>(transitionWeighting, threshold);  // ( range -100;+100,   range 0;+100 )
	transitionWeighting /= 100.; // range -1., +1.
	threshold /= 100.;      // range  0., +1.

	// Initial point
	points.push_back(0.); points.push_back(0.);

	double p2 = 3.0*threshold/4.0;                 // same one than in ipvibrance.cc
	double s0 = threshold + (1.0-threshold)/4.0;   // same one than in ipvibrance.cc

	// point at the beginning of the first linear transition
	points.push_back(p2); points.push_back(0.);

	// Y value of the chroma mean point, calculated to get a straight line between p2 and s0
	double chromaMean = (threshold/4.0) / (s0-p2);
	// move chromaMean up or down depending on transitionWeighting
	if (transitionWeighting > 0.0) {
		// positive values -> give more weight to Saturated
		chromaMean = (1.0-chromaMean) * transitionWeighting + chromaMean;
	}
	else if (transitionWeighting < 0.0) {
		// negative values -> give more weight to Pastels
		chromaMean =      chromaMean  * transitionWeighting + chromaMean;
	}

	// point at the location of the Top cursor, at the end of the first linear transition and the beginning of the second one
	points.push_back(threshold); points.push_back(chromaMean);

	if (threshold < 1.0) {

		// point at the end of the second linear transition
		points.push_back(s0); points.push_back(1.0);

		// end point
		points.push_back(1.0); points.push_back(1.0);
	}

	return points;
}
