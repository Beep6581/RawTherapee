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
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include <fstream>
#include <string>
#include "guiutils.h"
#include "multilangmgr.h"
#include "../rtengine/LUT.h"

#include <cstring>

DiagonalCurveEditor::DiagonalCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup) : CurveEditor::CurveEditor(text, static_cast<CurveEditorGroup*>(ceGroup), ceSubGroup) {

    // Order set in the same order than "enum DiagonalCurveType". Shouldn't change, for compatibility reason
    curveType->addEntry("curveType-linear.png",     M("CURVEEDITOR_LINEAR"));		// 0 Linear
    curveType->addEntry("curveType-spline.png",     M("CURVEEDITOR_CUSTOM"));		// 1 Spline
    curveType->addEntry("curveType-parametric.png", M("CURVEEDITOR_PARAMETRIC"));	// 2 Parametric
    curveType->addEntry("curveType-NURBS.png",      M("CURVEEDITOR_NURBS"));		// 3 NURBS
    curveType->setSelected(DCT_Linear);
    curveType->show();

    rangeLabels[0] = M("CURVEEDITOR_SHADOWS");
    rangeLabels[1] = M("CURVEEDITOR_DARKS");
    rangeLabels[2] = M("CURVEEDITOR_LIGHTS");
    rangeLabels[3] = M("CURVEEDITOR_HIGHLIGHTS");

    rangeMilestones[0] = 0.25;
    rangeMilestones[1] = 0.50;
    rangeMilestones[2] = 0.75;
}

std::vector<double> DiagonalCurveEditor::getCurve () {
	std::vector<double> curve;

	switch (selected) {
	case (DCT_Spline):
        return curve = customCurveEd;
	case (DCT_Parametric):
        return curve = paramCurveEd;
	case (DCT_NURBS):
        return curve = NURBSCurveEd;
	default:
		// returning Linear or Unchanged
		curve.push_back((double)(selected));
		return curve;
	}
}

void DiagonalCurveEditor::setRangeLabels(Glib::ustring r1, Glib::ustring r2, Glib::ustring r3, Glib::ustring r4) {
    rangeLabels[0] = r1;
    rangeLabels[1] = r2;
    rangeLabels[2] = r3;
    rangeLabels[3] = r4;
}

void DiagonalCurveEditor::getRangeLabels(Glib::ustring &r1, Glib::ustring &r2, Glib::ustring &r3, Glib::ustring &r4) {
    r1 = rangeLabels[0];
    r2 = rangeLabels[1];
    r3 = rangeLabels[2];
    r4 = rangeLabels[3];
}

/*
 * Admittedly that this method is called just after the instantiation of this class, we set the shcselector's default values
 */
void DiagonalCurveEditor::setRangeDefaultMilestones(double m1, double m2, double m3) {
    rangeMilestones[0] = m1;
    rangeMilestones[1] = m2;
    rangeMilestones[2] = m3;

    paramCurveEd.at(1) = m1;
    paramCurveEd.at(2) = m2;
    paramCurveEd.at(3) = m3;
}

void DiagonalCurveEditor::getRangeDefaultMilestones(double &m1, double &m2, double &m3) {
    m1 = rangeMilestones[0];
    m2 = rangeMilestones[1];
    m3 = rangeMilestones[2];
}

FlatCurveEditor::FlatCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup, bool isPeriodic) : CurveEditor::CurveEditor(text, static_cast<CurveEditorGroup*>(ceGroup), ceSubGroup) {

	periodic = isPeriodic;

    // Order set in the same order than "enum FlatCurveType". Shouldn't change, for compatibility reason
    curveType->addEntry("curveType-flatLinear.png",    M("CURVEEDITOR_LINEAR"));			// 0 Linear
    curveType->addEntry("curveType-controlPoints.png", M("CURVEEDITOR_MINMAXCPOINTS"));		// 1 Min/Max ControlPoints
    curveType->setSelected(FCT_Linear);
    curveType->show();
}

std::vector<double> FlatCurveEditor::getCurve () {
	std::vector<double> curve;

	switch (selected) {
	//case (Parametric):
    //    return curve = paramCurveEd;
	case (FCT_MinMaxCPoints):
        return curve = controlPointsCurveEd;
	default:
		// returning Linear or Unchanged
		curve.push_back((double)(selected));
		return curve;
	}
}

/*
 * CurveEditor (CurveEditorGroup* ceGroup, Glib::ustring text)
 *
 * parameters:
 * 		ceGroup = NULL or the address of the Widget that will receive the CurveTypeToggleButton
 * 		text    = (optional) label of the curve, displayed in the CurveTypeToggleButton, next to the image
 */
CurveEditor::CurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup) {

	bgHistValid = false;
	selected = DCT_Linear;
	bottomBarCP = NULL;
	leftBarCP = NULL;
	curveCP = NULL;

	group = ceGroup;
	subGroup = ceSubGroup;

    if (group && text.size())
    	curveType = new PopUpToggleButton(text + ":");
    else
    	curveType = new PopUpToggleButton();

    curveType->set_tooltip_text(M("CURVEEDITOR_TYPE"));
    // TODO: Does this signal have to be blocked when on curve type change ?
    curveType->signal_toggled().connect ( sigc::mem_fun(*this, &CurveEditor::curveTypeToggled) );
	typeconn  = curveType->signal_changed().connect (sigc::mem_fun(*this, &CurveEditor::typeSelectionChanged) );
}

void CurveEditor::setCurve (const std::vector<double>& p) {
	tempCurve = p;
	group->setCurveExternal(this, p);
}

CurveEditor::~CurveEditor () {
    delete curveType;
}

void CurveEditor::typeSelectionChanged (int n) {
	group->typeSelectionChanged(this, n);
}

void CurveEditor::curveTypeToggled() {
	group->curveTypeToggled(this);
}

bool CurveEditor::isUnChanged () {
    return curveType->getSelected()==subGroup->getValUnchanged();
}

void CurveEditor::setUnChanged (bool uc) {
	group->setUnChanged(uc, this);
}

/*
 * Update the backgrounds histograms
 */
void CurveEditor::updateBackgroundHistogram (LUTu & hist) {
	// Copy the histogram in the curve editor cache
	if (hist) {
		histogram=hist;
		bgHistValid = true;
	}
	else
		bgHistValid = false;
	
	// Then call the curve editor group to eventually update the histogram
	subGroup->updateBackgroundHistogram (this);
}

// Open up the curve if it has modifications and it's not already opened
// Returns: is it non-linear?
bool CurveEditor::openIfNonlinear() {
    bool nonLinear = tempCurve.size() && (tempCurve[0] > subGroup->getValLinear()) && (tempCurve[0] < subGroup->getValUnchanged());

    if (nonLinear && !curveType->get_active()) {
        // Will toggle the event doing the display
        curveType->set_active( true );
    }

    return nonLinear;
}

// Handles markup tooltips
void CurveEditor::setTooltip(Glib::ustring ttip) {
	curveType->set_tooltip_text(ttip.empty() ?
			Glib::ustring::compose("<b>%1</b> ", M("CURVEEDITOR_TYPE")) :
			Glib::ustring::compose("%1\n<b>%2</b>", ttip, M("CURVEEDITOR_TYPE")));
}

void CurveEditor::setLeftBarColorProvider(ColorProvider* cp) {
	leftBarCP = cp;
}

void CurveEditor::setBottomBarColorProvider(ColorProvider* cp) {
	bottomBarCP = cp;
}

void CurveEditor::setLeftBarBgGradient (const std::vector<GradientMilestone> &milestones) {
	leftBarBgGradient = milestones;
}

void CurveEditor::setBottomBarBgGradient (const std::vector<GradientMilestone> &milestones) {
	bottomBarBgGradient = milestones;
}

void CurveEditor::setCurveColorProvider(ColorProvider* cp) {
	curveCP = cp;
}

ColorProvider* CurveEditor::getLeftBarColorProvider() {
	return leftBarCP;
}

ColorProvider* CurveEditor::getBottomBarColorProvider() {
	return bottomBarCP;
}

ColorProvider* CurveEditor::getCurveColorProvider() {
	return curveCP;
}

std::vector<GradientMilestone> CurveEditor::getBottomBarBgGradient () const {
	return bottomBarBgGradient;
}

std::vector<GradientMilestone> CurveEditor::getLeftBarBgGradient () const {
	return leftBarBgGradient;
}

sigc::signal<void> CurveEditor::signal_curvegraph_enter() {
	return sig_curvegraph_enter;
}

sigc::signal<void> CurveEditor::signal_curvegraph_leave() {
	return sig_curvegraph_leave;
}

sigc::signal<void> CurveEditor::signal_curvepoint_click() {
	return sig_curvepoint_click;
}

sigc::signal<void> CurveEditor::signal_curvepoint_release() {
	return sig_curvepoint_release;
}
