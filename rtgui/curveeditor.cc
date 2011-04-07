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
#include <curveeditor.h>
#include <curveeditorgroup.h>
#include <fstream>
#include <string>
#include <guiutils.h>
#include <multilangmgr.h>

extern Glib::ustring argv0;

DiagonalCurveEditor::DiagonalCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup) : CurveEditor::CurveEditor(text, (CurveEditorGroup*) ceGroup, ceSubGroup) {

    // Order set in the same order than "enum DiagonalCurveType". Shouldn't change, for compatibility reason
    curveType->addEntry(argv0+"/images/curveType-linear.png", M("CURVEEDITOR_LINEAR"));			// 0 Linear
    curveType->addEntry(argv0+"/images/curveType-spline.png", M("CURVEEDITOR_CUSTOM"));			// 1 Spline
    curveType->addEntry(argv0+"/images/curveType-parametric.png", M("CURVEEDITOR_PARAMETRIC"));	// 2 Parametric
    curveType->addEntry(argv0+"/images/curveType-NURBS.png", M("CURVEEDITOR_NURBS"));			// 3 NURBS
    curveType->setSelected(DCT_Linear);
    curveType->show();
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

FlatCurveEditor::FlatCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup) : CurveEditor::CurveEditor(text, (CurveEditorGroup*) ceGroup, ceSubGroup) {

    // Order set in the same order than "enum FlatCurveType". Shouldn't change, for compatibility reason
    curveType->addEntry(argv0+"/images/curveType-flatLinear.png", M("CURVEEDITOR_LINEAR"));				// 0 Linear
    curveType->addEntry(argv0+"/images/curveType-controlPoints.png", M("CURVEEDITOR_MINMAXCPOINTS"));	// 1 Min/Max ControlPoints
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

    histogram = new unsigned int[256];	// histogram values

	group = ceGroup;
	subGroup = ceSubGroup;

    if (group && text.size())
    	curveType = Gtk::manage (new PopUpToggleButton(text + ":"));
    else
    	curveType = Gtk::manage (new PopUpToggleButton());

    curveType->set_tooltip_text(M("CURVEEDITOR_TYPE"));
    // TODO: Does this signal have to be blocked when on curve type change ?
    curveType->signal_toggled().connect ( sigc::mem_fun(*this, &CurveEditor::curveTypeToggled) );
	typeconn  = curveType->signal_changed().connect (sigc::mem_fun(*this, &CurveEditor::typeSelectionChanged) );
}

CurveEditor::~CurveEditor () {

	delete [] histogram;
}

void CurveEditor::setCurve (const std::vector<double>& p) {
	tempCurve = p;
	group->setCurveExternal(this, p);
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
void CurveEditor::updateBackgroundHistogram (unsigned int* hist) {
	// Copy the histogram in the curve editor cache
	if (hist!=NULL) {
		memcpy (histogram, hist, 256*sizeof(unsigned int));
		bgHistValid = true;
	}
	else
		bgHistValid = false;

	// Then call the curve editor group to eventually update the histogram
	subGroup->updateBackgroundHistogram (this);
}
