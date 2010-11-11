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

/*
 * CurveEditor (CurveEditorGroup* ceGroup, Glib::ustring text)
 *
 * parameters:
 * 		ceGroup = NULL or the address of the Widget that will receive the CurveTypeToggleButton
 * 		text    = (optional) label of the curve, displayed in the CurveTypeToggleButton, next to the image
 */
CurveEditor::CurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup) {

	bgHistValid = false;
	selected = Linear;

	group = ceGroup;

    if (group && text.size())
    	curveType = Gtk::manage (new PopUpToggleButton(text + ":", true));
    else
    	curveType = Gtk::manage (new PopUpToggleButton());

    curveType->set_image_position(Gtk::POS_RIGHT);
    // Order set in the same order than "enum CurveType". Shouldn't change, for compatibility reason
    curveType->addEntry(argv0+"/images/curveType-linear.png", M("CURVEEDITOR_LINEAR"));			// 0 Linear
    curveType->addEntry(argv0+"/images/curveType-spline.png", M("CURVEEDITOR_CUSTOM"));			// 1 Spline
    curveType->addEntry(argv0+"/images/curveType-parametric.png", M("CURVEEDITOR_PARAMETRIC"));	// 2 Parametric
    curveType->addEntry(argv0+"/images/curveType-NURBS.png", M("CURVEEDITOR_NURBS"));			// 3 NURBS
    curveType->setSelected(Linear);
    curveType->set_tooltip_text(M("CURVEEDITOR_TYPE"));
    // TODO: Does this signal have to be blocked when on curve type change ?
    curveType->signal_toggled().connect ( sigc::mem_fun(*this, &CurveEditor::curveTypeToggled) );
	typeconn  = curveType->signal_changed().connect (sigc::mem_fun(*this, &CurveEditor::typeSelectionChanged) );

    curveType->show();
}

/*
CurveEditor::~CurveEditor () {

}
*/

void CurveEditor::setCurve (const std::vector<double>& p) {
	tempCurve = p;
	group->setCurveExternal(this, p);
}

std::vector<double> CurveEditor::getCurve () {
	std::vector<double> curve;

	switch (selected) {
	case (Spline):
        return curve = customCurveEd;
	case (Parametric):
        return curve = paramCurveEd;
	case (NURBS):
        return curve = NURBSCurveEd;
	default:
		// returning Linear or Unchanged
		curve.push_back((double)(selected));
		return curve;
	}
}

void CurveEditor::typeSelectionChanged (int n) {
	group->typeSelectionChanged(this, n);
}

void CurveEditor::curveTypeToggled() {
	group->curveTypeToggled(this);
}

bool CurveEditor::isUnChanged () {
    return curveType->getSelected()==Unchanged;
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
	group->updateBackgroundHistogram (this);
}
