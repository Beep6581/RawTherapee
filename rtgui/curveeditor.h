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
#ifndef _CURVEEDITOR_
#define _CURVEEDITOR_

#include <gtkmm.h>
#include <popuptogglebutton.h>
#include <curveeditorgroup.h>
#include <mycurve.h>

class CurveEditorGroup;

/*
 * This class is an interface between RT and the curve editor group ; it handles the methods
 * related to a specific curve. It is created by CurveEditorGroup::addCurve
 */
class CurveEditor {

private:

	/*
	 * The curve editor contains only one widget (the curve type button) to receive the signals
	 * but it's co-handled by the CurveEditorGroup too
	*/

	// reflects the buttonType active selection ; used as a pre-'selectionChange' reminder value
    CurveType selected;

    PopUpToggleButton* curveType;
    LUTu histogram;	// histogram values
    bool bgHistValid;

	CurveEditorGroup* group;
	std::vector<double> tempCurve;
	std::vector<double> customCurveEd;
	std::vector<double> paramCurveEd;
	std::vector<double> NURBSCurveEd;
	sigc::connection typeconn;

public:

	friend class CurveEditorGroup;
	CurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup);
	~CurveEditor ();
	void typeSelectionChanged (int n);
	void curveTypeToggled();
	void setCurve (const std::vector<double>& p);
	std::vector<double> getCurve ();
	bool isUnChanged ();
	void setUnChanged (bool uc);
	void updateBackgroundHistogram (LUTu & hist);
};


#endif
