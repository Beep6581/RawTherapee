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

#include <popuptogglebutton.h>
#include <LUT.h>

class CurveEditorGroup;
class CurveEditorSubGroup;


/*
 *********************** Curve Editor ***********************
 */


/*
 * This class is an interface between RT and the curve editor group ; it handles the methods
 * related to a specific curve. It is created by CurveEditorGroup::addCurve
 */
class CurveEditor {

	friend class CurveEditorGroup;
	friend class CurveEditorSubGroup;
	friend class DiagonalCurveEditorSubGroup;
	friend class FlatCurveEditorSubGroup;
	friend class DiagonalCurveEditor;
	friend class FlatCurveEditor;

	protected:

		/*
		 * The curve editor contains only one widget (the curve type button) to receive the signals
		 * but it's co-handled by the CurveEditorGroup too
		*/

		PopUpToggleButton* curveType;
		LUTu histogram;	// histogram values
		bool bgHistValid;

		int selected;

		CurveEditorGroup* group;
		CurveEditorSubGroup* subGroup;

		std::vector<double> tempCurve;
		sigc::connection typeconn;

	public:

		CurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup);
		~CurveEditor ();
		void typeSelectionChanged (int n);
		void curveTypeToggled();
		bool isUnChanged ();
		void setUnChanged (bool uc);
		void updateBackgroundHistogram (LUTu & hist);
		void setCurve (const std::vector<double>& p);
		virtual std::vector<double> getCurve () = 0;
};


/*
 ******************** Diagonal Curve Editor ********************
 */


class DiagonalCurveEditor : public CurveEditor {

	friend class DiagonalCurveEditorSubGroup;

	protected:
		// reflects the buttonType active selection ; used as a pre-'selectionChange' reminder value
		std::vector<double> customCurveEd;
		std::vector<double> paramCurveEd;
		std::vector<double> NURBSCurveEd;

	public:
		DiagonalCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup);
		std::vector<double> getCurve ();
};


/*
 ********************** Flat Curve Editor **********************
 */


class FlatCurveEditor : public CurveEditor {

	friend class FlatCurveEditorSubGroup;

	protected:
		// reflects the buttonType active selection ; used as a pre-'selectionChange' reminder value
		std::vector<double> controlPointsCurveEd;

	public:
		FlatCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup);
		std::vector<double> getCurve ();
};

#endif
