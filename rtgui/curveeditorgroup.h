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
#ifndef _CURVEEDITORGROUP_
#define _CURVEEDITORGROUP_

#include <gtkmm.h>
#include <fstream>
#include <string>
#include <guiutils.h>
#include <mycurve.h>
#include <myflatcurve.h>
#include <mydiagonalcurve.h>
#include <shcselector.h>
#include <adjuster.h>

class CurveEditor;
class DiagonalCurveEditorSubGroup;
class FlatCurveEditorSubGroup;

/*
 * This class handle the curve widgets, shared between any number curve
 * - to add a curve to the list, use the 'addCurve' method
 * - to start a new line of curve button, use the 'newLine' method
 * - if you add more than one curve, you must add a "CurveEditor* ce" parameter to your listener
 */
class CurveEditorGroup : public Gtk::VBox, public CurveListener {

	friend class CurveEditor;
	friend class CurveEditorSubGroup;
	friend class DiagonalCurveEditorSubGroup;
	friend class FlatCurveEditorSubGroup;

protected:
	Gtk::Label* curveGroupLabel;
	Gtk::Button* curve_reset;
	std::vector<CurveEditor*> curveEditors;
	CurveEditor* displayedCurve;
	FlatCurveEditorSubGroup* flatSubGroup;
	DiagonalCurveEditorSubGroup* diagonalSubGroup;

	CurveListener* cl;
	ColorProvider* cp;

	unsigned int numberOfPackedCurve;

public:
	CurveEditorGroup(Glib::ustring groupLabel = "");
	~CurveEditorGroup();
	void newLine();
	void curveListComplete();
	void setBatchMode (bool batchMode);
	void setCurveExternal (CurveEditor* ce, const std::vector<double>& c);
	void setCurveListener (CurveListener* l) { cl = l; }
	void setColorProvider (ColorProvider* p) { cp = p; }
	CurveEditor* getDisplayedCurve () { return displayedCurve; }
	//void on_realize ();
	CurveEditor* addCurve(CurveType cType, Glib::ustring curveLabel);

protected:
	//void curveTypeToggled ();
	void curveTypeToggled (CurveEditor* ce);
	//void typeSelectionChanged (int n);
	void typeSelectionChanged (CurveEditor* ce, int n);
	void hideCurrentCurve ();
	void updateGUI (CurveEditor* ce);
	void curveResetPressed ();
	void curveChanged ();
	void setUnChanged (bool uc, CurveEditor* ce);
};

class CurveEditorSubGroup {

	friend class CurveEditorGroup;

protected:
	int valLinear;
	int valUnchanged;
	CurveEditorGroup *parent;

public:
	int getValUnchanged() { return valUnchanged; }
	virtual void updateBackgroundHistogram (CurveEditor* ce) {}
	virtual void setColorProvider (ColorProvider* p) = 0;

protected:
	Glib::ustring outputFile ();
	Glib::ustring inputFile ();

	virtual bool curveReset (int cType) = 0; // Reset a curve editor, return TRUE if successful (curve changed)
	virtual void storeCurveValues (CurveEditor* ce, const std::vector<double>& p) = 0;
	virtual void storeDisplayedCurve () = 0;
	virtual void restoreDisplayedHistogram() {};
	virtual void switchGUI() = 0;
	virtual void removeEditor () = 0;
	virtual const std::vector<double> getCurveFromGUI (int type) = 0;

};

#endif
