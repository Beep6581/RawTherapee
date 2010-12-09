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
#include <shcselector.h>
#include <adjuster.h>
#include <curveeditor.h>

/*
 * This class handle the curve widgets, shared between any number curve
 * - to add a curve to the list, use the 'addCurve' method
 * - to start a new line of curve button, use the 'newLine' method
 * - if you add more than one curve, you must add a "CurveEditor* ce" parameter to your listener
 */
class CurveEditorGroup : public Gtk::VBox, public CurveListener, public SHCListener, public AdjusterListener {

private:
	Gtk::Label* curveGroupLabel;

	Gtk::Button* curve_reset;
	Gtk::HBox* customCurveBox;
	Gtk::VBox* paramCurveBox;
	Gtk::HBox* NURBSCurveBox;

	MyCurve* customCurve;
	MyCurve* paramCurve;
	MyCurve* NURBSCurve;

	SHCSelector* shcSelector;
	Adjuster* highlights;
	Adjuster* lights;
	Adjuster* darks;
	Adjuster* shadows;

	Gtk::Button* saveCustom;
	Gtk::Button* loadCustom;
	Gtk::Button* saveNURBS;
	Gtk::Button* loadNURBS;

	CurveListener* cl;

	CurveType curveTypeIx;
	unsigned int numberOfPackedCurve;

	std::vector<CurveEditor*> curveEditors;
	CurveEditor* displayedCurve;

	int activeParamControl;

	void curveResetPressed ();
	void curveTypeToggled ();
	void typeSelectionChanged (CurveEditor* ce, int n);
	void curveTypeToggled (CurveEditor* ce);
	void savePressed ();
	void loadPressed ();
	void hideCurrentCurve ();
	void setUnChanged (bool uc, CurveEditor* ce);
	void storeDisplayedCurve ();
	void restoreDisplayedHistogram();
	void storeCurveValues (CurveEditor* ce, const std::vector<double>& p);
	void typeSelectionChanged (int n);
	void switchGUI();
	void updateGUI (CurveEditor* ce);
	void removeEditor ();
	void curveChanged ();
	void shcChanged ();
	void adjusterChanged (Adjuster* a, double newval);
	bool adjusterEntered (GdkEventCrossing* ev, int ac);
	bool adjusterLeft (GdkEventCrossing* ev, int ac);
	const std::vector<double> getCurveFromGUI (CurveType type);
	void updateBackgroundHistogram (CurveEditor* ce);

public:
	friend class CurveEditor;
	CurveEditorGroup(Glib::ustring groupLabel = "");
	~CurveEditorGroup();
	CurveEditor* addCurve(Glib::ustring curveLabel = "");
	void newLine();
	void curveListComplete();
	void setBatchMode (bool batchMode);
	void setCurveExternal (CurveEditor* ce, const std::vector<double>& c);
	//void on_realize ();
	void setCurveListener (CurveListener* l) { cl = l; }
};

#endif
