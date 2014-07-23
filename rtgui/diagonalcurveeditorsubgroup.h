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
#ifndef _DIAGONALCURVEEDITORSUBGROUP_
#define _DIAGONALCURVEEDITORSUBGROUP_

#include <gtkmm.h>
#include "curveeditorgroup.h"

class DiagonalCurveEditor;

class DiagonalCurveEditorSubGroup : public CurveEditorSubGroup, public SHCListener, public AdjusterListener {

	friend class DiagonalCurveEditor;

protected:
	Gtk::VBox* customCurveBox;
	Gtk::VBox* NURBSCurveBox;
	Gtk::VBox* paramCurveBox;

	MyDiagonalCurve* customCurve;
	MyDiagonalCurve* NURBSCurve;
	MyDiagonalCurve* paramCurve;

	SHCSelector* shcSelector;
	Adjuster* highlights;
	Adjuster* lights;
	Adjuster* darks;
	Adjuster* shadows;

	Adjuster *editedAdjuster;
	int editedAdjusterValue;

	Gtk::Button*       saveCustom;
	Gtk::Button*       loadCustom;
	Gtk::Button*       copyCustom;
	Gtk::Button*       pasteCustom;
	Gtk::ToggleButton* editCustom;
	sigc::connection   editCustomConn;
	Gtk::Button*       saveNURBS;
	Gtk::Button*       loadNURBS;
	Gtk::Button*       copyNURBS;
	Gtk::Button*       pasteNURBS;
	Gtk::ToggleButton* editNURBS;
	sigc::connection   editNURBSConn;
	Gtk::Button*       saveParam;
	Gtk::Button*       loadParam;
	Gtk::Button*       copyParam;
	Gtk::Button*       pasteParam;
	Gtk::ToggleButton* editParam;
	sigc::connection   editParamConn;

	int activeParamControl;

public:
	DiagonalCurveEditorSubGroup(CurveEditorGroup* prt, Glib::ustring& curveDir);
	virtual ~DiagonalCurveEditorSubGroup();

	DiagonalCurveEditor* addCurve(Glib::ustring curveLabel = "");
	virtual void updateBackgroundHistogram (CurveEditor* ce);
	void switchGUI();
	void refresh(CurveEditor *curveToRefresh);
	void editModeSwitchedOff ();
	void pipetteMouseOver(EditDataProvider *provider, int modifierKey);
	void pipetteButton1Pressed(EditDataProvider *provider, int modifierKey);
	void pipetteButton1Released(EditDataProvider *provider);
	void pipetteDrag(EditDataProvider *provider, int modifierKey);

	bool curveReset (CurveEditor *ce);

protected:
	void storeCurveValues (CurveEditor* ce, const std::vector<double>& p);
	void storeDisplayedCurve ();
	void restoreDisplayedHistogram ();
	void savePressed ();
	void loadPressed ();
	void copyPressed ();
	void pastePressed ();
	void editToggled (Gtk::ToggleButton *button);
	void removeEditor ();
	const std::vector<double> getCurveFromGUI (int type);
	void shcChanged ();
	void adjusterChanged (Adjuster* a, double newval);
	bool adjusterEntered (GdkEventCrossing* ev, int ac);
	bool adjusterLeft (GdkEventCrossing* ev, int ac);
	void setSubGroupRangeLabels(Glib::ustring r1, Glib::ustring r2, Glib::ustring r3, Glib::ustring r4);
	void setSubGroupBottomBarBgGradient();
};

#endif
