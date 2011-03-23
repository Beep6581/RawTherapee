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
#ifndef _FLATCURVEEDITORSUBGROUP_
#define _FLATCURVEEDITORSUBGROUP_

#include <gtkmm.h>
#include <curveeditorgroup.h>

class FlatCurveEditor;

class FlatCurveEditorSubGroup: public CurveEditorSubGroup  {

	friend class FlatCurveEditor;

private:
	Gtk::HBox* CPointsCurveBox;

	MyFlatCurve* CPointsCurve;

	Gtk::Button* saveCPoints;
	Gtk::Button* loadCPoints;

public:
	FlatCurveEditorSubGroup(CurveEditorGroup* prt);
	~FlatCurveEditorSubGroup();

	FlatCurveEditor* addCurve(Glib::ustring curveLabel = "");
	//virtual void updateBackgroundHistogram (CurveEditor* ce);
	virtual void setColorProvider (ColorProvider* p);

private:
	void storeCurveValues (CurveEditor* ce, const std::vector<double>& p);
	void storeDisplayedCurve ();
	void restoreDisplayedHistogram ();
	void savePressed ();
	void loadPressed ();
	void switchGUI();
	bool curveReset (int cType);
	void removeEditor ();
	const std::vector<double> getCurveFromGUI (int type);
};

#endif
