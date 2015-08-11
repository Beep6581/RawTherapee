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
#include "curveeditorgroup.h"

class FlatCurveEditor;

class FlatCurveEditorSubGroup: public CurveEditorSubGroup
{

    friend class FlatCurveEditor;

protected:
    Gtk::VBox* CPointsCurveBox;

    MyFlatCurve* CPointsCurve;

    CoordinateAdjuster *CPointsCoordAdjuster;

    Gtk::Button*       saveCPoints;
    Gtk::Button*       loadCPoints;
    Gtk::Button*       copyCPoints;
    Gtk::Button*       pasteCPoints;
    Gtk::ToggleButton* editPointCPoints;
    Gtk::ToggleButton* editCPoints;
    sigc::connection   editCPointsConn, editPointCPointsConn;

public:
    FlatCurveEditorSubGroup(CurveEditorGroup* prt, Glib::ustring& curveDir);
    virtual ~FlatCurveEditorSubGroup();

    FlatCurveEditor* addCurve(Glib::ustring curveLabel = "", bool periodic = true);
    //virtual void updateBackgroundHistogram (CurveEditor* ce);
    void switchGUI();
    void refresh(CurveEditor *curveToRefresh);
    void editModeSwitchedOff();
    void pipetteMouseOver(EditDataProvider *provider, int modifierKey);
    void pipetteButton1Pressed(EditDataProvider *provider, int modifierKey);
    void pipetteButton1Released(EditDataProvider *provider);
    void pipetteDrag(EditDataProvider *provider, int modifierKey);
    void showCoordinateAdjuster(CoordinateProvider *provider);
    void stopNumericalAdjustment();

    bool curveReset (CurveEditor *ce);

protected:
    void storeCurveValues (CurveEditor* ce, const std::vector<double>& p);
    void storeDisplayedCurve ();
    void restoreDisplayedHistogram ();
    void savePressed ();
    void loadPressed ();
    void copyPressed ();
    void pastePressed ();
    void removeEditor ();
    const std::vector<double> getCurveFromGUI (int type);
    void editPointToggled(Gtk::ToggleButton *button);
    void editToggled (Gtk::ToggleButton *button);
};

#endif
