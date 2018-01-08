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

#include "popuptogglebutton.h"
#include "../rtengine/LUT.h"
#include "coloredbar.h"
#include "edit.h"
#include "mydiagonalcurve.h"
#include "myflatcurve.h"

class CurveEditorGroup;
class CurveEditorSubGroup;


/*
 *********************** Curve Editor ***********************
 */


/** @brief This class is an interface between RT and the curve editor group
 * It handles the methods related to a specific curve. It is created by CurveEditorGroup::addCurve
 */
class CurveEditor : public rtedit::EditSubscriber
{

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
    LUTu histogram; // histogram values
    bool bgHistValid;

    bool remoteDrag;

    int selected;

    CurveEditorGroup* group;
    CurveEditorSubGroup* subGroup;
    Gtk::Widget* relatedWidget;
    bool expandRelatedWidget;

    std::vector<double> tempCurve;
    sigc::connection typeconn;

    ColorProvider* bottomBarCP;
    ColorProvider* leftBarCP;
    ColorProvider* curveCP;
    int bottomBarCId;
    int leftBarCId;
    int curveCId;
    std::vector<GradientMilestone> bottomBarBgGradient;
    std::vector<GradientMilestone> leftBarBgGradient;

    sigc::signal<void> sig_curvegraph_enter;
    sigc::signal<void> sig_curvegraph_leave;
    sigc::signal<void> sig_curvepoint_click;
    sigc::signal<void> sig_curvepoint_release;

public:

    CurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup);
    virtual ~CurveEditor ();
    void typeSelectionChanged (int n);
    void curveTypeToggled();
    bool isUnChanged ();
    void setUnChanged (bool uc);
    void updateBackgroundHistogram (LUTu & hist);

    void setLeftBarColorProvider(ColorProvider* cp, int callerId);
    void setBottomBarColorProvider(ColorProvider* cp, int callerId);
    void setCurveColorProvider(ColorProvider* cp, int callerId);
    void setBottomBarBgGradient (const std::vector<GradientMilestone> &milestones);
    void setLeftBarBgGradient (const std::vector<GradientMilestone> &milestones);
    ColorProvider* getLeftBarColorProvider();
    ColorProvider* getBottomBarColorProvider();
    ColorProvider* getCurveColorProvider();
    int getLeftBarCallerId();
    int getBottomBarCallerId();
    int getCurveCallerId();
    std::vector<GradientMilestone> getBottomBarBgGradient () const;
    std::vector<GradientMilestone> getLeftBarBgGradient () const;

    void refresh (); // refresh the display of the CurveEditor (e.g. when a ColoredBar has been changed from the outside)
    bool openIfNonlinear();  // Open up the curve if it has modifications and it's not already opened

    void setCurve (const std::vector<double>& p);
    virtual void setIdentityValue (const double iValue = 0.5) {};
    virtual double getIdentityValue ()
    {
        return 0.5;
    };
    virtual std::vector<double> getCurve () = 0;
    bool reset();
    void setTooltip(Glib::ustring ttip);

    sigc::signal<void> signal_curvegraph_enter();
    sigc::signal<void> signal_curvegraph_leave();
    sigc::signal<void> signal_curvepoint_click();
    sigc::signal<void> signal_curvepoint_release();

    void switchOffEditMode ();
    bool mouseOver(const int modifierKey);
    bool button1Pressed(const int modifierKey);
    bool button1Released();
    bool drag1(const int modifierKey);
    CursorShape getCursor(const int objectID);


};


/*
 ******************** Diagonal Curve Editor ********************
 */


class DiagonalCurveEditor : public CurveEditor
{

    friend class DiagonalCurveEditorSubGroup;

protected:
    // reflects the buttonType active selection ; used as a pre-'selectionChange' reminder value
    std::vector<double> customCurveEd;
    std::vector<double> customResetCurve;
    std::vector<double> paramCurveEd;
    std::vector<double> paramResetCurve;
    std::vector<double> NURBSCurveEd;
    std::vector<double> NURBSResetCurve;
    Glib::ustring rangeLabels[4];
    double rangeMilestones[3];

public:
    DiagonalCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup);
    std::vector<double> getCurve ();
    void setRangeLabels(Glib::ustring r1, Glib::ustring r2, Glib::ustring r3, Glib::ustring r4);
    void getRangeLabels(Glib::ustring &r1, Glib::ustring &r2, Glib::ustring &r3, Glib::ustring &r4);
    void setRangeDefaultMilestones(double m1, double m2, double m3);
    void getRangeDefaultMilestones(double &m1, double &m2, double &m3);

    // set the reset curve for a given curve type. This is optional; all curve type have a default reset curve
    void setResetCurve(DiagonalCurveType cType, const std::vector<double> &resetCurve);
};


/*
 ********************** Flat Curve Editor **********************
 */


class FlatCurveEditor : public CurveEditor
{

    friend class FlatCurveEditorSubGroup;

protected:
    // reflects the buttonType active selection ; used as a pre-'selectionChange' reminder value
    std::vector<double> controlPointsCurveEd;
    std::vector<double> controlPointsResetCurve;
    bool periodic;
    double identityValue;

public:
    FlatCurveEditor (Glib::ustring text, CurveEditorGroup* ceGroup, CurveEditorSubGroup* ceSubGroup, bool isPeriodic = true);
    virtual void setIdentityValue (const double iValue = 0.5)
    {
        identityValue = iValue;
    }
    virtual double getIdentityValue ()
    {
        return identityValue;
    };
    std::vector<double> getCurve ();

    // set the reset curve for a given curve type. This is optional; all curve type have a default reset curve
    void setResetCurve(FlatCurveType cType, const std::vector<double> &resetCurve);
};

#endif
