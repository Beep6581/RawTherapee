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
#include <mycurve.h>
#include <shcselector.h>
#include <adjuster.h>

class CurveEditor : public Gtk::VBox, public CurveListener, public SHCListener, public AdjusterListener {

        Gtk::ComboBoxText* curveType;
        Gtk::Button* curve_reset;
        Gtk::VBox* paramCurveBox;
        Gtk::VBox* customCurveBox;

        MyCurve* customCurve;
        MyCurve* paramCurve;
        SHCSelector* shcSelector;
        
        Adjuster* highlights;
        Adjuster* lights;
        Adjuster* darks;
        Adjuster* shadows;
        
        Gtk::Button* save;
        Gtk::Button* load;
        
        CurveListener* cl;
        
        bool realized;
        std::vector<double> tmpCurve;
        int curveTypeIx;
        
        int activeParamControl;
        
        sigc::connection typeconn;

    public:

        CurveEditor ();
        virtual ~CurveEditor ();
        void setBatchMode (bool batchMode);
        bool isUnChanged ();
        void setUnChanged (bool uc);
        
        void on_realize ();
        void setCurveListener (CurveListener* l) { cl = l; }
        void savePressed ();
        void loadPressed ();
        void typeSelectionChanged ();
        void setCurve (const std::vector<double>& c);
        std::vector<double> getCurve ();
        void curveChanged ();
        void curveResetPressed ();
        void shcChanged ();
        void adjusterChanged (Adjuster* a, double newval);
        bool adjusterEntered (GdkEventCrossing* ev, int ac);
        bool adjusterLeft (GdkEventCrossing* ev, int ac);
        void updateBackgroundHistogram (unsigned int* hist);
};


#endif
