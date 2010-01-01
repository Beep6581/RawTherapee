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

class CurveEditor : public Gtk::VBox {

        MyCurve* curve;
        Gtk::Button* linear;
        Gtk::Button* save;
        Gtk::Button* load;

    public:
        CurveEditor ();
        void setCurveListener (CurveListener* cl) { curve->setCurveListener (cl); }
        void linearPressed ();
        void savePressed ();
        void loadPressed ();
        void setCurve (const std::vector<double>& c);
        std::vector<double> getCurve ();
};


#endif
