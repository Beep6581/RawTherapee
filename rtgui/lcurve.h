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
#ifndef _LCURVE_H_
#define _LCURVE_H_

#include <gtkmm.h>
#include <adjuster.h>
#include <toolpanel.h>
#include <curveeditor.h>
#include <mycurve.h>

class LCurve : public Gtk::VBox, public AdjusterListener, public ToolPanel, public CurveListener {

  protected:
    Adjuster* brightness;
    Adjuster* black;
    Adjuster* contrast;
    Adjuster* hlcompr;
    Adjuster* shcompr;
    CurveEditor* shape;
    Gtk::Expander* curvexp;

  public:

    LCurve ();

    void read           (const rtengine::procparams::ProcParams* pp);
    void write          (rtengine::procparams::ProcParams* pp);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams);   

    void curveChanged ();
    void adjusterChanged (Adjuster* a, double newval);
    void expandCurve (bool isExpanded);
    bool isCurveExpanded ();
};

#endif
