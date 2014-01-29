/*
 *  This file is part of RawTherapee.
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
 *
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 */
 
#ifndef HSVEQUALIZER_H_INCLUDED
#define HSVEQUALIZER_H_INCLUDED

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "colorprovider.h"


class HSVEqualizer : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel, public CurveListener, public ColorProvider
{

protected:

    Gtk::CheckButton * enabled;

	CurveEditorGroup*  curveEditorG;
	FlatCurveEditor*   hshape;
	FlatCurveEditor*   sshape;
	FlatCurveEditor*   vshape;

public:

    HSVEqualizer ();
    virtual ~HSVEqualizer ();

    void read            (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
    void write           (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void curveChanged    (CurveEditor* ce);
    //void setDefaults     (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode    (bool batchMode);
    void setEditProvider (EditDataProvider *provider);
    void autoOpenCurve   ();
    virtual void colorForValue (double valX, double valY, int callerId, ColorCaller* caller);
   
    //void adjusterChanged (Adjuster* a, double newval);
};

#endif
