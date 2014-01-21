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
 *  © 2010 Emil Martinec <ejmartin@uchicago.edu>
 */
 
#ifndef DIRPYREQUALIZER_H_INCLUDED
#define DIRPYREQUALIZER_H_INCLUDED

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class DirPyrEqualizer : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel
{

protected:

    Gtk::CheckButton * enabled;
    Adjuster* multiplier[5]; 
    Adjuster* threshold;

    sigc::connection enaConn;
    sigc::connection lumaneutralPressedConn;
    sigc::connection lumacontrastPlusPressedConn;
    sigc::connection lumacontrastMinusPressedConn;
	
    bool lastEnabled;

public:

    DirPyrEqualizer ();
    virtual ~DirPyrEqualizer ();

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode        (bool batchMode);
    void setAdjusterBehavior (bool multiplieradd, bool thresholdadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledToggled ();
    void lumaneutralPressed ();
    void lumacontrastPlusPressed ();
    void lumacontrastMinusPressed ();
};

#endif
