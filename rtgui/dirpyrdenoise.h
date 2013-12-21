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
#ifndef _DIRPYRDENOISE_H_
#define _DIRPYRDENOISE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class DirPyrDenoise : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel {

  protected:
	Adjuster* luma;
	Adjuster* Ldetail;
	Adjuster* chroma;
	Adjuster* redchro;
	Adjuster* bluechro;	
	Adjuster* gamma;

    Gtk::CheckButton* enabled;
    bool lastEnabled;
    sigc::connection enaConn;
    Gtk::CheckButton* enhance;
    bool lastenhance;
    sigc::connection enhanConn;
	
//    Gtk::CheckButton* perform;
//    bool lastperform;
//    sigc::connection perfconn;
    MyComboBoxText*   dmethod;
    sigc::connection  dmethodconn;
	

  public:

    DirPyrDenoise ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged  ();
    void enhanceChanged  ();
//    void perform_toggled  ();
    void dmethodChanged      ();

    void setAdjusterBehavior (bool lumaadd, bool lumdetadd, bool chromaadd, bool chromaredadd, bool chromablueadd, bool gammaadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
};

#endif
