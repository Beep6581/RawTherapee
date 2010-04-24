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
#ifndef _SHARPENING_H_
#define _SHARPENING_H_

#include <gtkmm.h>
#include <adjuster.h>
#include <toolpanel.h>

class Sharpening : public Gtk::VBox, public AdjusterListener, public ToolPanel {

  protected:
    Gtk::ComboBoxText* method;
    Adjuster* dradius;
    Adjuster* damount;
    Adjuster* ddamping;
    Adjuster* diter;
    Gtk::VBox* usm;
    Gtk::VBox* rld;
    
    Adjuster* radius;
    Adjuster* amount;
    Adjuster* threshold;
    Adjuster* eradius;
    Adjuster* etolerance;
    Adjuster* hcamount;
    Gtk::VBox* edgebin;
    Gtk::VBox* hcbin;
    Gtk::VBox* edgebox;
    Gtk::VBox* hcbox;
    Gtk::CheckButton* enabled;  
    bool lastEnabled;
    sigc::connection enaConn;
    Gtk::CheckButton* edgesonly;  
    bool lastEdgesOnly;
    sigc::connection eonlyConn;
    Gtk::CheckButton* halocontrol;  
    bool lastHaloControl;
    sigc::connection hcConn;
	bool amountAdd;


	
  public:

    Sharpening ();
    virtual ~Sharpening ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged (Adjuster* a, double newval);
    void enabled_toggled ();
    void edgesonly_toggled ();
    void halocontrol_toggled ();
    void method_changed ();

    void setAdjusterBehavior (bool bamountadd);
};

#endif
