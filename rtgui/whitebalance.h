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
#ifndef _WB_H_
#define _WB_H_

#include <gtkmm.h>
#include <toolpanel.h>
#include <adjuster.h>
#include <wbprovider.h>

class SpotWBListener {

    public: 
        virtual void spotWBRequested (int size) {}
};

class WhiteBalance : public Gtk::VBox, public AdjusterListener, public ToolPanel {

  protected:
    Gtk::ComboBoxText* method;
    Gtk::ComboBoxText* spotsize;
    Adjuster* temp;
    Adjuster* green;
    Gtk::Button* spotbutton;
    int opt;
    double nextTemp;
    double nextGreen;
    WBProvider *wbp;
    SpotWBListener* wblistener;
    sigc::connection methconn;
    bool tempAdd, greenAdd;

  public:

    WhiteBalance ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);
    
    
    void optChanged ();
    void spotPressed ();
    void spotSizeChanged ();
    void adjusterChanged (Adjuster* a, double newval);
    int  getSize (); 
    void setWBProvider (WBProvider* p) { wbp = p; }
    void setSpotWBListener (SpotWBListener* l) { wblistener = l; }
    void setWB (int temp, double green);
    
    void setAdjusterBehavior (bool btempadd, bool bgreenadd);
};

#endif
