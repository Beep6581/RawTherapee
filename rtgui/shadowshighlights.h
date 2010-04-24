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
#ifndef _SHADOWSHIGHLIGHTS_H_
#define _SHADOWSHIGHLIGHTS_H_

#include <gtkmm.h>
#include <adjuster.h>
#include <toolpanel.h>

class ShadowsHighlights : public Gtk::VBox, public AdjusterListener, public ToolPanel {

  protected:
    Adjuster* highlights;
    Adjuster* h_tonalwidth;
    Adjuster* shadows;
    Adjuster* s_tonalwidth;
    Adjuster* lcontrast;
    Adjuster* radius;
    Gtk::CheckButton* enabled;
    Gtk::CheckButton* hq;
    bool                hAdd, sAdd, lcAdd;
    bool                lastEnabled, lastHQ;
    sigc::connection    enaConn, hqConn;

  public:

    ShadowsHighlights ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void hqChanged ();
    
    void setAdjusterBehavior (bool hadd, bool sadd, bool lcadd);
};

#endif
