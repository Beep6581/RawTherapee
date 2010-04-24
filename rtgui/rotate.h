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
#ifndef _ROTATE_H_
#define _ROTATE_H_

#include <gtkmm.h>
#include <adjuster.h>
#include <toolpanel.h>

class RotateListener {

    public:
        virtual void straightenRequested () {}
        virtual void autoCropRequested   () {}
};

class Rotate : public Gtk::VBox, public AdjusterListener, public ToolPanel {

  protected:
    Adjuster*           degree;
    Gtk::Button*        selectStraight;
    Gtk::Button*        autoCrop;
    RotateListener*     rlistener;
    Gtk::CheckButton*   fill;
    bool                degAdd;
    bool                lastFill;
    sigc::connection    fillConn;

  public:

    Rotate ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void straighten (double deg);

    void adjusterChanged        (Adjuster* a, double newval);
    void setAdjusterBehavior    (bool brotadd);
    void selectStraightPressed  ();
    void fillPressed            ();
    void autoCropPressed        ();
    void setRotateListener      (RotateListener* l) { rlistener = l; }
};

#endif
