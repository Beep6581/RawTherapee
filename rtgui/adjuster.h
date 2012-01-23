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
#ifndef _ADJUSTER_H_
#define _ADJUSTER_H_

#include <gtkmm.h>
#include "editedstate.h"
#include "guiutils.h"

class Adjuster;
class AdjusterListener {

  public:
    virtual void adjusterChanged (Adjuster* a, double newval) {}
};


class Adjuster : public Gtk::VBox {

  protected:
    Gtk::HBox* hbox;
    Gtk::Label* label;
    MyHScale* slider;
    MySpinButton* spin;
    Gtk::Button* reset;
    AdjusterListener* adjusterListener;
    sigc::connection delayConnection;
    sigc::connection spinChange;
    sigc::connection sliderChange;
    sigc::connection editedChange;
    bool listenerReady;
    double defaultVal;			// current default value (it can change when switching from ADD or SET mode)
    double initialDefaultVal;	// default value at construction time
    EditedState editedState;
    EditedState defEditedState;
    int digits;
    Gtk::CheckButton* editedCheckBox;
	bool afterReset;
    bool blocked;
    bool addMode;
    double vMin;
    double vMax;
    double vStep;

    double shapeValue (double a);
    void   refreshLabelStyle ();

  public:

    int delay;

    Adjuster (Glib::ustring label, double vmin, double vmax, double vstep, double vdefault, bool editedCheckBox=false);
    Adjuster (Gtk::Image *imgIcon, double vmin, double vmax, double vstep, double vdefault, bool editedCheckBox=false);
    virtual ~Adjuster ();
    void setAdjusterListener (AdjusterListener* alistener);

    double getValue ();
    int getIntValue ();
    Glib::ustring getTextValue ();
    void setValue (double a);   
    void setLimits (double vmin, double vmax, double vstep, double vdefault);
    void setEnabled (bool enabled);
    void setDefault (double def);
    void setNbDisplayedChars (int nbr) { spin->set_width_chars(nbr); }
    void setEditedState (EditedState eState);
    EditedState getEditedState ();
    void setDefaultEditedState (EditedState eState);
    void showEditedCB ();
    void block(bool isBlocked) { blocked = isBlocked; }
    
    void setAddMode(bool addM);
    bool getAddMode() { return addMode; };
    void spinChanged ();
    void sliderChanged ();
    bool notifyListener ();
    void resetPressed (GdkEventButton* event);
    void editedToggled ();
    double trimValue (double& val);
    float trimValue (float& val);
    int trimValue (int& val);
};

#endif
