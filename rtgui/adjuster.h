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
class AdjusterListener
{

public:
    virtual ~AdjusterListener() {};
    virtual void adjusterChanged (Adjuster* a, double newval) {}
    virtual void adjusterAutoToggled (Adjuster* a, bool newval) {}
};

typedef double(*double2double_fun)(double val);

class Adjuster : public Gtk::VBox
{

protected:
    Glib::ustring adjustmentName;
    Gtk::HBox* hbox;
    Gtk::Label* label;
    MyHScale* slider;
    MySpinButton* spin;
    Gtk::Button* reset;
    Gtk::CheckButton* automatic;
    AdjusterListener* adjusterListener;
    sigc::connection delayConnection;
    sigc::connection spinChange;
    sigc::connection sliderChange;
    sigc::connection editedChange;
    sigc::connection autoChange;
    sigc::connection buttonReleaseSlider;
    sigc::connection buttonReleaseSpin;
    bool listenerReady;
    double defaultVal;          // current default value (it can change when switching from ADD or SET mode)
    double ctorDefaultVal;      // default value at construction time
    EditedState editedState;
    EditedState defEditedState;
    EditedState autoState;
    EditedState defAutoState;
    int digits;
    Gtk::CheckButton* editedCheckBox;
    bool afterReset;
    bool blocked;
    bool addMode;
    bool eventPending;
    double vMin;
    double vMax;
    double vStep;

    double shapeValue (double a);
    void   refreshLabelStyle ();
    double2double_fun value2slider, slider2value;

public:

    int delay;

    Adjuster (Glib::ustring vlabel, double vmin, double vmax, double vstep, double vdefault, Gtk::Image *imgIcon1 = NULL, Gtk::Image *imgIcon2 = NULL, double2double_fun slider2value = NULL, double2double_fun value2slider = NULL);
    virtual ~Adjuster ();

    // Add an "Automatic" checkbox next to the reset button.
    void addAutoButton(Glib::ustring tooltip = "");
    // Remove the "Automatic" checkbox next to the reset button.
    void delAutoButton();
    // Send back the value of og the Auto checkbox
    bool getAutoValue ()
    {
        return automatic != NULL ? automatic->get_active () : false;
    }
    void setAutoValue (bool a);
    bool notifyListenerAutoToggled ();
    void autoToggled ();
    void setAutoInconsistent (bool i)
    {
        if (automatic) {
            automatic->set_inconsistent(i);
        }
    }
    bool getAutoInconsistent ()
    {
        return automatic ? automatic->get_inconsistent() : true /* we have to return something */;
    }

    void setAdjusterListener (AdjusterListener* alistener)
    {
        adjusterListener = alistener;
    }

    // return the value trimmed to the limits at construction time
    double getValue ()
    {
        return shapeValue(spin->get_value ());
    }
    // return the value trimmed to the limits at construction time
    int getIntValue ()
    {
        return spin->get_value_as_int ();
    }
    // return the value trimmed to the limits at construction time,
    // method only used by the history manager, so decoration is added if addMode=true
    Glib::ustring getTextValue ()
    {
        if (addMode) {
            return Glib::ustring::compose("<i>%1</i>", spin->get_text ());
        } else {
            return spin->get_text ();
        }
    }

    void setLabel (Glib::ustring lbl)
    {
        label->set_label(lbl);
    }
    void setValue (double a);
    void setLimits (double vmin, double vmax, double vstep, double vdefault);
    void setEnabled (bool enabled);
    void setDefault (double def);
    // will let the adjuster throw it's "changed" signal when the mouse button is released. Can work altogether with the delay value.
    void throwOnButtonRelease(bool throwOnBRelease = true);
    void setNbDisplayedChars (int nbr)
    {
        spin->set_width_chars(nbr);
    }
    void setEditedState (EditedState eState);
    EditedState getEditedState ();
    void setDefaultEditedState (EditedState eState);
    void showEditedCB ();
    bool block(bool isBlocked)
    {
        bool oldValue = blocked;
        blocked = isBlocked;
        return oldValue;
    }

    void setAddMode(bool addM);
    bool getAddMode()
    {
        return addMode;
    };
    void spinChanged ();
    void sliderChanged ();
    bool notifyListener ();
    void sliderReleased (GdkEventButton* event);
    void spinReleased (GdkEventButton* event);
    void resetValue (bool toInitial);
    void resetPressed (GdkEventButton* event);
    void editedToggled ();
    double trimValue (double val);
    float trimValue (float val);
    int trimValue (int val);
};

#endif
