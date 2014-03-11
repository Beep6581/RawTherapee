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
#include "adjuster.h"
#include <sigc++/class_slot.h>
#include <cmath>
#include "multilangmgr.h"
#include "../rtengine/rtengine.h"
#include "options.h"
#include "guiutils.h"
#include "rtimage.h"

#define MIN_RESET_BUTTON_HEIGHT 17

static double one2one(double val) { return val; }

Adjuster::Adjuster (Glib::ustring vlabel, double vmin, double vmax, double vstep, double vdefault, Gtk::Image *imgIcon1, Gtk::Image *imgIcon2, double2double_fun slider2value_, double2double_fun value2slider_) {

  Gtk::HBox *hbox2=NULL;
  label = NULL;
  adjusterListener = NULL;
  afterReset = false;
  blocked = false;
  automatic = NULL;
  eventPending = false;

  slider2value = slider2value_ ? slider2value_ : one2one;
  value2slider = value2slider_ ? value2slider_ : one2one;
  vMin = vmin;
  vMax = vmax;
  vStep = vstep;
  addMode = false;

  // TODO: let the user chose the default value of Adjuster::delay, for slow machines
  delay = options.adjusterDelay;		// delay is no more static, so we can set the delay individually (useful for the RAW editor tab)

  set_border_width (0);
  set_spacing (2);

  hbox = Gtk::manage (new Gtk::HBox ());
  hbox->set_border_width(0);
  hbox->set_spacing(2);

  editedCheckBox = NULL;

  if (!vlabel.empty()) {
    adjustmentName = vlabel;
    label = Gtk::manage (new Gtk::Label (adjustmentName, Gtk::ALIGN_LEFT));
  }

  reset = Gtk::manage (new Gtk::Button ());
  reset->add (*Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png")));
  reset->set_relief (Gtk::RELIEF_NONE);
  reset->set_border_width (0);
  reset->set_tooltip_text (M("ADJUSTER_RESET_TO_DEFAULT"));

  spin = Gtk::manage (new MySpinButton ());
  spin->set_has_frame(false);
  spin->set_name("FramelessSpinButton");

  reset->set_size_request (-1, spin->get_height() > MIN_RESET_BUTTON_HEIGHT ? spin->get_height(): MIN_RESET_BUTTON_HEIGHT);

  slider = Gtk::manage (new MyHScale ());
  slider->set_draw_value (false);

  pack_start (*hbox, true, true);

  if (vlabel.empty()) {
    // No label, everything goes in hbox
    if (imgIcon1) hbox->pack_start (*imgIcon1, Gtk::PACK_SHRINK, 0);
    hbox->pack_end (*reset, Gtk::PACK_SHRINK, 0); 
    hbox->pack_end (*spin, Gtk::PACK_SHRINK, 0);
    if (imgIcon2) hbox->pack_start (*imgIcon2, Gtk::PACK_SHRINK, 0);
    hbox->pack_start (*slider, true, true);
  }
  else {
    // A label is provided, spreading the widgets in 2 rows
    hbox->pack_start (*label);
    hbox->pack_end (*reset, Gtk::PACK_SHRINK, 0);
    hbox->pack_end (*spin, Gtk::PACK_SHRINK, 0);
    if (!imgIcon1 || !imgIcon2) {
      pack_start (*slider, true, true);
    }
    else {
      // A second HBox is necessary
      hbox2 = Gtk::manage (new Gtk::HBox());
      if (imgIcon1) hbox2->pack_start (*imgIcon1, Gtk::PACK_SHRINK, 0);
      hbox2->pack_start (*slider, true, true);
      if (imgIcon2) hbox2->pack_start (*imgIcon2, Gtk::PACK_SHRINK, 0);
      pack_start (*hbox2, true, true);
    }
  }

  setLimits (vmin, vmax, vstep, vdefault);
  
  defaultVal = shapeValue (vdefault);
  ctorDefaultVal = shapeValue (vdefault);
  editedState = defEditedState = Irrelevant;
  autoState = Irrelevant;

  sliderChange = slider->signal_value_changed().connect( sigc::mem_fun(*this, &Adjuster::sliderChanged) );
  spinChange = spin->signal_value_changed().connect ( sigc::mem_fun(*this, &Adjuster::spinChanged), true);
  reset->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &Adjuster::resetPressed) );
  slider->set_update_policy (Gtk::UPDATE_CONTINUOUS);
  
  show_all ();
}

Adjuster::~Adjuster () {

    sliderChange.block (true);
    spinChange.block (true);
    delayConnection.block (true);
    adjusterListener = NULL;
    if (automatic) delete automatic;
}

void Adjuster::addAutoButton (Glib::ustring tooltip) {
    if (!automatic) {
        automatic = new Gtk::CheckButton ();
        //automatic->add (*Gtk::manage (new RTImage ("processing.png")));
        automatic->set_border_width (0);
        automatic->set_tooltip_markup(tooltip.length() ? Glib::ustring::compose("<b>%1</b>\n\n%2", M("GENERAL_AUTO"), tooltip) : M("GENERAL_AUTO"));
        autoChange = automatic->signal_toggled().connect( sigc::mem_fun(*this, &Adjuster::autoToggled) );

        hbox->pack_end (*automatic, Gtk::PACK_SHRINK, 0);
        hbox->reorder_child (*automatic, 0);
    }
}

void Adjuster::delAutoButton () {
    if (automatic) {
        removeIfThere(hbox, automatic);
        delete automatic;
        automatic = NULL;
    }
}

void Adjuster::throwOnButtonRelease(bool throwOnBRelease) {

    if (throwOnBRelease) {
        if (!buttonReleaseSlider.connected())
            buttonReleaseSlider = slider->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &Adjuster::sliderReleased) );
        if (!buttonReleaseSpin.connected())
            buttonReleaseSpin = spin->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &Adjuster::spinReleased) ); // Use the same callback hook
    }
    else {
        if (buttonReleaseSlider.connected())
            buttonReleaseSlider.disconnect();
        if (buttonReleaseSpin.connected())
            buttonReleaseSpin.disconnect();
    }
    eventPending = false;
}

void Adjuster::setDefault (double def) {

    defaultVal = shapeValue (def);
}

void Adjuster::setDefaultEditedState (EditedState eState) {

    defEditedState = eState;
}

void Adjuster::autoToggled () {

    if (!editedCheckBox) {
        // If not used in the BatchEditor panel
        if (automatic->get_active()) {
            // Disable the slider and spin button
            spin->set_sensitive(false);
            slider->set_sensitive(false);
        }
        else {
            // Enable the slider and spin button
            spin->set_sensitive(true);
            slider->set_sensitive(true);
        }
    }

    if (adjusterListener!=NULL && !blocked) {
        adjusterListener->adjusterAutoToggled(this, automatic->get_active());
    }
}

void Adjuster::sliderReleased (GdkEventButton* event) {

    if ((event != NULL) && (event->button == 1)) {
        if (delayConnection.connected())
            delayConnection.disconnect ();
        notifyListener();
    }
}

void Adjuster::spinReleased (GdkEventButton* event) {

    if ((event != NULL) && delay==0) {
        if (delayConnection.connected())
            delayConnection.disconnect ();
        notifyListener();
    }
}

void Adjuster::resetValue (bool toInitial) {
    if (editedState!=Irrelevant) {
        editedState = defEditedState;
        if (editedCheckBox) {
            editedChange.block (true);
            editedCheckBox->set_active (defEditedState==Edited);
            editedChange.block (false);
        }
        refreshLabelStyle ();
    }
    afterReset = true;
    if (toInitial) {
        // resetting to the initial editing value, when the image has been loaded
	slider->set_value (addMode ? defaultVal : value2slider(defaultVal));
    }
    else {
        // resetting to the slider default value
        if (addMode)
            slider->set_value (0.);
        else
            slider->set_value (value2slider(ctorDefaultVal));
    }
}

// Please note that it won't change the "Auto" CheckBox's state, if there
void Adjuster::resetPressed (GdkEventButton* event) {

    if ((event != NULL) && (event->state & GDK_CONTROL_MASK) && (event->button == 1))
        resetValue(true);
    else
        resetValue(false);
}

double Adjuster::shapeValue (double a) {

  return round(a*pow(double(10), digits)) / pow(double(10), digits);
}

void Adjuster::setLimits (double vmin, double vmax, double vstep, double vdefault) {

  sliderChange.block (true);
  spinChange.block (true);
  for (digits=0; fabs(vstep*pow(double(10),digits)-floor(vstep*pow(double(10),digits)))>0.000000000001; digits++);
  spin->set_digits (digits);
  spin->set_increments (vstep, 2.0*vstep);
  spin->set_range (vmin, vmax);
  spin->updateSize();
  spin->set_value (shapeValue(vdefault));
  slider->set_digits (digits);
  slider->set_increments (vstep, 2.0*vstep);
  slider->set_range (addMode ? vmin : value2slider(vmin), addMode ? vmax : value2slider(vmax));
  slider->set_value (addMode ? shapeValue(vdefault) : value2slider(shapeValue(vdefault)));
  //defaultVal = shapeValue (vdefault);
  sliderChange.block (false);
  spinChange.block (false);
}

void Adjuster::setAddMode(bool addM) {
	if (addM != addMode) {
		// Switching the Adjuster to the new mode
		addMode = addM;
		if (addM) {
			// Switching to the relative mode
			double range = -vMin + vMax;
			if (range < 0.) range = -range;
			setLimits(-range, range, vStep, 0);
		}
		else {
			// Switching to the absolute mode
			setLimits(vMin, vMax, vStep, defaultVal);
		}
	}
}

void Adjuster::spinChanged () {

  if (delayConnection.connected())
    delayConnection.disconnect ();

  sliderChange.block (true);
  slider->set_value (addMode ? spin->get_value () : value2slider(spin->get_value ()));
  sliderChange.block (false);

  if (delay==0) {
    if (adjusterListener && !blocked) {
        if (!buttonReleaseSlider.connected() || afterReset) {
            eventPending = false;
            adjusterListener->adjusterChanged (this, spin->get_value ());
        }
        else eventPending = true;
    }
  }
  else {
    eventPending = true;
    delayConnection = Glib::signal_timeout().connect (sigc::mem_fun(*this, &Adjuster::notifyListener), delay);
  }

  if (editedState==UnEdited) {
    editedState = Edited;
    if (editedCheckBox) {
        editedChange.block (true);
        editedCheckBox->set_active (true);
        editedChange.block (false);
    }
    refreshLabelStyle ();
  }
  afterReset = false;
}

void Adjuster::sliderChanged () {
  
  if (delayConnection.connected())
    delayConnection.disconnect ();

  spinChange.block (true);
  spin->set_value (addMode ? slider->get_value () : slider2value(slider->get_value ()));
  spinChange.block (false);

  if (delay==0 || afterReset) {
    if (adjusterListener && !blocked) {
      if (!buttonReleaseSlider.connected() || afterReset) {
        eventPending = false;
        adjusterListener->adjusterChanged (this, spin->get_value ());
      }
      else eventPending = true;
    }
  }
  else {
    eventPending = true;
    delayConnection = Glib::signal_timeout().connect (sigc::mem_fun(*this, &Adjuster::notifyListener), delay);
  }

  if (!afterReset && editedState==UnEdited) {
    editedState = Edited;
    if (editedCheckBox) {
        editedChange.block (true);
        editedCheckBox->set_active (true);
        editedChange.block (false);
    }
    refreshLabelStyle ();
  }
  afterReset = false;
}

void Adjuster::setValue (double a) {

  spinChange.block (true);
  sliderChange.block (true);
  spin->set_value (shapeValue (a));
  slider->set_value (addMode ? shapeValue(a) : value2slider(shapeValue (a)));
  sliderChange.block (false);
  spinChange.block (false);
  afterReset = false;
}

void Adjuster::setAutoValue (bool a) {
    if (automatic) {
        bool oldVal = autoChange.block(true);
        automatic->set_active(a);
        autoChange.block(oldVal);
        if (!editedCheckBox) {
            // If not used in the BatchEditor panel
            if (a) {
                // Disable the slider and spin button
                spin->set_sensitive(false);
                slider->set_sensitive(false);
            }
            else {
                // Enable the slider and spin button
                spin->set_sensitive(true);
                slider->set_sensitive(true);
            }
        }
    }
}

bool Adjuster::notifyListener () {

  if (eventPending && adjusterListener!=NULL && !blocked) {
    adjusterListener->adjusterChanged (this, spin->get_value ());
  }
  eventPending = false;

  return false;
}

bool Adjuster::notifyListenerAutoToggled () {

  if (adjusterListener!=NULL && !blocked) {
    adjusterListener->adjusterAutoToggled(this, automatic->get_active());
  }
  return false;
}

void Adjuster::setEnabled (bool enabled) {

    bool autoVal = automatic && !editedCheckBox ? automatic->get_active() : true;
    spin->set_sensitive (enabled && autoVal);
    slider->set_sensitive (enabled && autoVal);
    if (automatic)
        automatic->set_sensitive (enabled);
}

void Adjuster::setEditedState (EditedState eState) {

    if (editedState!=eState) {
        if (editedCheckBox) {
            editedChange.block (true);
            editedCheckBox->set_active (eState==Edited);
            editedChange.block (false);
        }
        editedState = eState;
        refreshLabelStyle ();
    }
}

EditedState Adjuster::getEditedState () {

    if (editedState!=Irrelevant && editedCheckBox)
        editedState = editedCheckBox->get_active () ? Edited : UnEdited;
    return editedState;
}

void Adjuster::showEditedCB () {

    if (label)
        removeIfThere(hbox, label, false);

    if (!editedCheckBox) {
        editedCheckBox = Gtk::manage(new Gtk::CheckButton (adjustmentName));
        hbox->pack_start (*editedCheckBox, Gtk::PACK_SHRINK, 2);
        hbox->reorder_child (*editedCheckBox, 0);
        editedChange = editedCheckBox->signal_toggled().connect( sigc::mem_fun(*this, &Adjuster::editedToggled) );
        editedCheckBox->show();
    }
}

void Adjuster::refreshLabelStyle () {

/*	Glib::RefPtr<Gtk::Style> style = label->get_style ();
    Pango::FontDescription fd = style->get_font ();
    fd.set_weight (editedState==Edited ? Pango::WEIGHT_BOLD : Pango::WEIGHT_NORMAL);
    style->set_font (fd);
    label->set_style (style);
    label->queue_draw ();*/
}

void Adjuster::editedToggled () {

    if (adjusterListener && !blocked) {
        adjusterListener->adjusterChanged (this, spin->get_value ());
    }
    eventPending = false;
}

double Adjuster::trimValue (double& val) {

    if      (val > vMax) val = vMax;  // shapeValue(vMax) ?
    else if (val < vMin) val = vMin;  // shapeValue(vMin) ?
    return val;
}

int Adjuster::trimValue (int& val) {

    if      (val > (int)vMax) val = (int)vMax;  // shapeValue(vMax) ?
    else if (val < (int)vMin) val = (int)vMin;  // shapeValue(vMin) ?
    return val;
}

float Adjuster::trimValue (float& val) {

    if      (val > (float)vMax) val = (float)vMax;  // shapeValue(vMax) ?
    else if (val < (float)vMin) val = (float)vMin;  // shapeValue(vMin) ?
    return val;
}
