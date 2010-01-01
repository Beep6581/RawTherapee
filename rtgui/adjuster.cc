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
#include <adjuster.h>
#include <sigc++/class_slot.h>
#include <math.h>
#include <multilangmgr.h>

extern Glib::ustring argv0;

int Adjuster::delay = 1000;

Adjuster::Adjuster (Glib::ustring vlabel, double vmin, double vmax, double vstep, double vdefault, bool editedcb) {

  adjusterListener = NULL;
  afterReset = false;

  set_border_width (2);

  hbox = Gtk::manage (new Gtk::HBox ());

  label = Gtk::manage (new Gtk::Label (vlabel, Gtk::ALIGN_LEFT));

  if (editedcb) {
    editedCheckBox = new Gtk::CheckButton ();
    editedChange = editedCheckBox->signal_toggled().connect( sigc::mem_fun(*this, &Adjuster::editedToggled) );
	hbox->pack_start (*editedCheckBox);
  }
  else
    editedCheckBox = NULL;
    
  hbox->pack_start (*label);

  reset = Gtk::manage (new Gtk::Button ());
  reset->add (*Gtk::manage (new Gtk::Image (argv0+"/images/undo.png")));
  reset->set_relief (Gtk::RELIEF_NONE);
  reset->set_border_width (0);
  reset->set_tooltip_text (M("ADJUSTER_RESET_TO_DEFAULT"));

  hbox->pack_end (*reset, Gtk::PACK_SHRINK, 0); 
  
  spin = Gtk::manage (new Gtk::SpinButton ());
  spin->set_size_request (70, -1);

  hbox->pack_end (*spin, Gtk::PACK_SHRINK, 0);

  reset->set_size_request (-1, spin->get_height());
  
  slider = Gtk::manage (new Gtk::HScale ());
  slider->set_draw_value (false);

  pack_start (*hbox, false, false);
  pack_start (*slider, false, false);

  setLimits (vmin, vmax, vstep, vdefault);
  
  defaultVal = shapeValue (vdefault);
  editedState = defEditedState = Irrelevant;

  sliderChange = slider->signal_value_changed().connect( sigc::mem_fun(*this, &Adjuster::sliderChanged) );
  spinChange = spin->signal_value_changed().connect ( sigc::mem_fun(*this, &Adjuster::spinChanged), true);
  reset->signal_clicked().connect( sigc::mem_fun(*this, &Adjuster::resetPressed) );
  slider->set_update_policy (Gtk::UPDATE_CONTINUOUS);
  
  show_all ();
}

Adjuster::~Adjuster () {

    sliderChange.block (true);
    spinChange.block (true);
    delayConnection.block (true);
    adjusterListener = NULL;
}

void Adjuster::setDefault (double def) {

    defaultVal = shapeValue (def);
}

void Adjuster::setDefaultEditedState (EditedState eState) {

    defEditedState = eState;
}

void Adjuster::resetPressed () {

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
    slider->set_value (defaultVal);
}

double Adjuster::shapeValue (double a) {

  return round(a*pow(10, digits)) / pow(10, digits);
}

void Adjuster::setLimits (double vmin, double vmax, double vstep, double vdefault) {

  sliderChange.block (true);
  spinChange.block (true);
  for (digits=0; fabs(vstep*pow(10,digits)-floor(vstep*pow(10,digits)))>0.000000000001; digits++);
  spin->set_digits (digits);
  spin->set_increments (vstep, 2.0*vstep);
  spin->set_range (vmin, vmax);
  spin->set_value (shapeValue(vdefault));
  slider->set_digits (digits);
  slider->set_increments (vstep, 2.0*vstep);
  slider->set_range (vmin, vmax);
  slider->set_value (shapeValue(vdefault));
  defaultVal = shapeValue (vdefault);
  sliderChange.block (false);
  spinChange.block (false);
}

void Adjuster::setAdjusterListener (AdjusterListener* alistener) {

  adjusterListener = alistener;  
}

void Adjuster::spinChanged () {

  sliderChange.block (true);
  slider->set_value (spin->get_value ());
  sliderChange.block (false);

  if (delay==0) {
    if (adjusterListener!=NULL)
        adjusterListener->adjusterChanged (this, spin->get_value ());
  }
  else
    Glib::signal_idle().connect (sigc::mem_fun(*this, &Adjuster::notifyListener));
    
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
  spin->set_value (slider->get_value ());
  spinChange.block (false);

  if (delay==0) {
    if (adjusterListener)
        adjusterListener->adjusterChanged (this, spin->get_value ());
  }
  else
    delayConnection = Glib::signal_timeout().connect (sigc::mem_fun(*this, &Adjuster::notifyListener), delay);

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
  slider->set_value (shapeValue (a));
  sliderChange.block (false);
  spinChange.block (false);
  afterReset = false;
}

double Adjuster::getValue () {

  return spin->get_value ();
}

bool Adjuster::notifyListener () {

  gdk_threads_enter();
  
  if (adjusterListener!=NULL)
    adjusterListener->adjusterChanged (this, spin->get_value ());
  gdk_threads_leave();

  return false;
}

void Adjuster::setEnabled (bool enabled) {

    spin->set_sensitive (enabled);
    slider->set_sensitive (enabled);
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

    if (!editedCheckBox) {
        editedCheckBox = new Gtk::CheckButton ();
        hbox->pack_start (*editedCheckBox, Gtk::PACK_SHRINK, 2);
        hbox->reorder_child (*editedCheckBox, 0);
	    editedChange = editedCheckBox->signal_toggled().connect( sigc::mem_fun(*this, &Adjuster::editedToggled) );
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
	
    if (adjusterListener)
        adjusterListener->adjusterChanged (this, spin->get_value ());
}
