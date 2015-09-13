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
#include "thresholdadjuster.h"
#include <sigc++/class_slot.h>
#include <cmath>
#include "multilangmgr.h"
#include "../rtengine/rtengine.h"
#include "options.h"
#include "guiutils.h"
#include "rtimage.h"

#define MIN_RESET_BUTTON_HEIGHT 17

ThresholdAdjuster::ThresholdAdjuster (Glib::ustring label,
                                      double minValueBottom, double maxValueBottom, double defBottom, Glib::ustring labelBottom, unsigned int precisionBottom,
                                      double minValueTop,    double maxValueTop,    double defTop,    Glib::ustring labelTop,    unsigned int precisionTop,
                                      ThresholdCurveProvider* curveProvider, bool editedCheckBox)
    : tSelector(minValueBottom, maxValueBottom, defBottom, labelBottom, precisionBottom, minValueTop, maxValueTop, defTop, labelTop, precisionTop, curveProvider)

{
    initialDefaultVal[ThresholdSelector::TS_BOTTOMLEFT] = defBottom;
    initialDefaultVal[ThresholdSelector::TS_TOPLEFT] = defTop;
    initialDefaultVal[ThresholdSelector::TS_BOTTOMRIGHT] = 0.; // unused
    initialDefaultVal[ThresholdSelector::TS_TOPRIGHT] = 0.;    // unused

    initObject (label, editedCheckBox);
}

ThresholdAdjuster::ThresholdAdjuster (Glib::ustring label, double minValue, double maxValue, double defBottom,
                                      double defTop, unsigned int precision, bool startAtOne, bool editedCheckBox)
    : tSelector(minValue, maxValue, defBottom, defTop, precision, startAtOne)
{
    initialDefaultVal[ThresholdSelector::TS_BOTTOMLEFT] = defBottom;
    initialDefaultVal[ThresholdSelector::TS_TOPLEFT] = defTop;
    initialDefaultVal[ThresholdSelector::TS_BOTTOMRIGHT] = maxValue;
    initialDefaultVal[ThresholdSelector::TS_TOPRIGHT] = maxValue;

    initObject (label, editedCheckBox);
}

ThresholdAdjuster::ThresholdAdjuster (Glib::ustring label, double minValue, double maxValue,
                                      double defBottomLeft, double defTopLeft, double defBottomRight, double defTopRight,
                                      unsigned int precision, bool startAtOne, bool editedCheckBox)
    : tSelector(minValue, maxValue, defBottomLeft, defTopLeft,
                defBottomRight, defTopRight, precision, startAtOne)
{
    initialDefaultVal[ThresholdSelector::TS_BOTTOMLEFT] = defBottomLeft;
    initialDefaultVal[ThresholdSelector::TS_TOPLEFT] = defTopLeft;
    initialDefaultVal[ThresholdSelector::TS_BOTTOMRIGHT] = defBottomRight;
    initialDefaultVal[ThresholdSelector::TS_TOPRIGHT] = defTopRight;

    initObject (label, editedCheckBox);
}

void ThresholdAdjuster::initObject (Glib::ustring label, bool editedcb)
{

    adjusterListener = NULL;
    afterReset = false;
    blocked = false;

    addMode = false;

    delay = options.adjusterMinDelay;

    set_name("ThresholdAdjuster");

    hbox = Gtk::manage (new Gtk::HBox ());

    this->label = Gtk::manage (new Gtk::Label (label, Gtk::ALIGN_LEFT));

    if (editedcb) {
        editedCheckBox = Gtk::manage (new Gtk::CheckButton ());
        editedChange = editedCheckBox->signal_toggled().connect( sigc::mem_fun(*this, &ThresholdAdjuster::editedToggled) );
        hbox->pack_start (*editedCheckBox);
    } else {
        editedCheckBox = NULL;
    }

    hbox->pack_start (*this->label);

    reset = Gtk::manage (new Gtk::Button ());
    reset->add (*Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png")));
    reset->set_relief (Gtk::RELIEF_NONE);
    reset->set_border_width (0);
    reset->set_tooltip_text (M("ADJUSTER_RESET_TO_DEFAULT"));

    hbox->pack_end (*reset, Gtk::PACK_SHRINK, 0);

    reset->set_size_request (-1, this->label->get_height() > MIN_RESET_BUTTON_HEIGHT ? this->label->get_height() : MIN_RESET_BUTTON_HEIGHT);

    pack_start (*hbox, false, false);
    pack_start (tSelector, false, false);

    editedState = defEditedState = Irrelevant;

    selectorChange = tSelector.signal_value_changed().connect( sigc::mem_fun(*this, &ThresholdAdjuster::selectorChanged) );
    reset->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ThresholdAdjuster::resetPressed) );

    show_all ();
}

ThresholdAdjuster::~ThresholdAdjuster ()
{

    selectorChange.disconnect();
    delayConnection.block(true);
    adjusterListener = NULL;
}

void ThresholdAdjuster::setDefault (double bottom, double top)
{

    selectorChange.block (true);
    tSelector.setPositions(shapeValue(bottom), shapeValue(top));
    selectorChange.block (false);
}

void ThresholdAdjuster::setDefault (double bottomLeft, double topLeft, double bottomRight, double topRight)
{

    selectorChange.block (true);
    tSelector.setPositions(shapeValue(bottomLeft), shapeValue(topLeft), shapeValue(bottomRight), shapeValue(topRight));
    selectorChange.block (false);
}

void ThresholdAdjuster::setDefaultEditedState (EditedState eState)
{

    defEditedState = eState;
}

void ThresholdAdjuster::resetPressed (GdkEventButton* event)
{

    if (editedState != Irrelevant) {
        editedState = defEditedState;

        if (editedCheckBox) {
            editedChange.block (true);
            editedCheckBox->set_active (defEditedState == Edited);
            editedChange.block (false);
        }

        refreshLabelStyle ();
    }

    afterReset = true;

    if ((event != NULL) && (event->state & GDK_CONTROL_MASK) && (event->button == 1))
        // CTRL pressed : resetting to current default value
    {
        tSelector.reset();
    } else
        // no modifier key or addMode=true : resetting to initial default value
        tSelector.setPositions(initialDefaultVal[ThresholdSelector::TS_BOTTOMLEFT],
                               initialDefaultVal[ThresholdSelector::TS_TOPLEFT],
                               initialDefaultVal[ThresholdSelector::TS_BOTTOMRIGHT],
                               initialDefaultVal[ThresholdSelector::TS_TOPRIGHT]);
}

double ThresholdAdjuster::shapeValue (double a)
{

    unsigned int digit = tSelector.getPrecision();
    return round(a * pow(double(10), digit)) / pow(double(10), digit);
}

void ThresholdAdjuster::selectorChanged ()
{

    if (delayConnection.connected()) {
        delayConnection.disconnect ();
    }

    if (delay == 0) {
        if (adjusterListener && !blocked) {
            sendToListener ();
        }
    } else {
        delayConnection = Glib::signal_timeout().connect (sigc::mem_fun(*this, &ThresholdAdjuster::notifyListener), delay);
    }

    if (!afterReset && editedState == UnEdited) {
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

void ThresholdAdjuster::setValue (double bottom, double top)
{

    selectorChange.block (true);
    tSelector.setPositions(bottom, top);
    selectorChange.block (false);
    afterReset = false;
}

void ThresholdAdjuster::setValue (double bottomLeft, double topLeft, double bottomRight, double topRight)
{

    selectorChange.block (true);
    tSelector.setPositions(bottomLeft, topLeft, bottomRight, topRight);
    selectorChange.block (false);
    afterReset = false;
}

void ThresholdAdjuster::getValue (double& bottom, double& top)
{
    tSelector.getPositions<double> (bottom, top);
}
void ThresholdAdjuster::getValue (double& bottomLeft, double& topLeft, double& bottomRight, double& topRight)
{
    tSelector.getPositions<double> (bottomLeft, topLeft, bottomRight, topRight);
}
void ThresholdAdjuster::getValue (int& bottom, int& top)
{
    tSelector.getPositions<int> (bottom, top);
}
void ThresholdAdjuster::getValue (int& bottomLeft, int& topLeft, int& bottomRight, int& topRight)
{
    tSelector.getPositions<int> (bottomLeft, topLeft, bottomRight, topRight);
}

void ThresholdAdjuster::getValue (Glib::ustring& bottom, Glib::ustring& top)
{
    tSelector.getPositions (bottom, top);
}

void ThresholdAdjuster::getValue (Glib::ustring& bottomLeft, Glib::ustring& topLeft, Glib::ustring& bottomRight, Glib::ustring& topRight)
{
    tSelector.getPositions (bottomLeft, topLeft, bottomRight, topRight);
}

bool ThresholdAdjuster::notifyListener ()
{

    if (adjusterListener != NULL && !blocked) {
        GThreadLock lock;
        sendToListener();
    }

    return false;
}

void ThresholdAdjuster::setBgCurveProvider (ThresholdCurveProvider* provider)
{
    tSelector.setBgCurveProvider(provider);
}


void ThresholdAdjuster::setEnabled (bool enabled)
{

    tSelector.set_sensitive (enabled);
}

void ThresholdAdjuster::setEditedState (EditedState eState)
{

    if (editedState != eState) {
        if (editedCheckBox) {
            editedChange.block (true);
            editedCheckBox->set_active (eState == Edited);
            editedChange.block (false);
        }

        editedState = eState;
        refreshLabelStyle ();
    }
}

EditedState ThresholdAdjuster::getEditedState ()
{

    if (editedState != Irrelevant && editedCheckBox) {
        editedState = editedCheckBox->get_active () ? Edited : UnEdited;
    }

    return editedState;
}

void ThresholdAdjuster::showEditedCB ()
{

    if (!editedCheckBox) {
        editedCheckBox =  Gtk::manage(new Gtk::CheckButton ());
        hbox->pack_start (*editedCheckBox, Gtk::PACK_SHRINK, 2);
        hbox->reorder_child (*editedCheckBox, 0);
        editedChange = editedCheckBox->signal_toggled().connect( sigc::mem_fun(*this, &ThresholdAdjuster::editedToggled) );
    }
}

void ThresholdAdjuster::refreshLabelStyle ()
{

    /*  Glib::RefPtr<Gtk::Style> style = label->get_style ();
        Pango::FontDescription fd = style->get_font ();
        fd.set_weight (editedState==Edited ? Pango::WEIGHT_BOLD : Pango::WEIGHT_NORMAL);
        style->set_font (fd);
        label->set_style (style);
        label->queue_draw ();*/
}

void ThresholdAdjuster::editedToggled ()
{

    if (adjusterListener && !blocked) {
        sendToListener ();
    }
}

void ThresholdAdjuster::sendToListener ()
{
    if (tSelector.getPrecision() > 0) {
        // if precision is >0, then we assume that the listener is waiting for doubles
        rtengine::procparams::Threshold<double> t = tSelector.getPositions<double>();

        if (tSelector.isDouble()) {
            adjusterListener->adjusterChanged (this, t.value[0], t.value[1], t.value[2], t.value[3]);
            adjusterListener->adjusterChanged2 (this, t.value[0], t.value[1], t.value[2], t.value[3]);
        } else {
            adjusterListener->adjusterChanged (this, t.value[0], t.value[1]);
        }
    } else {
        // if precision is equal to 0, then we assume that the listener is waiting for integers
        rtengine::procparams::Threshold<int> t = tSelector.getPositions<int>();

        if (tSelector.isDouble()) {
            adjusterListener->adjusterChanged (this, t.value[0], t.value[1], t.value[2], t.value[3]);
            adjusterListener->adjusterChanged2 (this, t.value[0], t.value[1], t.value[2], t.value[3]);
        } else {
            adjusterListener->adjusterChanged (this, t.value[0], t.value[1]);
        }
    }
}

void ThresholdAdjuster::set_tooltip_markup(const Glib::ustring& markup)
{
    tSelector.set_tooltip_markup(markup);
}

void ThresholdAdjuster::set_tooltip_text(const Glib::ustring& text)
{
    tSelector.set_tooltip_text(text);
}

/* For better readability, this method create the history string of the parameter column,
 * so that the parameters list can be read in a more logical way (i.e. corresponding
 * to the startAtOne field)
 *
 * If separatedMode==true, the top slider is assumed to be the primary slider, then the bottom slider as the second one
 */
Glib::ustring ThresholdAdjuster::getHistoryString ()
{
    if (tSelector.isDouble()) {
        Glib::ustring bl, tl, br, tr;
        tSelector.getPositions(bl, tl, br, tr);
        return Glib::ustring::compose(tSelector.isStartAtOne() ? "%2, %1, %3, %4" : "%1, %2, %4, %3", bl, tl, br, tr);
    } else {
        Glib::ustring b, t;
        tSelector.getPositions(b, t);
        return Glib::ustring::compose(tSelector.isStartAtOne() || separatedMode ? "%2, %1" : "%1, %2", b, t);
    }
}
