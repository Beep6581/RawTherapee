/*
 *  This file is part of RawTherapee.
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
 *
 *  © 2010 Emil Martinec <ejmartin@uchicago.edu>
 */

#include "dirpyrequalizer.h"

using namespace rtengine;
using namespace rtengine::procparams;

DirPyrEqualizer::DirPyrEqualizer () : Gtk::VBox(), FoldableToolPanel(this) {

    set_border_width(4);

    enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
    enabled->set_active (true);
    pack_start(*enabled);
    enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrEqualizer::enabledToggled) );
    enabled->set_tooltip_markup (M("TP_SHARPENING_TOOLTIP"));

    Gtk::HSeparator *separator1 = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator1, Gtk::PACK_SHRINK, 2);
    
    Gtk::HBox * buttonBox1 = Gtk::manage (new Gtk::HBox());
    pack_start(*buttonBox1, Gtk::PACK_SHRINK, 2);

        Gtk::Button * lumacontrastMinusButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")));
        buttonBox1->pack_start(*lumacontrastMinusButton, Gtk::PACK_SHRINK, 2);
        lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumacontrastMinusPressed));

        Gtk::Button * lumaneutralButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")));
        buttonBox1->pack_start(*lumaneutralButton, Gtk::PACK_SHRINK, 2);
        lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumaneutralPressed));
        
        Gtk::Button * lumacontrastPlusButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")));
        buttonBox1->pack_start(*lumacontrastPlusButton, Gtk::PACK_SHRINK, 2);
        lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumacontrastPlusPressed));

    buttonBox1->show_all_children();

    Gtk::HSeparator *separator2 = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator2, Gtk::PACK_SHRINK, 2);

    for(int i = 0; i < 5; i++)
    {
        Glib::ustring ss;
        ss = Glib::ustring::format(i);
        if     (i == 0) ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        else if(i == 4) ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        multiplier[i] = Gtk::manage ( new Adjuster (ss, 0, 4, 0.01, 1.0) );
        multiplier[i]->setAdjusterListener(this);
        pack_start(*multiplier[i]);
    }

    Gtk::HSeparator *separator3 = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator3, Gtk::PACK_SHRINK, 2);

    threshold = Gtk::manage ( new Adjuster (M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 1, 0.01, 0.2) );
    threshold->setAdjusterListener(this);
    pack_start(*threshold);

    show_all_children ();

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}

DirPyrEqualizer::~DirPyrEqualizer () {

}

void DirPyrEqualizer::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {

        enabled->set_inconsistent (!pedited->dirpyrequalizer.enabled);

        for(int i = 0; i < 5; i++) {
            multiplier[i]->setEditedState (pedited->dirpyrequalizer.mult[i] ? Edited : UnEdited);
        }
        threshold->setEditedState (pedited->dirpyrequalizer.threshold ? Edited : UnEdited);
    }

    enaConn.block (true);
    enabled->set_active (pp->dirpyrequalizer.enabled);
    enaConn.block (false);
    lastEnabled = pp->dirpyrequalizer.enabled;
    
    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(pp->dirpyrequalizer.mult[i]);
    }
    threshold->setValue(pp->dirpyrequalizer.threshold);

    enableListener ();
}

void DirPyrEqualizer::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->dirpyrequalizer.enabled = enabled->get_active ();

    for (int i = 0; i < 5; i++) {
        pp->dirpyrequalizer.mult[i] = multiplier[i]->getValue();
    }
    pp->dirpyrequalizer.threshold = threshold->getValue();

    if (pedited) {

        pedited->dirpyrequalizer.enabled =  !enabled->get_inconsistent();

        for(int i = 0; i < 5; i++) {
            pedited->dirpyrequalizer.mult[i] = multiplier[i]->getEditedState();
        }
        pedited->dirpyrequalizer.threshold = threshold->getEditedState();
    }
}

void DirPyrEqualizer::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setDefault(defParams->dirpyrequalizer.mult[i]);
    }
    threshold->setDefault(defParams->dirpyrequalizer.threshold);
    
    if (pedited) {
        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(pedited->dirpyrequalizer.mult[i] ? Edited : UnEdited);
        }
        threshold->setDefaultEditedState(pedited->dirpyrequalizer.threshold ? Edited : UnEdited);
    }
    else {
        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(Irrelevant);
        }
        threshold->setDefaultEditedState(Irrelevant);
    }
}

void DirPyrEqualizer::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->showEditedCB();
    }
    threshold->showEditedCB();
}

void DirPyrEqualizer::adjusterChanged (Adjuster* a, double newval) {

    if (listener && enabled->get_active()) {
        if (a == threshold) {
            listener->panelChanged (EvDirPyrEqualizerThreshold,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(2), threshold->getValue()))
            );
        }
        else {
            listener->panelChanged (EvDirPyrEqualizer,
                         Glib::ustring::compose("%1, %2, %3, %4, %5",
                         Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[0]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[1]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[2]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[3]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[4]->getValue()))
            );
        }
    }    
}

void DirPyrEqualizer::enabledToggled () {

    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            enaConn.block (true);
            enabled->set_active (false);
            enaConn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);

        lastEnabled = enabled->get_active ();
    }

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_DISABLED"));
    }
}


void DirPyrEqualizer::lumaneutralPressed () {

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(1.0);
        adjusterChanged(multiplier[i], 1.0);
    }
}


void DirPyrEqualizer::lumacontrastPlusPressed () {

    for (int i = 0; i < 5; i++) {
        float inc = 0.05 * (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}


void DirPyrEqualizer::lumacontrastMinusPressed () {

    for (int i = 0; i < 5; i++) {
        float inc = -0.05 * (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}

void DirPyrEqualizer::setAdjusterBehavior (bool multiplieradd, bool thresholdadd) {

	for (int i=0; i<5; i++)
		multiplier[i]->setAddMode(multiplieradd);
	threshold->setAddMode(thresholdadd);
}

void DirPyrEqualizer::trimValues (rtengine::procparams::ProcParams* pp) {

	for (int i=0; i<5; i++)
		multiplier[i]->trimValue(pp->dirpyrequalizer.mult[i]);
	threshold->trimValue(pp->dirpyrequalizer.threshold);
}
