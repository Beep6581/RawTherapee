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
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 */

#include <equalizer.h>

using namespace rtengine;
using namespace rtengine::procparams;

Equalizer::Equalizer () : ToolPanel() {

    enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
    enabled->set_active (true);
    pack_start(*enabled);
    enabled->show();
    enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Equalizer::enabled_toggled) );

    Gtk::HSeparator *separator = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator, Gtk::PACK_SHRINK, 2);

    for(int i = 0; i < 8; i++)
    {
        std::stringstream ss;
        ss << i;
        if(i == 0)
            ss << " (" << M("TP_EQUALIZER_FINEST") << ")";
        if(i == 7)
            ss << " (" << M("TP_EQUALIZER_LARGEST") << ")";
        
        correction[i] = new Adjuster (ss.str(), -100, 100, 1, 0);
        correction[i]->setAdjusterListener(this);
        pack_start(*correction[i]);
    }

    show_all_children ();
}

Equalizer::~Equalizer () {

}

void Equalizer::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {

        enabled->set_inconsistent (!pedited->equalizer.enabled);

        for(int i = 0; i < 8; i++) {
            correction[i]->setEditedState (pedited->equalizer.c[i] ? Edited : UnEdited);
        }
    }

    enaConn.block (true);
    enabled->set_active (pp->equalizer.enabled);
    enaConn.block (false);
    lastEnabled = pp->equalizer.enabled;
    
    for (int i = 0; i < 8; i++) {
        correction[i]->setValue(pp->equalizer.c[i]);
    }

    enableListener ();
}

void Equalizer::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->equalizer.enabled = enabled->get_active ();

    for (int i = 0; i < 8; i++) {
        pp->equalizer.c[i] = (int) correction[i]->getValue();
    }

    if (pedited) {

        pedited->equalizer.enabled =  !enabled->get_inconsistent();

        for(int i = 0; i < 8; i++) {
            pedited->equalizer.c[i] = correction[i]->getEditedState();
        }
    }
}

void Equalizer::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    for (int i = 0; i < 8; i++) {
        correction[i]->setDefault(defParams->equalizer.c[i]);
    }
    
    if (pedited) {
        for (int i = 0; i < 8; i++) {
            correction[i]->setDefaultEditedState(pedited->equalizer.c[i] ? Edited : UnEdited);
        }
    }
    else {
        for (int i = 0; i < 8; i++) {
            correction[i]->setDefaultEditedState(Irrelevant);
        }
    }
}

void Equalizer::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    
    for (int i = 0; i < 8; i++) {
        correction[i]->showEditedCB();
    }
}

void Equalizer::adjusterChanged (Adjuster* a, double newval) {
    
    if (listener && enabled->get_active()) {
        std::stringstream ss;
        ss << "(";
        int i;
        for (i = 0; i < 8; i++) {
            if (i > 0) {
                ss << ", ";
            }
            ss << static_cast<int>(correction[i]->getValue());
        }
        ss << ")";
        listener->panelChanged (EvEqualizer, ss.str());
    }    
}

void Equalizer::enabled_toggled () {

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
            listener->panelChanged (EvEqlEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvEqlEnabled, M("GENERAL_DISABLED"));
    }
}


