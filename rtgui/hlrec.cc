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
#include <hlrec.h>
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

HLRecovery::HLRecovery () : ToolPanel() {

    enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLE")));
    enabled->set_active (false);
    pack_start (*enabled);

	method = Gtk::manage (new Gtk::ComboBoxText ());
	method->append_text (M("TP_HLREC_LUMINANCE"));
	method->append_text (M("TP_HLREC_CIELAB"));
	method->append_text (M("TP_HLREC_COLOR"));
	method->set_active (0);
	Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
	Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_HLREC_METHOD")));
	hb->pack_start (*lab, Gtk::PACK_SHRINK, 4);
	hb->pack_start (*method);
	pack_start (*hb);

    enaconn  = enabled->signal_toggled().connect( sigc::mem_fun(*this, &HLRecovery::enabledChanged) );
	methconn = method->signal_changed().connect ( sigc::mem_fun(*this, &HLRecovery::methodChanged) );

	show_all ();
}

void HLRecovery::read (const ProcParams* pp, const ParamsEdited* pedited) {
    
    disableListener ();


    if (pedited)
        enabled->set_inconsistent (!pedited->hlrecovery.enabled);
    enaconn.block (true);
    enabled->set_active  (pp->hlrecovery.enabled);
    enaconn.block (false);
    
    if (pedited && !pedited->hlrecovery.method) 
        method->set_active (3);
	else if (pp->hlrecovery.method=="Luminance")
		method->set_active (0);
	else if (pp->hlrecovery.method=="CIELab blending")
		method->set_active (1);
	else if (pp->hlrecovery.method=="Color")
		method->set_active (2);

    lastEnabled = pp->hlrecovery.enabled;

    enableListener ();
}

void HLRecovery::write (ProcParams* pp, ParamsEdited* pedited) {

    if (pedited) {
        pedited->hlrecovery.method    = method->get_active_row_number()!=3;
        pedited->hlrecovery.enabled   = !enabled->get_inconsistent();
    }

    pp->hlrecovery.enabled = enabled->get_active();
	if (method->get_active_row_number()==0)
		pp->hlrecovery.method = "Luminance";
	else if (method->get_active_row_number()==1)
		pp->hlrecovery.method = "CIELab blending";
	else if (method->get_active_row_number()==2)
		pp->hlrecovery.method = "Color";
}

void HLRecovery::enabledChanged () {

    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            enaconn.block (true);
            enabled->set_active (false);
            enaconn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);

        lastEnabled = enabled->get_active ();
    }

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvHREnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvHREnabled, M("GENERAL_DISABLED"));
    }  
}

void HLRecovery::methodChanged () {

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvHRMethod, method->get_active_text ());
    }  
}

void HLRecovery::setRaw (bool raw) {
    
    disableListener ();
    set_sensitive (raw);
    enableListener ();
}

void HLRecovery::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    method->append_text ("(Unchanged)");
}
