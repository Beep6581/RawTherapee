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

#include <hsvequalizer.h>

using namespace rtengine;
using namespace rtengine::procparams;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HSVEqualizer::HSVEqualizer () : Gtk::VBox(), FoldableToolPanel(this) {
	
	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
    enabled->set_active (true);
    pack_start(*enabled);
    enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &HSVEqualizer::enabledToggled) );	
	
	Gtk::HSeparator *hsvsepa = Gtk::manage (new  Gtk::HSeparator());
	pack_start(*hsvsepa, Gtk::PACK_SHRINK, 2);
	hsvsepa->show ();
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
	hb->set_border_width (4);
	hb->show ();
	Gtk::Label* hsvselect = Gtk::manage (new Gtk::Label (M("TP_HSVEQUALIZER_CHANNEL")+":"));
	hsvselect->show ();
	hsvchannel = Gtk::manage (new Gtk::ComboBoxText ());
	hsvchannel->append_text (M("TP_HSVEQUALIZER_SAT"));
	hsvchannel->append_text (M("TP_HSVEQUALIZER_VAL"));
	hsvchannel->append_text (M("TP_HSVEQUALIZER_HUE"));
	hsvchannel->show ();
	hb->pack_start(*hsvselect, Gtk::PACK_SHRINK, 4);
	hb->pack_start(*hsvchannel); 
	
	Gtk::Button * neutralButton = Gtk::manage (new Gtk::Button(M("TP_HSVEQUALIZER_NEUTRAL")));
	hb->pack_start(*neutralButton, Gtk::PACK_SHRINK, 2);
	neutralPressedConn = neutralButton->signal_pressed().connect( sigc::mem_fun(*this, &HSVEqualizer::neutralPressed));
	
	
	pack_start (*hb);
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*
	 Gtk::HBox * buttonBox17 = Gtk::manage (new Gtk::HBox());
	 pack_start(*buttonBox17, Gtk::PACK_SHRINK, 2);
	 
	 Gtk::Button * neutralButton = Gtk::manage (new Gtk::Button(M("TP_HSVEQUALIZER_NEUTRAL")));
	 buttonBox17->pack_start(*neutralButton, Gtk::PACK_SHRINK, 2);
	 neutralPressedConn = neutralButton->signal_pressed().connect( sigc::mem_fun(*this, &HSVEqualizer::neutralPressed));
	 
	 buttonBox17->show();
	 */
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	Gtk::HSeparator *hsvsepb = Gtk::manage (new  Gtk::HSeparator());
	pack_start(*hsvsepb, Gtk::PACK_SHRINK, 2);
	hsvsepb->show ();
	
	satbox = new Gtk::VBox ();
	sat[0] = new Adjuster (M("TP_HSVEQUALIZER1"), -100, 100, 0.1, 0);
	sat[1] = new Adjuster (M("TP_HSVEQUALIZER2"), -100, 100, 0.1, 0);
	sat[2] = new Adjuster (M("TP_HSVEQUALIZER3"), -100, 100, 0.1, 0);
	sat[3] = new Adjuster (M("TP_HSVEQUALIZER4"), -100, 100, 0.1, 0);
	sat[4] = new Adjuster (M("TP_HSVEQUALIZER5"), -100, 100, 0.1, 0);
	sat[5] = new Adjuster (M("TP_HSVEQUALIZER6"), -100, 100, 0.1, 0);
	sat[6] = new Adjuster (M("TP_HSVEQUALIZER7"), -100, 100, 0.1, 0);
	sat[7] = new Adjuster (M("TP_HSVEQUALIZER8"), -100, 100, 0.1, 0);
    for(int i = 0; i < 8; i++)
    {
        sat[i]->setAdjusterListener(this);
        satbox->pack_start(*sat[i]);
    }
	
    //show_all_children ();
	satbox->show ();
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	valbox = new Gtk::VBox ();
	val[0] = new Adjuster (M("TP_HSVEQUALIZER1"), -100, 100, 0.1, 0);
	val[1] = new Adjuster (M("TP_HSVEQUALIZER2"), -100, 100, 0.1, 0);
	val[2] = new Adjuster (M("TP_HSVEQUALIZER3"), -100, 100, 0.1, 0);
	val[3] = new Adjuster (M("TP_HSVEQUALIZER4"), -100, 100, 0.1, 0);
	val[4] = new Adjuster (M("TP_HSVEQUALIZER5"), -100, 100, 0.1, 0);
	val[5] = new Adjuster (M("TP_HSVEQUALIZER6"), -100, 100, 0.1, 0);
	val[6] = new Adjuster (M("TP_HSVEQUALIZER7"), -100, 100, 0.1, 0);
	val[7] = new Adjuster (M("TP_HSVEQUALIZER8"), -100, 100, 0.1, 0);
    for(int i = 0; i < 8; i++)
    {
        val[i]->setAdjusterListener(this);
        valbox->pack_start(*val[i]);
    }
	
    //show_all_children ();
	valbox->show ();	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	huebox = new Gtk::VBox ();
	
	
	hue[0] = new Adjuster (M("TP_HSVEQUALIZER1"), -100, 100, 0.1, 0);
	hue[1] = new Adjuster (M("TP_HSVEQUALIZER2"), -100, 100, 0.1, 0);
	hue[2] = new Adjuster (M("TP_HSVEQUALIZER3"), -100, 100, 0.1, 0);
	hue[3] = new Adjuster (M("TP_HSVEQUALIZER4"), -100, 100, 0.1, 0);
	hue[4] = new Adjuster (M("TP_HSVEQUALIZER5"), -100, 100, 0.1, 0);
	hue[5] = new Adjuster (M("TP_HSVEQUALIZER6"), -100, 100, 0.1, 0);
	hue[6] = new Adjuster (M("TP_HSVEQUALIZER7"), -100, 100, 0.1, 0);
	hue[7] = new Adjuster (M("TP_HSVEQUALIZER8"), -100, 100, 0.1, 0);
    for(int i = 0; i < 8; i++)
    {
        hue[i]->setAdjusterListener(this);
        huebox->pack_start(*hue[i]);
    }
	
    //huebox->show_all_children ();
	huebox->show ();
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    huebox->reference ();
    valbox->reference ();
    satbox->reference ();
	
	//enaConn 	= enabled->signal_toggled().connect( sigc::mem_fun(*this, &HSVEqualizer::enabled_toggled) );
	hsvchannel->signal_changed().connect( sigc::mem_fun(*this, &HSVEqualizer::hsvchannelChanged) );
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HSVEqualizer::~HSVEqualizer () {
	for (int i=0;i<8;i++)
	{
		delete hue[i];
		delete val[i];
		delete sat[i];
	}
}


void HSVEqualizer::read (const ProcParams* pp, const ParamsEdited* pedited) {
	
    disableListener ();
	
    if (pedited) {
		for (int i=0; i<8; i++) {
			sat[i]->setEditedState 		(pedited->hsvequalizer.sat[i] ? Edited : UnEdited);
			val[i]->setEditedState 		(pedited->hsvequalizer.val[i] ? Edited : UnEdited);
			hue[i]->setEditedState 		(pedited->hsvequalizer.hue[i] ? Edited : UnEdited);
		}
        enabled->set_inconsistent 		(!pedited->hsvequalizer.enabled);
    }
	
    enaConn.block (true);
    enabled->set_active (pp->hsvequalizer.enabled);
    enaConn.block (false);   
    lastEnabled = pp->hsvequalizer.enabled;
	
	for (int i=0; i<8; i++) {
		sat[i]->setValue        (pp->hsvequalizer.sat[i]);
		val[i]->setValue        (pp->hsvequalizer.val[i]);
		hue[i]->setValue        (pp->hsvequalizer.hue[i]);
	}
	
	if (pedited && !pedited->hsvequalizer.hsvchannel)
		hsvchannel->set_active (3);
	else if (pp->hsvequalizer.hsvchannel=="sat")
		hsvchannel->set_active (0);
	else if (pp->hsvequalizer.hsvchannel=="val")
		hsvchannel->set_active (1);
	else if (pp->hsvequalizer.hsvchannel=="hue")
		hsvchannel->set_active (2);
	
    enableListener ();
}

void HSVEqualizer::write (ProcParams* pp, ParamsEdited* pedited) {
	
	pp->hsvequalizer.enabled        = enabled->get_active ();
	for (int i=0; i<8; i++) {
		pp->hsvequalizer.sat[i]			= sat[i]->getValue();
		pp->hsvequalizer.val[i]			= val[i]->getValue();
		pp->hsvequalizer.hue[i]			= hue[i]->getValue();
	}
    
    if (hsvchannel->get_active_row_number()==0)
        pp->hsvequalizer.hsvchannel = "sat";
    else if (hsvchannel->get_active_row_number()==1)
        pp->hsvequalizer.hsvchannel = "val";
	else if (hsvchannel->get_active_row_number()==2)
        pp->hsvequalizer.hsvchannel = "hue";
	
    if (pedited) {
		pedited->hsvequalizer.enabled =  !enabled->get_inconsistent();//from dirpyreq
		for (int i=0; i<8; i++) {
			pedited->hsvequalizer.sat[i] 			= sat[i]->getEditedState ();
			pedited->hsvequalizer.val[i] 			= val[i]->getEditedState ();
			pedited->hsvequalizer.hue[i] 			= hue[i]->getEditedState ();
		}
    }
}

void HSVEqualizer::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {
	
    for (int i = 0; i < 8; i++) {
        sat[i]->setDefault(defParams->hsvequalizer.sat[i]);
		val[i]->setDefault(defParams->hsvequalizer.val[i]);
        hue[i]->setDefault(defParams->hsvequalizer.hue[i]);
    }
    
    if (pedited) {
        for (int i = 0; i < 8; i++) {
            sat[i]->setDefaultEditedState(pedited->hsvequalizer.sat[i] ? Edited : UnEdited);
			val[i]->setDefaultEditedState(pedited->hsvequalizer.val[i] ? Edited : UnEdited);
            hue[i]->setDefaultEditedState(pedited->hsvequalizer.hue[i] ? Edited : UnEdited);
        }
    }
    else {
        for (int i = 0; i < 8; i++) {
            sat[i]->setDefaultEditedState(Irrelevant);
			val[i]->setDefaultEditedState(Irrelevant);
            hue[i]->setDefaultEditedState(Irrelevant);
        }
    }
}

void HSVEqualizer::adjusterChanged (Adjuster* a, double newval) {
	
	if (listener && enabled->get_active()) {
        std::stringstream ss;
        ss << "(";
        int i;		
		if (hsvchannel->get_active_row_number()==0) {
			for (i = 0; i < 8; i++) {
				if (i > 0) {
					ss << ", ";
				}
				if (i == 4) {
					ss << "\n";
				}
				ss << static_cast<int>(sat[i]->getValue());
			}
			ss << ")";
			listener->panelChanged (EvHSVEqualizerS, ss.str());
		} 
		else if (hsvchannel->get_active_row_number()==1) {
			for (i = 0; i < 8; i++) {
				if (i > 0) {
					ss << ", ";
				}
				if (i == 4) {
					ss << "\n";
				}
				ss << static_cast<int>(val[i]->getValue());
			}
			ss << ")";
			listener->panelChanged (EvHSVEqualizerV, ss.str());
		}
		else if (hsvchannel->get_active_row_number()==2) {
			for (i = 0; i < 8; i++) {
				if (i > 0) {
					ss << ", ";
				}
				if (i == 4) {
					ss << "\n";
				}
				ss << static_cast<int>(hue[i]->getValue());
			}
			ss << ")";
			listener->panelChanged (EvHSVEqualizerH, ss.str());
		}
		
		//listener->panelChanged (EvHSVEqualizer, ss.str());
	}
}

void HSVEqualizer::enabledToggled () {
	
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
			listener->panelChanged (EvHSVEqEnabled, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvHSVEqEnabled, M("GENERAL_DISABLED"));
	}
}


void HSVEqualizer::hsvchannelChanged () {
	
	removeIfThere (this, satbox, false);
	removeIfThere (this, valbox, false);
	removeIfThere (this, huebox, false);
	
	if (hsvchannel->get_active_row_number()==0) 
		pack_start (*satbox);   
	else if (hsvchannel->get_active_row_number()==1) 
		pack_start (*valbox);
	else if (hsvchannel->get_active_row_number()==2) 
		pack_start (*huebox);
	
}

void HSVEqualizer::setBatchMode (bool batchMode) {
	
	ToolPanel::setBatchMode (batchMode);
	
	for (int i = 0; i < 8; i++) {
		sat[i]->showEditedCB();
		val[i]->showEditedCB();
		hue[i]->showEditedCB();
	}
	
	hsvchannel->append_text (M("GENERAL_UNCHANGED"));
}

void HSVEqualizer::neutralPressed () {
	if (hsvchannel->get_active_row_number()==0)
		for (int i = 0; i < 8; i++) {
			sat[i]->setValue(0);
			adjusterChanged(sat[i], 0);
		}
	else if (hsvchannel->get_active_row_number()==1)
		for (int i = 0; i < 8; i++) {
			val[i]->setValue(0);
			adjusterChanged(val[i], 0);
		}
	else if (hsvchannel->get_active_row_number()==2)
		for (int i = 0; i < 8; i++) {
			hue[i]->setValue(0);
			adjusterChanged(hue[i], 0);
		}
	
	
}
