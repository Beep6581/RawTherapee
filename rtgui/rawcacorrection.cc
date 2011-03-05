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
#include <rawcacorrection.h>
#include <guiutils.h>
#include <safegtk.h>
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

RAWCACorr::RAWCACorr () : Gtk::VBox(), FoldableToolPanel(this)
{
	caAutocorrect = Gtk::manage(new Gtk::CheckButton((M("TP_RAWCACORR_AUTO"))));
	caRed = Gtk::manage(new Adjuster (M("TP_RAWCACORR_CARED"),-4.0,4.0,0.1,0));
	caRed->setAdjusterListener (this);
	if (caRed->delay < 1000) caRed->delay = 1000;
	caRed->show();
	caBlue = Gtk::manage(new Adjuster (M("TP_RAWCACORR_CABLUE"),-4.0,4.0,0.1,0));
	caBlue->setAdjusterListener (this);
	if (caBlue->delay < 1000) caBlue->delay = 1000;
	caBlue->show();

	pack_start( *caAutocorrect, Gtk::PACK_SHRINK, 4);
	pack_start( *caRed, Gtk::PACK_SHRINK, 4);
	pack_start( *caBlue, Gtk::PACK_SHRINK, 4);

	caacsconn = caAutocorrect->signal_toggled().connect ( sigc::mem_fun(*this, &RAWCACorr::caCorrectionChanged), true);
}

void RAWCACorr::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();
	caacsconn.block (true);

	if(pedited ){
		caAutocorrect->set_inconsistent(!pedited->raw.caCorrection);
		caRed->setEditedState( pedited->raw.caRed ? Edited : UnEdited );
		caBlue->setEditedState( pedited->raw.caBlue ? Edited : UnEdited );
	}

	lastCA  = pp->raw.ca_autocorrect;

	caAutocorrect->set_active(pp->raw.ca_autocorrect);
	caRed->setValue (pp->raw.cared);
	caBlue->setValue (pp->raw.cablue);

	caacsconn.block (false);
	enableListener ();
}

void RAWCACorr::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.ca_autocorrect = caAutocorrect->get_active();
	pp->raw.cared = caRed->getValue();
	pp->raw.cablue = caBlue->getValue();

	if (pedited) {
		pedited->raw.caCorrection = !caAutocorrect->get_inconsistent();
		pedited->raw.caRed = caRed->getEditedState ();
		pedited->raw.caBlue = caBlue->getEditedState ();
	}

}

void RAWCACorr::adjusterChanged (Adjuster* a, double newval)
{
	if (listener) {

		Glib::ustring value = a->getTextValue();

		if (a == caRed)
			listener->panelChanged (EvPreProcessCARed,  value );
		else if (a == caBlue)
			listener->panelChanged (EvPreProcessCABlue,  value );
	}
}

void RAWCACorr::setBatchMode(bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	caRed->showEditedCB ();
	caBlue->showEditedCB ();
}

void RAWCACorr::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	caRed->setDefault( defParams->raw.cared);
	caBlue->setDefault( defParams->raw.cablue);

	if (pedited) {
		caRed->setDefaultEditedState( pedited->raw.caRed ? Edited : UnEdited);
		caBlue->setDefaultEditedState( pedited->raw.caBlue ? Edited : UnEdited);
	} else {
		caRed->setDefaultEditedState( Irrelevant );
		caBlue->setDefaultEditedState( Irrelevant );
	}
}

void RAWCACorr::caCorrectionChanged()
{
    if (batchMode) {
        if (caAutocorrect->get_inconsistent()) {
        	caAutocorrect->set_inconsistent (false);
        	caacsconn.block (true);
            caAutocorrect->set_active (false);
            caacsconn.block (false);
        }
        else if (lastCA)
        	caAutocorrect->set_inconsistent (true);

        lastCA = caAutocorrect->get_active ();

    }
    else {
    	// For non batch mode, we disable the red and blue slider if caAutocorrect is true
    	if (caAutocorrect->get_active ()) {
    		caRed->set_sensitive(false);
    		caBlue->set_sensitive(false);
    	}
    	else {
    		caRed->set_sensitive(true);
    		caBlue->set_sensitive(true);
    	}
    }
    if (listener)
        listener->panelChanged (EvPreProcessAutoCA, caAutocorrect->get_active()?M("GENERAL_ENABLED"):M("GENERAL_DISABLED"));
}
