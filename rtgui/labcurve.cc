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
#include "labcurve.h"
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

LCurve::LCurve () : Gtk::VBox(), FoldableToolPanel(this) {

	brightness = Gtk::manage (new Adjuster (M("TP_LABCURVE_BRIGHTNESS"), -100, 100, 1, 0));
	contrast   = Gtk::manage (new Adjuster (M("TP_LABCURVE_CONTRAST"), -100, 100, 1, 0));
	saturation   = Gtk::manage (new Adjuster (M("TP_LABCURVE_SATURATION"), -100, 100, 1, 0));

	pack_start (*brightness);
	brightness->show ();

	pack_start (*contrast);
	contrast->show ();

	pack_start (*saturation);
	saturation->show ();
		
	brightness->setAdjusterListener (this);
	contrast->setAdjusterListener (this);
	saturation->setAdjusterListener (this);
	
	//%%%%%%%%%%%%%%%%%%
	pack_start (*Gtk::manage (new  Gtk::HSeparator()));
	
	avoidclip = Gtk::manage (new Gtk::CheckButton (M("TP_LABCURVE_AVOIDCOLORCLIP")));
	
	pack_start (*avoidclip);
	pack_start (*Gtk::manage (new  Gtk::HSeparator()));
	
	enablelimiter = Gtk::manage (new Gtk::CheckButton (M("TP_LABCURVE_ENABLESATLIMITER")));
	pack_start (*enablelimiter);
	
	saturationlimiter = Gtk::manage ( new Adjuster (M("TP_LABCURVE_SATLIMIT"), 0, 100, 1.0, 40) );
	pack_start (*saturationlimiter);
	saturationlimiter->show ();
	saturationlimiter->reference ();  
	
	//saturation->setAdjusterListener (this);
	saturationlimiter->setAdjusterListener (this);
	acconn = avoidclip->signal_toggled().connect( sigc::mem_fun(*this, &LCurve::avoidclip_toggled) );
	elconn = enablelimiter->signal_toggled().connect( sigc::mem_fun(*this, &LCurve::enablelimiter_toggled) );
	//%%%%%%%%%%%%%%%%%%%

	Gtk::HSeparator *hsep3 = Gtk::manage (new  Gtk::HSeparator());
	hsep3->show ();
	pack_start (*hsep3);

	curveEditorG = new CurveEditorGroup (options.lastLabCurvesDir);
	curveEditorG->setCurveListener (this);
	curveEditorG->setColorProvider (this);

	lshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "L"));
	ashape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "a"));
	bshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "b"));

	// This will add the reset button at the end of the curveType buttons
	curveEditorG->curveListComplete();

	pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);
	
}

LCurve::~LCurve () {
	delete curveEditorG;
}

void LCurve::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        brightness->setEditedState (pedited->labCurve.brightness ? Edited : UnEdited);
        contrast->setEditedState (pedited->labCurve.contrast ? Edited : UnEdited);
		saturation->setEditedState (pedited->labCurve.saturation ? Edited : UnEdited);
		
		//%%%%%%%%%%%%%%%%%%%%%%
		saturationlimiter->setEditedState (pedited->labCurve.saturationlimit ? Edited : UnEdited);
        avoidclip->set_inconsistent (!pedited->labCurve.avoidclip);
        enablelimiter->set_inconsistent (!pedited->labCurve.enable_saturationlimiter);
		//%%%%%%%%%%%%%%%%%%%%%%

        lshape->setUnChanged (!pedited->labCurve.lcurve);
		ashape->setUnChanged (!pedited->labCurve.acurve);
        bshape->setUnChanged (!pedited->labCurve.bcurve);

    }

    brightness->setValue    (pp->labCurve.brightness);
    contrast->setValue      (pp->labCurve.contrast);
	saturation->setValue      (pp->labCurve.saturation);
	
	//%%%%%%%%%%%%%%%%%%%%%%
	saturationlimiter->setValue (pp->labCurve.saturationlimit);
    acconn.block (true);
    avoidclip->set_active (pp->labCurve.avoidclip);
    acconn.block (false);
    elconn.block (true);
    enablelimiter->set_active (pp->labCurve.enable_saturationlimiter);
    elconn.block (false);
	
    //removeIfThere (this, saturationlimiter, false);
    // if (enablelimiter->get_active () || enablelimiter->get_inconsistent())
    //    pack_start (*saturationlimiter);
	
    lastACVal = pp->labCurve.avoidclip;
    lastELVal = pp->labCurve.enable_saturationlimiter;
	//%%%%%%%%%%%%%%%%%%%%%%

    lshape->setCurve         (pp->labCurve.lcurve);
	ashape->setCurve         (pp->labCurve.acurve);
    bshape->setCurve         (pp->labCurve.bcurve);

    // Open up the first curve if selected
    bool active = lshape->openIfNonlinear();
    if (!active) ashape->openIfNonlinear();
    if (!active) bshape->openIfNonlinear();

    enableListener ();
}

void LCurve::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->labCurve.brightness    = brightness->getValue ();
    pp->labCurve.contrast      = (int)contrast->getValue ();
	pp->labCurve.saturation      = (int)saturation->getValue ();
	
	//%%%%%%%%%%%%%%%%%%%%%%
	pp->labCurve.avoidclip = avoidclip->get_active ();
    pp->labCurve.enable_saturationlimiter = enablelimiter->get_active ();
    pp->labCurve.saturationlimit = saturationlimiter->getValue ();
	//%%%%%%%%%%%%%%%%%%%%%%

    pp->labCurve.lcurve         = lshape->getCurve ();
	pp->labCurve.acurve         = ashape->getCurve ();
    pp->labCurve.bcurve         = bshape->getCurve ();

    if (pedited) {
        pedited->labCurve.brightness = brightness->getEditedState ();
        pedited->labCurve.contrast = contrast->getEditedState ();
		pedited->labCurve.saturation = saturation->getEditedState ();
		
		//%%%%%%%%%%%%%%%%%%%%%%
		pedited->labCurve.avoidclip = !avoidclip->get_inconsistent();
        pedited->labCurve.enable_saturationlimiter = !enablelimiter->get_inconsistent();
        pedited->labCurve.saturationlimit = saturationlimiter->getEditedState ();
		//%%%%%%%%%%%%%%%%%%%%%%

        pedited->labCurve.lcurve    = !lshape->isUnChanged ();
		pedited->labCurve.acurve    = !ashape->isUnChanged ();
        pedited->labCurve.bcurve    = !bshape->isUnChanged ();
    }
}

void LCurve::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    brightness->setDefault (defParams->labCurve.brightness);
    contrast->setDefault (defParams->labCurve.contrast);
	saturation->setDefault (defParams->labCurve.saturation);
    saturationlimiter->setDefault (defParams->labCurve.saturationlimit);

    if (pedited) {
        brightness->setDefaultEditedState (pedited->labCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->labCurve.contrast ? Edited : UnEdited);
		saturation->setDefaultEditedState (pedited->labCurve.saturation ? Edited : UnEdited);
        saturationlimiter->setDefaultEditedState (pedited->labCurve.saturationlimit ? Edited : UnEdited);

    }
    else {
        brightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
		saturation->setDefaultEditedState (Irrelevant);
		saturationlimiter->setDefaultEditedState (Irrelevant);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%
//Clipping control changed
void LCurve::avoidclip_toggled () {
	
    if (batchMode) {
        if (avoidclip->get_inconsistent()) {
            avoidclip->set_inconsistent (false);
            acconn.block (true);
            avoidclip->set_active (false);
            acconn.block (false);
        }
        else if (lastACVal)
            avoidclip->set_inconsistent (true);
		
        lastACVal = avoidclip->get_active ();
    }
	
    if (listener) {
        if (avoidclip->get_active ()) 
            listener->panelChanged (EvLAvoidClip, M("GENERAL_ENABLED"));
        else            
            listener->panelChanged (EvLAvoidClip, M("GENERAL_DISABLED"));
    }
}

void LCurve::enablelimiter_toggled () {
	
    if (batchMode) {
        if (enablelimiter->get_inconsistent()) {
            enablelimiter->set_inconsistent (false);
            elconn.block (true);
            enablelimiter->set_active (false);
            elconn.block (false);
        }
        else if (lastELVal)
            enablelimiter->set_inconsistent (true);
		
        lastELVal = enablelimiter->get_active ();
    }
	
    //removeIfThere (this, saturationlimiter, false);
    //if (enablelimiter->get_active () || enablelimiter->get_inconsistent())
    //    pack_start (*saturationlimiter);
	
    if (listener) {
        if (enablelimiter->get_active ()) 
            listener->panelChanged (EvLSatLimiter, M("GENERAL_ENABLED"));
        else            
            listener->panelChanged (EvLSatLimiter, M("GENERAL_DISABLED"));
    }
}
//%%%%%%%%%%%%%%%%%%%%%%


/*
 * Curve listener
 *
 * If more than one curve has been added, the curve listener is automatically
 * set to 'multi=true', and send a pointer of the modified curve in a parameter
 */
void LCurve::curveChanged (CurveEditor* ce) {

    if (listener) {
    	if (ce == lshape)
        	listener->panelChanged (EvLLCurve, M("HISTORY_CUSTOMCURVE"));
    	if (ce == ashape)
        	listener->panelChanged (EvLaCurve, M("HISTORY_CUSTOMCURVE"));
    	if (ce == bshape)
        	listener->panelChanged (EvLbCurve, M("HISTORY_CUSTOMCURVE"));
	}
}

void LCurve::adjusterChanged (Adjuster* a, double newval) {

    if (!listener)
        return;

    Glib::ustring costr;
    if (a==brightness)
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
	else if (a==saturationlimiter)
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(1), a->getValue());
    else
        costr = Glib::ustring::format ((int)a->getValue());

    if (a==brightness)
        listener->panelChanged (EvLBrightness, costr);
    else if (a==contrast)
        listener->panelChanged (EvLContrast, costr);
	else if (a==saturation)
        listener->panelChanged (EvLSaturation, costr);
	else if (a==saturationlimiter)
        listener->panelChanged (EvLSatLimit, costr);
}

void LCurve::colorForValue (double valX, double valY) {

	CurveEditor* ce = curveEditorG->getDisplayedCurve();

	if (ce == lshape) {         // L = f(L)
		red = (double)valY;
		green = (double)valY;
		blue = (double)valY;
	}
	else if (ce == ashape) {    // a = f(a)
		// TODO: To be implemented
		red = (double)valY;
		green = (double)valY;
		blue = (double)valY;
	}
	else if (ce == bshape) {    // b = f(b)
		red = (double)valY;
		green = (double)valY;
		blue = (double)valY;
	}
	else {
		printf("Error: no curve displayed!\n");
	}

}

void LCurve::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    brightness->showEditedCB ();
    contrast->showEditedCB ();
	saturation->showEditedCB ();
	saturationlimiter->showEditedCB ();

    curveEditorG->setBatchMode (batchMode);
}


void LCurve::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma){

    lshape->updateBackgroundHistogram (histLCurve);
}

void LCurve::setAdjusterBehavior (bool bradd, bool contradd, bool satadd) {

	brightness->setAddMode(bradd);
	contrast->setAddMode(contradd);
	saturation->setAddMode(satadd);
}

void LCurve::trimValues (rtengine::procparams::ProcParams* pp) {

	brightness->trimValue(pp->labCurve.brightness);
	contrast->trimValue(pp->labCurve.contrast);
	saturation->trimValue(pp->labCurve.saturation);
}
