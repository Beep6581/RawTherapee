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
#include <labcurve.h>
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

LCurve::LCurve () : Gtk::VBox(), FoldableToolPanel(this), brAdd(false), contrAdd(false), satAdd(false) {

	brightness = Gtk::manage (new Adjuster (M("TP_LABCURVE_BRIGHTNESS"), -100, 100, 0.01, 0));
	contrast   = Gtk::manage (new Adjuster (M("TP_LABCURVE_CONTRAST"), -100, 100, 1, 0));
	saturation   = Gtk::manage (new Adjuster (M("TP_LABCURVE_SATURATION"), -100, 100, 1, 0));

	pack_start (*brightness);
	brightness->show ();

	pack_start (*contrast);
	contrast->show ();

	pack_start (*saturation);
	saturation->show ();

	Gtk::HSeparator *hsep3 = Gtk::manage (new  Gtk::HSeparator());
	hsep3->show ();
	pack_start (*hsep3);

	curveEditorG = new CurveEditorGroup ();
	curveEditorG->setCurveListener (this);

	lshape = curveEditorG->addCurve("L");
	ashape = curveEditorG->addCurve("a");
	bshape = curveEditorG->addCurve("b");

	// This will add the reset button at the end of the curveType buttons
	curveEditorG->curveListComplete();

	pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);

	//curveEditorG->show();

	brightness->setAdjusterListener (this);
	contrast->setAdjusterListener (this);
	saturation->setAdjusterListener (this);

	//channel->signal_changed().connect( sigc::mem_fun(*this, &LCurve::channel) );


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

        lshape->setUnChanged (!pedited->labCurve.lcurve);
		ashape->setUnChanged (!pedited->labCurve.acurve);
        bshape->setUnChanged (!pedited->labCurve.bcurve);

    }

    brightness->setValue    (pp->labCurve.brightness);
    contrast->setValue      (pp->labCurve.contrast);
	saturation->setValue      (pp->labCurve.saturation);

    lshape->setCurve         (pp->labCurve.lcurve);
	ashape->setCurve         (pp->labCurve.acurve);
    bshape->setCurve         (pp->labCurve.bcurve);


    enableListener ();
}

void LCurve::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->labCurve.brightness    = brightness->getValue ();
    pp->labCurve.contrast      = (int)contrast->getValue ();
	pp->labCurve.saturation      = (int)saturation->getValue ();

    pp->labCurve.lcurve         = lshape->getCurve ();
	pp->labCurve.acurve         = ashape->getCurve ();
    pp->labCurve.bcurve         = bshape->getCurve ();

    if (pedited) {
        pedited->labCurve.brightness = brightness->getEditedState ();
        pedited->labCurve.contrast = contrast->getEditedState ();
		pedited->labCurve.saturation = saturation->getEditedState ();

        pedited->labCurve.lcurve    = !lshape->isUnChanged ();
		pedited->labCurve.acurve    = !ashape->isUnChanged ();
        pedited->labCurve.bcurve    = !bshape->isUnChanged ();
    }
}

void LCurve::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    brightness->setDefault (defParams->labCurve.brightness);
    contrast->setDefault (defParams->labCurve.contrast);
	saturation->setDefault (defParams->labCurve.saturation);

    if (pedited) {
        brightness->setDefaultEditedState (pedited->labCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->labCurve.contrast ? Edited : UnEdited);
		saturation->setDefaultEditedState (pedited->labCurve.saturation ? Edited : UnEdited);

    }
    else {
        brightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
		saturation->setDefaultEditedState (Irrelevant);
    }
}

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
    else
        costr = Glib::ustring::format ((int)a->getValue());

    if (a==brightness)
        listener->panelChanged (EvLBrightness, costr);
    else if (a==contrast)
        listener->panelChanged (EvLContrast, costr);
	else if (a==saturation)
        listener->panelChanged (EvLSaturation, costr);
}

//attempt to hide unused channels
/*void LCurve::channel_changed () {
	
    removeIfThere (this, lcurve, false);
    removeIfThere (this, acurve, false);
	removeIfThere (this, bcurve, false);

    if (channel->get_active_row_number()==0) 
        pack_start (*lcurve);   
    else if (channel->get_active_row_number()==1) 
        pack_start (*acurve);
	else if (method->get_active_row_number()==2) 
        pack_start (*bcurve);

    
    if (listener && enabled->get_active ())
        listener->panelChanged (EvLabCurvetype, channel->get_active_text ());
    
}*/

void LCurve::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    brightness->showEditedCB ();
    contrast->showEditedCB ();
	saturation->showEditedCB ();

    curveEditorG->setBatchMode (batchMode);
}

void LCurve::setAdjusterBehavior (bool bradd, bool contradd, bool satadd) {

    if ((!brAdd && bradd) || (brAdd && !bradd))
        brightness->setLimits (-100, 100, 1, 0);
    if ((!contrAdd && contradd) || (contrAdd && !contradd))
        contrast->setLimits (-100, 100, 1, 0);
	if ((!satAdd && satadd) || (satAdd && !satadd))
        saturation->setLimits (-100, 100, 1, 0);

    brAdd = bradd;
    contrAdd = contradd;
	satAdd = satadd;

}

void LCurve::updateCurveBackgroundHistogram (unsigned* hist) {
    
    lshape->updateBackgroundHistogram (hist);
}
