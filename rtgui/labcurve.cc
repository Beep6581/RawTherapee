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
#include "../rtengine/improcfun.h"

using namespace rtengine;
using namespace rtengine::procparams;

LCurve::LCurve () : Gtk::VBox(), FoldableToolPanel(this) {

	std::vector<GradientMilestone> bottomMilestones;

	brightness = Gtk::manage (new Adjuster (M("TP_LABCURVE_BRIGHTNESS"), -100., 100., 1., 0.));
	contrast   = Gtk::manage (new Adjuster (M("TP_LABCURVE_CONTRAST"), -100., 100., 1., 0.));
	chromaticity   = Gtk::manage (new Adjuster (M("TP_LABCURVE_CHROMATICITY"), -100., 100., 1., 0.));

	pack_start (*brightness);
	brightness->show ();

	pack_start (*contrast);
	contrast->show ();

	pack_start (*chromaticity);
	chromaticity->show ();
		
	brightness->setAdjusterListener (this);
	contrast->setAdjusterListener (this);
	chromaticity->setAdjusterListener (this);
	
	//%%%%%%%%%%%%%%%%%%
	pack_start (*Gtk::manage (new  Gtk::HSeparator()));
	
	bwtoning = Gtk::manage (new Gtk::CheckButton (M("TP_LABCURVE_BWTONING")));
	bwtoning->set_tooltip_markup (M("TP_LABCURVE_BWTONING_TIP"));
	pack_start (*bwtoning);

	avoidcolorshift = Gtk::manage (new Gtk::CheckButton (M("TP_LABCURVE_AVOIDCOLORSHIFT")));
	avoidcolorshift->set_tooltip_text (M("TP_LABCURVE_AVOIDCOLORSHIFT_TOOLTIP"));
	pack_start (*avoidcolorshift, Gtk::PACK_SHRINK, 4);

	rstprotection = Gtk::manage ( new Adjuster (M("TP_LABCURVE_RSTPROTECTION"), 0., 100., 0.1, 0.) );
	pack_start (*rstprotection);
	rstprotection->show ();
	
	rstprotection->setAdjusterListener (this);
	bwtconn= bwtoning->signal_toggled().connect( sigc::mem_fun(*this, &LCurve::bwtoning_toggled) );
	acconn = avoidcolorshift->signal_toggled().connect( sigc::mem_fun(*this, &LCurve::avoidcolorshift_toggled) );
	//%%%%%%%%%%%%%%%%%%%

	Gtk::HSeparator *hsep3 = Gtk::manage (new  Gtk::HSeparator());
	hsep3->show ();
	pack_start (*hsep3);

	curveEditorG = new CurveEditorGroup (options.lastLabCurvesDir);
	curveEditorG->setCurveListener (this);

	ccshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_LABCURVE_CURVEEDITOR_CC")));
	ccshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_CC_TOOLTIP"));
	ccshape->setRangeLabels(
			M("TP_LABCURVE_CURVEEDITOR_CC_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE2"),
			M("TP_LABCURVE_CURVEEDITOR_CC_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE4")
	);
	ccshape->setRangeDefaultMilestones(0.05, 0.2, 0.58);
	//cbgshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_LABCURVE_CURVEEDITOR_CBG"));
	//cbgshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_CBG_TOOLTIP"));
	// -0.1 rad < Hue < 1.6 rad
	for (int i=0; i<7; i++) {
		float R, G, B;
		float x = float(i)*(1.0f/6.0);
		Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
		bottomMilestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
	}
	chshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_LABCURVE_CURVEEDITOR_CH")));
	chshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_CH_TOOLTIP"));
	chshape->setBottomBarBgGradient(bottomMilestones);
	chshape->setCurveColorProvider(this);

	curveEditorG->newLine();

	bottomMilestones.clear();
	bottomMilestones.push_back( GradientMilestone(0., 0., 0., 0.) );
	bottomMilestones.push_back( GradientMilestone(1., 1., 1., 1.) );
	lshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "L"));
	lshape->setBottomBarBgGradient(bottomMilestones);
	lshape->setLeftBarBgGradient(bottomMilestones);
	ashape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "a"));
	ashape->setRangeLabels(
			M("TP_LABCURVE_CURVEEDITOR_A_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_A_RANGE2"),
			M("TP_LABCURVE_CURVEEDITOR_A_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_A_RANGE4")
	);
	bshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "b"));
	bshape->setRangeLabels(
			M("TP_LABCURVE_CURVEEDITOR_B_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_B_RANGE2"),
			M("TP_LABCURVE_CURVEEDITOR_B_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_B_RANGE4")
	);

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
		chromaticity->setEditedState (pedited->labCurve.chromaticity ? Edited : UnEdited);
		
		//%%%%%%%%%%%%%%%%%%%%%%
		rstprotection->setEditedState (pedited->labCurve.rstprotection ? Edited : UnEdited);
		bwtoning->set_inconsistent (!pedited->labCurve.bwtoning);
        avoidcolorshift->set_inconsistent (!pedited->labCurve.avoidcolorshift);
		//%%%%%%%%%%%%%%%%%%%%%%

        lshape->setUnChanged   (!pedited->labCurve.lcurve);
        ashape->setUnChanged   (!pedited->labCurve.acurve);
        bshape->setUnChanged   (!pedited->labCurve.bcurve);
        ccshape->setUnChanged  (!pedited->labCurve.cccurve);
        chshape->setUnChanged  (!pedited->labCurve.chcurve);
        //cbgshape->setUnChanged (!pedited->labCurve.cbgcurve);
    }
    else {
        //if bwtoning is enabled, chromaticity value, avoid color shift and rstprotection has no effect
        //ccshape->set_sensitive(!(!pp->labCurve.bwtoning));
        //chshape->set_sensitive(!(!pp->labCurve.bwtoning));
        chromaticity->set_sensitive(!pp->labCurve.bwtoning);
        rstprotection->set_sensitive( !pp->labCurve.bwtoning && pp->labCurve.chromaticity!=0 );
        avoidcolorshift->set_sensitive(!pp->labCurve.bwtoning);
    }

    brightness->setValue    (pp->labCurve.brightness);
    contrast->setValue      (pp->labCurve.contrast);
	chromaticity->setValue  (pp->labCurve.chromaticity);
	
	//%%%%%%%%%%%%%%%%%%%%%%
	rstprotection->setValue (pp->labCurve.rstprotection);

    bwtconn.block (true);
    acconn.block (true);
    bwtoning->set_active (pp->labCurve.bwtoning);
    avoidcolorshift->set_active (pp->labCurve.avoidcolorshift);
    bwtconn.block (false);
    acconn.block (false);

    lastBWTVal = pp->labCurve.bwtoning;
    lastACVal = pp->labCurve.avoidcolorshift;
	//%%%%%%%%%%%%%%%%%%%%%%

    lshape->setCurve   (pp->labCurve.lcurve);
    ashape->setCurve   (pp->labCurve.acurve);
    bshape->setCurve   (pp->labCurve.bcurve);
    ccshape->setCurve  (pp->labCurve.cccurve);
    chshape->setCurve  (pp->labCurve.chcurve);
    //cbgshape->setCurve (pp->labCurve.cbgcurve);

    // Open up the first curve if selected
    bool active = lshape->openIfNonlinear();
    if (!active) ashape->openIfNonlinear();
    if (!active) bshape->openIfNonlinear();
    if (!active) ccshape->openIfNonlinear();
    if (!active) chshape->openIfNonlinear();

    queue_draw();

    enableListener ();
}

void LCurve::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->labCurve.brightness    = brightness->getValue ();
    pp->labCurve.contrast      = (int)contrast->getValue ();
	pp->labCurve.chromaticity  = (int)chromaticity->getValue ();
	
	//%%%%%%%%%%%%%%%%%%%%%%
	pp->labCurve.bwtoning        = bwtoning->get_active ();
	pp->labCurve.avoidcolorshift = avoidcolorshift->get_active ();
    pp->labCurve.rstprotection   = rstprotection->getValue ();
	//%%%%%%%%%%%%%%%%%%%%%%

    pp->labCurve.lcurve  = lshape->getCurve ();
	pp->labCurve.acurve  = ashape->getCurve ();
    pp->labCurve.bcurve  = bshape->getCurve ();
    pp->labCurve.cccurve = ccshape->getCurve ();
    pp->labCurve.chcurve = chshape->getCurve ();
    //pp->labCurve.cbgcurve = cbgshape->getCurve ();

    if (pedited) {
        pedited->labCurve.brightness   = brightness->getEditedState ();
        pedited->labCurve.contrast     = contrast->getEditedState ();
		pedited->labCurve.chromaticity = chromaticity->getEditedState ();
		
		//%%%%%%%%%%%%%%%%%%%%%%
		pedited->labCurve.bwtoning        = !bwtoning->get_inconsistent();
		pedited->labCurve.avoidcolorshift = !avoidcolorshift->get_inconsistent();
        pedited->labCurve.rstprotection   = rstprotection->getEditedState ();
		//%%%%%%%%%%%%%%%%%%%%%%

        pedited->labCurve.lcurve    = !lshape->isUnChanged ();
		pedited->labCurve.acurve    = !ashape->isUnChanged ();
        pedited->labCurve.bcurve    = !bshape->isUnChanged ();
        pedited->labCurve.cccurve   = !ccshape->isUnChanged ();
        pedited->labCurve.chcurve   = !chshape->isUnChanged ();
        //pedited->labCurve.cbgcurve  = !bshape->isUnChanged ();
    }
}

void LCurve::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    brightness->setDefault (defParams->labCurve.brightness);
    contrast->setDefault (defParams->labCurve.contrast);
	chromaticity->setDefault (defParams->labCurve.chromaticity);
    rstprotection->setDefault (defParams->labCurve.rstprotection);

    if (pedited) {
        brightness->setDefaultEditedState (pedited->labCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->labCurve.contrast ? Edited : UnEdited);
		chromaticity->setDefaultEditedState (pedited->labCurve.chromaticity ? Edited : UnEdited);
        rstprotection->setDefaultEditedState (pedited->labCurve.rstprotection ? Edited : UnEdited);

    }
    else {
        brightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
		chromaticity->setDefaultEditedState (Irrelevant);
		rstprotection->setDefaultEditedState (Irrelevant);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%
//Color shift control changed
void LCurve::avoidcolorshift_toggled () {

    if (batchMode) {
        if (avoidcolorshift->get_inconsistent()) {
            avoidcolorshift->set_inconsistent (false);
            acconn.block (true);
            avoidcolorshift->set_active (false);
            acconn.block (false);
        }
        else if (lastACVal)
            avoidcolorshift->set_inconsistent (true);

        lastACVal = avoidcolorshift->get_active ();
    }

    if (listener) {
        if (avoidcolorshift->get_active ())
            listener->panelChanged (EvLAvoidColorShift, M("GENERAL_ENABLED"));
        else            
            listener->panelChanged (EvLAvoidColorShift, M("GENERAL_DISABLED"));
    }
}

//%%%%%%%%%%%%%%%%%%%%%%
//BW toning control changed
void LCurve::bwtoning_toggled () {

    if (batchMode) {
        if (bwtoning->get_inconsistent()) {
        	bwtoning->set_inconsistent (false);
        	bwtconn.block (true);
            bwtoning->set_active (false);
            bwtconn.block (false);
        }
        else if (lastBWTVal)
        	bwtoning->set_inconsistent (true);

        lastBWTVal = bwtoning->get_active ();
    }
    else {
        //ccshape->set_sensitive(!(!pp->labCurve.bwtoning));
        //chshape->set_sensitive(!(!pp->labCurve.bwtoning));
        chromaticity->set_sensitive( !(bwtoning->get_active ()) );
        rstprotection->set_sensitive( !(bwtoning->get_active ()) && chromaticity->getIntValue()!=0 );
        avoidcolorshift->set_sensitive( !(bwtoning->get_active ()) );
    }

    if (listener) {
    	if (bwtoning->get_active ())
        	listener->panelChanged (EvLBWtoning, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvLBWtoning, M("GENERAL_DISABLED"));
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
        if (ce == ccshape)
            listener->panelChanged (EvLCCCurve, M("HISTORY_CUSTOMCURVE"));
        if (ce == chshape)
            listener->panelChanged (EvLCHCurve, M("HISTORY_CUSTOMCURVE"));
        //if (ce == cbgshape)
        //    listener->panelChanged (EvLCBGCurve, M("HISTORY_CUSTOMCURVE"));
    }
}

void LCurve::adjusterChanged (Adjuster* a, double newval) {

    if (!listener)
        return;

    Glib::ustring costr;
    if (a==brightness)
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
	else if (a==rstprotection)
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(1), a->getValue());
    else
        costr = Glib::ustring::format ((int)a->getValue());

    if (a==brightness)
        listener->panelChanged (EvLBrightness, costr);
    else if (a==contrast)
        listener->panelChanged (EvLContrast, costr);
	else if (a==chromaticity) {
		if (!batchMode) {
	        rstprotection->set_sensitive( !(bwtoning->get_active ()) && chromaticity->getIntValue()!=0 );
		}
        listener->panelChanged (EvLSaturation, costr);
	}
	else if (a==rstprotection)
        listener->panelChanged (EvLRSTProtection, costr);
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
	else if (ce == ccshape) {    // c = f(c)
		red = (double)valY;
		green = (double)valY;
		blue = (double)valY;
	}
	else if (ce == chshape) {    // c = f(h)

		float r, g, b;

		Color::hsv2rgb01(float(valX), float(valY), 0.5f, r, g, b);

		red = double(r);
		green = double(g);
		blue = double(b);
	}
	
	else {
		printf("Error: no curve displayed!\n");
	}

}

void LCurve::setBatchMode (bool batchMode) {

	ToolPanel::setBatchMode (batchMode);
	brightness->showEditedCB ();
	contrast->showEditedCB ();
	chromaticity->showEditedCB ();
	rstprotection->showEditedCB ();

	curveEditorG->setBatchMode (batchMode);
}


void LCurve::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma){

	lshape->updateBackgroundHistogram (histLCurve);
}

void LCurve::setAdjusterBehavior (bool bradd, bool contradd, bool satadd) {

	brightness->setAddMode(bradd);
	contrast->setAddMode(contradd);
	chromaticity->setAddMode(satadd);
}

void LCurve::trimValues (rtengine::procparams::ProcParams* pp) {

	brightness->trimValue(pp->labCurve.brightness);
	contrast->trimValue(pp->labCurve.contrast);
	chromaticity->trimValue(pp->labCurve.chromaticity);
}
