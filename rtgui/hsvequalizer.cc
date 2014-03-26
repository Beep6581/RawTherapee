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

#include "hsvequalizer.h"
#include "../rtengine/color.h"

using namespace rtengine;
using namespace rtengine::procparams;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HSVEqualizer::HSVEqualizer () : FoldableToolPanel(this) {
	
	std::vector<GradientMilestone> bottomMilestones;
	float R, G, B;
	// -0.1 rad < Hue < 1.6 rad
	for (int i=0; i<7; i++) {
		float x = float(i)*(1.0f/6.0);
		Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
		bottomMilestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
	}

	curveEditorG = new CurveEditorGroup (options.lastHsvCurvesDir, M("TP_HSVEQUALIZER_CHANNEL"));
	curveEditorG->setCurveListener (this);

	hshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_HUE")));
	hshape->setEditID(EUID_HSV_H, BT_SINGLEPLANE_FLOAT);
	hshape->setBottomBarBgGradient(bottomMilestones);
	//hshape->setLeftBarColorProvider(this);  Not working yet
	hshape->setCurveColorProvider(this, 1);

	sshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_SAT")));
	sshape->setEditID(EUID_HSV_S, BT_SINGLEPLANE_FLOAT);
	sshape->setBottomBarBgGradient(bottomMilestones);
	//sshape->setLeftBarColorProvider(this);  Not working yet
	sshape->setCurveColorProvider(this, 2);

	vshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_VAL")));
	vshape->setEditID(EUID_HSV_V, BT_SINGLEPLANE_FLOAT);
	vshape->setBottomBarBgGradient(bottomMilestones);
	//vshape->setLeftBarColorProvider(this);  Not working yet
	vshape->setCurveColorProvider(this, 3);

	// This will add the reset button at the end of the curveType buttons
	curveEditorG->curveListComplete();

	pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);

	//curveEditorG->show();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HSVEqualizer::~HSVEqualizer () {
	delete curveEditorG;
}


void HSVEqualizer::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        hshape->setUnChanged (!pedited->hsvequalizer.hcurve);
        sshape->setUnChanged (!pedited->hsvequalizer.scurve);
        vshape->setUnChanged (!pedited->hsvequalizer.vcurve);
    }

    hshape->setCurve         (pp->hsvequalizer.hcurve);
	sshape->setCurve         (pp->hsvequalizer.scurve);
    vshape->setCurve         (pp->hsvequalizer.vcurve);
	
    enableListener ();
}

void HSVEqualizer::setEditProvider (EditDataProvider *provider) {
    hshape->setEditProvider(provider);
    sshape->setEditProvider(provider);
    vshape->setEditProvider(provider);
}

void HSVEqualizer::autoOpenCurve () {
    // Open up the first curve if selected
    bool active = hshape->openIfNonlinear();
    if (!active) sshape->openIfNonlinear();
    if (!active) vshape->openIfNonlinear();
}

void HSVEqualizer::write (ProcParams* pp, ParamsEdited* pedited) {
	
    pp->hsvequalizer.hcurve = hshape->getCurve ();
    pp->hsvequalizer.scurve = sshape->getCurve ();
    pp->hsvequalizer.vcurve = vshape->getCurve ();

    if (pedited) {

        pedited->hsvequalizer.hcurve = !hshape->isUnChanged ();
        pedited->hsvequalizer.scurve = !sshape->isUnChanged ();
        pedited->hsvequalizer.vcurve = !vshape->isUnChanged ();
    }
}

/*
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
*/

/*
 * Curve listener
 *
 * If more than one curve has been added, the curve listener is automatically
 * set to 'multi=true', and send a pointer of the modified curve in a parameter
 */
void HSVEqualizer::curveChanged (CurveEditor* ce) {

    if (listener) {
    	if (ce == hshape)
        	listener->panelChanged (EvHSVEqualizerH, M("HISTORY_CUSTOMCURVE"));
    	if (ce == sshape)
        	listener->panelChanged (EvHSVEqualizerS, M("HISTORY_CUSTOMCURVE"));
    	if (ce == vshape)
        	listener->panelChanged (EvHSVEqualizerV, M("HISTORY_CUSTOMCURVE"));
	}
}

/*
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
*/

void HSVEqualizer::colorForValue (double valX, double valY, int callerId, ColorCaller* caller) {

	float r, g, b;

	if (callerId == 1) {        // Hue = f(Hue)

		float h = float((valY - 0.5) * 2. + valX);
		if (h > 1.0f)
			h -= 1.0f;
		else if (h < 0.0f)
			h += 1.0f;
		Color::hsv2rgb01(h, 0.5f, 0.5f, r, g, b);
		caller->ccRed = double(r);
		caller->ccGreen = double(g);
		caller->ccBlue = double(b);
	}
	else if (callerId == 2) {   // Saturation = f(Hue)
		Color::hsv2rgb01(float(valX), float(valY), 0.5f, r, g, b);
		caller->ccRed = double(r);
		caller->ccGreen = double(g);
		caller->ccBlue = double(b);
	}
	else if (callerId == 3) {   // Value = f(Hue)
		Color::hsv2rgb01(float(valX), 0.5f, float(valY), r, g, b);
		caller->ccRed = double(r);
		caller->ccGreen = double(g);
		caller->ccBlue = double(b);
	}
	else {
		printf("Error: no curve displayed!\n");
	}

}

void HSVEqualizer::setBatchMode (bool batchMode) {
	
    ToolPanel::setBatchMode (batchMode);

    curveEditorG->setBatchMode (batchMode);
}
