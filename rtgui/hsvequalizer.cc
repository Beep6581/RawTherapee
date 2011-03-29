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
#include <utils.h>

using namespace rtengine;
using namespace rtengine::procparams;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HSVEqualizer::HSVEqualizer () : Gtk::VBox(), FoldableToolPanel(this) {
	
	curveEditorG = new CurveEditorGroup (M("TP_HSVEQUALIZER_CHANNEL"));
	curveEditorG->setCurveListener (this);
	curveEditorG->setColorProvider (this);

	hshape = (FlatCurveEditor*)curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_HUE"));
	sshape = (FlatCurveEditor*)curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_SAT"));
	vshape = (FlatCurveEditor*)curveEditorG->addCurve(CT_Flat, M("TP_HSVEQUALIZER_VAL"));

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

void HSVEqualizer::colorForValue (double valX, double valY) {

	float r, g, b;

	CurveEditor* ce = curveEditorG->getDisplayedCurve();

	if (ce == hshape) {        // Hue = f(Hue)

		float h = (float)((valY - 0.5) * 2. + valX);
		if (h > 1.0)
			h -= 1.0;
		else if (h < 0.0)
			h += 1.0;
		hsv2rgb(h, (float)0.5, (float)0.5, r, g, b);
		red = (double)r;
		green = (double)g;
		blue = (double)b;
	}
	else if (ce == sshape) {   // Saturation = f(Hue)
		hsv2rgb((float)valX, (float)valY, (float)0.5, r, g, b);
		red = (double)r;
		green = (double)g;
		blue = (double)b;
	}
	else if (ce == vshape) {   // Value = f(Hue)
		hsv2rgb((float)valX, (float)0.5, (float)valY, r, g, b);
		red = (double)r;
		green = (double)g;
		blue = (double)b;
	}
	else {
		printf("Error: no curve displayed!\n");
	}

}

void HSVEqualizer::setBatchMode (bool batchMode) {
	
    ToolPanel::setBatchMode (batchMode);

    curveEditorG->setBatchMode (batchMode);
}
