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

#include <gtkmm.h>
#include <fstream>
#include <string>
#include "guiutils.h"
#include "multilangmgr.h"
#include "guiutils.h"
#include "mycurve.h"
#include "shcselector.h"
#include "adjuster.h"
#include "mycurve.h"
#include "curveeditor.h"
#include "flatcurveeditorsubgroup.h"

FlatCurveEditorSubGroup::FlatCurveEditorSubGroup (CurveEditorGroup* prt) {

	valLinear = (int)FCT_Linear;
	valUnchanged = (int)FCT_Unchanged;
	parent = prt;

	// ControlPoints curve
	CPointsCurveBox = new Gtk::HBox ();
	CPointsCurveBox->set_spacing(4);
	CPointsCurve = Gtk::manage (new MyFlatCurve ());
	//CPointsCurve->set_size_request (GRAPH_SIZE+2*RADIUS+1, GRAPH_SIZE+2*RADIUS+1);
	CPointsCurve->setType (FCT_MinMaxCPoints);
	CPointsCurveBox->pack_start (*CPointsCurve, Gtk::PACK_EXPAND_WIDGET, 0);

	Gtk::VBox* CPointsbbox = Gtk::manage (new Gtk::VBox ());
	CPointsbbox->set_spacing(4);
	saveCPoints = Gtk::manage (new Gtk::Button ());
	saveCPoints->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
	loadCPoints = Gtk::manage (new Gtk::Button ());
	loadCPoints->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

	CPointsbbox->pack_end (*saveCPoints, Gtk::PACK_SHRINK, 0);
	CPointsbbox->pack_end (*loadCPoints, Gtk::PACK_SHRINK, 0);

	CPointsCurveBox->pack_end (*CPointsbbox, Gtk::PACK_SHRINK, 0);
	CPointsCurveBox->show_all ();

	saveCPoints->signal_clicked().connect( sigc::mem_fun(*this, &FlatCurveEditorSubGroup::savePressed) );
	loadCPoints->signal_clicked().connect( sigc::mem_fun(*this, &FlatCurveEditorSubGroup::loadPressed) );
	saveCPoints->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
	loadCPoints->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));

	CPointsCurve->setCurveListener (parent); // Send the message directly to the parent
}

FlatCurveEditorSubGroup::~FlatCurveEditorSubGroup() {
    delete CPointsCurveBox;
}

/*
 * Add a new curve to the curves list
 */
FlatCurveEditor* FlatCurveEditorSubGroup::addCurve(Glib::ustring curveLabel, bool isPeriodic) {
	FlatCurveEditor* newCE = new FlatCurveEditor(curveLabel, parent, this, isPeriodic);

	// Initialization of the new curve
	storeCurveValues(newCE, getCurveFromGUI(FCT_MinMaxCPoints));

	// We add it to the curve editor list
	parent->curveEditors.push_back(newCE);
	return newCE;
}

/*
 * Switch the editor widgets to the currently edited curve
 */
void FlatCurveEditorSubGroup::switchGUI() {

	removeEditor();

	FlatCurveEditor* dCurve = static_cast<FlatCurveEditor*>(parent->displayedCurve);

	if (dCurve) {

		// Initializing GUI values + repacking the appropriated widget
		//dCurve->typeconn.block(true);

		switch((FlatCurveType)(dCurve->curveType->getSelected())) {
		case (FCT_MinMaxCPoints):
			CPointsCurve->setPeriodicity(dCurve->periodic);		// Setting Periodicity before setting points
			CPointsCurve->setPoints (dCurve->controlPointsCurveEd);
			parent->pack_start (*CPointsCurveBox);
			CPointsCurveBox->check_resize();
			CPointsCurve->forceResize();
			break;
		default:	// (DCT_Linear, DCT_Unchanged)
			// ... do nothing
			break;
		}

		//dCurve->typeconn.block(false);
	}
}

void FlatCurveEditorSubGroup::savePressed () {

	Glib::ustring fname = outputFile();
	if (fname.size()) {
		std::ofstream f (fname.c_str());
		std::vector<double> p;
		//std::vector<double> p = customCurve->getPoints ();

		switch (parent->displayedCurve->selected) {
		case FCT_MinMaxCPoints:		// Control points
			p = CPointsCurve->getPoints ();
			break;
		default:
			break;
		}

		int ix = 0;
		if (p[ix]==(double)(FCT_Linear))
			f << "Linear\n";
		else if (p[ix]==(double)(FCT_MinMaxCPoints))
			f << "ControlPoints\n";
		ix++;
		for (unsigned int i=0; i<p.size()/2; i++, ix+=2)
			f << p[ix] << ' ' << p[ix+1] << std::endl;
		f.close ();
	}
}

void FlatCurveEditorSubGroup::loadPressed () {

	Glib::ustring fname = inputFile();
	if (fname.size()) {
		std::ifstream f (fname.c_str());
		if (f) {
			std::vector<double> p;
			std::string s;
			f >> s;
			if (s=="Linear")
				p.push_back ((double)(FCT_Linear));
			else if (s=="ControlPoints")
				p.push_back ((double)(FCT_MinMaxCPoints));
			else return;
			double x;
			while (f) {
				f >> x;
				if (f)
					p.push_back (x);
			}
			if (p[0] == (double)(FCT_MinMaxCPoints)) {
				CPointsCurve->setPoints (p);
				CPointsCurve->queue_draw ();
				CPointsCurve->notifyListener ();
			}
		}
    }
}

/*
 * Store the curves of the currently displayed type from the widgets to the CurveEditor object
 */
void FlatCurveEditorSubGroup::storeDisplayedCurve() {
	if (parent->displayedCurve) {
		switch (parent->displayedCurve->selected) {
		/*case (FCT_Parametric):
			storeCurveValues(parent->displayedCurve, getCurveFromGUI(FCT_Parametric));
			break;*/
		case (FCT_MinMaxCPoints):
			storeCurveValues(parent->displayedCurve, getCurveFromGUI(FCT_MinMaxCPoints));
			break;
		default:
			break;
		}
	}
}

/*
 * Restore the histogram to all types from the CurveEditor object to the widgets
 */
void FlatCurveEditorSubGroup::restoreDisplayedHistogram() {
	if (parent->displayedCurve) {
		//paramCurve->updateBackgroundHistogram (parent->displayedCurve->histogram);
		CPointsCurve->updateBackgroundHistogram (parent->displayedCurve->histogram);
	}

}

void FlatCurveEditorSubGroup::storeCurveValues (CurveEditor* ce, const std::vector<double>& p) {
	if (!p.empty()) {
		FlatCurveType t = static_cast<FlatCurveType>(p[0]);
		for (int i=0; i<(int)p.size(); i++)

		switch (t) {
		case (FCT_MinMaxCPoints):
			static_cast<FlatCurveEditor*>(ce)->controlPointsCurveEd = p;
			break;
		default:
			break;
		}
	}
}

/*
 * Called to update the parametric curve graph with new slider values
 */
const std::vector<double> FlatCurveEditorSubGroup::getCurveFromGUI (int type) {
	switch ((FlatCurveType)type) {
	case (FCT_MinMaxCPoints):
        return CPointsCurve->getPoints ();
	default: {
		// linear and other solutions
		std::vector<double> lcurve (1);
		lcurve[0] = (double)(FCT_Linear);
		return lcurve;
		}
	}
}

/*
 * Unlink the tree editor widgets from their parent box to hide them
 */
void FlatCurveEditorSubGroup::removeEditor () {
    removeIfThere (parent, CPointsCurveBox, false);
}

bool FlatCurveEditorSubGroup::curveReset(int cType) {
	switch ((FlatCurveType) cType) {
	case (FCT_MinMaxCPoints) :	// = Control cage
		CPointsCurve->reset ();
		return true;
		break;
	/*case (FCT_Parametric) :
		highlights->resetPressed();
		lights->resetPressed();
		darks->resetPressed();
		shadows->resetPressed();
		shcSelector->reset();
		paramCurve->reset ();
		return true;
		break;*/
	default:
		return false;
		break;
	}
}

void FlatCurveEditorSubGroup::setColorProvider (ColorProvider* p) {
	CPointsCurve->setColorProvider(p);
}

/*void FlatCurveEditorSubGroup::updateBackgroundHistogram (CurveEditor* ce) {
	CurveEditor* fce = (CurveEditor*)ce;
	if (fce==displayedCurve) {
		paramCurve->updateBackgroundHistogram (fce->bgHistValid ? fce->histogram : NULL);
		customCurve->updateBackgroundHistogram (fce->bgHistValid ? fce->histogram : NULL);
		NURBSCurve->updateBackgroundHistogram (fce->bgHistValid ? fce->histogram : NULL);
	}
}*/
