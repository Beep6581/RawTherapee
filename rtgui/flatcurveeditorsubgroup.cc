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

#include "clipboard.h"
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

FlatCurveEditorSubGroup::FlatCurveEditorSubGroup (CurveEditorGroup* prt, Glib::ustring& curveDir) : CurveEditorSubGroup(curveDir) {

	valLinear = (int)FCT_Linear;
	valUnchanged = (int)FCT_Unchanged;
	parent = prt;

	// ControlPoints curve
	CPointsCurveBox = new Gtk::VBox ();
	CPointsCurveBox->set_spacing(4);
	Gtk::HBox* CPointsCurveAndButtons = Gtk::manage (new Gtk::HBox ());
	CPointsCurveAndButtons->set_spacing(4);
	CPointsCurve = Gtk::manage (new MyFlatCurve ());
	CPointsCurve->set_size_request (GRAPH_SIZE+2*RADIUS+1, GRAPH_SIZE+2*RADIUS+1);
	CPointsCurve->setType (FCT_MinMaxCPoints);

	Gtk::Box* CPointsbbox; // curvebboxpos 0=above, 1=right, 2=below, 3=left
	if (options.curvebboxpos==1 || options.curvebboxpos==3) {
		CPointsbbox = Gtk::manage (new Gtk::VBox ());
	} else {
		CPointsbbox = Gtk::manage (new Gtk::HBox ());
	}
	CPointsbbox->set_spacing(4);

	pasteCPoints = Gtk::manage (new Gtk::Button ());
	pasteCPoints->add (*Gtk::manage (new RTImage ("edit-paste.png")));
	copyCPoints = Gtk::manage (new Gtk::Button ());
	copyCPoints->add (*Gtk::manage (new RTImage ("edit-copy.png")));
	saveCPoints = Gtk::manage (new Gtk::Button ());
	saveCPoints->add (*Gtk::manage (new RTImage ("gtk-save-large.png")));
	loadCPoints = Gtk::manage (new Gtk::Button ());
	loadCPoints->add (*Gtk::manage (new RTImage ("gtk-open.png")));
	editCPoints = Gtk::manage (new Gtk::ToggleButton());
	editCPoints->add (*Gtk::manage (new RTImage ("editmodehand.png")));
	editCPoints->set_tooltip_text(M("EDIT_PIPETTE_TOOLTIP"));
	editCPoints->hide();

	CPointsbbox->pack_end (*pasteCPoints, Gtk::PACK_SHRINK, 0);
	CPointsbbox->pack_end (*copyCPoints, Gtk::PACK_SHRINK, 0);
	CPointsbbox->pack_end (*saveCPoints, Gtk::PACK_SHRINK, 0);
	CPointsbbox->pack_end (*loadCPoints, Gtk::PACK_SHRINK, 0);
	CPointsbbox->pack_start(*editCPoints, Gtk::PACK_SHRINK, 0);

	CPointsCurveAndButtons->pack_start (*CPointsCurve, Gtk::PACK_EXPAND_WIDGET, 0);
	CPointsCurveAndButtons->pack_start (*CPointsbbox, Gtk::PACK_SHRINK, 0);
	CPointsCurveBox->pack_start (*CPointsCurveAndButtons, Gtk::PACK_EXPAND_WIDGET);
	if (options.curvebboxpos==0) {
		removeIfThere (CPointsCurveAndButtons, CPointsbbox, false);
		CPointsCurveBox->pack_start (*CPointsbbox);
		CPointsCurveBox->reorder_child(*CPointsbbox, 0);
	} else if (options.curvebboxpos==2) {
		removeIfThere (CPointsCurveAndButtons, CPointsbbox, false);
		CPointsCurveBox->pack_start (*CPointsbbox);
	} else if (options.curvebboxpos==3) {
		CPointsCurveAndButtons->reorder_child(*CPointsbbox, 0);
	}
	CPointsCurveBox->show_all ();

	saveCPoints->signal_clicked().connect( sigc::mem_fun(*this, &FlatCurveEditorSubGroup::savePressed) );
	loadCPoints->signal_clicked().connect( sigc::mem_fun(*this, &FlatCurveEditorSubGroup::loadPressed) );
	copyCPoints->signal_clicked().connect( sigc::mem_fun(*this, &FlatCurveEditorSubGroup::copyPressed) );
	pasteCPoints->signal_clicked().connect( sigc::mem_fun(*this, &FlatCurveEditorSubGroup::pastePressed) );
	editCPointsConn = editCPoints->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &FlatCurveEditorSubGroup::editToggled), editCPoints) );

	saveCPoints->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
	loadCPoints->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));
	copyCPoints->set_tooltip_text (M("CURVEEDITOR_TOOLTIPCOPY"));
	pasteCPoints->set_tooltip_text (M("CURVEEDITOR_TOOLTIPPASTE"));

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
	return newCE;
}

/*
 * Force the resize of the curve editor, if the displayed one is the requested one
 */
void FlatCurveEditorSubGroup::refresh(CurveEditor *curveToRefresh) {
	if (curveToRefresh != NULL && curveToRefresh == static_cast<FlatCurveEditor*>(parent->displayedCurve)) {
		switch(FlatCurveType(curveToRefresh->curveType->getSelected())) {
		case (FCT_MinMaxCPoints):
			CPointsCurve->refresh();
			break;
		default:	// (DCT_Linear, DCT_Unchanged)
			// ... do nothing
			break;
		}
	}
}

/*
 * Switch off the edit button
 */
void FlatCurveEditorSubGroup::editModeSwitchedOff () {
	// toggling off all edit buttons, even if only one is toggle on
	bool prevState = editCPointsConn.block(true);
	editCPoints->set_active(false);
	CPointsCurve->pipetteMouseOver(NULL, 0);
	CPointsCurve->setDirty(true);
	if (!prevState) editCPointsConn.block(false);
}

void FlatCurveEditorSubGroup::pipetteMouseOver(EditDataProvider *provider, int modifierKey) {
	CurveEditor *curveEditor = static_cast<FlatCurveEditor*>(parent->displayedCurve);
	switch((FlatCurveType)(curveEditor->curveType->getSelected())) {
	case (FCT_MinMaxCPoints):
		CPointsCurve->pipetteMouseOver(provider, modifierKey);
		CPointsCurve->setDirty(true);
		break;
	default:	// (DCT_Linear, DCT_Unchanged)
		// ... do nothing
		break;
	}
}

void FlatCurveEditorSubGroup::pipetteButton1Pressed(EditDataProvider *provider, int modifierKey) {
	CurveEditor *curveEditor = static_cast<FlatCurveEditor*>(parent->displayedCurve);
	switch((FlatCurveType)(curveEditor->curveType->getSelected())) {
	case (FCT_MinMaxCPoints):
		CPointsCurve->pipetteButton1Pressed(provider, modifierKey);
		CPointsCurve->setDirty(true);
		break;
	default:	// (DCT_Linear, DCT_Unchanged)
		// ... do nothing
		break;
	}
}

void FlatCurveEditorSubGroup::pipetteButton1Released(EditDataProvider *provider) {
	CurveEditor *curveEditor = static_cast<FlatCurveEditor*>(parent->displayedCurve);
	switch((FlatCurveType)(curveEditor->curveType->getSelected())) {
	case (FCT_MinMaxCPoints):
		CPointsCurve->pipetteButton1Released(provider);
		CPointsCurve->setDirty(true);
		break;
	default:	// (DCT_Linear, DCT_Unchanged)
		// ... do nothing
		break;
	}
}

void FlatCurveEditorSubGroup::pipetteDrag(EditDataProvider *provider, int modifierKey) {
	CurveEditor *curveEditor = static_cast<FlatCurveEditor*>(parent->displayedCurve);
	switch((FlatCurveType)(curveEditor->curveType->getSelected())) {
	case (FCT_MinMaxCPoints):
		CPointsCurve->pipetteDrag(provider, modifierKey);
		CPointsCurve->setDirty(true);
		break;
	default:	// (DCT_Linear, DCT_Unchanged)
		// ... do nothing
		break;
	}
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

		// first we update the colored bar

		ColorProvider *barColorProvider = dCurve->getLeftBarColorProvider();
		std::vector<GradientMilestone>  bgGradient = dCurve->getLeftBarBgGradient();
		if (barColorProvider == NULL && bgGradient.size() == 0) {
			// dCurve has no left colored bar, so we delete the object
			if (leftBar) {
				delete leftBar;
				leftBar = NULL;
			}
		}
		else {
			// dCurve has a ColorProvider or a background gradient defined, so we create/update the object
			if (!leftBar) {
				leftBar = new ColoredBar(RTO_Bottom2Top);
			}
			if (barColorProvider) {
				bgGradient.clear();
				leftBar->setColorProvider(barColorProvider, dCurve->getLeftBarCallerId());
				leftBar->setBgGradient (bgGradient);
			}
			else {
				leftBar->setColorProvider(NULL, -1);
				leftBar->setBgGradient (bgGradient);
			}
		}

		barColorProvider = dCurve->getBottomBarColorProvider();
		bgGradient = dCurve->getBottomBarBgGradient();
		if (barColorProvider == NULL && bgGradient.size() == 0) {
			// dCurve has no bottom colored bar, so we delete the object
			if (bottomBar) {
				delete bottomBar;
				bottomBar = NULL;
			}
		}
		else {
			// dCurve ave a ColorProvider or a background gradient defined, so we create/update the object
			if (!bottomBar) {
				bottomBar = new ColoredBar(RTO_Left2Right);
			}
			if (barColorProvider) {
				bgGradient.clear();
				bottomBar->setColorProvider(barColorProvider, dCurve->getBottomBarCallerId());
				bottomBar->setBgGradient (bgGradient);
			}
			else {
				bottomBar->setColorProvider(NULL, -1);
				bottomBar->setBgGradient (bgGradient);
			}
		}

		switch((FlatCurveType)(dCurve->curveType->getSelected())) {
		case (FCT_MinMaxCPoints):
			CPointsCurve->setPeriodicity(dCurve->periodic);		// Setting Periodicity before setting points
			CPointsCurve->setPoints (dCurve->controlPointsCurveEd);
			CPointsCurve->setColorProvider(dCurve->getCurveColorProvider(), dCurve->getCurveCallerId());
			CPointsCurve->setColoredBar(leftBar, bottomBar);
			CPointsCurve->forceResize();
			updateEditButton(dCurve, editCPoints, editCPointsConn);
			parent->pack_start (*CPointsCurveBox);
			CPointsCurveBox->check_resize();
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
			f << "Linear" << std::endl;
		else if (p[ix]==(double)(FCT_MinMaxCPoints))
			f << "ControlPoints" << std::endl;
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

void FlatCurveEditorSubGroup::copyPressed () {
// For compatibility use enum FlatCurveType here

	std::vector<double> curve;

	switch (parent->displayedCurve->selected) {
	case FCT_MinMaxCPoints:                // custom
		curve = CPointsCurve->getPoints ();
		clipboard.setFlatCurveData (curve,FCT_MinMaxCPoints);
		break;
	default:                       // (DCT_Linear, DCT_Unchanged)
		// ... do nothing
		break;
	}
}

void FlatCurveEditorSubGroup::pastePressed () {
// For compatibility use enum FlatCurveType here

	std::vector<double> curve;
	FlatCurveType type;

	type = clipboard.hasFlatCurveData();

	if (type == (FlatCurveType)parent->displayedCurve->selected) {
		curve = clipboard.getFlatCurveData ();
		switch (type) {
		case FCT_MinMaxCPoints:    // min/max control points
			CPointsCurve->setPoints (curve);
			CPointsCurve->queue_draw ();
			CPointsCurve->notifyListener ();
			break;
		default:                   // (FCT_Linear, FCT_Unchanged)
			// ... do nothing
			break;
		}
	}
	return;
}

void FlatCurveEditorSubGroup::editToggled (Gtk::ToggleButton *button) {
	FlatCurveEditor* dCurve = static_cast<FlatCurveEditor*>(parent->displayedCurve);
	if (!dCurve)
		// should never happen!
		return;

	if (button->get_active()) {
		dCurve->subscribe();
		CPointsCurve->notifyListener ();

	}
	else {
		dCurve->unsubscribe();
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

bool FlatCurveEditorSubGroup::curveReset(CurveEditor *ce) {
	if (!ce)
		return false;

	FlatCurveEditor *fce = static_cast<FlatCurveEditor*>(ce);

	switch (FlatCurveType(ce->selected)) {
	case (FCT_MinMaxCPoints) :	// = Control cage
		CPointsCurve->reset (fce->controlPointsResetCurve, fce->getIdentityValue());
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
	return true;
}

/*void FlatCurveEditorSubGroup::updateBackgroundHistogram (CurveEditor* ce) {
	CurveEditor* fce = (CurveEditor*)ce;
	if (fce==displayedCurve) {
		paramCurve->updateBackgroundHistogram (fce->bgHistValid ? fce->histogram : NULL);
		customCurve->updateBackgroundHistogram (fce->bgHistValid ? fce->histogram : NULL);
		NURBSCurve->updateBackgroundHistogram (fce->bgHistValid ? fce->histogram : NULL);
	}
}*/
