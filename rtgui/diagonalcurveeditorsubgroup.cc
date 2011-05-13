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
#include <guiutils.h>
#include <multilangmgr.h>
#include <guiutils.h>
#include <mycurve.h>
#include <shcselector.h>
#include <adjuster.h>
#include <mycurve.h>
#include <curveeditor.h>
#include <diagonalcurveeditorsubgroup.h>

DiagonalCurveEditorSubGroup::DiagonalCurveEditorSubGroup (CurveEditorGroup* prt) {

	valLinear = (int)DCT_Linear;
	valUnchanged = (int)DCT_Unchanged;
	parent = prt;

	activeParamControl = -1;

	// custom curve
	customCurveBox = new Gtk::HBox ();
	customCurveBox->set_spacing(4);
	customCurve = Gtk::manage (new MyDiagonalCurve ());
	customCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
	customCurve->setType (DCT_Spline);
	//customCurve->set_tooltip_text (M("CURVEEDITOR_TOOLTIPMOVESPEED"));
	customCurveBox->pack_start (*customCurve, Gtk::PACK_EXPAND_WIDGET, 0);

	Gtk::VBox* custombbox = Gtk::manage (new Gtk::VBox ());
	custombbox->set_spacing(4);
	saveCustom = Gtk::manage (new Gtk::Button ());
	saveCustom->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
	loadCustom = Gtk::manage (new Gtk::Button ());
	loadCustom->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

	custombbox->pack_end (*saveCustom, Gtk::PACK_SHRINK, 0);
	custombbox->pack_end (*loadCustom, Gtk::PACK_SHRINK, 0);

	customCurveBox->pack_end (*custombbox, Gtk::PACK_SHRINK, 0);
	customCurveBox->show_all ();

	saveCustom->signal_clicked().connect( sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::savePressed) );
	loadCustom->signal_clicked().connect( sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::loadPressed) );
	saveCustom->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
	loadCustom->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));

	// NURBS curve
	NURBSCurveBox = new Gtk::HBox ();
	NURBSCurveBox->set_spacing(4);
	NURBSCurve = Gtk::manage (new MyDiagonalCurve ());
	NURBSCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
	NURBSCurve->setType (DCT_NURBS);
	//customCurve->set_tooltip_text (M("CURVEEDITOR_TOOLTIPMOVESPEED"));
	NURBSCurveBox->pack_start (*NURBSCurve, Gtk::PACK_EXPAND_WIDGET, 0);

	Gtk::VBox* NURBSbbox = Gtk::manage (new Gtk::VBox ());
	NURBSbbox->set_spacing(4);
	saveNURBS = Gtk::manage (new Gtk::Button ());
	saveNURBS->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
	loadNURBS = Gtk::manage (new Gtk::Button ());
	loadNURBS->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

	NURBSbbox->pack_end (*saveNURBS, Gtk::PACK_SHRINK, 0);
	NURBSbbox->pack_end (*loadNURBS, Gtk::PACK_SHRINK, 0);

	NURBSCurveBox->pack_end (*NURBSbbox, Gtk::PACK_SHRINK, 0);
	NURBSCurveBox->show_all ();

	saveNURBS->signal_clicked().connect( sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::savePressed) );
	loadNURBS->signal_clicked().connect( sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::loadPressed) );
	saveNURBS->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
	loadNURBS->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));

	// parametric curve
	paramCurveBox = new Gtk::VBox ();
	paramCurveBox->set_spacing(4);
	paramCurve = Gtk::manage (new MyDiagonalCurve ());
	Gtk::Table* paramctab = Gtk::manage (new Gtk::Table (2,1));
	paramCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
	paramCurve->setType (DCT_Parametric);
	shcSelector = Gtk::manage (new SHCSelector ());
	shcSelector->set_size_request (GRAPH_SIZE, 20);

	paramctab->attach (*paramCurve, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 0, 0);
	paramctab->attach (*shcSelector, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, RADIUS, 0);

	paramCurveBox->pack_start (*paramctab, Gtk::PACK_EXPAND_WIDGET, 0);

	highlights = Gtk::manage (new Adjuster (M("CURVEEDITOR_HIGHLIGHTS"), -100, 100, 1, 0));
	lights     = Gtk::manage (new Adjuster (M("CURVEEDITOR_LIGHTS"), -100, 100, 1, 0));
	darks      = Gtk::manage (new Adjuster (M("CURVEEDITOR_DARKS"), -100, 100, 1, 0));
	shadows    = Gtk::manage (new Adjuster (M("CURVEEDITOR_SHADOWS"), -100, 100, 1, 0));

	Gtk::EventBox* evhighlights = Gtk::manage (new Gtk::EventBox ());
	Gtk::EventBox* evlights = Gtk::manage (new Gtk::EventBox ());
	Gtk::EventBox* evdarks = Gtk::manage (new Gtk::EventBox ());
	Gtk::EventBox* evshadows = Gtk::manage (new Gtk::EventBox ());

	evhighlights->add (*highlights);
	evlights->add (*lights);
	evdarks->add (*darks);
	evshadows->add (*shadows);

	paramCurveBox->pack_start (*Gtk::manage (new Gtk::HSeparator ()));
	paramCurveBox->pack_start (*evhighlights);
	paramCurveBox->pack_start (*evlights);
	paramCurveBox->pack_start (*evdarks);
	paramCurveBox->pack_start (*evshadows);
	paramCurveBox->show_all ();

	customCurveBox->reference ();
	paramCurveBox->reference ();

	customCurve->setCurveListener (parent); // Send the message directly to the parent
	NURBSCurve->setCurveListener (parent);
	paramCurve->setCurveListener (parent);
	shcSelector->setSHCListener (this);

	highlights->setAdjusterListener (this);
	lights->setAdjusterListener (this);
	darks->setAdjusterListener (this);
	shadows->setAdjusterListener (this);

	evhighlights->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evlights->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evdarks->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evshadows->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evhighlights->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterEntered), 4));
	evlights->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterEntered), 5));
	evdarks->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterEntered), 6));
	evshadows->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterEntered), 7));
	evhighlights->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterLeft), 4));
	evlights->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterLeft), 5));
	evdarks->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterLeft), 6));
	evshadows->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &DiagonalCurveEditorSubGroup::adjusterLeft), 7));
}

DiagonalCurveEditorSubGroup::~DiagonalCurveEditorSubGroup() {
    delete customCurveBox;
    delete paramCurveBox;
    delete NURBSCurveBox;
}

/*
 * Add a new curve to the curves list
 */
DiagonalCurveEditor* DiagonalCurveEditorSubGroup::addCurve(Glib::ustring curveLabel) {
	DiagonalCurveEditor* newCE = new DiagonalCurveEditor(curveLabel, parent, this);

	// Initialization of the new curve
	storeCurveValues(newCE, getCurveFromGUI(DCT_Spline));
	storeCurveValues(newCE, getCurveFromGUI(DCT_Parametric));
	storeCurveValues(newCE, getCurveFromGUI(DCT_NURBS));

	// We add it to the curve editor list
	parent->curveEditors.push_back(newCE);
	return newCE;
}

/*
 * Switch the editor widgets to the currently edited curve
 */
void DiagonalCurveEditorSubGroup::switchGUI() {

	removeEditor();

	DiagonalCurveEditor* dCurve = (DiagonalCurveEditor*)(parent->displayedCurve);

	if (dCurve) {

		// Initializing GUI values + repacking the appropriated widget
		//dCurve->typeconn.block(true);


		switch((DiagonalCurveType)(dCurve->curveType->getSelected())) {
		case (DCT_Spline):
			customCurve->setPoints (dCurve->customCurveEd);
			parent->pack_start (*customCurveBox);
			customCurveBox->check_resize();
			customCurve->forceResize();
			break;
		case (DCT_Parametric):
			paramCurve->setPoints (dCurve->paramCurveEd);
			shcSelector->setPositions (
					dCurve->paramCurveEd.at(1),
					dCurve->paramCurveEd.at(2),
					dCurve->paramCurveEd.at(3)
			);
			highlights->setValue (dCurve->paramCurveEd.at(4));
			lights->setValue (dCurve->paramCurveEd.at(5));
			darks->setValue (dCurve->paramCurveEd.at(6));
			shadows->setValue (dCurve->paramCurveEd.at(7));
			parent->pack_start (*paramCurveBox);
			paramCurve->forceResize();
			break;
		case (DCT_NURBS):
			NURBSCurve->setPoints (dCurve->NURBSCurveEd);
			parent->pack_start (*NURBSCurveBox);
			NURBSCurveBox->check_resize();
			NURBSCurve->forceResize();
			break;
		default:	// (DCT_Linear, DCT_Unchanged)
			// ... do nothing
			break;
		}

		//dCurve->typeconn.block(false);
	}
}

void DiagonalCurveEditorSubGroup::savePressed () {

	Glib::ustring fname = outputFile();
	if (fname.size()) {
		std::ofstream f (fname.c_str());
		std::vector<double> p;
		//std::vector<double> p = customCurve->getPoints ();

		switch (parent->displayedCurve->selected) {
		case DCT_Spline:		// custom
			p = customCurve->getPoints ();
			break;
		case DCT_NURBS:		// NURBS
			p = NURBSCurve->getPoints ();
			break;
		default:
			break;
		}

		int ix = 0;
		if (p[ix]==(double)(DCT_Linear))
			f << "Linear\n";
		else if (p[ix]==(double)(DCT_Spline))
			f << "Spline\n";
		else if (p[ix]==(double)(DCT_NURBS))
			f << "NURBS\n";
		ix++;
		for (unsigned int i=0; i<p.size()/2; i++, ix+=2)
			f << p[ix] << ' ' << p[ix+1] << std::endl;
		f.close ();
	}
}

void DiagonalCurveEditorSubGroup::loadPressed () {

	Glib::ustring fname = inputFile();
	if (fname.size()) {
		std::ifstream f (fname.c_str());
		if (f) {
			std::vector<double> p;
			std::string s;
			f >> s;
			if (s=="Linear")
				p.push_back ((double)(DCT_Linear));
			else if (s=="Spline")
				p.push_back ((double)(DCT_Spline));
			else if (s=="NURBS")
				p.push_back ((double)(DCT_NURBS));
			else return;
			double x;
			while (f) {
				f >> x;
				if (f)
					p.push_back (x);
			}
			if (p[0] == (double)(DCT_Spline)) {
				customCurve->setPoints (p);
				customCurve->queue_draw ();
				customCurve->notifyListener ();
			}
			else if (p[0] == (double)(DCT_NURBS)) {
				NURBSCurve->setPoints (p);
				NURBSCurve->queue_draw ();
				NURBSCurve->notifyListener ();
			}
		}
    }
}

/*
 * Store the curves of the currently displayed type from the widgets to the CurveEditor object
 */
void DiagonalCurveEditorSubGroup::storeDisplayedCurve() {
	if (parent->displayedCurve) {
		switch (parent->displayedCurve->selected) {
		case (DCT_Spline):
			storeCurveValues(parent->displayedCurve, getCurveFromGUI(DCT_Spline));
			break;
		case (DCT_Parametric):
			storeCurveValues(parent->displayedCurve, getCurveFromGUI(DCT_Parametric));
			break;
		case (DCT_NURBS):
			storeCurveValues(parent->displayedCurve, getCurveFromGUI(DCT_NURBS));
			break;
		default:
			break;
		}
	}
}

/*
 * Restore the histogram to all types from the CurveEditor object to the widgets
 */
void DiagonalCurveEditorSubGroup::restoreDisplayedHistogram() {
	if (parent->displayedCurve) {
		paramCurve->updateBackgroundHistogram (parent->displayedCurve->bgHistValid ? parent->displayedCurve->histogram : NULL);
		customCurve->updateBackgroundHistogram (parent->displayedCurve->bgHistValid ? parent->displayedCurve->histogram : NULL);
		NURBSCurve->updateBackgroundHistogram (parent->displayedCurve->bgHistValid ? parent->displayedCurve->histogram : NULL);
	}

}

void DiagonalCurveEditorSubGroup::storeCurveValues (CurveEditor* ce, const std::vector<double>& p) {
	if (p.size()) {
		DiagonalCurveType t = (DiagonalCurveType)p[0];
		for (int i=0; i<(int)p.size(); i++)

		switch (t) {
		case (DCT_Spline):
			((DiagonalCurveEditor*)ce)->customCurveEd = p;
			break;
		case (DCT_Parametric):
			((DiagonalCurveEditor*)ce)->paramCurveEd = p;
			break;
		case (DCT_NURBS):
			((DiagonalCurveEditor*)ce)->NURBSCurveEd = p;
			break;
		default:
			break;
		}
	}
}

/*
 * Called to update the parametric curve graph with new slider values
 */
const std::vector<double> DiagonalCurveEditorSubGroup::getCurveFromGUI (int type) {
	switch ((DiagonalCurveType)type) {
	case (DCT_Parametric): {
		std::vector<double> lcurve (8);
		lcurve[0] = (double)(DCT_Parametric);
		shcSelector->getPositions (lcurve[1], lcurve[2], lcurve[3]);
		lcurve[4] = highlights->getValue ();
		lcurve[5] = lights->getValue ();
		lcurve[6] = darks->getValue ();
		lcurve[7] = shadows->getValue ();
		return lcurve;
		}
	case (DCT_Spline):
        return customCurve->getPoints ();
	case (DCT_NURBS):
        return NURBSCurve->getPoints ();
	default: {
		// linear and other solutions
		std::vector<double> lcurve (1);
		lcurve[0] = (double)(DCT_Linear);
		return lcurve;
		}
	}
}

/*
 * Unlink the tree editor widgets from their parent box to hide them
 */
void DiagonalCurveEditorSubGroup::removeEditor () {
	removeIfThere (parent, customCurveBox, false);
	removeIfThere (parent, paramCurveBox, false);
	removeIfThere (parent, NURBSCurveBox, false);
}

bool DiagonalCurveEditorSubGroup::curveReset(int cType) {
	switch ((DiagonalCurveType) cType) {
	case (DCT_NURBS) :	// = Control cage
		NURBSCurve->reset ();
		return true;
		break;
	case (DCT_Spline) :	// = Custom
		customCurve->reset ();
		return true;
		break;
	case (DCT_Parametric) :
		highlights->resetPressed();
		lights->resetPressed();
		darks->resetPressed();
		shadows->resetPressed();
		shcSelector->reset();
		paramCurve->reset ();
		return true;
		break;
	default:
		return false;
		break;
	}
}

void DiagonalCurveEditorSubGroup::setColorProvider (ColorProvider* p) {
	customCurve->setColorProvider(p);
	paramCurve->setColorProvider(p);
	NURBSCurve->setColorProvider(p);
}

/*
 * Listener
 */
void DiagonalCurveEditorSubGroup::shcChanged () {

    paramCurve->setPoints (getCurveFromGUI(DCT_Parametric));
	storeDisplayedCurve();
	parent->curveChanged ();
}

/*
 * Listener
 */
void DiagonalCurveEditorSubGroup::adjusterChanged (Adjuster* a, double newval) {

    paramCurve->setPoints (getCurveFromGUI(DCT_Parametric));
	storeDisplayedCurve();
	parent->curveChanged ();
}

/*
 * Listener called when the mouse is over a parametric curve's slider
 */
bool DiagonalCurveEditorSubGroup::adjusterEntered (GdkEventCrossing* ev, int ac) {

    if (ev->detail != GDK_NOTIFY_INFERIOR) {
        activeParamControl = ac;
        paramCurve->setActiveParam (activeParamControl);
    }
    return true;
}

/*
 * Listener called when the mouse left the parametric curve's slider
 */
bool DiagonalCurveEditorSubGroup::adjusterLeft (GdkEventCrossing* ev, int ac) {

    if (ev->detail != GDK_NOTIFY_INFERIOR) {
        activeParamControl = -1;
        paramCurve->setActiveParam (activeParamControl);
    }
    return true;
}

void DiagonalCurveEditorSubGroup::updateBackgroundHistogram (CurveEditor* ce) {
	if (ce==parent->displayedCurve) {
		paramCurve->updateBackgroundHistogram (ce->bgHistValid ? ce->histogram : NULL);
		customCurve->updateBackgroundHistogram (ce->bgHistValid ? ce->histogram : NULL);
		NURBSCurve->updateBackgroundHistogram (ce->bgHistValid ? ce->histogram : NULL);
	}
}
