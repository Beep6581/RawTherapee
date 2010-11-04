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
 *
 *  Class created by Jean-Christophe FRISCH, aka 'Hombre'
 */
#include <curveeditorgroup.h>
#include <multilangmgr.h>

extern Glib::ustring argv0;

CurveEditorGroup::CurveEditorGroup (Glib::ustring groupLabel) : cl(NULL), activeParamControl(-1) {
	curveEditors.clear();
	displayedCurve = 0;
	numberOfPackedCurve = 0;

	// We set the label to the one provided as parameter, even if it's an empty string
	curveGroupLabel = Gtk::manage (new Gtk::Label (groupLabel+":", Gtk::ALIGN_LEFT));

	// custom curve
	customCurveBox = new Gtk::HBox ();
	Gtk::HBox* tmpa = Gtk::manage (new Gtk::HBox ());
	customCurve = Gtk::manage (new MyCurve ());
	//Gtk::AspectFrame* af = Gtk::manage (new Gtk::AspectFrame ("",Gtk::ALIGN_CENTER,Gtk::ALIGN_CENTER,1,false));
	//af->add (*customCurve);
	customCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
	customCurve->setType (Spline);
	//customCurve->set_tooltip_text (M("CURVEEDITOR_TOOLTIPMOVESPEED"));
	tmpa->pack_start (*customCurve, true, false, 4);
	customCurveBox->pack_start (*tmpa, true, true,4);

	Gtk::VBox* custombbox = Gtk::manage (new Gtk::VBox ());
	saveCustom = Gtk::manage (new Gtk::Button ());
	saveCustom->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
	loadCustom = Gtk::manage (new Gtk::Button ());
	loadCustom->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

	custombbox->pack_end (*saveCustom, Gtk::PACK_SHRINK, 4);
	custombbox->pack_end (*loadCustom, Gtk::PACK_SHRINK, 4);

	customCurveBox->pack_end (*custombbox, Gtk::PACK_SHRINK, 0);
	customCurveBox->show_all ();

	saveCustom->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditorGroup::savePressed) );
	loadCustom->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditorGroup::loadPressed) );
	saveCustom->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
	loadCustom->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));

	// NURBS curve
	NURBSCurveBox = new Gtk::HBox ();
	Gtk::HBox* tmpb = Gtk::manage (new Gtk::HBox ());
	NURBSCurve = Gtk::manage (new MyCurve ());
	//Gtk::AspectFrame* af = Gtk::manage (new Gtk::AspectFrame ("",Gtk::ALIGN_CENTER,Gtk::ALIGN_CENTER,1,false));
	//af->add (*customCurve);
	NURBSCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
	NURBSCurve->setType (NURBS);
	//customCurve->set_tooltip_text (M("CURVEEDITOR_TOOLTIPMOVESPEED"));
	tmpb->pack_start (*NURBSCurve, true, false, 4);
	NURBSCurveBox->pack_start (*tmpb, true, true,4);

	Gtk::VBox* NURBSbbox = Gtk::manage (new Gtk::VBox ());
	saveNURBS = Gtk::manage (new Gtk::Button ());
	saveNURBS->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
	loadNURBS = Gtk::manage (new Gtk::Button ());
	loadNURBS->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

	NURBSbbox->pack_end (*saveNURBS, Gtk::PACK_SHRINK, 4);
	NURBSbbox->pack_end (*loadNURBS, Gtk::PACK_SHRINK, 4);

	NURBSCurveBox->pack_end (*NURBSbbox, Gtk::PACK_SHRINK, 0);
	NURBSCurveBox->show_all ();

	saveNURBS->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditorGroup::savePressed) );
	loadNURBS->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditorGroup::loadPressed) );
	saveNURBS->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
	loadNURBS->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));

	// parametric curve
	paramCurveBox = new Gtk::VBox ();
	paramCurve = Gtk::manage (new MyCurve ());
	Gtk::Table* paramctab = Gtk::manage (new Gtk::Table (2,1));
	//Gtk::AspectFrame* afp = Gtk::manage (new Gtk::AspectFrame ("",Gtk::ALIGN_CENTER,Gtk::ALIGN_CENTER,1,false));
	//afp->add (*paramCurve);
	paramCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
	paramCurve->setType (Parametric);
	shcSelector = Gtk::manage (new SHCSelector ());
	shcSelector->set_size_request (GRAPH_SIZE, 20);

	paramctab->attach (*paramCurve, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 4, 4);
	paramctab->attach (*shcSelector, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, RADIUS+4, 0);

	Gtk::HBox* tmpc = Gtk::manage (new Gtk::HBox ());
	tmpc->pack_start (*paramctab, true, false);

	paramCurveBox->pack_start (*tmpc, true, true);

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

	customCurve->setCurveListener (this);
	NURBSCurve->setCurveListener (this);
	paramCurve->setCurveListener (this);
	shcSelector->setSHCListener (this);

	highlights->setAdjusterListener (this);
	lights->setAdjusterListener (this);
	darks->setAdjusterListener (this);
	shadows->setAdjusterListener (this);

	evhighlights->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evlights->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evdarks->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evshadows->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK);
	evhighlights->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterEntered), 4));
	evlights->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterEntered), 5));
	evdarks->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterEntered), 6));
	evshadows->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterEntered), 7));
	evhighlights->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterLeft), 4));
	evlights->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterLeft), 5));
	evdarks->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterLeft), 6));
	evshadows->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditorGroup::adjusterLeft), 7));

}

CurveEditorGroup::~CurveEditorGroup() {
    for (std::vector<CurveEditor*>::iterator i = curveEditors.begin(); i != curveEditors.end(); ++i)
    {
        delete *i;
    }
    delete customCurveBox;
    delete paramCurveBox;
    delete NURBSCurveBox;
}

void CurveEditorGroup::hideCurrentCurve() {
	// Setting the curve type to 'Unchanged' hide the CurveEditor
	if (displayedCurve)
		displayedCurve->curveType->set_active(false);
}

/*
 * Add a new curve to the curves list
 */
CurveEditor* CurveEditorGroup::addCurve(Glib::ustring curveLabel) {
	CurveEditor* newCE = new CurveEditor(curveLabel, this);

	// Initialization of the new curve
	storeCurveValues(newCE, getCurveFromGUI(Spline));
	storeCurveValues(newCE, getCurveFromGUI(Parametric));
	storeCurveValues(newCE, getCurveFromGUI(NURBS));

	// We add it to the curve editor list
	curveEditors.push_back(newCE);
	return newCE;
}

/*
 * Use this method to start a new line of button
 */
void CurveEditorGroup::newLine() {
	Gtk::HBox* headerBox;

	if (curveEditors.size() > numberOfPackedCurve) {
		headerBox = Gtk::manage (new Gtk::HBox ());

		if (!numberOfPackedCurve) {
			headerBox->pack_start(*curveGroupLabel, Gtk::PACK_SHRINK, 2);

			curve_reset = Gtk::manage (new Gtk::Button ());
			curve_reset->add (*Gtk::manage (new Gtk::Image (argv0+"/images/undo.png")));
			curve_reset->set_relief (Gtk::RELIEF_NONE);
			curve_reset->set_border_width (0);
			curve_reset->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLINEAR"));
			curve_reset->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditorGroup::curveResetPressed) );

			headerBox->pack_end (*curve_reset, Gtk::PACK_SHRINK, 0);
			curve_reset->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditorGroup::curveResetPressed) );
		}

		int j = numberOfPackedCurve;
		for (int i = (int)(curveEditors.size())-1; i >= j; i--)
		{
			headerBox->pack_end (*curveEditors[i]->curveType, Gtk::PACK_EXPAND_WIDGET, 2);
			numberOfPackedCurve++;
		}

		pack_start (*headerBox, Gtk::PACK_SHRINK, 2);
	}


}

/*
 * Create all the widgets now that the curve list is complete
 * This method should handle all curve number correctly, i.e. eventually display the curve type buttons
 * in a grid (or table)
 */
void CurveEditorGroup::curveListComplete() {
	newLine();

	// We check the length of the label ; if it contains only one char (':'), we set it to the right default string
	if (curveGroupLabel->get_label().size()==1)
		curveGroupLabel->set_label(M(curveEditors.size() > 1 ? "CURVEEDITOR_CURVES" : "CURVEEDITOR_CURVE") + ":");

	if (curveEditors.size() > 1)
		cl->setMulti(true);
}

/*
 * Callback method used when a curve type button has changed ;
 * it will activate the button, and so emit 'signal_toggled' (-> curveTypeToggled here under)
 */
void CurveEditorGroup::typeSelectionChanged (CurveEditor* ce, int n) {
	// Same type : do nothing
	if (ce==displayedCurve && (CurveType)n==ce->selected)
		return;

	if ((CurveType)(n)<Unchanged)
		ce->selected = (CurveType)n;

	// The user selected a new type from a toggled off button
	if (ce!=displayedCurve)
		// We toggle off the other curve: it will emit the toggle off signal
		hideCurrentCurve();

	// If the button was not pressed before
	if (!ce->curveType->get_active()) {
		storeDisplayedCurve();
		// We set it pressed : it will emit the toggle on signal and update the GUI
		ce->curveType->set_active( n>Linear && n<Unchanged );
		if (n==Linear || n==Unchanged) {
			// Since we do not activate the curve when the user switch the a toggled off button to 'Linear', we have to
			// to call the curve listener manually, because 'curveChanged' uses displayedCurve...
		    if (cl) {
		    	if (cl->isMulti())
		    		cl->curveChanged (ce);
		    	else
		    		cl->curveChanged ();
		    }
		}
		else
			curveChanged ();
	}
	else {
		// The button is already pressed so we switch the GUI ourselves
	   	switchGUI();
	   	curveChanged ();
	}
}

/*
 * Callback method used when a button has been toggled on/off
 * It then hide any other displayed curve and display it's curve
 */
void CurveEditorGroup::curveTypeToggled(CurveEditor* ce) {
	bool curveRestored = false;

	// Looking for the button state
	if (ce->curveType->get_active()) {
		// The button is now pressed, so we have to first hide all other CurveEditor
		hideCurrentCurve();

		displayedCurve = ce;

		if (ce->curveType->getSelected()==Unchanged) {
			curveRestored = true;
			ce->curveType->setSelected(ce->selected);
		}

		// then show this CurveEditor
		int ct = ce->curveType->getSelected();
		if (ct < Unchanged)
			restoreDisplayedHistogram();
	}
	else {
		// The button is now released, so we have to hide this CurveEditor
		displayedCurve = 0;
	}
   	switchGUI();

   	if (curveRestored)
   	    curveChanged ();

}

/*
 * Switch the editor widgets to the currently edited curve
 */
void CurveEditorGroup::switchGUI() {

	removeEditor();

	if (displayedCurve) {

		// Initializing GUI values + repacking the appropriated widget
		//displayedCurve->typeconn.block(true);

		switch((CurveType)(displayedCurve->curveType->getSelected())) {
		case (Spline):
			customCurve->setPoints (displayedCurve->customCurveEd);
			pack_start (*customCurveBox);
			break;
		case (Parametric):
			paramCurve->setPoints (displayedCurve->paramCurveEd);
			shcSelector->setPositions (
					displayedCurve->paramCurveEd.at(1),
					displayedCurve->paramCurveEd.at(2),
					displayedCurve->paramCurveEd.at(3)
			);
			highlights->setValue (displayedCurve->paramCurveEd.at(4));
			lights->setValue (displayedCurve->paramCurveEd.at(5));
			darks->setValue (displayedCurve->paramCurveEd.at(6));
			shadows->setValue (displayedCurve->paramCurveEd.at(7));
			pack_start (*paramCurveBox);
			break;
		case (NURBS):
			NURBSCurve->setPoints (displayedCurve->NURBSCurveEd);
			pack_start (*NURBSCurveBox);
			break;
		default:	// (Linear, Unchanged)
			// ... do nothing
			break;
		}

		//displayedCurve->typeconn.block(false);
	}
}

void CurveEditorGroup::savePressed () {

    Gtk::FileChooserDialog dialog(M("CURVEEDITOR_SAVEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
//    if (options.multiUser)
//        dialog.set_current_folder (Options::rtdir + "/" + options.profilePath);
//    else
//        dialog.set_current_folder (argv0 + "/" + options.profilePath);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("CURVEEDITOR_FILEDLGFILTERCURVE"));
    filter_pp.add_pattern("*.rtc");
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("CURVEEDITOR_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

    dialog.set_do_overwrite_confirmation (true);

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK) {

        std::string fname = dialog.get_filename();

        if (getExtension (fname)!="rtc")
            fname = fname + ".rtc";

        if (Glib::file_test (fname, Glib::FILE_TEST_EXISTS)) {
            Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
            Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
            int response = msgd.run ();
            if (response==Gtk::RESPONSE_NO)
                return;
        }

        std::ofstream f (fname.c_str());
        std::vector<double> p = customCurve->getPoints ();

        switch (displayedCurve->selected) {
        case Spline:		// custom
        	p = customCurve->getPoints ();
        	break;
        case NURBS:		// NURBS
        	p = NURBSCurve->getPoints ();
        	break;
        default:
        	break;
        }

        int ix = 0;
        if (p[ix]==(double)(Linear))
            f << "Linear\n";
        else if (p[ix]==(double)(Spline))
            f << "Spline\n";
        else if (p[ix]==(double)(NURBS))
            f << "NURBS\n";
        ix++;
        for (unsigned int i=0; i<p.size()/2; i++, ix+=2)
            f << p[ix] << ' ' << p[ix+1] << std::endl;
        f.close ();
    }
}

void CurveEditorGroup::loadPressed () {

    Gtk::FileChooserDialog dialog(M("CURVEEDITOR_LOADDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-open"), Gtk::RESPONSE_OK);

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("CURVEEDITOR_FILEDLGFILTERCURVE"));
    filter_pp.add_pattern("*.rtc");
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("CURVEEDITOR_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK) {
        std::ifstream f (dialog.get_filename().c_str());
        if (f) {
            std::vector<double> p;
            std::string s;
            f >> s;
            if (s=="Linear")
                p.push_back ((double)(Linear));
            else if (s=="Spline")
                p.push_back ((double)(Spline));
            else if (s=="NURBS")
                p.push_back ((double)(NURBS));
            else return;
            double x;
            while (f) {
                f >> x;
                if (f)
                    p.push_back (x);
            }
            if (p[0] == (double)(Spline)) {
				customCurve->setPoints (p);
				customCurve->queue_draw ();
				customCurve->notifyListener ();
            }
            else if (p[0] == (double)(NURBS)) {
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
void CurveEditorGroup::storeDisplayedCurve() {
	if (displayedCurve) {
		switch (displayedCurve->selected) {
		case (Spline):
			storeCurveValues(displayedCurve, getCurveFromGUI(Spline));
			break;
		case (Parametric):
			storeCurveValues(displayedCurve, getCurveFromGUI(Parametric));
			break;
		case (NURBS):
			storeCurveValues(displayedCurve, getCurveFromGUI(NURBS));
			break;
		default:
			break;
		}
	}
}

/*
 * Restore the histogram to all types from the CurveEditor object to the widgets
 */
void CurveEditorGroup::restoreDisplayedHistogram() {
	if (displayedCurve) {
		paramCurve->updateBackgroundHistogram (displayedCurve->bgHistValid ? displayedCurve->histogram : NULL);
		customCurve->updateBackgroundHistogram (displayedCurve->bgHistValid ? displayedCurve->histogram : NULL);
		NURBSCurve->updateBackgroundHistogram (displayedCurve->bgHistValid ? displayedCurve->histogram : NULL);
	}

}

void CurveEditorGroup::storeCurveValues (CurveEditor* ce, const std::vector<double>& p) {
	if (p.size()) {
		CurveType t = (CurveType)p[0];
		for (int i=0; i<(int)p.size(); i++)

		switch (t) {
		case (Spline):
			ce->customCurveEd = p;
			break;
		case (Parametric):
			ce->paramCurveEd = p;
			break;
		case (NURBS):
			ce->NURBSCurveEd = p;
			break;
		default:
			break;
		}
	}
}

/*
 * Update the GUI if the given curveEditor is currently displayed
 */
void CurveEditorGroup::updateGUI (CurveEditor* ce) {
	if (!ce) {
    	return;
    }

	// we update the curve type button to the corresponding curve type, only if it is not currently set to 'Unchanged'
	if (ce->curveType->getSelected()<Unchanged)
		ce->curveType->setSelected(ce->selected);

    // if not displayed or "unchanged" is selected, do not change gui
    if (ce==displayedCurve && ce->curveType->getSelected()<Unchanged) {
    	switchGUI();
    }
}

/*
 * Called from the outside to set the curve type & values
 */
void CurveEditorGroup::setCurveExternal (CurveEditor* ce, const std::vector<double>& c) {
	if (c.size()) {
		storeCurveValues(ce, c);			// The new curve is saved in the CurveEditor
		ce->selected = (CurveType)(c[0]);	// We set the selected curve type in the CurveEditor to the one of the specified curve
	}
	updateGUI(ce);							// And we update the GUI if necessary
}

/*
 * Called to update the parametric curve graph with new slider values
 */
const std::vector<double> CurveEditorGroup::getCurveFromGUI (CurveType type) {
	switch (type) {
	case (Parametric): {
		std::vector<double> lcurve (8);
		lcurve[0] = (double)(Parametric);
		shcSelector->getPositions (lcurve[1], lcurve[2], lcurve[3]);
		lcurve[4] = highlights->getValue ();
		lcurve[5] = lights->getValue ();
		lcurve[6] = darks->getValue ();
		lcurve[7] = shadows->getValue ();
		return lcurve;
		}
	case (Spline):
        return customCurve->getPoints ();
	case (NURBS):
        return NURBSCurve->getPoints ();
	default: {
		// linear and other solutions
		std::vector<double> lcurve (1);
		lcurve[0] = (double)(Linear);
		return lcurve;
		}
	}
}

/*
 * Unlink the tree editor widgets from their parent box to hide them
 */
void CurveEditorGroup::removeEditor () {
    removeIfThere (this, customCurveBox, false);
    removeIfThere (this, paramCurveBox, false);
    removeIfThere (this, NURBSCurveBox, false);
}

/*
 * Listener called when the user has modified the curve
 */
void CurveEditorGroup::curveChanged () {

	storeDisplayedCurve();
    if (cl) {
    	if (cl->isMulti())
    		cl->curveChanged (displayedCurve);
    	else
    		cl->curveChanged ();
    }
}

/*
 * Listener
 */
void CurveEditorGroup::shcChanged () {

    paramCurve->setPoints (getCurveFromGUI(Parametric));
	storeDisplayedCurve();
	if (cl->isMulti())
		cl->curveChanged (displayedCurve);
	else
		cl->curveChanged ();
}

/*
 * Listener
 */
void CurveEditorGroup::adjusterChanged (Adjuster* a, double newval) {

    paramCurve->setPoints (getCurveFromGUI(Parametric));
	storeDisplayedCurve();
	if (cl->isMulti())
		cl->curveChanged (displayedCurve);
	else
		cl->curveChanged ();
}

/*
 * Listener called when the mouse is over a parametric curve's slider
 */
bool CurveEditorGroup::adjusterEntered (GdkEventCrossing* ev, int ac) {

    if (ev->detail != GDK_NOTIFY_INFERIOR) {
        activeParamControl = ac;
        paramCurve->setActiveParam (activeParamControl);
    }
    return true;
}

/*
 * Listener called when the mouse left the parametric curve's slider
 */
bool CurveEditorGroup::adjusterLeft (GdkEventCrossing* ev, int ac) {

    if (ev->detail != GDK_NOTIFY_INFERIOR) {
        activeParamControl = -1;
        paramCurve->setActiveParam (activeParamControl);
    }
    return true;
}

/*
 * Call back method when the reset button is pressed :
 * reset the currently toggled on curve editor
 */
void CurveEditorGroup::curveResetPressed() {
	if (displayedCurve) {
		switch (displayedCurve->selected) {
		case (NURBS) :	// = Control cage
			NURBSCurve->reset ();
			curveChanged ();
			break;
		case (Spline) :	// = Custom
			customCurve->reset ();
			curveChanged ();
			break;
		case (Parametric) :
			highlights->resetPressed();
			lights->resetPressed();
			darks->resetPressed();
			shadows->resetPressed();
			shcSelector->reset();
			paramCurve->reset ();
			curveChanged ();
			break;
		default:
			break;
		}
	}
}

void CurveEditorGroup::setBatchMode (bool batchMode) {
	for (std::vector<CurveEditor*>::iterator i = curveEditors.begin(); i != curveEditors.end(); ++i) {
		(*i)->curveType->addEntry(argv0+"/images/curveType-unchanged.png", M("GENERAL_UNCHANGED"));
		(*i)->curveType->show();
	}
}

void CurveEditorGroup::setUnChanged (bool uc, CurveEditor* ce) {
	if (uc) {
		// the user selected several thumbnails, so we hide the editors and set the curveEditor selection to 'Unchanged'
		//ce->typeconn.block(true);
		// we hide the editor widgets
		hideCurrentCurve();
		// the curve type selected option is set to unchanged
		ce->curveType->setSelected(Unchanged);
		//ce->typeconn.block(false);
	}
	else {
		// we want it to use back the 'CurveEditor::setCurve' memorized in CurveEditor::tempCurve
		//ce->typeconn.block(true);
		// we switch back the curve type selected option to the one of the used curve
		ce->curveType->setSelected(ce->selected);
		updateGUI (ce);
		//ce->typeconn.block(false);
	}
}

void CurveEditorGroup::updateBackgroundHistogram (CurveEditor* ce) {
	if (ce==displayedCurve) {
		paramCurve->updateBackgroundHistogram (ce->bgHistValid ? ce->histogram : NULL);
		customCurve->updateBackgroundHistogram (ce->bgHistValid ? ce->histogram : NULL);
		NURBSCurve->updateBackgroundHistogram (ce->bgHistValid ? ce->histogram : NULL);
		printf(" - fait! (ce==displayedCurve)");
	}
}
