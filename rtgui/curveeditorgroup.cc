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

#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "diagonalcurveeditorsubgroup.h"
#include "flatcurveeditorsubgroup.h"
#include "multilangmgr.h"
#include "../rtengine/safegtk.h"
#include "rtimage.h"

extern Glib::ustring argv0;

CurveEditorGroup::CurveEditorGroup (Glib::ustring groupLabel) : cl(NULL), cp(NULL) {
	curveEditors.clear();
	displayedCurve = 0;
	numberOfPackedCurve = 0;
	flatSubGroup = 0;
	diagonalSubGroup = 0;

	// We set the label to the one provided as parameter, even if it's an empty string
	curveGroupLabel = Gtk::manage (new Gtk::Label (groupLabel+":", Gtk::ALIGN_LEFT));
}

CurveEditorGroup::~CurveEditorGroup() {
    for (std::vector<CurveEditor*>::iterator i = curveEditors.begin(); i != curveEditors.end(); ++i)
    {
        delete *i;
    }
    delete flatSubGroup;
    delete diagonalSubGroup;
}

void CurveEditorGroup::hideCurrentCurve() {
	// Setting the curve type to 'Unchanged' hide the CurveEditor
	if (displayedCurve)
		displayedCurve->curveType->set_active(false);
}

/*
 * Add a new curve to the curves list
 *
 * The "periodic" parameter is only used by flat curve editors
 */
CurveEditor* CurveEditorGroup::addCurve(CurveType cType, Glib::ustring curveLabel, bool periodic) {
	switch (cType) {
	case (CT_Diagonal):
		if (!diagonalSubGroup) {
			diagonalSubGroup = new DiagonalCurveEditorSubGroup(this);
		}
		return (CurveEditor*)diagonalSubGroup->addCurve(curveLabel);
	case (CT_Flat):
		if (!flatSubGroup) {
			flatSubGroup = new FlatCurveEditorSubGroup(this);
		}
		return (CurveEditor*)flatSubGroup->addCurve(curveLabel, periodic);
	default:
		return (CurveEditor*)NULL;
		break;
	}
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
			curve_reset->add (*Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png")));
			curve_reset->set_relief (Gtk::RELIEF_NONE);
			curve_reset->set_border_width (0);
			curve_reset->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLINEAR"));
			curve_reset->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditorGroup::curveResetPressed) );

			headerBox->pack_end (*curve_reset, Gtk::PACK_SHRINK, 0);
		}

		int j = numberOfPackedCurve;
		for (int i = (int)(curveEditors.size())-1; i >= j; i--)
		{
			headerBox->pack_end (*curveEditors[i]->curveType->buttonGroup, Gtk::PACK_EXPAND_WIDGET, 2);
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

	// Set the color provider
	if (cp) {
		if (flatSubGroup) flatSubGroup->setColorProvider(cp);
		if (diagonalSubGroup) diagonalSubGroup->setColorProvider(cp);
	}

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
	if (ce==displayedCurve && n==(int)ce->selected)
		return;

	if (n<ce->subGroup->valUnchanged)
		ce->selected = n;

	// The user selected a new type from a toggled off button
	if (ce!=displayedCurve)
		// We toggle off the other curve: it will emit the toggle off signal
		hideCurrentCurve();

	// If the button was not pressed before
	if (!ce->curveType->get_active()) {
		ce->subGroup->storeDisplayedCurve();
		// We set it pressed : it will emit the toggle on signal and update the GUI
		ce->curveType->set_active( n>ce->subGroup->valLinear && n<ce->subGroup->valUnchanged );
		if (n==ce->subGroup->valLinear || n==ce->subGroup->valUnchanged) {
			// Since we do not activate the curve when the user switch the toggled off button to 'Linear', we have to
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
		ce->subGroup->switchGUI();
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

		if (ce->curveType->getSelected()==ce->subGroup->valUnchanged) {
			curveRestored = true;
			ce->curveType->setSelected(ce->selected);
		}

		// then show this CurveEditor
		int ct = ce->curveType->getSelected();
		if (ct < ce->subGroup->valUnchanged)
			ce->subGroup->restoreDisplayedHistogram();
	}
	else {
		// The button is now released, so we have to hide this CurveEditor
		displayedCurve = 0;
	}
	ce->subGroup->switchGUI();

   	if (curveRestored)
   	    curveChanged ();

}

/*
 * Update the GUI if the given curveEditor is currently displayed
 */
void CurveEditorGroup::updateGUI (CurveEditor* ce) {
	if (!ce) {
    	return;
    }

	// we update the curve type button to the corresponding curve type, only if it is not currently set to 'Unchanged'
	if (ce->curveType->getSelected()<ce->subGroup->valUnchanged)
		ce->curveType->setSelected(ce->selected);

    // if not displayed or "unchanged" is selected, do not change gui
    if (ce==displayedCurve && ce->curveType->getSelected()<ce->subGroup->valUnchanged) {
    	ce->subGroup->switchGUI();
    }
}

/*
 * Called from the outside to set the curve type & values
 */
void CurveEditorGroup::setCurveExternal (CurveEditor* ce, const std::vector<double>& c) {
	if (!c.empty()) {
		ce->subGroup->storeCurveValues(ce, c);	// The new curve is saved in the CurveEditor
		(ce)->selected = c[0];		// We set the selected curve type in the CurveEditor to the one of the specified curve
	}
	updateGUI((CurveEditor*)ce);				// And we update the GUI if necessary
}

/*
 * Listener called when the user has modified the curve
 */
void CurveEditorGroup::curveChanged () {

	displayedCurve->subGroup->storeDisplayedCurve();
    if (cl) {
    	if (cl->isMulti())
    		cl->curveChanged (displayedCurve);
    	else
    		cl->curveChanged ();
    }
}

/*
 * Call back method when the reset button is pressed :
 * reset the currently toggled on curve editor
 */
void CurveEditorGroup::curveResetPressed () {
	if (displayedCurve) {
		if (displayedCurve->subGroup->curveReset(displayedCurve->selected)) {
			curveChanged();
		}
	}
}

void CurveEditorGroup::setBatchMode (bool batchMode) {
	for (std::vector<CurveEditor*>::iterator i = curveEditors.begin(); i != curveEditors.end(); ++i) {
		(*i)->curveType->addEntry("curveType-unchanged.png", M("GENERAL_UNCHANGED"));
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
		ce->curveType->setSelected(ce->subGroup->valUnchanged);
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

Glib::ustring CurveEditorSubGroup::outputFile () {

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

	Glib::ustring fname;
    do {
		int result = dialog.run();

		fname = dialog.get_filename();

		if (result==Gtk::RESPONSE_OK) {

			if (getExtension (fname)!="rtc")
				fname = fname + ".rtc";

			if (safe_file_test (fname, Glib::FILE_TEST_EXISTS)) {
				Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
				Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
				int response = msgd.run ();
				if (response==Gtk::RESPONSE_YES)
					break;
			}
			else
				break;
		}
		else {
			fname = "";
			break;
		}
    } while (1);

    return fname;
}

Glib::ustring CurveEditorSubGroup::inputFile () {

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

   	Glib::ustring fname;
    if (result==Gtk::RESPONSE_OK) {
   		fname = dialog.get_filename();
    	if (Glib::file_test (fname, Glib::FILE_TEST_EXISTS))
    		return fname;
    }
    fname = "";
    return fname;
}
