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
#include "saveasdlg.h"
#include <multilangmgr.h>
#include <guiutils.h>
#include <safegtk.h>
#include <rtimage.h>

extern Options options;
SaveAsDialog::SaveAsDialog (Glib::ustring initialDir) {

	set_title(M("GENERAL_SAVE"));

    Gtk::VBox* vbox = get_vbox ();

    fchooser = Gtk::manage( new Gtk::FileChooserWidget (Gtk::FILE_CHOOSER_ACTION_SAVE) );
    fchooser->set_current_folder (initialDir);

    filter_jpg.set_name(M("SAVEDLG_JPGFILTER"));
    filter_jpg.add_pattern("*.jpg");

    filter_tif.set_name(M("SAVEDLG_JPGFILTER"));
    filter_tif.add_pattern("*.tif");

    filter_png.set_name(M("SAVEDLG_JPGFILTER"));
    filter_png.add_pattern("*.png");

    vbox->pack_start (*fchooser);

    Gtk::HSeparator* hsep1 = Gtk::manage( new Gtk::HSeparator () );
    vbox->pack_start (*hsep1, Gtk::PACK_SHRINK, 2);

// Unique filename option
// ~~~~~~~~~~~~~~~~~~~~~~
    autoSuffix = Gtk::manage( new Gtk::CheckButton (M("SAVEDLG_AUTOSUFFIX")) );
    autoSuffix->set_active(options.autoSuffix);

    vbox->pack_start (*autoSuffix, Gtk::PACK_SHRINK, 4);

    Gtk::HSeparator* hsep2 = Gtk::manage( new Gtk::HSeparator () );
    vbox->pack_start (*hsep2, Gtk::PACK_SHRINK, 2);

// Output Options
// ~~~~~~~~~~~~~~
    formatOpts = Gtk::manage( new SaveFormatPanel () );
    formatOpts->init (options.saveFormat);
    formatOpts->setListener (this);

    vbox->pack_start (*formatOpts, Gtk::PACK_SHRINK, 4);

    Gtk::HSeparator* hsep3 = Gtk::manage( new Gtk::HSeparator () );
    vbox->pack_start (*hsep3, Gtk::PACK_SHRINK, 2);

// queue/immediate
// ~~~~~~~~~~~~~
    immediately    = Gtk::manage( new Gtk::RadioButton (M("SAVEDLG_SAVEIMMEDIATELY")) );
    putToQueueHead = Gtk::manage( new Gtk::RadioButton (M("SAVEDLG_PUTTOQUEUEHEAD")) );
    putToQueueTail = Gtk::manage( new Gtk::RadioButton (M("SAVEDLG_PUTTOQUEUETAIL")) );

    // There is no queue in simple mode, so no need to choose
    if (!simpleEditor) {
    vbox->pack_start (*immediately, Gtk::PACK_SHRINK, 4);
    vbox->pack_start (*putToQueueHead, Gtk::PACK_SHRINK, 4);
    vbox->pack_start (*putToQueueTail, Gtk::PACK_SHRINK, 4);
    }

    immediately->set_active (true);
    Gtk::RadioButton::Group g = immediately->get_group();
    putToQueueHead->set_group (g);
    putToQueueTail->set_group (g);
        
// buttons
// ~~~~~~    
    Gtk::Button* ok     = Gtk::manage( new Gtk::Button (M("GENERAL_OK")) );
    Gtk::Button* cancel = Gtk::manage( new Gtk::Button (M("GENERAL_CANCEL")) );

    ok->set_image (*Gtk::manage(new RTImage ("addtags.png")));
    cancel->set_image (*Gtk::manage(new RTImage ("gtk-cancel.png")));

    ok->signal_clicked().connect( sigc::mem_fun(*this, &SaveAsDialog::okPressed) );
    cancel->signal_clicked().connect( sigc::mem_fun(*this, &SaveAsDialog::cancelPressed) );

    get_action_area()->pack_end (*ok, Gtk::PACK_SHRINK, 4);
    get_action_area()->pack_end (*cancel, Gtk::PACK_SHRINK, 4);

    set_border_width (4);
    show_all_children ();
}

bool SaveAsDialog::getAutoSuffix () {

    return autoSuffix->get_active();
}

bool SaveAsDialog::getImmediately () {

    return immediately->get_active ();
}

bool SaveAsDialog::getToHeadOfQueue () {

    return putToQueueHead->get_active ();
}

bool SaveAsDialog::getToTailOfQueue () {

    return putToQueueTail->get_active ();
}

Glib::ustring SaveAsDialog::getFileName () {

	// fname is empty if the dialog has been cancelled
	if (fname.length())
		return removeExtension(fname) + Glib::ustring(".") + formatOpts->getFormat().format;
	else
		return "";
}

Glib::ustring SaveAsDialog::getDirectory () {

    return fchooser->get_current_folder ();
}

SaveFormat SaveAsDialog::getFormat () {

    return formatOpts->getFormat ();
}

void SaveAsDialog::okPressed () {

    fname = fchooser->get_filename();

    // checking if the filename field is empty. The user have to click Cancel if he don't want to specify a filename
    // NB: There seem to be a bug in Gtkmm2.22 / FileChooserWidget : if you suppress the filename entry and
    //     click on a folder in the list, the filename field is empty but get_filename will return the folder's path :/
    if (!fname.length() || safe_file_test (fname, Glib::FILE_TEST_IS_DIR)) {
        Glib::ustring msg_ = Glib::ustring("<b>") + M("MAIN_MSG_EMPTYFILENAME") + "</b>";
        Gtk::MessageDialog msgd (*this, msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_OK, true);
        msgd.run ();
        return;
    }
    response = Gtk::RESPONSE_OK;
    hide ();
}

void SaveAsDialog::cancelPressed () {

    fname = "";
    response = Gtk::RESPONSE_CANCEL;
    hide ();
}

void SaveAsDialog::formatChanged (Glib::ustring f) {

    if (f=="jpg")
        fchooser->set_filter (filter_jpg);
    else if (f=="png") 
        fchooser->set_filter (filter_png);
    else if (f=="tif")
        fchooser->set_filter (filter_tif);
}

void SaveAsDialog::setInitialFileName (Glib::ustring fname) {

    fchooser->set_current_name(fname);
}
