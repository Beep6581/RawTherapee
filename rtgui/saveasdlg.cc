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

extern Options options;
SaveAsDialog::SaveAsDialog (Glib::ustring initialDir) {

    Gtk::VBox* vbox = get_vbox ();

    fchooser = new Gtk::FileChooserWidget (Gtk::FILE_CHOOSER_ACTION_SAVE);
    fchooser->set_current_folder (initialDir);

    filter_jpg.set_name(M("SAVEDLG_JPGFILTER"));
    filter_jpg.add_pattern("*.jpg");

    filter_tif.set_name(M("SAVEDLG_JPGFILTER"));
    filter_tif.add_pattern("*.tif");

    filter_png.set_name(M("SAVEDLG_JPGFILTER"));
    filter_png.add_pattern("*.png");

    vbox->pack_start (*fchooser);

    Gtk::HSeparator* hsep1 = new Gtk::HSeparator ();
    vbox->pack_start (*hsep1, Gtk::PACK_SHRINK, 2);

// Output Options
// ~~~~~~~~~~~~~~
    formatOpts = new SaveFormatPanel ();
    formatOpts->init (options.saveFormat);
    formatOpts->setListener (this);

    vbox->pack_start (*formatOpts, Gtk::PACK_SHRINK, 4);

    Gtk::HSeparator* hsep2 = new Gtk::HSeparator ();
    vbox->pack_start (*hsep2, Gtk::PACK_SHRINK, 2);

// queue/immediate
// ~~~~~~~~~~~~~
    immediately    = new Gtk::RadioButton (M("SAVEDLG_SAVEIMMEDIATELY"));
    putToQueueHead = new Gtk::RadioButton (M("SAVEDLG_PUTTOQUEUEHEAD"));
    putToQueueTail = new Gtk::RadioButton (M("SAVEDLG_PUTTOQUEUETAIL"));
    vbox->pack_start (*immediately, Gtk::PACK_SHRINK, 4);
    vbox->pack_start (*putToQueueHead, Gtk::PACK_SHRINK, 4);
    vbox->pack_start (*putToQueueTail, Gtk::PACK_SHRINK, 4);
    immediately->set_active (true);
    Gtk::RadioButton::Group g = immediately->get_group();
    putToQueueHead->set_group (g);
    putToQueueTail->set_group (g);
        
// buttons
// ~~~~~~    
    Gtk::Button* ok     = new Gtk::Button (M("GENERAL_OK"));
    Gtk::Button* cancel = new Gtk::Button (M("GENERAL_CANCEL"));

    ok->set_image (*(new Gtk::Image (Gtk::StockID("gtk-ok"), Gtk::ICON_SIZE_BUTTON)));
    cancel->set_image (*(new Gtk::Image (Gtk::StockID("gtk-cancel"), Gtk::ICON_SIZE_BUTTON)));

    ok->signal_clicked().connect( sigc::mem_fun(*this, &SaveAsDialog::okPressed) );
    cancel->signal_clicked().connect( sigc::mem_fun(*this, &SaveAsDialog::cancelPressed) );

    get_action_area()->pack_end (*ok, Gtk::PACK_SHRINK, 4);
    get_action_area()->pack_end (*cancel, Gtk::PACK_SHRINK, 4);

    set_border_width (4);
    show_all_children ();
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

    return fname;
}

Glib::ustring SaveAsDialog::getDirectory () {

    return fchooser->get_current_folder ();
}

SaveFormat SaveAsDialog::getFormat () {

    return formatOpts->getFormat ();
}

void SaveAsDialog::okPressed () {

    fname = fchooser->get_filename();
    hide ();
}

void SaveAsDialog::cancelPressed () {

    fname = "";
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
