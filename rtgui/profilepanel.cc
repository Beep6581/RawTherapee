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
#include "profilepanel.h"
#include "options.h"
#include "profilestore.h"
#include "clipboard.h"
#include "multilangmgr.h"
#include "../rtengine/safegtk.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

extern Glib::ustring argv0;

PartialPasteDlg* ProfilePanel::partialProfileDlg;


void ProfilePanel::init () {
    partialProfileDlg = new PartialPasteDlg("Foo");
}

void ProfilePanel::cleanup () {
    delete partialProfileDlg;
}

ProfilePanel::ProfilePanel (bool readOnly) {

    tpc = NULL;
  
    profiles = Gtk::manage (new MyComboBoxText ());
    Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
    hbox->show ();
//    pack_start (*profiles, Gtk::PACK_SHRINK, 4);

    pack_start (*hbox, Gtk::PACK_SHRINK, 4);
    
    load = Gtk::manage (new Gtk::Button ());
    load->add (*Gtk::manage (new RTImage ("gtk-open.png")));
    if (!readOnly) save = Gtk::manage (new Gtk::Button ());
    if (!readOnly) save->add (*Gtk::manage (new RTImage ("gtk-save-large.png")));
    if (!readOnly) copy = Gtk::manage (new Gtk::Button ());
    if (!readOnly) copy->add (*Gtk::manage (new RTImage ("edit-copy.png")));
    paste = Gtk::manage (new Gtk::Button ());
    paste->add (*Gtk::manage (new RTImage ("edit-paste.png")));

    hbox->pack_start (*profiles);
    hbox->pack_start (*load, Gtk::PACK_SHRINK, 1);
    if (!readOnly) hbox->pack_start (*save, Gtk::PACK_SHRINK, 1);
    hbox->pack_start (*copy, Gtk::PACK_SHRINK, 1);
    if (!readOnly) hbox->pack_start (*paste, Gtk::PACK_SHRINK, 1);

    load->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::load_clicked) );
    if (!readOnly) save->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::save_clicked) );
    if (!readOnly) copy->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::copy_clicked) );
    paste->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::paste_clicked) );

    custom = NULL;
    lastsaved = NULL;
    dontupdate = false;

    refreshProfileList ();

    profiles->set_active (0);
    old = profiles->get_active_text();
    changeconn = profiles->signal_changed().connect( sigc::mem_fun(*this, &ProfilePanel::selection_changed) );

    load->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPLOAD"));
    if (!readOnly) save->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPSAVE"));
    if (!readOnly) copy->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPCOPY"));
    paste->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPPASTE"));
    
    show_all_children ();
}

ProfilePanel::~ProfilePanel () {

    if (custom)    { custom->deleteInstance();    delete custom;    }
    if (lastsaved) { lastsaved->deleteInstance(); delete lastsaved; }
}

void ProfilePanel::refreshProfileList () {

    Glib::ustring oldsel = profiles->get_active_text ();
    changeconn.block (true);

    // clear items
    profiles->clear_items ();
    pparams.clear ();

    // re-parse profile directories (deletes old ones)
    profileStore.parseProfiles ();
    pparams = profileStore.getProfileNames ();
    for (unsigned int i=0; i<pparams.size(); i++)
        profiles->append_text (pparams[i]);

    if (custom)
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
    if (lastsaved)
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")");

    profiles->set_active_text (oldsel);
    changeconn.block (false);
}

void ProfilePanel::save_clicked (GdkEventButton* event) {

    if (event->button != 1)
        return;

    Gtk::FileChooserDialog dialog(M("PROFILEPANEL_SAVEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    FileChooserLastFolderPersister persister( &dialog, options.loadSaveProfilePath );

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("PROFILEPANEL_FILEDLGFILTERPP"));
    filter_pp.add_pattern("*"+paramFileExtension);
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("PROFILEPANEL_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

//    dialog.set_do_overwrite_confirmation (true);

    savedialog = &dialog;

    bool done = false;
    do {
        int result = dialog.run();
        dialog.hide();

        if (result==Gtk::RESPONSE_OK) {

            std::string fname = dialog.get_filename();

            bool hasext = true;
	    size_t dotpos = fname.find_last_of ('.');
            if (dotpos==Glib::ustring::npos)
                hasext = false;
	    size_t dirpos1 = fname.find_last_of ('/');
            if (dirpos1!=Glib::ustring::npos && dirpos1>dotpos)
                hasext = false;
	    size_t dirpos2 = fname.find_last_of ('\\');
            if (dirpos2!=Glib::ustring::npos && dirpos2>dotpos)
                hasext = false;

            if (!hasext)
                fname = fname + paramFileExtension;

            if (safe_file_test (fname, Glib::FILE_TEST_EXISTS)) {
                Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
                Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
                int response = msgd.run ();
                if (response==Gtk::RESPONSE_NO)
                    // make another try
                    continue;
            }

            PartialProfile* toSave = NULL;
            if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")")
                toSave = custom;
            else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")")
                toSave = lastsaved;
            else
                toSave = profileStore.getProfile (profiles->get_active_text());

            if (toSave) {
                if (event->state & Gdk::CONTROL_MASK) {
                    // opening the partial paste dialog window
                    partialProfileDlg->set_title(M("PROFILEPANEL_SAVEPPASTE"));
                    int i = partialProfileDlg->run();
                    partialProfileDlg->hide();
                    if (i != Gtk::RESPONSE_OK)
                        return;

                    // saving the partial profile
                    PartialProfile ppTemp(true);
                    partialProfileDlg->applyPaste (ppTemp.pparams, ppTemp.pedited, toSave->pparams, toSave->pedited);
                    ppTemp.pparams->save (fname, "", ppTemp.pedited);
                    ppTemp.deleteInstance();
                }
                else {
                    // saving a full profile
                    toSave->pparams->save (fname);
                }
                refreshProfileList ();
            }
        }
        done = true;
    } while (!done);
    return;
}

/*
 * Copy the actual full profile to the clipboard
 */
void ProfilePanel::copy_clicked (GdkEventButton* event) {

    if (event->button != 1)
        return;

    PartialProfile* toSave = NULL;
    if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")") 
        toSave = custom;
    else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")") 
        toSave = lastsaved; 
    else
        toSave = profileStore.getProfile (profiles->get_active_text());

    // toSave has to be a complete procparams
    if (toSave) {
        if (event->state & Gdk::CONTROL_MASK) {
            // opening the partial paste dialog window
            partialProfileDlg->set_title(M("PROFILEPANEL_COPYPPASTE"));
            int i = partialProfileDlg->run();
            partialProfileDlg->hide();
            if (i != Gtk::RESPONSE_OK)
                return;

            // saving a partial profile
            PartialProfile ppTemp(true);
            partialProfileDlg->applyPaste (ppTemp.pparams, ppTemp.pedited, toSave->pparams, toSave->pedited);
            clipboard.setPartialProfile(ppTemp);
            ppTemp.deleteInstance();
        }
        else
            clipboard.setProcParams (*toSave->pparams);
    }
    return;
}

/*
 * Load a potentially partial profile
 */
void ProfilePanel::load_clicked (GdkEventButton* event) {

    if (event->button != 1)
        return;

    Gtk::FileChooserDialog dialog(M("PROFILEPANEL_LOADDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN);
    FileChooserLastFolderPersister persister( &dialog, options.loadSaveProfilePath );

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-open"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("PROFILEPANEL_FILEDLGFILTERPP"));
    filter_pp.add_pattern("*"+paramFileExtension);
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("PROFILEPANEL_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

    int result = dialog.run();
    dialog.hide();

    if (result==Gtk::RESPONSE_OK) {
        Glib::ustring fname = dialog.get_filename();

        if (event->state & Gdk::CONTROL_MASK) {
            // opening the partial paste dialog window
            partialProfileDlg->set_title(M("PROFILEPANEL_LOADPPASTE"));
            int i = partialProfileDlg->run();
            partialProfileDlg->hide();
            if (i != Gtk::RESPONSE_OK)
                return;
        }
        bool customCreated = false;
        if (!custom) {
            custom = new PartialProfile (true);
            custom->set(true);
            customCreated = true;
        }
        int err = custom->load (fname);
        if (!err) {
            bool prevState = changeconn.block(true);
            Glib::ustring newEntry = Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")";
            profiles->append_text     (newEntry);
            profiles->set_active_text (newEntry);
            old = profiles->get_active_text();
            changeconn.block(prevState);

            if (event->state & Gdk::CONTROL_MASK) {
                // applying partial profile
                PartialProfile ppTemp(true);
                // the 2 next line modify custom->pedited without modifying custom->pparams
                partialProfileDlg->applyPaste (ppTemp.pparams, ppTemp.pedited, custom->pparams, custom->pedited);
                *custom->pedited = *ppTemp.pedited;
                ppTemp.deleteInstance();
            }

            changeTo (custom, M("PROFILEPANEL_PFILE"));
        }
        else if (customCreated) {
            // we delete custom
            custom->deleteInstance();
            delete custom;
        }
    }
    return;
}

/*
 * Paste a full profile from the clipboard
 */
void ProfilePanel::paste_clicked (GdkEventButton* event) {

    if (event->button != 1)
        return;
    if (!clipboard.hasProcParams())
        return;

    if (event->state & Gdk::CONTROL_MASK) {
        partialProfileDlg->set_title(M("PROFILEPANEL_PASTEPPASTE"));
        int i = partialProfileDlg->run();
        partialProfileDlg->hide();
        if (i != Gtk::RESPONSE_OK)
            return;
    }

    bool prevState = changeconn.block(true);
    Glib::ustring newEntry = Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")";

    if (!custom) {
        custom = new PartialProfile (true);
        custom->pedited->set(true);
        profiles->append_text (newEntry);
    }
    ProcParams pp = clipboard.getProcParams ();
    *custom->pparams = pp;

    profiles->set_active_text (newEntry);
    old = profiles->get_active_text();

    changeconn.block(prevState);

    if (event->state & Gdk::CONTROL_MASK) {
        // applying partial profile
        PartialProfile ppTemp(true);
        // the 2 next line modify custom->pedited without modifying custom->pparams
        partialProfileDlg->applyPaste (ppTemp.pparams, ppTemp.pedited, custom->pparams, custom->pedited);
        *custom->pedited = *ppTemp.pedited;
        ppTemp.deleteInstance();
    }

    changeTo (custom, M("HISTORY_FROMCLIPBOARD"));
    return;
}

void ProfilePanel::changeTo (PartialProfile* newpp, Glib::ustring profname) {

    if (!newpp)
        return;

    if (tpc)
        tpc->profileChange (newpp, EvProfileChanged, profname);  
}

void ProfilePanel::selection_changed () {

    Glib::ustring entry;
    if (profiles->get_active_text() == (entry = Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")")) {
        if (!dontupdate) 
            changeTo (custom, entry);
    }
    else if (profiles->get_active_text() == (entry = Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")"))
            changeTo (lastsaved, entry);
    else {
    	PartialProfile* s = profileStore.getProfile (profiles->get_active_text());
        if (s)
            changeTo (s, profiles->get_active_text());
    }
    old = profiles->get_active_text ();
    dontupdate = false;
}

void ProfilePanel::procParamsChanged (rtengine::procparams::ProcParams* p, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited) {

    // to prevent recursion, filter out the events caused by the profilepanel
    if (ev==EvProfileChanged || ev==EvPhotoLoaded)
        return;

    Glib::ustring entry = Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")";
    if (profiles->get_active_text() != entry) {
        dontupdate = true;
        if (!custom) {
            custom = new PartialProfile (true);
            custom->set(true);
            profiles->append_text (entry);
        }
        profiles->set_active_text (entry);
        old = profiles->get_active_text();
    }
    *custom->pparams = *p;
}

void ProfilePanel::initProfile (const Glib::ustring& profname, ProcParams* lastSaved) {

    changeconn.block (true);

    profiles->clear_items ();
    pparams.clear ();

    pparams = profileStore.getProfileNames ();
    for (unsigned int i=0; i<pparams.size(); i++)
        profiles->append_text (pparams[i]);

    if (custom) {
        custom->deleteInstance();
        delete custom; custom = NULL;
    }

    if (lastsaved) {
        lastsaved->deleteInstance();
        delete lastsaved; lastsaved = NULL;
    }
    if (lastSaved) {
        ParamsEdited* pe = new ParamsEdited();
        pe->set(true);
        lastsaved = new PartialProfile(lastSaved, pe);
    }

    Glib::ustring defline = profname;
    PartialProfile* defprofile = profileStore.getProfile (profname);

    if (lastsaved) {
        defline = Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")";
        defprofile = lastsaved;
        profiles->append_text (defline);
    }

    if (tpc) {
        if (lastsaved)
            tpc->setDefaults (lastsaved->pparams);
        else
            tpc->setDefaults (profileStore.getProfile (profname)->pparams);
    }
    if (defprofile) {
        old = defline;
        profiles->set_active_text (defline);
        changeconn.block (false);
        if (tpc)
            tpc->profileChange (defprofile, EvPhotoLoaded, defline);  
    }
    else {
        bool dels = false;
        // select first valid profile
        old = "";
        profiles->set_active (0);
        PartialProfile* s = profileStore.getProfile (profiles->get_active_text());
        if (!s) {
            changeconn.block (false);
            s = new PartialProfile (true);
            s->set(true);
            dels = true;  // we've created a temporary PartialProfile, so we set a flag to destroy it
            if (tpc)
                tpc->profileChange (s, EvPhotoLoaded, DEFPROFILE_INTERNAL);
        }
        else {
            Glib::ustring cProfile = profiles->get_active_text();
            changeconn.block (false);
            if (tpc)
                tpc->profileChange (s, EvPhotoLoaded, cProfile);
        }

        if (dels) {
            s->deleteInstance();
            delete s;
        }
    }
}



