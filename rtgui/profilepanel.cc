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

PartialPasteDlg* ProfilePanel::partialProfileDlg;


void ProfilePanel::init () {
    partialProfileDlg = new PartialPasteDlg("Foo");
}

void ProfilePanel::cleanup () {
    delete partialProfileDlg;
}

ProfilePanel::ProfilePanel (bool readOnly) : storedPProfile(NULL), lastFilename(""), imagePath("") {

    tpc = NULL;

    profileFillModeOnImage  = new RTImage("profile-filled.png");
    profileFillModeOffImage = new RTImage("profile-partial.png");
    fillMode = Gtk::manage (new Gtk::ToggleButton());
    fillMode->set_active(options.filledProfile);
    fillMode->add( options.filledProfile ? *profileFillModeOnImage : *profileFillModeOffImage );
    fillMode->signal_toggled().connect ( sigc::mem_fun(*this, &ProfilePanel::profileFillModeToggled) );
    fillMode->set_tooltip_text(M("PROFILEPANEL_MODE_TIP"));

    // Create the Combobox
    profiles = Gtk::manage (new ProfileStoreComboBox ());

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

    hbox->pack_start (*fillMode, Gtk::PACK_SHRINK, 1);
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

    profileStore.addListener(this);

    changeconn = profiles->signal_changed().connect( sigc::mem_fun(*this, &ProfilePanel::selection_changed) );

    load->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPLOAD"));
    if (!readOnly) save->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPSAVE"));
    if (!readOnly) copy->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPCOPY"));
    paste->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPPASTE"));
    
    show_all_children ();
}

ProfilePanel::~ProfilePanel () {

    profileStore.removeListener(this);
    if (custom)    { custom->deleteInstance();    delete custom;    }
    if (lastsaved) { lastsaved->deleteInstance(); delete lastsaved; }
    delete profileFillModeOnImage;
    delete profileFillModeOffImage;
}

bool ProfilePanel::isCustomSelected() {
    if (profiles->getCurrentLabel() == Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")"))
        return true;
    return false;
}

bool ProfilePanel::isLastSavedSelected() {
    if (profiles->getCurrentLabel() == Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")"))
        return true;
    return false;
}

Gtk::TreeIter ProfilePanel::getCustomRow() {
    Gtk::TreeIter row;
    if (custom)
        row = profiles->getRowFromLabel(Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")"));
    return row;
}

Gtk::TreeIter ProfilePanel::getLastSavedRow() {
    Gtk::TreeIter row;
    if (lastsaved) {
        row = profiles->getRowFromLabel(Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")"));
    }
    return row;
}

Gtk::TreeIter ProfilePanel::addCustomRow() {
    const ProfileStoreEntry *customPSE = new ProfileStoreEntry(Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")"), PSET_FILE, 0, 0);
    Gtk::TreeIter newEntry = profiles->addRow(customPSE);
    return newEntry;
}

Gtk::TreeIter ProfilePanel::addLastSavedRow() {
    const ProfileStoreEntry *lastSavedPSE = new ProfileStoreEntry(Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")"), PSET_FILE, 0, 0);
    Gtk::TreeIter newEntry = profiles->addRow(lastSavedPSE);
    return newEntry;
}

void ProfilePanel::storeCurrentValue () {
    // TODO: Find a way to get and restore the current selection; the following line can't work anymore
    storedValue = profiles->getFullPathFromActiveRow();
    if (!isCustomSelected() && !isLastSavedSelected()) {
        // storing the current entry's procparams, if not "Custom" or "LastSaved"

        // for now, the storedPProfile has default internal values
        const ProfileStoreEntry *entry = profiles->getSelectedEntry();
        const PartialProfile *currProfile;
        if (entry && (currProfile = profileStore.getProfile(entry))!=NULL) {
            // now storedPProfile has the current entry's values
            storedPProfile = new PartialProfile(currProfile->pparams, currProfile->pedited, true);
        }
        else
            storedPProfile = new PartialProfile(true);
    }
}

/* Get the ProfileStore's entry list and recreate the combobox entries
 * If you want want to update the ProfileStore list itself (rescan the dir tree), use its "parseProfiles" method instead
 */
void ProfilePanel::updateProfileList () {

    bool ccPrevState = changeconn.block(true);

    // rescan file tree
    profiles->updateProfileList();

    if (custom)
        addCustomRow();

    if (lastsaved)
        addLastSavedRow();

    changeconn.block (ccPrevState);
}

void ProfilePanel::restoreValue () {
    bool ccPrevState = changeconn.block(true);

    if (!profiles->setActiveRowFromFullPath(storedValue) && storedPProfile) {
        if (custom) delete custom;
        custom = new PartialProfile (storedPProfile->pparams, storedPProfile->pedited, true);
        Gtk::TreeIter custRow = getCustomRow();
        if (custRow)
            profiles->set_active(custRow);
        else
            profiles->set_active (addCustomRow());
    }

    currRow = profiles->get_active();

    changeconn.block (ccPrevState);

    storedValue = "";

    if (storedPProfile) {
        storedPProfile->deleteInstance();
        delete storedPProfile;
        storedPProfile = NULL;
    }
}

void ProfilePanel::save_clicked (GdkEventButton* event) {

    if (event->button != 1)
        return;

    Gtk::FileChooserDialog dialog(M("PROFILEPANEL_SAVEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    FileChooserLastFolderPersister persister( &dialog, options.loadSaveProfilePath );
    dialog.set_current_name (lastFilename);

    //Add the user's default (or global if multiuser=false) profile path to the Shortcut list
#ifdef WIN32
    // Dirty workaround, waiting for a clean solution by using exceptions!
    if (!safe_is_shortcut_dir(options.getPreferredProfilePath()))
#endif
    try {
        dialog.add_shortcut_folder(options.getPreferredProfilePath());
    }
    catch (Glib::Error &err) {}
    //Add the image's path to the Shortcut list
#ifdef WIN32
    // Dirty workaround, waiting for a clean solution by using exceptions!
    if (!safe_is_shortcut_dir(imagePath))
#endif
    try {
        dialog.add_shortcut_folder(imagePath);
    }
    catch (Glib::Error &err) {}

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("FILECHOOSER_FILTER_PP"));
    filter_pp.add_pattern("*"+paramFileExtension);
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

//    dialog.set_do_overwrite_confirmation (true);

    bool done = false;
    do {
        if (dialog.run()==Gtk::RESPONSE_OK) {

            std::string fname = dialog.get_filename();
            Glib::ustring ext = getExtension (fname);

            if (("." + ext) != paramFileExtension)
                fname += paramFileExtension;

            if (!confirmOverwrite (dialog, fname))
                continue;

            lastFilename = Glib::path_get_basename (fname);

            const PartialProfile* toSave;
            if (isCustomSelected())
                toSave = custom;
            else if (isLastSavedSelected())
                toSave = lastsaved;
            else {
                const ProfileStoreEntry* entry = profiles->getSelectedEntry();
                toSave = entry ? profileStore.getProfile (profiles->getSelectedEntry()) : NULL;
            }

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
                    int retCode = ppTemp.pparams->save (fname, "", true, ppTemp.pedited);
                    ppTemp.deleteInstance();
                    if (retCode)
                        writeFailed(dialog, fname);
                    else {
                        done=true;
                        bool ccPrevState = changeconn.block(true);
                        profileStore.parseProfiles();
                        changeconn.block (ccPrevState);
                    }
                }
                else {
                    // saving a full profile
                    int retCode = toSave->pparams->save (fname);
                    if (retCode)
                        writeFailed(dialog, fname);
                    else {
                        done=true;
                        bool ccPrevState = changeconn.block(true);
                        profileStore.parseProfiles();
                        changeconn.block (ccPrevState);
                    }
                }
            }
            else done = true;
        }
        else done = true;
    } while (!done);
    return;
}

/*
 * Copy the actual full profile to the clipboard
 */
void ProfilePanel::copy_clicked (GdkEventButton* event) {

    if (event->button != 1)
        return;

    const PartialProfile* toSave;
    if (isCustomSelected())
        toSave = custom;
    else if (isLastSavedSelected())
        toSave = lastsaved; 
    else {
        const ProfileStoreEntry* entry = profiles->getSelectedEntry();
        toSave = entry ? profileStore.getProfile (entry) : NULL;
    }

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

    //Add the user's default (or global if multiuser=false) profile path to the Shortcut list
#ifdef WIN32
    // Dirty workaround, waiting for a clean solution by using exceptions!
    if (!safe_is_shortcut_dir(options.getPreferredProfilePath()))
#endif
    try {
        dialog.add_shortcut_folder(options.getPreferredProfilePath());
    }
    catch (Glib::Error &err) {}

    //Add the image's path to the Shortcut list
#ifdef WIN32
    // Dirty workaround, waiting for a clean solution by using exceptions!
    if (!safe_is_shortcut_dir(imagePath))
#endif
    try {
        dialog.add_shortcut_folder(imagePath);
    }
    catch (Glib::Error &err) {}

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-open"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("FILECHOOSER_FILTER_PP"));
    filter_pp.add_pattern("*"+paramFileExtension);
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("FILECHOOSER_FILTER_ANY"));
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
            customCreated = true;
        }

        ProcParams pp;
        ParamsEdited pe;
        int err = pp.load (fname, &pe);
        if (!err) {
            if (!customCreated && fillMode->get_active())
                custom->pparams->setDefaults();
            custom->set(true);

            bool prevState = changeconn.block(true);
            Gtk::TreeIter newEntry = addCustomRow();
            profiles->set_active (newEntry);
            currRow = profiles->get_active();
            changeconn.block(prevState);

            // Now we have procparams initialized to default if fillMode is on
            // and paramsedited initialized to default in all cases

            if (event->state & Gdk::CONTROL_MASK)
                // custom.pparams = loadedFile.pparams filtered by ( loadedFile.pedited & partialPaste.pedited )
                partialProfileDlg->applyPaste (custom->pparams, !fillMode->get_active()?custom->pedited:NULL, &pp, &pe);
            else {
                // custom.pparams = loadedFile.pparams filtered by ( loadedFile.pedited )
                pe.combine(*custom->pparams, pp, true);
                if (!fillMode->get_active())
                    *custom->pedited = pe;
            }

            changeTo (custom, M("PROFILEPANEL_PFILE"));
        }
        else if (customCreated) {
            // we delete custom
            custom->deleteInstance();
            delete custom; custom = NULL;
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

    if (!custom) {
        custom = new PartialProfile (true);
        if (isLastSavedSelected()) {
            *custom->pparams = *lastsaved->pparams;
        }
        else {
            const ProfileStoreEntry* entry = profiles->getSelectedEntry();
            if (entry) {
                const PartialProfile* partProfile = profileStore.getProfile (entry);
                *custom->pparams = *partProfile->pparams;
            }
        }
        profiles->set_active (addCustomRow());
        currRow = profiles->get_active();
    }
    else {
        if (fillMode->get_active())
            custom->pparams->setDefaults();
        profiles->set_active(getCustomRow());
        currRow = profiles->get_active();
    }
    custom->pedited->set(true);

    changeconn.block(prevState);

    // Now we have procparams initialized to default if fillMode is on
    // and paramsedited initialized to default in all cases

    ProcParams pp = clipboard.getProcParams ();
    if (clipboard.hasPEdited()) {
        ParamsEdited pe = clipboard.getParamsEdited();
        if (event->state & Gdk::CONTROL_MASK)
            // custom.pparams = clipboard.pparams filtered by ( clipboard.pedited & partialPaste.pedited )
            partialProfileDlg->applyPaste (custom->pparams, !fillMode->get_active()?custom->pedited:NULL, &pp, &pe);
        else {
            // custom.pparams = clipboard.pparams filtered by ( clipboard.pedited )
            pe.combine(*custom->pparams, pp, true);
            if (!fillMode->get_active())
                *custom->pedited = pe;
        }
    }
    else {
        if (event->state & Gdk::CONTROL_MASK)
            // custom.pparams = clipboard.pparams filtered by ( partialPaste.pedited )
            partialProfileDlg->applyPaste (custom->pparams, NULL, &pp, NULL);
        else {
            // custom.pparams = clipboard.pparams non filtered
            *custom->pparams = pp;
        }
    }

    changeTo (custom, M("HISTORY_FROMCLIPBOARD"));
    return;
}

void ProfilePanel::changeTo (const PartialProfile* newpp, Glib::ustring profname) {

    if (!newpp)
        return;

    if (tpc)
        tpc->profileChange (newpp, EvProfileChanged, profname);  
}

void ProfilePanel::selection_changed () {

    if (isCustomSelected()) {
        if (!dontupdate) 
            changeTo (custom, Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")"));
    }
    else if (isLastSavedSelected())
            changeTo (lastsaved, Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")"));
    else {
        const ProfileStoreEntry *pse = profiles->getSelectedEntry();
        if (pse->type == PSET_FOLDER) {
            // this entry is invalid, restoring the old value
            bool ccPrevState = changeconn.block(true);
            profiles->set_active(currRow);
            changeconn.block(ccPrevState);
            dontupdate = false;
            return;
        }
        else
            currRow = profiles->get_active();

        const PartialProfile* s = profileStore.getProfile (pse);
        if (s) {
            if (fillMode->get_active() && s->pedited) {
                ParamsEdited pe(true);
                PartialProfile s2(s->pparams, &pe, false);
                changeTo (&s2, pse->label+"+");
            }
            else
                changeTo (s, pse->label);
        }
    }
    dontupdate = false;
}

void ProfilePanel::procParamsChanged (rtengine::procparams::ProcParams* p, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited) {

    // to prevent recursion, filter out the events caused by the profilepanel
    if (ev==EvProfileChanged || ev==EvPhotoLoaded)
        return;

    if (!isCustomSelected()) {
        dontupdate = true;
        if (!custom) {
            custom = new PartialProfile (true);
            custom->set(true);
            profiles->set_active (addCustomRow());
            currRow = profiles->get_active();
        }
        else {
            profiles->set_active(getCustomRow());
            currRow = profiles->get_active();
        }
    }
    *custom->pparams = *p;
}

/** @brief Initialize the Profile panel with a default profile, overridden by the last saved profile if provided
 *
 * The file tree has already been created on object's construction. We add here the Custom, LastSaved and/or Internal item.
 *
 * @param profileFullPath   full path of the profile; must start by the virtual root (${G} or ${U}, and without suffix
 * @param lastSaved         pointer to the last saved ProcParam; may be NULL
 */
void ProfilePanel::initProfile (const Glib::ustring& profileFullPath, ProcParams* lastSaved) {

    const ProfileStoreEntry *pse = NULL;
    const PartialProfile *defprofile = NULL;

    bool ccPrevState = changeconn.block(true);

    if (custom) {
        custom->deleteInstance();
        delete custom; custom = NULL;
    }

    if (lastsaved) {
        lastsaved->deleteInstance();
        delete lastsaved; lastsaved = NULL;
    }
    if (lastSaved) {
        ParamsEdited* pe = new ParamsEdited(true);
        // copying the provided last saved profile to ProfilePanel::lastsaved
        lastsaved = new PartialProfile(lastSaved, pe);
    }

    // update the content of the combobox; will add 'custom' and 'lastSaved' if necessary
    updateProfileList();

    Gtk::TreeIter lasSavedEntry;
    // adding the Last Saved combobox entry, if needed
    if (lastsaved) {
        defprofile = lastsaved;
        lasSavedEntry = getLastSavedRow();
    }

    if (!(pse = profileStore.findEntryFromFullPath(profileFullPath))) {
        // entry not found, pse = the Internal ProfileStoreEntry
        pse = profileStore.getInternalDefaultPSE();
    }

    defprofile = profileStore.getProfile (pse);

    // selecting the "Internal" entry
    profiles->setInternalEntry ();
    currRow = profiles->get_active();

    if (lastsaved) {
        if (lasSavedEntry) profiles->set_active (lasSavedEntry);
        currRow = profiles->get_active();
        if (tpc) {
            tpc->setDefaults   (lastsaved->pparams);
            tpc->profileChange (lastsaved, EvPhotoLoaded, profiles->getSelectedEntry()->label);
        }
    }
    else {
        if (pse) {
            profiles->setActiveRowFromEntry(pse);
            currRow = profiles->get_active();
        }
        if (tpc) {
            tpc->setDefaults   (defprofile->pparams);
            tpc->profileChange (defprofile, EvPhotoLoaded, profiles->getSelectedEntry()->label);
        }
    }
    changeconn.block (ccPrevState);
}

void ProfilePanel::setInitialFileName (const Glib::ustring& filename) {
    lastFilename = Glib::path_get_basename(filename) + paramFileExtension;
    imagePath = Glib::path_get_dirname(filename);
}

void ProfilePanel::profileFillModeToggled() {
    if (fillMode->get_active()) {
        // The button is pressed, we'll use the profileFillModeOnImage
        fillMode->set_image(*profileFillModeOnImage);
    }
    else {
        // The button is released, we'll use the profileFillModeOffImage
        fillMode->set_image(*profileFillModeOffImage);
    }
}

void ProfilePanel::writeOptions() {
    options.filledProfile = fillMode->get_active();
}

