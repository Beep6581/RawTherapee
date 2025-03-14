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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "profilepanel.h"

#include "clipboard.h"
#include "multilangmgr.h"
#include "options.h"
#include "profilestorecombobox.h"
#include "paramsedited.h"
#include "pathutils.h"
#include "rtimage.h"

#include "rtengine/procparams.h"
#include "rtengine/procevents.h"

using namespace rtengine;
using namespace rtengine::procparams;

PartialPasteDlg* ProfilePanel::partialProfileDlg = nullptr;
Gtk::Window* ProfilePanel::parent;

void ProfilePanel::init (Gtk::Window* parentWindow)
{
    parent = parentWindow;
}

void ProfilePanel::cleanup ()
{
    delete partialProfileDlg;
}

ProfilePanel::ProfilePanel () : storedPProfile(nullptr),
    modeOn("profile-filled"), modeOff("profile-partial"),
    profileFillImage(Gtk::manage(new RTImage(options.filledProfile ? modeOn : modeOff, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    lastSavedPSE(nullptr), customPSE(nullptr)
{
    tpc = nullptr;

    fillMode = Gtk::manage (new Gtk::ToggleButton());
    fillMode->set_active(options.filledProfile);
    fillMode->add(*profileFillImage);
    fillMode->signal_toggled().connect ( sigc::mem_fun(*this, &ProfilePanel::profileFillModeToggled) );
    fillMode->set_tooltip_text(M("PROFILEPANEL_MODE_TOOLTIP"));
//GTK318
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
    fillMode->set_margin_right(2);
#endif
//GTK318
    setExpandAlignProperties(fillMode, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);

    // Create the Combobox
    profiles = Gtk::manage (new ProfileStoreComboBox ());
    setExpandAlignProperties(profiles, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    load = Gtk::manage (new Gtk::Button ());
    load->add (*Gtk::manage (new RTImage ("folder-open", Gtk::ICON_SIZE_LARGE_TOOLBAR)));
    load->get_style_context()->add_class("Left");
    load->set_margin_left(2);
    setExpandAlignProperties(load, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    save = Gtk::manage (new Gtk::Button ());
    save->add (*Gtk::manage (new RTImage ("save", Gtk::ICON_SIZE_LARGE_TOOLBAR)));
    save->get_style_context()->add_class("MiddleH");
    setExpandAlignProperties(save, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    copy = Gtk::manage (new Gtk::Button ());
    copy->add (*Gtk::manage (new RTImage ("copy", Gtk::ICON_SIZE_LARGE_TOOLBAR)));
    copy->get_style_context()->add_class("MiddleH");
    setExpandAlignProperties(copy, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    paste = Gtk::manage (new Gtk::Button ());
    paste->add (*Gtk::manage (new RTImage ("paste", Gtk::ICON_SIZE_LARGE_TOOLBAR)));
    paste->get_style_context()->add_class("Right");
    setExpandAlignProperties(paste, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);

    attach_next_to (*fillMode, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*profiles, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*load, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*save, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*copy, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*paste, Gtk::POS_RIGHT, 1, 1);

    setExpandAlignProperties(this, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    show ();

    load->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::load_clicked) );
    save->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::save_clicked) );
    copy->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::copy_clicked) );
    paste->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &ProfilePanel::paste_clicked) );

    custom = nullptr;
    lastsaved = nullptr;
    dontupdate = false;

    ProfileStore::getInstance()->addListener(this);

    changeconn = profiles->signal_changed().connect( sigc::mem_fun(*this, &ProfilePanel::selection_changed) );

    load->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPLOAD"));
    save->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPSAVE"));
    copy->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPCOPY"));
    paste->set_tooltip_markup (M("PROFILEPANEL_TOOLTIPPASTE"));

    show_all_children ();
}

ProfilePanel::~ProfilePanel ()
{

    ProfileStore::getInstance()->removeListener(this);

    if (custom)    {
        custom->deleteInstance();
        delete custom;
    }

    if (lastsaved) {
        lastsaved->deleteInstance();
        delete lastsaved;
    }

    delete lastSavedPSE;
    delete customPSE;
}

bool ProfilePanel::isCustomSelected()
{
    return profiles->getCurrentLabel().collate_key() == Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")").collate_key();
}

bool ProfilePanel::isLastSavedSelected()
{
    return profiles->getCurrentLabel().collate_key() == Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")").collate_key();
}

Gtk::TreeIter ProfilePanel::getCustomRow()
{
    Gtk::TreeIter row;

    if (custom) {
        row = profiles->getRowFromLabel(Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")"));
    }

    return row;
}

Gtk::TreeIter ProfilePanel::getLastSavedRow()
{
    Gtk::TreeIter row;

    if (lastsaved) {
        row = profiles->getRowFromLabel(Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")"));
    }

    return row;
}

Gtk::TreeIter ProfilePanel::addCustomRow()
{
    if(customPSE) {
        profiles->deleteRow(customPSE);
        delete customPSE;
        customPSE = nullptr;
    }

    customPSE = new ProfileStoreEntry(Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")"), PSET_FILE, 0, 0);
    Gtk::TreeIter newEntry = profiles->addRow(customPSE);
    return newEntry;
}

Gtk::TreeIter ProfilePanel::addLastSavedRow()
{
    if(lastSavedPSE) {
        profiles->deleteRow(lastSavedPSE);
        delete lastSavedPSE;
        lastSavedPSE = nullptr;
    }

    lastSavedPSE = new ProfileStoreEntry(Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")"), PSET_FILE, 0, 0);
    Gtk::TreeIter newEntry = profiles->addRow(lastSavedPSE);
    return newEntry;
}

void ProfilePanel::storeCurrentValue ()
{
    // TODO: Find a way to get and restore the current selection; the following line can't work anymore
    storedValue = profiles->getFullPathFromActiveRow();

    if (!isCustomSelected() && !isLastSavedSelected()) {
        // storing the current entry's procparams, if not "Custom" or "LastSaved"

        // for now, the storedPProfile has default internal values
        const ProfileStoreEntry *entry = profiles->getSelectedEntry();
        const PartialProfile *currProfile;

        if (entry && (currProfile = ProfileStore::getInstance()->getProfile(entry)) != nullptr) {
            // now storedPProfile has the current entry's values
            storedPProfile = new PartialProfile(currProfile->pparams, currProfile->pedited, true);
        } else {
            storedPProfile = new PartialProfile(true);
        }
    }
}

/* Get the ProfileStore's entry list and recreate the combobox entries
 * If you want want to update the ProfileStore list itself (rescan the dir tree), use its "parseProfiles" method instead
 */
void ProfilePanel::updateProfileList ()
{

    bool ccPrevState = changeconn.block(true);

    // rescan file tree
    profiles->updateProfileList();

    if (custom) {
        addCustomRow();
    }

    if (lastsaved) {
        addLastSavedRow();
    }

    changeconn.block (ccPrevState);
}

void ProfilePanel::restoreValue ()
{
    bool ccPrevState = changeconn.block(true);

    if (!profiles->setActiveRowFromFullPath(storedValue) && storedPProfile) {
        if (custom) {
            delete custom;
        }

        custom = new PartialProfile (storedPProfile->pparams, storedPProfile->pedited, true);
        Gtk::TreeIter custRow = getCustomRow();

        if (custRow) {
            profiles->set_active(custRow);
        } else {
            profiles->set_active (addCustomRow());
        }
    }

    currRow = profiles->get_active();

    changeconn.block (ccPrevState);

    storedValue = "";

    if (storedPProfile) {
        storedPProfile->deleteInstance();
        delete storedPProfile;
        storedPProfile = nullptr;
    }
}

void ProfilePanel::save_clicked (GdkEventButton* event)
{
    if (event->button != 1) {
        return;
    }

    const PartialProfile* toSave;

    if (isCustomSelected()) {
        toSave = custom;
    } else if (isLastSavedSelected()) {
        toSave = lastsaved;
    } else {
        const auto entry = profiles->getSelectedEntry();
        toSave = entry ? ProfileStore::getInstance()->getProfile(entry) : nullptr;
    }

    // If no entry has been selected or anything unpredictable happened, toSave
    // can be nullptr.
    if (toSave == nullptr) {
        return;
    }

    // If it's a partial profile, it's more intuitive to first allow the user
    // to choose which parameters to save before showing the Save As dialog
    // #5491
    const auto isPartial = event->state & Gdk::CONTROL_MASK;
    if (isPartial) {
        if (!partialProfileDlg) {
            partialProfileDlg = new PartialPasteDlg(Glib::ustring(), parent);
        }

        partialProfileDlg->set_title(M("PROFILEPANEL_SAVEPPASTE"));
        partialProfileDlg->updateSpotWidget(toSave->pparams);
        const auto response = partialProfileDlg->run();
        partialProfileDlg->hide();

        if (response != Gtk::RESPONSE_OK) {
            return;
        }
    }

    Gtk::FileChooserDialog dialog(getToplevelWindow(this), M("PROFILEPANEL_SAVEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    bindCurrentFolder(dialog, options.loadSaveProfilePath);
    dialog.set_current_name(lastFilename);

    //Add the user's default (or global if multiuser=false) profile path to the Shortcut list
    try {
        dialog.add_shortcut_folder(options.getPreferredProfilePath());
    } catch (Glib::Error&) {}

    //Add the image's path to the Shortcut list
    try {
        dialog.add_shortcut_folder(imagePath);
    } catch (Glib::Error&) {}

    //Add response buttons to the dialog:
    dialog.add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(M("GENERAL_SAVE"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:
    auto filter_pp = Gtk::FileFilter::create();
    filter_pp->set_name(M("FILECHOOSER_FILTER_PP"));
    filter_pp->add_pattern("*" + paramFileExtension);
    dialog.add_filter(filter_pp);

    auto filter_any = Gtk::FileFilter::create();
    filter_any->set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any->add_pattern("*");
    dialog.add_filter(filter_any);

    while (true) {

        // Run the saving dialog and let the user select a path and filename.
        const auto response = dialog.run();

        if (response != Gtk::RESPONSE_OK) {
            // Just exit the loop, cause the user cancels the dialog.
            
            break;
        } else {
            // Go on with saving the the profile.

            auto fname = dialog.get_filename();

            if (("." + getExtension(fname)) != paramFileExtension) {
                fname += paramFileExtension;
            }

            if (!confirmOverwrite(dialog, fname)) {

                // The user doesn't want to override the existing file. So, just restart the loop,
                // so the user can select a different path or file name.
                continue;
            }

            lastFilename = Glib::path_get_basename(fname);

            auto retCode = -1;

            if (isPartial) {
                // Build partial profile
                PartialProfile ppTemp(true);
                partialProfileDlg->applyPaste(ppTemp.pparams, ppTemp.pedited, toSave->pparams, nullptr);
                
                // Save partial profile
                retCode = ppTemp.pparams->save(fname, "", true, ppTemp.pedited);
                
                // Cleanup
                ppTemp.deleteInstance();
            } else {
                // Save full profile
                retCode = toSave->pparams->save(fname);
            }

            if (retCode == 0) {
                // Saving the profile was successfull.

                const auto ccPrevState = changeconn.block(true);
                ProfileStore::getInstance()->parseProfiles();
                changeconn.block(ccPrevState);

                // Because saving has been successfull, just leave the loop;
                break;
            } else {
                // Saving the profile was not successfull.

                writeFailed(dialog, fname);

                // In case the saving process was not successfull (missing permissions, ...)
                // reopen the dialog and try again.
                continue;
            }
        }
    }
}

/*
 * Copy the actual full profile to the clipboard
 */
void ProfilePanel::copy_clicked (GdkEventButton* event)
{

    if (event->button != 1) {
        return;
    }

    const PartialProfile* toSave;

    if (isCustomSelected()) {
        toSave = custom;
    } else if (isLastSavedSelected()) {
        toSave = lastsaved;
    } else {
        const ProfileStoreEntry* entry = profiles->getSelectedEntry();
        toSave = entry ? ProfileStore::getInstance()->getProfile (entry) : nullptr;
    }

    // toSave has to be a complete procparams
    if (toSave) {
        if (event->state & Gdk::CONTROL_MASK) {
            // opening the partial paste dialog window
            if(!partialProfileDlg) {
                partialProfileDlg = new PartialPasteDlg (Glib::ustring (), parent);
            }
            partialProfileDlg->set_title(M("PROFILEPANEL_COPYPPASTE"));
            partialProfileDlg->updateSpotWidget(toSave->pparams);
            int i = partialProfileDlg->run();
            partialProfileDlg->hide();

            if (i != Gtk::RESPONSE_OK) {
                return;
            }

            // saving a partial profile
            PartialProfile ppTemp(true);
            partialProfileDlg->applyPaste (ppTemp.pparams, ppTemp.pedited, toSave->pparams, toSave->pedited);
            clipboard.setPartialProfile(ppTemp);
            ppTemp.deleteInstance();
        } else {
            clipboard.setProcParams (*toSave->pparams);
        }
    }

    return;
}

/*
 * Load a potentially partial profile
 */
void ProfilePanel::load_clicked (GdkEventButton* event)
{

    if (event->button != 1) {
        return;
    }

    Gtk::FileChooserDialog dialog (getToplevelWindow (this), M("PROFILEPANEL_LOADDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN);
    bindCurrentFolder (dialog, options.loadSaveProfilePath);

    //Add the user's default (or global if multiuser=false) profile path to the Shortcut list
    try {
        dialog.add_shortcut_folder(options.getPreferredProfilePath());
    } catch (Glib::Error&) {}

    //Add the image's path to the Shortcut list
    try {
        dialog.add_shortcut_folder(imagePath);
    } catch (Glib::Error&) {}

    //Add response buttons to the dialog:
    dialog.add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(M("GENERAL_OPEN"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:
    Glib::RefPtr<Gtk::FileFilter> filter_pp = Gtk::FileFilter::create();
    filter_pp->set_name(M("FILECHOOSER_FILTER_PP"));
    filter_pp->add_pattern("*" + paramFileExtension);
    dialog.add_filter(filter_pp);

    Glib::RefPtr<Gtk::FileFilter> filter_any = Gtk::FileFilter::create();
    filter_any->set_name(M("FILECHOOSER_FILTER_ANY"));
    filter_any->add_pattern("*");
    dialog.add_filter(filter_any);

    int result = dialog.run();
    dialog.hide();

    if (result == Gtk::RESPONSE_OK) {
        Glib::ustring fname = dialog.get_filename();
		printf("fname=%s\n", fname.c_str());

        bool customCreated = false;

        if (!custom) {
            custom = new PartialProfile (true);
            customCreated = true;
        }

        ProcParams pp;
        ParamsEdited pe;
        int err = pp.load (fname, &pe);

        if (!err) {
            if (!customCreated && fillMode->get_active()) {
                custom->pparams->setDefaults();

                // Clearing all LocallabSpotEdited to be compliant with default pparams
                custom->pedited->locallab.spots.clear();
            }

            // For each Locallab spot, loaded profile pp only contains activated tools params
            // Missing tool params in pe shall be also set to true to avoid a "spot merge" issue
            for (int i = 0; i < (int)pe.locallab.spots.size(); i++) {
                pe.locallab.spots.at(i).set(true);
            }

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
            {
                if(!partialProfileDlg) {
                    partialProfileDlg = new PartialPasteDlg (Glib::ustring (), parent);
                }

                // opening the partial paste dialog window
                partialProfileDlg->set_title(M("PROFILEPANEL_LOADPPASTE"));
                partialProfileDlg->updateSpotWidget(&pp);
                int i = partialProfileDlg->run();
                partialProfileDlg->hide();

                if (i != Gtk::RESPONSE_OK) {
                    return;
                }

                partialProfileDlg->applyPaste (custom->pparams, !fillMode->get_active() ? custom->pedited : nullptr, &pp, &pe);
            } else {
                // custom.pparams = loadedFile.pparams filtered by ( loadedFile.pedited )
                pe.combine(*custom->pparams, pp, true);

                if (!fillMode->get_active()) {
                    *custom->pedited = pe;
                } else {
                    // Resize custom->pedited to be compliant with pe spot size
                    custom->pedited->locallab.spots.resize(pe.locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
                }
            }

            changeTo (custom, M("PROFILEPANEL_PFILE"));
        } else if (customCreated) {
            // we delete custom
            custom->deleteInstance();
            delete custom;
            custom = nullptr;
        }
    }

    return;
}

/*
 * Paste a full profile from the clipboard
 */
void ProfilePanel::paste_clicked (GdkEventButton* event)
{

    if (event->button != 1) {
        return;
    }

    if (!clipboard.hasProcParams()) {
        return;
    }

    bool prevState = changeconn.block(true);

    if (!custom) {
        custom = new PartialProfile (true); // custom pedited is initialized to false

        if (isLastSavedSelected()) {
            *custom->pparams = *lastsaved->pparams;

            // Setting LocallabSpotEdited number coherent with spots number in lastsaved->pparams
            custom->pedited->locallab.spots.clear();
            custom->pedited->locallab.spots.resize(custom->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(false));
        } else {
            const ProfileStoreEntry* entry = profiles->getSelectedEntry();

            if (entry) {
                const PartialProfile* partProfile = ProfileStore::getInstance()->getProfile (entry);
                *custom->pparams = *partProfile->pparams;

                // Setting LocallabSpotEdited number coherent with spots number in partProfile->pparams
                custom->pedited->locallab.spots.clear();
                custom->pedited->locallab.spots.resize(custom->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(false));
            }
        }

        profiles->set_active (addCustomRow());
        currRow = profiles->get_active();
    } else {
        if (fillMode->get_active()) {
            custom->pparams->setDefaults();

            // Clear all LocallabSpotEdited to be compliant with default pparams
            custom->pedited->locallab.spots.clear();
        } else if (!isCustomSelected ()) {
            if (isLastSavedSelected()) {
                *custom->pparams = *lastsaved->pparams;

                // Setting LocallabSpotEdited number coherent with spots number in lastsaved->pparams
                custom->pedited->locallab.spots.clear();
                custom->pedited->locallab.spots.resize(custom->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
            } else {
                const ProfileStoreEntry* entry = profiles->getSelectedEntry();

                if (entry) {
                    const PartialProfile* partProfile = ProfileStore::getInstance()->getProfile (entry);
                    *custom->pparams = *partProfile->pparams;

                    // Setting LocallabSpotEdited number coherent with spots number in partProfile->pparams
                    custom->pedited->locallab.spots.clear();
                    custom->pedited->locallab.spots.resize(custom->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
                }
            }
        }

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
        {
            if(!partialProfileDlg) {
                partialProfileDlg = new PartialPasteDlg (Glib::ustring (), parent);
            }

            partialProfileDlg->set_title(M("PROFILEPANEL_PASTEPPASTE"));
            partialProfileDlg->updateSpotWidget(&pp);
            int i = partialProfileDlg->run();
            partialProfileDlg->hide();

            if (i != Gtk::RESPONSE_OK) {
                return;
            }

            partialProfileDlg->applyPaste (custom->pparams, !fillMode->get_active() ? custom->pedited : nullptr, &pp, &pe);
        } else {
            // custom.pparams = clipboard.pparams filtered by ( clipboard.pedited )
            pe.combine(*custom->pparams, pp, true);

            if (!fillMode->get_active()) {
                *custom->pedited = pe;
            } else {
                // Setting LocallabSpotEdited number coherent with spots number in custom->pparams
                custom->pedited->locallab.spots.clear();
                custom->pedited->locallab.spots.resize(custom->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
            }
        }
    } else {
        if (event->state & Gdk::CONTROL_MASK)
            // custom.pparams = clipboard.pparams filtered by ( partialPaste.pedited )
        {
            if(!partialProfileDlg) {
                partialProfileDlg = new PartialPasteDlg (Glib::ustring (), parent);
            }

            partialProfileDlg->set_title(M("PROFILEPANEL_PASTEPPASTE"));
            partialProfileDlg->updateSpotWidget(&pp);
            int i = partialProfileDlg->run();
            partialProfileDlg->hide();

            if (i != Gtk::RESPONSE_OK) {
                return;
            }

            partialProfileDlg->applyPaste (custom->pparams, nullptr, &pp, nullptr);

            // Setting LocallabSpotEdited number coherent with spots number in custom->pparams
            custom->pedited->locallab.spots.clear();
            custom->pedited->locallab.spots.resize(custom->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
        } else {
            // custom.pparams = clipboard.pparams non filtered
            *custom->pparams = pp;

            // Setting LocallabSpotEdited number coherent with spots number in custom->pparams
            custom->pedited->locallab.spots.clear();
            custom->pedited->locallab.spots.resize(custom->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
        }
    }

    changeTo (custom, M("HISTORY_FROMCLIPBOARD"));
    return;
}

void ProfilePanel::changeTo (const PartialProfile* newpp, Glib::ustring profname)
{

    if (!newpp) {
        return;
    }

    if (tpc) {
        tpc->profileChange (newpp, EvProfileChanged, profname);
    }
}

void ProfilePanel::selection_changed ()
{

    if (isCustomSelected()) {
        if (!dontupdate) {
            changeTo (custom, Glib::ustring ("(" + M("PROFILEPANEL_PCUSTOM") + ")"));
        }
    } else if (isLastSavedSelected()) {
        changeTo (lastsaved, Glib::ustring ("(" + M("PROFILEPANEL_PLASTSAVED") + ")"));
    } else {
        const ProfileStoreEntry *pse = profiles->getSelectedEntry();

        if (pse->type == PSET_FOLDER) {
            // this entry is invalid, restoring the old value
            bool ccPrevState = changeconn.block(true);
            profiles->set_active(currRow);
            changeconn.block(ccPrevState);
            dontupdate = false;
            return;
        } else {
            currRow = profiles->get_active();
        }

        const PartialProfile* s = ProfileStore::getInstance()->getProfile (pse);

        if (s) {
            if (fillMode->get_active() && s->pedited) {
                ParamsEdited pe(true);

                // Setting LocallabSpotEdited number coherent with spots number in s->pparams
                pe.locallab.spots.resize(s->pparams->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));

                PartialProfile s2(s->pparams, &pe, false);
                changeTo (&s2, pse->label + "+");
            } else {
                changeTo (s, pse->label);
            }
        }
    }

    dontupdate = false;
}

void ProfilePanel::procParamsChanged(
    const rtengine::procparams::ProcParams* p,
    const rtengine::ProcEvent& ev,
    const Glib::ustring& descr,
    const ParamsEdited* paramsEdited
)
{
    // to prevent recursion, filter out the events caused by the profilepanel
    if (ev == EvProfileChanged || ev == EvPhotoLoaded) {
        return;
    }

    if (!isCustomSelected()) {
        dontupdate = true;

        if (!custom) {
            custom = new PartialProfile (true);
            custom->set(true);
            profiles->set_active (addCustomRow());
            currRow = profiles->get_active();
        } else {
            profiles->set_active(getCustomRow());
            currRow = profiles->get_active();
        }
    }

    *custom->pparams = *p;

    // Setting LocallabSpotEdited number coherent with spots number in p
    custom->pedited->locallab.spots.clear();
    custom->pedited->locallab.spots.resize(p->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
}

void ProfilePanel::clearParamChanges()
{
}

/** @brief Initialize the Profile panel with a default profile, overridden by the last saved profile if provided
 *
 * The file tree has already been created on object's construction. We add here the Custom, LastSaved and/or Internal item.
 *
 * @param profileFullPath   full path of the profile; must start by the virtual root (${G} or ${U}, and without suffix
 * @param lastSaved         pointer to the last saved ProcParam; may be NULL
 */
void ProfilePanel::initProfile (const Glib::ustring& profileFullPath, ProcParams* lastSaved)
{

    const ProfileStoreEntry *pse = nullptr;
    const PartialProfile *defprofile = nullptr;

    bool ccPrevState = changeconn.block(true);

    if (custom) {
        custom->deleteInstance();
        delete custom;
        custom = nullptr;
    }

    if (lastsaved) {
        lastsaved->deleteInstance();
        delete lastsaved;
        lastsaved = nullptr;
    }

    if (lastSaved) {
        ParamsEdited* pe = new ParamsEdited(true);
        // Setting LocallabSpotEdited number coherent with lastSaved->locallab spots number (initialized at true such as pe)
        pe->locallab.spots.resize(lastSaved->locallab.spots.size(), LocallabParamsEdited::LocallabSpotEdited(true));
        // copying the provided last saved profile to ProfilePanel::lastsaved
        lastsaved = new PartialProfile(lastSaved, pe);
    }

    // update the content of the combobox; will add 'custom' and 'lastSaved' if necessary
    updateProfileList();

    Gtk::TreeIter lasSavedEntry;

    // adding the Last Saved combobox entry, if needed
    if (lastsaved) {
        lasSavedEntry = getLastSavedRow();
    }

    if (!(pse = ProfileStore::getInstance()->findEntryFromFullPath(profileFullPath))) {
        // entry not found, pse = the Internal ProfileStoreEntry
        pse = ProfileStore::getInstance()->getInternalDefaultPSE();
    }

    defprofile = ProfileStore::getInstance()->getProfile (pse);

    // selecting the "Internal" entry
    profiles->setInternalEntry ();
    currRow = profiles->get_active();

    if (lastsaved) {
        if (lasSavedEntry) {
            profiles->set_active (lasSavedEntry);
        }

        currRow = profiles->get_active();

        if (tpc) {
            tpc->setDefaults   (lastsaved->pparams);
            tpc->profileChange (lastsaved, EvPhotoLoaded, profiles->getSelectedEntry()->label, nullptr, true);
        }
    } else {
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

void ProfilePanel::setInitialFileName (const Glib::ustring& filename)
{
    lastFilename = Glib::path_get_basename(filename) + paramFileExtension;
    imagePath = Glib::path_get_dirname(filename);
}

void ProfilePanel::profileFillModeToggled()
{
    if (fillMode->get_active()) {
        // The button is pressed, we'll use the profileFillModeOnImage
        profileFillImage->set_from_icon_name(modeOn);
    } else {
        // The button is released, we'll use the profileFillModeOffImage
        profileFillImage->set_from_icon_name(modeOff);
    }
}

void ProfilePanel::writeOptions()
{
    options.filledProfile = fillMode->get_active();
}

