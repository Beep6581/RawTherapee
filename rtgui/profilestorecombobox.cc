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
#include "rtengine/profilestore.h"
#include "profilestorecombobox.h"

#include "rtengine/dynamicprofile.h"
#include "options.h"
#include "toolpanel.h"
#include "guiutils.h"

ProfileStoreLabel::ProfileStoreLabel (const ProfileStoreEntry *entry) : Gtk::Label (entry->label), entry (entry)
{
    set_alignment (0, 0.5);
    set_ellipsize (Pango::ELLIPSIZE_END);
    show();
}

ProfileStoreComboBox::ProfileStoreComboBox ()
{
    updateProfileList();
    setPreferredWidth (50, 120);
}

Glib::ustring ProfileStoreComboBox::getCurrentLabel() const
{
    const Gtk::TreeModel::const_iterator currRow = get_active();

    if (currRow) {
        const ProfileStoreEntry *currEntry = (*currRow)[methodColumns.profileStoreEntry];
        return currEntry->label;
    }

    return {};
}

const ProfileStoreEntry* ProfileStoreComboBox::getSelectedEntry() const
{
    const Gtk::TreeModel::const_iterator currRow_ = get_active();
    const Gtk::TreeModel::Row currRow = *currRow_;

    if (currRow) {
        return currRow[methodColumns.profileStoreEntry];
    } else {
        return nullptr;
    }
}

/** @brief Recursive method to update the combobox entries */
void ProfileStoreComboBox::refreshProfileList_ (Gtk::TreeModel::Row *parentRow, int parentFolderId, bool initial, const std::vector<const ProfileStoreEntry*> *entryList)
{
    for (auto entry : *entryList) {
        if (entry->parentFolderId == parentFolderId) {  // filtering the entry of the same folder
            if (entry->type == PSET_FOLDER) {
                Glib::ustring folderPath ( ProfileStore::getInstance()->getPathFromId (entry->folderId) );

                if (options.useBundledProfiles || ((folderPath != "${G}" ) && (folderPath != "${U}" ))) {
                    // creating the new submenu
                    Gtk::TreeModel::Row newSubMenu;

                    if (initial) {
                        newSubMenu = * (refTreeModel->append());
                    } else {
                        newSubMenu = * (refTreeModel->append (parentRow->children()));
                    }

                    // creating and assigning the custom Label object
                    newSubMenu[methodColumns.label] = entry->label;
                    newSubMenu[methodColumns.profileStoreEntry] = entry;
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 18
                    // HACK: Workaround for bug in Gtk+ 3.18...
                    Gtk::TreeModel::Row menuHeader = * (refTreeModel->append (newSubMenu->children()));
                    menuHeader[methodColumns.label] = "-";
                    menuHeader[methodColumns.profileStoreEntry] = entry;
#endif
                    refreshProfileList_ (&newSubMenu, entry->folderId, false, entryList);
                } else {
                    refreshProfileList_ (parentRow, entry->folderId, true, entryList);
                }
            } else {
                Gtk::TreeModel::Row newItem;

                // creating a menu entry
                if (initial) {
                    newItem = * (refTreeModel->append());
                } else {
                    newItem = * (refTreeModel->append (parentRow->children()));
                }

                newItem[methodColumns.label] = entry->label;
                newItem[methodColumns.profileStoreEntry] = entry;
            }
        }
    }
}
/** @brief Get the ProfileStore's entry list and recreate the combobox entries.
  * If you want to update the ProfileStore list itself (rescan the dir tree), use the "ProfileStore::parseProfiles" method instead
  *
  * This method has to be called by the ProfileStoreListener having a ProfileStoreComboBox.
  */
void ProfileStoreComboBox::updateProfileList ()
{
    // clear items
    clear();
    refTreeModel.clear();
    // Create the Tree model
    refTreeModel = Gtk::TreeStore::create (methodColumns);
    // Assign the model to the Combobox
    set_model (refTreeModel);

    // this will lock the profilestore's entry list too
    const std::vector<const ProfileStoreEntry*> *entryList = ProfileStore::getInstance()->getFileList();

    //profileStore.dumpFolderList();
    refreshProfileList_ (NULL, entryList->at (0)->parentFolderId, true, entryList);

    if (entryList->at (0)->parentFolderId != 0) {
        // special case for the Internal default entry
        addRow (ProfileStore::getInstance()->getInternalDefaultPSE());
    }

    // releasing the profilestore's entry list mutex
    ProfileStore::getInstance()->releaseFileList();

    pack_start (methodColumns.label, false);

    Gtk::CellRendererText* cellRenderer = dynamic_cast<Gtk::CellRendererText*> (get_first_cell());
    cellRenderer->property_ellipsize() = Pango::ELLIPSIZE_MIDDLE;
    cellRenderer->property_ellipsize_set() = true;
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromEntry_ (Gtk::TreeModel::Children childs, const ProfileStoreEntry *pse) const
{

    for (const auto& iter : childs) {
        const Gtk::TreeModel::Row row = *iter;
        // Hombre: is there a smarter way of knowing if this row has childs?
        const ProfileStoreEntry *pse_ = row[methodColumns.profileStoreEntry];

        if (pse_->type == PSET_FOLDER) {
            const Gtk::TreeIter rowInSubLevel = findRowFromEntry_ (iter->children(), pse);

            if (rowInSubLevel) {
                // entry found
                return rowInSubLevel;
            }
        } else if (pse_ == pse) {
            // entry found
            return iter;
        }
    }

    return childs.end();
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromEntry (const ProfileStoreEntry *pse) const
{
    Gtk::TreeModel::Children childs = refTreeModel->children();

    if (pse) {
        return findRowFromEntry_ (childs, pse);
    }

    return childs.end();
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromFullPath_ (Gtk::TreeModel::Children childs, int parentFolderId, const Glib::ustring &name) const
{

    for (const auto& iter : childs) {
        const Gtk::TreeModel::Row row = *iter;
        // Hombre: is there a smarter way of knowing if this row has childs?
        const ProfileStoreEntry *pse = row[methodColumns.profileStoreEntry];

        if (pse->type == PSET_FOLDER) {
            const Gtk::TreeIter rowInSubLevel = findRowFromFullPath_ (iter->children(), parentFolderId, name);

            if (rowInSubLevel) {
                // entry found
                return rowInSubLevel;
            }
        } else if (parentFolderId == pse->parentFolderId && name == pse->label) {
            // entry found
            return iter;
        }
    }

    return childs.end();
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromFullPath (const Glib::ustring &path) const
{

    if (path.empty()) {
        return {};
    }

    const ProfileStore *profileStore = ProfileStore::getInstance();

    if (path == DEFPROFILE_INTERNAL) {
        return findRowFromEntry (profileStore->getInternalDefaultPSE());
    }

    if (path == DEFPROFILE_DYNAMIC) {
        return findRowFromEntry (profileStore->getInternalDynamicPSE());
    }

    // removing the filename
    const Glib::ustring fName = Glib::path_get_basename(path);

    if (!fName.empty()) {
        int parentFolderId = profileStore->findFolderId (Glib::path_get_dirname (path.substr (0, path.length() - fName.length())));
        // 1. find the path in the folder list
        if (parentFolderId != -1) {
            return findRowFromFullPath_ (refTreeModel->children(), parentFolderId, fName);
        }
    }

    return {};
}

/** @brief Get the absolute full path of the active row entry.
  * @return The absolute full path of the active row entry, or the "Internal" keyword,
  *         or an empty string if the ComboBox is in an invalid state
  */
Glib::ustring ProfileStoreComboBox::getFullPathFromActiveRow() const
{
    const Gtk::TreeModel::const_iterator currRowI = get_active();

    if (!currRowI) {
        return {};
    }

    Gtk::TreeModel::Row currRow = *currRowI;

    if (currRow) {
        const ProfileStoreEntry *currEntry = currRow[methodColumns.profileStoreEntry];

        if (!currEntry) {
            return {};
        }

        const ProfileStore *profileStore = ProfileStore::getInstance();
        if (currEntry == profileStore->getInternalDefaultPSE()) {
            return DEFPROFILE_INTERNAL;
        }

        if (currEntry == profileStore->getInternalDynamicPSE()) {
            return DEFPROFILE_DYNAMIC;
        }

        return Glib::build_filename (profileStore->getPathFromId (currEntry->parentFolderId), currEntry->label);
    }

    return {};
}

bool ProfileStoreComboBox::setActiveRowFromFullPath (const Glib::ustring &path)
{
    if (!path.empty()) {
        const Gtk::TreeIter row = findRowFromFullPath (path);

        if (row) {
            set_active (row);
            return true;
        }
    }

    return false;
}

bool ProfileStoreComboBox::setActiveRowFromEntry (const ProfileStoreEntry *pse)
{
    if (pse) {
        const Gtk::TreeIter row = findRowFromEntry (pse);

        if (row) {
            set_active (row);
            return true;
        }
    }

    return false;
}

bool ProfileStoreComboBox::setInternalEntry ()
{
    return setActiveRowFromEntry (ProfileStore::getInstance()->getInternalDefaultPSE());
}

/** @brief Get the row from the first level of the tree that match the provided name */
Gtk::TreeIter ProfileStoreComboBox::getRowFromLabel (const Glib::ustring &name) const
{
    const Gtk::TreeModel::Children childs = refTreeModel->children();

    if (!name.empty()) {
        for (const auto& iter : childs) {
            const Gtk::TreeModel::Row currRow = *iter;
            const ProfileStoreEntry *pse = currRow[methodColumns.profileStoreEntry];

            if (pse->label == name) {
                return currRow;
            }
        }
    }

    return childs.end();
    //return refTreeModel->get_iter(""); // is this fast? We want to send back a null, anvalid or end() iterator object here
}

/** @brief Add a new row to the first level of the tree */
Gtk::TreeIter ProfileStoreComboBox::addRow (const ProfileStoreEntry *profileStoreEntry)
{
    Gtk::TreeIter newEntry = refTreeModel->append();
    Gtk::TreeModel::Row row = *newEntry;
    row[methodColumns.label] = profileStoreEntry->label;
    row[methodColumns.profileStoreEntry] = profileStoreEntry;
    return newEntry;
}

/** @brief Delete a row from the first level of the tree */
void ProfileStoreComboBox::deleteRow (const ProfileStoreEntry *profileStoreEntry)
{
    Gtk::TreeIter entry = findRowFromEntry(profileStoreEntry);
    if (entry) {
        refTreeModel->erase(entry);
    }
}
