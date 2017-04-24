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
#ifndef _PROFILESTORECOMBOBOX_
#define _PROFILESTORECOMBOBOX_

#include <map>
#include <vector>
#include <glibmm.h>

#include "../rtengine/rtengine.h"
#include "../rtengine/noncopyable.h"
#include "../rtengine/profilestore.h"

#include "threadutils.h"
#include "paramsedited.h"
#include "guiutils.h"


/**
 * @brief subclass of Gtk::Label with extra fields for Combobox and Menu, to link with a ProfileStoreEntry
 */
class ProfileStoreLabel : public Gtk::Label
{

public:
    const ProfileStoreEntry *entry;

#ifndef NDEBUG
    ProfileStoreLabel() : Gtk::Label("*** error ***"), entry(nullptr) {}
#else
    ProfileStoreLabel() : Gtk::Label(""), entry(NULL) {}
#endif

    /** @brief Create a new ProfileStoreLabel
      *
      * @param entry      Pointer to the ProfileStoreEntry object, be it a directory or a file
      */
    explicit ProfileStoreLabel(const ProfileStoreEntry *entry);
    ProfileStoreLabel (const ProfileStoreLabel &other);
};

class ProfileStoreComboBox : public MyComboBox
{

protected:
    class MethodColumns : public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::ustring> label;
        Gtk::TreeModelColumn<const ProfileStoreEntry*> profileStoreEntry;
        MethodColumns()
        {
            add(label);
            add(profileStoreEntry);
        }
    };

    Glib::RefPtr<Gtk::TreeStore> refTreeModel;
    MethodColumns methodColumns;
    void refreshProfileList_ (Gtk::TreeModel::Row *parentRow, int parentFolderId, bool initial, const std::vector<const ProfileStoreEntry*> *entryList);
    Gtk::TreeIter findRowFromEntry_ (Gtk::TreeModel::Children childs, const ProfileStoreEntry *pse);
    Gtk::TreeIter findRowFromFullPath_(Gtk::TreeModel::Children childs, int parentFolderId, Glib::ustring &name);

public:
    ProfileStoreComboBox();
    void updateProfileList();
    Glib::ustring getCurrentLabel();
    const ProfileStoreEntry* getSelectedEntry();
    Gtk::TreeIter findRowFromEntry (const ProfileStoreEntry *pse);
    Gtk::TreeIter findRowFromFullPath (Glib::ustring path);
    Glib::ustring getFullPathFromActiveRow ();
    bool setActiveRowFromFullPath (Glib::ustring oldPath);
    bool setActiveRowFromEntry (const ProfileStoreEntry *pse);
    bool setInternalEntry ();
    Gtk::TreeIter getRowFromLabel(Glib::ustring name);
    Gtk::TreeIter addRow(const ProfileStoreEntry *profileStoreEntry);
};

#endif
