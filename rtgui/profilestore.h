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
#ifndef _PROFILESTORE_
#define _PROFILESTORE_

#include <map>
#include <vector>
#include <glibmm.h>

#include "../rtengine/rtengine.h"
#include "../rtengine/noncopyable.h"

#include "threadutils.h"
#include "paramsedited.h"
#include "guiutils.h"


/** @brief This will implement callback functions for the ProfileStore
  *
  */
class ProfileStoreListener
{

public:
    virtual ~ProfileStoreListener() {}

    /** @brief Called whenever the current value has to be stored before update. */
    virtual void storeCurrentValue() {}
    /** @brief Called whenever the file list has been updated and the content of the listener has to be updated. */
    virtual void updateProfileList() = 0;
    /** @brief Called whenever the profile list has changed and the old value have to be restored (if possible). */
    virtual void restoreValue() {}
};

/// @brief ProfileStoreEntry type (folder or file)
typedef enum PSE_TYPE {
    PSET_FOLDER,
    PSET_FILE
} PSEType;


/**
 * @brief Entry of the profiles store, consisting of a type (folder/file) & label.
 *
 * Will be used as key element in the name / PartialProfile mapping
 */
class ProfileStoreEntry
{

public:

    Glib::ustring label;                  /// Label to be used in menu or combobox = profile's filename
    PSEType type;                         /// either PSET_FOLDER or PSET_FILE
    unsigned short parentFolderId;        /// index of the element's path in the folder list; id == 0 is reserved
    unsigned short folderId;              /// index of the folder's own path in the folder list; will be null for file entries

    /** @brief Create a new ProfileStoreLabel with null values, that will have to be set later with setValues or the copy operator
      */
    ProfileStoreEntry();

    /** @brief Create a new ProfileStoreLabel with values
      * @param label         Label to be used in menu or combobox; also used as profile's filename
      * @param type          either PSET_FOLDER or PSET_FILE
      * @param parentFolder  index of the elements's path in the folder list
      * @param folder        index of the folder's own path in the folder list
      */
    ProfileStoreEntry(Glib::ustring label, PSEType type, unsigned short parentFolder, unsigned short folder);

    /** @brief Set the values of the object after its instantiation
      * @param label         Label to be used in menu or combobox; also used as profile's filename
      * @param type          either PSET_FOLDER or PSET_FILE
      * @param parentFolder  index of the elements's path in the folder list
      * @param folder        index of the folder's own path in the folder list
      */
    void setValues(Glib::ustring label, PSEType type, unsigned short parentFolder, unsigned short folder);
};


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


/** @brief Store the profiles bundled with RawTharapee and created by the user.
  *
  * This store can be queried by the GUI to display a Tree of the profiles available
  * in the user's and system's profile directory and subdirectories.
  */
class ProfileStore :
    public rtengine::NonCopyable
{

    typedef enum {
        STORESTATE_NOTINITIALIZED,
        STORESTATE_BEINGINITIALIZED,
        STORESTATE_INITIALIZED,
        STORESTATE_DELETED
    } StoreState;

private:
    struct SortProfiles {
        bool operator ()(const ProfileStoreEntry* const a1, const ProfileStoreEntry* const a2)
        {
            return a1->parentFolderId == a2->parentFolderId ? a1->label < a2->label :  a1->parentFolderId < a2->parentFolderId;
        }
    };

    MyMutex *parseMutex;
    StoreState storeState;
    rtengine::procparams::AutoPartialProfile *internalDefaultProfile;
    ProfileStoreEntry *internalDefaultEntry;

    /** Alphabetically ordered list of folder and files through Gtk::Label sub-class;
      * ready to be used in Menu and Combobox
      * The first element (#0) will be a fake path so ProfileStoreEntry can be attached to the root container */
    std::vector<Glib::ustring> folders;

    /** Alphabetically ordered list of folder and files through Gtk::Label derived class;
      * ready to be used in Menu and Combobox */
    std::vector<const ProfileStoreEntry*> entries;

    /** List of PartialProfiles from the indexed files */
    std::map<const ProfileStoreEntry*, rtengine::procparams::AutoPartialProfile*> partProfiles;

    /** List of the client of this store */
    std::list<ProfileStoreListener*> listeners;

    /** @brief Method to recursively parse a profile folder with a level depth arbitrarily limited to 3
      *
      * @param realPath       current full path of the scanned directory ; e.g.:  ~/MyProfiles/
      * @param virtualPath    current full path that will be saved in "options" ; must start with either ${U} or ${G},
      *                       standing for User's and Global's (RT) profile folder, respectively
      * @param currDir        name of the directory to scan; it's the last element of the virtualPath
      * @param parentId       path entry of the parent folder
      * @param level          current level of the directory tree
      * @param displayLevel0  if true, level 0 is created in order to have a User's and Bundled profiles separation (i.e. 2 root directories are expected)
      *                       if false, only one root directory is expected
      */
    bool parseDir (Glib::ustring& realPath, Glib::ustring& virtualPath, Glib::ustring& currDir, unsigned int parentId, unsigned char level, bool displayLevel0);
    void _parseProfiles ();
    void clearFileList ();
    void clearProfileList ();
    const ProfileStoreEntry* findEntryFromFullPathU(Glib::ustring path);

public:

    ProfileStore();
    ~ProfileStore();

    bool init ();
    void parseProfiles ();
    int findFolderId(const Glib::ustring &path);
    const ProfileStoreEntry*                     findEntryFromFullPath(Glib::ustring path);
    const rtengine::procparams::PartialProfile*  getProfile (Glib::ustring path);
    const rtengine::procparams::PartialProfile*  getProfile (const ProfileStoreEntry* entry);
    const std::vector<const ProfileStoreEntry*>* getFileList ();
    void                                         releaseFileList ();
    const rtengine::procparams::ProcParams*      getDefaultProcParams (bool isRaw);
    const rtengine::procparams::PartialProfile*  getDefaultPartialProfile (bool isRaw);
    const Glib::ustring                          getPathFromId(int folderId);
    const ProfileStoreEntry*                     getInternalDefaultPSE()
    {
        return internalDefaultEntry;
    }

    void addListener(ProfileStoreListener *listener);
    void removeListener(ProfileStoreListener *listener);

    void dumpFolderList();

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

extern ProfileStore profileStore;

#endif
