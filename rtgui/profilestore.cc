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
#include "profilestore.h"
#include "options.h"
#include "toolpanel.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"

ProfileStore profileStore;

using namespace rtengine;
using namespace rtengine::procparams;

ProfileStore::ProfileStore () : parseMutex(NULL), storeState(STORESTATE_NOTINITIALIZED), internalDefaultProfile(NULL), internalDefaultEntry(NULL) {
    internalDefaultProfile = new AutoPartialProfile();
    internalDefaultProfile->set(true);
}

bool ProfileStore::init () {
    if (storeState == STORESTATE_DELETED)
        return false;
    if (storeState == STORESTATE_NOTINITIALIZED) {
        storeState = STORESTATE_BEINGINITIALIZED;
        parseMutex = new MyMutex();
        _parseProfiles ();
        storeState = STORESTATE_INITIALIZED;
    }
    return true;
}

ProfileStore::~ProfileStore () {

    // This lock prevent object's suppression while scanning the directories
    storeState = STORESTATE_DELETED;
    MyMutex::MyLock lock(*parseMutex);

    clearProfileList ();
    partProfiles.clear ();
    clearFileList();
    delete internalDefaultProfile;
    lock.release();
    delete parseMutex;
    parseMutex = NULL;
}

/*
 * Public method to parse the profiles' directories
 * Since there's a race condition in the multithreaded environment on this object,
 * parseProfiles may need to ask for initialization of this object, and then will
 * ask a mutex lock on it, has it been initialized by this call or not
 *
 * This method will scan the directory tree again and update the profile list. When finished,
 * the listeners will be called in order to update with the new list
 */
void ProfileStore::parseProfiles () {

    if (!init())
        // I don't even know if this situation can occur
        return;

    for (std::list<ProfileStoreListener*>::iterator i=listeners.begin(); i!=listeners.end(); ++i)
        (*i)->storeCurrentValue();

    {
    MyMutex::MyLock lock(*parseMutex);

    _parseProfiles ();
    }

    for (std::list<ProfileStoreListener*>::iterator i=listeners.begin(); i!=listeners.end(); ++i) {
        (*i)->updateProfileList();
        (*i)->restoreValue();
    }
}

void ProfileStore::_parseProfiles () {

    // Acquire the GUI, since the tree model can interact with combobox
    GThreadLock threadLock;

    // clear loaded profiles
    folders.clear();
    clearFileList();
    clearProfileList ();

    folders.push_back("<<< ROOT >>>");  // Fake path, so parentFolderId == 0 will be used to attach a ProfileStoreEntry to the root container, not sub-menu

    Glib::ustring p1 = options.getUserProfilePath();
    Glib::ustring p2 = options.getGlobalProfilePath();
    bool displayLevel0 = options.useBundledProfiles && !p1.empty() && !p2.empty() && p1 != p2;

    Glib::ustring virtualPath("${U}");
    Glib::ustring currDir("${U}");
    parseDir (p1, virtualPath, currDir, 0, 0, displayLevel0);
    if (displayLevel0) {
        virtualPath = "${G}";
        currDir = "${G}";
        parseDir (p2, virtualPath, currDir, 0, 0, displayLevel0);
    }
    // entries and partProfiles are empty, but the entry and profiles already exist (they have survived to clearFileList and clearProfileList)
    if (!internalDefaultEntry)
        internalDefaultEntry = new ProfileStoreEntry(Glib::ustring("(") + M("PROFILEPANEL_PINTERNAL") + Glib::ustring(")"), PSET_FILE, 0, 0);
    entries.push_back(internalDefaultEntry);
    partProfiles[internalDefaultEntry] = internalDefaultProfile;


    // Check if the default profiles has been found.
    if (findEntryFromFullPathU(options.defProfRaw) == NULL) {
        options.setDefProfRawMissing(true);
        if (options.rtSettings.verbose)
            printf("WARNING: Default profile \"%s\" for raw images not found!\n", options.defProfRaw.c_str());
    }
    if (findEntryFromFullPathU(options.defProfImg) == NULL) {
        options.setDefProfImgMissing(true);
        if (options.rtSettings.verbose)
            printf("WARNING: Default profile \"%s\" for standard images not found!\n", options.defProfImg.c_str());
    }
}

/// @return Returns true if some files has been found (directories are ignored)
bool ProfileStore::parseDir (Glib::ustring& realPath, Glib::ustring& virtualPath, Glib::ustring& currDir, unsigned int parentId, unsigned char level, bool displayLevel0) {

    bool fileFound = false;
    unsigned int folder = 0; // folder's own Id

    // reload the available profiles from the profile dir
    if (!realPath.empty() && safe_file_test(realPath, Glib::FILE_TEST_EXISTS) && safe_file_test (realPath, Glib::FILE_TEST_IS_DIR)) {

        // add this entry to the folder list
        folders.push_back(virtualPath);
        folder = (unsigned int)(folders.size())-1;

        if (level>0 || displayLevel0) {
            // replace the virtual folder name by a localized text
            if (currDir == "${U}")
                currDir = M("PROFILEPANEL_MYPROFILES");
            else if (currDir == "${G}")
                currDir = M("PROFILEPANEL_GLOBALPROFILES");

            // add this localized text to the file list
            entries.push_back( new ProfileStoreEntry(currDir, PSET_FOLDER, parentId, folder) );
        }

        // walking through the directory
        Glib::Dir* dir = NULL;
        dir = new Glib::Dir (realPath);
        for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
            currDir = *i;
            if (currDir == "." || currDir == "..")
                continue;

            Glib::ustring fname = Glib::build_filename(realPath, currDir);
            if (safe_file_test (fname, Glib::FILE_TEST_IS_DIR)) {
                Glib::ustring vp(Glib::build_filename(virtualPath, currDir));
                Glib::ustring rp(Glib::build_filename(realPath,    currDir));
                parseDir (rp, vp, currDir, folder, level+1, 0);
            }
            else {
                size_t lastdot = currDir.find_last_of ('.');
                if (lastdot!=Glib::ustring::npos && lastdot<=currDir.size()-4 && !currDir.casefold().compare (lastdot, 4, paramFileExtension)) {
                    // file found
                    if( options.rtSettings.verbose )
                        printf ("Processing file %s...", fname.c_str());

                    Glib::ustring name = currDir.substr(0,lastdot);

                    // create the partial profile
                    AutoPartialProfile *pProf = new AutoPartialProfile();
                    int res = pProf->load (fname);
                    if (!res && pProf->pparams->ppVersion>=220) {
                        fileFound = true;

                        if( options.rtSettings.verbose )
                            printf ("OK\n");

                        // adding this file to the list
                        ProfileStoreEntry* filePSE = new ProfileStoreEntry(name, PSET_FILE, folder, 0);
                        entries.push_back(filePSE);

                        // map the partial profile
                        partProfiles[filePSE] = pProf;
                        //partProfiles.insert( std::pair<ProfileStoreEntry*, rtengine::procparams::AutoPartialProfile*> (folderPSE, pProf) );
                    }
                    else if( options.rtSettings.verbose ) {
                        printf ("failed!\n");
                    }
                }
            }
        }
        delete dir;
    }

    if (!fileFound && (level>0 || displayLevel0)) {
        // no files found in this level, we delete the subdirectory entry
        folders.pop_back();
        entries.pop_back();
    }

    return fileFound;
}

int ProfileStore::findFolderId(const Glib::ustring &path) {
    for (std::vector<Glib::ustring>::iterator i=folders.begin(); i!=folders.end(); i++) {
        if (*i == path) {
            return i - folders.begin();
        }
    }
    return -1;
}

/** @brief Return the ProfileStoreEntry object that match the given file and path
  *
  * @param fullPath  Path of the file; the filename may end by the standard extension,
  *                  but have to begin with a virtual location ( ${G} or ${U} )
  *                  Will return null on invalid path or if the entry can't be found
  */
const ProfileStoreEntry* ProfileStore::findEntryFromFullPathU(Glib::ustring path) {

    if (path.empty())
        return NULL;

    if (path == DEFPROFILE_INTERNAL)
        return internalDefaultEntry;

    size_t lastdot = path.find_last_of ('.');
    if (lastdot!=Glib::ustring::npos && lastdot<=path.size()-4 && !path.casefold().compare (lastdot, 4, paramFileExtension))
        // removing the extension
        path = path.substr(0,lastdot);

    // removing the filename
    Glib::ustring fName = Glib::path_get_basename(path);
    if (!fName.empty()) {
        path = path.substr(0, path.length()-fName.length());
    }
    else {
        // path is malformed, returning NULL;
        return NULL;
    }

    path = Glib::path_get_dirname(path);

    // 1. find the path in the folder list
    int parentFolderId = findFolderId(path);

    if (parentFolderId == -1) {
        return NULL;
    }

    // 2. find the entry that match the given filename and parentFolderId
    for (std::vector<const ProfileStoreEntry*>::iterator i=entries.begin(); i!=entries.end(); i++) {
        if (((*i)->parentFolderId)==parentFolderId  &&  (*i)->label==fName)
            return *i;
    }
    return NULL;
}

/** Protected version of findEntryFromFullPathU */
const ProfileStoreEntry* ProfileStore::findEntryFromFullPath(Glib::ustring path) {
    MyMutex::MyLock lock(*parseMutex);
    return findEntryFromFullPathU(path);
}

const PartialProfile* ProfileStore::getProfile (Glib::ustring path) {

    if (!init())
        // I don't even know if this situation can occur
        return NULL;

    const ProfileStoreEntry *pse = findEntryFromFullPath(path);
    if (!pse)
        return NULL;

    return getProfile(pse);
}

const PartialProfile* ProfileStore::getProfile (const ProfileStoreEntry* entry) {

    if (!init())
        // I don't even know if this situation can occur
        return NULL;

    MyMutex::MyLock lock(*parseMutex);

    if (entry == internalDefaultEntry)
        return internalDefaultProfile;

    std::map<const ProfileStoreEntry*, rtengine::procparams::AutoPartialProfile*>::iterator iter = partProfiles.find(entry);
    if (iter != partProfiles.end()) {
        return iter->second;
    }
    else {
        // This shouldn't happen!
        #ifndef NDEBUG
        printf("WARNING! Profile not found!\n");
        #endif
        return NULL;
    }
}

/** @brief Get a pointer to the profile's vector list
 *
 * This method grants you unique access to the vector list through Mutex locking.
 * When you're done with the file list, you MUST call the releaseFileList method to release the lock.
 */
const std::vector<const ProfileStoreEntry*>* ProfileStore::getFileList () {
    /*if (!init()) {
        // I don't even know if this situation can occur
        return NULL;
    }*/

    parseMutex->lock();

    return &entries;
}

void ProfileStore::releaseFileList() {
    parseMutex->unlock();
}

/*
 * Send back a pointer to the default procparams for raw or standard images.
 * If the profile doesn't already exist in the profile list,
 * it will add it with default internal values, so this method never fails
 */
const ProcParams* ProfileStore::getDefaultProcParams (bool isRaw) {

    if (!init())
        // I don't even know if this situation can occur
        return NULL;
    //Note: the mutex is locked in getProfile, called below

    const PartialProfile* pProf = getProfile (isRaw ? options.defProfRaw : options.defProfImg);

    if (!pProf) pProf = internalDefaultProfile;

    return pProf->pparams;
}

/*
 * Send back a pointer to the default partial profile for raw or standard images.
 * If it doesn't already exist in the profile list, it will add it with default internal values,
 * so this method will never fails
 */
const PartialProfile* ProfileStore::getDefaultPartialProfile (bool isRaw) {

    if (!init())
        // I don't even know if this situation can occur
        return NULL;
    //Note: the mutex is locked in getProfile, called below

    const PartialProfile* pProf = getProfile (isRaw ? options.defProfRaw : options.defProfImg);

    if (!pProf) pProf = internalDefaultProfile;

    return pProf;
}

const Glib::ustring ProfileStore::getPathFromId(int folderId) {
	return folders.at(folderId);
}


void ProfileStore::clearFileList() {
    for (std::vector<const ProfileStoreEntry*>::iterator i=entries.begin(); i!=entries.end(); ++i)
        if (*i != internalDefaultEntry) delete *i;
    entries.clear();
}

void ProfileStore::clearProfileList() {
    for (std::map<const ProfileStoreEntry*, rtengine::procparams::AutoPartialProfile*>::iterator i=partProfiles.begin(); i!=partProfiles.end(); ++i)
        if (i->second != internalDefaultProfile) delete i->second;
    partProfiles.clear();
}

void ProfileStore::addListener(ProfileStoreListener *listener) {
    listeners.push_back(listener);
}

void ProfileStore::removeListener(ProfileStoreListener *listener) {
    listeners.remove(listener);
}

ProfileStoreEntry::ProfileStoreEntry() : label(""), type(PSET_FOLDER), parentFolderId(0), folderId(0) {}

ProfileStoreEntry::ProfileStoreEntry(Glib::ustring label, PSEType type, unsigned short parentFolder, unsigned short folder) : label(label), type(type), parentFolderId(parentFolder), folderId(folder) {}

void ProfileStoreEntry::setValues(Glib::ustring label, PSEType type, unsigned short parentFolder, unsigned short folder) {
    this->label = label;
    this->type = type;
    parentFolderId = parentFolder;
    folderId = folder;
}

ProfileStoreLabel::ProfileStoreLabel(const ProfileStoreEntry *entry) : Gtk::Label(entry->label), entry(entry) {
    set_alignment(0, 0.5);
    show();
}

ProfileStoreComboBox::ProfileStoreComboBox () {
    updateProfileList();
}

Glib::ustring ProfileStoreComboBox::getCurrentLabel() {
    Glib::ustring currLabel;
    Gtk::TreeModel::iterator currRow = get_active();

    if (currRow) {
        const ProfileStoreEntry *currEntry = (*currRow)[methodColumns.profileStoreEntry];
        return currEntry->label;
    }
    return currLabel;
}

const ProfileStoreEntry* ProfileStoreComboBox::getSelectedEntry() {
    Gtk::TreeModel::iterator currRow_ = get_active();
    Gtk::TreeModel::Row currRow = *currRow_;
    if (currRow)
        return currRow[methodColumns.profileStoreEntry];
    else
        return NULL;
}

/** @brief Recursive method to update the combobox entries */
void ProfileStoreComboBox::refreshProfileList_ (Gtk::TreeModel::Row *parentRow, int parentFolderId, bool initial, const std::vector<const ProfileStoreEntry*> *entryList) {
    for (std::vector<const ProfileStoreEntry*>::const_iterator i=entryList->begin(); i!=entryList->end(); i++) {
        if ((*i)->parentFolderId == parentFolderId) {  // filtering the entry of the same folder
            if ((*i)->type == PSET_FOLDER) {
                Glib::ustring folderPath( profileStore.getPathFromId((*i)->folderId) );
                if (options.useBundledProfiles || ((folderPath != "${G}" ) && (folderPath != "${U}" ))) {
                    // creating the new submenu
                    Gtk::TreeModel::Row newSubMenu;
                    if (initial) {
                        newSubMenu = *(refTreeModel->append());
                    }
                    else {
                        newSubMenu = *(refTreeModel->append(parentRow->children()));
                    }

                    // creating and assigning the custom Label object
                    newSubMenu[methodColumns.label] = (*i)->label;
                    newSubMenu[methodColumns.profileStoreEntry] = *i;

                    refreshProfileList_ (&newSubMenu, (*i)->folderId, false, entryList);
                }
                else
                    refreshProfileList_ (parentRow, (*i)->folderId, true, entryList);
            }
            else {
                Gtk::TreeModel::Row newItem;
                // creating a menu entry
                if (initial)
                    newItem = *(refTreeModel->append());
                else
                    newItem = *(refTreeModel->append(parentRow->children()));
                newItem[methodColumns.label] = (*i)->label;
                newItem[methodColumns.profileStoreEntry] = *i;
            }
        }
    }
}

/** @brief Get the ProfileStore's entry list and recreate the combobox entries.
  * If you want to update the ProfileStore list itself (rescan the dir tree), use the "ProfileStore::parseProfiles" method instead
  *
  * This method has to be called by the ProfileStoreListener having a ProfileStoreComboBox.
  */
void ProfileStoreComboBox::updateProfileList () {

    // clear items
    clear();
    refTreeModel.clear();
    // Create the Tree model
    refTreeModel = Gtk::TreeStore::create(methodColumns);
    // Assign the model to the Combobox
    set_model(refTreeModel);


    // this will lock the profilestore's entry list too
    const std::vector<const ProfileStoreEntry*> *entryList = profileStore.getFileList();

    Gtk::TreeModel::Row root;
    refreshProfileList_ (&root, entryList->at(0)->parentFolderId, true, entryList);

    if (entryList->at(0)->parentFolderId != 0) {
        // special case for the Internal default entry
        addRow(profileStore.getInternalDefaultPSE());
    }

    // releasing the profilestore's entry list mutex
    profileStore.releaseFileList();

    pack_start(methodColumns.label, false);
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromEntry_ (Gtk::TreeModel::Children childs, const ProfileStoreEntry *pse) {
    Gtk::TreeModel::Row row;
    Gtk::TreeIter rowInSubLevel;
    for(Gtk::TreeModel::Children::iterator iter = childs.begin(); iter != childs.end(); ++iter) {
        row = *iter;
        // Hombre: is there a smarter way of knowing if this row has childs?
        const ProfileStoreEntry *pse_ = row[methodColumns.profileStoreEntry];
        if (pse_->type == PSET_FOLDER) {
            rowInSubLevel = findRowFromEntry_ (iter->children(), pse);
            if (rowInSubLevel) {
                // entry found
                return rowInSubLevel;
            }
        }
        else if (pse_ == pse) {
            // entry found
            return iter;
        }
    }
    return childs.end();
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromEntry (const ProfileStoreEntry *pse) {
    Gtk::TreeModel::Children childs = refTreeModel->children();
    if (pse) {
        Gtk::TreeIter row = findRowFromEntry_ (childs, pse);
        return row;
    }
    return childs.end();
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromFullPath_ (Gtk::TreeModel::Children childs, int parentFolderId, Glib::ustring &name) {
    Gtk::TreeModel::Row row;
    Gtk::TreeIter rowInSubLevel;
    for(Gtk::TreeModel::Children::iterator iter = childs.begin(); iter != childs.end(); ++iter) {
        row = *iter;
        // Hombre: is there a smarter way of knowing if this row has childs?
        const ProfileStoreEntry *pse = row[methodColumns.profileStoreEntry];
        if (pse->type == PSET_FOLDER) {
            rowInSubLevel = findRowFromFullPath_ (iter->children(), parentFolderId, name);
            if (rowInSubLevel) {
                // entry found
                return rowInSubLevel;
            }
        }
        else if (parentFolderId==pse->parentFolderId && name==pse->label) {
            // entry found
            return iter;
        }
    }
    return childs.end();
}

Gtk::TreeIter ProfileStoreComboBox::findRowFromFullPath (Glib::ustring path) {
    Gtk::TreeIter row;

    if (path.empty())
        return row;

    if (path == DEFPROFILE_INTERNAL) {
        row = findRowFromEntry(profileStore.getInternalDefaultPSE());
        return row;
    }

    // removing the filename
    Glib::ustring fName = Glib::path_get_basename(path);
    if (!fName.empty()) {
        path = path.substr(0, path.length()-fName.length());
    }
    else {
        // path is malformed;
        return row;
    }

    path = Glib::path_get_dirname(path);
    int parentFolderId = profileStore.findFolderId(path);

    // 1. find the path in the folder list
    if (parentFolderId != -1)
        row = findRowFromFullPath_ (refTreeModel->children(), parentFolderId, fName);

    return row;
}

/** @brief Get the absolute full path of the active row entry.
  * @return The absolute full path of the active row entry, or the "Internal" keyword,
  *         or an empty string if the ComboBox is in an invalid state
  */
Glib::ustring ProfileStoreComboBox::getFullPathFromActiveRow() {
    Glib::ustring path;
    Gtk::TreeModel::iterator currRowI = get_active();
    if (!currRowI)
        return path;

    Gtk::TreeModel::Row currRow = *currRowI;

    if (currRow) {

        const ProfileStoreEntry *currEntry = currRow[methodColumns.profileStoreEntry];
        if (!currEntry)
            return path;

        if (currEntry == profileStore.getInternalDefaultPSE())
            return Glib::ustring(DEFPROFILE_INTERNAL);

        path = Glib::build_filename(profileStore.getPathFromId(currEntry->parentFolderId), currEntry->label);
    }
    return path;
}

bool ProfileStoreComboBox::setActiveRowFromFullPath(Glib::ustring path) {
    if (!path.empty()) {
        Gtk::TreeIter row = findRowFromFullPath(path);
        if (row) {
            set_active(row);
            return true;
        }
    }
    return false;
}

bool ProfileStoreComboBox::setActiveRowFromEntry(const ProfileStoreEntry *pse) {
    if (pse) {
        Gtk::TreeIter row = findRowFromEntry(pse);
        if (row) {
            set_active(row);
            return true;
        }
    }
    return false;
}

bool ProfileStoreComboBox::setInternalEntry () {
    return setActiveRowFromEntry(profileStore.getInternalDefaultPSE());
}

/** @brief Get the row from the first level of the tree that match the provided name */
Gtk::TreeIter ProfileStoreComboBox::getRowFromLabel(Glib::ustring name) {
    Gtk::TreeIter row;
    Gtk::TreeModel::Children childs = refTreeModel->children();
    if (!name.empty()) {
        Gtk::TreeModel::Row currRow;
        for(Gtk::TreeModel::Children::iterator iter = childs.begin(); iter != childs.end(); ++iter) {
            currRow = *iter;
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
Gtk::TreeIter ProfileStoreComboBox::addRow(const ProfileStoreEntry *profileStoreEntry) {
    Gtk::TreeIter newEntry = refTreeModel->append();
    Gtk::TreeModel::Row row = *newEntry;
    row[methodColumns.label] = profileStoreEntry->label;
    row[methodColumns.profileStoreEntry] = profileStoreEntry;
    return newEntry;
}

