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

#include "dynamicprofile.h"
#include "../rtgui/options.h"
#include "../rtgui/multilangmgr.h"

using namespace rtengine;
using namespace rtengine::procparams;

ProfileStore::ProfileStore () : storeState (STORESTATE_NOTINITIALIZED), internalDefaultProfile (nullptr), internalDefaultEntry (nullptr), internalDynamicEntry (nullptr), loadAll (true)
{
    internalDefaultProfile = new AutoPartialProfile();
    internalDefaultProfile->set (true);
}

ProfileStore* ProfileStore::getInstance()
{
    static ProfileStore instance;
    return &instance;
}


bool ProfileStore::init (bool loadAll)
{
    if (storeState == STORESTATE_DELETED) {
        return false;
    }

    this->loadAll = loadAll;

    if ((storeState == STORESTATE_NOTINITIALIZED || storeState == STORESTATE_DIRTY) && loadAll) {
        storeState = STORESTATE_BEINGINITIALIZED;
        _parseProfiles ();
        storeState = STORESTATE_INITIALIZED;
    }

    return storeState == STORESTATE_INITIALIZED;
}

ProfileStore::~ProfileStore ()
{
    if (storeState == STORESTATE_NOTINITIALIZED) {
        return;
    }

    // This lock prevent object's suppression while scanning the directories
    storeState = STORESTATE_DELETED;

    {
        MyMutex::MyLock lock (parseMutex);

        clearProfileList ();
        partProfiles.clear ();
        clearFileList();
        delete internalDefaultProfile;
        delete internalDefaultEntry;
        delete internalDynamicEntry;
    }
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
void ProfileStore::parseProfilesOnce ()
{

    for (auto listener : listeners) {
        listener->storeCurrentValue();
    }

    init (true);  // safe even if already initialized

    for (auto listener : listeners) {
        listener->updateProfileList();
        listener->restoreValue();
    }
}

void ProfileStore::parseProfiles ()
{

    storeState = STORESTATE_DIRTY;
    parseProfilesOnce ();
}

void ProfileStore::_parseProfiles ()
{
    // clear loaded profiles
    folders.clear();
    clearFileList();
    clearProfileList ();

    folders.push_back ("<<< ROOT >>>"); // Fake path, so parentFolderId == 0 will be used to attach a ProfileStoreEntry to the root container, not sub-menu

    Glib::ustring p1 = options.getUserProfilePath();
    Glib::ustring p2 = options.getGlobalProfilePath();
    bool displayLevel0 = options.useBundledProfiles && !p1.empty() && !p2.empty() && p1 != p2;

    Glib::ustring virtualPath ("${U}");
    Glib::ustring currDir ("${U}");
    parseDir (p1, virtualPath, currDir, 0, 0, displayLevel0);

    if (displayLevel0) {
        virtualPath = "${G}";
        currDir = "${G}";
        parseDir (p2, virtualPath, currDir, 0, 0, displayLevel0);
    }

    // sort profiles
    std::sort (entries.begin(), entries.end(), SortProfiles() );

    // entries and partProfiles are empty, but the entry and profiles already exist (they have survived to clearFileList and clearProfileList)
    if (!internalDefaultEntry) {
        internalDefaultEntry = new ProfileStoreEntry (Glib::ustring ("(") + M ("PROFILEPANEL_PINTERNAL") + Glib::ustring (")"), PSET_FILE, 0, 0);
    }

    entries.push_back (internalDefaultEntry);
    partProfiles[internalDefaultEntry] = internalDefaultProfile;

    if (!internalDynamicEntry) {
        internalDynamicEntry = new ProfileStoreEntry (Glib::ustring ("(") + M ("PROFILEPANEL_PDYNAMIC") + Glib::ustring (")"), PSET_FILE, 0, 0);
        // do not add it to the entries. This is here only for the preferences dialog
    }

    // Check if the default profiles has been found.
    if (findEntryFromFullPathU (options.defProfRaw) == nullptr) {
        options.setDefProfRawMissing (true);

        if (options.rtSettings.verbose) {
            printf ("WARNING: Default profile \"%s\" for raw images not found!\n", options.defProfRaw.c_str());
        }
    }

    if (findEntryFromFullPathU (options.defProfImg) == nullptr) {
        options.setDefProfImgMissing (true);

        if (options.rtSettings.verbose) {
            printf ("WARNING: Default profile \"%s\" for standard images not found!\n", options.defProfImg.c_str());
        }
    }
}

/// @return Returns true if some files has been found (directories are ignored)
bool ProfileStore::parseDir (Glib::ustring& realPath, Glib::ustring& virtualPath, Glib::ustring& currDir, unsigned int parentId, unsigned char level, bool displayLevel0)
{
    bool fileFound = false;

    // reload the available profiles from the profile dir
    if (!realPath.empty() && Glib::file_test (realPath, Glib::FILE_TEST_EXISTS) && Glib::file_test (realPath, Glib::FILE_TEST_IS_DIR)) {
        unsigned int folder = 0; // folder's own Id

        // add this entry to the folder list
        folders.push_back (virtualPath);
        folder = (unsigned int) (folders.size()) - 1;

        if (level > 0 || displayLevel0) {
            // replace the virtual folder name by a localized text
            if (currDir == "${U}") {
                currDir = M ("PROFILEPANEL_MYPROFILES");
            } else if (currDir == "${G}") {
                currDir = M ("PROFILEPANEL_GLOBALPROFILES");
            }

            // add this localized text to the file list
            entries.push_back ( new ProfileStoreEntry (currDir, PSET_FOLDER, parentId, folder) );
        }

        // walking through the directory
        Glib::Dir* dir = nullptr;
        dir = new Glib::Dir (realPath);

        for (Glib::DirIterator i = dir->begin(); i != dir->end(); ++i) {
            currDir = *i;

            if (currDir == "." || currDir == "..") {
                continue;
            }

            Glib::ustring fname = Glib::build_filename (realPath, currDir);

            if (Glib::file_test (fname, Glib::FILE_TEST_IS_DIR)) {
                Glib::ustring vp (Glib::build_filename (virtualPath, currDir));
                Glib::ustring rp (Glib::build_filename (realPath,    currDir));
                fileFound = parseDir (rp, vp, currDir, folder, level + 1, 0);
            } else {
                size_t lastdot = currDir.find_last_of ('.');

                if (lastdot != Glib::ustring::npos && lastdot == currDir.length() - 4 && currDir.substr (lastdot).casefold() == paramFileExtension) {
                    // file found
                    if ( options.rtSettings.verbose ) {
                        printf ("Processing file %s...", fname.c_str());
                    }

                    Glib::ustring name = currDir.substr (0, lastdot);

                    // create the partial profile
                    AutoPartialProfile *pProf = new AutoPartialProfile();
                    int res = pProf->load (fname);

                    if (!res && pProf->pparams->ppVersion >= 220) {
                        fileFound = true;

                        if ( options.rtSettings.verbose ) {
                            printf ("OK\n");
                        }

                        // adding this file to the list
                        ProfileStoreEntry* filePSE = new ProfileStoreEntry (name, PSET_FILE, folder, 0);
                        entries.push_back (filePSE);

                        // map the partial profile
                        partProfiles[filePSE] = pProf;
                        //partProfiles.insert( std::pair<ProfileStoreEntry*, rtengine::procparams::AutoPartialProfile*> (filePSE, pProf) );
                    } else if ( options.rtSettings.verbose ) {
                        printf ("failed!\n");
                    }
                }
            }
        }

        delete dir;
    }

    if (!fileFound && (level > 0 || displayLevel0)) {
        // no files found in this level, we delete the subdirectory entry
        folders.pop_back();

        delete entries.back();
        entries.pop_back();
    }

    return fileFound;
}

int ProfileStore::findFolderId (const Glib::ustring &path)
{
    // initialization must have been done when calling this
    for (std::vector<Glib::ustring>::iterator i = folders.begin(); i != folders.end(); ++i) {
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
const ProfileStoreEntry* ProfileStore::findEntryFromFullPathU (Glib::ustring path)
{
    if (path.empty()) {
        return nullptr;
    }

    if (storeState == STORESTATE_NOTINITIALIZED) {
        parseProfilesOnce();
    }

    if (path == DEFPROFILE_INTERNAL || path == DEFPROFILE_DYNAMIC) {
        return internalDefaultEntry;
    }

    // consistently apply casefold() to make sure dot position is correct
    const Glib::ustring casefolded_path = path.casefold();
    const Glib::ustring::size_type lastdot_pos = casefolded_path.find_last_of ('.');

    if (
        lastdot_pos != Glib::ustring::npos
        && lastdot_pos <= casefolded_path.size() - 4
        && !casefolded_path.compare (lastdot_pos, 4, paramFileExtension)
    ) {
        // removing the extension
        // now use dot position without casefold()
        path = path.substr (0, path.find_last_of ('.'));
    }

    // dir separator may come from options file and may be \ or /, we convert them to G_DIR_SEPARATOR_S
    if (path.size() > 4 && (path[4] == '/' || path[4] == '\\')) {
        path = path.substr (0, 4) + G_DIR_SEPARATOR_S + path.substr (5);
    }

    // removing the filename
    Glib::ustring fName = Glib::path_get_basename (path);

    if (!fName.empty()) {
        path = path.substr (0, path.length() - fName.length());
    } else {
        // path is malformed, returning NULL;
        return nullptr;
    }

    path = Glib::path_get_dirname (path);

    // 1. find the path in the folder list
    int parentFolderId = findFolderId (path);

    if (parentFolderId == -1) {
        return nullptr;
    }

    // 2. find the entry that match the given filename and parentFolderId
    if (parentFolderId >= 0) {
        for (auto entry : entries) {
            if (entry->parentFolderId == parentFolderId  &&  entry->label == fName) {
                return entry;
            }
        }
    }

    return nullptr;
}

/** Protected version of findEntryFromFullPathU */
const ProfileStoreEntry* ProfileStore::findEntryFromFullPath (Glib::ustring path)
{
    MyMutex::MyLock lock (parseMutex);
    return findEntryFromFullPathU (path);
}

const PartialProfile* ProfileStore::getProfile (Glib::ustring path)
{

    if (storeState == STORESTATE_NOTINITIALIZED) {
        parseProfilesOnce();
    }

    const ProfileStoreEntry *pse = findEntryFromFullPath (path);

    if (!pse) {
        return nullptr;
    }

    return getProfile (pse);
}

const PartialProfile* ProfileStore::getProfile (const ProfileStoreEntry* entry)
{

    if (storeState == STORESTATE_NOTINITIALIZED) {
        parseProfilesOnce();
    }

    MyMutex::MyLock lock (parseMutex);

    if (entry == internalDefaultEntry) {
        return internalDefaultProfile;
    }

    std::map<const ProfileStoreEntry*, rtengine::procparams::AutoPartialProfile*>::iterator iter = partProfiles.find (entry);

    if (iter != partProfiles.end()) {
        return iter->second;
    } else {
        // This shouldn't happen!
#ifndef NDEBUG
        printf ("WARNING! Profile not found!\n");
#endif
        return nullptr;
    }
}

/** @brief Get a pointer to the profile's vector list
 *
 * This method grants you unique access to the vector list through Mutex locking.
 * When you're done with the file list, you MUST call the releaseFileList method to release the lock.
 */
const std::vector<const ProfileStoreEntry*>* ProfileStore::getFileList ()
{

    if (storeState == STORESTATE_NOTINITIALIZED) {
        parseProfilesOnce();
    }

    parseMutex.lock();

    return &entries;
}

void ProfileStore::releaseFileList()
{
    parseMutex.unlock();
}

/*
 * Send back a pointer to the default procparams for raw or standard images.
 * If the profile doesn't already exist in the profile list,
 * it will add it with default internal values, so this method never fails
 */
const ProcParams* ProfileStore::getDefaultProcParams (bool isRaw)
{

    //Note: the mutex is locked in getProfile, called below
    //      eventual initialization is done there too

    const PartialProfile* pProf = getProfile (isRaw ? options.defProfRaw : options.defProfImg);

    if (!pProf) {
        pProf = internalDefaultProfile;
    }

    return pProf->pparams;
}

/*
 * Send back a pointer to the default partial profile for raw or standard images.
 * If it doesn't already exist in the profile list, it will add it with default internal values,
 * so this method will never fails
 */
const PartialProfile* ProfileStore::getDefaultPartialProfile (bool isRaw)
{

    //Note: the mutex is locked in getProfile, called below
    //      eventual initialization is done there too

    const PartialProfile* pProf = getProfile (isRaw ? options.defProfRaw : options.defProfImg);

    if (!pProf) {
        pProf = internalDefaultProfile;
    }

    return pProf;
}

const Glib::ustring ProfileStore::getPathFromId (int folderId)
{
    // initialization must have been done when calling this
    return folders.at (folderId);
}


void ProfileStore::clearFileList()
{
    for (auto entry : entries) {
        if (entry != internalDefaultEntry) {
            delete entry;
        }
    }

    entries.clear();
}

void ProfileStore::clearProfileList()
{
    for (auto partProfile : partProfiles) {
        if (partProfile.second != internalDefaultProfile) {
            delete partProfile.second;
        }
    }

    partProfiles.clear();
}

void ProfileStore::addListener (ProfileStoreListener *listener)
{
    listeners.push_back (listener);
}

void ProfileStore::removeListener (ProfileStoreListener *listener)
{
    listeners.remove (listener);
}

void ProfileStore::dumpFolderList()
{
    printf ("Folder list:\n------------\n");

    for (unsigned int i = 0; i < folders.size(); i++) {
        printf (" #%3ud - %s\n", i, folders.at (i).c_str());
    }

    printf ("\n");
}

PartialProfile *ProfileStore::loadDynamicProfile (const ImageMetaData *im)
{
    if (storeState == STORESTATE_NOTINITIALIZED) {
        parseProfilesOnce();
    }

    PartialProfile *ret = new PartialProfile (true, true);

    if (!rulesLoaded) {
        loadRules();
    }

    for (auto rule : dynamicRules) {
        if (rule.matches (im)) {
            if (options.rtSettings.verbose) {
                printf ("found matching profile %s\n", rule.profilepath.c_str());
            }

            const PartialProfile *p = getProfile (rule.profilepath);

            if (p != nullptr) {
                p->applyTo (ret->pparams);
            } else {
                printf ("ERROR loading matching profile from: %s\n", rule.profilepath.c_str());
            }
        }
    }

    return ret;
}

ProfileStoreEntry::ProfileStoreEntry() : label (""), type (PSET_FOLDER), parentFolderId (0), folderId (0) {}

ProfileStoreEntry::ProfileStoreEntry (Glib::ustring label, PSEType type, unsigned short parentFolder, unsigned short folder) : label (label), type (type), parentFolderId (parentFolder), folderId (folder) {}

void ProfileStoreEntry::setValues (Glib::ustring label, PSEType type, unsigned short parentFolder, unsigned short folder)
{
    this->label = label;
    this->type = type;
    parentFolderId = parentFolder;
    folderId = folder;
}

