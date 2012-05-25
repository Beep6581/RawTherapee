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
#include "../rtengine/safegtk.h"

ProfileStore profileStore;

using namespace rtengine;
using namespace rtengine::procparams;

ProfileStore::ProfileStore () {
    storeState = STORESTATE_NOTINITIALIZED;
    parseMutex = NULL;
}

bool ProfileStore::init () {
    if (storeState == STORESTATE_DELETED)
        return false;
    if (storeState == STORESTATE_NOTINITIALIZED) {
        storeState = STORESTATE_BEINGINITIALIZED;
        parseMutex = new Glib::Mutex();
        _parseProfiles ();
        storeState = STORESTATE_INITIALIZED;
    }
	return true;
}

ProfileStore::~ProfileStore () {

    // This lock prevent object's suppression while scanning the directories
    storeState = STORESTATE_DELETED;
    Glib::Mutex::Lock lock(*parseMutex);

    for (std::map<Glib::ustring,PartialProfile*>::iterator i = partProfiles.begin(); i!=partProfiles.end(); i++) {
        if (i->second->pparams) delete i->second->pparams;
        if (i->second->pedited) delete i->second->pedited;
        delete i->second;
    }
    partProfiles.clear ();
    lock.release();
    delete parseMutex;
    parseMutex = NULL;
}

/*
 * Public method to parse the profiles' directories
 * Since there's a race condition in the multithreaded environment on this object,
 * parseProfiles may need to ask for initialization of this object, and then will
 * ask a mutex lock on it, has it been initialized by this call or not
 */
void ProfileStore::parseProfiles () {

    if (!init())
        // I don't even know if this situation can occur
        return;
    Glib::Mutex::Lock lock(*parseMutex);

    _parseProfiles ();
}

void ProfileStore::_parseProfiles () {

    // clear loaded profiles
    for (std::map<Glib::ustring,PartialProfile*>::iterator i = partProfiles.begin(); i!=partProfiles.end(); i++) {
        delete i->second->pparams;
        delete i->second->pedited;
        delete i->second;
    }
    partProfiles.clear ();

    parseDir (options.getUserProfilePath());
    parseDir (options.getGlobalProfilePath());
}

void ProfileStore::parseDir (const Glib::ustring& pdir) {

  // reload the available profiles from the profile dir
  if (pdir!="" && safe_file_test(pdir, Glib::FILE_TEST_EXISTS) && safe_file_test(pdir, Glib::FILE_TEST_IS_DIR)) {
    // process directory
    Glib::ustring dirname = pdir;
    Glib::Dir* dir = NULL;
    dir = new Glib::Dir (dirname);
    dirname = dirname + "/";
    for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
      Glib::ustring fname = dirname + *i;
      Glib::ustring sname = *i;
      // ignore directories
      if (!safe_file_test (fname, Glib::FILE_TEST_IS_DIR)) {
        size_t lastdot = sname.find_last_of ('.');
        if (lastdot!=Glib::ustring::npos && lastdot<=sname.size()-4 && !sname.casefold().compare (lastdot, 4, paramFileExtension)) {
          if( options.rtSettings.verbose )
            printf ("Processing file %s...\n", fname.c_str());
          Glib::ustring name = sname.substr(0,lastdot);
          std::map<Glib::ustring,PartialProfile*>::iterator j = partProfiles.find(name);
          if (j!=partProfiles.end()) {
            j->second->deleteInstance();
            delete j->second;
            partProfiles.erase (j);
          }
          PartialProfile* pProf = new PartialProfile (true);
          int res = pProf->load (fname);
          if (!res && pProf->pparams->ppVersion>=220) {
            partProfiles[name] = pProf;
          }
          else {
            pProf->deleteInstance();
            delete pProf;
          }
        }
      }
    }
    delete dir;
  }
  // Check if the default profiles has been found. If no, create default instance
  // This operation is safe: if the profile is finally found in another directory, the profile will be updated
  if (partProfiles.find(options.defProfRaw) == partProfiles.end()) {
    PartialProfile* pProf = new PartialProfile (true);
    pProf->set(true);
    partProfiles[options.defProfRaw] = pProf;
  }
  if (partProfiles.find(options.defProfImg) == partProfiles.end()) {
    PartialProfile* pProf = new PartialProfile (true);
    pProf->set(true);
    partProfiles[options.defProfImg] = pProf;
  }
}

PartialProfile* ProfileStore::getProfile (const Glib::ustring& profname) {

    if (!init())
        // I don't even know if this situation can occur
        return NULL;
    Glib::Mutex::Lock lock(*parseMutex);

    if (partProfiles.find(profname) != partProfiles.end()) {
        return partProfiles[profname];
    }
    else {
        return NULL;
    }
}

std::vector<Glib::ustring> ProfileStore::getProfileNames () {

    std::vector<Glib::ustring> ret;

    if (!init())
        // I don't even know if this situation can occur
        return ret;
    Glib::Mutex::Lock lock(*parseMutex);

    for (std::map<Glib::ustring,PartialProfile*>::iterator i = partProfiles.begin(); i!=partProfiles.end(); i++)
        ret.push_back (i->first);
    return ret;
}

/*
 * Send back a pointer to the default procparams for raw or standard images.
 * If the profile doesn't already exist in the profile list,
 * it will add it with default internal values, so this method never fails
 */
ProcParams* ProfileStore::getDefaultProcParams (bool isRaw) {

    if (!init())
        // I don't even know if this situation can occur
        return NULL;
    //Note: the mutex is locked in getProfile, called below

    PartialProfile* pProf = getProfile (isRaw ? options.defProfRaw : options.defProfImg);
    // NOTE: pProf should not be NULL anymore, since init() should have created the default profiles already, but the code is left as is
    if (!pProf) {
        pProf = new PartialProfile (true);
        pProf->set(true);
        partProfiles[DEFPROFILE_INTERNAL] = pProf;
    }
    return pProf->pparams;
}

