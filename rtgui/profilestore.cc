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
#include <profilestore.h>
#include <options.h>

ProfileStore profileStore;

using namespace rtengine;
using namespace rtengine::procparams;

extern Glib::ustring argv0;

void ProfileStore::parseProfiles () {

    // clear loaded profiles
    for (std::map<Glib::ustring,ProcParams*>::iterator i = pparams.begin(); i!=pparams.end(); i++)
        delete i->second;
    pparams.clear ();

    if (options.multiUser) {
        Glib::ustring userPD = options.rtdir + "/" + options.profilePath;
        if (!Glib::file_test (userPD, Glib::FILE_TEST_IS_DIR))
            g_mkdir_with_parents (userPD.c_str(), 511);
        parseDir (userPD);
    }
    parseDir (argv0 + "/" + options.profilePath);
}

void ProfileStore::parseDir (const Glib::ustring& pdir) {

  // reload the available profiles from the profile dir
  if (pdir!="") {
    // process directory
    Glib::ustring dirname = pdir;
    Glib::Dir* dir = NULL;
    try {
        dir = new Glib::Dir (dirname);
    }
    catch (const Glib::FileError& fe) {
        return;
    }
    dirname = dirname + "/";
    for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
      Glib::ustring fname = dirname + *i;
      Glib::ustring sname = *i;
      // ignore directories
      if (!Glib::file_test (fname, Glib::FILE_TEST_IS_DIR)) {
        int lastdot = sname.find_last_of ('.');
        if (lastdot!=Glib::ustring::npos && lastdot<=sname.size()-4 && !sname.casefold().compare (lastdot, 4, ".pp2")) {
          printf ("processing file %s...\n", fname.c_str());
          Glib::ustring name = sname.substr(0,lastdot);
          if (pparams.find(name)!=pparams.end()) {
            delete pparams[name];
            pparams.erase (pparams.find(name));
          }
          ProcParams* pp = new ProcParams ();
          int res = pp->load (fname);
          if (!res && pp->version>=220) 
            pparams[name] = pp;
          else
            delete pp;
        }
      }
    }
    delete dir;
  }
}

rtengine::procparams::ProcParams* ProfileStore::getProfile (const Glib::ustring& profname) {

    return pparams[profname];
}

std::vector<Glib::ustring> ProfileStore::getProfileNames () {

    std::vector<Glib::ustring> ret;
    for (std::map<Glib::ustring,ProcParams*>::iterator i = pparams.begin(); i!=pparams.end(); i++)
        ret.push_back (i->first);
    return ret;
}

rtengine::procparams::ProcParams* ProfileStore::getDefaultProcParams (bool isRaw) {

    if (!isRaw)
        return getProfile (options.defProfImg);
    else
        return getProfile (options.defProfRaw);
}

