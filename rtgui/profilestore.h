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
#include "../rtengine/rtengine.h"
#include "threadutils.h"
#include "paramsedited.h"
#include <glibmm.h>

class ProfileStore {

        typedef enum {
            STORESTATE_NOTINITIALIZED,
            STORESTATE_BEINGINITIALIZED,
            STORESTATE_INITIALIZED,
            STORESTATE_DELETED
        } StoreState;

        MyMutex *parseMutex;
        StoreState storeState;
        std::map<Glib::ustring, rtengine::procparams::PartialProfile*> partProfiles;
        void parseDir (const Glib::ustring& pdir);
        void _parseProfiles ();

    public:

        ProfileStore();
        ~ProfileStore();
        bool init ();
        void parseProfiles ();
        const rtengine::procparams::PartialProfile* getProfile (const Glib::ustring& profname);
        std::vector<Glib::ustring>                  getProfileNames ();
        const rtengine::procparams::ProcParams*     getDefaultProcParams (bool isRaw);
        const rtengine::procparams::PartialProfile* getDefaultPartialProfile (bool isRaw);
};

extern ProfileStore profileStore;

#endif
