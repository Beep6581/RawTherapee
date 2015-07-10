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
#ifndef _CACHEMANAGER_
#define _CACHEMANAGER_

#include <string>
#include <map>
#include <glibmm.h>
#include "thumbnail.h"
#include <cstdio>
#include "../rtengine/procparams.h"
#include "threadutils.h"

class Thumbnail;

class CacheManager {

        typedef std::pair<std::string, Thumbnail*> string_thumb_pair;
        typedef std::map<std::string, Thumbnail*> string_thumb_map;

        string_thumb_map openEntries;
        Glib::ustring    baseDir;
        MyMutex          mutex_;

        void deleteDir (const Glib::ustring& dirName);

        CacheManager () {}

    public:

        static CacheManager* getInstance(void);

        void        init        ();
        Thumbnail*  getEntry    (const Glib::ustring& fname);
        void        deleteEntry (const Glib::ustring& fname);
        void        renameEntry (const std::string& oldfilename, const std::string& oldmd5, const std::string& newfilename);

        void        closeThumbnail (Thumbnail* t);

        const Glib::ustring& getBaseDir     ()       { MyMutex::MyLock lock(mutex_); return baseDir; }
        void  closeCache ();

        static std::string getMD5 (const Glib::ustring& fname);

        void clearAll ();
        void clearThumbImages ();
        void clearProfiles ();
        void clearFromCache(const Glib::ustring& fname, bool leavenotrace);

        void applyCacheSizeLimitation ();

        Glib::ustring    getCacheFileName (const Glib::ustring& subdir, const Glib::ustring& fname, const Glib::ustring& md5);
};

#define cacheMgr CacheManager::getInstance()

#endif

