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
#include <cachemanager.h>
#include <options.h>
#include <glib/gstdio.h>
#include <giomm.h>
#include <guiutils.h>
#include <procparamchangers.h>
#include <safegtk.h>

CacheManager*
CacheManager::getInstance(void)
{
	static CacheManager* instance_ = 0;
	if ( instance_ == 0 )
	{
		static Glib::Mutex smutex_;
		Glib::Mutex::Lock lock(smutex_);
		if ( instance_ == 0 )
		{
			instance_ = new CacheManager();
		}
	}
	return instance_;
}

void CacheManager::init () {

    Glib::Mutex::Lock lock(mutex_);

    openEntries.clear ();
    baseDir = options.cacheBaseDir;

    if (!safe_file_test (baseDir, Glib::FILE_TEST_IS_DIR))
        safe_g_mkdir_with_parents (baseDir, 511);
    if (!safe_file_test (Glib::build_filename (baseDir, "profiles"), Glib::FILE_TEST_IS_DIR))
        safe_g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "profiles")), 511);
    if (!safe_file_test (Glib::build_filename (baseDir, "images"), Glib::FILE_TEST_IS_DIR))
        safe_g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "images")), 511);
    if (!safe_file_test (Glib::build_filename (baseDir, "aehistograms"), Glib::FILE_TEST_IS_DIR))
        safe_g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "aehistograms")), 511);
    if (!safe_file_test (Glib::build_filename (baseDir, "embprofiles"), Glib::FILE_TEST_IS_DIR))
        safe_g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "embprofiles")), 511);
    if (!safe_file_test (Glib::build_filename (baseDir, "data"), Glib::FILE_TEST_IS_DIR))
        safe_g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "data")), 511);
}

Thumbnail* CacheManager::getEntry (const Glib::ustring& fname) {

    Thumbnail* res = NULL;
    
    // take manager lock and search for entry, if found return it else release
    // lock and create it
	{
        Glib::Mutex::Lock lock(mutex_);

		string_thumb_map::iterator r = openEntries.find (fname);
		// if it is open, return it
		if (r!=openEntries.end()) {
			r->second->increaseRef ();
			return r->second;
		}
	}

    // compute the md5
    std::string md5 = getMD5 (fname);
    if (md5=="")
        return NULL;

    // build path name
    Glib::ustring cfname = getCacheFileName ("data", fname, md5) + ".txt";

    // let's see if we have it in the cache
    if (safe_file_test (cfname, Glib::FILE_TEST_EXISTS)) {
        CacheImageData* cfs = new CacheImageData ();
        int e = cfs->load (cfname);
        if (!e && cfs->supported==true) 
            res = new Thumbnail (this, fname, cfs);
        if (res && !res->isSupported ()) {
            delete res;
            res = NULL;
        }
        delete cfs;
    }

	// if not, create a new one
    if (!res) {
        res = new Thumbnail (this, fname, md5);
        if (!res->isSupported ()) {
            delete res;
            res = NULL;
        }
    }

    // retake the lock and see if it was added while we we're unlocked, if it
    // was use it over our version. if not added we create the cache entry
    if (res)
	{
		Glib::Mutex::Lock lock(mutex_);

		string_thumb_map::iterator r = openEntries.find (fname);
		if (r!=openEntries.end()) {
			delete res;
			r->second->increaseRef ();
			return r->second;
		}

		// it wasn't, create a new entry
		openEntries[fname] = res;
	}

    return res;
}


void CacheManager::deleteEntry (const Glib::ustring& fname) {

    Glib::Mutex::Lock lock(mutex_);

    // check if it is opened
    string_thumb_map::iterator r = openEntries.find (fname);
    // if it is open, dont delete it
    if (r!=openEntries.end()) {
		std::string md5 = r->second->getMD5 ();

        // decrease reference count; this will call back into CacheManager so
        // we release the lock for it.
		{
            lock.release();
			r->second->decreaseRef ();
            lock.acquire();
		}

		// if in the editor, the thumbnail still exists. If not, delete it:
		r = openEntries.find (fname);
	    if (r==openEntries.end() && md5!="") {
			safe_g_remove (getCacheFileName ("data", fname, md5) + ".txt");
			safe_g_remove (getCacheFileName ("profiles", fname, md5) + paramFileExtension);
			safe_g_remove (getCacheFileName ("images", fname, md5) + ".cust16");
			safe_g_remove (getCacheFileName ("images", fname, md5) + ".cust");
			safe_g_remove (getCacheFileName ("images", fname, md5) + ".jpg");
			safe_g_remove (getCacheFileName ("aehistograms", fname, md5));
			safe_g_remove (getCacheFileName ("embprofiles", fname, md5) + ".icc");
		}
	}
	else {
	    std::string md5 = getMD5 (fname);
	    if (md5!="") {
	        safe_g_remove (getCacheFileName ("data", fname, md5) + ".txt");
	        safe_g_remove (getCacheFileName ("profiles", fname, md5) + paramFileExtension);
	        safe_g_remove (getCacheFileName ("images", fname, md5) + ".cust16");
	        safe_g_remove (getCacheFileName ("images", fname, md5) + ".cust");
	        safe_g_remove (getCacheFileName ("images", fname, md5) + ".jpg");
	        safe_g_remove (getCacheFileName ("aehistograms", fname, md5));
	        safe_g_remove (getCacheFileName ("embprofiles", fname, md5) + ".icc");
	    }
	}
}


void CacheManager::renameEntry (const std::string& oldfilename, const std::string& oldmd5, const std::string& newfilename) {

    Glib::Mutex::Lock lock(mutex_);

    std::string newmd5 = getMD5 (newfilename);

    safe_g_rename (getCacheFileName ("profiles", oldfilename, oldmd5) + paramFileExtension, (getCacheFileName ("profiles", newfilename, newmd5) + paramFileExtension).c_str());
    safe_g_rename (getCacheFileName ("images", oldfilename, oldmd5) + ".cust16", getCacheFileName ("images", newfilename, newmd5) + ".cust16");
    safe_g_rename (getCacheFileName ("images", oldfilename, oldmd5) + ".cust", getCacheFileName ("images", newfilename, newmd5) + ".cust");
    safe_g_rename (getCacheFileName ("images", oldfilename, oldmd5) + ".jpg", getCacheFileName ("images", newfilename, newmd5) + ".jpg");
    safe_g_rename (getCacheFileName ("aehistograms", oldfilename, oldmd5), getCacheFileName ("aehistograms", newfilename, newmd5));
    safe_g_rename (getCacheFileName ("embprofiles", oldfilename, oldmd5) + ".icc", getCacheFileName ("embprofiles", newfilename, newmd5) + ".icc");
    safe_g_rename (getCacheFileName ("data", oldfilename, oldmd5) + ".txt", getCacheFileName ("data", newfilename, newmd5) + ".txt");

    // check if it is opened
    string_thumb_map::iterator r = openEntries.find (oldfilename);
    // if it is open, update md5
    if (r!=openEntries.end()) {
        Thumbnail* t = r->second;
        openEntries.erase (r);
        t->setFileName (newfilename);
        openEntries[newfilename] = t;
        t->updateCache ();
        t->saveThumbnail ();
    }
}

void CacheManager::closeThumbnail (Thumbnail* t) {

    Glib::Mutex::Lock lock(mutex_);

    t->updateCache ();
    string_thumb_map::iterator r = openEntries.find (t->getFileName());
    if (r!=openEntries.end()) 
        openEntries.erase (r);
    delete t;
}

void CacheManager::closeCache () {

    Glib::Mutex::Lock lock(mutex_);

    applyCacheSizeLimitation ();
}

void CacheManager::clearAll () {

    Glib::Mutex::Lock lock(mutex_);

    deleteDir ("images");
    deleteDir ("aehistograms");
    deleteDir ("embprofiles");
    deleteDir ("profiles");
    deleteDir ("data");
    
    // re-generate thumbnail images and clear profiles of open thumbnails
    //string_thumb_map::iterator i;
    //for (i=openEntries.begin(); i!=openEntries.end(); i++) {
    //    i->second->clearProcParams (CACHEMGR);
    //    i->second->generateThumbnailImage ();
    //    i->second->updateCache ();
    //}
}
void CacheManager::clearThumbImages () {

    Glib::Mutex::Lock lock(mutex_);

    deleteDir ("images");
    deleteDir ("aehistograms");
    deleteDir ("embprofiles");
}

void CacheManager::clearProfiles () {
    Glib::Mutex::Lock lock(mutex_);

    deleteDir ("profiles");
}

void CacheManager::deleteDir (const Glib::ustring& dirName) {

    try {
        Glib::Dir* dir = new Glib::Dir (Glib::build_filename (baseDir, dirName));
        for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) 
            safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, dirName), *i));
        delete dir;
    }
    catch (const Glib::FileError& fe) {
    }
}

std::string CacheManager::getMD5 (const Glib::ustring& fname) {

    Glib::RefPtr<Gio::File> file = Gio::File::create_for_path (fname);
    if (file && file->query_exists())	{
        Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info (file);
        if (info)
            return Glib::Checksum::compute_checksum (Glib::Checksum::CHECKSUM_MD5, Glib::ustring::compose ("%1%2", fname, info->get_size()));
    }
    return "";
}

Glib::ustring CacheManager::getCacheFileName (const Glib::ustring& subdir, const Glib::ustring& fname, const Glib::ustring& md5) {

    Glib::ustring cfn = Glib::build_filename (baseDir, subdir);
    Glib::ustring cname = Glib::path_get_basename (fname) + "." + md5;
    return Glib::build_filename (cfn, cname);
}

void CacheManager::applyCacheSizeLimitation () {

    std::vector<FileMTimeInfo> flist;
    Glib::ustring dataDir = Glib::build_filename (baseDir, "data");
    Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (dataDir);

    safe_build_file_list (dir, flist);

    if (flist.size() > options.maxCacheEntries) {
        std::sort (flist.begin(), flist.end());
        while (flist.size() > options.maxCacheEntries) {
            safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, "data"), flist.front().fname) + ".txt");
            safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, "images"), flist.front().fname) + ".cust16");
            safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, "images"), flist.front().fname) + ".cust");
            safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, "images"), flist.front().fname) + ".jpg");
            safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, "aehistograms"), flist.front().fname));
            safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, "embprofiles"), flist.front().fname) + ".icc");
            //                safe_g_remove (Glib::build_filename (Glib::build_filename (baseDir, "profiles"), flist.front().fname) + paramFileExtension);
            flist.erase (flist.begin());
        }
    }
}

