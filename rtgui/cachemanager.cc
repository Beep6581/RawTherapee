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

CacheManager cacheMgr;

void CacheManager::init () {

    openEntries.clear ();
    baseDir = options.cacheBaseDir;

    if (!Glib::file_test (baseDir, Glib::FILE_TEST_IS_DIR))
        g_mkdir_with_parents (baseDir.c_str(), 511);
    if (!Glib::file_test (Glib::build_filename (baseDir, "profiles"), Glib::FILE_TEST_IS_DIR))
        g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "profiles")).c_str(), 511);
    if (!Glib::file_test (Glib::build_filename (baseDir, "images"), Glib::FILE_TEST_IS_DIR))
        g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "images")).c_str(), 511);
    if (!Glib::file_test (Glib::build_filename (baseDir, "aehistograms"), Glib::FILE_TEST_IS_DIR))
        g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "aehistograms")).c_str(), 511);
    if (!Glib::file_test (Glib::build_filename (baseDir, "embprofiles"), Glib::FILE_TEST_IS_DIR))
        g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "embprofiles")).c_str(), 511);
    if (!Glib::file_test (Glib::build_filename (baseDir, "data"), Glib::FILE_TEST_IS_DIR))
        g_mkdir_with_parents (Glib::ustring(Glib::build_filename (baseDir, "data")).c_str(), 511);
}

Thumbnail* CacheManager::getEntry (const Glib::ustring& fname) {

    Thumbnail* res = NULL;
    
    std::map<std::string, Thumbnail*>::iterator r = openEntries.find (fname);
    // if it is open, return it
    if (r!=openEntries.end()) {
        r->second->increaseRef ();
        return r->second;
    }

    // compute the md5
    std::string md5 = getMD5 (fname);
    if (md5=="")
        return NULL;

    // build path name
    Glib::ustring cfname = getCacheFileName ("data", fname, md5) + ".txt";

    // let's see if we have it in the cache
    if (Glib::file_test (cfname, Glib::FILE_TEST_EXISTS)) {
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

    if (res)
        openEntries[fname] = res;
    return res;
}


void CacheManager::deleteEntry (const Glib::ustring& fname) {

    // check if it is opened
    std::map<std::string, Thumbnail*>::iterator r = openEntries.find (fname);
    // if it is open, dont delete it
    if (r!=openEntries.end()) {
		std::string md5 = r->second->getMD5 ();
		r->second->decreaseRef ();
		// if in the editor, the thumbnail still exists. If not, delete it:
		r = openEntries.find (fname);
	    if (r==openEntries.end() && md5!="") {
			::g_remove ((getCacheFileName ("data", fname, md5) + ".txt").c_str());
			::g_remove ((getCacheFileName ("profiles", fname, md5) + ".pp2").c_str());
			::g_remove ((getCacheFileName ("images", fname, md5) + ".cust").c_str());
			::g_remove ((getCacheFileName ("images", fname, md5) + ".jpg").c_str());
			::g_remove ((getCacheFileName ("aehistograms", fname, md5)).c_str());
			::g_remove ((getCacheFileName ("embprofiles", fname, md5) + ".icc").c_str());
		}
	}
	else {
	    std::string md5 = getMD5 (fname);
	    if (md5!="") {
	        ::g_remove ((getCacheFileName ("data", fname, md5) + ".txt").c_str());
	        ::g_remove ((getCacheFileName ("profiles", fname, md5) + ".pp2").c_str());
	        ::g_remove ((getCacheFileName ("images", fname, md5) + ".cust").c_str());
	        ::g_remove ((getCacheFileName ("images", fname, md5) + ".jpg").c_str());
	        ::g_remove ((getCacheFileName ("aehistograms", fname, md5)).c_str());
	        ::g_remove ((getCacheFileName ("embprofiles", fname, md5) + ".icc").c_str());
	    }
	}
}


void CacheManager::renameEntry (const std::string& oldfilename, const std::string& oldmd5, const std::string& newfilename) {

    std::string newmd5 = getMD5 (newfilename);

    ::g_rename ((getCacheFileName ("profiles", oldfilename, oldmd5) + ".pp2").c_str(), (getCacheFileName ("profiles", newfilename, newmd5) + ".pp2").c_str());
    ::g_rename ((getCacheFileName ("images", oldfilename, oldmd5) + ".cust").c_str(), (getCacheFileName ("images", newfilename, newmd5) + ".cust").c_str());
    ::g_rename ((getCacheFileName ("images", oldfilename, oldmd5) + ".jpg").c_str(), (getCacheFileName ("images", newfilename, newmd5) + ".jpg").c_str());
    ::g_rename ((getCacheFileName ("aehistograms", oldfilename, oldmd5)).c_str(), (getCacheFileName ("aehistograms", newfilename, newmd5)).c_str());
    ::g_rename ((getCacheFileName ("embprofiles", oldfilename, oldmd5) + ".icc").c_str(), (getCacheFileName ("embprofiles", newfilename, newmd5) + ".icc").c_str());
    ::g_rename ((getCacheFileName ("data", oldfilename, oldmd5) + ".txt").c_str(), (getCacheFileName ("data", newfilename, newmd5) + ".txt").c_str());

    // check if it is opened
    std::map<std::string, Thumbnail*>::iterator r = openEntries.find (oldfilename);
    // if it is open, update md5
    if (r!=openEntries.end()) {
        Thumbnail* t = r->second;
        openEntries.erase (r);
        t->setFileName (newfilename);
        openEntries[newfilename] = t;
        t->updateCache ();
        t->reSaveThumbnail ();
    }
}

void CacheManager::closeThumbnail (Thumbnail* t) {

    t->updateCache ();
    std::map<std::string, Thumbnail*>::iterator r = openEntries.find (t->getFileName());
    if (r!=openEntries.end()) 
        openEntries.erase (r);
    delete t;
}

void CacheManager::closeCache () {

    applyCacheSizeLimitation ();
}

void CacheManager::clearAll () {

    deleteDir ("images");
    deleteDir ("aehistograms");
    deleteDir ("embprofiles");
    deleteDir ("profiles");
    deleteDir ("data");
    
    // re-generate thumbnail images and clear profiles of open thumbnails
    std::map<std::string, Thumbnail*>::iterator i;
    for (i=openEntries.begin(); i!=openEntries.end(); i++) {
        i->second->clearProcParams (CACHEMGR);
        i->second->generateThumbnailImage ();
        i->second->updateCache ();
    }
}
void CacheManager::clearThumbImages () {

    deleteDir ("images");
    deleteDir ("aehistograms");
    deleteDir ("embprofiles");
    
    // re-generate thumbnail images of open thumbnails
    std::map<std::string, Thumbnail*>::iterator i;
    for (i=openEntries.begin(); i!=openEntries.end(); i++)
        i->second->generateThumbnailImage ();
}

void CacheManager::clearProfiles () {

    deleteDir ("profiles");
    // clear profiles of open thumbnails
    std::map<std::string, Thumbnail*>::iterator i;
    for (i=openEntries.begin(); i!=openEntries.end(); i++)
        i->second->clearProcParams (CACHEMGR);
}

void CacheManager::deleteDir (const Glib::ustring& dirName) {

    try {
        Glib::Dir* dir = new Glib::Dir (Glib::build_filename (baseDir, dirName));
        for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) 
            ::g_remove (Glib::build_filename (Glib::build_filename (baseDir, dirName), *i).c_str());
        delete dir;
    }
    catch (const Glib::FileError& fe) {
    }
}

std::string CacheManager::getMD5 (const Glib::ustring& fname) {

    Glib::RefPtr<Gio::File> file = Gio::File::create_for_path (fname);
    if (file)	{
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
            ::g_remove ((Glib::build_filename (Glib::build_filename (baseDir, "data"), flist.front().fname) + ".txt").c_str());
            ::g_remove ((Glib::build_filename (Glib::build_filename (baseDir, "images"), flist.front().fname) + ".cust").c_str());
            ::g_remove ((Glib::build_filename (Glib::build_filename (baseDir, "images"), flist.front().fname) + ".jpg").c_str());
            ::g_remove ((Glib::build_filename (Glib::build_filename (baseDir, "aehistograms"), flist.front().fname)).c_str());
            ::g_remove ((Glib::build_filename (Glib::build_filename (baseDir, "embprofiles"), flist.front().fname) + ".icc").c_str());
//                ::g_remove ((Glib::build_filename (Glib::build_filename (baseDir, "profiles"), flist.front().fname) + ".pp2").c_str());
            flist.erase (flist.begin());
        }
    }
}

