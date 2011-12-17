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
#include "cacheimagedata.h"
#include <vector>
#include <glib/gstdio.h>
#include "../rtengine/safekeyfile.h"
#include "../rtengine/safegtk.h"
#include "version.h"

CacheImageData::CacheImageData () 
    : md5(""), supported(false), format(FT_Invalid), rankOld(-1), inTrashOld(false), recentlySaved(false),
    timeValid(false), exifValid(false), thumbImgType(0) {
}

int CacheImageData::load (const Glib::ustring& fname) {

    rtengine::SafeKeyFile keyFile;
    
    try {
        if (!keyFile.load_from_file (fname)) 
            return 1;

        if (keyFile.has_group ("General")) { 
            if (keyFile.has_key ("General", "MD5"))             md5         = keyFile.get_string ("General", "MD5");
            if (keyFile.has_key ("General", "Version"))         version     = keyFile.get_string ("General", "Version");
            if (keyFile.has_key ("General", "Supported"))       supported   = keyFile.get_boolean ("General", "Supported");
            if (keyFile.has_key ("General", "Format"))          format      = (ThFileType)keyFile.get_integer ("General", "Format");
            if (keyFile.has_key ("General", "Rank"))            rankOld     = keyFile.get_integer ("General", "Rank");
            if (keyFile.has_key ("General", "InTrash"))         inTrashOld  = keyFile.get_boolean ("General", "InTrash");
            if (keyFile.has_key ("General", "RecentlySaved"))   recentlySaved = keyFile.get_boolean ("General", "RecentlySaved");
        }

        timeValid = keyFile.has_group ("DateTime");

        if (timeValid) { 
            if (keyFile.has_key ("DateTime", "Year"))   year    = keyFile.get_integer ("DateTime", "Year");
            if (keyFile.has_key ("DateTime", "Month"))  month   = keyFile.get_integer ("DateTime", "Month");
            if (keyFile.has_key ("DateTime", "Day"))    day     = keyFile.get_integer ("DateTime", "Day");
            if (keyFile.has_key ("DateTime", "Hour"))   hour    = keyFile.get_integer ("DateTime", "Hour");
            if (keyFile.has_key ("DateTime", "Min"))    min     = keyFile.get_integer ("DateTime", "Min");
            if (keyFile.has_key ("DateTime", "Sec"))    sec     = keyFile.get_integer ("DateTime", "Sec");
            if (keyFile.has_key ("DateTime", "MSec"))   msec    = keyFile.get_integer ("DateTime", "MSec");
        }

		exifValid = false;
        
		if (keyFile.has_group ("ExifInfo")) {
			exifValid = true;
			if (keyFile.has_key ("ExifInfo", "Valid")) exifValid = keyFile.get_boolean ("ExifInfo", "Valid");

			if (exifValid) { 
				if (keyFile.has_key ("ExifInfo", "FNumber"))    fnumber     = keyFile.get_double ("ExifInfo", "FNumber");
				if (keyFile.has_key ("ExifInfo", "Shutter"))    shutter     = keyFile.get_double ("ExifInfo", "Shutter");
				if (keyFile.has_key ("ExifInfo", "FocalLen"))   focalLen    = keyFile.get_double ("ExifInfo", "FocalLen");
				if (keyFile.has_key ("ExifInfo", "ISO"))        iso         = keyFile.get_integer ("ExifInfo", "ISO");
				if (keyFile.has_key ("ExifInfo", "ExpComp"))    expcomp     = keyFile.get_string ("ExifInfo", "ExpComp");
			}
			if (keyFile.has_key ("ExifInfo", "Lens"))       lens        = keyFile.get_string ("ExifInfo", "Lens");
			if (keyFile.has_key ("ExifInfo", "Camera"))     camera      = keyFile.get_string ("ExifInfo", "Camera");
		}
		
		if (keyFile.has_group ("FileInfo")) {
			if (keyFile.has_key ("FileInfo", "Filetype"))   filetype    = keyFile.get_string ("FileInfo", "Filetype");
		}
			
			
        if (format==FT_Raw && keyFile.has_group ("ExtraRawInfo")) {
            if (keyFile.has_key ("ExtraRawInfo", "ThumbImageType"))     thumbImgType    = keyFile.get_integer ("ExtraRawInfo", "ThumbImageType");
            if (keyFile.has_key ("ExtraRawInfo", "ThumbImageOffset"))   thumbOffset     = keyFile.get_integer ("ExtraRawInfo", "ThumbImageOffset");
        }
        else {
            rotate = 0;
            thumbImgType = 0;
        }
        return 0;
    }
    catch (Glib::Error) {
        return 1;
    }
}

int CacheImageData::save (const Glib::ustring& fname) {

    rtengine::SafeKeyFile keyFile;
    
    if (safe_file_test(fname,Glib::FILE_TEST_EXISTS)) keyFile.load_from_file (fname); 

    keyFile.set_string  ("General", "MD5", md5);
    keyFile.set_string  ("General", "Version", VERSION); // Application's version
    keyFile.set_boolean ("General", "Supported", supported);
    keyFile.set_integer ("General", "Format", format);
    keyFile.set_boolean ("General", "RecentlySaved", recentlySaved);

    // remove the old implementation of Rank and InTrash from cache
    if (keyFile.has_key ("General", "Rank")) keyFile.remove_key("General", "Rank");
    if (keyFile.has_key ("General", "InTrash")) keyFile.remove_key("General", "InTrash");

    if (timeValid) { 
        keyFile.set_integer ("DateTime", "Year", year);
        keyFile.set_integer ("DateTime", "Month", month);
        keyFile.set_integer ("DateTime", "Day", day);
        keyFile.set_integer ("DateTime", "Hour", hour);
        keyFile.set_integer ("DateTime", "Min", min);
        keyFile.set_integer ("DateTime", "Sec", sec);
        keyFile.set_integer ("DateTime", "MSec", msec);
    }

    keyFile.set_boolean  ("ExifInfo", "Valid", exifValid);
    if (exifValid) { 
        keyFile.set_double  ("ExifInfo", "FNumber", fnumber);
        keyFile.set_double  ("ExifInfo", "Shutter", shutter);
        keyFile.set_double  ("ExifInfo", "FocalLen", focalLen);
        keyFile.set_integer ("ExifInfo", "ISO", iso);
        keyFile.set_string  ("ExifInfo", "ExpComp", expcomp);
    }
    keyFile.set_string  ("ExifInfo", "Lens", lens);
    keyFile.set_string  ("ExifInfo", "Camera", camera);
    keyFile.set_string  ("FileInfo", "Filetype", filetype);

    if (format==FT_Raw) {
        keyFile.set_integer ("ExtraRawInfo", "ThumbImageType", thumbImgType);
        keyFile.set_integer ("ExtraRawInfo", "ThumbImageOffset", thumbOffset);
    }

    FILE *f = safe_g_fopen (fname, "wt");
    if (!f)
        return 1;
    else {
        fprintf (f, "%s", keyFile.to_data().c_str());
        fclose (f);
        return 0;
    }}

