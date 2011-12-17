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
#ifndef _CACHEIMAGEDATA_
#define _CACHEIMAGEDATA_

#include <glibmm.h>
#include "options.h"

class CacheImageData {

    public:

        // basic informations
        Glib::ustring  md5;
        Glib::ustring  version;
        bool  supported;
        ThFileType  format;
        char  rankOld; // old implementation of rank
        bool  inTrashOld; // old implementation of inTrash
        bool  recentlySaved;

        // time/date info
        bool  timeValid;
        short year;
        char  month;
        char  day;
        char  hour;
        char  min;
        char  sec;
        char  msec;

        // exif info
        bool  exifValid;
        double fnumber;
        double shutter;
        double focalLen;
        unsigned iso;
        Glib::ustring lens;
        Glib::ustring camera;
        Glib::ustring filetype;
        Glib::ustring expcomp;

        // additional info on raw images
        int   rotate;
        int   thumbImgType;
        int   thumbOffset;

		enum
		{
			FULL_THUMBNAIL = 0,  // was the thumbnail generated from whole file
			QUICK_THUMBNAIL = 1, // was rhe thumbnail generated from embedded jpeg
		};

        CacheImageData ();
        
        int load (const Glib::ustring& fname);
        int save (const Glib::ustring& fname);
};
#endif
