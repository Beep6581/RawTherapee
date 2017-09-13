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
#include "../rtengine/rtengine.h"

class CacheImageData: public rtengine::ImageMetaData
{

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
    // exif info
    bool  exifValid;
    double fnumber;
    double shutter;
    double focalLen, focalLen35mm;
    float focusDist;
    unsigned iso;
    Glib::ustring lens;
    Glib::ustring camMake;
    Glib::ustring camModel;
    Glib::ustring filetype;
    Glib::ustring expcomp;

    // store a copy of the autoWB's multipliers computed in Thumbnail::_generateThumbnailImage
    // they are not stored in the cache file by this class, but by rtengine::Thumbnail
    // -1 = Unknown
    double redAWBMul, greenAWBMul, blueAWBMul;

    // additional info on raw images
    int   rotate;
    int   thumbImgType;

    enum {
        FULL_THUMBNAIL = 0,  // was the thumbnail generated from whole file
        QUICK_THUMBNAIL = 1  // was the thumbnail generated from embedded jpeg
    };

    CacheImageData ();

    int load (const Glib::ustring& fname);
    int save (const Glib::ustring& fname);

    //-------------------------------------------------------------------------
    // ImageMetaData interface
    //-------------------------------------------------------------------------

    bool hasExif() const { return false; }
    const rtexif::TagDirectory *getExifData() const { return NULL; }
    bool hasIPTC() const { return false; }
    const rtengine::procparams::IPTCPairs getIPTCData () const { return rtengine::procparams::IPTCPairs(); }
    struct tm getDateTime () const { struct tm ret; return ret; }
    time_t getDateTimeAsTS() const { return time_t(-1); }
    int getISOSpeed() const { return iso; }
    double getFNumber() const { return fnumber; }
    double getFocalLen() const { return focalLen; }
    double getFocalLen35mm() const { return focalLen35mm; }
    float getFocusDist() const { return focusDist; }
    double getShutterSpeed() const { return shutter; }
    double getExpComp() const { return atof(expcomp.c_str()); }
    std::string getMake() const { return camMake; }
    std::string getModel() const { return camModel; }
    std::string getLens() const { return lens; }
    std::string getOrientation() const { return ""; } // TODO
};
#endif
