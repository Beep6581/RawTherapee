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
#include "../rtengine/imageformat.h"

class CacheImageData: public rtengine::FramesMetaData
{

public:

    // basic information
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
    unsigned short frameCount;
    double fnumber;
    double shutter;
    double focalLen, focalLen35mm;
    float focusDist;
    unsigned iso;
    bool isHDR;
    bool isPixelShift;
    int sensortype;
    rtengine::IIO_Sample_Format sampleFormat;
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
    // FramesMetaData interface
    //-------------------------------------------------------------------------

    unsigned int getRootCount () const { return -1; }
    unsigned int getFrameCount () const { return frameCount; }
    bool hasExif (unsigned int frame = 0) const  { return false; }
    rtexif::TagDirectory* getRootExifData (unsigned int root = 0) const { return nullptr; }
    rtexif::TagDirectory* getFrameExifData (unsigned int frame = 0) const { return nullptr; }
    rtexif::TagDirectory* getBestExifData (rtengine::ImageSource *imgSource, rtengine::procparams::RAWParams *rawParams) const { return nullptr; }
    bool hasIPTC (unsigned int frame = 0) const { return false; }
    rtengine::procparams::IPTCPairs getIPTCData (unsigned int frame = 0) const { return rtengine::procparams::IPTCPairs(); }
    tm getDateTime (unsigned int frame = 0) const { return tm{}; }
    time_t getDateTimeAsTS(unsigned int frame = 0) const { return time_t(-1); }
    int getISOSpeed (unsigned int frame = 0) const { return iso; }
    double getFNumber  (unsigned int frame = 0) const { return fnumber; }
    double getFocalLen (unsigned int frame = 0) const { return focalLen; }
    double getFocalLen35mm (unsigned int frame = 0) const { return focalLen35mm; }
    float getFocusDist (unsigned int frame = 0) const { return focusDist; }
    double getShutterSpeed (unsigned int frame = 0) const { return shutter; }
    double getExpComp (unsigned int frame = 0) const { return atof(expcomp.c_str()); }
    std::string getMake     (unsigned int frame = 0) const { return camMake; }
    std::string getModel    (unsigned int frame = 0) const { return camModel; }
    std::string getLens     (unsigned int frame = 0) const { return lens; }
    std::string getOrientation (unsigned int frame = 0) const { return ""; } // TODO
    bool getPixelShift (unsigned int frame = 0) const { return isPixelShift; }
    bool getHDR (unsigned int frame = 0) const { return isHDR; }
    std::string getRawType (unsigned int frame) const { return isPixelShift ? "PS" : isHDR ? "HDR" : "STD"; }
    rtengine::IIOSampleFormat getSampleFormat (unsigned int frame = 0) const { return sampleFormat; }
};
#endif
