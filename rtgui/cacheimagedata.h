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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <glibmm/ustring.h>

#include "options.h"

#include "../rtengine/imageformat.h"
#include "../rtengine/rtengine.h"

class CacheImageData :
    public rtengine::FramesMetaData
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
    int rating;
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

    unsigned int getRootCount () const override { return -1; }
    unsigned int getFrameCount () const override { return frameCount; }
    bool hasExif (unsigned int frame = 0) const override  { return false; }
    rtexif::TagDirectory* getRootExifData (unsigned int root = 0) const override { return nullptr; }
    rtexif::TagDirectory* getFrameExifData (unsigned int frame = 0) const override { return nullptr; }
    rtexif::TagDirectory* getBestExifData (rtengine::ImageSource *imgSource, rtengine::procparams::RAWParams *rawParams) const override { return nullptr; }
    bool hasIPTC (unsigned int frame = 0) const override { return false; }
    rtengine::procparams::IPTCPairs getIPTCData (unsigned int frame = 0) const override;
    tm getDateTime (unsigned int frame = 0) const override { return tm{}; }
    time_t getDateTimeAsTS(unsigned int frame = 0) const override { return time_t(-1); }
    int getISOSpeed (unsigned int frame = 0) const override { return iso; }
    double getFNumber  (unsigned int frame = 0) const override { return fnumber; }
    double getFocalLen (unsigned int frame = 0) const override { return focalLen; }
    double getFocalLen35mm (unsigned int frame = 0) const override { return focalLen35mm; }
    float getFocusDist (unsigned int frame = 0) const override { return focusDist; }
    double getShutterSpeed (unsigned int frame = 0) const override { return shutter; }
    double getExpComp (unsigned int frame = 0) const override { return atof(expcomp.c_str()); }
    std::string getMake     (unsigned int frame = 0) const override { return camMake; }
    std::string getModel    (unsigned int frame = 0) const override { return camModel; }
    std::string getLens     (unsigned int frame = 0) const override { return lens; }
    std::string getOrientation (unsigned int frame = 0) const override { return ""; } // TODO
    int getRating (unsigned int frame = 0) const override { return rating; } // FIXME-piotr : missing rating
    bool getPixelShift () const override { return isPixelShift; }
    bool getHDR (unsigned int frame = 0) const override { return isHDR; }
    std::string getImageType (unsigned int frame) const override { return isPixelShift ? "PS" : isHDR ? "HDR" : "STD"; }
    rtengine::IIOSampleFormat getSampleFormat (unsigned int frame = 0) const override { return sampleFormat; }
};
