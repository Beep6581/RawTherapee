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
#ifndef __IMAGEDATA_H__
#define __IMAGEDATA_H__

#include <cstdio>
#include <memory>
#include <string>

#include <glibmm.h>

#include <libiptcdata/iptc-data.h>

#include "../rtexif/rtexif.h"

#include "procparams.h"
#include "rawimage.h"
#include "rtengine.h"

namespace rtengine
{

class FrameData
{

protected:
    rtexif::TagDirectory* root;
    IptcData* iptc;

    struct tm time;
    time_t timeStamp;
    int iso_speed;
    double aperture;
    double focal_len, focal_len35mm;
    float focus_dist;  // dist: 0=unknown, 10000=infinity
    double shutter;
    double expcomp;
    std::string make, model, serial;
    std::string orientation;
    std::string lens;
    IIOSampleFormat sampleFormat;

    // each frame has the knowledge of "being an"
    // or "being part of an" HDR or PS image
    bool isPixelShift;
    bool isHDR;

    void extractInfo ();

public:

    FrameData ();
    FrameData (rtexif::ExifManager &exifManager);
    virtual ~FrameData ();

    bool getPixelShift () const;
    bool getHDR () const;
    IIOSampleFormat getSampleFormat () const;
    rtexif::TagDirectory* getExifData () const;
    procparams::IPTCPairs getIPTCData () const;
    bool hasExif () const;
    bool hasIPTC () const;
    tm getDateTime () const;
    time_t getDateTimeAsTS () const;
    int getISOSpeed () const;
    double getFNumber () const;
    double getFocalLen () const;
    double getFocalLen35mm () const;
    float getFocusDist () const;
    double getShutterSpeed () const;
    double getExpComp  () const;
    std::string getMake () const;
    std::string getModel () const;
    std::string getLens () const;
    std::string getSerialNumber () const;
    std::string getOrientation () const;
};

class RawFrameData : public FrameData
{
public:
    RawFrameData (rtexif::ExifManager &exifManager);
};

class JpegFrameData : public FrameData
{
public:
    JpegFrameData (rtexif::ExifManager &exifManager);
};

class TiffFrameData : public FrameData
{
public:
    TiffFrameData (rtexif::ExifManager &exifManager);
};

class FramesData : public FramesMetaData {
private:
    std::vector<FrameData*> frames;
    unsigned int dcrawFrameCount;

public:
    FramesData (const Glib::ustring& fname, std::unique_ptr<RawMetaDataLocation> rml = nullptr, bool firstFrameOnly = false, bool loadAll = false);
    ~FramesData ();

    void setDCRawFrameCount (unsigned int frameCount);
    unsigned int getFrameCount () const;
    FrameData *getFrameData (int frame) const;
    bool getPixelShift (unsigned int frame = 0) const;
    bool getHDR (unsigned int frame = 0) const;
    IIOSampleFormat getSampleFormat (unsigned int frame = 0) const;
    rtexif::TagDirectory* getExifData (unsigned int frame = 0) const;
    procparams::IPTCPairs getIPTCData (unsigned int frame = 0) const;
    bool hasExif (unsigned int frame = 0) const;
    bool hasIPTC (unsigned int frame = 0) const;
    tm getDateTime (unsigned int frame = 0) const;
    time_t getDateTimeAsTS (unsigned int frame = 0) const;
    int getISOSpeed (unsigned int frame = 0) const;
    double getFNumber (unsigned int frame = 0) const;
    double getFocalLen (unsigned int frame = 0) const;
    double getFocalLen35mm (unsigned int frame = 0) const;
    float getFocusDist (unsigned int frame = 0) const;
    double getShutterSpeed (unsigned int frame = 0) const;
    double getExpComp (unsigned int frame = 0) const;
    std::string getMake (unsigned int frame = 0) const;
    std::string getModel (unsigned int frame = 0) const;
    std::string getLens (unsigned int frame = 0) const;
    std::string getSerialNumber (unsigned int frame = 0) const;
    std::string getOrientation (unsigned int frame = 0) const;
};


}
#endif
