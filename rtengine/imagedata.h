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
#include "rawimage.h"
#include <string>
#include <glibmm.h>
#include "../rtexif/rtexif.h"
#include "procparams.h"
#include <libiptcdata/iptc-data.h>
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
    int isHDR; // Number of frame

    void extractInfo ();

public:

    FrameData ();
    FrameData (rtexif::ExifManager &exifManager);
    virtual ~FrameData ();

    bool getPixelShift () const
    {
        return isPixelShift;
    }
    int getHDR () const
    {
        return isHDR;
    }

    IIOSampleFormat getSampleFormat () const
    {
        return sampleFormat;
    }

    const rtexif::TagDirectory*   getExifData () const
    {
        return root;
    }
    const procparams::IPTCPairs   getIPTCData () const;

    bool hasExif () const
    {
        return root && root->getCount();
    }
    bool hasIPTC () const
    {
        return iptc;
    }

    struct tm   getDateTime () const {
        return time;
    }
    time_t      getDateTimeAsTS() const
    {
        return timeStamp;
    }
    int         getISOSpeed () const
    {
        return iso_speed;
    }
    double      getFNumber  () const
    {
        return aperture;
    }
    double      getFocalLen () const
    {
        return focal_len;
    }
    double      getFocalLen35mm () const
    {
        return focal_len35mm;
    }
    float       getFocusDist () const
    {
        return focus_dist;
    }
    double      getShutterSpeed () const
    {
        return shutter;
    }
    double      getExpComp  () const
    {
        return expcomp;
    }
    std::string getMake     () const
    {
        return make;
    }
    std::string getModel    () const
    {
        return model;
    }
    std::string getLens     () const
    {
        return lens;
    }
    std::string getSerialNumber () const
    {
        return serial;
    }
    std::string getOrientation () const
    {
        return orientation;
    }
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
    int dcrawFrameCount;

public:
    FramesData (Glib::ustring fname, RawMetaDataLocation* rml = nullptr, bool firstFrameOnly = false, bool loadAll = false);
    ~FramesData ();

    void setDCRawFrameCount (int frameCount)
    {
        dcrawFrameCount = frameCount;
    }

    int getFrameCount () const
    {
        return dcrawFrameCount ? dcrawFrameCount : frames.size();
    }
    FrameData *getFrameData (int frame) const
    {
        return frames.at(frame);
    }

    bool getPixelShift (int frame = 0) const
    {
        // So far only Pentax provide multi-frame HDR file.
        // Only the first frame contains the HDR tag
        // If more brand have to be supported, this rule may need
        // to evolve

        //return frames.at(frame)->getPixelShift ();
        return frames.at(0)->getPixelShift ();
    }
    int getHDR (int frame = 0) const
    {
        // So far only Pentax provide multi-frame HDR file.
        // Only the first frame contains the HDR tag
        // If more brand have to be supported, this rule may need
        // to evolve

        //return frames.at(frame)->getPixelShift ();
        if (frames.size()) {
            return frames.at(frame)->getHDR ();
        } else {
            return 0;
        }
    }

    IIOSampleFormat getSampleFormat (int frame = 0) const
    {
        return frames.at(frame)->getSampleFormat ();
    }

    const rtexif::TagDirectory*   getExifData (int frame = 0) const
    {
        return frames.at(frame)->getExifData ();
    }
    const procparams::IPTCPairs   getIPTCData (int frame = 0) const
    {
        return frames.at(frame)->getIPTCData ();
    }

    bool hasExif (int frame = 0) const
    {
        return frames.at(frame)->hasExif ();
    }
    bool hasIPTC (int frame = 0) const
    {
        return frames.at(frame)->hasIPTC ();
    }

    struct tm   getDateTime (int frame = 0) const {
        return frames.at(frame)->getDateTime ();
    }
    time_t      getDateTimeAsTS(int frame = 0) const
    {
        return frames.at(frame)->getDateTimeAsTS ();
    }
    int         getISOSpeed (int frame = 0) const
    {
        return frames.at(frame)->getISOSpeed ();
    }
    double      getFNumber  (int frame = 0) const
    {
        return frames.at(frame)->getFNumber ();
    }
    double      getFocalLen (int frame = 0) const
    {
        return frames.at(frame)->getFocalLen ();
    }
    double      getFocalLen35mm (int frame = 0) const
    {
        return frames.at(frame)->getFocalLen35mm ();
    }
    float       getFocusDist (int frame = 0) const
    {
        return frames.at(frame)->getFocusDist ();
    }
    double      getShutterSpeed (int frame = 0) const
    {
        return frames.at(frame)->getShutterSpeed ();
    }
    double      getExpComp  (int frame = 0) const
    {
        return frames.at(frame)->getExpComp ();
    }
    std::string getMake     (int frame = 0) const
    {
        return frames.at(frame)->getMake ();
    }
    std::string getModel    (int frame = 0) const
    {
        return frames.at(frame)->getModel ();
    }
    std::string getLens     (int frame = 0) const
    {
        return frames.at(frame)->getLens ();
    }
    std::string getSerialNumber (int frame = 0) const
    {
        return frames.at(frame)->getSerialNumber ();
    }
    std::string getOrientation (int frame = 0) const
    {
        return frames.at(frame)->getOrientation ();
    }

};


}
#endif
