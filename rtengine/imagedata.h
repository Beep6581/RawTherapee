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

#include <cstdio>
#include <memory>
#include <string>
#include <vector>


#include <libiptcdata/iptc-data.h>

#include "imageio.h"

namespace Glib
{

class ustring;

}

namespace rtexif
{

class TagDirectory;
}

namespace rtengine
{

class FrameData final
{

protected:
    rtexif::TagDirectory* frameRootDir;
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
    int rating;
    std::string lens;
    IIOSampleFormat sampleFormat;

    // each frame has the knowledge of "being an"
    // or "being part of an" HDR or PS image
    bool isPixelShift;
    bool isHDR;

public:

    FrameData (rtexif::TagDirectory* frameRootDir, rtexif::TagDirectory* rootDir, rtexif::TagDirectory* firstRootDir);
    virtual ~FrameData ();

    bool getPixelShift () const;
    bool getHDR () const;
    std::string getImageType () const;
    IIOSampleFormat getSampleFormat () const;
    rtexif::TagDirectory* getExifData () const;
    procparams::IPTCPairs getIPTCData () const;
    static procparams::IPTCPairs getIPTCData (IptcData* iptc_);
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
    int getRating () const;
};

class FramesData final : public FramesMetaData {
private:
    // frame's root IFD, can be a file root IFD or a SUB-IFD
    std::vector<std::unique_ptr<FrameData>> frames;
    // root IFD in the file
    std::vector<rtexif::TagDirectory*> roots;
    IptcData* iptc;
    unsigned int dcrawFrameCount;

public:
    explicit FramesData (const Glib::ustring& fname, std::unique_ptr<RawMetaDataLocation> rml = nullptr, bool firstFrameOnly = false);
    ~FramesData () override;

    void setDCRawFrameCount (unsigned int frameCount);
    unsigned int getRootCount () const override;
    unsigned int getFrameCount () const override;
    bool getPixelShift () const override;
    bool getHDR (unsigned int frame = 0) const override;
    std::string getImageType (unsigned int frame) const override;
    IIOSampleFormat getSampleFormat (unsigned int frame = 0) const override;
    rtexif::TagDirectory* getFrameExifData (unsigned int frame = 0) const override;
    rtexif::TagDirectory* getRootExifData (unsigned int root = 0) const override;
    rtexif::TagDirectory* getBestExifData (ImageSource *imgSource, procparams::RAWParams *rawParams) const override;
    procparams::IPTCPairs getIPTCData (unsigned int frame = 0) const override;
    bool hasExif (unsigned int frame = 0) const override;
    bool hasIPTC (unsigned int frame = 0) const override;
    tm getDateTime (unsigned int frame = 0) const override;
    time_t getDateTimeAsTS (unsigned int frame = 0) const override;
    int getISOSpeed (unsigned int frame = 0) const override;
    double getFNumber (unsigned int frame = 0) const override;
    double getFocalLen (unsigned int frame = 0) const override;
    double getFocalLen35mm (unsigned int frame = 0) const override;
    float getFocusDist (unsigned int frame = 0) const override;
    double getShutterSpeed (unsigned int frame = 0) const override;
    double getExpComp (unsigned int frame = 0) const override;
    std::string getMake (unsigned int frame = 0) const override;
    std::string getModel (unsigned int frame = 0) const override;
    std::string getLens (unsigned int frame = 0) const override;
    std::string getSerialNumber (unsigned int frame = 0) const;
    std::string getOrientation (unsigned int frame = 0) const override;
    int getRating (unsigned int frame = 0) const override;
};


}
