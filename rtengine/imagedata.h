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

#include "dnggainmap.h"
#include "imageio.h"
#include "metadata.h"

namespace Glib
{

class ustring;

}

namespace rtengine
{

class FramesData final :
    public FramesMetaData
{
private:
    bool ok_;
    Glib::ustring fname_;
    unsigned int dcrawFrameCount;
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
    struct tm modTime;
    time_t modTimeStamp;
    bool isPixelShift;
    bool isHDR;
    bool isDNG;
    std::uint32_t fixBadPixelsConstant;
    bool hasFixBadPixelsConstant_{false};
    std::vector<GainMap> gain_maps_;
    int w_;
    int h_;

public:
    explicit FramesData(const Glib::ustring& fname, time_t ts = 0);

    void setDCRawFrameCount(unsigned int frameCount);
    unsigned int getFrameCount() const override;
    bool getPixelShift() const override;
    bool getHDR() const override;
    bool getDNG() const override;
    std::string getImageType() const override;
    IIOSampleFormat getSampleFormat() const override;
    bool hasExif() const override;
    tm getDateTime() const override;
    time_t getDateTimeAsTS() const override;
    int getISOSpeed() const override;
    double getFNumber() const override;
    double getFocalLen() const override;
    double getFocalLen35mm() const override;
    float getFocusDist() const override;
    double getShutterSpeed() const override;
    double getExpComp() const override;
    std::string getMake() const override;
    std::string getModel() const override;
    std::string getLens() const override;
    std::string getSerialNumber() const;
    std::string getOrientation() const override;
    Glib::ustring getFileName() const override;
    int getRating() const override;
    std::uint32_t getFixBadPixelsConstant() const override;
    bool hasFixBadPixelsConstant() const override;
    std::vector<GainMap> getGainMaps() const override;
    void getDimensions(int &w, int &h) const override;

    void fillBasicTags(Exiv2::ExifData &exif) const;

    void setDimensions(int w, int h);
};

}
