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

class ImageData : public ImageMetaData
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

    void extractInfo ();

public:

    ImageData (Glib::ustring fname, RawMetaDataLocation* rml = NULL);
    virtual ~ImageData ();

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
}
#endif
