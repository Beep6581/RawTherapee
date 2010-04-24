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

#include <stdio.h>
#include <common.h>
#include <string>
#include <glibmm.h>
#include <rtexif.h>
#include <procparams.h>
#include <libiptcdata/iptc-data.h>
#include <rtengine.h>

namespace rtengine {

class ImageData : public ImageMetaData {

  protected:
    rtexif::TagDirectory* root;
    IptcData* iptc;

    struct tm time;
    int iso_speed;
    double aperture;
    double focal_len;
    double shutter;
    std::string make, model;
    std::string lens;

    void extractInfo ();
    
  public:

    ImageData (Glib::ustring fname, RawMetaDataLocation* rml=NULL);
    ~ImageData ();

    const rtexif::TagDirectory*   getExifData () const { return root; }
    const std::vector<procparams::IPTCPair> getIPTCData () const;

    bool hasExif () const { return root && root->getCount(); }
    bool hasIPTC () const { return iptc; }

    struct tm   getDateTime () const { return time;      }
    int         getISOSpeed () const { return iso_speed; }
    double      getFNumber  () const { return aperture;  }
    double      getFocalLen () const { return focal_len;   }
    double      getShutterSpeed () const { return shutter;   }
    std::string getMake     () const { return make;      }
    std::string getModel    () const { return model;     }
    std::string getLens     () const { return lens;      }
};
};
#endif
