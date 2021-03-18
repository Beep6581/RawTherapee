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

#include <memory>

#include <glibmm/ustring.h>

#include <libiptcdata/iptc-data.h>

#include "iimage.h"
#include "imagedimensions.h"
#include "imageformat.h"
#include "rtengine.h"

enum {
    IMIO_SUCCESS,
    IMIO_CANNOTREADFILE,
    IMIO_INVALIDHEADER,
    IMIO_HEADERERROR,
    IMIO_READERROR,
    IMIO_VARIANTNOTSUPPORTED,
    IMIO_FILETYPENOTSUPPORTED,
    IMIO_CANNOTWRITEFILE
};

namespace rtexif
{

class TagDirectory;

}

namespace rtengine
{

class ColorTemp;
class ProgressListener;
class Imagefloat;

namespace procparams
{

class ExifPairs;

}

class ImageIO : virtual public ImageDatas
{

protected:
    ProgressListener* pl;
    cmsHPROFILE embProfile;
    char* profileData;
    int profileLength;
    char* loadedProfileData;
    int loadedProfileLength;
    const std::unique_ptr<procparams::ExifPairs> exifChange;
    IptcData* iptc;
    const rtexif::TagDirectory* exifRoot;
    MyMutex imutex;
    IIOSampleFormat sampleFormat;
    IIOSampleArrangement sampleArrangement;

private:
    void deleteLoadedProfileData( );

public:
    static Glib::ustring errorMsg[6];

    ImageIO();
    ~ImageIO() override;

    void setProgressListener (ProgressListener* l);
    void setSampleFormat(IIOSampleFormat sFormat);
    IIOSampleFormat getSampleFormat() const;
    void setSampleArrangement(IIOSampleArrangement sArrangement);
    IIOSampleArrangement getSampleArrangement() const;

    virtual void getStdImage (const ColorTemp &ctemp, int tran, Imagefloat* image, PreviewProps pp) const = 0;
    virtual int getBPS () const = 0;
    virtual void getScanline (int row, unsigned char* buffer, int bps, bool isFloat = false) const = 0;
    virtual void setScanline (int row, const unsigned char* buffer, int bps, unsigned int numSamples = 3) = 0;
    virtual const char* getType () const = 0;

    int load (const Glib::ustring &fname);
    int save (const Glib::ustring &fname) const;

    int loadPNG (const Glib::ustring &fname);
    int loadJPEG (const Glib::ustring &fname);
    int loadTIFF (const Glib::ustring &fname);
    static int getPNGSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement);
    static int getTIFFSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement);

    int loadJPEGFromMemory (const char* buffer, int bufsize);
    int loadPPMFromMemory(const char* buffer, int width, int height, bool swap, int bps);

    int savePNG (const Glib::ustring &fname, int bps = -1) const;
    int saveJPEG (const Glib::ustring &fname, int quality = 100, int subSamp = 3) const;
    int saveTIFF (const Glib::ustring &fname, int bps = -1, bool isFloat = false, bool uncompressed = false) const;

    cmsHPROFILE getEmbeddedProfile () const;
    void getEmbeddedProfileData (int& length, unsigned char*& pdata) const;

    void setMetadata (const rtexif::TagDirectory* eroot);
    void setMetadata (const rtexif::TagDirectory* eroot, const rtengine::procparams::ExifPairs& exif, const rtengine::procparams::IPTCPairs& iptcc);
    void setOutputProfile (const char* pdata, int plen);

    MyMutex& mutex ();
};

}
