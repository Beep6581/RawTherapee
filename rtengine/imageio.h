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
#ifndef _IMAGEIO_
#define _IMAGEIO_

#define IMIO_SUCCESS               0
#define IMIO_CANNOTREADFILE        1
#define IMIO_INVALIDHEADER         2
#define IMIO_HEADERERROR           3
#define IMIO_READERROR             4
#define IMIO_VARIANTNOTSUPPORTED   5
#define IMIO_FILETYPENOTSUPPORTED  6
#define IMIO_CANNOTWRITEFILE       7

#include <memory>

#include <glibmm.h>
#include "rtengine.h"
#include "imageformat.h"
#include "imagedimensions.h"
#include "iimage.h"
#include "colortemp.h"
#include "procparams.h"

namespace rtengine
{

class ProgressListener;
class Imagefloat;

class MetadataInfo {
public:
    explicit MetadataInfo(const Glib::ustring &src=Glib::ustring()):
        src_(src) {}

    const Glib::ustring &filename() const { return src_; }

    const rtengine::procparams::ExifPairs &exif() const { return exif_; }
    const rtengine::procparams::IPTCPairs &iptc() const { return iptc_; }
    void setExif(const rtengine::procparams::ExifPairs &exif) { exif_ = exif; }
    void setIptc(const rtengine::procparams::IPTCPairs &iptc) { iptc_ = iptc; }

private:
    Glib::ustring src_;
    rtengine::procparams::ExifPairs exif_;
    rtengine::procparams::IPTCPairs iptc_;
};

class ImageIO : virtual public ImageDatas
{

protected:
    ProgressListener* pl;
    cmsHPROFILE embProfile;
    char* profileData;
    int profileLength;
    char* loadedProfileData;
    bool loadedProfileDataJpg;
    int loadedProfileLength;
    MyMutex imutex;
    IIOSampleFormat sampleFormat;
    IIOSampleArrangement sampleArrangement;
    MetadataInfo metadataInfo;

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
    virtual void setScanline (int row, unsigned char* buffer, int bps, unsigned int numSamples = 3) = 0;
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

    void setMetadata(const MetadataInfo &info) { metadataInfo = info; }
    void setOutputProfile (const char* pdata, int plen);

    bool saveMetadata(const Glib::ustring &fname) const;

    MyMutex& mutex ();
};

}
#endif
