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

#include "rtengine.h"
#include <glibmm.h>
#include "procparams.h"
#include <libiptcdata/iptc-data.h>
#include "../rtexif/rtexif.h"
#include "imagedimensions.h"
#include "iimage.h"
#include "../rtgui/threadutils.h"

namespace rtengine {

    class ProgressListener;
    class Imagefloat;

    typedef enum IIO_Sample_Format {
        IIOSF_UNKNOWN        = 0,       // Unknown or Unsupported file type; Has to remain 0
        //IIOSF_SIGNED_INT         ,    // Not yet supported
        IIOSF_UNSIGNED_CHAR  = 1<<0,
        IIOSF_UNSIGNED_SHORT = 1<<1,
        //IIOSF_HALF               ,    // OpenEXR & NVidia's Half Float, not yet supported
        IIOSF_LOGLUV24       = 1<<2,
        IIOSF_LOGLUV32       = 1<<3,
        IIOSF_FLOAT          = 1<<4
    } IIOSampleFormat;

    typedef enum IIO_Sample_Arrangement {
        IIOSA_UNKNOWN,       // Unknown or Unsupported file type
        IIOSA_CHUNKY,
        IIOSA_PLANAR
    } IIOSampleArrangement;

    typedef enum SensorType {
        ST_NONE,   // use this value if the image is already demosaiced (i.e. not a raw file)
        ST_BAYER,
        ST_FUJI_XTRANS,
        //ST_FOVEON,
        //ST_FUJI_EXR
    } eSensorType;

class ImageIO : virtual public ImageDatas {

    protected:
        ProgressListener* pl;
        cmsHPROFILE embProfile;
        char* profileData;
        int profileLength;
        char* loadedProfileData;
        bool loadedProfileDataJpg;
        int loadedProfileLength;
        procparams::ExifPairs exifChange;
        IptcData* iptc;
        const rtexif::TagDirectory* exifRoot;
        MyMutex imutex;
        IIOSampleFormat sampleFormat;
        IIOSampleArrangement sampleArrangement;

	private:
		void deleteLoadedProfileData( ) { if(loadedProfileData) {if(loadedProfileDataJpg) free(loadedProfileData); else delete[] loadedProfileData;} loadedProfileData = NULL; }
    public:
        static Glib::ustring errorMsg[6];

        ImageIO () : pl (NULL), embProfile(NULL), profileData(NULL), profileLength(0), loadedProfileData(NULL),loadedProfileDataJpg(false),
                     loadedProfileLength(0), iptc(NULL), exifRoot (NULL), sampleFormat(IIOSF_UNKNOWN),
                     sampleArrangement(IIOSA_UNKNOWN) {}

        virtual ~ImageIO ();

        void setProgressListener (ProgressListener* l) { pl = l; }

        void                 setSampleFormat(IIOSampleFormat sFormat) { sampleFormat = sFormat; }
        IIOSampleFormat      getSampleFormat() { return sampleFormat; }
        void                 setSampleArrangement(IIOSampleArrangement sArrangement) { sampleArrangement = sArrangement; }
        IIOSampleArrangement getSampleArrangement() { return sampleArrangement; }

        virtual void    getStdImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, procparams::ToneCurveParams hrp) {
            printf("getStdImage NULL!\n");
        }

        virtual int     getBPS      () =0;
        virtual void    getScanline (int row, unsigned char* buffer, int bps) {}
        virtual void    setScanline (int row, unsigned char* buffer, int bps, float minValue[3]=NULL, float  maxValue[3]=NULL) {}

        virtual bool    readImage   (Glib::ustring &fname, FILE *fh) { return false; };
        virtual bool    writeImage  (Glib::ustring &fname, FILE *fh) { return false; };

        int load (Glib::ustring fname);
        int save (Glib::ustring fname);

        int loadPNG  (Glib::ustring fname);
        int loadJPEG (Glib::ustring fname);
        int loadTIFF (Glib::ustring fname);
        static int getPNGSampleFormat  (Glib::ustring fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement);
        static int getTIFFSampleFormat (Glib::ustring fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement);

        int loadJPEGFromMemory (const char* buffer, int bufsize);
        int loadPPMFromMemory(const char* buffer,int width,int height, bool swap, int bps);

        int savePNG  (Glib::ustring fname, int compression = -1, volatile int bps = -1);
        int saveJPEG (Glib::ustring fname, int quality = 100, int subSamp=3);
        int saveTIFF (Glib::ustring fname, int bps = -1, bool uncompressed = false);

        cmsHPROFILE getEmbeddedProfile () { return embProfile; }
        void        getEmbeddedProfileData (int& length, unsigned char*& pdata) { length = loadedProfileLength; pdata = (unsigned char*)loadedProfileData; }

        void setMetadata (const rtexif::TagDirectory* eroot);
        void setMetadata (const rtexif::TagDirectory* eroot, const rtengine::procparams::ExifPairs& exif, const rtengine::procparams::IPTCPairs& iptcc);
        void setOutputProfile  (char* pdata, int plen);
        MyMutex& mutex () { return imutex; }
};

}
#endif
