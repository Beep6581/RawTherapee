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

#include "rtengine.h"
#include <glibmm.h>
#include "procparams.h"
#include <libiptcdata/iptc-data.h>
#include "../rtexif/rtexif.h"

namespace rtengine {

class ImageIO {

    protected:
        ProgressListener* pl;
        cmsHPROFILE embProfile;
        char* profileData;
        int profileLength;
        char* loadedProfileData;
        int loadedProfileLength;
        procparams::ExifPairs exifChange;
        IptcData* iptc;
        const rtexif::TagDirectory* exifRoot;
        Glib::Mutex imutex;

    public:
        static Glib::ustring errorMsg[6];

        ImageIO () : pl (NULL), embProfile(NULL), profileData(NULL), loadedProfileData(NULL), loadedProfileLength(0), iptc(NULL), exifRoot (NULL) {}

        virtual ~ImageIO ();

        void setProgressListener (ProgressListener* l) { pl = l; }

        virtual int     getW            () =0;
        virtual int     getH            () =0;
        virtual void    allocate        (int width, int height) =0;
        virtual int     getBPS          () =0;
        virtual void    getScanline     (int row, unsigned char* buffer, int bps) {}
        virtual void    setScanline     (int row, unsigned char* buffer, int bps) {}

        int load (Glib::ustring fname);
        int save (Glib::ustring fname);

        int loadPNG  (Glib::ustring fname);
        int loadJPEG (Glib::ustring fname);
        int loadTIFF (Glib::ustring fname);

        int loadJPEGFromMemory (const char* buffer, int bufsize);
        int loadPPMFromMemory(const char* buffer,int width,int height, bool swap, int bps);

        int savePNG  (Glib::ustring fname, int compression = -1, volatile int bps = -1);
        int saveJPEG (Glib::ustring fname, int quality = 100);
        int saveTIFF (Glib::ustring fname, int bps = -1, bool uncompressed = false);

        cmsHPROFILE getEmbeddedProfile () { return embProfile; }
        void        getEmbeddedProfileData (int& length, unsigned char*& pdata) { length = loadedProfileLength; pdata = (unsigned char*)loadedProfileData; }

        void setMetadata (const rtexif::TagDirectory* eroot);
        void setMetadata (const rtexif::TagDirectory* eroot, const rtengine::procparams::ExifPairs& exif, const rtengine::procparams::IPTCPairs& iptcc);
        void setOutputProfile  (char* pdata, int plen);
        Glib::Mutex& mutex () { return imutex; }
};

};
#endif
