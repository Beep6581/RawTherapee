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

#include <rtengine.h>
#include <procparams.h>
#include <libiptcdata/iptc-data.h>
#include <rtexif.h>
#include <png.h>

namespace rtengine {

class ImageIO {

    protected:
        ProgressListener* pl;
        cmsHPROFILE embProfile;
        char* profileData;
        int profileLength;
        char* loadedProfileData;
        int loadedProfileLength;
        std::vector<std::pair<std::string,std::string> > exifChange;
        IptcData* iptc;
        const rtexif::TagDirectory* exifRoot;

    public:
        static Glib::ustring errorMsg[6];
    	
        ImageIO () : pl (NULL), embProfile(NULL), profileData(NULL), exifRoot (NULL), iptc(NULL), loadedProfileData(NULL), loadedProfileLength(0) {}
        
        virtual ~ImageIO ();

        void setProgressListener (ProgressListener* l) { pl = l; }

        virtual int     getW            () =0;
		virtual int     getH            () =0;
		virtual void    allocate        (int width, int height) =0;
        virtual int     getBPS          () =0;
        virtual void    getScanline     (int row, unsigned char* buffer, int bps) {}
        virtual void    setScanline     (int row, unsigned char* buffer, int bps) {}

        int load (const String& fname);
        int save (const String& fname);

        int loadPNG  (const String& fname);
        int loadJPEG (const String& fname);
        int loadTIFF (const String& fname);

        int savePNG  (const String& fname, int compression = -1, int bps = -1);
        int saveJPEG (const String& fname, int quality = 100);
        int saveTIFF (const String& fname, int bps = -1, bool uncompressed = false);
        
        cmsHPROFILE getEmbeddedProfile () { return embProfile; }
        void        getEmbeddedProfileData (int& length, unsigned char*& pdata) { length = loadedProfileLength; pdata = (unsigned char*)loadedProfileData; }

        void setMetadata (const rtexif::TagDirectory* eroot, const std::vector<ExifPair>& exif, const std::vector<IPTCPair>& iptcc);
        void setOutputProfile  (char* pdata, int plen);


        static void png_read_data  (png_structp png_ptr, png_bytep data, png_size_t length);
        static void png_write_data (png_structp png_ptr, png_bytep data, png_size_t length);
        static void png_flush      (png_structp png_ptr);
};

};
#endif
