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
//
// A class representing a 16 bit rgb image with separate planes and 16 byte aligned data
//
#ifndef _IMAGEFLOAT_
#define _IMAGEFLOAT_

#include "imageio.h"
#include "rtengine.h"

namespace rtengine {
    using namespace procparams;

class Image8;
class Image16;

/*
 * Image type used by most tools; expected range: [0.0 ; 65535.0]
 */
class Imagefloat : public IImagefloat, public ImageIO {

    public:

        Imagefloat ();
        Imagefloat (int width, int height);
        ~Imagefloat ();

        Imagefloat*          copy ();

        Image8*              to8();
        Image16*             to16();

        virtual void         getStdImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, procparams::ToneCurveParams hrp);

        virtual const char*  getType     () const { return sImagefloat; }
        virtual int          getBPS      () { return 8*sizeof(float); }
        virtual void         getScanline (int row, unsigned char* buffer, int bps);
        virtual void         setScanline (int row, unsigned char* buffer, int bps, float *minValue=NULL, float *maxValue=NULL);

        // functions inherited from IImagefloat:
        virtual MyMutex&     getMutex () { return mutex (); }
        virtual cmsHPROFILE  getProfile () { return getEmbeddedProfile (); }
        virtual int          getBitsPerPixel () { return 8*sizeof(float); }
        virtual int          saveToFile (Glib::ustring fname) { return save (fname); }
        virtual int          saveAsPNG  (Glib::ustring fname, int compression = -1, int bps = -1) { return savePNG (fname, compression, bps); }
        virtual int          saveAsJPEG (Glib::ustring fname, int quality = 100, int subSamp = 3) { return saveJPEG (fname, quality, subSamp); }
        virtual int          saveAsTIFF (Glib::ustring fname, int bps = -1, bool uncompressed = false) { return saveTIFF (fname, bps, uncompressed); }
        virtual void         setSaveProgressListener (ProgressListener* pl) { setProgressListener (pl); }
        virtual void         free () { delete this; }

        virtual void         normalizeFloat(float srcMinVal, float srcMaxVal);
        void                 normalizeFloatTo1();
        void                 normalizeFloatTo65535();
        void                 calcCroppedHistogram(const ProcParams &params, float scale, LUTu & hist);

        void                 ExecCMSTransform(cmsHTRANSFORM hTransform);
};

}
#endif
