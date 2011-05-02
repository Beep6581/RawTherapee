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

#include <imageio.h>
#include <rtengine.h>

namespace rtengine {
    using namespace procparams;
//enum TypeInterpolation { TI_Nearest, TI_Bilinear };

class Image8;
class Image16;

class Imagefloat : public ImageIO, public IImagefloat {

    private:
        unsigned char* unaligned;

    public:
        int rowstride;
        int planestride;

        int width;
        int height;

        float* data;
        
        float** r;
        float** g;
        float** b;
  
  
        Imagefloat ();
        Imagefloat (int width, int height);
        ~Imagefloat ();

        Imagefloat* copy ();
	
        Image8* to8() const;
        Image16* to16() const;


        Imagefloat* rotate (int deg);
        Imagefloat* hflip ();
        Imagefloat* vflip ();
        //Imagefloat* resize (int nw, int nh, TypeInterpolation interp);

        virtual int     getW            () { return width;  }
        virtual int     getH            () { return height; }
        virtual void    allocate        (int width, int height);
        virtual int     getBPS          () { return 8*sizeof(float); }
        //virtual void    getScanline     (int row, unsigned char* buffer, int bps);
        //virtual void    setScanline     (int row, unsigned char* buffer, int bps);

        // functions inherited from IImagefloat:
        virtual int getWidth ()  { return width; }
        virtual int getHeight () { return height; }
        virtual Glib::Mutex& getMutex () { return mutex (); }
        virtual cmsHPROFILE getProfile () { return getEmbeddedProfile (); }
        virtual int getBitsPerPixel () { return 16; }
        virtual int saveToFile (Glib::ustring fname) { return save (fname); }
        virtual int saveAsPNG  (Glib::ustring fname, int compression = -1, int bps = -1) { return savePNG (fname, compression, bps); }
        virtual int saveAsJPEG (Glib::ustring fname, int quality = 100) { return saveJPEG (fname, quality); }
        virtual int saveAsTIFF (Glib::ustring fname, int bps = -1, bool uncompressed = false) { return saveTIFF (fname, bps, uncompressed); }
        virtual void setSaveProgressListener (ProgressListener* pl) { return setProgressListener (pl); } 
        virtual void free () { delete this; }
        virtual float** getRPlane () { return r; }
        virtual float** getGPlane () { return g; }
        virtual float** getBPlane () { return b; }

        void calcCroppedHistogram(const ProcParams &params, float scale, LUTu & hist);

        void ExecCMSTransform(cmsHTRANSFORM hTransform, bool safe);
    };
};
#endif
