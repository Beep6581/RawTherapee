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
// A class representing a 8 bit rgb image without alpha channel
//
#ifndef _IMAGE8_
#define _IMAGE8_

#include "imageio.h"
#include "rtengine.h"
#include "imagefloat.h"

namespace rtengine
{

class Image8 : public IImage8, public ImageIO
{

public:

    Image8 ();
    Image8 (int width, int height);
    ~Image8 ();

    Image8* copy () const;

    virtual void getStdImage (const ColorTemp &ctemp, int tran, Imagefloat* image, PreviewProps pp) const;

    virtual const char* getType () const
    {
        return sImage8;
    }
    virtual int getBPS () const
    {
        return 8 * sizeof(unsigned char);
    }
    virtual void getScanline (int row, unsigned char* buffer, int bps, bool isFloat = false) const;
    virtual void setScanline (int row, unsigned char* buffer, int bps, unsigned int numSamples);

    // functions inherited from IImage*:
    virtual MyMutex& getMutex ()
    {
        return mutex ();
    }

    virtual cmsHPROFILE getProfile () const
    {
        return getEmbeddedProfile ();
    }

    virtual int getBitsPerPixel () const
    {
        return 8 * sizeof(unsigned char);
    }

    virtual int saveToFile (const Glib::ustring &fname) const
    {
        return save (fname);
    }

    virtual int saveAsPNG (const Glib::ustring &fname, int bps = -1) const
    {
        return savePNG (fname, bps);
    }

    virtual int saveAsJPEG (const Glib::ustring &fname, int quality = 100, int subSamp = 3) const
    {
        return saveJPEG (fname, quality, subSamp);
    }

    virtual int saveAsTIFF (const Glib::ustring &fname, int bps = -1, bool isFloat = false, bool uncompressed = false) const
    {
        return saveTIFF (fname, bps, isFloat, uncompressed);
    }

    virtual void setSaveProgressListener (ProgressListener* pl)
    {
        setProgressListener (pl);
    }

    virtual void free ()
    {
        delete this;
    }

};

}
#endif
