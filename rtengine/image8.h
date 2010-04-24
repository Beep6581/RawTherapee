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

#include <imageio.h>
#include <rtengine.h>

namespace rtengine {

class Image8 : public ImageIO, public IImage8 {

    public:
        unsigned char* data;
        int width;
        int height;
        
        Image8 ();
        Image8 (int width, int height);
       ~Image8 ();
     
        unsigned char r (int row, int col);
        unsigned char g (int row, int col);
        unsigned char b (int row, int col);
        void r (int row, int col, unsigned char val);
        void g (int row, int col, unsigned char val);
        void b (int row, int col, unsigned char val);
    
        virtual int     getW            () { return width;  }
        virtual int     getH            () { return height; }
        virtual void    allocate        (int width, int height);
        virtual int     getBPS          () { return 8; }
        virtual void    getScanline     (int row, unsigned char* buffer, int bps);
        virtual void    setScanline     (int row, unsigned char* buffer, int bps);
    
        // functions inherited from IImage*:
        virtual Glib::Mutex& getMutex () { return mutex (); }
        virtual cmsHPROFILE getProfile () { return getEmbeddedProfile (); }
        virtual int getWidth ()  { return width; }
        virtual int getHeight () { return height; }
        virtual int getBitsPerPixel () { return 16; }
        virtual int saveToFile (Glib::ustring fname) { return save (fname); }
        virtual int saveAsPNG  (Glib::ustring fname, int compression = -1, int bps = -1) { return savePNG (fname, compression, bps); }
        virtual int saveAsJPEG (Glib::ustring fname, int quality = 100) { return saveJPEG (fname, quality); }
        virtual int saveAsTIFF (Glib::ustring fname, int bps = -1) { return saveTIFF (fname, bps); }
        virtual void setSaveProgressListener (ProgressListener* pl) { setProgressListener (pl); } 
        virtual void free () { delete this; }
        virtual const unsigned char* getData () { return data; }

};
};
#endif
