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
#ifndef _IIMAGE_
#define _IIMAGE_

#include <lcms.h>
#include <glibmm.h>

namespace rtengine {

    class ProgressListener;
    
  /** This class represents an image (the result of the image processing) */
    class IImage {
        public:
        /** Returns a mutex that can is useful in many situations. No image operations shuold be performed without locking this mutex.
          * @return The mutex */
            virtual Glib::Mutex& getMutex ();
            virtual cmsHPROFILE getProfile ();
        /** Returns the width of the image.
          * @return The width of the image */
            virtual int getWidth ();
        /** Returns the height of the image.
          * @return The height of the image */
            virtual int getHeight ();
        /** Returns the bits per pixel of the image.
          * @return The bits per pixel of the image */
            virtual int getBitsPerPixel ();
        /** Saves the image to file. It autodetects the format (jpg, tif, png are supported).
          * @param fname is the name of the file 
            @return the error code, 0 if none */
            virtual int saveToFile (Glib::ustring fname);
        /** Saves the image to file in a png format.
          * @param fname is the name of the file 
          * @param compression is the amount of compression (0-6), -1 corresponds to the default
          * @param bps can be 8 or 16 depending on the bits per pixels the output file will have
            @return the error code, 0 if none */
            virtual int saveAsPNG  (Glib::ustring fname, int compression = -1, int bps = -1);
        /** Saves the image to file in a jpg format.
          * @param fname is the name of the file 
          * @param quality is the quality of the jpeg (0...100), set it to -1 to use default
            @return the error code, 0 if none */
            virtual int saveAsJPEG (Glib::ustring fname, int quality = 100);
        /** Saves the image to file in a tif format.
          * @param fname is the name of the file 
          * @param bps can be 8 or 16 depending on the bits per pixels the output file will have
            @return the error code, 0 if none */
            virtual int saveAsTIFF (Glib::ustring fname, int bps = -1);
        /** Sets the progress listener if you want to follow the progress of the image saving operations (optional).
          * @param pl is the pointer to the class implementing the ProgressListener interface */
            virtual void setSaveProgressListener (ProgressListener* pl);
        /** Free the image */
            virtual void free ();
    };

  /** This class represents an image having a classical 8 bits/pixel representation */
    class IImage8 : public IImage {
        public:
        /** Returns the pixel data, in r/g/b order from top left to bottom right continously.
          * @return a pointer to the pixel data */
            virtual const unsigned char* getData ();
    };

  /** This class represents an image having a 16 bits/pixel planar representation. 
      The planes are stored as two dimensional arrays. All the rows have a 8-byte alignment. */
    class IImage16 : public IImage {
        public:
        /** Returns the "red" plane data.
          * @return the two dimensional array of the red plane */
            virtual unsigned short** getRPlane ();
        /** Returns the "green" plane data.
          * @return the two dimensional array of the green plane */
            virtual unsigned short** getGPlane ();
        /** Returns the "blue" plane data.
          * @return the two dimensional array of the blue plane */
            virtual unsigned short** getBPlane ();
    };
}

#endif
