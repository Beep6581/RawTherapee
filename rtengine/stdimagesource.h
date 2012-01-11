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
#ifndef _STDIMAGESOURCE_
#define _STDIMAGESOURCE_

#include "imagesource.h"

namespace rtengine {

class StdImageSource : public ImageSource {

    protected:
        Image16* img;
        ColorTemp wb;
        ProgressListener* plistener;
        bool full;
        float** hrmap[3];
		char** needhr;
        int max[3];

        void transform           (PreviewProps pp, int tran, int &sx1, int &sy1, int &sx2, int &sy2);
        void transformPixel      (int x, int y, int tran, int& tx, int& ty);
        bool rgbSourceModified;
    public:
        StdImageSource ();
        ~StdImageSource ();
    
        int         load        (Glib::ustring fname, bool batch = false);
        void        getImage    (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, HRecParams hrp, ColorManagementParams cmp, RAWParams raw);
        ColorTemp   getWB       () { return wb; }
        ColorTemp   getAutoWB   ();
        ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue, int tran);
        bool        isWBProviderReady () { return true; };

        void        getAutoExpHistogram (LUTu &histogram, int& histcompr);
        
        double      getDefGain  () { return 0.0; }
        double      getGamma    () { return 0.0; }
        
        void        getFullSize (int& w, int& h, int tr = TR_NONE);
        void        getSize     (int tran, PreviewProps pp, int& w, int& h);

        ImageData*  getImageData () { return idata; }
        void        setProgressListener (ProgressListener* pl) { plistener = pl; }
        static void colorSpaceConversion (Imagefloat* im, ColorManagementParams cmp, cmsHPROFILE embedded);
		static void colorSpaceConversion16 (Image16* im, ColorManagementParams cmp, cmsHPROFILE embedded);


        static inline double intpow (double a, int b) { double r = 1.0; for (int i=0; i<b; i++) r *= a; return r; }
        bool        IsrgbSourceModified() {return rgbSourceModified;}
    protected:
        void    getImage_   (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, bool first, HRecParams hrp);
        void    hflip       (Imagefloat* im);
        void    vflip       (Imagefloat* im);
};
};
#endif
