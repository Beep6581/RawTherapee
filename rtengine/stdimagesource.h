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

namespace rtengine
{

class StdImageSource : public ImageSource
{

protected:
    ImageIO* img;
    ColorTemp wb;
    ProgressListener* plistener;
    bool full;
    float** hrmap[3];
    char** needhr;
    int max[3];
    bool rgbSourceModified;

    //void transformPixel             (int x, int y, int tran, int& tx, int& ty);
    void getSampleFormat (Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement);

public:
    StdImageSource ();
    ~StdImageSource ();

    int         load        (Glib::ustring fname, bool batch = false);
    void        getImage    (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, ToneCurveParams hrp, ColorManagementParams cmp, RAWParams raw);
    ColorTemp   getWB       ()
    {
        return wb;
    }
    void        getAutoWBMultipliers (double &rm, double &gm, double &bm);
    ColorTemp   getSpotWB   (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D> &blue, int tran, double equal);

    bool        isWBProviderReady ()
    {
        return true;
    };

    void        getAutoExpHistogram (LUTu &histogram, int& histcompr);

    double      getDefGain  ()
    {
        return 0.0;
    }

    void        getFullSize (int& w, int& h, int tr = TR_NONE);
    void        getSize     (int tran, PreviewProps pp, int& w, int& h);

    ImageData*  getImageData ()
    {
        return idata;
    }
    ImageIO*    getImageIO   ()
    {
        return img;
    }
    ImageMatrices* getImageMatrices ()
    {
        return (ImageMatrices*)NULL;
    }
    bool        isRAW() const
    {
        return false;
    }

    void        setProgressListener (ProgressListener* pl)
    {
        plistener = pl;
    }

    void        convertColorSpace(Imagefloat* image, ColorManagementParams cmp, ColorTemp &wb);// RAWParams raw will not be used for non-raw files (see imagesource.h)
    static void colorSpaceConversion (Imagefloat* im, ColorManagementParams cmp, cmsHPROFILE embedded, IIOSampleFormat sampleFormat);
    //static void colorSpaceConversion16 (Image16* im, ColorManagementParams cmp, cmsHPROFILE embedded);

    static inline double intpow (double a, int b)
    {
        double r = 1.0;

        for (int i = 0; i < b; i++) {
            r *= a;
        }

        return r;
    }
    bool        IsrgbSourceModified()
    {
        return rgbSourceModified;
    }
};
}
#endif
