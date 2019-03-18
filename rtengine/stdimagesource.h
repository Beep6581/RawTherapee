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
    int max[3];
    bool rgbSourceModified;

    //void transformPixel             (int x, int y, int tran, int& tx, int& ty);
    void getSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement);

public:
    StdImageSource ();
    ~StdImageSource () override;

    int         load        (const Glib::ustring &fname) override;
    void        getImage    (const ColorTemp &ctemp, int tran, Imagefloat* image, const PreviewProps &pp, const ToneCurveParams &hrp, const RAWParams &raw) override;
    ColorTemp   getWB       () const override
    {
        return wb;
    }
    void        getAutoWBMultipliers (double &rm, double &gm, double &bm) override;
    ColorTemp   getSpotWB   (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D> &blue, int tran, double equal) override;

    eSensorType getSensorType() const override {return ST_NONE;}
    bool isMono() const override {return false;}

    bool        isWBProviderReady () override
    {
        return true;
    };

    void        getAutoExpHistogram (LUTu &histogram, int& histcompr) override;

    double      getDefGain  () const override
    {
        return 0.0;
    }

    void        getFullSize (int& w, int& h, int tr = TR_NONE) override;
    void        getSize     (const PreviewProps &pp, int& w, int& h) override;

    ImageIO*    getImageIO   ()
    {
        return img;
    }
    ImageMatrices* getImageMatrices () override
    {
        return (ImageMatrices*)nullptr;
    }
    bool        isRAW() const override
    {
        return false;
    }

    void        setProgressListener (ProgressListener* pl) override
    {
        plistener = pl;
    }

    void        convertColorSpace(Imagefloat* image, const ColorManagementParams &cmp, const ColorTemp &wb) override;// RAWParams raw will not be used for non-raw files (see imagesource.h)
    static void colorSpaceConversion (Imagefloat* im, const ColorManagementParams &cmp, cmsHPROFILE embedded, IIOSampleFormat sampleFormat);

    bool        isRGBSourceModified() const override
    {
        return rgbSourceModified;
    }
    void setCurrentFrame(unsigned int frameNum) override {}
    int getFrameCount() override {return 1;}
    int getFlatFieldAutoClipValue() override {return 0;}


    void getRawValues(int x, int y, int rotate, int &R, int &G, int &B) override { R = G = B = 0;}

    void        flushRGB          () override;

};
}
#endif
