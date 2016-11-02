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
#ifndef _IMAGESOURCE_
#define _IMAGESOURCE_

#include <glibmm.h>
#include <vector>
#include "rtengine.h"
#include "colortemp.h"
#include "procparams.h"
#include "coord2d.h"
#include "dcp.h"
#include "LUT.h"
#include "imagedata.h"
#include "image8.h"
#include "image16.h"
#include "imagefloat.h"

namespace rtengine
{

using namespace procparams;

class ImageMatrices
{

public:
    double rgb_cam[3][3] = {};
    double cam_rgb[3][3] = {};
    double xyz_cam[3][3] = {};
    double cam_xyz[3][3] = {};
};

class ImageSource : public InitialImage
{

private:
    int references;

protected:
    double redAWBMul, greenAWBMul, blueAWBMul; // local copy of the multipliers, to avoid recomputing the values
    cmsHPROFILE embProfile;
    Glib::ustring fileName;
    ImageData* idata;
    ImageMatrices imatrices;
    double dirpyrdenoiseExpComp;

public:
    ImageSource () : references (1), redAWBMul(-1.), greenAWBMul(-1.), blueAWBMul(-1.),
        embProfile(nullptr), idata(nullptr), dirpyrdenoiseExpComp(INFINITY) {}

    virtual ~ImageSource            () {}
    virtual int         load        (const Glib::ustring &fname, int imageNum = 0, bool batch = false) = 0;
    virtual void        preprocess  (const RAWParams &raw, const LensProfParams &lensProf, const CoarseTransformParams& coarse, bool prepareDenoise = true) {};
    virtual void        demosaic    (const RAWParams &raw) {};
    virtual void        retinex       (ColorManagementParams cmp, RetinexParams  deh, ToneCurveParams Tc, LUTf & cdcurve, LUTf & mapcurve, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, multi_array2D<float, 4> &conversionBuffer, bool dehacontlutili, bool mapcontlutili, bool useHsl, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax, LUTu &histLRETI) {};
    virtual void        retinexPrepareCurves       (RetinexParams retinexParams, LUTf &cdcurve, LUTf &mapcurve, RetinextransmissionCurve &retinextransmissionCurve, RetinexgaintransmissionCurve &retinexgaintransmissionCurve, bool &retinexcontlutili, bool &mapcontlutili, bool &useHsl, LUTu & lhist16RETI, LUTu & histLRETI) {};
    virtual void        retinexPrepareBuffers      (ColorManagementParams cmp, RetinexParams retinexParams, multi_array2D<float, 4> &conversionBuffer, LUTu &lhist16RETI) {};
    virtual void        flushRawData       () {};
    virtual void        flushRGB           () {};
    virtual void        HLRecovery_Global  (ToneCurveParams hrp) {};
    virtual void        HLRecovery_inpaint (float** red, float** green, float** blue) {};
    virtual void        MSR(LabImage* lab, LUTf & mapcurve, bool &mapcontlutili, int width, int height, int skip, RetinexParams deh, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax) {};

    virtual bool        IsrgbSourceModified() const = 0; // tracks whether cached rgb output of demosaic has been modified

    virtual void setCurrentFrame(unsigned int frameNum) = 0;


    // use right after demosaicing image, add coarse transformation and put the result in the provided Imagefloat*
    virtual void        getImage    (const ColorTemp &ctemp, int tran, Imagefloat* image, const PreviewProps &pp, const ToneCurveParams &hlp, const ColorManagementParams &cmp, const RAWParams &raw) = 0;
    virtual eSensorType getSensorType ()
    {
        return ST_NONE;
    }
    // true is ready to provide the AutoWB, i.e. when the image has been demosaiced for RawImageSource
    virtual bool        isWBProviderReady () = 0;

    virtual void        convertColorSpace    (Imagefloat* image, const ColorManagementParams &cmp, const ColorTemp &wb) = 0; // DIRTY HACK: this method is derived in rawimagesource and strimagesource, but (...,RAWParams raw) will be used ONLY for raw images
    virtual void        getAutoWBMultipliers (double &rm, double &gm, double &bm) = 0;
    virtual ColorTemp   getWB       () const = 0;
    virtual ColorTemp   getSpotWB   (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D> &blue, int tran, double equal) = 0;

    virtual double      getDefGain  () const
    {
        return 1.0;
    }

    virtual void        getFullSize (int& w, int& h, int tr = TR_NONE) {}
    virtual void        getSize     (PreviewProps pp, int& w, int& h) = 0;
    virtual int         getRotateDegree() const
    {
        return 0;
    }

    virtual ImageData*     getImageData () = 0;
    virtual ImageMatrices* getImageMatrices () = 0;
    virtual bool        isRAW() const = 0;
    virtual DCPProfile* getDCP(const ColorManagementParams &cmp, ColorTemp &wb, DCPProfile::ApplyState &as)
    {
        return nullptr;
    };

    virtual void        setProgressListener (ProgressListener* pl) {}

    void        increaseRef ()
    {
        references++;
    }
    void        decreaseRef ()
    {
        references--;

        if (!references) {
            delete this;
        }
    }

    virtual void        getAutoExpHistogram (LUTu & histogram, int& histcompr) = 0;
    virtual void        getRAWHistogram (LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw)
    {
        histRedRaw.clear();
        histGreenRaw.clear();
        histBlueRaw.clear(); // only some sources will supply this
    }
    double getDirPyrDenoiseExpComp ( )
    {
        return dirpyrdenoiseExpComp;
    }
    // functions inherited from the InitialImage interface
    virtual Glib::ustring getFileName ()
    {
        return fileName;
    }
    virtual cmsHPROFILE getEmbeddedProfile ()
    {
        return embProfile;
    }
    virtual const ImageMetaData* getMetaData ()
    {
        return idata;
    }
    virtual ImageSource* getImageSource ()
    {
        return this;
    }
};
}
#endif
