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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <array>
#include <vector>

#include <glibmm.h>

#include "colortemp.h"
#include "coord2d.h"
#include "dcp.h"
#include "image16.h"
#include "image8.h"
#include "imagedata.h"
#include "imagefloat.h"
#include "LUT.h"
#include "rtengine.h"

namespace rtengine
{

namespace procparams
{

struct CoarseTransformParams;
struct ColorManagementParams;
struct FilmNegativeParams;
struct LensProfParams;
struct RAWParams;
struct RetinexParams;
struct ToneCurveParams;
struct CaptureSharpeningParams;
}

class ImageMatrices
{

public:
    double rgb_cam[3][3] = {};
    double cam_rgb[3][3] = {};
    double xyz_cam[3][3] = {};
    double cam_xyz[3][3] = {};
};

// TODO: Move implementation to .cc (Flössie)
class ImageSource : public InitialImage
{

private:
    int references;

protected:
    double redAWBMul, greenAWBMul, blueAWBMul; // local copy of the multipliers, to avoid recomputing the values
    cmsHPROFILE embProfile;
    Glib::ustring fileName;
    FramesData* idata;
    ImageMatrices imatrices;
    double dirpyrdenoiseExpComp;

public:
    ImageSource () : references (1), redAWBMul(-1.), greenAWBMul(-1.), blueAWBMul(-1.),
        embProfile(nullptr), idata(nullptr), dirpyrdenoiseExpComp(INFINITY) {}

    ~ImageSource            () override {}
    virtual int         load        (const Glib::ustring &fname) = 0;
    virtual void        preprocess  (const procparams::RAWParams &raw, const procparams::LensProfParams &lensProf, const procparams::CoarseTransformParams& coarse, bool prepareDenoise = true) {};
    virtual void        filmNegativeProcess (const procparams::FilmNegativeParams &params) {};
    virtual bool        getFilmNegativeExponents (Coord2D spotA, Coord2D spotB, int tran, const FilmNegativeParams& currentParams, std::array<float, 3>& newExps) { return false; };
    virtual bool        getFilmNegativeMedians (Coord2D topLeft, Coord2D bottomRight, int tran, const FilmNegativeParams &currentParams, std::array<float, 3>& medians) { return false; };
    virtual void        demosaic    (const procparams::RAWParams &raw, bool autoContrast, double &contrastThreshold, bool cache = false) {};
    virtual void        retinex       (const procparams::ColorManagementParams& cmp, const procparams::RetinexParams &deh, const procparams::ToneCurveParams& Tc, LUTf & cdcurve, LUTf & mapcurve, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, multi_array2D<float, 4> &conversionBuffer, bool dehacontlutili, bool mapcontlutili, bool useHsl, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax, LUTu &histLRETI) {};
    virtual void        retinexPrepareCurves       (const procparams::RetinexParams &retinexParams, LUTf &cdcurve, LUTf &mapcurve, RetinextransmissionCurve &retinextransmissionCurve, RetinexgaintransmissionCurve &retinexgaintransmissionCurve, bool &retinexcontlutili, bool &mapcontlutili, bool &useHsl, LUTu & lhist16RETI, LUTu & histLRETI) {};
    virtual void        retinexPrepareBuffers      (const procparams::ColorManagementParams& cmp, const procparams::RetinexParams &retinexParams, multi_array2D<float, 4> &conversionBuffer, LUTu &lhist16RETI) {};
    virtual void        flushRawData       () {};
    virtual void        flushRGB           () {};
    virtual void        HLRecovery_Global  (const procparams::ToneCurveParams &hrp) {};
    virtual void        HLRecovery_inpaint (float** red, float** green, float** blue) {};

    virtual bool        isRGBSourceModified () const = 0; // tracks whether cached rgb output of demosaic has been modified

    virtual void        setBorder (unsigned int border) {}
    virtual void        setCurrentFrame (unsigned int frameNum) = 0;
    virtual int         getFrameCount () = 0;
    virtual int         getFlatFieldAutoClipValue () = 0;


    // use right after demosaicing image, add coarse transformation and put the result in the provided Imagefloat*
    virtual void        getImage    (const ColorTemp &ctemp, int tran, Imagefloat* image, const PreviewProps &pp, const procparams::ToneCurveParams &hlp, const RAWParams &raw) = 0;
    virtual eSensorType getSensorType () const = 0;
    virtual bool        isMono () const = 0;
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
    virtual void        getSize     (const PreviewProps &pp, int& w, int& h) = 0;
    virtual int         getRotateDegree() const
    {
        return 0;
    }

    virtual ImageMatrices* getImageMatrices () = 0;
    virtual bool           isRAW () const = 0;
    virtual DCPProfile*    getDCP (const procparams::ColorManagementParams &cmp, DCPProfile::ApplyState &as)
    {
        return nullptr;
    };

    virtual void        setProgressListener (ProgressListener* pl) {}

    void        increaseRef () override
    {
        references++;
    }
    void        decreaseRef () override
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

    // for RAW files, compute a tone curve using histogram matching on the embedded thumbnail
    virtual void getAutoMatchedToneCurve(const procparams::ColorManagementParams &cp, std::vector<double> &outCurve)
    {
        outCurve = { 0.0 };
    }
    
    double getDirPyrDenoiseExpComp ( )
    {
        return dirpyrdenoiseExpComp;
    }
    // functions inherited from the InitialImage interface
    Glib::ustring getFileName () override
    {
        return fileName;
    }
    cmsHPROFILE getEmbeddedProfile () override
    {
        return embProfile;
    }
    const FramesMetaData* getMetaData () override
    {
        return idata;
    }
    ImageSource* getImageSource () override
    {
        return this;
    }
    virtual void getRawValues(int x, int y, int rotate, int &R, int &G, int &B) = 0;
    virtual void captureSharpening(const procparams::CaptureSharpeningParams &sharpeningParams, bool showMask, double &conrastThreshold, double &radius) = 0;
};

}
