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

#include <glibmm/ustring.h>

#include "coord2d.h"
#include "imagedata.h"
#include "rtengine.h"
#include "colortemp.h"
#include "array2D.h"

template<typename T>
class LUT;

using LUTf = LUT<float>;

template<typename T, const size_t num>
class multi_array2D;

namespace rtengine
{

class ColorTemp;
class DCPProfile;
class DCPProfileApplyState;
class Imagefloat;
class RetinexgaintransmissionCurve;
class RetinextransmissionCurve;

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
};

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
        embProfile(nullptr), idata(nullptr), dirpyrdenoiseExpComp(RT_INFINITY) {}

    ~ImageSource            () override {}
    virtual int         load        (const Glib::ustring &fname) = 0;
    virtual void        preprocess  (const procparams::RAWParams &raw, const procparams::LensProfParams &lensProf, const procparams::CoarseTransformParams& coarse, float &reddeha, float &greendeha, float &bluedeha, bool prepareDenoise = true) {};
    virtual void        demosaic    (const procparams::RAWParams &raw, bool autoContrast, double &contrastThreshold, bool cache = false) {};
    virtual void        retinex       (const procparams::ColorManagementParams& cmp, const procparams::RetinexParams &deh, const procparams::ToneCurveParams& Tc, LUTf & cdcurve, LUTf & mapcurve, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, multi_array2D<float, 4> &conversionBuffer, bool dehacontlutili, bool mapcontlutili, bool useHsl, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax, LUTu &histLRETI) {};
    virtual void        retinexPrepareCurves       (const procparams::RetinexParams &retinexParams, LUTf &cdcurve, LUTf &mapcurve, RetinextransmissionCurve &retinextransmissionCurve, RetinexgaintransmissionCurve &retinexgaintransmissionCurve, bool &retinexcontlutili, bool &mapcontlutili, bool &useHsl, LUTu & lhist16RETI, LUTu & histLRETI) {};
    virtual void        retinexPrepareBuffers      (const procparams::ColorManagementParams& cmp, const procparams::RetinexParams &retinexParams, multi_array2D<float, 4> &conversionBuffer, LUTu &lhist16RETI) {};
    virtual void        flush           () = 0;
    virtual void        HLRecovery_Global  (const procparams::ToneCurveParams &hrp) {};

    virtual bool        isRGBSourceModified () const = 0; // tracks whether cached rgb output of demosaic has been modified

    virtual void        setBorder (unsigned int border) {}
    virtual void        setCurrentFrame (unsigned int frameNum) = 0;
    virtual int         getFrameCount () = 0;
    virtual int         getFlatFieldAutoClipValue () = 0;

    virtual void        getWBMults  (const ColorTemp &ctemp, const procparams::RAWParams &raw, std::array<float, 4>& scale_mul, float &autoGainComp, float &rm, float &gm, float &bm) const = 0;

    // use right after demosaicing image, add coarse transformation and put the result in the provided Imagefloat*
    virtual void        getImage    (const ColorTemp &ctemp, int tran, Imagefloat* image, const PreviewProps &pp, const procparams::ToneCurveParams &hlp, const procparams::RAWParams &raw) = 0;
    virtual eSensorType getSensorType () const = 0;
    virtual bool        isMono () const = 0;
    // true is ready to provide the AutoWB, i.e. when the image has been demosaiced for RawImageSource
    virtual bool        isWBProviderReady () = 0;

    virtual void        convertColorSpace    (Imagefloat* image, const procparams::ColorManagementParams &cmp, const ColorTemp &wb) = 0; // DIRTY HACK: this method is derived in rawimagesource and strimagesource, but (...,RAWParams raw) will be used ONLY for raw images
    virtual void        getAutoWBMultipliers (double &rm, double &gm, double &bm) = 0;
    virtual void        getAutoWBMultipliersitc(bool extra, double &tempref, double &greenref, double &tempitc, double & greenitc, float &temp0, float &delta, int &bia, int &dread, int &kcam, int &nocam, float &studgood, float &minchrom, int &kmin, float &minhist, float &maxhist, int begx, int begy, int yEn, int xEn, int cx, int cy, int bf_h, int bf_w, double &rm, double &gm, double &bm, const procparams::WBParams & wbpar, const procparams::ColorManagementParams &cmp, const procparams::RAWParams &raw, const procparams::ToneCurveParams &hrp) = 0;
    virtual void getAutoWBMultipliersItcGreen(
        procparams::ProcParams &params,
        bool forcewbgrey,
        int kcam,
        double greenitc,
        bool extra,
        float &temp0,
        float &delta,
        int &bia,
        int &dread,
        int nocam,
        float &studgood,
        float &minchrom,
        int &kmin,
        float &minhist,
        float &maxhist,
        int fh,
        int fw,
        ColorTemp &currWB,
        int tempnotisraw,
        double greennotisraw,
        bool skipRecalculate,
        ColorTemp &autoWB,
        double &rm,
        double &gm,
        double &bm
        );
    virtual ColorTemp   getWB       () const = 0;
    virtual ColorTemp   getSpotWB   (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D> &blue, int tran, double equal, StandardObserver observer) = 0;
    virtual void        WBauto(bool extra, double &tempref, double &greenref, array2D<float> &redloc, array2D<float> &greenloc, array2D<float> &blueloc, int bfw, int bfh, double &avg_rm, double &avg_gm, double &avg_bm, double &tempitc, double &greenitc, float &temp0, float &delta,  int &bia, int &dread, int &kcam, int &nocam, float &studgood, float &minchrom, int &kmin, float &minhist, float &maxhist, bool &twotimes, const procparams::WBParams & wbpar, int begx, int begy, int yEn, int xEn, int cx, int cy, const procparams::ColorManagementParams &cmp, const procparams::RAWParams &raw, const procparams::ToneCurveParams &hrp) = 0;
    virtual void        getrgbloc(int begx, int begy, int yEn, int xEn, int cx, int cy, int bf_h, int bf_w, const procparams::WBParams & wbpar) = 0;

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
    virtual bool           isGainMapSupported () const = 0;
    virtual DCPProfile*    getDCP (const procparams::ColorManagementParams &cmp, DCPProfileApplyState &as)
    {
        return nullptr;
    };

    virtual void        setProgressListener (ProgressListener* pl) {}

    void        increaseRef () final
    {
        references++;
    }
    void        decreaseRef () final
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
    virtual void getAutoMatchedToneCurve(const procparams::ColorManagementParams &cp, const procparams::RAWParams &rawParams, StandardObserver observer, std::vector<double> &outCurve)
    {
        outCurve = { 0.0 };
    }

    double getDirPyrDenoiseExpComp () const
    {
        return dirpyrdenoiseExpComp;
    }
    // functions inherited from the InitialImage interface
    Glib::ustring getFileName() const final override
    {
        return fileName;
    }
    cmsHPROFILE getEmbeddedProfile() const final override
    {
        return embProfile;
    }
    const FramesMetaData *getMetaData() const final override
    {
        return idata;
    }
    ImageSource* getImageSource () override
    {
        return this;
    }
    virtual void getRawValues(int x, int y, int rotate, int &R, int &G, int &B) = 0;
    virtual void captureSharpening(const procparams::CaptureSharpeningParams &sharpeningParams, bool showMask, double &conrastThreshold, double &radius) = 0;
    virtual void wbMul2Camera(double &rm, double &gm, double &bm) = 0;
    virtual void wbCamera2Mul(double &rm, double &gm, double &bm) = 0;

};

}
