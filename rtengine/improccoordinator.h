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
#ifndef _IMPROCCOORDINATOR_H_
#define _IMPROCCOORDINATOR_H_

#include <memory>

#include "rtengine.h"
#include "improcfun.h"
#include "image8.h"
#include "image16.h"
#include "imagesource.h"
#include "procevents.h"
#include "dcrop.h"
#include "LUT.h"
#include "../rtgui/threadutils.h"

namespace rtengine
{

using namespace procparams;

class Crop;

/** @brief Manages the image processing, espc. of the preview windows
  *
  * There is one ImProcCoordinator per edit panel.
  *
  * The ImProcCoordinator handle an sized down image representation of the full image, that is used when paning
  * and in the Navigator object.
  *
  * Each ImProcCoordinator handles an rtengine::Crop list, which process images too with their own pipeline,
  * but using this class' LUT and other precomputed parameters. The main preview area is displaying a non framed Crop object,
  * while detail windows are framed Crop objects.
  */
class ImProcCoordinator : public StagedImageProcessor
{

    friend class Crop;

protected:
    Imagefloat *orig_prev;
    Imagefloat *oprevi;
    LabImage *oprevl;
    LabImage *nprevl;
    Imagefloat *fattal_11_dcrop_cache; // global cache for ToneMapFattal02 used in 1:1 detail windows (except when denoise is active)
    Image8 *previmg;  // displayed image in monitor color space, showing the output profile as well (soft-proofing enabled, which then correspond to workimg) or not
    Image8 *workimg;  // internal image in output color space for analysis
    CieImage *ncie;

    ImageSource* imgsrc;

    ColorTemp currWB;
    ColorTemp autoWB;

    double lastAwbEqual;
    double lastAwbTempBias;

    Glib::ustring monitorProfile;
    RenderingIntent monitorIntent;
    bool softProof;
    bool gamutCheck;
    bool sharpMask;

    int scale;
    bool highDetailPreprocessComputed;
    bool highDetailRawComputed;
    bool allocated;

    void freeAll ();

    // Precomputed values used by DetailedCrop ----------------------------------------------

    float bwAutoR, bwAutoG, bwAutoB;
    float CAMMean;
    LUTf hltonecurve;
    LUTf shtonecurve;
    LUTf tonecurve;

    LUTf lumacurve;
    LUTf chroma_acurve;
    LUTf chroma_bcurve;
    LUTf satcurve;
    LUTf lhskcurve;
    LUTf clcurve;
//    multi_array2D<float, 3> conversionBuffer;
    multi_array2D<float, 4> conversionBuffer;
    LUTf wavclCurve;
    LUTf clToningcurve;
    LUTf cl2Toningcurve;
    LUTf Noisecurve;
    LUTf NoiseCCcurve;

    LUTu vhist16, vhist16bw;
    LUTu lhist16CAM;
    LUTu lhist16CCAM;
    LUTu lhist16RETI;
    LUTu lhist16CLlad, lhist16LClad;
    LUTu histRed, histRedRaw;
    LUTu histGreen, histGreenRaw;
    LUTu histBlue, histBlueRaw;
    LUTu histLuma, histToneCurve, histToneCurveBW, histLCurve, histCCurve;
    LUTu histLLCurve, histLCAM, histCCAM, histClad, bcabhist, histChroma, histLRETI;

    LUTf CAMBrightCurveJ, CAMBrightCurveQ;

    LUTf rCurve;
    LUTf gCurve;
    LUTf bCurve;
    ToneCurve customToneCurve1;
    ToneCurve customToneCurve2;
    ColorGradientCurve ctColorCurve;
    OpacityCurve ctOpacityCurve;
    NoiseCurve noiseLCurve;
    NoiseCurve noiseCCurve;
    WavCurve wavCLVCurve;
    WavOpacityCurveRG waOpacityCurveRG;
    WavOpacityCurveBY waOpacityCurveBY;
    WavOpacityCurveW waOpacityCurveW;
    WavOpacityCurveWL waOpacityCurveWL;
    RetinextransmissionCurve dehatransmissionCurve;
    RetinexgaintransmissionCurve dehagaintransmissionCurve;

    ColorAppearance customColCurve1;
    ColorAppearance customColCurve2;
    ColorAppearance customColCurve3;
    ToneCurve beforeToneCurveBW;
    ToneCurve afterToneCurveBW;

    LUTu rcurvehist, rcurvehistCropped, rbeforehist;
    LUTu gcurvehist, gcurvehistCropped, gbeforehist;
    LUTu bcurvehist, bcurvehistCropped, bbeforehist;

    // ------------------------------------------------------------------------------------

    int fw, fh, tr, fullw, fullh;
    int pW, pH;

    ProgressListener* plistener;
    PreviewImageListener* imageListener;
    AutoExpListener* aeListener;
    AutoCamListener* acListener;
    AutoBWListener* abwListener;
    AutoWBListener* awbListener;
    FlatFieldAutoClipListener *flatFieldAutoClipListener;
    AutoContrastListener *bayerAutoContrastListener;
    AutoContrastListener *xtransAutoContrastListener;
    FrameCountListener *frameCountListener;
    ImageTypeListener *imageTypeListener;

    AutoColorTonListener* actListener;
    AutoChromaListener* adnListener;
    WaveletListener* awavListener;
    RetinexListener* dehaListener;

    HistogramListener* hListener;
    std::vector<SizeListener*> sizeListeners;

    std::vector<Crop*> crops;

    bool resultValid;

    MyMutex minit;  // to gain mutually exclusive access to ... to what exactly?

    void progress (Glib::ustring str, int pr);
    void reallocAll ();
    void updateLRGBHistograms ();
    void setScale (int prevscale);
    void updatePreviewImage (int todo, bool panningRelatedChange);

    MyMutex mProcessing;
    const std::unique_ptr<ProcParams> params;

    // for optimization purpose, the output profile, output rendering intent and
    // output BPC will trigger a regeneration of the profile on parameter change only
    // and automatically
    Glib::ustring lastOutputProfile;
    RenderingIntent lastOutputIntent;
    bool lastOutputBPC;

    // members of the updater:
    Glib::Thread* thread;
    MyMutex updaterThreadStart;
    MyMutex paramsUpdateMutex;
    int  changeSinceLast;
    bool updaterRunning;
    const std::unique_ptr<ProcParams> nextParams;
    bool destroying;
    bool utili;
    bool autili;
    bool butili;
    bool ccutili;
    bool cclutili;
    bool clcutili;
    bool opautili;
    bool wavcontlutili;
    void startProcessing ();
    void process ();
    float colourToningSatLimit;
    float colourToningSatLimitOpacity;
    bool highQualityComputed;
    cmsHTRANSFORM customTransformIn;
    cmsHTRANSFORM customTransformOut;

    ImProcFunctions ipf;

public:

    ImProcCoordinator ();
    ~ImProcCoordinator () override;
    void assign     (ImageSource* imgsrc);

    void        getParams (procparams::ProcParams* dst) override;

    void        startProcessing (int changeCode) override;
    ProcParams* beginUpdateParams () override;
    void        endUpdateParams (ProcEvent change) override;  // must be called after beginUpdateParams, triggers update
    void        endUpdateParams (int changeFlags) override;
    void        stopProcessing () override;


    void setPreviewScale    (int scale) override
    {
        setScale (scale);
    }
    int  getPreviewScale    () override
    {
        return scale;
    }

    //void fullUpdatePreviewImage  ();

    int getFullWidth () override
    {
        return fullw;
    }
    int getFullHeight () override
    {
        return fullh;
    }

    int getPreviewWidth () override
    {
        return pW;
    }
    int getPreviewHeight () override
    {
        return pH;
    }

    DetailedCrop* createCrop  (::EditDataProvider *editDataProvider, bool isDetailWindow) override;

    bool getAutoWB   (double& temp, double& green, double equal, double tempBias) override;
    void getCamWB    (double& temp, double& green) override;
    void getSpotWB   (int x, int y, int rectSize, double& temp, double& green) override;
    void getAutoCrop (double ratio, int &x, int &y, int &w, int &h) override;
    bool getHighQualComputed() override;
    void setHighQualComputed() override;
    void setMonitorProfile (const Glib::ustring& profile, RenderingIntent intent) override;
    void getMonitorProfile (Glib::ustring& profile, RenderingIntent& intent) const override;
    void setSoftProofing   (bool softProof, bool gamutCheck) override;
    void getSoftProofing   (bool &softProof, bool &gamutCheck) override;
    void setSharpMask      (bool sharpMask) override;
    bool updateTryLock () override
    {
        return updaterThreadStart.trylock();
    }
    void updateUnLock () override
    {
        updaterThreadStart.unlock();
    }

    void setProgressListener (ProgressListener* pl) override
    {
        plistener = pl;
    }
    void setPreviewImageListener    (PreviewImageListener* il) override
    {
        imageListener = il;
    }
    void setSizeListener     (SizeListener* il) override
    {
        sizeListeners.push_back (il);
    }
    void delSizeListener     (SizeListener* il) override
    {
        std::vector<SizeListener*>::iterator it = std::find (sizeListeners.begin(), sizeListeners.end(), il);

        if (it != sizeListeners.end()) {
            sizeListeners.erase (it);
        }
    }
    void setAutoExpListener  (AutoExpListener* ael) override
    {
        aeListener = ael;
    }
    void setHistogramListener (HistogramListener *h) override
    {
        hListener = h;
    }
    void setAutoCamListener  (AutoCamListener* acl) override
    {
        acListener = acl;
    }
    void setAutoBWListener   (AutoBWListener* abw) override
    {
        abwListener = abw;
    }
    void setAutoWBListener   (AutoWBListener* awb) override
    {
        awbListener = awb;
    }
    void setAutoColorTonListener   (AutoColorTonListener* bwct) override
    {
        actListener = bwct;
    }
    void setAutoChromaListener  (AutoChromaListener* adn) override
    {
        adnListener = adn;
    }
    void setRetinexListener  (RetinexListener* adh) override
    {
        dehaListener = adh;
    }
    void setWaveletListener  (WaveletListener* awa) override
    {
        awavListener = awa;
    }

    void setFrameCountListener  (FrameCountListener* fcl) override
    {
        frameCountListener = fcl;
    }

    void setFlatFieldAutoClipListener  (FlatFieldAutoClipListener* ffacl) override
    {
        flatFieldAutoClipListener = ffacl;
    }
    void setBayerAutoContrastListener  (AutoContrastListener* acl) override
    {
        bayerAutoContrastListener = acl;
    }

    void setXtransAutoContrastListener  (AutoContrastListener* acl) override
    {
        xtransAutoContrastListener = acl;
    }

    void setImageTypeListener  (ImageTypeListener* itl) override
    {
        imageTypeListener = itl;
    }

    void saveInputICCReference (const Glib::ustring& fname, bool apply_wb) override;

    InitialImage*  getInitialImage () override
    {
        return imgsrc;
    }

    cmsHTRANSFORM& getCustomTransformIn ()
    {
        return customTransformIn;
    }

    cmsHTRANSFORM& getCustomTransformOut ()
    {
        return customTransformOut;
    }

    struct DenoiseInfoStore {
        DenoiseInfoStore () : chM (0), max_r{}, max_b{}, ch_M{}, valid (false)  {}
        float chM;
        float max_r[9];
        float max_b[9];
        float ch_M[9];
        bool valid;

    } denoiseInfoStore;

};
}
#endif
