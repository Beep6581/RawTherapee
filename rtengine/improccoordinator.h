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
private:


protected:
    Imagefloat *orig_prev;
    Imagefloat *oprevi;
    LabImage *oprevl;
    LabImage *nprevl;
    Image8 *previmg;  // displayed image in monitor color space, showing the output profile as well (soft-proofing enabled, which then correspond to workimg) or not
    Image8 *workimg;  // internal image in output color space for analysis
    CieImage *ncie;

    ImageSource* imgsrc;

    SHMap* shmap;

    ColorTemp currWB;
    ColorTemp autoWB;

    double lastAwbEqual;
    double lastAwbTempBias;

    ImProcFunctions ipf;

    Glib::ustring monitorProfile;
    RenderingIntent monitorIntent;
    bool softProof;
    bool gamutCheck;

    int scale;
    bool highDetailPreprocessComputed;
    bool highDetailRawComputed;
    bool allocated;

    void freeAll ();

    // Precomputed values used by DetailedCrop ----------------------------------------------

    float bwAutoR, bwAutoG, bwAutoB;
    float CAMMean;
    int coordX, coordY, localX, localY;

    LUTf hltonecurve;
    LUTf shtonecurve;
    LUTf tonecurve;

    LUTf lumacurve;
//    LUTf localcurve;
    LUTf chroma_acurve;
    LUTf chroma_bcurve;
    LUTf satcurve;
    LUTf lhskcurve;
    LUTf clcurve;
//    multi_array2D<float, 3> conversionBuffer;
    multi_array2D<float, 4> conversionBuffer;
    LUTf wavclCurve;
    LUTf clToningcurve;
    LUTf lllocalcurve;
    LUTf cclocalcurve;
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
    LocretigainCurve locRETgainCurve;
    LocretigainCurverab locRETgainCurverab;
    LocLHCurve loclhCurve;

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
    localListener* aloListener;
    AutoWBListener* awbListener;
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
    void updatePreviewImage (int todo, Crop* cropCall = nullptr);

    MyMutex mProcessing;
    ProcParams params;

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
    ProcParams nextParams;
    ProcParams nextParams2;
    bool destroying;
    bool utili;
    bool locutili;
    bool autili;
    bool butili;
    bool ccutili;
    bool cclutili;
    bool clcutili;
    bool opautili;
    bool wavcontlutili;
    bool locallutili;
    bool localcutili;

    int **dataspot;
    std::string *retistr;
    std::string *llstr;
    std::string *lhstr;
    std::string *ccstr;

    LUTi circrads;
    LUTi centerx;
    LUTi centery;
    LUTi locx;
    LUTi locy;
    LUTi locxl;
    LUTi locyt;
    LUTi lights;
    LUTi contrs;
    LUTi chroms;
    LUTi sensis;
    LUTi transits;
    LUTi inverss;
    LUTi curvactivs;
    LUTi smeths;
    LUTi curens;
    LUTi radiuss;
    LUTi strengths;
    LUTi sensibns;
    LUTi inversrads;
    LUTi strs;
    LUTi chrrts;
    LUTi neighs;
    LUTi varts;
    LUTi sensihs;
    LUTi inversrets;
    LUTi retinexs;
    LUTi sps;
    LUTi sharradiuss;
    LUTi sharamounts;
    LUTi shardampings;
    LUTi inversshas;
    LUTi shariters;
    LUTi sensishas;
    LUTi qualitys;
    LUTi thress;
    LUTi proxis;
    LUTi noiselumfs;
    LUTi noiselumcs;
    LUTi noisechrofs;
    LUTi noisechrocs;
    LUTi mult0s;
    LUTi mult1s;
    LUTi mult2s;
    LUTi mult3s;
    LUTi mult4s;
    LUTi thresholds;
    LUTi sensicbs;
    LUTi activlums;
    LUTi strens;
    LUTi gammas;
    LUTi estops;
    LUTi scaltms;
    LUTi reweis;
    LUTi sensitms;
    LUTi qualitycurves;

    int versionmip;
    LUTi sizeretics;
    LUTi reticurvs;
    LUTi retrabs;
    LUTi llcurvs;
    LUTi sizellcs;
    LUTi lhcurvs;
    LUTi sizelhcs;
    LUTi cccurvs;
    LUTi sizecccs;

    LUTf huerefs;
    LUTf chromarefs;
    LUTf lumarefs;

    double huer, chromar, lumar;
    void startProcessing ();
    void process ();
    float colourToningSatLimit;
    float colourToningSatLimitOpacity;

public:

    ImProcCoordinator ();
    ~ImProcCoordinator ();
    void assign     (ImageSource* imgsrc);

    void        getParams (procparams::ProcParams* dst)
    {
        *dst = params;
    }

    void        startProcessing (int changeCode);
    ProcParams* beginUpdateParams ();
    void        endUpdateParams (ProcEvent change);  // must be called after beginUpdateParams, triggers update
    void        endUpdateParams (int changeFlags);
    void        stopProcessing ();
//    void updatePreviewImage (int todo, Crop* cropCall = NULL);

    std::string *retistrsav;

    void setPreviewScale    (int scale)
    {
        setScale (scale);
    }
    int  getPreviewScale    ()
    {
        return scale;
    }

    //void fullUpdatePreviewImage  ();

    int getFullWidth ()
    {
        return fullw;
    }
    int getFullHeight ()
    {
        return fullh;
    }

    int getPreviewWidth ()
    {
        return pW;
    }
    int getPreviewHeight ()
    {
        return pH;
    }

    DetailedCrop* createCrop  (::EditDataProvider *editDataProvider, bool isDetailWindow);

    bool getAutoWB   (double& temp, double& green, double equal, double tempBias);
    void getCamWB    (double& temp, double& green);
    void getSpotWB   (int x, int y, int rectSize, double& temp, double& green);
    void getAutoCrop (double ratio, int &x, int &y, int &w, int &h);

    void setMonitorProfile (const Glib::ustring& profile, RenderingIntent intent);
    void getMonitorProfile (Glib::ustring& profile, RenderingIntent& intent) const;
    void setSoftProofing   (bool softProof, bool gamutCheck);
    void getSoftProofing   (bool &softProof, bool &gamutCheck);

    bool updateTryLock ()
    {
        return updaterThreadStart.trylock();
    }
    void updateUnLock ()
    {
        updaterThreadStart.unlock();
    }

    void setProgressListener (ProgressListener* pl)
    {
        plistener = pl;
    }
    void setPreviewImageListener    (PreviewImageListener* il)
    {
        imageListener = il;
    }
    void setSizeListener     (SizeListener* il)
    {
        sizeListeners.push_back (il);
    }
    void delSizeListener     (SizeListener* il)
    {
        std::vector<SizeListener*>::iterator it = std::find (sizeListeners.begin(), sizeListeners.end(), il);

        if (it != sizeListeners.end()) {
            sizeListeners.erase (it);
        }
    }
    void setAutoExpListener  (AutoExpListener* ael)
    {
        aeListener = ael;
    }
    void setHistogramListener (HistogramListener *h)
    {
        hListener = h;
    }
    void setAutoCamListener  (AutoCamListener* acl)
    {
        acListener = acl;
    }
    void setAutoBWListener   (AutoBWListener* abw)
    {
        abwListener = abw;
    }
    void setlocalListener   (localListener* alo)
    {
        aloListener = alo;
    }

    void setAutoWBListener   (AutoWBListener* awb)
    {
        awbListener = awb;
    }
    void setAutoColorTonListener   (AutoColorTonListener* bwct)
    {
        actListener = bwct;
    }
    void setAutoChromaListener  (AutoChromaListener* adn)
    {
        adnListener = adn;
    }
    void setRetinexListener  (RetinexListener* adh)
    {
        dehaListener = adh;
    }
    void setWaveletListener  (WaveletListener* awa)
    {
        awavListener = awa;
    }

    void setFrameCountListener  (FrameCountListener* fcl)
    {
        frameCountListener = fcl;
    }

    void setImageTypeListener  (ImageTypeListener* itl)
    {
        imageTypeListener = itl;
    }

    void saveInputICCReference (const Glib::ustring& fname, bool apply_wb);

    InitialImage*  getInitialImage ()
    {
        return imgsrc;
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
