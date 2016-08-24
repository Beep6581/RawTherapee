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

protected:
    Imagefloat *orig_prev;
    Imagefloat *oprevi;
    LabImage *oprevl;
    LabImage *nprevl;
    Image8 *previmg;
    Image8 *workimg;
    CieImage *ncie;

    ImageSource* imgsrc;

    SHMap* shmap;

    ColorTemp currWB;
    ColorTemp autoWB;

    double lastAwbEqual;

    ImProcFunctions ipf;

    Glib::ustring monitorProfile;

    RenderingIntent monitorIntent;

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
    float chaut, redaut, blueaut, maxredaut, maxblueaut,  minredaut, minblueaut, nresi, highresi, chromina, sigma, lumema;

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
    void updatePreviewImage (int todo, Crop* cropCall = NULL);

    MyMutex mProcessing;
    ProcParams params;

    // members of the updater:
    Glib::Thread* thread;
    MyMutex updaterThreadStart;
    MyMutex paramsUpdateMutex;
    int  changeSinceLast;
    bool updaterRunning;
    ProcParams nextParams;
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

public:

    ImProcCoordinator ();
    ~ImProcCoordinator ();
    void assign     (ImageSource* imgsrc);

    void        getParams (procparams::ProcParams* dst)
    {
        *dst = params;
    }

    void        startProcessing(int changeCode);
    ProcParams* beginUpdateParams ();
    void        endUpdateParams (ProcEvent change);  // must be called after beginUpdateParams, triggers update
    void        endUpdateParams (int changeFlags);
    void        stopProcessing ();


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

    bool getAutoWB   (double& temp, double& green, double equal);
    void getCamWB    (double& temp, double& green);
    void getSpotWB   (int x, int y, int rectSize, double& temp, double& green);
    void getAutoCrop (double ratio, int &x, int &y, int &w, int &h);

    void setMonitorProfile (const Glib::ustring& profile, RenderingIntent intent);
    void getMonitorProfile (Glib::ustring& profile, RenderingIntent& intent) const;

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
    void setHistogramListener(HistogramListener *h)
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

    void saveInputICCReference (const Glib::ustring& fname, bool apply_wb);

    InitialImage*  getInitialImage ()
    {
        return imgsrc;
    }
};
}
#endif
