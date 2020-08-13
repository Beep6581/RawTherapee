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
#include <fstream>
#include <iostream>
#include <string>
#include <glibmm/thread.h>

#include "improccoordinator.h"

#include "cieimage.h"
#include "color.h"
#include "colortemp.h"
#include "jaggedarray.h"
#include "curves.h"
#include "dcp.h"
#include "iccstore.h"
#include "image8.h"
#include "imagefloat.h"
#include "improcfun.h"
#include "labimage.h"
#include "lcp.h"
#include "procparams.h"
#include "refreshmap.h"
#include "guidedfilter.h"

#include "../rtgui/options.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
using rtengine::Coord2D;
Coord2D translateCoord(const rtengine::ImProcFunctions& ipf, int fw, int fh, int x, int y) {

    const std::vector<Coord2D> points = {Coord2D(x, y)};

    std::vector<Coord2D> red;
    std::vector<Coord2D> green;
    std::vector<Coord2D> blue;
    ipf.transCoord(fw, fh, points, red, green, blue);

    return green[0];
}

}

namespace rtengine
{

ImProcCoordinator::ImProcCoordinator() :
    orig_prev(nullptr),
    oprevi(nullptr),
    oprevl(nullptr),
    nprevl(nullptr),
    fattal_11_dcrop_cache(nullptr),
    previmg(nullptr),
    workimg(nullptr),
    ncie(nullptr),
    imgsrc(nullptr),
    lastAwbEqual(0.),
    lastAwbTempBias(0.0),
    lastAwbauto(""),
    monitorIntent(RI_RELATIVE),
    softProof(false),
    gamutCheck(false),
    sharpMask(false),
    sharpMaskChanged(false),
    scale(10),
    highDetailPreprocessComputed(false),
    highDetailRawComputed(false),
    allocated(false),
    bwAutoR(-9000.f),
    bwAutoG(-9000.f),
    bwAutoB(-9000.f),
    CAMMean(NAN),
    hltonecurve(65536),
    shtonecurve(65536),
    tonecurve(65536, 0),  //,1);
    lumacurve(32770, 0),  // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
    chroma_acurve(65536, 0),
    chroma_bcurve(65536, 0),
    satcurve(65536, 0),
    lhskcurve(65536, 0),
    clcurve(65536, 0),
    conversionBuffer(1, 1),
    wavclCurve(65536, 0),
    clToningcurve(65536, 0),
    cl2Toningcurve(65536, 0),
    Noisecurve(65536, 0),
    NoiseCCcurve(65536, 0),
    vhist16(65536), vhist16bw(65536),
    lhist16CAM(65536),
    lhist16CCAM(65536),
    lhist16RETI(),
    lhist16LClad(65536),
    histRed(256), histRedRaw(256),
    histGreen(256), histGreenRaw(256),
    histBlue(256), histBlueRaw(256),
    histLuma(256),
    histToneCurve(256),
    histToneCurveBW(256),
    histLCurve(256),
    histCCurve(256),
    histLLCurve(256),

    histLCAM(256),
    histCCAM(256),
    histClad(256),
    bcabhist(256),
    histChroma(256),

    histLRETI(256),

    CAMBrightCurveJ(), CAMBrightCurveQ(),

    rCurve(),
    gCurve(),
    bCurve(),
    ctColorCurve(),
    rcurvehist(256), rcurvehistCropped(256), rbeforehist(256),
    gcurvehist(256), gcurvehistCropped(256), gbeforehist(256),
    bcurvehist(256), bcurvehistCropped(256), bbeforehist(256),
    fw(0), fh(0), tr(0),
    fullw(1), fullh(1),
    pW(-1), pH(-1),
    plistener(nullptr),
    imageListener(nullptr),
    aeListener(nullptr),
    acListener(nullptr),
    abwListener(nullptr),
    awbListener(nullptr),
    flatFieldAutoClipListener(nullptr),
    bayerAutoContrastListener(nullptr),
    xtransAutoContrastListener(nullptr),
    pdSharpenAutoContrastListener(nullptr),
    pdSharpenAutoRadiusListener(nullptr),
    frameCountListener(nullptr),
    imageTypeListener(nullptr),
    filmNegListener(nullptr),
    actListener(nullptr),
    adnListener(nullptr),
    awavListener(nullptr),
    dehaListener(nullptr),
    hListener(nullptr),
    resultValid(false),
    params(new procparams::ProcParams),
    lastOutputProfile("BADFOOD"),
    lastOutputIntent(RI__COUNT),
    lastOutputBPC(false),
    thread(nullptr),
    changeSinceLast(0),
    updaterRunning(false),
    nextParams(new procparams::ProcParams),
    destroying(false),
    utili(false),
    autili(false),
    butili(false),
    ccutili(false),
    cclutili(false),
    clcutili(false),
    opautili(false),
    wavcontlutili(false),
    colourToningSatLimit(0.f),
    colourToningSatLimitOpacity(0.f),
    highQualityComputed(false),
    customTransformIn(nullptr),
    customTransformOut(nullptr),
    ipf(params.get(), true),

    // Locallab
    locallListener(nullptr),
    lllocalcurve(65536, LUT_CLIP_OFF),
    cllocalcurve(65536, LUT_CLIP_OFF),
    lclocalcurve(65536, LUT_CLIP_OFF),
    cclocalcurve(65536, LUT_CLIP_OFF),
    rgblocalcurve(65536, LUT_CLIP_OFF),
    exlocalcurve(65536, LUT_CLIP_OFF),
    hltonecurveloc(65536, LUT_CLIP_OFF), //32768
    shtonecurveloc(65536, LUT_CLIP_OFF),
    tonecurveloc(65536, LUT_CLIP_OFF),
    lightCurveloc(32770, LUT_CLIP_OFF),
    lmasklocalcurve(65536, LUT_CLIP_OFF),
    lmaskexplocalcurve(65536, LUT_CLIP_OFF),
    lmaskSHlocalcurve(65536, LUT_CLIP_OFF),
    lmaskviblocalcurve(65536, LUT_CLIP_OFF),
    lmasktmlocalcurve(65536, LUT_CLIP_OFF),
    lmaskretilocalcurve(65536, LUT_CLIP_OFF),
    lmaskcblocalcurve(65536, LUT_CLIP_OFF),
    lmaskbllocalcurve(65536, LUT_CLIP_OFF),
    lmasklclocalcurve(65536, LUT_CLIP_OFF),
    lmasklocal_curve(65536, LUT_CLIP_OFF),
    lastspotdup(false),
    previewDeltaE(false),
    locallColorMask(0),
    locallColorMaskinv(0),
    locallExpMask(0),
    locallExpMaskinv(0),
    locallSHMask(0),
    locallSHMaskinv(0),
    locallvibMask(0),
    localllcMask(0),
    locallcbMask(0),
    locallretiMask(0),
    locallsoftMask(0),
    localltmMask(0),
    locallblMask(0),
    locallsharMask(0),
    locall_Mask(0),
    retistrsav(nullptr)
{
}

ImProcCoordinator::~ImProcCoordinator()
{

    destroying = true;
    updaterThreadStart.lock();

    if (updaterRunning && thread) {
        thread->join();
    }

    mProcessing.lock();
    mProcessing.unlock();
    freeAll();

    if (fattal_11_dcrop_cache) {
        delete fattal_11_dcrop_cache;
        fattal_11_dcrop_cache = nullptr;
    }

    std::vector<Crop*> toDel = crops;

    for (size_t i = 0; i < toDel.size(); i++) {
        delete toDel[i];
    }

    imgsrc->decreaseRef();

    if (customTransformIn) {
        cmsDeleteTransform(customTransformIn);
        customTransformIn = nullptr;
    }

    if (customTransformOut) {
        cmsDeleteTransform(customTransformOut);
        customTransformOut = nullptr;
    }

    updaterThreadStart.unlock();
}

void ImProcCoordinator::assign(ImageSource* imgsrc)
{
    this->imgsrc = imgsrc;
}

void ImProcCoordinator::getParams(procparams::ProcParams* dst)
{
    *dst = *params;
}

DetailedCrop* ImProcCoordinator::createCrop(::EditDataProvider *editDataProvider, bool isDetailWindow)
{

    return new Crop(this, editDataProvider, isDetailWindow);
}


// todo: bitmask containing desired actions, taken from changesSinceLast
void ImProcCoordinator::updatePreviewImage(int todo, bool panningRelatedChange)
{
    // TODO Locallab printf

    MyMutex::MyLock processingLock(mProcessing);

    bool highDetailNeeded = options.prevdemo == PD_Sidecar ? true : (todo & M_HIGHQUAL);
                //    printf("metwb=%s \n", params->wb.method.c_str());

    // Check if any detail crops need high detail. If not, take a fast path short cut
    if (!highDetailNeeded) {
        for (size_t i = 0; i < crops.size(); i++) {
            if (crops[i]->get_skip() == 1) {   // skip=1 -> full  resolution
                highDetailNeeded = true;
                break;
            }
        }
    }

    if (((todo & ALL) == ALL) || (todo & M_MONITOR) || panningRelatedChange || (highDetailNeeded && options.prevdemo != PD_Sidecar)) {
        bwAutoR = bwAutoG = bwAutoB = -9000.f;

        if (todo == CROP && ipf.needsPCVignetting()) {
            todo |= TRANSFORM;    // Change about Crop does affect TRANSFORM
        }

        RAWParams rp = params->raw;
        ColorManagementParams cmp = params->icm;
        LCurveParams  lcur = params->labCurve;
        
        if (!highDetailNeeded) {
            // if below 100% magnification, take a fast path
            if (rp.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::NONE) && rp.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::MONO)) {
                rp.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
            }

            //bayerrp.all_enhance = false;

            if (rp.xtranssensor.method != RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::NONE) && rp.xtranssensor.method != RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::MONO)) {
                rp.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);
            }

            rp.bayersensor.ccSteps = 0;
            rp.xtranssensor.ccSteps = 0;
            //rp.deadPixelFilter = rp.hotPixelFilter = false;
        }

        if (frameCountListener) {
            frameCountListener->FrameCountChanged(imgsrc->getFrameCount(), params->raw.bayersensor.imageNum);
        }

        // raw auto CA is bypassed if no high detail is needed, so we have to compute it when high detail is needed
        if ((todo & M_PREPROC) || (!highDetailPreprocessComputed && highDetailNeeded)) {
            imgsrc->setCurrentFrame(params->raw.bayersensor.imageNum);

            imgsrc->preprocess(rp, params->lensProf, params->coarse);

            if (flatFieldAutoClipListener && rp.ff_AutoClipControl) {
                flatFieldAutoClipListener->flatFieldAutoClipValueChanged(imgsrc->getFlatFieldAutoClipValue());
            }

            imgsrc->getRAWHistogram(histRedRaw, histGreenRaw, histBlueRaw);

            highDetailPreprocessComputed = highDetailNeeded;

            // After preprocess, run film negative processing if enabled
            if (
                (todo & M_RAW)
                && (
                    imgsrc->getSensorType() == ST_BAYER
                    || imgsrc->getSensorType() == ST_FUJI_XTRANS
                )
                && params->filmNegative.enabled
            ) {
                std::array<float, 3> filmBaseValues = {
                    static_cast<float>(params->filmNegative.redBase),
                    static_cast<float>(params->filmNegative.greenBase),
                    static_cast<float>(params->filmNegative.blueBase)
                };
                imgsrc->filmNegativeProcess(params->filmNegative, filmBaseValues);
                if (filmNegListener && params->filmNegative.redBase <= 0.f) {
                    filmNegListener->filmBaseValuesChanged(filmBaseValues);
                }
            }
        }

        /*
        Demosaic is kicked off only when
        Detail considerations:
            accurate detail is not displayed yet needed based on preview specifics (driven via highDetailNeeded flag)
        OR
        HLR considerations:
            Color HLR alters rgb output of demosaic, so re-demosaic is needed when Color HLR is being turned off;
            if HLR is enabled and changing method *from* Color to any other method
            OR HLR gets disabled when Color method was selected
        */
        // If high detail (=100%) is newly selected, do a demosaic update, since the last was just with FAST

        if (imageTypeListener) {
            imageTypeListener->imageTypeChanged(imgsrc->isRAW(), imgsrc->getSensorType() == ST_BAYER, imgsrc->getSensorType() == ST_FUJI_XTRANS, imgsrc->isMono());
        }

        if ((todo & M_RAW)
                || (!highDetailRawComputed && highDetailNeeded)
                || (params->toneCurve.hrenabled && params->toneCurve.method != "Color" && imgsrc->isRGBSourceModified())
                || (!params->toneCurve.hrenabled && params->toneCurve.method == "Color" && imgsrc->isRGBSourceModified())) {

            if (settings->verbose) {
                if (imgsrc->getSensorType() == ST_BAYER) {
                    printf("Demosaic Bayer image n.%d using method: %s\n", rp.bayersensor.imageNum + 1, rp.bayersensor.method.c_str());
                } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                    printf("Demosaic X-Trans image with using method: %s\n", rp.xtranssensor.method.c_str());
                }
            }

            if (imgsrc->getSensorType() == ST_BAYER) {
                if (params->raw.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::PIXELSHIFT)) {
                    imgsrc->setBorder(params->raw.bayersensor.border);
                } else {
                    imgsrc->setBorder(std::max(params->raw.bayersensor.border, 2));
                }
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                imgsrc->setBorder(params->raw.xtranssensor.border);
            }

            bool autoContrast = imgsrc->getSensorType() == ST_BAYER ? params->raw.bayersensor.dualDemosaicAutoContrast : params->raw.xtranssensor.dualDemosaicAutoContrast;
            double contrastThreshold = imgsrc->getSensorType() == ST_BAYER ? params->raw.bayersensor.dualDemosaicContrast : params->raw.xtranssensor.dualDemosaicContrast;
            imgsrc->demosaic(rp, autoContrast, contrastThreshold, params->pdsharpening.enabled);

            if (imgsrc->getSensorType() == ST_BAYER && bayerAutoContrastListener && autoContrast) {
                bayerAutoContrastListener->autoContrastChanged(contrastThreshold);
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS && xtransAutoContrastListener && autoContrast) {

                xtransAutoContrastListener->autoContrastChanged(contrastThreshold);
            }

            // if a demosaic happened we should also call getimage later, so we need to set the M_INIT flag
            todo |= (M_INIT | M_CSHARP);

        }

        if ((todo & (M_RAW | M_CSHARP)) && params->pdsharpening.enabled) {
            double pdSharpencontrastThreshold = params->pdsharpening.contrast;
            double pdSharpenRadius = params->pdsharpening.deconvradius;
            imgsrc->captureSharpening(params->pdsharpening, sharpMask, pdSharpencontrastThreshold, pdSharpenRadius);

            if (pdSharpenAutoContrastListener && params->pdsharpening.autoContrast) {
                pdSharpenAutoContrastListener->autoContrastChanged(pdSharpencontrastThreshold);
            }

            if (pdSharpenAutoRadiusListener && params->pdsharpening.autoRadius) {
                pdSharpenAutoRadiusListener->autoRadiusChanged(pdSharpenRadius);
            }
        }


        if ((todo & M_RAW)
                || (!highDetailRawComputed && highDetailNeeded)
                || (params->toneCurve.hrenabled && params->toneCurve.method != "Color" && imgsrc->isRGBSourceModified())
                || (!params->toneCurve.hrenabled && params->toneCurve.method == "Color" && imgsrc->isRGBSourceModified())) {
            if (highDetailNeeded) {
                highDetailRawComputed = true;
            } else {
                highDetailRawComputed = false;
            }

            if (params->retinex.enabled) {
                lhist16RETI(32768);
                lhist16RETI.clear();

                imgsrc->retinexPrepareBuffers(params->icm, params->retinex, conversionBuffer, lhist16RETI);
            }
        }

        if (todo & (M_INIT | M_LINDENOISE | M_HDR)) {
            if (params->wb.method == "autitcgreen") {
                imgsrc->getrgbloc(0, 0, fh, fw, 0, 0, fh, fw);
            }
        }

        if ((todo & (M_RETINEX | M_INIT)) && params->retinex.enabled) {
            bool dehacontlutili = false;
            bool mapcontlutili = false;
            bool useHsl = false;
            LUTf cdcurve(65536, 0);
            LUTf mapcurve(65536, 0);

            imgsrc->retinexPrepareCurves(params->retinex, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, dehacontlutili, mapcontlutili, useHsl, lhist16RETI, histLRETI);
            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            imgsrc->retinex(params->icm, params->retinex,  params->toneCurve, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, conversionBuffer, dehacontlutili, mapcontlutili, useHsl, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax, histLRETI);   //enabled Retinex

            if (dehaListener) {
                dehaListener->minmaxChanged(maxCD, minCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
            }
        }

        const bool autowb = (params->wb.method == "autold" || params->wb.method == "autitcgreen");
        if (settings->verbose) {
            printf("automethod=%s \n", params->wb.method.c_str());
        }
        if (todo & (M_INIT | M_LINDENOISE | M_HDR)) {
            MyMutex::MyLock initLock(minit);  // Also used in crop window

            imgsrc->HLRecovery_Global(params->toneCurve);   // this handles Color HLRecovery


            if (settings->verbose) {
                printf("Applying white balance, color correction & sRBG conversion...\n");
            }

            currWB = ColorTemp(params->wb.temperature, params->wb.green, params->wb.equal, params->wb.method);
            float studgood = 1000.f;

            if (!params->wb.enabled) {
                currWB = ColorTemp();
            } else if (params->wb.method == "Camera") {
                currWB = imgsrc->getWB();
                lastAwbauto = ""; //reinitialize auto
            } else if (autowb) {
                if (params->wb.method == "autitcgreen" || lastAwbEqual != params->wb.equal || lastAwbTempBias != params->wb.tempBias || lastAwbauto != params->wb.method) {
                    double rm, gm, bm;
                    double tempitc = 5000.f;
                    double greenitc = 1.;
                    currWBitc = imgsrc->getWB();
                    double tempref = currWBitc.getTemp() * (1. + params->wb.tempBias);
                    double greenref = currWBitc.getGreen();
                    if (settings->verbose && params->wb.method ==  "autitcgreen") {
                        printf("tempref=%f greref=%f\n", tempref, greenref);
                    }

                    imgsrc->getAutoWBMultipliersitc(tempref, greenref, tempitc, greenitc, studgood, 0, 0, fh, fw, 0, 0, fh, fw, rm, gm, bm,  params->wb, params->icm, params->raw);

                    if (params->wb.method ==  "autitcgreen") {
                        params->wb.temperature = tempitc;
                        params->wb.green = greenitc;
                        currWB = ColorTemp(params->wb.temperature, params->wb.green, 1., params->wb.method);
                        currWB.getMultipliers(rm, gm, bm);
                    }

                    if (rm != -1.) {
                        double bias = params->wb.tempBias;

                        if (params->wb.method ==  "autitcgreen") {
                            bias = 0.;
                        }

                        autoWB.update(rm, gm, bm, params->wb.equal, bias);
                        lastAwbEqual = params->wb.equal;
                        lastAwbTempBias = params->wb.tempBias;
                        lastAwbauto = params->wb.method;
                    } else {
                        lastAwbEqual = -1.;
                        lastAwbTempBias = 0.0;
                        lastAwbauto = "";
                        autoWB.useDefaults(params->wb.equal);
                    }
                    
                    
                }

                currWB = autoWB;
            }

            if (params->wb.enabled) {
                params->wb.temperature = currWB.getTemp();
                params->wb.green = currWB.getGreen();
            }

            if (autowb && awbListener && params->wb.method ==  "autitcgreen") {
                awbListener->WBChanged(params->wb.temperature, params->wb.green, studgood);
            } 

            if (autowb && awbListener && params->wb.method ==  "autold") {
                awbListener->WBChanged(params->wb.temperature, params->wb.green, -1.f);
            } 

            /*
                    GammaValues g_a;
                    double pwr = 1.0 / params->icm.gampos;
                    double ts = params->icm.slpos;


                    int mode = 0;
                    Color::calcGamma(pwr, ts, mode, g_a); // call to calcGamma with selected gamma and slope
                        printf("ga[0]=%f ga[1]=%f ga[2]=%f ga[3]=%f ga[4]=%f\n", g_a[0],g_a[1],g_a[2],g_a[3],g_a[4]);

                        Glib::ustring datal;
                        datal = "lutsrgb.txt";
                                ofstream fou(datal, ios::out | ios::trunc);

                    for(int i=0; i < 212; i++) {
                        //printf("igamma2=%i\n", (int) 65535.f*Color::igamma2(i/212.0));
                                float gam = Color::igamma2(i/211.0);
                                int lutga = nearbyint(65535.f* gam);
                              //  fou << 65535*(int)Color::igamma2(i/212.0) << endl;
                                fou << i << " " << lutga << endl;

                    }
                            fou.close();
            */
            int tr = getCoarseBitMask(params->coarse);

            imgsrc->getFullSize(fw, fh, tr);

            // Will (re)allocate the preview's buffers
            setScale(scale);
            PreviewProps pp(0, 0, fw, fh, scale);
            // Tells to the ImProcFunctions' tools what is the preview scale, which may lead to some simplifications
            ipf.setScale(scale);

            imgsrc->getImage(currWB, tr, orig_prev, pp, params->toneCurve, params->raw);
            denoiseInfoStore.valid = false;
            //ColorTemp::CAT02 (orig_prev, &params) ;
            //   printf("orig_prevW=%d\n  scale=%d",orig_prev->width, scale);
            /* Issue 2785, disabled some 1:1 tools
                    if (todo & M_LINDENOISE) {
                        DirPyrDenoiseParams denoiseParams = params->dirpyrDenoise;
                        if (denoiseParams.enabled && (scale==1)) {
                            Imagefloat *calclum = NULL ;

                            denoiseParams.getCurves(noiseLCurve,noiseCCurve);
                            int nbw=6;//nb tile W
                            int nbh=4;//

                            float ch_M[nbw*nbh];
                            float max_r[nbw*nbh];
                            float max_b[nbw*nbh];

                            if(denoiseParams.Lmethod == "CUR") {
                                if(noiseLCurve)
                                    denoiseParams.luma = 0.5f;
                                else
                                    denoiseParams.luma = 0.0f;
                            } else if(denoiseParams.Lmethod == "SLI")
                                noiseLCurve.Reset();


                            if(noiseLCurve || noiseCCurve){//only allocate memory if enabled and scale=1
                                // we only need image reduced to 1/4 here
                                calclum = new Imagefloat ((pW+1)/2, (pH+1)/2);//for luminance denoise curve
                                for(int ii=0;ii<pH;ii+=2){
                                    for(int jj=0;jj<pW;jj+=2){
                                        calclum->r(ii>>1,jj>>1) = orig_prev->r(ii,jj);
                                        calclum->g(ii>>1,jj>>1) = orig_prev->g(ii,jj);
                                        calclum->b(ii>>1,jj>>1) = orig_prev->b(ii,jj);
                                    }
                                }
                                imgsrc->convertColorSpace(calclum, params->icm, currWB);//calculate values after colorspace conversion
                            }

                            int kall=1;
                            ipf.RGB_denoise(kall, orig_prev, orig_prev, calclum, ch_M, max_r, max_b, imgsrc->isRAW(), denoiseParams, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, nresi, highresi);
                        }
                    }
            */
            imgsrc->convertColorSpace(orig_prev, params->icm, currWB);

            ipf.firstAnalysis(orig_prev, *params, vhist16);
        }

        if ((todo & M_HDR) && (params->fattal.enabled || params->dehaze.enabled)) {
            if (fattal_11_dcrop_cache) {
                delete fattal_11_dcrop_cache;
                fattal_11_dcrop_cache = nullptr;
            }

            ipf.dehaze(orig_prev, params->dehaze);
            ipf.ToneMapFattal02(orig_prev, params->fattal, 3, 0, nullptr, 0, 0, 0);

            if (oprevi != orig_prev) {
                delete oprevi;
            }
        }

        oprevi = orig_prev;

        // Remove transformation if unneeded
        bool needstransform = ipf.needsTransform(fw, fh, imgsrc->getRotateDegree(), imgsrc->getMetaData());

        if ((needstransform || ((todo & (M_TRANSFORM | M_RGBCURVE))  && params->dirpyrequalizer.cbdlMethod == "bef" && params->dirpyrequalizer.enabled && !params->colorappearance.enabled))) {
            assert(oprevi);
            Imagefloat *op = oprevi;
            oprevi = new Imagefloat(pW, pH);

            if (needstransform)
                ipf.transform(op, oprevi, 0, 0, 0, 0, pW, pH, fw, fh,
                              imgsrc->getMetaData(), imgsrc->getRotateDegree(), false);
            else {
                op->copyData(oprevi);
            }
        }

        if ((todo & (M_TRANSFORM | M_RGBCURVE))  && params->dirpyrequalizer.cbdlMethod == "bef" && params->dirpyrequalizer.enabled && !params->colorappearance.enabled) {
            const int W = oprevi->getWidth();
            const int H = oprevi->getHeight();
            LabImage labcbdl(W, H);
            ipf.rgb2lab(*oprevi, labcbdl, params->icm.workingProfile);
            ipf.dirpyrequalizer(&labcbdl, scale);
            ipf.lab2rgb(labcbdl, *oprevi, params->icm.workingProfile);
        }

        if (todo & M_AUTOEXP) {
            if (params->toneCurve.autoexp) {
                LUTu aehist;
                int aehistcompr;
                imgsrc->getAutoExpHistogram(aehist, aehistcompr);
                ipf.getAutoExp(aehist, aehistcompr, params->toneCurve.clip, params->toneCurve.expcomp,
                               params->toneCurve.brightness, params->toneCurve.contrast, params->toneCurve.black, params->toneCurve.hlcompr, params->toneCurve.hlcomprthresh);

                if (aeListener)
                    aeListener->autoExpChanged(params->toneCurve.expcomp, params->toneCurve.brightness, params->toneCurve.contrast,
                                               params->toneCurve.black, params->toneCurve.hlcompr, params->toneCurve.hlcomprthresh, params->toneCurve.hrenabled);
            }

            if (params->toneCurve.histmatching) {
                if (!params->toneCurve.fromHistMatching) {
                    imgsrc->getAutoMatchedToneCurve(params->icm, params->toneCurve.curve);
                }

                if (params->toneCurve.autoexp) {
                    params->toneCurve.expcomp = 0.0;
                }

                params->toneCurve.autoexp = false;
                params->toneCurve.curveMode = ToneCurveMode::FILMLIKE;
                params->toneCurve.curve2 = { 0 };
                params->toneCurve.brightness = 0;
                params->toneCurve.contrast = 0;
                params->toneCurve.black = 0;
                params->toneCurve.fromHistMatching = true;

                if (aeListener) {
                    aeListener->autoMatchedToneCurveChanged(params->toneCurve.curveMode, params->toneCurve.curve);
                }
            }

            // Encoding log with locallab
            if (params->locallab.enabled && !params->locallab.spots.empty()) {
                const int sizespot = (int)params->locallab.spots.size();

                float *sourceg = nullptr;
                sourceg = new float[sizespot];
                float *targetg = nullptr;
                targetg = new float[sizespot];
                bool *log = nullptr;
                log = new bool[sizespot];
                bool *autocomput = nullptr;
                autocomput = new bool[sizespot];
                float *blackev = nullptr;
                blackev = new float[sizespot];
                float *whiteev = nullptr;
                whiteev = new float[sizespot];
                bool *Autogr = nullptr;
                Autogr = new bool[sizespot];
                
                float *locx = nullptr;
                locx = new float[sizespot];
                float *locy = nullptr;
                locy = new float[sizespot];
                float *locxL = nullptr;
                locxL = new float[sizespot];
                float *locyT = nullptr;
                locyT = new float[sizespot];
                float *centx = nullptr;
                centx = new float[sizespot];
                float *centy = nullptr;
                centy = new float[sizespot];
                
                for (int sp = 0; sp < sizespot; sp++) {
                    log[sp] = params->locallab.spots.at(sp).explog;
                    autocomput[sp] = params->locallab.spots.at(sp).autocompute;
                    blackev[sp] = params->locallab.spots.at(sp).blackEv;
                    whiteev[sp] = params->locallab.spots.at(sp).whiteEv;
                    sourceg[sp] = params->locallab.spots.at(sp).sourceGray;
                    Autogr[sp] = params->locallab.spots.at(sp).Autogray;
                    targetg[sp] = params->locallab.spots.at(sp).targetGray;
                    locx[sp] = params->locallab.spots.at(sp).loc.at(0) / 2000.0;
                    locy[sp] = params->locallab.spots.at(sp).loc.at(2) / 2000.0;
                    locxL[sp] = params->locallab.spots.at(sp).loc.at(1) / 2000.0;
                    locyT[sp] = params->locallab.spots.at(sp).loc.at(3) / 2000.0;
                    centx[sp] = params->locallab.spots.at(sp).centerX / 2000.0 + 0.5;
                    centy[sp] = params->locallab.spots.at(sp).centerY / 2000.0 + 0.5;

                    const bool fullim = params->locallab.spots.at(sp).fullimage;

                    if (log[sp] && autocomput[sp]) {
                        constexpr int SCALE = 10;
                        int fw, fh, tr = TR_NONE;
                        imgsrc->getFullSize(fw, fh, tr);
                        PreviewProps pp(0, 0, fw, fh, SCALE);

                        float ysta = std::max(static_cast<float>(centy[sp] - locyT[sp]), 0.f);
                        float yend = std::min(static_cast<float>(centy[sp] + locy[sp]), 1.f);
                        float xsta = std::max(static_cast<float>(centx[sp] - locxL[sp]), 0.f);
                        float xend = std::min(static_cast<float>(centx[sp] + locx[sp]), 1.f);

                        if (fullim) {
                            ysta = 0.f;
                            yend = 1.f;
                            xsta = 0.f;
                            xend = 1.f;
                        }

                        ipf.getAutoLogloc(sp, imgsrc, sourceg, blackev, whiteev, Autogr, fw, fh, xsta, xend, ysta, yend, SCALE);

                        params->locallab.spots.at(sp).blackEv = blackev[sp];
                        params->locallab.spots.at(sp).whiteEv = whiteev[sp];
                        params->locallab.spots.at(sp).sourceGray = sourceg[sp];

                        if (locallListener) {
                            locallListener->logencodChanged(blackev[sp], whiteev[sp], sourceg[sp], targetg[sp]);
                        }
                    }
                }

                delete [] locx;
                delete [] locy;
                delete [] locxL;
                delete [] locyT;
                delete [] centx;
                delete [] centy;

                delete [] Autogr;
                delete [] whiteev;
                delete [] blackev;
                delete [] targetg;
                delete [] sourceg;
                delete [] log;
                delete [] autocomput;
            }
        }

        if (todo & (M_AUTOEXP | M_RGBCURVE)) {
            if (params->icm.workingTRC == "Custom") { //exec TRC IN free
                if (oprevi == orig_prev) {
                    oprevi = new Imagefloat(pW, pH);
                    orig_prev->copyData(oprevi);
                }

                const Glib::ustring profile = params->icm.workingProfile;

                if (profile == "sRGB" || profile == "Adobe RGB" || profile == "ProPhoto" || profile == "WideGamut" || profile == "BruceRGB" || profile == "Beta RGB" || profile == "BestRGB" || profile == "Rec2020" || profile == "ACESp0" || profile == "ACESp1") {
                    const int cw = oprevi->getWidth();
                    const int ch = oprevi->getHeight();

                    // put gamma TRC to 1
                    if (customTransformIn) {
                        cmsDeleteTransform(customTransformIn);
                        customTransformIn = nullptr;
                    }

                    ipf.workingtrc(oprevi, oprevi, cw, ch, -5, params->icm.workingProfile, 2.4, 12.92310, customTransformIn, true, false, true);

                    //adjust TRC
                    if (customTransformOut) {
                        cmsDeleteTransform(customTransformOut);
                        customTransformOut = nullptr;
                    }

                    ipf.workingtrc(oprevi, oprevi, cw, ch, 5, params->icm.workingProfile, params->icm.workingTRCGamma, params->icm.workingTRCSlope, customTransformOut, false, true, true);
                }
            }
        }


        if ((todo & M_RGBCURVE) || (todo & M_CROP)) {
            //complexCurve also calculated pre-curves histogram depending on crop
            CurveFactory::complexCurve(params->toneCurve.expcomp, params->toneCurve.black / 65535.0,
                                       params->toneCurve.hlcompr, params->toneCurve.hlcomprthresh,
                                       params->toneCurve.shcompr, params->toneCurve.brightness, params->toneCurve.contrast,
                                       params->toneCurve.curve, params->toneCurve.curve2,
                                       vhist16, hltonecurve, shtonecurve, tonecurve, histToneCurve, customToneCurve1, customToneCurve2, 1);

            CurveFactory::RGBCurve(params->rgbCurves.rcurve, rCurve, 1);
            CurveFactory::RGBCurve(params->rgbCurves.gcurve, gCurve, 1);
            CurveFactory::RGBCurve(params->rgbCurves.bcurve, bCurve, 1);


            opautili = false;

            if (params->colorToning.enabled) {
                TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
                double wp[3][3] = {
                    {wprof[0][0], wprof[0][1], wprof[0][2]},
                    {wprof[1][0], wprof[1][1], wprof[1][2]},
                    {wprof[2][0], wprof[2][1], wprof[2][2]}
                };
                params->colorToning.getCurves(ctColorCurve, ctOpacityCurve, wp, opautili);
                CurveFactory::diagonalCurve2Lut(params->colorToning.clcurve, clToningcurve, scale == 1 ? 1 : 16);
                CurveFactory::diagonalCurve2Lut(params->colorToning.cl2curve, cl2Toningcurve, scale == 1 ? 1 : 16);
            }

            if (params->blackwhite.enabled) {
                CurveFactory::curveBW(params->blackwhite.beforeCurve, params->blackwhite.afterCurve, vhist16bw, histToneCurveBW, beforeToneCurveBW, afterToneCurveBW, 1);
            }

            colourToningSatLimit = float (params->colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
            colourToningSatLimitOpacity = 1.f - (float (params->colorToning.saturatedOpacity) / 100.f);

            int satTH = 80;
            int satPR = 30;
            int indi = 0;

            if (params->colorToning.enabled  && params->colorToning.autosat && params->colorToning.method != "LabGrid") { //for colortoning evaluation of saturation settings
                float moyS = 0.f;
                float eqty = 0.f;
                ipf.moyeqt(oprevi, moyS, eqty); //return image : mean saturation and standard dev of saturation
                //printf("moy=%f ET=%f\n", moyS,eqty);
                float satp = ((moyS + 1.5f * eqty) - 0.3f) / 0.7f; //1.5 sigma ==> 93% pixels with high saturation -0.3 / 0.7 convert to Hombre scale

                if (satp >= 0.92f) {
                    satp = 0.92f;    //avoid values too high (out of gamut)
                }

                if (satp <= 0.15f) {
                    satp = 0.15f;    //avoid too low values
                }

                //satTH=(int) 100.f*satp;
                //satPR=(int) 100.f*(moyS-0.85f*eqty);//-0.85 sigma==>20% pixels with low saturation
                colourToningSatLimit = 100.f * satp;
                satTH = (int) 100.f * satp;

                colourToningSatLimitOpacity = 100.f * (moyS - 0.85f * eqty); //-0.85 sigma==>20% pixels with low saturation
                satPR = (int) 100.f * (moyS - 0.85f * eqty);
            }

            if (actListener && params->colorToning.enabled) {
                if (params->blackwhite.enabled && params->colorToning.autosat) {
                    actListener->autoColorTonChanged(0, satTH, satPR);    //hide sliders only if autosat
                    indi = 0;
                } else {
                    if (params->colorToning.autosat) {
                        if (params->colorToning.method == "Lab") {
                            indi = 1;
                        } else if (params->colorToning.method == "RGBCurves") {
                            indi = 1;
                        } else if (params->colorToning.method == "RGBSliders") {
                            indi = 1;
                        } else if (params->colorToning.method == "Splico") {
                            indi = 2;
                        } else if (params->colorToning.method == "Splitlr") {
                            indi = 2;
                        }
                    }
                }
            }

            // if it's just crop we just need the histogram, no image updates
            if (todo & M_RGBCURVE) {
                //initialize rrm bbm ggm different from zero to avoid black screen in some cases
                double rrm = 33.;
                double ggm = 33.;
                double bbm = 33.;

                DCPProfileApplyState as;
                DCPProfile *dcpProf = imgsrc->getDCP(params->icm, as);

                ipf.rgbProc(oprevi, oprevl, nullptr, hltonecurve, shtonecurve, tonecurve, params->toneCurve.saturation,
                            rCurve, gCurve, bCurve, colourToningSatLimit, colourToningSatLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, beforeToneCurveBW, afterToneCurveBW, rrm, ggm, bbm, bwAutoR, bwAutoG, bwAutoB, params->toneCurve.expcomp, params->toneCurve.hlcompr, params->toneCurve.hlcomprthresh, dcpProf, as, histToneCurve);

                if (params->blackwhite.enabled && params->blackwhite.autoc && abwListener) {
                    if (settings->verbose) {
                        printf("ImProcCoordinator / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", static_cast<double>(bwAutoR), static_cast<double>(bwAutoG), static_cast<double>(bwAutoB));
                    }

                    abwListener->BWChanged((float) rrm, (float) ggm, (float) bbm);
                }

                if (params->colorToning.enabled && params->colorToning.autosat && actListener) {
                    actListener->autoColorTonChanged(indi, (int) colourToningSatLimit, (int)colourToningSatLimitOpacity);  //change sliders autosat
                }

                // correct GUI black and white with value
            }

            //  ipf.Lab_Tile(oprevl, oprevl, scale);

            // compute L channel histogram
            int x1, y1, x2, y2;
            params->crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
        }

//    lhist16(32768);
        if (todo & (M_LUMACURVE | M_CROP)) {
            LUTu lhist16(32768);
            lhist16.clear();
#ifdef _OPENMP
            const int numThreads = min(max(pW * pH / (int)lhist16.getSize(), 1), omp_get_max_threads());
            #pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
            {
                LUTu lhist16thr(lhist16.getSize());
                lhist16thr.clear();
#ifdef _OPENMP
                #pragma omp for nowait
#endif

                for (int x = 0; x < pH; x++)
                    for (int y = 0; y < pW; y++) {
                        int pos = (int)(oprevl->L[x][y]);
                        lhist16thr[pos]++;
                    }

#ifdef _OPENMP
                #pragma omp critical
#endif
                lhist16 += lhist16thr;
            }
#ifdef _OPENMP
            static_cast<void>(numThreads);  // to silence cppcheck warning
#endif
            CurveFactory::complexLCurve(params->labCurve.brightness, params->labCurve.contrast, params->labCurve.lcurve, lhist16, lumacurve, histLCurve, scale == 1 ? 1 : 16, utili);
        }

        if (todo & M_LUMACURVE) {

            clcutili = CurveFactory::diagonalCurve2Lut(params->labCurve.clcurve, clcurve, scale == 1 ? 1 : 16);

            CurveFactory::complexsgnCurve(autili, butili, ccutili, cclutili, params->labCurve.acurve, params->labCurve.bcurve, params->labCurve.cccurve,
                                          params->labCurve.lccurve, chroma_acurve, chroma_bcurve, satcurve, lhskcurve, scale == 1 ? 1 : 16);
        }

        //scale = 1;

        if ((todo & (M_LUMINANCE + M_COLOR)) || (todo & M_AUTOEXP)) {
            nprevl->CopyFrom(oprevl);

            //  int maxspot = 1;
            //*************************************************************
            // locallab
            //*************************************************************

            if (params->locallab.enabled && !params->locallab.spots.empty()) {
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
                 *  2017 2018 Jacques Desmis <jdesmis@gmail.com>
                 *  2019 Pierre Cabrera <pierre.cab@gmail.com>
                 */
                const std::unique_ptr<LabImage> reserv(new LabImage(*oprevl, true));
                const std::unique_ptr<LabImage> lastorigimp(new LabImage(*oprevl, true));
                float **shbuffer = nullptr;
                int sca = 1;
                double huere, chromare, lumare, huerefblu, chromarefblu, lumarefblu, sobelre;
                float avge;
                std::vector<LocallabListener::locallabRef> locallref;
                std::vector<LocallabListener::locallabRetiMinMax> locallretiminmax;
                huerefs.resize(params->locallab.spots.size());
                huerefblurs.resize(params->locallab.spots.size());
                chromarefblurs.resize(params->locallab.spots.size());
                lumarefblurs.resize(params->locallab.spots.size());
                chromarefs.resize(params->locallab.spots.size());
                lumarefs.resize(params->locallab.spots.size());
                sobelrefs.resize(params->locallab.spots.size());
                avgs.resize(params->locallab.spots.size());

                for (int sp = 0; sp < (int)params->locallab.spots.size(); sp++) {
                    // Set local curves of current spot to LUT
                    locRETgainCurve.Set(params->locallab.spots.at(sp).localTgaincurve);
                    locRETtransCurve.Set(params->locallab.spots.at(sp).localTtranscurve);
                    const bool LHutili = loclhCurve.Set(params->locallab.spots.at(sp).LHcurve);
                    const bool HHutili = lochhCurve.Set(params->locallab.spots.at(sp).HHcurve);
                    const bool CHutili = locchCurve.Set(params->locallab.spots.at(sp).CHcurve);
                    const bool lcmasutili = locccmasCurve.Set(params->locallab.spots.at(sp).CCmaskcurve);
                    const bool llmasutili = locllmasCurve.Set(params->locallab.spots.at(sp).LLmaskcurve);
                    const bool lhmasutili = lochhmasCurve.Set(params->locallab.spots.at(sp).HHmaskcurve);
                    const bool lhhmasutili = lochhhmasCurve.Set(params->locallab.spots.at(sp).HHhmaskcurve);
                    const bool llmasexputili = locllmasexpCurve.Set(params->locallab.spots.at(sp).LLmaskexpcurve);
                    const bool lcmasexputili = locccmasexpCurve.Set(params->locallab.spots.at(sp).CCmaskexpcurve);
                    const bool lhmasexputili = lochhmasexpCurve.Set(params->locallab.spots.at(sp).HHmaskexpcurve);
                    const bool llmasSHutili = locllmasSHCurve.Set(params->locallab.spots.at(sp).LLmaskSHcurve);
                    const bool lcmasSHutili = locccmasSHCurve.Set(params->locallab.spots.at(sp).CCmaskSHcurve);
                    const bool lhmasSHutili = lochhmasSHCurve.Set(params->locallab.spots.at(sp).HHmaskSHcurve);
                    const bool llmasvibutili = locllmasvibCurve.Set(params->locallab.spots.at(sp).LLmaskvibcurve);
                    const bool lcmasvibutili = locccmasvibCurve.Set(params->locallab.spots.at(sp).CCmaskvibcurve);
                    const bool lhmasvibutili = lochhmasvibCurve.Set(params->locallab.spots.at(sp).HHmaskvibcurve);
                    const bool llmascbutili = locllmascbCurve.Set(params->locallab.spots.at(sp).LLmaskcbcurve);
                    const bool lcmascbutili = locccmascbCurve.Set(params->locallab.spots.at(sp).CCmaskcbcurve);
                    const bool lhmascbutili = lochhmascbCurve.Set(params->locallab.spots.at(sp).HHmaskcbcurve);
                    const bool llmaslcutili = locllmaslcCurve.Set(params->locallab.spots.at(sp).LLmasklccurve);
                    const bool lcmaslcutili = locccmaslcCurve.Set(params->locallab.spots.at(sp).CCmasklccurve);
                    const bool lhmaslcutili = lochhmaslcCurve.Set(params->locallab.spots.at(sp).HHmasklccurve);
                    const bool llmasretiutili = locllmasretiCurve.Set(params->locallab.spots.at(sp).LLmaskreticurve);
                    const bool lcmasretiutili = locccmasretiCurve.Set(params->locallab.spots.at(sp).CCmaskreticurve);
                    const bool lhmasretiutili = lochhmasretiCurve.Set(params->locallab.spots.at(sp).HHmaskreticurve);
                    const bool llmastmutili = locllmastmCurve.Set(params->locallab.spots.at(sp).LLmasktmcurve);
                    const bool lcmastmutili = locccmastmCurve.Set(params->locallab.spots.at(sp).CCmasktmcurve);
                    const bool lhmastmutili = lochhmastmCurve.Set(params->locallab.spots.at(sp).HHmasktmcurve);
                    const bool llmasblutili = locllmasblCurve.Set(params->locallab.spots.at(sp).LLmaskblcurve);
                    const bool lcmasblutili = locccmasblCurve.Set(params->locallab.spots.at(sp).CCmaskblcurve);
                    const bool lhmasblutili = lochhmasblCurve.Set(params->locallab.spots.at(sp).HHmaskblcurve);
                    const bool lcmas_utili = locccmas_Curve.Set(params->locallab.spots.at(sp).CCmask_curve);
                    const bool llmas_utili = locllmas_Curve.Set(params->locallab.spots.at(sp).LLmask_curve);
                    const bool lhmas_utili = lochhmas_Curve.Set(params->locallab.spots.at(sp).HHmask_curve);
                    const bool lhhmas_utili = lochhhmas_Curve.Set(params->locallab.spots.at(sp).HHhmask_curve);
                    const bool lmasutiliblwav = loclmasCurveblwav.Set(params->locallab.spots.at(sp).LLmaskblcurvewav);
                    const bool lmasutilicolwav = loclmasCurvecolwav.Set(params->locallab.spots.at(sp).LLmaskcolcurvewav);
                    const bool locwavutili = locwavCurve.Set(params->locallab.spots.at(sp).locwavcurve);
                    const bool loclevwavutili = loclevwavCurve.Set(params->locallab.spots.at(sp).loclevwavcurve);
                    const bool locconwavutili = locconwavCurve.Set(params->locallab.spots.at(sp).locconwavcurve);
                    const bool loccompwavutili = loccompwavCurve.Set(params->locallab.spots.at(sp).loccompwavcurve);
                    const bool loccomprewavutili = loccomprewavCurve.Set(params->locallab.spots.at(sp).loccomprewavcurve);
                    const bool locwavdenutili = locwavCurveden.Set(params->locallab.spots.at(sp).locwavcurveden);
                    const bool locedgwavutili = locedgwavCurve.Set(params->locallab.spots.at(sp).locedgwavcurve);
                    const bool lmasutili_wav = loclmasCurve_wav.Set(params->locallab.spots.at(sp).LLmask_curvewav);
                    const bool locallutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).llcurve, lllocalcurve, sca);
                    const bool localclutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).clcurve, cllocalcurve, sca);
                    const bool locallcutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).lccurve, lclocalcurve, sca);
                    const bool localcutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).cccurve, cclocalcurve, sca);
                    const bool localrgbutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).rgbcurve, rgblocalcurve, sca);
                    const bool localexutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).excurve, exlocalcurve, sca);
                    const bool localmaskutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmaskcurve, lmasklocalcurve, sca);
                    const bool localmaskexputili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmaskexpcurve, lmaskexplocalcurve, sca);
                    const bool localmaskSHutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).LmaskSHcurve, lmaskSHlocalcurve, sca);
                    const bool localmaskvibutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmaskvibcurve, lmaskviblocalcurve, sca);
                    const bool localmasktmutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmasktmcurve, lmasktmlocalcurve, sca);
                    const bool localmaskretiutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmaskreticurve, lmaskretilocalcurve, sca);
                    const bool localmaskcbutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmaskcbcurve, lmaskcblocalcurve, sca);
                    const bool localmaskblutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmaskblcurve, lmaskbllocalcurve, sca);
                    const bool localmasklcutili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmasklccurve, lmasklclocalcurve, sca);
                    const bool localmask_utili = CurveFactory::diagonalCurve2Lut(params->locallab.spots.at(sp).Lmask_curve, lmasklocal_curve, sca);
                    double ecomp = params->locallab.spots.at(sp).expcomp;
                    double black = params->locallab.spots.at(sp).black;
                    double hlcompr = params->locallab.spots.at(sp).hlcompr;
                    double hlcomprthresh = params->locallab.spots.at(sp).hlcomprthresh;
                    double shcompr = params->locallab.spots.at(sp).shcompr;
                    double br = params->locallab.spots.at(sp).lightness;
                    double cont = params->locallab.spots.at(sp).contrast;

                    if (black < 0. && params->locallab.spots.at(sp).expMethod == "pde") {
                        black *= 1.5;
                    }

                    // Reference parameters computation
                    if (params->locallab.spots.at(sp).spotMethod == "exc") {
                        ipf.calc_ref(sp, reserv.get(), reserv.get(), 0, 0, pW, pH, scale, huerefblu, chromarefblu, lumarefblu, huere, chromare, lumare, sobelre, avge, locwavCurveden, locwavdenutili);
                    } else {
                        ipf.calc_ref(sp, nprevl, nprevl, 0, 0, pW, pH, scale, huerefblu, chromarefblu, lumarefblu, huere, chromare, lumare, sobelre, avge, locwavCurveden, locwavdenutili);
                    }

                    double huerblu = huerefblurs[sp] = huerefblu;
                    double chromarblu = chromarefblurs[sp] = chromarefblu;
                    double lumarblu = lumarefblurs[sp] = lumarefblu;
                    double huer = huerefs[sp] = huere;
                    double chromar = chromarefs[sp] = chromare;
                    double lumar = lumarefs[sp] = lumare ;
                    double sobeler = sobelrefs[sp] = sobelre;
                    float avg = avgs[sp] = avge;
                    CurveFactory::complexCurvelocal(ecomp, black / 65535., hlcompr, hlcomprthresh, shcompr, br, cont, lumar,
                                                    hltonecurveloc, shtonecurveloc, tonecurveloc, lightCurveloc, avg,
                                                    sca);

                    // Save Locallab mask curve references for current spot
                    LocallabListener::locallabRef spotref;
                    spotref.huer = huer;
                    spotref.lumar = lumar;
                    spotref.chromar = chromar;
                    locallref.push_back(spotref);

                    // Locallab tools computation
                    /* Notes:
                     * - shbuffer is used as nullptr
                     */
                    // Locallab mask is only showed in detailed image
                    float minCD;
                    float maxCD;
                    float mini;
                    float maxi;
                    float Tmean;
                    float Tsigma;
                    float Tmin;
                    float Tmax;
                    int lastsav;
                    ipf.Lab_Local(3, sp, (float**)shbuffer, nprevl, nprevl, reserv.get(), lastorigimp.get(), 0, 0, pW, pH, scale, locRETgainCurve, locRETtransCurve,
                                  lllocalcurve, locallutili,
                                  cllocalcurve, localclutili,
                                  lclocalcurve, locallcutili,
                                  loclhCurve,  lochhCurve, locchCurve,
                                  lmasklocalcurve, localmaskutili,
                                  lmaskexplocalcurve, localmaskexputili,
                                  lmaskSHlocalcurve, localmaskSHutili,
                                  lmaskviblocalcurve, localmaskvibutili,
                                  lmasktmlocalcurve, localmasktmutili,
                                  lmaskretilocalcurve, localmaskretiutili,
                                  lmaskcblocalcurve, localmaskcbutili,
                                  lmaskbllocalcurve, localmaskblutili,
                                  lmasklclocalcurve, localmasklcutili,
                                  lmasklocal_curve, localmask_utili,

                                  locccmasCurve, lcmasutili, locllmasCurve, llmasutili, lochhmasCurve, lhmasutili, lochhhmasCurve, lhhmasutili, locccmasexpCurve, lcmasexputili, locllmasexpCurve, llmasexputili, lochhmasexpCurve, lhmasexputili,
                                  locccmasSHCurve, lcmasSHutili, locllmasSHCurve, llmasSHutili, lochhmasSHCurve, lhmasSHutili,
                                  locccmasvibCurve, lcmasvibutili, locllmasvibCurve, llmasvibutili, lochhmasvibCurve, lhmasvibutili,
                                  locccmascbCurve, lcmascbutili, locllmascbCurve, llmascbutili, lochhmascbCurve, lhmascbutili,
                                  locccmasretiCurve, lcmasretiutili, locllmasretiCurve, llmasretiutili, lochhmasretiCurve, lhmasretiutili,
                                  locccmastmCurve, lcmastmutili, locllmastmCurve, llmastmutili, lochhmastmCurve, lhmastmutili,
                                  locccmasblCurve, lcmasblutili, locllmasblCurve, llmasblutili, lochhmasblCurve, lhmasblutili,
                                  locccmaslcCurve, lcmaslcutili, locllmaslcCurve, llmaslcutili, lochhmaslcCurve, lhmaslcutili,
                                  locccmas_Curve, lcmas_utili, locllmas_Curve, llmas_utili, lochhmas_Curve, lhmas_utili,
                                  lochhhmas_Curve, lhhmas_utili,
                                  loclmasCurveblwav, lmasutiliblwav,
                                  loclmasCurvecolwav, lmasutilicolwav,
                                  locwavCurve, locwavutili,
                                  loclevwavCurve, loclevwavutili,
                                  locconwavCurve, locconwavutili,
                                  loccompwavCurve, loccompwavutili,
                                  loccomprewavCurve, loccomprewavutili,
                                  locwavCurveden, locwavdenutili,
                                  locedgwavCurve, locedgwavutili,
                                  loclmasCurve_wav, lmasutili_wav,
                                  LHutili, HHutili, CHutili, cclocalcurve, localcutili, rgblocalcurve, localrgbutili, localexutili, exlocalcurve, hltonecurveloc, shtonecurveloc, tonecurveloc, lightCurveloc,
                                  huerblu, chromarblu, lumarblu, huer, chromar, lumar, sobeler, lastsav, false, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);

                    if (sp + 1u < params->locallab.spots.size()) {
                        // do not copy for last spot as it is not needed anymore
                        lastorigimp->CopyFrom(nprevl);
                    }

                    // Save Locallab Retinex min/max for current spot
                    LocallabListener::locallabRetiMinMax retiMinMax;
                    retiMinMax.cdma = maxCD;
                    retiMinMax.cdmin = minCD;
                    retiMinMax.mini = mini;
                    retiMinMax.maxi = maxi;
                    retiMinMax.Tmean = Tmean;
                    retiMinMax.Tsigma = Tsigma;
                    retiMinMax.Tmin = Tmin;
                    retiMinMax.Tmax = Tmax;
                    locallretiminmax.push_back(retiMinMax);

                    // Recalculate references after
                    if (params->locallab.spots.at(sp).spotMethod == "exc") {
                        ipf.calc_ref(sp, reserv.get(), reserv.get(), 0, 0, pW, pH, scale, huerefblu, chromarefblu, lumarefblu, huer, chromar, lumar, sobeler, avg, locwavCurveden, locwavdenutili);
                    } else {
                        ipf.calc_ref(sp, nprevl, nprevl, 0, 0, pW, pH, scale, huerefblu, chromarefblu, lumarefblu, huer, chromar, lumar, sobeler, avg, locwavCurveden, locwavdenutili);
                    }

                    // Update Locallab reference values according to recurs parameter
                    if (params->locallab.spots.at(sp).recurs) {
                        locallref.at(sp).chromar = chromar;
                        locallref.at(sp).lumar = lumar;
                        locallref.at(sp).huer = huer;
                    }
                }

                // Transmit Locallab reference values and Locallab Retinex min/max to LocallabListener
                if (locallListener) {
                    locallListener->refChanged(locallref, params->locallab.selspot);
                    locallListener->minmaxChanged(locallretiminmax, params->locallab.selspot);
                }
            }

            //*************************************************************
            // end locallab
            //*************************************************************

            histCCurve.clear();
            histLCurve.clear();
            ipf.chromiLuminanceCurve(nullptr, pW, nprevl, nprevl, chroma_acurve, chroma_bcurve, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, histCCurve, histLCurve);
            ipf.vibrance(nprevl, params->vibrance, params->toneCurve.hrenabled, params->icm.workingProfile);
            ipf.labColorCorrectionRegions(nprevl);

            if ((params->colorappearance.enabled && !params->colorappearance.tonecie) || (!params->colorappearance.enabled)) {
                ipf.EPDToneMap(nprevl, 0, scale);
            }

            if (params->dirpyrequalizer.cbdlMethod == "aft") {
                if (((params->colorappearance.enabled && !settings->autocielab) || (!params->colorappearance.enabled))) {
                    ipf.dirpyrequalizer(nprevl, scale);
                }
            }

            wavcontlutili = CurveFactory::diagonalCurve2Lut(params->wavelet.wavclCurve, wavclCurve, scale == 1 ? 1 : 16);

            if ((params->wavelet.enabled)) {
                WaveletParams WaveParams = params->wavelet;
                WaveParams.getCurves(wavCLVCurve, wavblcurve, waOpacityCurveRG, waOpacityCurveSH, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL);
                int kall = 0;
                LabImage *unshar = nullptr;
                Glib::ustring provis;
                LabImage *provradius = nullptr;
                bool procont = WaveParams.expcontrast;
                bool prochro = WaveParams.expchroma;
                bool proedge = WaveParams.expedge;
                bool profin = WaveParams.expfinal;
                bool proton = WaveParams.exptoning;
                bool pronois = WaveParams.expnoise; 

                if(WaveParams.showmask) {
                 //   WaveParams.showmask = false;
                 //   WaveParams.expclari = true;
                }

                if (WaveParams.softrad > 0.f) {
                    provradius = new LabImage(*nprevl, true);
                }

                if ((WaveParams.ushamethod == "sharp" || WaveParams.ushamethod == "clari") && WaveParams.expclari && WaveParams.CLmethod != "all") {
                    provis = params->wavelet.CLmethod;
                    params->wavelet.CLmethod = "all";
                    ipf.ip_wavelet(nprevl, nprevl, kall, WaveParams, wavCLVCurve, wavblcurve, waOpacityCurveRG, waOpacityCurveSH, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL, wavclCurve, scale);
                    unshar = new LabImage(*nprevl, true);

                    params->wavelet.CLmethod = provis;

                    WaveParams.expcontrast = false;
                    WaveParams.expchroma = false;
                    WaveParams.expedge = false;
                    WaveParams.expfinal = false;
                    WaveParams.exptoning = false;
                    WaveParams.expnoise = false; 
                }

                ipf.ip_wavelet(nprevl, nprevl, kall, WaveParams, wavCLVCurve, wavblcurve, waOpacityCurveRG, waOpacityCurveSH, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL, wavclCurve, scale);


                if ((WaveParams.ushamethod == "sharp" || WaveParams.ushamethod == "clari") && WaveParams.expclari && WaveParams.CLmethod != "all") {
                    WaveParams.expcontrast = procont;
                    WaveParams.expchroma = prochro;
                    WaveParams.expedge = proedge;
                    WaveParams.expfinal = profin;
                    WaveParams.exptoning = proton;
                    WaveParams.expnoise = pronois;
                    
                    if (WaveParams.softrad > 0.f) {

                        array2D<float> ble(pW, pH);
                        array2D<float> guid(pW, pH);
                        Imagefloat *tmpImage = nullptr;
                        tmpImage = new Imagefloat(pW, pH);

#ifdef _OPENMP
                        #pragma omp parallel for
#endif

                        for (int ir = 0; ir < pH; ir++)
                            for (int jr = 0; jr < pW; jr++) {
                                float X, Y, Z;
                                float L = provradius->L[ir][jr];
                                float a = provradius->a[ir][jr];
                                float b = provradius->b[ir][jr];
                                Color::Lab2XYZ(L, a, b, X, Y, Z);

                                guid[ir][jr] = Y / 32768.f;
                                float La = nprevl->L[ir][jr];
                                float aa = nprevl->a[ir][jr];
                                float ba = nprevl->b[ir][jr];
                                Color::Lab2XYZ(La, aa, ba, X, Y, Z);
                                tmpImage->r(ir, jr) = X;
                                tmpImage->g(ir, jr) = Y;
                                tmpImage->b(ir, jr) = Z;
                                ble[ir][jr] = Y / 32768.f;
                            }
    
                        double epsilmax = 0.0001;
                        double epsilmin = 0.00001;
                        double aepsil = (epsilmax - epsilmin) / 100.f;
                        double bepsil = epsilmin; //epsilmax - 100.f * aepsil;
                        double epsil = aepsil * WaveParams.softrad + bepsil;

                        float blur = 10.f / scale * (0.5f + 0.8f * WaveParams.softrad);
                        // rtengine::guidedFilter(guid, ble, ble, blur, 0.001, multiTh);
                        rtengine::guidedFilter(guid, ble, ble, blur, epsil, false);



#ifdef _OPENMP
                        #pragma omp parallel for
#endif

                        for (int ir = 0; ir < pH; ir++)
                            for (int jr = 0; jr < pW; jr++) {
                                float X = tmpImage->r(ir, jr);
                                float Y = 32768.f * ble[ir][jr];
                                float Z = tmpImage->b(ir, jr);
                                float L, a, b;
                                Color::XYZ2Lab(X, Y, Z, L, a, b);
                                nprevl->L[ir][jr] =  L;
                            }
      
                    delete tmpImage;

                    }
                    
                }

                if ((WaveParams.ushamethod == "sharp" || WaveParams.ushamethod == "clari")  && WaveParams.expclari && WaveParams.CLmethod != "all") {
                    float mL = (float)(WaveParams.mergeL / 100.f);
                    float mC = (float)(WaveParams.mergeC / 100.f);
                    float mL0;
                    float mC0;
                    float background = 0.f;
                    int show = 0; 



                    if ((WaveParams.CLmethod == "one" || WaveParams.CLmethod == "inf")  && WaveParams.Backmethod == "black") {
                        mL0 = mC0 = 0.f;
                        mL = - 1.5f * mL;
                        mC = -mC;
                        background = 12000.f;
                        show = 0;
                    } else if (WaveParams.CLmethod == "sup" && WaveParams.Backmethod == "resid") {
                        mL0 = mL;
                        mC0 = mC;
                        background = 0.f;
                        show = 0;
                    } else {
                        mL0 = mL = mC0 = mC = 0.f;
                        background = 0.f;
                        show = 0;
                    }
                float indic = 1.f;

                if(WaveParams.showmask){
                    mL0 = mC0 = -1.f;
                    indic = -1.f;
                    mL = fabs(mL);
                    mC = fabs(mC);
                    show = 1;
                }
#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int x = 0; x < pH; x++)
                        for (int y = 0; y < pW; y++) {
                            nprevl->L[x][y] = LIM((1.f + mL0) * (unshar->L[x][y]) + show * background - mL * indic * nprevl->L[x][y], 0.f, 32768.f);
                            nprevl->a[x][y] = (1.f + mC0) * (unshar->a[x][y]) - mC * indic * nprevl->a[x][y];
                            nprevl->b[x][y] = (1.f + mC0) * (unshar->b[x][y]) - mC * indic * nprevl->b[x][y];
                        }

                    delete unshar;
                    unshar    = NULL;
/*
                    if (WaveParams.softrad > 0.f) {
                        array2D<float> ble(pW, pH);
                        array2D<float> guid(pW, pH);
#ifdef _OPENMP
                        #pragma omp parallel for
#endif

                        for (int ir = 0; ir < pH; ir++)
                            for (int jr = 0; jr < pW; jr++) {
                                ble[ir][jr] = (nprevl->L[ir][jr]  - provradius->L[ir][jr]) / 32768.f;
                                guid[ir][jr] = provradius->L[ir][jr] / 32768.f;
                            }
                        double epsilmax = 0.001;
                        double epsilmin = 0.0001;
                        double aepsil = (epsilmax - epsilmin) / 90.f;
                        double bepsil = epsilmax - 100.f * aepsil;
                        double epsil = aepsil * WaveParams.softrad + bepsil;

                        float blur = 10.f / scale * (0.001f + 0.8f * WaveParams.softrad);
                        // rtengine::guidedFilter(guid, ble, ble, blur, 0.001, multiTh);
                        rtengine::guidedFilter(guid, ble, ble, blur, epsil, false);



#ifdef _OPENMP
                        #pragma omp parallel for
#endif

                        for (int ir = 0; ir < pH; ir++)
                            for (int jr = 0; jr < pW; jr++) {
                                nprevl->L[ir][jr] =  provradius->L[ir][jr] + 32768.f * ble[ir][jr];
                            }
                    }
*/
                    if (WaveParams.softrad > 0.f) {

                        delete provradius;
                        provradius    = NULL;

                    }


                }
               
            }

            ipf.softLight(nprevl, params->softlight);

            if (params->colorappearance.enabled) {
                // L histo  and Chroma histo for ciecam
                // histogram well be for Lab (Lch) values, because very difficult to do with J,Q, M, s, C
                int x1, y1, x2, y2;
                params->crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
                lhist16CAM.clear();
                lhist16CCAM.clear();

                if (!params->colorappearance.datacie) {
                    for (int x = 0; x < pH; x++)
                        for (int y = 0; y < pW; y++) {
                            int pos = CLIP((int)(nprevl->L[x][y]));
                            int posc = CLIP((int)sqrt(nprevl->a[x][y] * nprevl->a[x][y] + nprevl->b[x][y] * nprevl->b[x][y]));
                            lhist16CAM[pos]++;
                            lhist16CCAM[posc]++;
                        }
                }

                CurveFactory::curveLightBrightColor(params->colorappearance.curve, params->colorappearance.curve2, params->colorappearance.curve3,
                                                    lhist16CAM, histLCAM, lhist16CCAM, histCCAM,
                                                    customColCurve1, customColCurve2, customColCurve3, 1);

                const FramesMetaData* metaData = imgsrc->getMetaData();
                int imgNum = 0;

                if (imgsrc->isRAW()) {
                    if (imgsrc->getSensorType() == ST_BAYER) {
                        imgNum = rtengine::LIM<unsigned int>(params->raw.bayersensor.imageNum, 0, metaData->getFrameCount() - 1);
                    } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                        //imgNum = rtengine::LIM<unsigned int>(params->raw.xtranssensor.imageNum, 0, metaData->getFrameCount() - 1);
                    }
                }

                float fnum = metaData->getFNumber(imgNum);          // F number
                float fiso = metaData->getISOSpeed(imgNum) ;        // ISO
                float fspeed = metaData->getShutterSpeed(imgNum) ;  // Speed
                double fcomp = metaData->getExpComp(imgNum);        // Compensation +/-
                double adap;

                if (fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) { //if no exif data or wrong
                    adap = 2000.;
                } else {
                    double E_V = fcomp + log2(double ((fnum * fnum) / fspeed / (fiso / 100.f)));
                    E_V += params->toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                    E_V += log2(params->raw.expos);  // exposure raw white point ; log2 ==> linear to EV
                    adap = pow(2.0, E_V - 3.0);  // cd / m2
                    // end calculation adaptation scene luminosity
                }

                float d, dj, yb;
                bool execsharp = false;

                if (!ncie) {
                    ncie = new CieImage(pW, pH);
                }

                if (!CAMBrightCurveJ && (params->colorappearance.algo == "JC" || params->colorappearance.algo == "JS" || params->colorappearance.algo == "ALL")) {
                    CAMBrightCurveJ(32768, 0);
                }

                if (!CAMBrightCurveQ && (params->colorappearance.algo == "QM" || params->colorappearance.algo == "ALL")) {
                    CAMBrightCurveQ(32768, 0);
                }

                // Issue 2785, only float version of ciecam02 for navigator and pan background
                CAMMean = NAN;
                CAMBrightCurveJ.dirty = true;
                CAMBrightCurveQ.dirty = true;

                ipf.ciecam_02float(ncie, float (adap), pW, 2, nprevl, params.get(), customColCurve1, customColCurve2, customColCurve3, histLCAM, histCCAM, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 0, scale, execsharp, d, dj, yb, 1);

                if ((params->colorappearance.autodegree || params->colorappearance.autodegreeout) && acListener && params->colorappearance.enabled && !params->colorappearance.presetcat02) {
                    acListener->autoCamChanged(100.* (double)d, 100.* (double)dj);
                }

                if (params->colorappearance.autoadapscen && acListener && params->colorappearance.enabled && !params->colorappearance.presetcat02) {
                    acListener->adapCamChanged(adap);    //real value of adapt scene
                }

                if (params->colorappearance.autoybscen && acListener && params->colorappearance.enabled && !params->colorappearance.presetcat02) {
                    acListener->ybCamChanged((int) yb);    //real value Yb scene
                }

                if (params->colorappearance.enabled && params->colorappearance.presetcat02  && params->colorappearance.autotempout) {
              //      acListener->wbCamChanged(params->wb.temperature, params->wb.green);    //real temp and tint
                    acListener->wbCamChanged(params->wb.temperature, 1.f);    //real temp and tint = 1.
                }
                
            } else {
                // CIECAM is disabled, we free up its image buffer to save some space
                if (ncie) {
                    delete ncie;
                }

                ncie = nullptr;

                if (CAMBrightCurveJ) {
                    CAMBrightCurveJ.reset();
                }

                if (CAMBrightCurveQ) {
                    CAMBrightCurveQ.reset();
                }
            }
        }

        // Update the monitor color transform if necessary
        if ((todo & M_MONITOR) || (lastOutputProfile != params->icm.outputProfile) || lastOutputIntent != params->icm.outputIntent || lastOutputBPC != params->icm.outputBPC) {
            lastOutputProfile = params->icm.outputProfile;
            lastOutputIntent = params->icm.outputIntent;
            lastOutputBPC = params->icm.outputBPC;
            ipf.updateColorProfiles(monitorProfile, monitorIntent, softProof, gamutCheck);
        }
    }

// process crop, if needed
    for (size_t i = 0; i < crops.size(); i++)
        if (crops[i]->hasListener() && (panningRelatedChange || (highDetailNeeded && options.prevdemo != PD_Sidecar) || (todo & (M_MONITOR | M_RGBCURVE | M_LUMACURVE)) || crops[i]->get_skip() == 1)) {
            crops[i]->update(todo);     // may call ourselves
        }

    if (panningRelatedChange || (todo & M_MONITOR)) {
        if ((todo != CROP && todo != MINUPDATE) || (todo & M_MONITOR)) {
            MyMutex::MyLock prevImgLock(previmg->getMutex());

            try {
                // Computing the preview image, i.e. converting from WCS->Monitor color space (soft-proofing disabled) or WCS->Printer profile->Monitor color space (soft-proofing enabled)
                ipf.lab2monitorRgb(nprevl, previmg);

                // Computing the internal image for analysis, i.e. conversion from WCS->Output profile
                delete workimg;
                workimg = ipf.lab2rgb(nprevl, 0, 0, pW, pH, params->icm);
            } catch (char * str) {
                return;
            }
        }

        if (!resultValid) {
            resultValid = true;

            if (imageListener) {
                imageListener->setImage(previmg, scale, params->crop);
            }
        }

        if (imageListener)
            // TODO: The WB tool should be advertised too in order to get the AutoWB's temp and green values
        {
            imageListener->imageReady(params->crop);
        }

        if (hListener) {
            updateLRGBHistograms();
            hListener->histogramChanged(histRed, histGreen, histBlue, histLuma, histToneCurve, histLCurve, histCCurve, /*histCLurve, histLLCurve,*/ histLCAM, histCCAM, histRedRaw, histGreenRaw, histBlueRaw, histChroma, histLRETI);
        }
    }

    if (orig_prev != oprevi) {
        delete oprevi;
        oprevi = nullptr;
    }


}


void ImProcCoordinator::freeAll()
{

    if (allocated) {
        if (orig_prev != oprevi) {
            delete oprevi;
        }

        oprevi    = nullptr;
        delete orig_prev;
        orig_prev = nullptr;
        delete oprevl;
        oprevl    = nullptr;
        delete nprevl;
        nprevl    = nullptr;

        if (ncie) {
            delete ncie;
        }

        ncie      = nullptr;

        if (imageListener) {
            imageListener->delImage(previmg);
        } else {
            delete previmg;
        }

        delete workimg;

    }

    allocated = false;
}

/** @brief Handles image buffer (re)allocation and trigger sizeChanged of SizeListener[s]
 * If the scale change, this method will free all buffers and reallocate ones of the new size.
 * It will then tell to the SizeListener that size has changed (sizeChanged)
 *
 * @param prevscale New Preview's scale.
 */
void ImProcCoordinator::setScale(int prevscale)
{

    tr = getCoarseBitMask(params->coarse);

    int nW, nH;
    imgsrc->getFullSize(fw, fh, tr);

    prevscale++;

    do {
        prevscale--;
        PreviewProps pp(0, 0, fw, fh, prevscale);
        imgsrc->getSize(pp, nW, nH);
    } while (nH < 400 && prevscale > 1 && (nW * nH < 1000000));  // sctually hardcoded values, perhaps a better choice is possible

    if (nW != pW || nH != pH) {

        freeAll();

        pW = nW;
        pH = nH;

        orig_prev = new Imagefloat(pW, pH);
        oprevi = orig_prev;
        oprevl = new LabImage(pW, pH);
        nprevl = new LabImage(pW, pH);

        //ncie is only used in ImProcCoordinator::updatePreviewImage, it will be allocated on first use and deleted if not used anymore
        previmg = new Image8(pW, pH);
        workimg = new Image8(pW, pH);

        allocated = true;
    }

    scale = prevscale;
    resultValid = false;
    fullw = fw;
    fullh = fh;

    if (!sizeListeners.empty())
        for (size_t i = 0; i < sizeListeners.size(); i++) {
            sizeListeners[i]->sizeChanged(fullw, fullh, fw, fh);
        }
}


void ImProcCoordinator::updateLRGBHistograms()
{

    int x1, y1, x2, y2;
    params->crop.mapToResized(pW, pH, scale, x1, x2, y1, y2);

#ifdef _OPENMP
    #pragma omp parallel sections
#endif
    {
#ifdef _OPENMP
        #pragma omp section
#endif
        {
            histChroma.clear();

            for (int i = y1; i < y2; i++)
                for (int j = x1; j < x2; j++)
                {
                    histChroma[(int)(sqrtf(SQR(nprevl->a[i][j]) + SQR(nprevl->b[i][j])) / 188.f)]++;      //188 = 48000/256
                }
        }
#ifdef _OPENMP
        #pragma omp section
#endif
        {
            histLuma.clear();

            for (int i = y1; i < y2; i++)
                for (int j = x1; j < x2; j++)
                {
                    histLuma[(int)(nprevl->L[i][j] / 128.f)]++;
                }
        }
#ifdef _OPENMP
        #pragma omp section
#endif
        {
            histRed.clear();
            histGreen.clear();
            histBlue.clear();

            for (int i = y1; i < y2; i++)
            {
                int ofs = (i * pW + x1) * 3;

                for (int j = x1; j < x2; j++) {
                    int r = workimg->data[ofs++];
                    int g = workimg->data[ofs++];
                    int b = workimg->data[ofs++];

                    histRed[r]++;
                    histGreen[g]++;
                    histBlue[b]++;
                }
            }
        }
    }

}

bool ImProcCoordinator::getAutoWB(double& temp, double& green, double equal, double tempBias)
{

    if (imgsrc) {
        if (lastAwbEqual != equal || lastAwbTempBias != tempBias || lastAwbauto != params->wb.method) {
// Issue 2500            MyMutex::MyLock lock(minit);  // Also used in crop window
            double rm, gm, bm;
            params->wb.method = "autold";//same result as before muliple Auto WB
            
           // imgsrc->getAutoWBMultipliers(rm, gm, bm);
            double tempitc = 5000.;
            double greenitc = 1.;
            float studgood = 1000.f;
            double tempref, greenref;
            imgsrc->getAutoWBMultipliersitc(tempref, greenref, tempitc, greenitc, studgood,  0, 0, fh, fw, 0, 0, fh, fw, rm, gm, bm,  params->wb, params->icm, params->raw);

            if (rm != -1) {
                autoWB.update(rm, gm, bm, equal, tempBias);
                lastAwbEqual = equal;
                lastAwbTempBias = tempBias;
                lastAwbauto = params->wb.method;
            } else {
                lastAwbEqual = -1.;
                autoWB.useDefaults(equal);
                lastAwbauto = "";
                lastAwbTempBias = 0.0;
            }
        }

        temp = autoWB.getTemp();
        green = autoWB.getGreen();
        return true;
    } else {
        //temp = autoWB.getTemp();
        temp = -1.0;
        green = -1.0;
        return false;
    }
}

void ImProcCoordinator::getCamWB(double& temp, double& green)
{

    if (imgsrc) {
        temp = imgsrc->getWB().getTemp();
        green = imgsrc->getWB().getGreen();
    }
}

void ImProcCoordinator::getSpotWB(int x, int y, int rect, double& temp, double& tgreen)
{

    ColorTemp ret;

    {
        MyMutex::MyLock lock(mProcessing);
        std::vector<Coord2D> points, red, green, blue;

        for (int i = y - rect; i <= y + rect; i++)
            for (int j = x - rect; j <= x + rect; j++) {
                points.push_back(Coord2D(j, i));
            }

        ipf.transCoord(fw, fh, points, red, green, blue);

        int tr = getCoarseBitMask(params->coarse);

        ret = imgsrc->getSpotWB(red, green, blue, tr, params->wb.equal);
        currWB = ColorTemp(params->wb.temperature, params->wb.green, params->wb.equal, params->wb.method);
        //double rr,gg,bb;
        //currWB.getMultipliers(rr,gg,bb);

    } // end of mutex lockong

    if (ret.getTemp() > 0) {
        temp = ret.getTemp();
        tgreen = ret.getGreen();
    } else {
        temp = currWB.getTemp();
        tgreen = currWB.getGreen();
    }
}

bool ImProcCoordinator::getFilmNegativeExponents(int xA, int yA, int xB, int yB, std::array<float, 3>& newExps)
{
    MyMutex::MyLock lock(mProcessing);

    const int tr = getCoarseBitMask(params->coarse);

    const Coord2D p1 = translateCoord(ipf, fw, fh, xA, yA);
    const Coord2D p2 = translateCoord(ipf, fw, fh, xB, yB);

    return imgsrc->getFilmNegativeExponents(p1, p2, tr, params->filmNegative, newExps);
}

bool ImProcCoordinator::getRawSpotValues(int x, int y, int spotSize, std::array<float, 3>& rawValues)
{
    MyMutex::MyLock lock(mProcessing);

    return imgsrc->getRawSpotValues(translateCoord(ipf, fw, fh, x, y), spotSize,
        getCoarseBitMask(params->coarse), params->filmNegative, rawValues);
}

void ImProcCoordinator::getAutoCrop(double ratio, int &x, int &y, int &w, int &h)
{

    MyMutex::MyLock lock(mProcessing);

    LensCorrection *pLCPMap = nullptr;

    if (params->lensProf.useLcp() && imgsrc->getMetaData()->getFocalLen() > 0) {
        const std::shared_ptr<LCPProfile> pLCPProf = LCPStore::getInstance()->getProfile(params->lensProf.lcpFile);

        if (pLCPProf) pLCPMap = new LCPMapper(pLCPProf, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(), imgsrc->getMetaData()->getFocusDist(),
                                                  0, false, params->lensProf.useDist, fullw, fullh, params->coarse, imgsrc->getRotateDegree());
    }

    double fillscale = ipf.getTransformAutoFill(fullw, fullh, pLCPMap);

    if (ratio > 0) {
        w = fullw * fillscale;
        h = w / ratio;

        if (h > fullh * fillscale) {
            h = fullh * fillscale;
            w = h * ratio;
        }
    } else {
        w = fullw * fillscale;
        h = fullh * fillscale;
    }

    x = (fullw - w) / 2;
    y = (fullh - h) / 2;
}

void ImProcCoordinator::setMonitorProfile(const Glib::ustring& profile, RenderingIntent intent)
{
    monitorProfile = profile;
    monitorIntent = intent;
}

void ImProcCoordinator::getMonitorProfile(Glib::ustring& profile, RenderingIntent& intent) const
{
    profile = monitorProfile;
    intent = monitorIntent;
}

void ImProcCoordinator::setSoftProofing(bool softProof, bool gamutCheck)
{
    this->softProof = softProof;
    this->gamutCheck = gamutCheck;
}

void ImProcCoordinator::getSoftProofing(bool &softProof, bool &gamutCheck)
{
    softProof = this->softProof;
    gamutCheck = this->gamutCheck;
}

ProcEvent ImProcCoordinator::setSharpMask(bool sharpMask)
{
    if (this->sharpMask != sharpMask) {
        sharpMaskChanged = true;
        this->sharpMask = sharpMask;
        return params->pdsharpening.enabled ? rtengine::EvPdShrMaskToggled : rtengine::EvShrEnabled;
    } else {
        sharpMaskChanged = false;
        return rtengine::EvShrEnabled;
    }
}

void ImProcCoordinator::saveInputICCReference(const Glib::ustring& fname, bool apply_wb)
{

    MyMutex::MyLock lock(mProcessing);

    int fW, fH;

    int tr = getCoarseBitMask(params->coarse);

    imgsrc->getFullSize(fW, fH, tr);
    PreviewProps pp(0, 0, fW, fH, 1);
    ProcParams ppar = *params;
    ppar.toneCurve.hrenabled = false;
    ppar.icm.inputProfile = "(none)";
    Imagefloat* im = new Imagefloat(fW, fH);
    imgsrc->preprocess(ppar.raw, ppar.lensProf, ppar.coarse);
    double dummy = 0.0;
    imgsrc->demosaic(ppar.raw, false, dummy);
    ColorTemp currWB = ColorTemp(params->wb.temperature, params->wb.green, params->wb.equal, params->wb.method);

    if (params->wb.method == "Camera") {
        currWB = imgsrc->getWB();
    } else if (params->wb.method == "autold") {
        if (lastAwbEqual != params->wb.equal || lastAwbTempBias != params->wb.tempBias) {
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers(rm, gm, bm);

            if (rm != -1.) {
                autoWB.update(rm, gm, bm, params->wb.equal, params->wb.tempBias);
                lastAwbEqual = params->wb.equal;
                lastAwbTempBias = params->wb.tempBias;
            } else {
                lastAwbEqual = -1.;
                lastAwbTempBias = 0.0;
                autoWB.useDefaults(params->wb.equal);
            }
        }

        currWB = autoWB;
    }

    if (!apply_wb) {
        currWB = ColorTemp(); // = no white balance
    }

    imgsrc->getImage(currWB, tr, im, pp, ppar.toneCurve, ppar.raw);
    ImProcFunctions ipf(&ppar, true);

    if (ipf.needsTransform(fW, fH, imgsrc->getRotateDegree(), imgsrc->getMetaData())) {
        Imagefloat* trImg = new Imagefloat(fW, fH);
        ipf.transform(im, trImg, 0, 0, 0, 0, fW, fH, fW, fH,
                      imgsrc->getMetaData(), imgsrc->getRotateDegree(), true);
        delete im;
        im = trImg;
    }

    if (params->crop.enabled) {
        Imagefloat *tmpim = new Imagefloat(params->crop.w, params->crop.h);
        int cx = params->crop.x;
        int cy = params->crop.y;
        int cw = params->crop.w;
        int ch = params->crop.h;
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = cy; i < cy + ch; i++) {
            for (int j = cx; j < cx + cw; j++) {
                tmpim->r(i - cy, j - cx) = im->r(i, j);
                tmpim->g(i - cy, j - cx) = im->g(i, j);
                tmpim->b(i - cy, j - cx) = im->b(i, j);
            }
        }

        delete im;
        im = tmpim;
    }

    // image may contain out of range samples, clip them to avoid wrap-arounds
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < im->getHeight(); i++) {
        for (int j = 0; j < im->getWidth(); j++) {
            im->r(i, j) = CLIP(im->r(i, j));
            im->g(i, j) = CLIP(im->g(i, j));
            im->b(i, j) = CLIP(im->b(i, j));
        }
    }

    int imw, imh;
    double tmpScale = ipf.resizeScale(params.get(), fW, fH, imw, imh);

    if (tmpScale != 1.0) {
        Imagefloat* tempImage = new Imagefloat(imw, imh);
        ipf.resize(im, tempImage, tmpScale);
        delete im;
        im = tempImage;
    }

    im->setMetadata(imgsrc->getMetaData()->getRootExifData());

    im->saveTIFF(fname, 16, false, true);
    delete im;

    if (plistener) {
        plistener->setProgressState(false);
    }

    //im->saveJPEG (fname, 85);
}

void ImProcCoordinator::stopProcessing()
{

    updaterThreadStart.lock();

    if (updaterRunning && thread) {
        changeSinceLast = 0;
        thread->join();
    }

    updaterThreadStart.unlock();
}

void ImProcCoordinator::startProcessing()
{

#undef THREAD_PRIORITY_NORMAL

    if (!destroying) {
        if (!updaterRunning) {
            updaterThreadStart.lock();
            thread = nullptr;
            updaterRunning = true;
            updaterThreadStart.unlock();

            //batchThread->yield(); //the running batch should wait other threads to avoid conflict

            thread = Glib::Thread::create(sigc::mem_fun(*this, &ImProcCoordinator::process), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);

        }
    }
}

void ImProcCoordinator::startProcessing(int changeCode)
{
    paramsUpdateMutex.lock();
    changeSinceLast |= changeCode;
    paramsUpdateMutex.unlock();

    startProcessing();
}

void ImProcCoordinator::process()
{
    if (plistener) {
        plistener->setProgressState(true);
    }

    paramsUpdateMutex.lock();

    while (changeSinceLast) {
        const bool panningRelatedChange =
            params->toneCurve.isPanningRelatedChange(nextParams->toneCurve)
            || params->labCurve != nextParams->labCurve
            || params->locallab != nextParams->locallab
            || params->localContrast != nextParams->localContrast
            || params->rgbCurves != nextParams->rgbCurves
            || params->colorToning != nextParams->colorToning
            || params->vibrance != nextParams->vibrance
            || params->wb.isPanningRelatedChange(nextParams->wb)
            || params->colorappearance != nextParams->colorappearance
            || params->epd != nextParams->epd
            || params->fattal != nextParams->fattal
            || params->sh != nextParams->sh
            || params->crop != nextParams->crop
            || params->coarse != nextParams->coarse
            || params->commonTrans != nextParams->commonTrans
            || params->rotate != nextParams->rotate
            || params->distortion != nextParams->distortion
            || params->lensProf != nextParams->lensProf
            || params->perspective != nextParams->perspective
            || params->gradient != nextParams->gradient
            || params->pcvignette != nextParams->pcvignette
            || params->cacorrection != nextParams->cacorrection
            || params->vignetting != nextParams->vignetting
            || params->chmixer != nextParams->chmixer
            || params->blackwhite != nextParams->blackwhite
            || params->icm != nextParams->icm
            || params->hsvequalizer != nextParams->hsvequalizer
            || params->filmSimulation != nextParams->filmSimulation
            || params->softlight != nextParams->softlight
            || params->raw != nextParams->raw
            || params->retinex != nextParams->retinex
            || params->wavelet != nextParams->wavelet
            || params->dirpyrequalizer != nextParams->dirpyrequalizer
            || params->dehaze != nextParams->dehaze
            || params->pdsharpening != nextParams->pdsharpening
            || sharpMaskChanged;

        sharpMaskChanged = false;
        *params = *nextParams;
        int change = changeSinceLast;
        changeSinceLast = 0;
        paramsUpdateMutex.unlock();

        // M_VOID means no update, and is a bit higher that the rest
        if (change & (M_VOID - 1)) {
            updatePreviewImage(change, panningRelatedChange);
        }

        paramsUpdateMutex.lock();
    }

    paramsUpdateMutex.unlock();
    updaterRunning = false;

    if (plistener) {
        plistener->setProgressState(false);
    }
}

ProcParams* ImProcCoordinator::beginUpdateParams()
{
    paramsUpdateMutex.lock();

    return nextParams.get();
}

void ImProcCoordinator::endUpdateParams(ProcEvent change)
{
    int action = RefreshMapper::getInstance()->getAction(change);
    endUpdateParams(action);
}

void ImProcCoordinator::endUpdateParams(int changeFlags)
{
    changeSinceLast |= changeFlags;

    paramsUpdateMutex.unlock();
    startProcessing();
}

bool ImProcCoordinator::getHighQualComputed()
{
    // this function may only be called from detail windows
    if (!highQualityComputed) {
        if (options.prevdemo == PD_Sidecar) {
            // we already have high quality preview
            setHighQualComputed();
        } else {
            for (size_t i = 0; i < crops.size() - 1; ++i) { // -1, because last entry is the freshly created detail window
                if (crops[i]->get_skip() == 1) {   // there is at least one crop with skip == 1 => we already have high quality preview
                    setHighQualComputed();
                    break;
                }
            }
        }
    }

    return highQualityComputed;
}

void ImProcCoordinator::setHighQualComputed()
{
    highQualityComputed = true;
}

}
