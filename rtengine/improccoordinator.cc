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
#include "improccoordinator.h"
#include "curves.h"
#include "mytime.h"
#include "refreshmap.h"
#include "../rtgui/ppversion.h"
#include "colortemp.h"
#include "improcfun.h"
#include <iostream>
#include <fstream>
#include <string>
#include "../rtgui/md5helper.h"
#include "../rtgui/thresholdselector.h"
#include <unistd.h>

#include "iccstore.h"
#ifdef _OPENMP
#include <omp.h>
#endif
namespace rtengine
{

extern const Settings* settings;

ImProcCoordinator::ImProcCoordinator()
    : orig_prev(nullptr), oprevi(nullptr), oprevl(nullptr), nprevl(nullptr), reserv(nullptr), fattal_11_dcrop_cache(nullptr), previmg(nullptr), workimg(nullptr),
      ncie(nullptr), imgsrc(nullptr), lastAwbEqual(0.), lastAwbTempBias(0.0), ipf(&params, true), monitorIntent(RI_RELATIVE),
      softProof(false), gamutCheck(false), sharpMask(false), scale(10), highDetailPreprocessComputed(false), highDetailRawComputed(false),
      allocated(false), bwAutoR(-9000.f), bwAutoG(-9000.f), bwAutoB(-9000.f), CAMMean(NAN), coordX(0), coordY(0), localX(0), localY(0),
      ctColorCurve(),
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
      lllocalcurve(65536, 0),
      cclocalcurve(65536, 0),
      sklocalcurve(65536, 0),
      exlocalcurve(65536, 0),
      hltonecurveloc(65536, 0), //32768
      shtonecurveloc(65536, 0),
      tonecurveloc(65536, 0),
      lightCurveloc(32770, 0),
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
      rcurvehist(256), rcurvehistCropped(256), rbeforehist(256),
      gcurvehist(256), gcurvehistCropped(256), gbeforehist(256),
      bcurvehist(256), bcurvehistCropped(256), bbeforehist(256),
      fw(0), fh(0), tr(0),
      fullw(1), fullh(1),
      pW(-1), pH(-1),
      plistener(nullptr), awbListener(nullptr), imageListener(nullptr), aeListener(nullptr), acListener(nullptr), abwListener(nullptr), actListener(nullptr), adnListener(nullptr), awavListener(nullptr), dehaListener(nullptr), frameCountListener(nullptr), imageTypeListener(nullptr), hListener(nullptr),
      resultValid(false), lastOutputProfile("BADFOOD"), lastOutputIntent(RI__COUNT), lastOutputBPC(false), thread(nullptr), changeSinceLast(0), updaterRunning(false), destroying(false), utili(false), autili(false),
      butili(false), ccutili(false), cclutili(false), clcutili(false), opautili(false),  wavcontlutili(false),
      dataspot(nullptr), maxdata(0), retistr(nullptr), llstr(nullptr), lhstr(nullptr), ccstr(nullptr), hhstr(nullptr), skinstr(nullptr), pthstr(nullptr), exstr(nullptr),
      circrads(500, -10000),
      centerx(500, -10000),
      centery(500, -10000),
      centerxbufs(500, -10000),
      centerybufs(500, -10000),
      adjblurs(500, -10000),
      cutpasts(500, -10000),
      lastdusts(500, -10000),
      blurmets(500, -10000),
      dustmets(500, -10000),

      locx(500, -10000),
      locy(500, -10000),
      locxl(500, -10000),
      locyt(500, -10000),
      lights(500, -100000),
      contrs(500, -10000),
      chroms(500, -10000),
      sensis(500, -10000),
      expcomps(500, -10000),
      blacks(500, -10000),
      hlcomprs(500, -10000),
      hlcomprthreshs(500, -10000),
      shcomprs(500, -10000),
      sensiexs(500, -10000),

      transits(500, -10000),

      inverss(500, -10000),
      curvactivs(500, -10000),
      smeths(500, -10000),
      curens(500, -10000),
      radiuss(500, -10000),
      strengths(500, -10000),
      sensibns(500, -10000),
      inversrads(500, -10000),
      strs(500, 10000),
      chrrts(500, -10000),
      neighs(500, -10000),
      varts(500, -10000),
      sensihs(500, -10000),
      inversrets(500, -10000),
      retinexs(500, -10000),
      sps(500, -10000),
      sharradiuss(500, -10000),
      sharamounts(500, -10000),
      shardampings(500, -10000),
      inversshas(500, -10000),
      shariters(500, -10000),
      sensishas(500, -10000),
      qualitys(500, -10000),
      thress(500, -10000),
      proxis(500, -10000),
      noiselumfs(500, -10000),
      noiselumcs(500, -10000),
      noiselumdetails(500, -10000),
      noiselequals(500, -10000),
      noisechrodetails(500, -10000),
      bilaterals(500, -10000),
      sensidens(500, -10000),
      noisechrofs(500, -10000),
      noisechrocs(500, -10000),
      mult0s(500, -10000),
      mult1s(500, -10000),
      mult2s(500, -10000),
      mult3s(500, -10000),
      mult4s(500, -10000),
      chromacbdls(500, -10000),
      thresholds(500, -10000),
      sensicbs(500, -10000),
      activlums(500, -10000),
      versionmip(0),
      mipver(0),
      strens(500, -10000),
      gammas(500, -10000),
      estops(500, -10000),
      scaltms(500, -10000),
      reweis(500, -10000),
      sensitms(500, -10000),
      qualitycurves(500, -10000),
      sizeretics(500, -10000),
      reticurvs(25000, -10000),  //allow 500 values for each control point * 500
      retrabs(500, -10000),
      llcurvs(25000, -10000),  //allow 500 values for each control point * 500
      sizellcs(500, -10000),
      lhcurvs(25000, -10000),  //allow 500 values for each control point * 500
      hhcurvs(25000, -10000),  //allow 500 values for each control point * 500
      sizelhcs(500, -10000),
      sizehhcs(500, -10000),
      cccurvs(25000, -10000),  //allow 500 values for each control point * 500
      sizecccs(500, -10000),
      sensivs(500, -10000),
      saturateds(500, -10000),
      pastels(500, -10000),
      psthresholds(500, -10000),
      protectskinss(500, -10000),
      avoidcolorshifts(500, -10000),
      pastsattogs(500, -10000),
      skintonescurves(25000, -10000),
      sizeskintonecurves(500, -10000),
      excurves(25000, -10000),
      sizeexcurves(500, -10000),
      shapemets(500, -1000),
      exclumets(500, -1000),
      sensiexclus(500, -1000),
      strucs(500, -1000),
      warms(500, -1000),
      expdenois(500, -10000),
      expcolors(500, -10000),
      expvibrances(500, -10000),
      expblurs(500, -10000),
      exptonemaps(500, -10000),
      expretis(500, -10000),
      expsharps(500, -10000),
      expcbdls(500, -10000),
      expexposes(500, -10000),


      huerefs(500, -100000.f),
      huerefblurs(500, -100000.f),
      chromarefs(500, -100000.f),
      lumarefs(500, -100000.f),
      sobelrefs(500, -100000.f),
      huer(0),
      huerblu(0),
      chromar(0),
      lumar(0),
      sobeler(0),
      colourToningSatLimit(0.f), colourToningSatLimitOpacity(0.f), lastspotdup(false), highQualityComputed(false),

      retistrsav(nullptr)
{}

void ImProcCoordinator::assign(ImageSource* imgsrc)
{
    this->imgsrc = imgsrc;
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
    updaterThreadStart.unlock();
}

DetailedCrop* ImProcCoordinator::createCrop(::EditDataProvider *editDataProvider, bool isDetailWindow)
{

    return new Crop(this, editDataProvider, isDetailWindow);
}


// todo: bitmask containing desired actions, taken from changesSinceLast
// cropCall: calling crop, used to prevent self-updates  ...doesn't seem to be used
void ImProcCoordinator::updatePreviewImage(int todo, Crop* cropCall)
{
    // TODO Locallab printf
    printf("updatePreviewImage\n");

    MyMutex::MyLock processingLock(mProcessing);
    int numofphases = 14;
    int readyphase = 0;

    bwAutoR = bwAutoG = bwAutoB = -9000.f;

    if (todo == CROP && ipf.needsPCVignetting()) {
        todo |= TRANSFORM;    // Change about Crop does affect TRANSFORM
    }

    bool highDetailNeeded = false;

    if (options.prevdemo == PD_Sidecar) {
        highDetailNeeded = true;    //i#2664
    } else {
        highDetailNeeded = (todo & M_HIGHQUAL);
    }

    // Check if any detail crops need high detail. If not, take a fast path short cut
    if (!highDetailNeeded) {
        for (size_t i = 0; i < crops.size(); i++)
            if (crops[i]->get_skip() == 1) {   // skip=1 -> full resolution
                highDetailNeeded = true;
                break;
            }
    }

    RAWParams rp = params.raw;
    ColorManagementParams cmp = params.icm;
    LCurveParams  lcur = params.labCurve;

    if (!highDetailNeeded) {
        // if below 100% magnification, take a fast path
        if (rp.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::NONE) && rp.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::NONE)) {
            rp.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
        }

        //bayerrp.all_enhance = false;

        if (rp.xtranssensor.method != RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::NONE) && rp.xtranssensor.method != RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::NONE)) {
            rp.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);
        }

        rp.bayersensor.ccSteps = 0;
        rp.xtranssensor.ccSteps = 0;
        //rp.deadPixelFilter = rp.hotPixelFilter = false;
    }

    progress("Applying white balance, color correction & sRGB conversion...", 100 * readyphase / numofphases);

    if (frameCountListener) {
        frameCountListener->FrameCountChanged(imgsrc->getFrameCount(), params.raw.bayersensor.imageNum);
    }

    // raw auto CA is bypassed if no high detail is needed, so we have to compute it when high detail is needed
    if ((todo & M_PREPROC) || (!highDetailPreprocessComputed && highDetailNeeded)) {
        imgsrc->setCurrentFrame(params.raw.bayersensor.imageNum);

        imgsrc->preprocess(rp, params.lensProf, params.coarse);
        imgsrc->getRAWHistogram(histRedRaw, histGreenRaw, histBlueRaw);

        highDetailPreprocessComputed = highDetailNeeded;
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
            || (params.toneCurve.hrenabled && params.toneCurve.method != "Color" && imgsrc->isRGBSourceModified())
            || (!params.toneCurve.hrenabled && params.toneCurve.method == "Color" && imgsrc->isRGBSourceModified())) {

        if (settings->verbose) {
            if (imgsrc->getSensorType() == ST_BAYER) {
                printf("Demosaic Bayer image n.%d using method: %s\n", rp.bayersensor.imageNum + 1, rp.bayersensor.method.c_str());
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                printf("Demosaic X-Trans image with using method: %s\n", rp.xtranssensor.method.c_str());
            }
        }

        imgsrc->demosaic(rp);   //enabled demosaic
        // if a demosaic happened we should also call getimage later, so we need to set the M_INIT flag
        todo |= M_INIT;

        if (highDetailNeeded) {
            highDetailRawComputed = true;
        } else {
            highDetailRawComputed = false;
        }

        if (params.retinex.enabled) {
            lhist16RETI(32768);
            lhist16RETI.clear();

            imgsrc->retinexPrepareBuffers(params.icm, params.retinex, conversionBuffer, lhist16RETI);
        }
    }

    if ((todo & (M_RETINEX | M_INIT)) && params.retinex.enabled) {
        bool dehacontlutili = false;
        bool mapcontlutili = false;
        bool useHsl = false;
        LUTf cdcurve(65536, 0);
        LUTf mapcurve(65536, 0);

        imgsrc->retinexPrepareCurves(params.retinex, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, dehacontlutili, mapcontlutili, useHsl, lhist16RETI, histLRETI);
        float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
        imgsrc->retinex(params.icm, params.retinex,  params.toneCurve, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, conversionBuffer, dehacontlutili, mapcontlutili, useHsl, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax, histLRETI);   //enabled Retinex

        if (dehaListener) {
            dehaListener->minmaxChanged(maxCD, minCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
        }
    }

    if (todo & (M_INIT | M_LINDENOISE | M_HDR)) {
        MyMutex::MyLock initLock(minit);  // Also used in crop window

        imgsrc->HLRecovery_Global(params.toneCurve);   // this handles Color HLRecovery


        if (settings->verbose) {
            printf("Applying white balance, color correction & sRBG conversion...\n");
        }

        currWB = ColorTemp(params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

        if (!params.wb.enabled) {
            currWB = ColorTemp();
        } else if (params.wb.method == "Camera") {
            currWB = imgsrc->getWB();
        } else if (params.wb.method == "Auto") {
            if (lastAwbEqual != params.wb.equal || lastAwbTempBias != params.wb.tempBias) {
                double rm, gm, bm;
                imgsrc->getAutoWBMultipliers(rm, gm, bm);

                if (rm != -1.) {
                    autoWB.update(rm, gm, bm, params.wb.equal, params.wb.tempBias);
                    lastAwbEqual = params.wb.equal;
                    lastAwbTempBias = params.wb.tempBias;
                } else {
                    lastAwbEqual = -1.;
                    lastAwbTempBias = 0.0;
                    autoWB.useDefaults(params.wb.equal);
                }

                //double rr,gg,bb;
                //autoWB.getMultipliers(rr,gg,bb);
            }

            currWB = autoWB;
        }

        if (params.wb.enabled) {
            params.wb.temperature = currWB.getTemp();
            params.wb.green = currWB.getGreen();
        }

        if (params.wb.method == "Auto" && awbListener && params.wb.enabled) {
            awbListener->WBChanged(params.wb.temperature, params.wb.green);
        }

        int tr = getCoarseBitMask(params.coarse);

        imgsrc->getFullSize(fw, fh, tr);

        // Will (re)allocate the preview's buffers
        setScale(scale);
        PreviewProps pp(0, 0, fw, fh, scale);
        // Tells to the ImProcFunctions' tools what is the preview scale, which may lead to some simplifications
        ipf.setScale(scale);

        imgsrc->getImage(currWB, tr, orig_prev, pp, params.toneCurve, params.raw);
        denoiseInfoStore.valid = false;
        //ColorTemp::CAT02 (orig_prev, &params) ;
        //   printf("orig_prevW=%d\n  scale=%d",orig_prev->width, scale);
        /* Issue 2785, disabled some 1:1 tools
                if (todo & M_LINDENOISE) {
                    DirPyrDenoiseParams denoiseParams = params.dirpyrDenoise;
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
                            imgsrc->convertColorSpace(calclum, params.icm, currWB);//calculate values after colorspace conversion
                        }

                        int kall=1;
                        ipf.RGB_denoise(kall, orig_prev, orig_prev, calclum, ch_M, max_r, max_b, imgsrc->isRAW(), denoiseParams, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, nresi, highresi);
                    }
                }
        */
        imgsrc->convertColorSpace(orig_prev, params.icm, currWB);

        ipf.firstAnalysis(orig_prev, params, vhist16);
    }

    readyphase++;

    if ((todo & M_HDR) && params.fattal.enabled) {
        if (fattal_11_dcrop_cache) {
            delete fattal_11_dcrop_cache;
            fattal_11_dcrop_cache = nullptr;
        }

        ipf.ToneMapFattal02(orig_prev);

        if (oprevi != orig_prev) {
            delete oprevi;
        }
    }

    oprevi = orig_prev;

    progress("Rotate / Distortion...", 100 * readyphase / numofphases);
    // Remove transformation if unneeded
    bool needstransform = ipf.needsTransform();

    if ((needstransform || ((todo & (M_TRANSFORM | M_RGBCURVE))  && params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled))) {
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

    if ((todo & (M_TRANSFORM | M_RGBCURVE))  && params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled) {
        const int W = oprevi->getWidth();
        const int H = oprevi->getHeight();
        LabImage labcbdl(W, H);
        ipf.rgb2lab(*oprevi, labcbdl, params.icm.working);
        ipf.dirpyrequalizer(&labcbdl, scale);
        ipf.lab2rgb(labcbdl, *oprevi, params.icm.working);
    }

    readyphase++;
    progress("Preparing shadow/highlight map...", 100 * readyphase / numofphases);

    readyphase++;

    if (todo & M_AUTOEXP) {
        if (params.toneCurve.autoexp) {
            LUTu aehist;
            int aehistcompr;
            imgsrc->getAutoExpHistogram(aehist, aehistcompr);
            ipf.getAutoExp(aehist, aehistcompr, params.toneCurve.clip, params.toneCurve.expcomp,
                           params.toneCurve.brightness, params.toneCurve.contrast, params.toneCurve.black, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh);

            if (aeListener)
                aeListener->autoExpChanged(params.toneCurve.expcomp, params.toneCurve.brightness, params.toneCurve.contrast,
                                           params.toneCurve.black, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, params.toneCurve.hrenabled);
        }

        if (params.toneCurve.histmatching) {
            imgsrc->getAutoMatchedToneCurve(params.icm, params.toneCurve.curve);

            if (params.toneCurve.autoexp) {
                params.toneCurve.expcomp = 0.0;
            }

            params.toneCurve.autoexp = false;
            params.toneCurve.curveMode = ToneCurveParams::TcMode::FILMLIKE;
            params.toneCurve.curve2 = { 0 };
            params.toneCurve.brightness = 0;
            params.toneCurve.contrast = 0;
            params.toneCurve.black = 0;

            if (aeListener) {
                aeListener->autoMatchedToneCurveChanged(params.toneCurve.curveMode, params.toneCurve.curve);
            }
        }
    }

    progress("Exposure curve & CIELAB conversion...", 100 * readyphase / numofphases);

    if ((todo & M_RGBCURVE) || (todo & M_CROP)) {
//        if (hListener) oprevi->calcCroppedHistogram(params, scale, histCropped);

        //complexCurve also calculated pre-curves histogram depending on crop
        CurveFactory::complexCurve(params.toneCurve.expcomp, params.toneCurve.black / 65535.0,
                                   params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh,
                                   params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast,
                                   params.toneCurve.curve, params.toneCurve.curve2,
                                   vhist16, hltonecurve, shtonecurve, tonecurve, histToneCurve, customToneCurve1, customToneCurve2, 1);

        CurveFactory::RGBCurve(params.rgbCurves.rcurve, rCurve, 1);
        CurveFactory::RGBCurve(params.rgbCurves.gcurve, gCurve, 1);
        CurveFactory::RGBCurve(params.rgbCurves.bcurve, bCurve, 1);


        opautili = false;

        if (params.colorToning.enabled) {
            TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params.icm.working);
            double wp[3][3] = {
                {wprof[0][0], wprof[0][1], wprof[0][2]},
                {wprof[1][0], wprof[1][1], wprof[1][2]},
                {wprof[2][0], wprof[2][1], wprof[2][2]}
            };
            params.colorToning.getCurves(ctColorCurve, ctOpacityCurve, wp, opautili);
            CurveFactory::curveToning(params.colorToning.clcurve, clToningcurve, scale == 1 ? 1 : 16);
            CurveFactory::curveToning(params.colorToning.cl2curve, cl2Toningcurve, scale == 1 ? 1 : 16);
        }

        if (params.blackwhite.enabled) {
            CurveFactory::curveBW(params.blackwhite.beforeCurve, params.blackwhite.afterCurve, vhist16bw, histToneCurveBW, beforeToneCurveBW, afterToneCurveBW, 1);
        }

        colourToningSatLimit = float (params.colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
        colourToningSatLimitOpacity = 1.f - (float (params.colorToning.saturatedOpacity) / 100.f);

        int satTH = 80;
        int satPR = 30;
        int indi = 0;

        if (params.colorToning.enabled  && params.colorToning.autosat && params.colorToning.method != "LabGrid") { //for colortoning evaluation of saturation settings
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

        if (actListener) {
            //if(params.blackwhite.enabled) {actListener->autoColorTonChanged(0, satTH, satPR);}
            if (params.blackwhite.enabled && params.colorToning.autosat) {
                actListener->autoColorTonChanged(0, satTH, satPR);    //hide sliders only if autosat
                indi = 0;
            } else {
                if (params.colorToning.autosat) {
                    if (params.colorToning.method == "Lab") {
                        indi = 1;
                    } else if (params.colorToning.method == "RGBCurves") {
                        indi = 1;
                    } else if (params.colorToning.method == "RGBSliders") {
                        indi = 1;
                    } else if (params.colorToning.method == "Splico") {
                        indi = 2;
                    } else if (params.colorToning.method == "Splitlr") {
                        indi = 2;
                    }

                    //actListener->autoColorTonChanged(indi, satTH, satPR);
                }
            }
        }

        // if it's just crop we just need the histogram, no image updates
        if (todo & M_RGBCURVE) {
            //initialize rrm bbm ggm different from zero to avoid black screen in some cases
            double rrm = 33.;
            double ggm = 33.;
            double bbm = 33.;

            DCPProfile::ApplyState as;
            DCPProfile *dcpProf = imgsrc->getDCP(params.icm, as);

            ipf.rgbProc(oprevi, oprevl, nullptr, hltonecurve, shtonecurve, tonecurve, params.toneCurve.saturation,
                        rCurve, gCurve, bCurve, colourToningSatLimit, colourToningSatLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, beforeToneCurveBW, afterToneCurveBW, rrm, ggm, bbm, bwAutoR, bwAutoG, bwAutoB, params.toneCurve.expcomp, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, dcpProf, as, histToneCurve);

            if (params.blackwhite.enabled && params.blackwhite.autoc && abwListener) {
                if (settings->verbose) {
                    printf("ImProcCoordinator / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", bwAutoR, bwAutoG, bwAutoB);
                }

                abwListener->BWChanged((float) rrm, (float) ggm, (float) bbm);
            }

            if (params.colorToning.autosat && actListener) {
                if (settings->verbose) {
                    printf("ImProcCoordinator / Auto CT:  indi=%d   satH=%d  satPR=%d\n", indi, (int)colourToningSatLimit, (int) colourToningSatLimitOpacity);
                }

                actListener->autoColorTonChanged(indi, (int) colourToningSatLimit, (int)colourToningSatLimitOpacity);  //change sliders autosat
            }

            // correct GUI black and white with value
        }

        //  ipf.Lab_Tile(oprevl, oprevl, scale);

        // compute L channel histogram
        int x1, y1, x2, y2;
        params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
    }

    readyphase++;
    lhist16(32768);

    if (todo & (M_LUMACURVE | M_CROP)) {
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
        CurveFactory::complexLCurve(params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, lhist16, lumacurve, histLCurve, scale == 1 ? 1 : 16, utili);
    }

    if (todo & M_LUMACURVE) {

        CurveFactory::curveCL(clcutili, params.labCurve.clcurve, clcurve, scale == 1 ? 1 : 16);

        CurveFactory::complexsgnCurve(autili, butili, ccutili, cclutili, params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,
                                      params.labCurve.lccurve, chroma_acurve, chroma_bcurve, satcurve, lhskcurve, scale == 1 ? 1 : 16);
    }

    //scale = 1;
    if (todo & (M_LUMINANCE + M_COLOR)) {
        nprevl->CopyFrom(oprevl);
        reserv->CopyFrom(oprevl);

        int maxspot = 1;
        progress("Applying Color Boost...", 100 * readyphase / numofphases);

        //*************************************************************
        // locallab
        //*************************************************************

        if (params.locallab.enabled) {
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
             *  2018 Pierre Cabrera <pierre.cab@gmail.com>
             */

            float **shbuffer = nullptr;
            int sca = 1;

            for (int sp = 0; sp < params.locallab.nbspot; sp++) {
                // Set local curves of current spot to LUT
                LHutili = false;
                HHutili = false;
                locallutili = false;
                localexutili = false;
                localcutili = false;
                localskutili = false;
                locRETgainCurve.Set(params.locallab.localTgaincurve.at(sp));
                loclhCurve.Set(params.locallab.LHcurve.at(sp), LHutili);
                lochhCurve.Set(params.locallab.HHcurve.at(sp), HHutili);
                CurveFactory::curveLocal(locallutili, params.locallab.llcurve.at(sp), lllocalcurve, sca);
                CurveFactory::curveCCLocal(localcutili, params.locallab.cccurve.at(sp), cclocalcurve, sca);
                CurveFactory::curveskLocal(localskutili, params.locallab.skintonescurve.at(sp), sklocalcurve, sca);
                CurveFactory::curveexLocal(localexutili, params.locallab.excurve.at(sp), exlocalcurve, sca);
                double ecomp = params.locallab.expcomp.at(sp);
                double black = params.locallab.black.at(sp);
                double hlcompr = params.locallab.hlcompr.at(sp);
                double hlcomprthresh = params.locallab.hlcomprthresh.at(sp);
                double shcompr = params.locallab.shcompr.at(sp);
                double br = params.locallab.lightness.at(sp);
                CurveFactory::complexCurvelocal(ecomp, black / 65535., hlcompr, hlcomprthresh, shcompr, br,
                                                hltonecurveloc, shtonecurveloc, tonecurveloc, lightCurveloc,
                                                sca);

                // Reference parameters computation
                double huere, chromare, lumare, huerefblu, sobelre;

                if (params.locallab.spotMethod.at(sp) == "exc") {
                    ipf.calc_ref(sp, reserv, reserv, 0, 0, pW, pH, scale, huerefblu, huere, chromare, lumare, sobelre);
                } else {
                    ipf.calc_ref(sp, nprevl, nprevl, 0, 0, pW, pH, scale, huerefblu, huere, chromare, lumare, sobelre);
                }

                huerblu = huerefblurs[sp] = huerefblu;
                huer = huerefs[sp] = huere;
                chromar = chromarefs[sp] = chromare;
                lumar = lumarefs[sp] = lumare ;
                sobeler = sobelrefs[sp] = sobelre;

                printf("huerblu=%f, huer=%f, chromar=%f, lumar=%f, sobeler=%f\n", huerblu, huer, chromar, lumar, sobeler);

                // Locallab tools computation
                /* Notes:
                 * - maxspot, huerefs, centerx and centery aren't used in Lab_Local (only for printf) so values aren't important
                 * - shbuffer is used as nullptr
                 */
                ipf.Lab_Local(3, maxspot, sp, huerefs, sobelrefs, centerx, centery, (float**)shbuffer, nprevl, nprevl, reserv, 0, 0, pW, pH, scale, locRETgainCurve, lllocalcurve, loclhCurve,  lochhCurve,
                              LHutili, HHutili, cclocalcurve, localskutili, sklocalcurve, localexutili, exlocalcurve, hltonecurveloc, shtonecurveloc, tonecurveloc, lightCurveloc, huerblu, huer, chromar, lumar, sobeler);

                // Clear local curves
                lllocalcurve.clear();
                cclocalcurve.clear();
                sklocalcurve.clear();
                exlocalcurve.clear();
            }
        }

        //*************************************************************
        // end locallab
        //*************************************************************

        histCCurve.clear();
        histLCurve.clear();
        ipf.chromiLuminanceCurve(nullptr, pW, nprevl, nprevl, chroma_acurve, chroma_bcurve, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, histCCurve, histLCurve);
        ipf.vibrance(nprevl);

        if ((params.colorappearance.enabled && !params.colorappearance.tonecie) || (!params.colorappearance.enabled)) {
            ipf.EPDToneMap(nprevl, 5, scale);
        }

        // for all treatments Defringe, Sharpening, Contrast detail , Microcontrast they are activated if "CIECAM" function are disabled
        readyphase++;

        /* Issue 2785, disabled some 1:1 tools
                if (scale==1) {
                    if((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)){
                        progress ("Denoising luminance impulse...",100*readyphase/numofphases);
                        ipf.impulsedenoise (nprevl);
                        readyphase++;
                    }
                    if((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)){
                        progress ("Defringing...",100*readyphase/numofphases);
                        ipf.defringe (nprevl);
                        readyphase++;
                    }
                    if (params.sharpenEdge.enabled) {
                        progress ("Edge sharpening...",100*readyphase/numofphases);
                        ipf.MLsharpen (nprevl);
                        readyphase++;
                    }
                    if (params.sharpenMicro.enabled) {
                        if(( params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)){
                            progress ("Microcontrast...",100*readyphase/numofphases);
                            ipf.MLmicrocontrast (nprevl);
                            readyphase++;
                        }
                    }
                    if(((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) && params.sharpening.enabled) {
                        progress ("Sharpening...",100*readyphase/numofphases);

                        float **buffer = new float*[pH];
                        for (int i=0; i<pH; i++)
                            buffer[i] = new float[pW];

                        ipf.sharpening (nprevl, (float**)buffer);

                        for (int i=0; i<pH; i++)
                            delete [] buffer[i];
                        delete [] buffer;
                        readyphase++;
                    }
                }
        */
        if (params.dirpyrequalizer.cbdlMethod == "aft") {
            if (((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled))) {
                progress("Pyramid wavelet...", 100 * readyphase / numofphases);
                ipf.dirpyrequalizer(nprevl, scale);
                //ipf.Lanczoslab (ip_wavelet(LabImage * lab, LabImage * dst, const procparams::EqualizerParams & eqparams), nprevl, 1.f/scale);
                readyphase++;
            }
        }


        wavcontlutili = false;
        //CurveFactory::curveWavContL ( wavcontlutili,params.wavelet.lcurve, wavclCurve, LUTu & histogramwavcl, LUTu & outBeforeWavCLurveHistogram,int skip);
        CurveFactory::curveWavContL(wavcontlutili, params.wavelet.wavclCurve, wavclCurve, scale == 1 ? 1 : 16);


        if ((params.wavelet.enabled)) {
            WaveletParams WaveParams = params.wavelet;
            //      WaveParams.getCurves(wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY);
            WaveParams.getCurves(wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL);

            int kall = 0;
            progress("Wavelet...", 100 * readyphase / numofphases);
            //  ipf.ip_wavelet(nprevl, nprevl, kall, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, scale);
            ipf.ip_wavelet(nprevl, nprevl, kall, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL, wavclCurve, scale);

        }


        if (params.colorappearance.enabled) {
            //L histo  and Chroma histo for ciecam
            // histogram well be for Lab (Lch) values, because very difficult to do with J,Q, M, s, C
            int x1, y1, x2, y2;
            params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
            lhist16CAM.clear();
            lhist16CCAM.clear();

            if (!params.colorappearance.datacie) {
                for (int x = 0; x < pH; x++)
                    for (int y = 0; y < pW; y++) {
                        int pos = CLIP((int)(nprevl->L[x][y]));
                        int posc = CLIP((int)sqrt(nprevl->a[x][y] * nprevl->a[x][y] + nprevl->b[x][y] * nprevl->b[x][y]));
                        lhist16CAM[pos]++;
                        lhist16CCAM[posc]++;
                    }
            }

            CurveFactory::curveLightBrightColor(params.colorappearance.curve, params.colorappearance.curve2, params.colorappearance.curve3,
                                                lhist16CAM, histLCAM, lhist16CCAM, histCCAM,
                                                customColCurve1, customColCurve2, customColCurve3, 1);

            const FramesMetaData* metaData = imgsrc->getMetaData();
            int imgNum = 0;

            if (imgsrc->isRAW()) {
                if (imgsrc->getSensorType() == ST_BAYER) {
                    imgNum = rtengine::LIM<unsigned int> (params.raw.bayersensor.imageNum, 0, metaData->getFrameCount() - 1);
                } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                    //imgNum = rtengine::LIM<unsigned int>(params.raw.xtranssensor.imageNum, 0, metaData->getFrameCount() - 1);
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
                E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                E_V += log2(params.raw.expos);  // exposure raw white point ; log2 ==> linear to EV
                adap = powf(2.f, E_V - 3.f);  // cd / m2
                // end calculation adaptation scene luminosity
            }

            float d, dj, yb;
            bool execsharp = false;

            if (!ncie) {
                ncie = new CieImage(pW, pH);
            }

            if (!CAMBrightCurveJ && (params.colorappearance.algo == "JC" || params.colorappearance.algo == "JS" || params.colorappearance.algo == "ALL")) {
                CAMBrightCurveJ(32768, 0);
            }

            if (!CAMBrightCurveQ && (params.colorappearance.algo == "QM" || params.colorappearance.algo == "ALL")) {
                CAMBrightCurveQ(32768, 0);
            }

            // Issue 2785, only float version of ciecam02 for navigator and pan background
            CAMMean = NAN;
            CAMBrightCurveJ.dirty = true;
            CAMBrightCurveQ.dirty = true;

            ipf.ciecam_02float(ncie, float (adap), pW, 2, nprevl, &params, customColCurve1, customColCurve2, customColCurve3, histLCAM, histCCAM, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, scale, execsharp, d, dj, yb, 1);

            if ((params.colorappearance.autodegree || params.colorappearance.autodegreeout) && acListener && params.colorappearance.enabled) {
                acListener->autoCamChanged(100.* (double)d, 100.* (double)dj);
            }

            if (params.colorappearance.autoadapscen && acListener && params.colorappearance.enabled) {
                acListener->adapCamChanged(adap);    //real value of adapt scene
            }

            if (params.colorappearance.autoybscen && acListener && params.colorappearance.enabled) {
                acListener->ybCamChanged((int) yb);    //real value Yb scene
            }

            readyphase++;
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
    if ((todo & M_MONITOR) || (lastOutputProfile != params.icm.output) || lastOutputIntent != params.icm.outputIntent || lastOutputBPC != params.icm.outputBPC) {
        lastOutputProfile = params.icm.output;
        lastOutputIntent = params.icm.outputIntent;
        lastOutputBPC = params.icm.outputBPC;
        ipf.updateColorProfiles(monitorProfile, monitorIntent, softProof, gamutCheck);
    }

// process crop, if needed
    for (size_t i = 0; i < crops.size(); i++)
        if (crops[i]->hasListener() && cropCall != crops[i]) {
            crops[i]->update(todo);     // may call ourselves
        }

    progress("Conversion to RGB...", 100 * readyphase / numofphases);

    if ((todo != CROP && todo != MINUPDATE) || (todo & M_MONITOR)) {
        MyMutex::MyLock prevImgLock(previmg->getMutex());

        try {
            // Computing the preview image, i.e. converting from WCS->Monitor color space (soft-proofing disabled) or WCS->Printer profile->Monitor color space (soft-proofing enabled)
            ipf.lab2monitorRgb(nprevl, previmg);

            // Computing the internal image for analysis, i.e. conversion from WCS->Output profile
            delete workimg;
            workimg = ipf.lab2rgb(nprevl, 0, 0, pW, pH, params.icm);
        } catch (char * str) {
            progress("Error converting file...", 0);
            return;
        }
    }

    if (!resultValid) {
        resultValid = true;

        if (imageListener) {
            imageListener->setImage(previmg, scale, params.crop);
        }
    }

    if (imageListener)
        // TODO: The WB tool should be advertised too in order to get the AutoWB's temp and green values
    {
        imageListener->imageReady(params.crop);
    }

    readyphase++;

    if (hListener) {
        updateLRGBHistograms();
        hListener->histogramChanged(histRed, histGreen, histBlue, histLuma, histToneCurve, histLCurve, histCCurve, /*histCLurve, histLLCurve,*/ histLCAM, histCCAM, histRedRaw, histGreenRaw, histBlueRaw, histChroma, histLRETI);
    }

}


void ImProcCoordinator::freeAll()
{

    if (settings->verbose) {
        printf("freeall starts %d\n", (int)allocated);
    }

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
        delete reserv;
        reserv    = nullptr;

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

    if (settings->verbose) {
        printf("setscale before lock\n");
    }

    tr = getCoarseBitMask(params.coarse);

    int nW, nH;
    imgsrc->getFullSize(fw, fh, tr);

    prevscale++;

    do {
        prevscale--;
        PreviewProps pp(0, 0, fw, fh, prevscale);
        imgsrc->getSize(pp, nW, nH);
    } while (nH < 400 && prevscale > 1 && (nW * nH < 1000000));  // sctually hardcoded values, perhaps a better choice is possible

    if (settings->verbose) {
        printf("setscale starts (%d, %d)\n", nW, nH);
    }

    if (nW != pW || nH != pH) {

        freeAll();

        pW = nW;
        pH = nH;

        orig_prev = new Imagefloat(pW, pH);
        oprevi = orig_prev;
        oprevl = new LabImage(pW, pH);
        nprevl = new LabImage(pW, pH);
        reserv = new LabImage(pW, pH);

        //  nprevloc = new LabImage (pW, pH);
        //ncie is only used in ImProcCoordinator::updatePreviewImage, it will be allocated on first use and deleted if not used anymore
        previmg = new Image8(pW, pH);
        workimg = new Image8(pW, pH);

        allocated = true;
    }

    scale = prevscale;
    resultValid = false;
    fullw = fw;
    fullh = fh;

    if (settings->verbose) {
        printf("setscale ends\n");
    }

    if (!sizeListeners.empty())
        for (size_t i = 0; i < sizeListeners.size(); i++) {
            sizeListeners[i]->sizeChanged(fullw, fullh, fw, fh);
        }

    if (settings->verbose) {
        printf("setscale ends2\n");
    }

}


void ImProcCoordinator::updateLRGBHistograms()
{

    int x1, y1, x2, y2;
    params.crop.mapToResized(pW, pH, scale, x1, x2, y1, y2);

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            histChroma.clear();

            for (int i = y1; i < y2; i++)
                for (int j = x1; j < x2; j++)
                {
                    histChroma[(int)(sqrtf(SQR(nprevl->a[i][j]) + SQR(nprevl->b[i][j])) / 188.f)]++;      //188 = 48000/256
                }
        }
        #pragma omp section
        {
            histLuma.clear();

            for (int i = y1; i < y2; i++)
                for (int j = x1; j < x2; j++)
                {
                    histLuma[(int)(nprevl->L[i][j] / 128.f)]++;
                }
        }
        #pragma omp section
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

void ImProcCoordinator::progress(Glib::ustring str, int pr)
{

    /*  if (plistener) {
        plistener->setProgressStr (str);
        plistener->setProgress ((double)pr / 100.0);
      }*/
}

bool ImProcCoordinator::getAutoWB(double& temp, double& green, double equal, double tempBias)
{

    if (imgsrc) {
        if (lastAwbEqual != equal || lastAwbTempBias != tempBias) {
// Issue 2500            MyMutex::MyLock lock(minit);  // Also used in crop window
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers(rm, gm, bm);

            if (rm != -1) {
                autoWB.update(rm, gm, bm, equal, tempBias);
                lastAwbEqual = equal;
                lastAwbTempBias = tempBias;
            } else {
                lastAwbEqual = -1.;
                autoWB.useDefaults(equal);
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

void ImProcCoordinator::getCamWB(double & temp, double & green)
{

    if (imgsrc) {
        temp = imgsrc->getWB().getTemp();
        green = imgsrc->getWB().getGreen();
    }
}

void ImProcCoordinator::getSpotWB(int x, int y, int rect, double & temp, double & tgreen)
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

        int tr = getCoarseBitMask(params.coarse);

        ret = imgsrc->getSpotWB(red, green, blue, tr, params.wb.equal);
        currWB = ColorTemp(params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);
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

void ImProcCoordinator::getAutoCrop(double ratio, int &x, int &y, int &w, int &h)
{

    MyMutex::MyLock lock(mProcessing);

    LensCorrection *pLCPMap = nullptr;

    if (params.lensProf.useLcp() && imgsrc->getMetaData()->getFocalLen() > 0) {
        const std::shared_ptr<LCPProfile> pLCPProf = LCPStore::getInstance()->getProfile(params.lensProf.lcpFile);

        if (pLCPProf) pLCPMap = new LCPMapper(pLCPProf, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(), imgsrc->getMetaData()->getFocusDist(),
                                                  0, false, params.lensProf.useDist, fullw, fullh, params.coarse, imgsrc->getRotateDegree());
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

void ImProcCoordinator::setMonitorProfile(const Glib::ustring & profile, RenderingIntent intent)
{
    monitorProfile = profile;
    monitorIntent = intent;
}

void ImProcCoordinator::getMonitorProfile(Glib::ustring & profile, RenderingIntent & intent) const
{
    profile = monitorProfile;
    intent = monitorIntent;
}

void ImProcCoordinator::setSoftProofing(bool softProof, bool gamutCheck)
{
    this->softProof = softProof;
    this->gamutCheck = gamutCheck;
}

void ImProcCoordinator::getSoftProofing(bool & softProof, bool & gamutCheck)
{
    softProof = this->softProof;
    gamutCheck = this->gamutCheck;
}

void ImProcCoordinator::setSharpMask(bool sharpMask)
{
    this->sharpMask = sharpMask;
}

void ImProcCoordinator::saveInputICCReference(const Glib::ustring & fname, bool apply_wb)
{

    MyMutex::MyLock lock(mProcessing);

    int fW, fH;

    int tr = getCoarseBitMask(params.coarse);

    imgsrc->getFullSize(fW, fH, tr);
    PreviewProps pp(0, 0, fW, fH, 1);
    ProcParams ppar = params;
    ppar.toneCurve.hrenabled = false;
    ppar.icm.input = "(none)";
    Imagefloat* im = new Imagefloat(fW, fH);
    imgsrc->preprocess(ppar.raw, ppar.lensProf, ppar.coarse);
    imgsrc->demosaic(ppar.raw);
    ColorTemp currWB = ColorTemp(params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

    if (params.wb.method == "Camera") {
        currWB = imgsrc->getWB();
    } else if (params.wb.method == "Auto") {
        if (lastAwbEqual != params.wb.equal || lastAwbTempBias != params.wb.tempBias) {
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers(rm, gm, bm);

            if (rm != -1.) {
                autoWB.update(rm, gm, bm, params.wb.equal, params.wb.tempBias);
                lastAwbEqual = params.wb.equal;
                lastAwbTempBias = params.wb.tempBias;
            } else {
                lastAwbEqual = -1.;
                lastAwbTempBias = 0.0;
                autoWB.useDefaults(params.wb.equal);
            }
        }

        currWB = autoWB;
    }

    if (!apply_wb) {
        currWB = ColorTemp(); // = no white balance
    }

    imgsrc->getImage(currWB, tr, im, pp, ppar.toneCurve, ppar.raw);
    ImProcFunctions ipf(&ppar, true);

    if (ipf.needsTransform()) {
        Imagefloat* trImg = new Imagefloat(fW, fH);
        ipf.transform(im, trImg, 0, 0, 0, 0, fW, fH, fW, fH,
                      imgsrc->getMetaData(), imgsrc->getRotateDegree(), true);
        delete im;
        im = trImg;
    }

    if (params.crop.enabled) {
        Imagefloat *tmpim = new Imagefloat(params.crop.w, params.crop.h);
        int cx = params.crop.x;
        int cy = params.crop.y;
        int cw = params.crop.w;
        int ch = params.crop.h;
        #pragma omp parallel for

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
    #pragma omp parallel for

    for (int i = 0; i < im->getHeight(); i++) {
        for (int j = 0; j < im->getWidth(); j++) {
            im->r(i, j) = CLIP(im->r(i, j));
            im->g(i, j) = CLIP(im->g(i, j));
            im->b(i, j) = CLIP(im->b(i, j));
        }
    }

    int imw, imh;
    double tmpScale = ipf.resizeScale(&params, fW, fH, imw, imh);

    if (tmpScale != 1.0) {
        Imagefloat* tempImage = new Imagefloat(imw, imh);
        ipf.resize(im, tempImage, tmpScale);
        delete im;
        im = tempImage;
    }

    im->setMetadata(imgsrc->getMetaData()->getRootExifData());

    im->saveTIFF(fname, 16, true);
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
        params = nextParams;
        int change = changeSinceLast;
        changeSinceLast = 0;
        paramsUpdateMutex.unlock();

        // M_VOID means no update, and is a bit higher that the rest
        if (change & (M_VOID - 1)) {
            updatePreviewImage(change);
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

    return &nextParams;
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

void ImProcCoordinator::spotduplic(int **dataspot, int spottodupli, int maxdata)
{
    //perhaps some datas are redondant..to verify
    circrads[0] = circrads[spottodupli] = dataspot[2][0] = dataspot[2][spottodupli] = dataspot[2][spottodupli - 1];
    locx[0] = locx[spottodupli] = dataspot[3][0] = dataspot[3][spottodupli] = dataspot[3][spottodupli - 1];
    locy[0] = locy[spottodupli] = dataspot[4][0] =  dataspot[4][spottodupli] = dataspot[4][spottodupli - 1];
    locyt[0] = locyt[spottodupli] = dataspot[5][0] = dataspot[5][spottodupli] = dataspot[5][spottodupli - 1];
    locxl[0] =  locxl[spottodupli] = dataspot[6][0] = dataspot[6][spottodupli] = dataspot[6][spottodupli - 1];
    //change only center position to 200 200 to see changes
    centerx[0] = centerx[spottodupli]  = dataspot[7][0] = dataspot[7][spottodupli] = 200; //not to center
    centery[0] = centery[spottodupli]  = dataspot[8][0] = dataspot[8][spottodupli] = 200;//not to center to see it is a duplicated spot
    //
    lights[0] = lights[spottodupli] = dataspot[9][0] = dataspot[9][spottodupli] = dataspot[9][spottodupli - 1];
    contrs[0] = contrs[spottodupli] = dataspot[10][0] = dataspot[10][spottodupli] = dataspot[10][spottodupli - 1];
    chroms[0] = chroms[spottodupli] = dataspot[11][0] = dataspot[11][spottodupli] = dataspot[11][spottodupli - 1];
    sensis[0] = sensis[spottodupli] = dataspot[12][0] = dataspot[12][spottodupli] = dataspot[12][spottodupli - 1];
    transits[0] = transits[spottodupli] = dataspot[13][0] = dataspot[13][spottodupli] = dataspot[13][spottodupli - 1];
    inverss[0] = inverss[spottodupli] = dataspot[14][0] = dataspot[14][spottodupli] = dataspot[14][spottodupli - 1];
    smeths[0]  = smeths[spottodupli] = dataspot[15][0] = dataspot[15][spottodupli] = dataspot[15][spottodupli - 1];
    //no chnage to spot current value  16
    radiuss[0]  = radiuss[spottodupli] = dataspot[17][0] = dataspot[17][spottodupli] = dataspot[17][spottodupli - 1];
    strengths[0]  = strengths[spottodupli] = dataspot[18][0] = dataspot[18][spottodupli] = dataspot[18][spottodupli - 1];
    sensibns[0]  = sensibns[spottodupli] = dataspot[19][0] = dataspot[19][spottodupli] = dataspot[19][spottodupli - 1];
    inversrads[0]  = inversrads[spottodupli] = dataspot[20][0] = dataspot[20][spottodupli] = dataspot[20][spottodupli - 1];
    strs[0]  = strs[spottodupli] = dataspot[21][0] = dataspot[21][spottodupli] = dataspot[21][spottodupli - 1];
    chrrts[0]  = chrrts[spottodupli] = dataspot[22][0] = dataspot[22][spottodupli] = dataspot[22][spottodupli - 1];
    neighs[0]  = neighs[spottodupli] = dataspot[23][0] = dataspot[23][spottodupli] = dataspot[23][spottodupli - 1];
    varts[0]  = varts[spottodupli] = dataspot[24][0] = dataspot[24][spottodupli] = dataspot[24][spottodupli - 1];
    sensihs[0]  = sensihs[spottodupli] = dataspot[25][0] = dataspot[25][spottodupli] = dataspot[25][spottodupli - 1];
    inversrets[0]  = inversrets[spottodupli] = dataspot[26][0] = dataspot[26][spottodupli] = dataspot[26][spottodupli - 1];
    retinexs[0]  = retinexs[spottodupli] = dataspot[27][0] = dataspot[27][spottodupli] = dataspot[27][spottodupli - 1];
    sharradiuss[0]  = sharradiuss[spottodupli] = dataspot[28][0] = dataspot[28][spottodupli] = dataspot[28][spottodupli - 1];
    sharamounts[0]  = sharamounts[spottodupli] = dataspot[29][0] = dataspot[29][spottodupli] = dataspot[29][spottodupli - 1];
    shardampings[0]  = shardampings[spottodupli] = dataspot[30][0] = dataspot[30][spottodupli] = dataspot[30][spottodupli - 1];
    shariters[0]  = shariters[spottodupli] = dataspot[31][0] = dataspot[31][spottodupli] = dataspot[31][spottodupli - 1];
    sensishas[0]  = sensishas[spottodupli] = dataspot[32][0] = dataspot[32][spottodupli] = dataspot[32][spottodupli - 1];
    inversshas[0]  = inversshas[spottodupli] = dataspot[33][0] = dataspot[33][spottodupli] = dataspot[33][spottodupli - 1];
    qualitys[0]  = qualitys[spottodupli] = dataspot[34][0] = dataspot[34][spottodupli] = dataspot[34][spottodupli - 1];
    thress[0]  = thress[spottodupli] = dataspot[35][0] = dataspot[35][spottodupli] = dataspot[35][spottodupli - 1];
    proxis[0]  = proxis[spottodupli] = dataspot[36][0] = dataspot[36][spottodupli] = dataspot[36][spottodupli - 1];
    noiselumfs[0]  = noiselumfs[spottodupli] = dataspot[37][0] = dataspot[37][spottodupli] = dataspot[37][spottodupli - 1];
    noiselumcs[0]  = noiselumcs[spottodupli] = dataspot[38][0] = dataspot[38][spottodupli] = dataspot[38][spottodupli - 1];
    noisechrofs[0]  = noisechrofs[spottodupli] = dataspot[39][0] = dataspot[39][spottodupli] = dataspot[39][spottodupli - 1];
    noisechrocs[0]  = noisechrocs[spottodupli] = dataspot[40][0] = dataspot[40][spottodupli] = dataspot[40][spottodupli - 1];
    mult0s[0]  = mult0s[spottodupli] = dataspot[41][0] = dataspot[41][spottodupli] = dataspot[41][spottodupli - 1];
    mult1s[0]  = mult1s[spottodupli] = dataspot[42][0] = dataspot[42][spottodupli] = dataspot[42][spottodupli - 1];
    mult2s[0]  = mult2s[spottodupli] = dataspot[43][0] = dataspot[43][spottodupli] = dataspot[43][spottodupli - 1];
    mult3s[0]  = mult3s[spottodupli] = dataspot[44][0] = dataspot[44][spottodupli] = dataspot[44][spottodupli - 1];
    mult4s[0]  = mult4s[spottodupli] = dataspot[45][0] = dataspot[45][spottodupli] = dataspot[45][spottodupli - 1];
    thresholds[0]  = thresholds[spottodupli] = dataspot[46][0] = dataspot[46][spottodupli] = dataspot[46][spottodupli - 1];
    sensicbs[0]  = sensicbs[spottodupli] = dataspot[47][0] = dataspot[47][spottodupli] = dataspot[47][spottodupli - 1];
    activlums[0]  = activlums[spottodupli] = dataspot[48][0] = dataspot[48][spottodupli] = dataspot[48][spottodupli - 1];
    strens[0]  = strens[spottodupli] = dataspot[49][0] = dataspot[49][spottodupli] = dataspot[49][spottodupli - 1];
    gammas[0]  = gammas[spottodupli] = dataspot[50][0] = dataspot[50][spottodupli] = dataspot[50][spottodupli - 1];
    estops[0]  = estops[spottodupli] = dataspot[51][0] = dataspot[51][spottodupli] = dataspot[51][spottodupli - 1];
    scaltms[0]  = scaltms[spottodupli] = dataspot[52][0] = dataspot[52][spottodupli] = dataspot[52][spottodupli - 1];
    reweis[0]  = reweis[spottodupli] = dataspot[53][0] = dataspot[53][spottodupli] = dataspot[53][spottodupli - 1];
    sensitms[0]  = sensitms[spottodupli] = dataspot[54][0] = dataspot[54][spottodupli] = dataspot[54][spottodupli - 1];
    retrabs[0]  = retrabs[spottodupli] = dataspot[55][0] = dataspot[55][spottodupli] = dataspot[55][spottodupli - 1];
    curvactivs[0]  = curvactivs[spottodupli] = dataspot[56][0] = dataspot[56][spottodupli] = dataspot[56][spottodupli - 1];
    qualitycurves[0]  = qualitycurves[spottodupli] = dataspot[57][0] = dataspot[57][spottodupli] = dataspot[57][spottodupli - 1];
    sensivs[0]  = sensivs[spottodupli] = dataspot[58][0] = dataspot[58][spottodupli] = dataspot[58][spottodupli - 1];
    pastels[0]  = pastels[spottodupli] = dataspot[59][0] = dataspot[59][spottodupli] = dataspot[59][spottodupli - 1];
    saturateds[0]  = saturateds[spottodupli] = dataspot[60][0] = dataspot[60][spottodupli] = dataspot[60][spottodupli - 1];
    protectskinss[0]  = protectskinss[spottodupli] = dataspot[61][0] = dataspot[61][spottodupli] = dataspot[61][spottodupli - 1];
    avoidcolorshifts[0]  = avoidcolorshifts[spottodupli] = dataspot[62][0] = dataspot[62][spottodupli] = dataspot[62][spottodupli - 1];
    pastsattogs[0]  = pastsattogs[spottodupli] = dataspot[63][0] = dataspot[63][spottodupli] = dataspot[63][spottodupli - 1];
    expcomps[0]  = expcomps[spottodupli] = dataspot[64][0] = dataspot[64][spottodupli] = dataspot[64][spottodupli - 1];
    blacks[0]  = blacks[spottodupli] = dataspot[65][0] = dataspot[65][spottodupli] = dataspot[65][spottodupli - 1];
    hlcomprs[0]  = hlcomprs[spottodupli] = dataspot[66][0] = dataspot[66][spottodupli] = dataspot[66][spottodupli - 1];
    hlcomprthreshs[0]  = hlcomprthreshs[spottodupli] = dataspot[67][0] = dataspot[67][spottodupli] = dataspot[67][spottodupli - 1];
    shcomprs[0]  = shcomprs[spottodupli] = dataspot[68][0] = dataspot[68][spottodupli] = dataspot[68][spottodupli - 1];
    sensiexs[0]  = sensiexs[spottodupli] = dataspot[69][0] = dataspot[69][spottodupli] = dataspot[69][spottodupli - 1];
    centerxbufs[0]  = centerxbufs[spottodupli] = dataspot[70][0] = dataspot[70][spottodupli] = dataspot[70][spottodupli - 1];
    centerybufs[0]  = centerybufs[spottodupli] = dataspot[71][0] = dataspot[71][spottodupli] = dataspot[71][spottodupli - 1];
    adjblurs[0]  = adjblurs[spottodupli] = dataspot[72][0] = dataspot[72][spottodupli] = dataspot[72][spottodupli - 1];
    cutpasts[0]  = cutpasts[spottodupli] = dataspot[73][0] = dataspot[73][spottodupli] = dataspot[73][spottodupli - 1];
    chromacbdls[0]  = chromacbdls[spottodupli] = dataspot[74][0] = dataspot[74][spottodupli] = dataspot[74][spottodupli - 1];
    lastdusts[0]  = lastdusts[spottodupli] = dataspot[75][0] = dataspot[75][spottodupli] = dataspot[75][spottodupli - 1];
    blurmets[0]  = blurmets[spottodupli] = dataspot[76][0] = dataspot[76][spottodupli] = dataspot[76][spottodupli - 1];
    dustmets[0]  = dustmets[spottodupli] = dataspot[77][0] = dataspot[77][spottodupli] = dataspot[77][spottodupli - 1];
    exclumets[0]  = exclumets[spottodupli] = dataspot[78][0] = dataspot[78][spottodupli] = dataspot[78][spottodupli - 1];
    sensiexclus[0]  = sensiexclus[spottodupli] = dataspot[79][0] = dataspot[79][spottodupli] = dataspot[79][spottodupli - 1];
    strucs[0]  = strucs[spottodupli] = dataspot[80][0] = dataspot[80][spottodupli] = dataspot[80][spottodupli - 1];
    warms[0]  = warms[spottodupli] = dataspot[81][0] = dataspot[81][spottodupli] = dataspot[81][spottodupli - 1];
    noiselumdetails[0]  = noiselumdetails[spottodupli] = dataspot[82][0] = dataspot[82][spottodupli] = dataspot[82][spottodupli - 1];
    noisechrodetails[0]  = noisechrodetails[spottodupli] = dataspot[83][0] = dataspot[83][spottodupli] = dataspot[83][spottodupli - 1];
    sensidens[0]  = sensidens[spottodupli] = dataspot[84][0] = dataspot[84][spottodupli] = dataspot[84][spottodupli - 1];
    expdenois[0]  = expdenois[spottodupli] = dataspot[85][0] = dataspot[85][spottodupli] = dataspot[85][spottodupli - 1];
    expcolors[0]  = expcolors[spottodupli] = dataspot[86][0] = dataspot[86][spottodupli] = dataspot[86][spottodupli - 1];
    expvibrances[0]  = expvibrances[spottodupli] = dataspot[87][0] = dataspot[87][spottodupli] = dataspot[87][spottodupli - 1];
    expblurs[0]  = expblurs[spottodupli] = dataspot[88][0] = dataspot[88][spottodupli] = dataspot[88][spottodupli - 1];
    exptonemaps[0]  = exptonemaps[spottodupli] = dataspot[89][0] = dataspot[89][spottodupli] = dataspot[89][spottodupli - 1];
    expretis[0]  = expretis[spottodupli] = dataspot[90][0] = dataspot[90][spottodupli] = dataspot[90][spottodupli - 1];
    expsharps[0]  = expsharps[spottodupli] = dataspot[91][0] = dataspot[91][spottodupli] = dataspot[91][spottodupli - 1];
    expcbdls[0]  = expcbdls[spottodupli] = dataspot[92][0] = dataspot[92][spottodupli] = dataspot[92][spottodupli - 1];
    expexposes[0]  = expexposes[spottodupli] = dataspot[93][0] = dataspot[93][spottodupli] = dataspot[93][spottodupli - 1];
    bilaterals[0]  = bilaterals[spottodupli] = dataspot[94][0] = dataspot[94][spottodupli] = dataspot[94][spottodupli - 1];
    noiselequals[0]  = noiselequals[spottodupli] = dataspot[95][0] = dataspot[95][spottodupli] = dataspot[95][spottodupli - 1];
    shapemets[0]  = shapemets[spottodupli] = dataspot[96][0] = dataspot[96][spottodupli] = dataspot[96][spottodupli - 1];

    //datas for end ... references hue, etc.
    huerefblurs[0] = huerefblurs[spottodupli] = dataspot[maxdata - 5][0] = dataspot[maxdata - 5][spottodupli] = dataspot[maxdata - 5][spottodupli - 1];
    huerefs[0] = huerefs[spottodupli] = dataspot[maxdata - 4][0] = dataspot[maxdata - 4][spottodupli] = dataspot[maxdata - 4][spottodupli - 1];
    chromarefs[0] = chromarefs[spottodupli] = dataspot[maxdata - 3][0] = dataspot[maxdata - 3][spottodupli] = dataspot[maxdata - 3][spottodupli - 1];
    lumarefs[0] = lumarefs[spottodupli] = dataspot[maxdata - 2][0] = dataspot[maxdata - 2][spottodupli] = dataspot[maxdata - 2][spottodupli - 1];
    sobelrefs[0] = sobelrefs[spottodupli] = dataspot[maxdata - 1][0] = dataspot[maxdata - 1][spottodupli] = dataspot[maxdata - 1][spottodupli - 1];

    //perhaps not good after ?? to verify and to complete ?? difficult but "only" curves
    retistr[spottodupli] = retistr[spottodupli - 1];
    llstr[spottodupli] = llstr[spottodupli - 1];
    lhstr[spottodupli] = lhstr[spottodupli - 1];
    ccstr[spottodupli] = ccstr[spottodupli - 1];
    hhstr[spottodupli] = hhstr[spottodupli - 1];
    skinstr[spottodupli] = skinstr[spottodupli - 1];
    pthstr[spottodupli] = pthstr[spottodupli - 1];
    exstr[spottodupli] = exstr[spottodupli - 1];

}

void ImProcCoordinator::changenumberofspot(int **dataspot, int maxdata, int maxspot, int ns, Glib::ustring datal, int versionmip)
{
    ofstream fic(datal, ios::out | ios::app);  // ouverture en écriture avec effacement du fichier ouvert


    for (int sp = ns + 1 ; sp < maxspot; sp++) { // spots default
        int t_sp = sp;
        int t_mipversion = versionmip;
        int t_circrad = 18;
        int t_locX = 250;
        int t_locY = 250;
        int t_locYT = 250;
        int t_locXL = 250;
        int t_centerX = 0;
        int t_centerY = 0;
        int t_lightness = 0;
        int t_contrast = 0;
        int t_chroma = 0;
        int t_sensi = 19;
        int t_transit = 60;
        int t_invers = 0;
        int t_Smeth = 0;
        int t_currentspot = 1;
        int t_radius = 1;
        int t_strength = 0;
        int t_sensibn = 40;
        int t_inversrad = 0;
        int t_str = 0;
        int t_chrrt = 0;
        int t_neigh = 50;
        int t_vart = 200;
        int t_sensih = 19;
        int t_inversret = 0;
        int t_retinexMethod = 2;
        int t_sharradius = 40;
        int t_sharamount = 75;
        int t_shardamping = 75;
        int t_shariter = 30;
        int t_sensisha = 19;
        int t_inverssha = 0;
        int t_qualityMethod = 1;
        int t_thres = 18;
        int t_proxi = 0;
        int t_noiselumf = 0;
        int t_noiselumc = 0;
        int t_noisechrof = 0;
        int t_noisechroc = 0;
        int t_mult0 = 100;
        int t_mult1 = 100;
        int t_mult2 = 100;
        int t_mult3 = 100;
        int t_mult4 = 100;
        int t_threshold = 20;
        int t_sensicb = 19;
        int t_activlum = 0;
        //10001 TM
        int t_stren = 0;
        int t_gamma = 100;
        int t_estop = 140;
        int t_scaltm = 10;
        int t_rewei = 0;
        int t_sensitm = 19;

        //10002 curve
        int t_retrab = 500;

        std::string t_curvret = "1000A0B120C350D350E700F500G350H350I1000J120K350L350M";//12 points
        //10003
        std::string t_curvll = "3000A0B0C1000D1000E"; //"3000A0B0C499D501E1000F1000G"; //"3000A0B0C1000D1000E";//0 points with marks
        //10004
        std::string t_curvlh = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y";
        //10005
        int t_curvactiv = 0;
        //10006
        std::string t_curvcc = "3000A0B0C1000D1000E";
        //10007
        int t_qualitycurveMethod = 0;
        //10008
        std::string t_curvhh = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y";

        //10009
        int t_sensiv = 19;
        int t_pastel = 0;
        int t_saturated = 0;
        std::string t_psthres = "0A75B";
        int t_proskin = 0;
        int t_avoidcsh = 0;
        int t_pastsat = 0;
        std::string t_curvskin = "3000A0B0C1000D1000E"; //"3000A0B0C499D501E1000F1000G"; //"3000A0B0C1000D1000E";//0 points with marks

        int t_expcomp       = 0;
        int t_black         = 0;
        int t_hlcompr       = 20;
        int t_hlcomprthresh = 33;
        int t_shcompr       = 50;
        int t_sensiex = 19;
        //10010
        std::string t_curvex = "3000A0B0C1000D1000E";

        //10012
        int t_centerXbuf = 0;
        int t_centerYbuf = 0;
        int t_adjblur = 0;
        int t_cutpast = 0;

        //10013
        int t_chromacbdl = 0;

        //10014
        int t_lastdust = 0;
        int t_blurMethod = 0;
        int t_dustMethod = 1;

        //10016
        int t_excludemeth = 0;
        int t_sensiexclu = 19;
        int t_struc = 0;

        //10017
        int t_warm = 0;
        //10018
        int t_noiselumdetail = 0;
        //10019
        int t_noisechrodetail = 0;
        //10020
        int t_sensiden = 30;

        //10021
        int t_expdenoi = 0;

        int t_expcolor = 0;
        int t_expvibrance = 0;
        int t_expblur = 0;
        int t_exptonemap = 0;
        int t_expreti = 0;
        int t_expsharp = 0;
        int t_expcbdl = 0;
        int t_expexpose = 0;

        //10022
        int t_bilateral = 0;

        //10023
        int t_noiselequal = 7;
        //10024
        int t_shapemeth = 0;

        fic << "Mipversion=" << t_mipversion << '@' << endl;
        fic << "Spot=" << t_sp << '@' << endl;
        fic << "Circrad=" << t_circrad << '@' << endl;
        fic << "LocX=" << t_locX << '@' << endl;
        fic << "LocY=" << t_locY << '@' << endl;
        fic << "LocYT=" << t_locYT << '@' << endl;
        fic << "LocXL=" << t_locXL << '@' << endl ;
        fic << "CenterX=" << t_centerX << '@' << endl;
        fic << "CenterY=" << t_centerY << '@' << endl;
        fic << "Lightness=" << t_lightness << '@' << endl;
        fic << "Contrast=" << t_contrast << '@' <<  endl;
        fic << "Chroma=" << t_chroma << '@' << endl;
        fic << "Sensi=" << t_sensi << '@' << endl;
        fic << "Transit=" << t_transit << '@' << endl;
        fic << "Invers=" << t_invers << '@' << endl;
        fic << "Smethod=" << t_Smeth << '@' << endl;
        fic << "Currentspot=" << t_currentspot << '@' << endl;
        fic << "Radius=" << t_radius << '@' << endl;
        fic << "Strength=" << t_strength << '@' << endl;
        fic << "Sensibn=" << t_sensibn << '@' << endl;
        fic << "Inversrad=" << t_inversrad << '@' << endl;
        fic << "Str=" << t_str << '@' << endl;
        fic << "Chroma=" << t_chrrt << '@' << endl;
        fic << "Neigh=" << t_neigh << '@' << endl;
        fic << "Vart=" << t_vart << '@' << endl;
        fic << "Sensih=" << t_sensih << '@' << endl;
        fic << "Inversret=" << t_inversret << '@' << endl;
        fic << "retinexMethod=" << t_retinexMethod << '@' << endl;
        fic << "Sharradius=" << t_sharradius << '@' << endl;
        fic << "Sharamount=" << t_sharamount << '@' << endl;
        fic << "Shardamping=" << t_shardamping << '@' << endl;
        fic << "Shariter=" << t_shariter << '@' << endl;
        fic << "Sensisha=" << t_sensisha << '@' << endl;
        fic << "Inverssha=" << t_inverssha << '@' << endl;
        fic << "qualityMethod=" << t_qualityMethod << '@' << endl;
        fic << "Thres=" << t_thres << '@' << endl;
        fic << "Proxi=" << t_proxi << '@' << endl;
        fic << "Noiselumf=" << t_noiselumf << '@' << endl;
        fic << "Noiselumc=" << t_noiselumc << '@' << endl;
        fic << "Noisechrof=" << t_noisechrof << '@' << endl;
        fic << "Noisechroc=" << t_noisechroc << '@' << endl;
        fic << "Mult0=" << t_mult0 << '@' << endl;
        fic << "Mult1=" << t_mult1 << '@' << endl;
        fic << "Mult2=" << t_mult2 << '@' << endl;
        fic << "Mult3=" << t_mult3 << '@' << endl;
        fic << "Mult4=" << t_mult4 << '@' << endl;
        fic << "Threshold=" << t_threshold << '@' << endl;
        fic << "Sensicb=" << t_sensicb << '@' << endl;
        fic << "Activblurlum=" << t_activlum << '@' << endl;

        fic << "Stren=" << t_stren << '@' << endl;
        fic << "Gamma=" << t_gamma << '@' << endl;
        fic << "Estop=" << t_estop << '@' << endl;
        fic << "Scaltm=" << t_scaltm << '@' << endl;
        fic << "Rewei=" << t_rewei << '@' << endl;
        fic << "Sensitm=" << t_sensitm << '@' << endl;
        fic << "Retrab=" << t_retrab << '@' << endl;
        fic << "Curvactiv=" << t_curvactiv << '@' << endl;
        fic << "qualitycurveMethod=" << t_qualitycurveMethod << '@' << endl;

        fic << "Sensiv=" << t_sensiv << '@' << endl;
        fic << "Pastel=" << t_pastel << '@' << endl;
        fic << "Saturated=" << t_saturated << '@' << endl;
        fic << "Proskin=" << t_proskin << '@' << endl;
        fic << "Avoidcsh=" << t_avoidcsh << '@' << endl;
        fic << "Pastsat=" << t_pastsat << '@' << endl;

        fic << "Expcomp=" << t_expcomp << '@' << endl;
        fic << "Black=" << t_black << '@' << endl;
        fic << "Hlcompr=" << t_hlcompr << '@' << endl;
        fic << "Hlcomprthresh=" << t_hlcomprthresh << '@' << endl;
        fic << "Shcompr=" << t_shcompr  << '@' << endl;
        fic << "Sensiex=" << t_sensiex << '@' << endl;

        fic << "CenterXbuf=" << t_centerXbuf << '@' << endl;
        fic << "CenterYbuf=" << t_centerYbuf << '@' << endl;
        fic << "Adjblur=" << t_adjblur << '@' << endl;
        fic << "Cutpast=" << t_cutpast << '@' <<  endl;

        fic << "Chromacbdl=" << t_chromacbdl << '@' <<  endl;
        fic << "Lastdust=" << t_lastdust << '@' <<  endl;
        fic << "BlurMethod=" << t_blurMethod << '@' <<  endl;
        fic << "DustMethod=" << t_dustMethod << '@' <<  endl;

        fic << "ExcludeMethod=" << t_excludemeth << '@' <<  endl;
        fic << "Sensiexclu=" << t_sensiexclu << '@' << endl;
        fic << "Struc=" << t_struc << '@' << endl;
        fic << "Warm=" << t_warm << '@' << endl;
        fic << "Noiselumdetail=" << t_noiselumdetail << '@' << endl;
        fic << "Noisechrodetail=" << t_noisechrodetail << '@' << endl;

        fic << "Sensiden=" << t_sensiden << '@' << endl;
        fic << "Expdenoi=" << t_expdenoi << '@' << endl;
        fic << "Expcolor=" << t_expcolor << '@' << endl;
        fic << "Expvibrance=" << t_expvibrance << '@' << endl;
        fic << "Expblur=" << t_expblur << '@' << endl;
        fic << "Exptonemap=" << t_exptonemap << '@' << endl;
        fic << "Expreti=" << t_expreti << '@' << endl;
        fic << "Expsharp=" << t_expsharp << '@' << endl;
        fic << "Expcbdl=" << t_expcbdl << '@' << endl;
        fic << "Expexpose=" << t_expexpose << '@' << endl;

        fic << "Bilateral=" << t_bilateral << '@' << endl;
        fic << "Noiselequal=" << t_noiselequal << '@' << endl;
        fic << "ShapeMethod=" << t_shapemeth << '@' <<  endl;

        fic << "curveReti=" << t_curvret << '@' << endl;
        fic << "curveLL=" << t_curvll << '@' << endl;
        fic << "curveLH=" << t_curvlh << '@' << endl;
        fic << "curveCC=" << t_curvcc << '@' << endl;
        fic << "curveHH=" << t_curvhh << '@' << endl;
        fic << "curveskin=" << t_curvskin << '@' << endl;
        fic << "pthres=" << t_psthres << '@' << endl;
        fic << "curveex=" << t_curvex << '@' << endl;



        fic << endl;
    }

    fic.close();

    ifstream fich2(datal, ios::in);

    if (fich2) {

        std::string line2;
        std::string spotline2;
        int cont2 = 0;
        int ns2 = 0;
        int maxin = maxdata - 5; //70 ;//64

        while (getline(fich2, line2)) {
            spotline2 = line2;
            std::size_t pos2 = spotline2.find("=");
            std::size_t posend2 = spotline2.find("@");  //in case of for futur use

            if (spotline2.substr(0, pos2) == "Mipversion") {
                std::string strversion = spotline2.substr(pos2 + 1, (posend2 - pos2));
                versionmip = std::stoi(strversion.c_str());
            }

            if (spotline2.substr(0, pos2) == "Spot") {
                cont2 = 0;
            }

            cont2++;
            std::string str32 = spotline2.substr(pos2 + 1, (posend2 - pos2));

            if (cont2 == 1) {
                ns2 =  std::stoi(str32.c_str());
            }

            if (cont2 >= 2  && cont2 < 16) {
                dataspot[cont2][ns2] = std::stoi(str32.c_str());
            }

            if (spotline2.substr(0, pos2) == "Currentspot") {
                dataspot[16][0] = std::stoi(str32.c_str());
            }

            if (cont2 > 16  && cont2 < maxin) {
                dataspot[cont2][ns2] = std::stoi(str32.c_str());
            }

            if (spotline2.substr(0, pos2) == "curveReti") {
                retistr[ns2] = str32;
            }

            if (spotline2.substr(0, pos2) == "curveLL") {
                llstr[ns2] = str32;
            }

            if (spotline2.substr(0, pos2) == "curveLH") {
                lhstr[ns2] = str32;
            }

            if (spotline2.substr(0, pos2) == "curveCC") {
                ccstr[ns2] = str32;
            }

            if (spotline2.substr(0, pos2) == "curveHH") {
                hhstr[ns2] = str32;
            }

            if (spotline2.substr(0, pos2) == "curveskin") {
                skinstr[ns2] = str32;
            }

            if (spotline2.substr(0, pos2) == "pthres") {
                pthstr[ns2] = str32;
            }

            if (spotline2.substr(0, pos2) == "curveex") {
                exstr[ns2] = str32;
            }

        }

        fich2.close() ;
    }
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
