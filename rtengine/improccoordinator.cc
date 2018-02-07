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
      ncie(nullptr), imgsrc(nullptr), shmap(nullptr), lastAwbEqual(0.), lastAwbTempBias(0.0), ipf(&params, true), monitorIntent(RI_RELATIVE),
      softProof(false), gamutCheck(false), scale(10), highDetailPreprocessComputed(false), highDetailRawComputed(false),
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
      hltonecurveloc(32768, 0), //32768
      shtonecurveloc(32768, 0),
      tonecurveloc(32768, 0),
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
      plistener(nullptr), awbListener(nullptr), imageListener(nullptr), aeListener(nullptr), acListener(nullptr), abwListener(nullptr),  aloListener(nullptr), actListener(nullptr), adnListener(nullptr), awavListener(nullptr), dehaListener(nullptr), frameCountListener(nullptr), imageTypeListener(nullptr), hListener(nullptr),
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
      colourToningSatLimit(0.f), colourToningSatLimitOpacity(0.f), lastspotdup(false), highQualityComputed (false), 

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
        imageTypeListener->imageTypeChanged(imgsrc->isRAW(), imgsrc->getSensorType() == ST_BAYER, imgsrc->getSensorType() == ST_FUJI_XTRANS);
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
                            imgsrc->convertColorSpace(calclum, params.icm, currWB);//claculate values after colorspace conversion
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

    if ((todo & M_BLURMAP) && params.sh.enabled) {
        double radius = sqrt(double (pW * pW + pH * pH)) / 2.0;
        double shradius = params.sh.radius;

        if (!params.sh.hq) {
            shradius *= radius / 1800.0;
        }

        if (!shmap) {
            shmap = new SHMap(pW, pH, true);
        }

        shmap->update(oprevi, shradius, ipf.lumimul, params.sh.hq, scale);
    }



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

            ipf.rgbProc(oprevi, oprevl, nullptr, hltonecurve, shtonecurve, tonecurve, shmap, params.toneCurve.saturation,
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

        int maxspot = settings->nspot + 1;
        progress("Applying Color Boost...", 100 * readyphase / numofphases);


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
             */


            //*********************************************************
            //advertissment
            //we can probably put all these function outside main process
            // but for now, I think it s mode readable
            // we have similar process in dcrop. and simpleprocess.cc
            // see rawpedia for all "fantaisies"
            // all this code is probably not optimal...but actually it run :)
            //there are probably errors...
            //***********************************************************

            bool isascii = true;
            std::string mdfive = getMD5(imgsrc->getFileName());


            Glib::ustring datainterm = imgsrc->getFileName() + ".ii";//extansion ii arbitrary to test if mip file is possible

            ofstream finterm(datainterm, ios::out);

            if (finterm.fail()) {
                printf("Non ascii Mip file possible..switch to Profiles\n");
                isascii = false;
            } else {
                printf("ascii Mip file possible!\n");
            }

            finterm.close();

            if (isascii == true) {
                if (std::remove(datainterm.c_str()) != 0) {
                    perror("Error deleting test ii file");
                } else {
                    puts("Test ii file successfully deleted");
                }
            }

            Glib::ustring pop = options.cacheBaseDir + "/mip/";

            Glib::ustring datal;

            if (options.mip == MI_opt || !isascii) {
                datal = pop + Glib::path_get_basename(imgsrc->getFileName() + "." + mdfive + ".mip");
            }

            if (options.mip == MI_prev && isascii) {//&& isascii
                datal = imgsrc->getFileName() + ".mip";
            }

            /*
            //test to see if wofstream and wifstream works with NON ASCII, but it's bad
                        wofstream test(datal, ios::out);
                        if(test.fail()) printf("ca va pas\n");
                        else ("ca va bien\n");
                        test.close();
            */
            ifstream fic0(datal, ios::in);

            printf("mip files in=%s\n", datal.c_str());
            //    if(! fic0.fail())    {
            float **shbuffer = nullptr;
            versionmip = 0;
            //   int maxdat;
            int sca = 1;
            //string delim ==> delimiter to separate integer in a string, 70 is largely enough for curves : noramlly 3 to 21 must be suffisant
            //curious method, but I am a poor informatic scientist...
            std::string delim[69] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
                                     "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
                                     "&", "#", "{", "[", "]", "}", "$", "*", "?", ">", "!", ";", "<", "(", ")", "+", "-"
                                    };


            maxdata = 102; //101 10023 //100 10022 //99 10021 // 90 10020 //88 10019//87 10018  //86 10017 //85 10016;// 82 10015//78;//73 for 10011
            //same value in simpleprocess.cc
            //same value in locallab.h for int nextdatasp[102];//102 = maxdata

            if (fic0) {
                //find current version mip
                std::string line;
                std::string spotline;
                //       int cont = 0;

                while (getline(fic0, line)) {
                    spotline = line;
                    std::size_t pos = spotline.find("=");
                    std::size_t posend = spotline.find("@");  //in case of for futur use

                    if (spotline.substr(0, pos) == "Mipversion") {
                        std::string strversion = spotline.substr(pos + 1, (posend - pos));
                        versionmip = std::stoi(strversion.c_str());
                    }


                }

                fic0.close();
            }

            printf("current mipvers=%i\n", versionmip);
            ifstream fic(datal, ios::in);



            if (fic.fail() || versionmip == 0  || params.locallab.nbspot == 0) { //initialize mip with default values if no file or old file to prevent crash

                ofstream fic(datal, ios::out | ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

                if (params.locallab.nbspot == 0) {
                    params.locallab.nbspot = 1;
                }


                if (fic) {
                    mipver = 10024;//to actualize  for each change, must be change also at the end just before save all datas mip files

                    //***************************************************************************************
                    //initialize new values when first utilisation of Locallab. Prepare creation of Mip files
                    //****************************************************************************************
                    for (int sp = 1; sp < maxspot; sp++) { // spots default
                        int t_sp = sp;
                        int t_mipversion = mipver;//new value for each change neads here, if it is first use
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

                        // end versionmip = 10000

                        //begin versionmip = 10001 Tone mapping
                        int t_stren = 0;
                        int t_gamma = 100;
                        int t_estop = 140;
                        int t_scaltm = 10;
                        int t_rewei = 0;
                        int t_sensitm = 19;

                        //versionmip = 10002 Reticurv
                        int t_retrab = 500;

                        std::string t_curvret = "1000A0B120C350D350E700F500G350H350I1000J120K350L350M";//12 points
                        //10003
                        //std::string t_curvll = "0A";
                        std::string t_curvll = "3000A0B0C1000D1000E"; //"3000A0B0C499D501E1000F1000G";// "3000A0B0C1000D1000E";//with that it works !

                        //versionmip = 10004 LHcurv
                        std::string t_curvlh = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y";
                        //10005
                        int t_curvactiv = 0;
                        //10006
                        std::string t_curvcc = "3000A0B0C1000D1000E"; //"3000A0B0C499D501E1000F1000G";// "3000A0B0C1000D1000E";//with that it works !

                        //10007
                        int t_qualitycurveMethod = 0;
                        //versionmip = 10008 HHcurv
                        std::string t_curvhh = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y";

                        //10009
                        int t_sensiv = 19;
                        int t_pastel = 0;
                        int t_saturated = 0;
                        std::string t_psthres = "0A75B";
                        int t_proskin = 0;
                        int t_avoidcsh = 0;
                        int t_pastsat = 0;
                        std::string t_curvskin = "3000A0B0C1000D1000E"; //"3000A0B0C499D501E1000F1000G";// "3000A0B0C1000D1000E";//with that it works !

                        //10010
                        int t_expcomp       = 0;
                        int t_black         = 0;
                        int t_hlcompr       = 20;
                        int t_hlcomprthresh = 33;
                        int t_shcompr       = 50;
                        int t_sensiex = 19;

                        //10011
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
                        //10019
                        int t_sensiden = 30;

                        //10021
                        int t_expdenoi = 0;

                        //10022
                        int t_bilateral = 0;

                        int t_expcolor = 0;
                        int t_expvibrance = 0;
                        int t_expblur = 0;
                        int t_exptonemap = 0;
                        int t_expreti = 0;
                        int t_expsharp = 0;
                        int t_expcbdl = 0;
                        int t_expexpose = 0;
                        //10023
                        int t_noiselequal = 7;

                        //10024
                        int t_shapemeth = 0;

                        //all variables except locRETgainCurve 'coomon for all)
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

                } else

                {
                    cerr << "can't open file !" << endl;
                }

            }

            //***************************************************************************************
            //End initialize new values when first utilisation of Locallab. Prepare creation of Mip files
            //****************************************************************************************


            int realspot = params.locallab.nbspot;

            if (realspot >= maxspot) {
                params.locallab.nbspot = realspot = 1;
            }

            std::string inser;

            //create data for mip files
            dataspot = new int*[maxdata];

            //initilize data "0" with params
            for (int i = 0; i < maxdata; i++) {
                dataspot[i] = new int[maxspot];
            }

            retistr = new std::string[maxspot];
            llstr = new std::string[maxspot];
            lhstr = new std::string[maxspot];
            ccstr = new std::string[maxspot];
            hhstr = new std::string[maxspot];
            skinstr = new std::string[maxspot];
            pthstr = new std::string[maxspot];
            exstr = new std::string[maxspot];

            //******************************************************************
            //initialize data[xx][0] and Lut cache with params
            //******************************************************************
            {
                sps[0] = 0;
                dataspot[2][0] =  circrads[0] = params.locallab.circrad;//copy params in dataspot and in LUTi
                dataspot[3][0] =  locx[0] = params.locallab.locX;
                dataspot[4][0] =  locy[0] = params.locallab.locY;
                dataspot[5][0] =  locyt[0] = params.locallab.locYT;
                dataspot[6][0] =  locxl[0] = params.locallab.locXL;
                dataspot[7][0] =  centerx[0] = params.locallab.centerX;
                dataspot[8][0] =  centery[0] = params.locallab.centerY;
                dataspot[9][0] =  lights[0] = params.locallab.lightness;
                dataspot[10][0] =  contrs[0] = params.locallab.contrast;
                dataspot[11][0] =  chroms[0] = params.locallab.chroma;
                dataspot[12][0] =  sensis[0] = params.locallab.sensi;
                dataspot[13][0] =  transits[0] = params. locallab.transit;

                if (!params.locallab.invers) {
                    dataspot[14][0] = inverss[0] = 0;
                } else {
                    dataspot[14][0] = inverss[0] = 1;
                }

                if (params.locallab.Smethod == "IND") {
                    dataspot[15][0] = smeths[0] =  0;
                } else if (params.locallab.Smethod == "SYM") {
                    dataspot[15][0] = smeths[0] =  1;
                } else if (params.locallab.Smethod == "INDSL") {
                    dataspot[15][0] = smeths[0] =  2;
                } else if (params.locallab.Smethod == "SYMSL") {
                    dataspot[15][0] =  smeths[0] = 3;
                }

                dataspot[16][0] = curens[0] = params.locallab.nbspot;
                dataspot[17][0] =  radiuss[0] = params.locallab.radius;
                dataspot[18][0] =  strengths[0] = params.locallab.strength;
                dataspot[19][0] =  sensibns[0] = params.locallab.sensibn;


                if (!params.locallab.inversrad) {
                    dataspot[20][0] =  inversrads[0] = 0;
                } else {
                    dataspot[20][0] =  inversrads[0] = 1;
                }

                dataspot[21][0] = strs[0] = params.locallab.str;
                dataspot[22][0] = chrrts[0] = params.locallab.chrrt;
                dataspot[23][0] = neighs[0] = params.locallab.neigh;
                dataspot[24][0] = varts[0] = params.locallab.vart;
                dataspot[25][0] = sensihs[0] = params.locallab.sensih;

                if (!params.locallab.inversret) {
                    dataspot[26][0] =  inversrets[0] = 0;
                } else {
                    dataspot[26][0] =  inversrets[0] = 1;
                }

                if (params.locallab.retinexMethod == "low") {
                    dataspot[27][0] =  retinexs[0] = 0;
                } else if (params.locallab.retinexMethod == "uni") {
                    dataspot[27][0] =  retinexs[0] = 1;
                } else if (params.locallab.retinexMethod == "high") {
                    dataspot[27][0] =  retinexs[0] = 2;
                }

                dataspot[28][0] = sharradiuss[0] = params.locallab.sharradius;
                dataspot[29][0] = sharamounts[0] = params.locallab.sharamount;
                dataspot[30][0] = shardampings[0] = params.locallab.shardamping;
                dataspot[31][0] = shariters[0] = params.locallab.shariter;
                dataspot[32][0] = sensishas[0] = params.locallab.sensisha;

                if (!params.locallab.inverssha) {
                    dataspot[33][0] =  inversshas[0] = 0;
                } else {
                    dataspot[33][0] =  inversshas[0] = 1;
                }

                if (params.locallab.qualityMethod == "std") {
                    dataspot[34][0] =  qualitys[0] = 0;
                } else if (params.locallab.qualityMethod == "enh") {
                    dataspot[34][0] =  qualitys[0] = 1;
                } else if (params.locallab.qualityMethod == "enhden") {
                    dataspot[34][0] =  qualitys[0] = 2;
                }

                dataspot[35][0] = thress[0] = params.locallab.thres;
                dataspot[36][0] = proxis[0] = params.locallab.proxi;
                dataspot[37][0] = noiselumfs[0] = params.locallab.noiselumf;
                dataspot[38][0] = noiselumcs[0] = params.locallab.noiselumc;
                dataspot[39][0] = noisechrofs[0] = params.locallab.noisechrof;
                dataspot[40][0] = noisechrocs[0] = params.locallab.noisechroc;

                dataspot[41][0] = mult0s[0] = params.locallab.mult[0];
                dataspot[42][0] = mult1s[0] = params.locallab.mult[1];
                dataspot[43][0] = mult2s[0] = params.locallab.mult[2];
                dataspot[44][0] = mult3s[0] = params.locallab.mult[3];
                dataspot[45][0] = mult4s[0] = params.locallab.mult[4];
                dataspot[46][0] = thresholds[0] = params.locallab.threshold;
                dataspot[47][0] = sensicbs[0] = params.locallab.sensicb;

                if (!params.locallab.activlum) {
                    dataspot[48][0] =  activlums[0] = 0;
                } else {
                    dataspot[48][0] =  activlums[0] = 1;
                }

                dataspot[49][0] = strens[0] = params.locallab.stren;
                dataspot[50][0] = gammas[0] = params.locallab.gamma;
                dataspot[51][0] = estops[0] = params.locallab.estop;
                dataspot[52][0] = scaltms[0] = params.locallab.scaltm;
                dataspot[53][0] = reweis[0] = params.locallab.rewei;
                dataspot[54][0] = sensitms[0] = params.locallab.sensitm;
                dataspot[55][0] = retrabs[0] = params.locallab.retrab;

                if (!params.locallab.curvactiv) {
                    dataspot[56][0] = curvactivs[0] = 0;
                } else {
                    dataspot[56][0] = curvactivs[0] = 1;
                }

                if (params.locallab.qualitycurveMethod == "none") {
                    dataspot[57][0] =  qualitycurves[0] = 0;
                } else if (params.locallab.qualitycurveMethod == "std") {
                    dataspot[57][0] =  qualitycurves[0] = 1;
                } else if (params.locallab.qualitycurveMethod == "enh") {
                    dataspot[57][0] =  qualitycurves[0] = 2;
                }


                dataspot[58][0] = sensivs[0] = params.locallab.sensiv;
                dataspot[59][0] = pastels[0] = params.locallab.pastels;
                dataspot[60][0] = saturateds[0] = params.locallab.saturated;

                if (!params.locallab.protectskins) {
                    dataspot[61][0] = protectskinss[0] = 0;
                } else {
                    dataspot[61][0] = protectskinss[0] = 1;
                }

                if (!params.locallab.avoidcolorshift) {
                    dataspot[62][0] = avoidcolorshifts[0] = 0;
                } else {
                    dataspot[62][0] = avoidcolorshifts[0] = 1;
                }

                if (!params.locallab.pastsattog) {
                    dataspot[63][0] = pastsattogs[0] = 0;
                } else {
                    dataspot[63][0] = pastsattogs[0] = 1;
                }

                dataspot[64][0] = expcomps[0] = params.locallab.expcomp;
                dataspot[65][0] = blacks[0] = params.locallab.black;
                dataspot[66][0] = hlcomprs[0] = params.locallab.hlcompr;
                dataspot[67][0] = hlcomprthreshs[0] = params.locallab.hlcomprthresh;
                dataspot[68][0] = shcomprs[0] = params.locallab.shcompr;
                dataspot[69][0] = sensiexs[0] = params.locallab.sensiex;

                dataspot[70][0] = centerxbufs[0] = params.locallab.centerXbuf;
                dataspot[71][0] = centerybufs[0] = params.locallab.centerYbuf;
                dataspot[72][0] = adjblurs[0] = params.locallab.adjblur;

                if (!params.locallab.cutpast) {
                    dataspot[73][0] = cutpasts[0] = 0;
                } else {
                    dataspot[73][0] = cutpasts[0] = 1;
                }


                dataspot[74][0] = chromacbdls[0] = params.locallab.chromacbdl;

                if (!params.locallab.lastdust) {
                    dataspot[75][0] = lastdusts[0] = 0;
                } else {
                    dataspot[75][0] = lastdusts[0] = 1;
                }

                if (params.locallab.blurMethod == "norm") {
                    dataspot[76][0] =  blurmets[0] = 0;
                } else if (params.locallab.blurMethod == "inv") {
                    dataspot[76][0] =  blurmets[0] = 1;
                } else if (params.locallab.blurMethod == "sym") {
                    dataspot[76][0] =  blurmets[0] = 2;
                }

                if (params.locallab.dustMethod == "cop") {
                    dataspot[77][0] =  dustmets[0] = 0;
                } else if (params.locallab.dustMethod == "mov") {
                    dataspot[77][0] =  dustmets[0] = 1;
                } else if (params.locallab.dustMethod == "pas") {
                    dataspot[77][0] =  dustmets[0] = 2;
                }

                if (params.locallab.Exclumethod == "norm") {
                    dataspot[78][0] =  exclumets[0] = 0;
                } else if (params.locallab.Exclumethod == "exc") {
                    dataspot[78][0] =  exclumets[0] = 1;
                }

                dataspot[79][0] = sensiexclus[0] = params.locallab.sensiexclu;
                dataspot[80][0] = strucs[0] = params.locallab.struc;
                dataspot[81][0] = warms[0] = params.locallab.warm;
                dataspot[82][0] = noiselumdetails[0] = params.locallab.noiselumdetail;
                dataspot[83][0] = noisechrodetails[0] = params.locallab.noisechrodetail;
                dataspot[84][0] = sensidens[0] = params.locallab.sensiden;

                if (!params.locallab.expdenoi) {
                    dataspot[85][0] = expdenois[0] = 0;
                } else {
                    dataspot[85][0] = expdenois[0] = 1;
                }

                if (!params.locallab.expcolor) {
                    dataspot[86][0] = expcolors[0] = 0;
                } else {
                    dataspot[86][0] = expcolors[0] = 1;
                }

                if (!params.locallab.expvibrance) {
                    dataspot[87][0] = expvibrances[0] = 0;
                } else {
                    dataspot[87][0] = expvibrances[0] = 1;
                }

                if (!params.locallab.expblur) {
                    dataspot[88][0] = expblurs[0] = 0;
                } else {
                    dataspot[88][0] = expblurs[0] = 1;
                }

                if (!params.locallab.exptonemap) {
                    dataspot[89][0] = exptonemaps[0] = 0;
                } else {
                    dataspot[89][0] = exptonemaps[0] = 1;
                }

                if (!params.locallab.expreti) {
                    dataspot[90][0] = expretis[0] = 0;
                } else {
                    dataspot[90][0] = expretis[0] = 1;
                }

                if (!params.locallab.expsharp) {
                    dataspot[91][0] = expsharps[0] = 0;
                } else {
                    dataspot[91][0] = expsharps[0] = 1;
                }

                if (!params.locallab.expcbdl) {
                    dataspot[92][0] = expcbdls[0] = 0;
                } else {
                    dataspot[92][0] = expcbdls[0] = 1;
                }

                if (!params.locallab.expexpose) {
                    dataspot[93][0] = expexposes[0] = 0;
                } else {
                    dataspot[93][0] = expexposes[0] = 1;
                }

                dataspot[94][0] = bilaterals[0] = params.locallab.bilateral;
                dataspot[95][0] = noiselequals[0] = params.locallab.noiselequal;

                if (params.locallab.shapemethod == "ELI") {
                    dataspot[96][0] =  shapemets[0] = 0;
                } else if (params.locallab.shapemethod == "RECT") {
                    dataspot[96][0] =  shapemets[0] = 1;
                }

                // for all curves work around - I do not know how to do with params curves...
                //curve Reti local
                int siz = params.locallab.localTgaincurve.size();

                if (siz > 69) {//max due to codage with strcurv_data ()
                    siz = 69;    //to avoid crash
                }

                int s_datcur[siz + 1];

                for (int j = 0; j < siz; j++) {
                    s_datcur[j] = reticurvs[0 + j] = (int)(1000. * params.locallab.localTgaincurve[j]);
                }

                std::string cur_str = "";

                for (int j = 0; j < siz; j++) {
                    cur_str = cur_str + std::to_string(s_datcur[j]) + delim[j];
                }

                inser = retistr[0] = cur_str + "@";
                //end retistr

                //curve local L Lum
                int sizl = params.locallab.llcurve.size();

                if (sizl > 69) {
                    sizl = 69;//to avoid crash
                }

                int s_datcurl[sizl + 1];

                for (int j = 0; j < sizl; j++) {
                    s_datcurl[j] = llcurvs[0 + j] = (int)(1000. * params.locallab.llcurve[j]);
                }

                std::string ll_str = "";

                for (int j = 0; j < sizl; j++) {
                    ll_str = ll_str + std::to_string(s_datcurl[j]) + delim[j];
                }

                llstr[0] = ll_str + "@";
                //end local L f(L)

                //curve local C chrom
                int sizc = params.locallab.cccurve.size();

                if (sizc > 69) {//max
                    sizc = 69;//to avoid crash
                }

                int s_datcurc[sizc + 1];

                for (int j = 0; j < sizc; j++) {
                    s_datcurc[j] = cccurvs[0 + j] = (int)(1000. * params.locallab.cccurve[j]);
                }

                std::string cc_str = "";

                for (int j = 0; j < sizc; j++) {
                    cc_str = cc_str + std::to_string(s_datcurc[j]) + delim[j];
                }

                ccstr[0] = cc_str + "@";
                //end local C f(C)


                //curve local L f(H)
                int sizh = params.locallab.LHcurve.size();

                if (sizh > 69) {
                    sizh = 69;//to avoid crash
                }

                //     int s_curh[sizh + 1];
                int s_datcurh[sizh + 1];

                for (int j = 0; j < sizh; j++) {
                    s_datcurh[j] = lhcurvs[0 + j] = (int)(1000. * params.locallab.LHcurve[j]);
                }

                std::string lh_str = "";

                for (int j = 0; j < sizh; j++) {
                    lh_str = lh_str + std::to_string(s_datcurh[j]) + delim[j];
                }

                lhstr[0] = lh_str + "@";


                //HH curve
                //curve local H f(H)
                int sizhh = params.locallab.HHcurve.size();

                if (sizhh > 69) {
                    sizhh = 69;//to avoid crash
                }


                int s_datcurhh[sizhh + 1];

                for (int j = 0; j < sizhh; j++) {
                    s_datcurhh[j] = hhcurvs[0 + j] = (int)(1000. * params.locallab.HHcurve[j]);
                }

                std::string hh_str = "";

                for (int j = 0; j < sizhh; j++) {
                    hh_str = hh_str + std::to_string(s_datcurhh[j]) + delim[j];
                }

                hhstr[0] = hh_str + "@";

                //end local L = f(H)

                //Skin curve
                int sizsk = params.locallab.skintonescurve.size();

                if (sizsk > 69) {
                    sizsk = 69;//to avoid crash
                }


                int s_datcursk[sizsk + 1];

                for (int j = 0; j < sizsk; j++) {
                    s_datcursk[j] = skintonescurves[0 + j] = (int)(1000. * params.locallab.skintonescurve[j]);
                }

                std::string sk_str = "";

                for (int j = 0; j < sizsk; j++) {
                    sk_str = sk_str + std::to_string(s_datcursk[j]) + delim[j];
                }

                skinstr[0] = sk_str + "@";

                //end local skin


                //PSThreshold
                int sizps = 2;
                int s_datps[sizps + 1];
                s_datps[1] =  psthresholds[1] =  static_cast<int>(params.locallab.psthreshold.getTopLeft());

                s_datps[0] =  psthresholds[0] = static_cast<int>(params.locallab.psthreshold.getBottomLeft());

                std::string ps_str = "";

                ps_str = ps_str  + std::to_string(s_datps[0]) +  delim[0] + std::to_string(s_datps[1]) +  delim[1];

                pthstr[0] = ps_str + "@";
                //end local ps
				
				
                //Exp curve
                int sizex = params.locallab.excurve.size();

                if (sizex > 69) {
                    sizex = 69;//to avoid crash
                }


                int s_datcurex[sizex + 1];

                for (int j = 0; j < sizex; j++) {
                    s_datcurex[j] = excurves[0 + j] = (int)(1000. * params.locallab.excurve[j]);
                }

                std::string ex_str = "";

                for (int j = 0; j < sizex; j++) {
                    ex_str = ex_str + std::to_string(s_datcurex[j]) + delim[j];
                }

                exstr[0] = ex_str + "@";

                //end local Exp



                if (params.locallab.anbspot == 0) {
                    //update GUI and MIP after current spot ==> params, shift with the other alolistener
                    if (aloListener  && params.locallab.anbspot == 0) {
                        aloListener->localretChanged(dataspot, retistr[0], llstr[0], lhstr[0], ccstr[0], hhstr[0], skinstr[0], pthstr[0], exstr[0], 0, 1);
                    }
                }

                locallutili = false;
				localexutili = false;
                localcutili = false;
                localskutili = false;

                LHutili = false;
                HHutili = false;

            }
            //******************************************************************
            //end initialize data[xx][0] and cache with params
            //******************************************************************


            //********************************************************************
            //read mip file
            //********************************************************************
            int ns = 0;
            //    int realsp = params.locallab.nbspot;
            bool excurvret = true;
            bool excurvll = true;
            bool excurvlh = true;
            bool excurvcc = true;
            bool excurvhh = true;
            bool excurvsk = true;
            bool excpth = true;
            bool excurvex = true;

            ifstream fich(datal, ios::in);

            //read mip file
            if (fich) {//may be a file with versionmip = 10000
                //we must add new fields at the good place
                std::string line;
                std::string spotline;
                int cont = 0;
                int maxind = maxdata - 4 ; //

                if (versionmip == 10000) {
                    maxind = 49;
                    excurvret = false;
                    excurvll = false;
                    excurvlh = false;
                    excurvhh = false;
                    excurvsk = false;
                    excpth = false;
                    excurvex = false;

                }

                if (versionmip == 10001) {
                    maxind = 55;
                    excurvret = false;
                    excurvll = false;
                    excurvlh = false;
                    excurvhh = false;
                    excurvsk = false;
                    excpth = false;
                    excurvex = false;

                }

                if (versionmip == 10004) {
                    maxind = 56;
                }

                if (versionmip == 10005) {
                    excurvcc = false;
                }

                if (versionmip == 10006) {
                    maxind = 57;
                }

                if (versionmip == 10008) {
                    maxind = 58;
                }

                if (versionmip == 10009) {
                    maxind = 69;
                }

// I have forgotten 10010  ==> probably crash...but now it's passed...
//enabled this code after...
                if (versionmip == 10011) {
                    maxind = 70;
                }

                if (versionmip == 10012) {
                    maxind = 74;
                }

                if (versionmip == 10013) {
                    maxind = 77;
                }

                if (versionmip == 10015) {
                    maxind = 77;
                }

                if (versionmip == 10016) {
                    maxind = 80;
                }

                if (versionmip == 10017) {
                    maxind = 81;
                }

                if (versionmip == 10018) {
                    maxind = 82;
                }

                if (versionmip == 10019) {
                    maxind = 83;
                }

                if (versionmip == 10020) {
                    maxind = 84;
                }

                if (versionmip == 10021) {
                    maxind = 92;
                }

                if (versionmip == 10022) {
                    maxind = 93;
                }

                if (versionmip == 10023) {
                    maxind = 94;
                }

                while (getline(fich, line)) {
                    spotline = line;
                    std::size_t pos = spotline.find("=");
                    std::size_t posend = spotline.find("@");  //in case of for futur use

                    if (spotline.substr(0, pos) == "Mipversion") {
                        std::string strversion = spotline.substr(pos + 1, (posend - pos));
                        versionmip = std::stoi(strversion.c_str());
                    }

                    if (spotline.substr(0, pos) == "Spot") {
                        cont = 0;
                    }

                    cont++;
                    std::string str3 = spotline.substr(pos + 1, (posend - pos));

                    if (cont == 1) {
                        ns =  std::stoi(str3.c_str());
                    }

                    if (ns < maxspot) {
                        if (cont >= 2  && cont < 16) {
                            dataspot[cont][ns] = std::stoi(str3.c_str());
                        }

                        if (spotline.substr(0, pos) == "Currentspot") {
                            dataspot[16][0] = std::stoi(str3.c_str());
                        }

                        if (cont > 16  && cont < maxind) {
                            dataspot[cont][ns] = std::stoi(str3.c_str());
                        }


                        if (excurvret && spotline.substr(0, pos) == "curveReti") {
                            retistr[ns] = str3;
                        }

                        if (excurvll && spotline.substr(0, pos) == "curveLL") {
                            llstr[ns] = str3;
                        }


                        if (excurvlh && spotline.substr(0, pos) == "curveLH") {

                            lhstr[ns] = str3;
                        }

                        if (excurvcc && spotline.substr(0, pos) == "curveCC") {
                            ccstr[ns] = str3;
                        }

                        if (excurvhh && spotline.substr(0, pos) == "curveHH") {
                            hhstr[ns] = str3;
                        }

                        if (excurvsk && spotline.substr(0, pos) == "curveskin") {
                            skinstr[ns] = str3;
                        }

                        if (excpth && spotline.substr(0, pos) == "pthres") {
                            pthstr[ns] = str3;
                        }

                        if (excurvex && spotline.substr(0, pos) == "curveex") {
                            exstr[ns] = str3;
                        }
                    }

                }

                fich.close();
            }

            //new filed for each update
            //new fields for TM
            if (versionmip == 10000) {
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[49][sp] = 0; //stren
                    dataspot[50][sp] = 100; //gamma
                    dataspot[51][sp] = 140; //estop
                    dataspot[52][sp] = 10; //scaltm
                    dataspot[53][sp] = 0; //rewei
                    dataspot[54][sp] = 40; //sensitm

                }
            }

            if (versionmip <= 10001) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[55][sp] = 500; //retrab
                    std::string cur_str = "1000A0B120C350D350E700F500G350H350I1000J120K350L350M";//12 points
                    retistr[sp] = cur_str + "@";
                }
            }

            if (versionmip <= 10002) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    std::string ll_str = "3000A0B0C1000D1000E"; //"3000A0B0C499D501E1000F1000G"; //"3000A0B0C1000D1000E"; //"3000A0B0C200D200E800F800G1000H1000I";//"0A"
                    llstr[sp] = ll_str + "@";
                }
            }

            if (versionmip <= 10003) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    std::string lh_str = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y";
                    lhstr[sp] = lh_str + "@";
                }
            }

            if (versionmip <= 10004) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[56][sp] = 0; //curvactiv
                }
            }

            if (versionmip <= 10005) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    std::string cc_str = "3000A0B0C1000D1000E";
                    ccstr[sp] = cc_str + "@";
                }
            }

            if (versionmip <= 10006) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[57][sp] = 0; //qualitycurveMethod
                }
            }

            if (versionmip <= 10007) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    std::string hh_str = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y";
                    hhstr[sp] = hh_str + "@";
                }
            }

            if (versionmip <= 10008) {
                //vibrance
                for (int sp = 1; sp < maxspot; sp++) { // spots default

                    dataspot[58][sp] = 19;
                    dataspot[59][sp] = 0;
                    dataspot[60][sp] = 0;
                    dataspot[61][sp] = 0;
                    dataspot[62][sp] = 0;
                    dataspot[63][sp] = 0;
                    std::string sk_str = "3000A0B0C1000D1000E";
                    skinstr[sp] = sk_str + "@";
                    pthstr[sp] = "0A75B@";
                }
            }

            if (versionmip <= 10009) {//exposure
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[64][sp] = 0;
                    dataspot[65][sp] = 0;
                    dataspot[66][sp] = 20;
                    dataspot[67][sp] = 33;
                    dataspot[68][sp] = 50;
                    dataspot[69][sp] = 19;

                }
            }

            if (versionmip <= 10010) {

                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    std::string ex_str = "3000A0B0C1000D1000E";
                    exstr[sp] = ex_str + "@";
                }
            }

            if (versionmip <= 10011) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[70][sp] = 0;
                    dataspot[71][sp] = 0;
                    dataspot[72][sp] = 0;
                    dataspot[73][sp] = 0;

                }
            }

            if (versionmip <= 10012) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[74][sp] = 0;

                }
            }

            if (versionmip <= 10013) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[75][sp] = 0;
                    dataspot[76][sp] = 0;
                    dataspot[77][sp] = 1;

                }
            }

            if (versionmip <= 10015) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[78][sp] = 0;
                    dataspot[79][sp] = 19;
                    dataspot[80][sp] = 0;

                }
            }

            if (versionmip <= 10016) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[81][sp] = 0;
                }
            }

            if (versionmip <= 10017) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[82][sp] = 0;
                }
            }

            if (versionmip <= 10018) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[83][sp] = 0;
                }
            }

            if (versionmip <= 10019) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[84][sp] = 30;
                }
            }

            if (versionmip <= 10020) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[85][sp] = 0;
                    dataspot[86][sp] = 0;
                    dataspot[87][sp] = 0;
                    dataspot[88][sp] = 0;
                    dataspot[89][sp] = 0;
                    dataspot[90][sp] = 0;
                    dataspot[91][sp] = 0;
                    dataspot[92][sp] = 0;
                    dataspot[93][sp] = 0;
                }
            }

            if (versionmip <= 10021) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[94][sp] = 0;

                }
            }

            if (versionmip <= 10022) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[95][sp] = 7;

                }
            }

            if (versionmip <= 10023) {//
                for (int sp = 1; sp < maxspot; sp++) { // spots default
                    dataspot[96][sp] = 0;

                }
            }

            //**************************************************************
            //here we change the number of spot if change in options
            //**************************************************************
            if (ns < (maxspot - 1)) { //only for increasing, in case of decreasing, datas are "forgoten"
                changenumberofspot(dataspot, maxdata, maxspot, ns, datal, versionmip);

            }

            //**************************************************************************
            //end change number of spots
            //**************************************************************************

            //************************************************************************************************
            //duplicated spot
            //************************************************************************************************
            int spottodupli = dataspot[16][0];//current value

            if (params.locallab.spotduplicated) {
                lastspotdup = true;//probably unused
            }

            if (params.locallab.spotduplicated && spottodupli >= 1) {

                spotduplic(dataspot, spottodupli, maxdata);

                if (aloListener && params.locallab.spotduplicated) {
                    //update GUI and MIP
                    int sp = spottodupli;
                    int maxreal = maxdata;
                    aloListener->localChanged(dataspot, retistr[sp], llstr[sp], lhstr[sp], ccstr[sp], hhstr[sp], skinstr[sp], pthstr[sp], exstr[sp], sp, maxreal);
                    aloListener->spotdupChanged(false);//put checkbox to false, and spotduplicated to false
                    params.locallab.spotduplicated = false;

                }

            }

            //************************************************
            //end duplicated spot
            //************************************************

            //*************************************************************************
            //main algorithm for all spots
            //*************************************************************************

            for (int sp = 1; sp < maxspot; sp++) { //spots default
                params.locallab.huerefblur = dataspot[maxdata - 5][sp] / 100.;
                params.locallab.hueref = dataspot[maxdata - 4][sp] / 100.;
                params.locallab.chromaref = dataspot[maxdata - 3][sp];
                params.locallab.lumaref = dataspot[maxdata - 2][sp];
                params.locallab.sobelref = dataspot[maxdata - 1][sp];
                params.locallab.circrad = circrads[sp] = dataspot[2][sp];
                params.locallab.locX = locx[sp] = dataspot[3][sp];
                params.locallab.locY = locy[sp] = dataspot[4][sp];
                params.locallab.locYT = locyt[sp] = dataspot[5][sp];
                params.locallab.locXL = locxl[sp] = dataspot[6][sp];
                params.locallab.centerX = centerx[sp] = dataspot[7][sp];
                params.locallab.centerY = centery[sp] = dataspot[8][sp];
                params.locallab.lightness = lights[sp] = dataspot[9][sp];
                params.locallab.contrast = contrs[sp] = dataspot[10][sp];
                params.locallab.chroma = chroms[sp] = dataspot[11][sp];
                params.locallab.sensi = sensis[sp] = dataspot[12][sp];
                params.locallab.transit = transits[sp] = dataspot[13][sp];
                sps[sp] = sp;

                if (dataspot[14][sp] ==  0) {
                    inverss[sp] = 0;
                    params.locallab.invers = false;
                } else {
                    inverss[sp] = 1;
                    params.locallab.invers = true;
                }

                if (dataspot[15][sp] ==  0) {
                    smeths[sp] = 0;
                    params.locallab.Smethod = "IND" ;
                } else if (dataspot[15][sp] ==  1) {
                    smeths[sp] = 1;
                    params.locallab.Smethod = "SYM" ;
                } else if (dataspot[15][sp] ==  2) {
                    smeths[sp] = 2;
                    params.locallab.Smethod = "INDSL";
                } else if (dataspot[15][sp] ==  3) {
                    smeths[sp] = 3;
                    params.locallab.Smethod = "SYMSL";
                }

                radiuss[sp] = dataspot[17][sp];
                strengths[sp] = dataspot[18][sp];
                params.locallab.radius = dataspot[17][sp];
                params.locallab.strength = dataspot[18][sp];
                params.locallab.sensibn = sensibns[sp] = dataspot[19][sp];

                if (dataspot[20][sp] ==  0) {
                    inversrads[sp] = 0;
                    params.locallab.inversrad = false;
                } else {
                    inversrads[sp] = 1;
                    params.locallab.inversrad = true;
                }


                params.locallab.str = strs[sp] = dataspot[21][sp];
                params.locallab.chrrt = chrrts[sp] = dataspot[22][sp];
                params.locallab.neigh = neighs[sp] = dataspot[23][sp];
                params.locallab.vart = varts[sp] = dataspot[24][sp];
                params.locallab.sensih = sensihs[sp] = dataspot[25][sp];

                if (dataspot[26][sp] ==  0) {
                    inversrets[sp] = 0;
                    params.locallab.inversret = false;
                } else {
                    inversrets[sp] = 1;
                    params.locallab.inversret = true;
                }

                if (dataspot[27][sp] ==  0) {
                    retinexs[sp] = 0;
                    params.locallab.retinexMethod = "low" ;
                } else if (dataspot[27][sp] ==  1) {
                    retinexs[sp] = 1;
                    params.locallab.retinexMethod = "uni" ;
                } else if (dataspot[27][sp] ==  2) {
                    retinexs[sp] = 2;
                    params.locallab.retinexMethod = "high";
                }

                sharradiuss[sp] = dataspot[28][sp];
                params.locallab.sharradius = dataspot[28][sp];

                params.locallab.sharamount = sharamounts[sp] = dataspot[29][sp];
                params.locallab.shardamping = shardampings[sp] = dataspot[30][sp];
                params.locallab.shariter = shariters[sp] = dataspot[31][sp];
                params.locallab.sensisha = sensishas[sp] = dataspot[32][sp];

                if (dataspot[33][sp] ==  0) {
                    inversshas[sp] = 0;
                    params.locallab.inverssha = false;
                } else {
                    inversshas[sp] = 1;
                    params.locallab.inverssha = true;
                }

                if (dataspot[34][sp] ==  0) {
                    qualitys[sp] = 0;
                    params.locallab.qualityMethod = "std" ;
                } else if (dataspot[34][sp] ==  1) {
                    qualitys[sp] = 1;
                    params.locallab.qualityMethod = "enh" ;
                } else if (dataspot[34][sp] ==  2) {
                    qualitys[sp] = 2;
                    params.locallab.qualityMethod = "enhden" ;
                }

                params.locallab.thres = thress[sp] = dataspot[35][sp];
                params.locallab.proxi = proxis[sp] = dataspot[36][sp];
                params.locallab.noiselumf = noiselumfs[sp] = dataspot[37][sp];
                params.locallab.noiselumc = noiselumcs[sp] = dataspot[38][sp];
                params.locallab.noisechrof = noisechrofs[sp] = dataspot[39][sp];
                params.locallab.noisechroc = noisechrocs[sp] = dataspot[40][sp];
                params.locallab.mult[0] = mult0s[sp] = dataspot[41][sp];
                params.locallab.mult[1] = mult1s[sp] = dataspot[42][sp];
                params.locallab.mult[2] = mult2s[sp] = dataspot[43][sp];
                params.locallab.mult[3] = mult3s[sp] = dataspot[44][sp];
                params.locallab.mult[4] = mult4s[sp] = dataspot[45][sp];
                params.locallab.threshold = thresholds[sp] = dataspot[46][sp];
                params.locallab.sensicb = sensicbs[sp] = dataspot[47][sp];

                if (dataspot[48][sp] ==  0) {
                    activlums[sp] = 0;
                    params.locallab.activlum = false;
                } else {
                    activlums[sp] = 1;
                    params.locallab.activlum = true;
                }

                params.locallab.stren = strens[sp] = dataspot[49][sp];
                params.locallab.gamma = gammas[sp] = dataspot[50][sp];
                params.locallab.estop = estops[sp] = dataspot[51][sp];
                params.locallab.scaltm = scaltms[sp] = dataspot[52][sp];
                params.locallab.rewei = reweis[sp] = dataspot[53][sp];
                params.locallab.sensitm = sensitms[sp] = dataspot[54][sp];
                params.locallab.retrab = retrabs[sp] = dataspot[55][sp];

                if (dataspot[56][sp] ==  0) {
                    curvactivs[sp] = 0;
                    params.locallab.curvactiv = false;
                } else {
                    curvactivs[sp] = 1;
                    params.locallab.curvactiv = true;
                }

                if (dataspot[57][sp] ==  0) {
                    qualitycurves[sp] = 0;
                    params.locallab.qualitycurveMethod = "none" ;
                } else if (dataspot[57][sp] ==  1) {
                    qualitycurves[sp] = 1;
                    params.locallab.qualitycurveMethod = "std" ;
                } else if (dataspot[57][sp] ==  2) {
                    qualitycurves[sp] = 2;
                    params.locallab.qualitycurveMethod = "enh" ;
                }

                params.locallab.sensiv = sensivs[sp] = dataspot[58][sp];
                params.locallab.pastels = pastels[sp] =  dataspot[59][sp];
                params.locallab.saturated = saturateds[sp] = dataspot[60][sp];

                if (dataspot[61][sp] ==  0) {
                    protectskinss[sp] = 0;
                    params.locallab.protectskins = false;
                } else {
                    protectskinss[sp] = 1;
                    params.locallab.protectskins  = true;
                }

                if (dataspot[62][sp] ==  0) {
                    avoidcolorshifts[sp] = 0;
                    params.locallab.avoidcolorshift = false;
                } else {
                    avoidcolorshifts[sp] = 1;
                    params.locallab.avoidcolorshift  = true;
                }

                if (dataspot[63][sp] ==  0) {
                    pastsattogs[sp] = 0;
                    params.locallab.pastsattog = false;
                } else {
                    pastsattogs[sp] = 1;
                    params.locallab.pastsattog  = true;
                }



                params.locallab.expcomp = expcomps[sp] = dataspot[64][sp];
                params.locallab.black = blacks[sp] = dataspot[65][sp];
                params.locallab.hlcompr = hlcomprs[sp] = dataspot[66][sp];
                params.locallab.hlcomprthresh = hlcomprthreshs[sp] = dataspot[67][sp];
                params.locallab.shcompr = shcomprs[sp] = dataspot[68][sp];
                params.locallab.sensiex = sensiexs[sp] = dataspot[69][sp];

                params.locallab.centerXbuf = centerxbufs[sp] = dataspot[70][sp];
                params.locallab.centerYbuf = centerybufs[sp] = dataspot[71][sp];
                params.locallab.adjblur = adjblurs[sp] = dataspot[72][sp];

                if (dataspot[73][sp] ==  0) {
                    cutpasts[sp] = 0;
                    params.locallab.cutpast = false;
                } else {
                    cutpasts[sp] = 1;
                    params.locallab.cutpast  = true;
                }

                params.locallab.chromacbdl = chromacbdls[sp] = dataspot[74][sp];

                if (dataspot[75][sp] ==  0) {
                    lastdusts[sp] = 0;
                    params.locallab.lastdust = false;
                } else {
                    lastdusts[sp] = 1;
                    params.locallab.lastdust  = true;
                }

                if (dataspot[76][sp] ==  0) {
                    blurmets[sp] = 0;
                    params.locallab.blurMethod = "norm" ;
                } else if (dataspot[76][sp] ==  1) {
                    blurmets[sp] = 1;
                    params.locallab.blurMethod = "inv" ;
                } else if (dataspot[76][sp] ==  2) {
                    blurmets[sp] = 2;
                    params.locallab.blurMethod = "sym" ;
                }

                if (dataspot[77][sp] ==  0) {
                    dustmets[sp] = 0;
                    params.locallab.dustMethod = "cop" ;
                } else if (dataspot[77][sp] ==  1) {
                    dustmets[sp] = 1;
                    params.locallab.dustMethod = "mov" ;
                } else if (dataspot[77][sp] ==  2) {
                    dustmets[sp] = 2;
                    params.locallab.dustMethod = "pas" ;
                }

                if (dataspot[78][sp] ==  0) {
                    exclumets[sp] = 0;
                    params.locallab.Exclumethod = "norm" ;
                } else if (dataspot[78][sp] ==  1) {
                    exclumets[sp] = 1;
                    params.locallab.Exclumethod = "exc" ;
                }

                params.locallab.sensiexclu = sensiexclus[sp] = dataspot[79][sp];
                params.locallab.struc = strucs[sp] = dataspot[80][sp];
                params.locallab.warm = warms[sp] = dataspot[81][sp];
                params.locallab.noiselumdetail = noiselumdetails[sp] = dataspot[82][sp];
                params.locallab.noisechrodetail = noisechrodetails[sp] = dataspot[83][sp];
                params.locallab.sensiden = sensidens[sp] = dataspot[84][sp];

                if (dataspot[85][sp] ==  0) {
                    expdenois[sp] = 0;
                    params.locallab.expdenoi = false;
                } else {
                    expdenois[sp] = 1;
                    params.locallab.expdenoi = true;
                }

                if (dataspot[86][sp] ==  0) {
                    expcolors[sp] = 0;
                    params.locallab.expcolor = false;
                } else {
                    expcolors[sp] = 1;
                    params.locallab.expcolor = true;
                }

                if (dataspot[87][sp] ==  0) {
                    expvibrances[sp] = 0;
                    params.locallab.expvibrance = false;
                } else {
                    expvibrances[sp] = 1;
                    params.locallab.expvibrance = true;
                }

                if (dataspot[88][sp] ==  0) {
                    expblurs[sp] = 0;
                    params.locallab.expblur = false;
                } else {
                    expblurs[sp] = 1;
                    params.locallab.expblur = true;
                }

                if (dataspot[89][sp] ==  0) {
                    exptonemaps[sp] = 0;
                    params.locallab.exptonemap = false;
                } else {
                    exptonemaps[sp] = 1;
                    params.locallab.exptonemap = true;
                }

                if (dataspot[90][sp] ==  0) {
                    expretis[sp] = 0;
                    params.locallab.expreti = false;
                } else {
                    expretis[sp] = 1;
                    params.locallab.expreti = true;
                }

                if (dataspot[91][sp] ==  0) {
                    expsharps[sp] = 0;
                    params.locallab.expsharp = false;
                } else {
                    expsharps[sp] = 1;
                    params.locallab.expsharp = true;
                }

                if (dataspot[92][sp] ==  0) {
                    expcbdls[sp] = 0;
                    params.locallab.expcbdl = false;
                } else {
                    expcbdls[sp] = 1;
                    params.locallab.expcbdl = true;
                }

                if (dataspot[93][sp] ==  0) {
                    expexposes[sp] = 0;
                    params.locallab.expexpose = false;
                } else {
                    expcbdls[sp] = 1;
                    params.locallab.expexpose = true;
                }

                params.locallab.bilateral = bilaterals[sp] = dataspot[94][sp];
                params.locallab.noiselequal = noiselequals[sp] = dataspot[95][sp];

                if (dataspot[96][sp] ==  0) {
                    shapemets[sp] = 0;
                    params.locallab.shapemethod = "ELI" ;
                } else if (dataspot[96][sp] ==  1) {
                    shapemets[sp] = 1;
                    params.locallab.shapemethod = "RECT" ;
                }

                int *s_datc;
                s_datc = new int[70];
                int siz;

                ipf.strcurv_data(retistr[sp], s_datc, siz); //convert data in int string with strcurv_data () - it is a work around !

                sizeretics[sp] = siz;

                std::vector<double>   cretiend;

                for (int j = 0; j < siz; j++) {
                    reticurvs[sp * 500 + j] =  s_datc[j];
                    cretiend.push_back((double)(s_datc[j]) / 1000.);
                }

                delete [] s_datc;

                int *s_datcl;
                s_datcl = new int[70];
                int sizl;

                ipf.strcurv_data(llstr[sp], s_datcl, sizl);

                sizellcs[sp] = sizl;

                std::vector<double>   cllend;

                for (int j = 0; j < sizl; j++) {
                    llcurvs[sp * 500 + j] =  s_datcl[j];
                    cllend.push_back((double)(s_datcl[j]) / 1000.);
                }

                delete [] s_datcl;


                int *s_datcc;
                s_datcc = new int[70];
                int sizc;

                ipf.strcurv_data(ccstr[sp], s_datcc, sizc);

                sizecccs[sp] = sizc;

                std::vector<double>   cccend;

                for (int j = 0; j < sizc; j++) {
                    cccurvs[sp * 500 + j] =  s_datcc[j];
                    cccend.push_back((double)(s_datcc[j]) / 1000.);
                }

                delete [] s_datcc;

                int *s_datch;
                s_datch = new int[70];
                int sizh;

                ipf.strcurv_data(lhstr[sp], s_datch, sizh);

                sizelhcs[sp] = sizh;

                std::vector<double>   clhend;

                for (int j = 0; j < sizh; j++) {
                    lhcurvs[sp * 500 + j] =  s_datch[j];
                    clhend.push_back((double)(s_datch[j]) / 1000.);
                }

                delete [] s_datch;


                int *s_datchh;
                s_datchh = new int[70];
                int sizhh;

                ipf.strcurv_data(hhstr[sp], s_datchh, sizhh);

                sizehhcs[sp] = sizhh;

                std::vector<double>   chhend;

                for (int j = 0; j < sizhh; j++) {
                    hhcurvs[sp * 500 + j] =  s_datchh[j];
                    chhend.push_back((double)(s_datchh[j]) / 1000.);
                }

                delete [] s_datchh;


                int *s_datcsk;
                s_datcsk = new int[70];
                int sizsk;

                ipf.strcurv_data(skinstr[sp], s_datcsk, sizsk);

                sizeskintonecurves[sp] = sizsk;

                std::vector<double>   cskend;

                for (int j = 0; j < sizsk; j++) {
                    skintonescurves[sp * 500 + j] =  s_datcsk[j];
                    cskend.push_back((double)(s_datcsk[j]) / 1000.);
                }

                delete [] s_datcsk;

                //PSThreshold + 1
                int sizps = 2;
                int s_datcps[sizps + 1];
                ipf.strcurv_data(pthstr[sp], s_datcps, sizps);

                psthresholds[sp * 500] = s_datcps[0];
                psthresholds[sp * 500 + 1] = s_datcps[1];
                //  printf("A 0=%i 1=%i\n", s_datcps[0], s_datcps[1]);
                params.locallab.psthreshold.setValues(s_datcps[0], s_datcps[1]);

                //end local PS

                //exposure
                int *s_datcexx;
                s_datcexx = new int[70];
                int sizexx;

                ipf.strcurv_data(exstr[sp], s_datcexx, sizexx);

                sizeexcurves[sp] = sizexx;

                std::vector<double>   cexend;

                for (int j = 0; j < sizexx; j++) {
                    excurves[sp * 500 + j] =  s_datcexx[j];
                    cexend.push_back((double)(s_datcexx[j]) / 1000.);
                }

                delete [] s_datcexx;


                params.locallab.localTgaincurve.clear();
                params.locallab.localTgaincurve = cretiend;

                //       int lenc = params.locallab.localTgaincurve.size();

                params.locallab.llcurve.clear();
                params.locallab.llcurve = cllend;

                params.locallab.LHcurve.clear();
                params.locallab.LHcurve = clhend;

                params.locallab.cccurve.clear();
                params.locallab.cccurve = cccend;

                params.locallab.HHcurve.clear();
                params.locallab.HHcurve = chhend;

                params.locallab.skintonescurve.clear();
                params.locallab.skintonescurve = cskend;


                params.locallab.excurve.clear();
                params.locallab.excurve = cexend;

                locallutili = false;
                localcutili = false;
                localskutili = false;
                localexutili = false;

                LHutili = false;
                HHutili = false;
                std::string t_curvhhref = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y@";

                if (lhstr[sp].c_str() != t_curvhhref) {
                    //LHutili = true;
                }

                if (hhstr[sp].c_str() != t_curvhhref) {
                  //  HHutili = true;
                }

                std::string t_curvskinref = "3000A0B0C1000D1000E@";
                std::string t_none = "0A@";

                if (skinstr[sp].c_str() != t_curvskinref  && skinstr[sp].c_str() != t_none) {
                 //   localskutili = true;
                }

                std::string t_curvexref = "3000A0B0C1000D1000E@";

                if (exstr[sp].c_str() != t_curvexref  && exstr[sp].c_str() != t_none) {
                  //  localexutili = true;
                }

                params.locallab.getCurves(locRETgainCurve, locRETgainCurverab, loclhCurve, lochhCurve, LHutili, HHutili);
                CurveFactory::curveLocal(locallutili, params.locallab.llcurve, lllocalcurve, sca);
                CurveFactory::curveCCLocal(localcutili, params.locallab.cccurve, cclocalcurve, sca);
                CurveFactory::curveskLocal(localskutili, params.locallab.skintonescurve, sklocalcurve, sca);
                CurveFactory::curveexLocal(localexutili, params.locallab.excurve, exlocalcurve, sca);
                //provisory
                double br = 0.;
                double contr = 0.;
                double ecomp = params.locallab.expcomp;
                double black = params.locallab.black;
                double hlcompr = params.locallab.hlcompr;
                double hlcomprthresh = params.locallab.hlcomprthresh;
                double shcompr = params.locallab.shcompr;

                CurveFactory::complexCurvelocal(ecomp, black / 65535., hlcompr, hlcomprthresh, shcompr, br, contr,
                                                lhist16, hltonecurveloc, shtonecurveloc, tonecurveloc,
                                                sca);

                double huere, chromare, lumare, huerefblu;
                double sobelre;

                ipf.calc_ref(nprevl, nprevl, 0, 0, pW, pH, scale, huerefblu, huere, chromare, lumare, sobelre);
                huerblu = huerefblu;
                huer = huere;
                chromar = chromare;
                lumar = lumare ;
                sobeler = sobelre;
                params.locallab.huerefblur = huerblu;
                params.locallab.hueref = huer;
                params.locallab.chromaref = chromar;
                params.locallab.lumaref = lumar;
                params.locallab.sobelref = sobeler;

                dataspot[maxdata - 5][sp] = huerefblurs[sp] = 100.f * params.locallab.huerefblur;
                dataspot[maxdata - 4][sp] = huerefs[sp] = 100.f * params.locallab.hueref;
                dataspot[maxdata - 3][sp] = chromarefs[sp] = params.locallab.chromaref;
                dataspot[maxdata - 2][sp] = lumarefs[sp] = params.locallab.lumaref;
                dataspot[maxdata - 1][sp] = sobelrefs[sp] = params.locallab.sobelref;
                //printf("sp=%i huerefsp=%f\n", sp, huerefs[sp]);
                ipf.Lab_Local(3, maxspot, sp, huerefs, sobelrefs, centerx, centery, (float**)shbuffer, nprevl, nprevl, reserv, 0, 0, pW, pH, scale, locRETgainCurve, lllocalcurve, loclhCurve,  lochhCurve,
                              LHutili, HHutili, cclocalcurve, localskutili, sklocalcurve, localexutili, exlocalcurve, hltonecurveloc, shtonecurveloc, tonecurveloc, params.locallab.huerefblur, params.locallab.hueref, params.locallab.chromaref, params.locallab.lumaref, params.locallab.sobelref);
                lllocalcurve.clear();
                cclocalcurve.clear();
                sklocalcurve.clear();
                exlocalcurve.clear();

            }

            //***********************************************************
            //end main algoritm
            //***********************************************************

            int sp ;
            sp = realspot;
            //now for current spot
            int maxreal = maxdata;

            //*************************************************************
            //update GUI and Mip files
            //*************************************************************
            if (aloListener && realspot != dataspot[16][0]) {
                //update GUI and MIP
                aloListener->localChanged(dataspot, retistr[sp], llstr[sp], lhstr[sp], ccstr[sp], hhstr[sp], skinstr[sp], pthstr[sp], exstr[sp], sp, maxreal);
            }


            //****************************************************************
            //now works on current spot
            //****************************************************************
            params.locallab.huerefblur = INFINITY;
            params.locallab.hueref = INFINITY;
            params.locallab.chromaref = INFINITY;
            params.locallab.lumaref = INFINITY;
            params.locallab.sobelref = INFINITY;
            locallutili = false;
            localexutili = false;
              //  locallutili = false;
                localcutili = false;
                localskutili = false;
             //   localexutili = false;

                LHutili = false;
                HHutili = false;

            sps[sp] = sp;
            dataspot[2][sp] = circrads[sp] = params.locallab.circrad = dataspot[2][0];
            dataspot[3][sp] = locx[sp] = params.locallab.locX = dataspot[3][0];
            dataspot[4][sp] = locy[sp] = params.locallab.locY = dataspot[4][0];
            dataspot[5][sp] = locyt[sp] = params.locallab.locYT = dataspot[5][0];
            dataspot[6][sp] = locxl[sp] = params.locallab.locXL = dataspot[6][0];
            dataspot[7][sp] = centerx[sp] = params.locallab.centerX = dataspot[7][0];
            dataspot[8][sp] = centery[sp] = params.locallab.centerY = dataspot[8][0];
            dataspot[9][sp] = lights[sp] = params.locallab.lightness = dataspot[9][0];
            dataspot[10][sp] = contrs[sp] = params.locallab.contrast = dataspot[10][0];
            dataspot[11][sp] = chroms[sp] = params.locallab.chroma = dataspot[11][0];
            dataspot[12][sp] = sensis[sp] = params.locallab.sensi = dataspot[12][0];
            dataspot[13][sp] = transits[sp] = params.locallab.transit = dataspot[13][0];

            if (dataspot[14][0] == 0) {
                params.locallab.invers = false;
                dataspot[14][sp] = 0;
                inverss[sp] = 0;
            } else {
                params.locallab.invers = true;
                dataspot[14][sp] = 1;
                inverss[sp] = 1;
            }

            if (dataspot[15][0] == 0) {
                params.locallab.Smethod = "IND" ;
                smeths[sp] = 0;
                dataspot[15][sp] = 0;
            } else if (dataspot[15][0] == 1) {
                params.locallab.Smethod = "SYM" ;
                smeths[sp] = 1;
                dataspot[15][sp] = 1;

            } else if (dataspot[15][0] == 2) {
                params.locallab.Smethod = "INDSL" ;
                smeths[sp] = 2;
                dataspot[15][sp] = 2;
            } else if (dataspot[15][0] == 3) {
                params.locallab.Smethod = "SYMSL" ;
                smeths[sp] = 3;
                dataspot[15][sp] = 3;
            }

            params.locallab.radius = dataspot[17][0];
            params.locallab.strength = dataspot[18][0];
            params.locallab.sensibn = dataspot[19][0];

            dataspot[17][sp] = radiuss[sp] = params.locallab.radius;
            dataspot[18][sp] = strengths[sp] = params.locallab.strength;
            dataspot[19][sp] = sensibns[sp] = params.locallab.sensibn;

            if (dataspot[20][0] == 0) {
                params.locallab.inversrad = false;
                dataspot[20][sp] = 0;
                inversrads[sp] = 0;
            } else {
                params.locallab.inversrad = true;
                dataspot[20][sp] = 1;
                inversrads[sp] = 1;
            }

            dataspot[21][sp] = strs[sp] = params.locallab.str = dataspot[21][0];
            dataspot[22][sp] = chrrts[sp] = params.locallab.chrrt = dataspot[22][0];
            dataspot[23][sp] = neighs[sp] = params.locallab.neigh = dataspot[23][0];
            dataspot[24][sp] = varts[sp] = params.locallab.vart = dataspot[24][0];
            dataspot[25][sp] = sensihs[sp] = params.locallab.sensih = dataspot[25][0];

            if (dataspot[26][0] == 0) {
                params.locallab.inversret = false;
                inversrets[sp] = 0;
                dataspot[26][sp] = 0;
            } else {
                params.locallab.inversret = true;
                inversrets[sp] = 1;
                dataspot[26][sp] = 1;
            }

            if (dataspot[27][0] == 0) {
                params.locallab.retinexMethod = "low" ;
                retinexs[sp] = 0;
                dataspot[27][sp] = 0;
            } else if (dataspot[27][0] == 1) {
                params.locallab.retinexMethod = "uni" ;
                retinexs[sp] = 1;
                dataspot[27][sp] = 1;
            } else if (dataspot[27][0] == 2) {
                params.locallab.Smethod = "high" ;
                retinexs[sp] = 2;
                dataspot[27][sp] = 2;
            }

            dataspot[28][sp] = sharradiuss[sp] = params.locallab.sharradius = dataspot[28][0];

            dataspot[29][sp] = sharamounts[sp] = params.locallab.sharamount = dataspot[29][0];
            dataspot[30][sp] = shardampings[sp] = params.locallab.shardamping = dataspot[30][0];
            dataspot[31][sp] = shariters[sp] = params.locallab.shariter = dataspot[31][0];
            dataspot[32][sp] = sensishas[sp] = params.locallab.sensisha = dataspot[32][0];

            if (dataspot[33][0] == 0) {
                params.locallab.inverssha = 0;
                inversshas[sp] = 0;
                dataspot[33][sp] = 0;
            } else {
                params.locallab.inverssha = 1;
                inversshas[sp] = 1;
                dataspot[33][sp] = 1;
            }

            if (dataspot[34][0] == 0) {
                params.locallab.qualityMethod = "std" ;
                qualitys[sp] = 0;
                dataspot[34][sp] = 0;
            } else if (dataspot[34][0] == 1) {
                params.locallab.qualityMethod = "enh" ;
                qualitys[sp] = 1;
                dataspot[34][sp] = 1;
            } else if (dataspot[34][0] == 2) {
                params.locallab.qualityMethod = "enhden" ;
                qualitys[sp] = 2;
                dataspot[34][sp] = 2;
            }

            dataspot[35][sp] = thress[sp] = params.locallab.thres = dataspot[35][0];
            dataspot[36][sp] = proxis[sp] = params.locallab.proxi = dataspot[36][0];
            dataspot[37][sp] = noiselumfs[sp] = params.locallab.noiselumf = dataspot[37][0];
            dataspot[38][sp] = noiselumcs[sp] = params.locallab.noiselumc = dataspot[38][0];
            dataspot[39][sp] = noisechrofs[sp] = params.locallab.noisechrof = dataspot[39][0];
            dataspot[40][sp] = noisechrocs[sp] = params.locallab.noisechroc = dataspot[40][0];
            dataspot[41][sp] = mult0s[sp] = params.locallab.mult[0] = dataspot[41][0];
            dataspot[42][sp] = mult1s[sp] = params.locallab.mult[1] = dataspot[42][0];
            dataspot[43][sp] = mult2s[sp] = params.locallab.mult[2] = dataspot[43][0];
            dataspot[44][sp] = mult3s[sp] = params.locallab.mult[3] = dataspot[44][0];
            dataspot[45][sp] = mult4s[sp] = params.locallab.mult[4] = dataspot[45][0];
            dataspot[46][sp] = thresholds[sp] = params.locallab.threshold = dataspot[46][0];
            dataspot[47][sp] = sensicbs[sp] = params.locallab.sensicb = dataspot[47][0];

            if (dataspot[48][0] == 0) {
                params.locallab.activlum = 0;
                activlums[sp] = 0;
                dataspot[48][sp] = 0;
            } else {
                params.locallab.activlum = 1;
                activlums[sp] = 1;
                dataspot[48][sp] = 1;
            }

            dataspot[49][sp] = strens[sp] = params.locallab.stren = dataspot[49][0];
            dataspot[50][sp] = gammas[sp] = params.locallab.gamma = dataspot[50][0];
            dataspot[51][sp] = estops[sp] = params.locallab.estop = dataspot[51][0];
            dataspot[52][sp] = scaltms[sp] = params.locallab.scaltm = dataspot[52][0];
            dataspot[53][sp] = reweis[sp] = params.locallab.rewei = dataspot[53][0];
            dataspot[54][sp] = sensitms[sp] = params.locallab.sensitm = dataspot[54][0];
            dataspot[55][sp] = retrabs[sp] = params.locallab.retrab = dataspot[55][0];

            if (dataspot[56][0] == 0) {
                params.locallab.curvactiv = false;
                dataspot[56][sp] = 0;
                curvactivs[sp] = 0;
            } else {
                params.locallab.curvactiv = true;
                dataspot[56][sp] = 1;
                curvactivs[sp] = 1;
            }

            if (dataspot[57][0] == 0) {
                params.locallab.qualitycurveMethod = "none" ;
                qualitycurves[sp] = 0;
                dataspot[57][sp] = 0;
            } else if (dataspot[57][0] == 1) {
                params.locallab.qualitycurveMethod = "std" ;
                qualitycurves[sp] = 1;
                dataspot[57][sp] = 1;
            } else if (dataspot[57][0] == 2) {
                params.locallab.qualitycurveMethod = "enh" ;
                qualitycurves[sp] = 2;
                dataspot[57][sp] = 2;
            }

            dataspot[58][sp] = sensivs[sp] = params.locallab.sensiv = dataspot[58][0];
            dataspot[59][sp] = pastels[sp] = params.locallab.pastels = dataspot[59][0];
            dataspot[60][sp] = saturateds[sp] = params.locallab.saturated = dataspot[60][0];

            if (dataspot[61][0] == 0) {
                params.locallab.protectskins = false;
                dataspot[61][sp] = 0;
                protectskinss[sp] = 0;
            } else {
                params.locallab.protectskins = true;
                dataspot[61][sp] = 1;
                protectskinss[sp] = 1;
            }

            if (dataspot[62][0] == 0) {
                params.locallab.avoidcolorshift = false;
                dataspot[62][sp] = 0;
                avoidcolorshifts[sp] = 0;
            } else {
                params.locallab.avoidcolorshift = true;
                dataspot[62][sp] = 1;
                avoidcolorshifts[sp] = 1;
            }

            if (dataspot[63][0] == 0) {
                params.locallab.pastsattog = false;
                dataspot[63][sp] = 0;
                pastsattogs[sp] = 0;
            } else {
                params.locallab.pastsattog = true;
                dataspot[63][sp] = 1;
                pastsattogs[sp] = 1;
            }

            dataspot[64][sp] = expcomps[sp] = params.locallab.expcomp = dataspot[64][0];
            dataspot[65][sp] = blacks[sp] = params.locallab.black = dataspot[65][0];
            dataspot[66][sp] = hlcomprs[sp] = params.locallab.hlcompr = dataspot[66][0];
            dataspot[67][sp] = hlcomprthreshs[sp] = params.locallab.hlcomprthresh = dataspot[67][0];
            dataspot[68][sp] = shcomprs[sp] = params.locallab.shcompr = dataspot[68][0];
            dataspot[69][sp] = sensiexs[sp] = params.locallab.sensiex = dataspot[69][0];

            dataspot[70][sp] = centerxbufs[sp] = params.locallab.centerXbuf = dataspot[70][0];
            dataspot[71][sp] = centerybufs[sp] = params.locallab.centerYbuf = dataspot[71][0];
            dataspot[72][sp] = adjblurs[sp] = params.locallab.adjblur = dataspot[72][0];

            if (dataspot[73][0] == 0) {
                params.locallab.cutpast = false;
                dataspot[73][sp] = 0;
                cutpasts[sp] = 0;
            } else {
                params.locallab.cutpast = true;
                dataspot[73][sp] = 1;
                cutpasts[sp] = 1;
            }

            dataspot[74][sp] = chromacbdls[sp] = params.locallab.chromacbdl = dataspot[74][0];

            if (dataspot[75][0] == 0) {
                params.locallab.lastdust = false;
                dataspot[75][sp] = 0;
                lastdusts[sp] = 0;
            } else {
                params.locallab.lastdust = true;
                dataspot[75][sp] = 1;
                lastdusts[sp] = 1;
            }

            if (dataspot[76][0] == 0) {
                params.locallab.blurMethod = "norm" ;
                blurmets[sp] = 0;
                dataspot[76][sp] = 0;
            } else if (dataspot[76][0] == 1) {
                params.locallab.blurMethod = "inv" ;
                blurmets[sp] = 1;
                dataspot[76][sp] = 1;
            } else if (dataspot[76][0] == 2) {
                params.locallab.blurMethod = "sym" ;
                blurmets[sp] = 2;
                dataspot[76][sp] = 2;
            }

            if (dataspot[77][0] == 0) {
                params.locallab.dustMethod = "cop" ;
                dustmets[sp] = 0;
                dataspot[77][sp] = 0;
            } else if (dataspot[77][0] == 1) {
                params.locallab.dustMethod = "mov" ;
                dustmets[sp] = 1;
                dataspot[77][sp] = 1;
            } else if (dataspot[77][0] == 2) {
                params.locallab.dustMethod = "pas" ;
                dustmets[sp] = 2;
                dataspot[77][sp] = 2;
            }

            if (dataspot[78][0] == 0) {
                params.locallab.Exclumethod = "norm" ;
                exclumets[sp] = 0;
                dataspot[78][sp] = 0;
            } else if (dataspot[78][0] == 1) {
                params.locallab.Exclumethod = "exc" ;
                exclumets[sp] = 1;
                dataspot[78][sp] = 1;
            }

            dataspot[79][sp] = sensiexclus[sp] = params.locallab.sensiexclu = dataspot[79][0];
            dataspot[80][sp] = strucs[sp] = params.locallab.struc = dataspot[80][0];
            dataspot[81][sp] = warms[sp] = params.locallab.warm = dataspot[81][0];
            dataspot[82][sp] = noiselumdetails[sp] = params.locallab.noiselumdetail = dataspot[82][0];
            dataspot[83][sp] = noisechrodetails[sp] = params.locallab.noisechrodetail = dataspot[83][0];
            dataspot[84][sp] = sensidens[sp] = params.locallab.sensiden = dataspot[84][0];

            if (dataspot[85][0] == 0) {
                params.locallab.expdenoi = false;
                dataspot[85][sp] = 0;
                expdenois[sp] = 0;
            } else {
                params.locallab.expdenoi = true;
                dataspot[85][sp] = 1;
                expdenois[sp] = 1;
            }

            if (dataspot[86][0] == 0) {
                params.locallab.expcolor = false;
                dataspot[86][sp] = 0;
                expcolors[sp] = 0;
            } else {
                params.locallab.expcolor = true;
                dataspot[86][sp] = 1;
                expcolors[sp] = 1;
            }

            if (dataspot[87][0] == 0) {
                params.locallab.expvibrance = false;
                dataspot[87][sp] = 0;
                expvibrances[sp] = 0;
            } else {
                params.locallab.expvibrance = true;
                dataspot[87][sp] = 1;
                expvibrances[sp] = 1;
            }

            if (dataspot[88][0] == 0) {
                params.locallab.expblur = false;
                dataspot[88][sp] = 0;
                expblurs[sp] = 0;
            } else {
                params.locallab.expblur = true;
                dataspot[88][sp] = 1;
                expblurs[sp] = 1;
            }

            if (dataspot[89][0] == 0) {
                params.locallab.exptonemap = false;
                dataspot[89][sp] = 0;
                exptonemaps[sp] = 0;
            } else {
                params.locallab.exptonemap = true;
                dataspot[89][sp] = 1;
                exptonemaps[sp] = 1;
            }

            if (dataspot[90][0] == 0) {
                params.locallab.expreti = false;
                dataspot[90][sp] = 0;
                expretis[sp] = 0;
            } else {
                params.locallab.expreti = true;
                dataspot[90][sp] = 1;
                expretis[sp] = 1;
            }

            if (dataspot[91][0] == 0) {
                params.locallab.expsharp = false;
                dataspot[91][sp] = 0;
                expsharps[sp] = 0;
            } else {
                params.locallab.expsharp = true;
                dataspot[91][sp] = 1;
                expsharps[sp] = 1;
            }

            if (dataspot[92][0] == 0) {
                params.locallab.expcbdl = false;
                dataspot[92][sp] = 0;
                expcbdls[sp] = 0;
            } else {
                params.locallab.expcbdl = true;
                dataspot[92][sp] = 1;
                expcbdls[sp] = 1;
            }

            if (dataspot[93][0] == 0) {
                params.locallab.expexpose = false;
                dataspot[93][sp] = 0;
                expexposes[sp] = 0;
            } else {
                params.locallab.expexpose = true;
                dataspot[93][sp] = 1;
                expexposes[sp] = 1;
            }

            dataspot[94][sp] = bilaterals[sp] = params.locallab.bilateral = dataspot[94][0];
            dataspot[95][sp] = noiselequals[sp] = params.locallab.noiselequal = dataspot[95][0];

            if (dataspot[96][0] == 0) {
                params.locallab.shapemethod = "ELI" ;
                shapemets[sp] = 0;
                dataspot[96][sp] = 0;
            } else if (dataspot[96][0] == 1) {
                params.locallab.shapemethod = "RECT" ;
                shapemets[sp] = 1;
                dataspot[96][sp] = 1;
            }

            int *s_datc;
            s_datc = new int[70];
            int siz;

            ipf.strcurv_data(retistr[0], s_datc, siz);
            sizeretics[sp] = siz;
            std::vector<double>   cretiend;

            retistr[sp] = retistr[0];

            for (int j = 0; j < siz; j++) {
                reticurvs[sp * 500 + j] = s_datc[j];
                cretiend.push_back((double)(s_datc[j]) / 1000.);

            }

            params.locallab.localTgaincurve.clear();
            params.locallab.localTgaincurve = cretiend;

            delete [] s_datc;

            int *s_datcl;
            s_datcl = new int[70];
            int sizl;

            ipf.strcurv_data(llstr[0], s_datcl, sizl);
            sizellcs[sp] = sizl;
            std::vector<double>   cllend;

            llstr[sp] = llstr[0];

            for (int j = 0; j < sizl; j++) {
                llcurvs[sp * 500 + j] = s_datcl[j];
                cllend.push_back((double)(s_datcl[j]) / 1000.);

            }

            params.locallab.llcurve.clear();
            params.locallab.llcurve = cllend;

            delete [] s_datcl;


            int *s_datcc;
            s_datcc = new int[70];
            int sizc;

            ipf.strcurv_data(ccstr[0], s_datcc, sizc);
            sizecccs[sp] = sizc;
            std::vector<double>   cccend;

            ccstr[sp] = ccstr[0];

            for (int j = 0; j < sizc; j++) {
                cccurvs[sp * 500 + j] = s_datcc[j];
                cccend.push_back((double)(s_datcc[j]) / 1000.);

            }

            params.locallab.cccurve.clear();
            params.locallab.cccurve = cccend;

            delete [] s_datcc;

            int *s_datch;
            s_datch = new int[70];
            int sizh;

            ipf.strcurv_data(lhstr[0], s_datch, sizh);
            sizelhcs[sp] = sizh;
            std::vector<double>   clhend;

            lhstr[sp] = lhstr[0];

            for (int j = 0; j < sizh; j++) {
                lhcurvs[sp * 500 + j] = s_datch[j];
                clhend.push_back((double)(s_datch[j]) / 1000.);

            }

            params.locallab.LHcurve.clear();
            params.locallab.LHcurve = clhend;

            delete [] s_datch;

            int *s_datchh;
            s_datchh = new int[70];
            int sizhh;

            ipf.strcurv_data(hhstr[0], s_datchh, sizhh);
            sizehhcs[sp] = sizhh;
            std::vector<double>   chhend;

            hhstr[sp] = hhstr[0];

            for (int j = 0; j < sizhh; j++) {
                hhcurvs[sp * 500 + j] = s_datchh[j];
                chhend.push_back((double)(s_datchh[j]) / 1000.);

            }

            params.locallab.HHcurve.clear();
            params.locallab.HHcurve = chhend;

            delete [] s_datchh;

            int *s_datcsk;
            s_datcsk = new int[70];
            int sizsk;

            ipf.strcurv_data(skinstr[0], s_datcsk, sizsk);
            sizeskintonecurves[sp] = sizsk;
            std::vector<double>   cskend;

            skinstr[sp] = skinstr[0];

            for (int j = 0; j < sizsk; j++) {
                skintonescurves[sp * 500 + j] = s_datcsk[j];
                cskend.push_back((double)(s_datcsk[j]) / 1000.);

            }

            params.locallab.skintonescurve.clear();
            params.locallab.skintonescurve = cskend;

            delete [] s_datcsk;



            //PSThreshold + 1
            int sizps = 2;
            int s_datcps[sizps + 1];
            ipf.strcurv_data(pthstr[0], s_datcps, sizps);

            psthresholds[sp * 500] = s_datcps[0];
            psthresholds[sp * 500 + 1] = s_datcps[1];
//          printf("B 0=%i 1=%i\n", s_datcps[0], s_datcps[1]);
            std::string ps_str2 = "";

            ps_str2 = ps_str2  + std::to_string(s_datcps[0]) +  delim[0] + std::to_string(s_datcps[1]) +  delim[1];
            pthstr[0] = ps_str2 + "@";

            pthstr[sp] = pthstr[0];
            params.locallab.psthreshold.setValues(s_datcps[0], s_datcps[1]);

            //end local PS

			
            //expos
            int *s_datcex;
            s_datcex = new int[70];
            int sizex;
            //   printf ("ex0=%s \n", exstr[0].c_str());
            ipf.strcurv_data(exstr[0], s_datcex, sizex);
            sizeexcurves[sp] = sizex;
            std::vector<double>   cexend;

            exstr[sp] = exstr[0];

            for (int j = 0; j < sizex; j++) {
                excurves[sp * 500 + j] = s_datcex[j];
                cexend.push_back((double)(s_datcex[j]) / 1000.);

            }

            params.locallab.excurve.clear();
            params.locallab.excurve = cexend;

            delete [] s_datcex;

            LHutili = false;
            HHutili = false;

            std::string t_curvhhref2 = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y@";

            if (hhstr[sp].c_str() != t_curvhhref2) {
              //  HHutili = true;
            }

            if (lhstr[sp].c_str() != t_curvhhref2) {
             //   LHutili = true;
            }

            params.locallab.getCurves(locRETgainCurve, locRETgainCurverab, loclhCurve, lochhCurve, LHutili, HHutili);
            locallutili = false;
            localcutili = false;
            localskutili = false;
            localexutili = false;

            std::string t_curvskinref2 = "3000A0B0C1000D1000E@";
            std::string t_none2 = "0A@";

            if (skinstr[sp].c_str() != t_curvskinref2 && skinstr[sp].c_str() != t_none2) {
              // localskutili = true;
            }

            std::string t_curvexref2 = "3000A0B0C1000D1000E@";

            if (exstr[sp].c_str() != t_curvexref2 && exstr[sp].c_str() != t_none2) {
              //  localexutili = true;
            }

            double br = 0.;
            double contr = 0.;
            double ecomp = params.locallab.expcomp;
            double black = params.locallab.black;
            double hlcompr = params.locallab.hlcompr;
            double hlcomprthresh = params.locallab.hlcomprthresh;
            double shcompr = params.locallab.shcompr;

            CurveFactory::complexCurvelocal(ecomp, black / 65535., hlcompr, hlcomprthresh, shcompr, br, contr,
                                            lhist16, hltonecurveloc, shtonecurveloc, tonecurveloc,
                                            sca);

            CurveFactory::curveLocal(locallutili, params.locallab.llcurve, lllocalcurve, sca);
            CurveFactory::curveCCLocal(localcutili, params.locallab.cccurve, cclocalcurve, sca);
            CurveFactory::curveskLocal(localskutili, params.locallab.skintonescurve, sklocalcurve, sca);
            CurveFactory::curveexLocal(localexutili, params.locallab.excurve, exlocalcurve, sca);

            params.locallab.huerefblur = huerefblurs[sp] / 100.;
            params.locallab.hueref = huerefs[sp] / 100.;
            params.locallab.chromaref = chromarefs[sp];
            params.locallab.lumaref = lumarefs[sp];
            params.locallab.sobelref = sobelrefs[sp];
            ipf.Lab_Local(3, maxspot, sp, huerefs, sobelrefs, centerx, centery, (float**)shbuffer, nprevl, nprevl, reserv, 0, 0, pW, pH, scale, locRETgainCurve, lllocalcurve, loclhCurve, lochhCurve, LHutili, HHutili, cclocalcurve,
                          localskutili, sklocalcurve, localexutili, exlocalcurve, hltonecurveloc, shtonecurveloc, tonecurveloc, params.locallab.huerefblur, params.locallab.hueref, params.locallab.chromaref, params.locallab.lumaref, params.locallab.sobelref);
            lllocalcurve.clear();
            cclocalcurve.clear();
            sklocalcurve.clear();
            exlocalcurve.clear();

            //*******************************************************
            //end current spot
            //*******************************************************

            //*******************************************************
            //write mip file in real time
            //*******************************************************

            ofstream fou(datal, ios::out | ios::trunc);

            if (fou) {
                mipver = 10024; //confirm new versionmip - last - same value as at beginning,

                for (int spe = 1; spe < maxspot; spe++) {
                    int t_sp = spe;
                    int t_mipversion = mipver;//actualize for current usage
                    int t_circrad  = dataspot[2][spe];
                    int t_locX  = dataspot[3][spe];
                    int t_locY  = dataspot[4][spe];
                    int t_locYT  = dataspot[5][spe];
                    int t_locXL  = dataspot[6][spe];
                    int t_centerX  = dataspot[7][spe];
                    int t_centerY  = dataspot[8][spe];
                    int t_lightness  = dataspot[9][spe];
                    int t_contrast  = dataspot[10][spe];
                    int t_chroma  = dataspot[11][spe];
                    int t_sensi  = dataspot[12][spe];
                    int t_transit  = dataspot[13][spe];
                    int t_invers = dataspot[14][spe];
                    int t_Smeth = dataspot[15][spe];
                    int t_currentspot  = realspot;
                    int t_radius = dataspot[17][spe];
                    int t_strength = dataspot[18][spe];
                    int t_sensibn = dataspot[19][spe];
                    int t_inversrad = dataspot[20][spe];
                    int t_str = dataspot[21][spe];
                    int t_chrrt = dataspot[22][spe];
                    int t_neigh = dataspot[23][spe];
                    int t_vart = dataspot[24][spe];
                    int t_sensih = dataspot[25][spe];
                    int t_inversret = dataspot[26][spe];
                    int t_retinexMethod = dataspot[27][spe];
                    int t_sharradius = dataspot[28][spe];
                    int t_sharamount = dataspot[29][spe];
                    int t_shardamping = dataspot[30][spe];
                    int t_shariter = dataspot[31][spe];
                    int t_sensisha = dataspot[32][spe];
                    int t_inverssha = dataspot[33][spe];
                    int t_qualityMethod =  dataspot[34][spe];
                    int t_thres =  dataspot[35][spe];
                    int t_proxi =  dataspot[36][spe];
                    int t_noiselumf = dataspot[37][spe];
                    int t_noiselumc = dataspot[38][spe];
                    int t_noisechrof = dataspot[39][spe];
                    int t_noisechroc = dataspot[40][spe];
                    int t_mult0 = dataspot[41][spe];
                    int t_mult1 = dataspot[42][spe];
                    int t_mult2 = dataspot[43][spe];
                    int t_mult3 = dataspot[44][spe];
                    int t_mult4 = dataspot[45][spe];
                    int t_threshold = dataspot[46][spe];
                    int t_sensicb = dataspot[47][spe];
                    int t_activlum = dataspot[48][spe];

                    int t_stren = dataspot[49][spe];
                    int t_gamma = dataspot[50][spe];
                    int t_estop = dataspot[51][spe];
                    int t_scaltm = dataspot[52][spe];
                    int t_rewei = dataspot[53][spe];
                    int t_sensitm = dataspot[54][spe];
                    int t_retrab = dataspot[55][spe];
                    int t_curvactiv = dataspot[56][spe];
                    int t_qualitycurveMethod =  dataspot[57][spe];

                    int t_sensiv = dataspot[58][spe];
                    int t_pastel = dataspot[59][spe];
                    int t_saturated = dataspot[60][spe];
                    int t_proskin = dataspot[61][spe];
                    int t_avoidcsh = dataspot[62][spe];
                    int t_pastsat = dataspot[63][spe];

                    int t_expcomp = dataspot[64][spe];
                    int t_black = dataspot[65][spe];
                    int t_hlcompr = dataspot[66][spe];
                    int t_hlcomprthresh = dataspot[67][spe];
                    int t_shcompr = dataspot[68][spe];
                    int t_sensiex = dataspot[69][spe];

                    int t_centerXbuf = dataspot[70][spe];
                    int t_centerYbuf = dataspot[71][spe];
                    int t_adjblur = dataspot[72][spe];
                    int t_cutpast = dataspot[73][spe];

                    int t_chromacbdl = dataspot[74][spe];

                    int t_lastdust = dataspot[75][spe];
                    int t_blurMethod = dataspot[76][spe];
                    int t_dustMethod = dataspot[77][spe];

                    int t_excludemeth = dataspot[78][spe];
                    int t_sensiexclu = dataspot[79][spe];
                    int t_struc = dataspot[80][spe];
                    int t_warm = dataspot[81][spe];
                    int t_noiselumdetail = dataspot[82][spe];
                    int t_noisechrodetail = dataspot[83][spe];
                    int t_sensiden = dataspot[84][spe];
                    int t_expdenoi = dataspot[85][spe];

                    int t_expcolor = dataspot[86][spe];
                    int t_expvibrance = dataspot[87][spe];
                    int t_expblur = dataspot[88][spe];
                    int t_exptonemap = dataspot[89][spe];
                    int t_expreti = dataspot[90][spe];
                    int t_expsharp = dataspot[91][spe];
                    int t_expcbdl = dataspot[92][spe];
                    int t_expexpose = dataspot[93][spe];

                    int t_bilateral = dataspot[94][spe];
                    int t_noiselequal = dataspot[95][spe];
                    int t_shapemeth = dataspot[96][spe];

                    int t_huerefblur = dataspot[maxdata - 5][spe];
                    int t_hueref = dataspot[maxdata - 4][spe];
                    int t_chromaref = dataspot[maxdata - 3][spe];
                    int t_lumaref = dataspot[maxdata - 2][spe];
                    int t_sobelref = dataspot[maxdata - 1][spe];



                    std::string t_curvret = retistr[spe];
                    std::string t_curvll = llstr[spe];
                    std::string t_curvlh = lhstr[spe];
                    std::string t_curvcc = ccstr[spe];
                    std::string t_curvhh = hhstr[spe];
                    std::string t_curvskin = skinstr[spe];
                    std::string t_psthres = pthstr[spe];
                    std::string t_curvex = exstr[spe];

                    fou << "Mipversion=" << t_mipversion << '@' << endl;
                    fou << "Spot=" << t_sp << '@' << endl;
                    fou << "Circrad=" << t_circrad << '@' << endl;
                    fou << "LocX=" << t_locX << '@' << endl;
                    fou << "LocY=" << t_locY << '@' << endl;
                    fou << "LocYT=" << t_locYT << '@' << endl;
                    fou << "LocXL=" << t_locXL << '@' << endl ;
                    fou << "CenterX=" << t_centerX << '@' << endl;
                    fou << "CenterY=" << t_centerY << '@' << endl;
                    fou << "Lightness=" << t_lightness << '@' << endl;
                    fou << "Contrast=" << t_contrast << '@' <<  endl;
                    fou << "Chroma=" << t_chroma << '@' << endl;
                    fou << "Sensi=" << t_sensi << '@' << endl;
                    fou << "Transit=" << t_transit << '@' << endl;
                    fou << "Invers=" << t_invers << '@' << endl;
                    fou << "Smethod=" << t_Smeth << '@' << endl;
                    fou << "Currentspot=" << t_currentspot << '@' << endl;
                    fou << "Radius=" << t_radius << '@' << endl;
                    fou << "Strength=" << t_strength << '@' << endl;
                    fou << "Sesibn=" << t_sensibn << '@' << endl;
                    fou << "Inversrad=" << t_inversrad << '@' << endl;
                    fou << "Str=" << t_str << '@' << endl;
                    fou << "Chroma=" << t_chrrt << '@' << endl;
                    fou << "Neigh=" << t_neigh << '@' << endl;
                    fou << "Vart=" << t_vart << '@' << endl;
                    fou << "Sensih=" << t_sensih << '@' << endl;
                    fou << "Inversret=" << t_inversret << '@' << endl;
                    fou << "retinexMethod=" << t_retinexMethod << '@' << endl;
                    fou << "Sharradius=" << t_sharradius << '@' << endl;
                    fou << "Sharamount=" << t_sharamount << '@' << endl;
                    fou << "Shardamping=" << t_shardamping << '@' << endl;
                    fou << "Shariter=" << t_shariter << '@' << endl;
                    fou << "Sensisha=" << t_sensisha << '@' << endl;
                    fou << "Inverssha=" << t_inverssha << '@' << endl;
                    fou << "qualityMethod=" << t_qualityMethod << '@' << endl;
                    fou << "Thres=" << t_thres << '@' << endl;
                    fou << "Proxi=" << t_proxi << '@' << endl;
                    fou << "Noiselumf=" << t_noiselumf << '@' << endl;
                    fou << "Noiselumc=" << t_noiselumc << '@' << endl;
                    fou << "Noisechrof=" << t_noisechrof << '@' << endl;
                    fou << "Noisechroc=" << t_noisechroc << '@' << endl;
                    fou << "Mult0=" << t_mult0 << '@' << endl;
                    fou << "Mult1=" << t_mult1 << '@' << endl;
                    fou << "Mult2=" << t_mult2 << '@' << endl;
                    fou << "Mult3=" << t_mult3 << '@' << endl;
                    fou << "Mult4=" << t_mult4 << '@' << endl;
                    fou << "Threshold=" << t_threshold << '@' << endl;
                    fou << "Sensicb=" << t_sensicb << '@' << endl;
                    fou << "Activblurlum=" << t_activlum << '@' << endl;

                    fou << "Stren=" << t_stren << '@' << endl;
                    fou << "Gamma=" << t_gamma << '@' << endl;
                    fou << "Estop=" << t_estop << '@' << endl;
                    fou << "Scaltm=" << t_scaltm << '@' << endl;
                    fou << "Rewei=" << t_rewei << '@' << endl;
                    fou << "Sensitm=" << t_sensitm << '@' << endl;

                    fou << "Retrab=" << t_retrab << '@' << endl;
                    fou << "Curvactiv=" << t_curvactiv << '@' << endl;
                    fou << "qualitycurveMethod=" << t_qualitycurveMethod << '@' << endl;

                    fou << "Sensiv=" << t_sensiv << '@' << endl;
                    fou << "Pastel=" << t_pastel << '@' << endl;
                    fou << "Saturated=" << t_saturated << '@' << endl;
                    fou << "Proskin=" << t_proskin << '@' << endl;
                    fou << "Avoidcsh=" << t_avoidcsh << '@' << endl;
                    fou << "Pastsat=" << t_pastsat << '@' << endl;

                    fou << "Expcomp=" << t_expcomp << '@' << endl;
                    fou << "Black=" << t_black << '@' << endl;
                    fou << "Hlcompr=" << t_hlcompr << '@' << endl;
                    fou << "Hlcomprthresh=" << t_hlcomprthresh << '@' << endl;
                    fou << "Shcompr=" << t_shcompr  << '@' << endl;
                    fou << "Sensiex=" << t_sensiex << '@' << endl;

                    fou << "CenterXbuf=" << t_centerXbuf << '@' << endl;
                    fou << "CenterYbuf=" << t_centerYbuf << '@' << endl;
                    fou << "Adjblur=" << t_adjblur << '@' << endl;
                    fou << "Cutpast=" << t_cutpast << '@' <<  endl;

                    fou << "Chromacbdl=" << t_chromacbdl << '@' <<  endl;
                    fou << "Lastdust=" << t_lastdust << '@' <<  endl;
                    fou << "BlurMethod=" << t_blurMethod << '@' <<  endl;
                    fou << "DustMethod=" << t_dustMethod << '@' <<  endl;

                    fou << "ExcludeMethod=" << t_excludemeth << '@' <<  endl;
                    fou << "Sensiexclu=" << t_sensiexclu << '@' << endl;
                    fou << "Struc=" << t_struc << '@' << endl;
                    fou << "Warm=" << t_warm << '@' << endl;
                    fou << "Noiselumdetail=" << t_noiselumdetail << '@' << endl;
                    fou << "Noisechrodetail=" << t_noisechrodetail << '@' << endl;
                    fou << "Sensiden=" << t_sensiden << '@' << endl;
                    fou << "Expdenoi=" << t_expdenoi << '@' << endl;
                    fou << "Expcolor=" << t_expcolor << '@' << endl;
                    fou << "Expvibrance=" << t_expvibrance << '@' << endl;
                    fou << "Expblur=" << t_expblur << '@' << endl;
                    fou << "Exptonemap=" << t_exptonemap << '@' << endl;
                    fou << "Expreti=" << t_expreti << '@' << endl;
                    fou << "Expsharp=" << t_expsharp << '@' << endl;
                    fou << "Expcbdl=" << t_expcbdl << '@' << endl;
                    fou << "Expexpose=" << t_expexpose << '@' << endl;

                    fou << "Bilateral=" << t_bilateral << '@' << endl;
                    fou << "Noiselequal=" << t_noiselequal << '@' << endl;
                    fou << "ShapeMethod=" << t_shapemeth << '@' <<  endl;

                    fou << "huerefblur=" << t_huerefblur << '@' << endl;
                    fou << "hueref=" << t_hueref << '@' << endl;
                    fou << "chromaref=" << t_chromaref << '@' << endl;
                    fou << "lumaref=" << t_lumaref << '@' << endl;
                    fou << "sobelref=" << t_sobelref << '@' << endl;
                    fou << "curveReti=" << t_curvret << endl;
                    fou << "curveLL=" << t_curvll  << endl;
                    fou << "curveLH=" << t_curvlh  << endl;
                    fou << "curveCC=" << t_curvcc  << endl;
                    fou << "curveHH=" << t_curvhh   << endl;
                    fou << "curveskin=" << t_curvskin  << endl;
                    fou << "pthres=" << t_psthres  << endl;
                    fou << "curveex=" << t_curvex  << endl;

                    fou << endl;
                }

                fou.close();
            }

            //********************************************************
            //end write mip file
            //*********************************************************

            //clean all
            for (int i = 0; i < maxdata; i++) {
                delete [] dataspot[i];
            }

            delete [] dataspot;

            delete [] retistr;
            delete [] llstr;
            delete [] lhstr;
            delete [] ccstr;
            delete [] hhstr;
            delete [] skinstr;
            delete [] pthstr;
            delete [] exstr;

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

        if (shmap) {
            delete shmap;
        }

        shmap = nullptr;

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

        if (params.sh.enabled) {
            shmap = new SHMap(pW, pH, true);
        }

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
bool ImProcCoordinator::getHighQualComputed() {
    // this function may only be called from detail windows
    if(!highQualityComputed) {
        if(options.prevdemo == PD_Sidecar) {
            // we already have high quality preview
            setHighQualComputed();
        } else {
            for (size_t i = 0; i < crops.size() - 1; ++i) { // -1, because last entry is the freshly created detail window
                if (crops[i]->get_skip() == 1 ) {  // there is at least one crop with skip == 1 => we already have high quality preview
                    setHighQualComputed();
                    break;
                }
            }
        }
    }
    return highQualityComputed;
}

void ImProcCoordinator::setHighQualComputed() {
    highQualityComputed = true;
}

}
