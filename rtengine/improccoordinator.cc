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
#include "iccstore.h"
#ifdef _OPENMP
#include <omp.h>
#endif
namespace rtengine
{

extern const Settings* settings;

ImProcCoordinator::ImProcCoordinator ()
    : orig_prev(nullptr), oprevi(nullptr), oprevl(nullptr), nprevl(nullptr), previmg(nullptr), workimg(nullptr),
      ncie(nullptr), imgsrc(nullptr), shmap(nullptr), lastAwbEqual(0.), lastAwbTempBias(0.0), ipf(&params, true), monitorIntent(RI_RELATIVE),
      softProof(false), gamutCheck(false), scale(10), highDetailPreprocessComputed(false), highDetailRawComputed(false),
      allocated(false), bwAutoR(-9000.f), bwAutoG(-9000.f), bwAutoB(-9000.f), CAMMean(NAN),

      ctColorCurve(),
      hltonecurve(65536),
      shtonecurve(65536),
      tonecurve(65536, 0), //,1);
      lumacurve(32770, 0), // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
      chroma_acurve(65536, 0),
      chroma_bcurve(65536, 0),
      satcurve(65536, 0),
      lhskcurve(65536, 0),
      clcurve(65536, 0),
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
      rcurvehist(256), rcurvehistCropped(256), rbeforehist(256),
      gcurvehist(256), gcurvehistCropped(256), gbeforehist(256),
      bcurvehist(256), bcurvehistCropped(256), bbeforehist(256),
      fw(0), fh(0), tr(0),
      fullw(1), fullh(1),
      pW(-1), pH(-1),
      plistener(nullptr), imageListener(nullptr), aeListener(nullptr), acListener(nullptr), abwListener(nullptr), awbListener(nullptr), actListener(nullptr), adnListener(nullptr), awavListener(nullptr), dehaListener(nullptr), hListener(nullptr),
      resultValid(false), lastOutputProfile("BADFOOD"), lastOutputIntent(RI__COUNT), lastOutputBPC(false), thread(nullptr), changeSinceLast(0), updaterRunning(false), destroying(false), utili(false), autili(false), wavcontlutili(false),
      butili(false), ccutili(false), cclutili(false), clcutili(false), opautili(false), conversionBuffer(1, 1), colourToningSatLimit(0.f), colourToningSatLimitOpacity(0.f)
{}

void ImProcCoordinator::assign (ImageSource* imgsrc)
{
    this->imgsrc = imgsrc;
}

ImProcCoordinator::~ImProcCoordinator ()
{

    destroying = true;
    updaterThreadStart.lock ();

    if (updaterRunning && thread) {
        thread->join ();
    }

    mProcessing.lock();
    mProcessing.unlock();
    freeAll ();

    std::vector<Crop*> toDel = crops;

    for (size_t i = 0; i < toDel.size(); i++) {
        delete toDel[i];
    }

    imgsrc->decreaseRef ();
    updaterThreadStart.unlock ();
}

DetailedCrop* ImProcCoordinator::createCrop  (::EditDataProvider *editDataProvider, bool isDetailWindow)
{

    return new Crop (this, editDataProvider, isDetailWindow);
}


// todo: bitmask containing desired actions, taken from changesSinceLast
// cropCall: calling crop, used to prevent self-updates  ...doesn't seem to be used
void ImProcCoordinator::updatePreviewImage (int todo, Crop* cropCall)
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
            if (crops[i]->get_skip() == 1 ) {  // skip=1 -> full resolution
                highDetailNeeded = true;
                break;
            }
    }

    RAWParams rp = params.raw;
    ColorManagementParams cmp = params.icm;
    LCurveParams  lcur = params.labCurve;

    if( !highDetailNeeded ) {
        // if below 100% magnification, take a fast path
        if(rp.bayersensor.method != RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::none] && rp.bayersensor.method != RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::mono]) {
            rp.bayersensor.method = RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::fast];
        }

        //bayerrp.all_enhance = false;

        if(rp.xtranssensor.method != RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::none] && rp.xtranssensor.method != RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::mono]) {
            rp.xtranssensor.method = RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::fast];
        }

        rp.bayersensor.ccSteps = 0;
        rp.xtranssensor.ccSteps = 0;
        //rp.deadPixelFilter = rp.hotPixelFilter = false;
    }

    progress ("Applying white balance, color correction & sRGB conversion...", 100 * readyphase / numofphases);

    // raw auto CA is bypassed if no high detail is needed, so we have to compute it when high detail is needed
    if ( (todo & M_PREPROC) || (!highDetailPreprocessComputed && highDetailNeeded)) {
        imgsrc->preprocess( rp, params.lensProf, params.coarse );
        imgsrc->getRAWHistogram( histRedRaw, histGreenRaw, histBlueRaw );

        if (highDetailNeeded) {
            highDetailPreprocessComputed = true;
        } else {
            highDetailPreprocessComputed = false;
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

    if (   (todo & M_RAW)
            || (!highDetailRawComputed && highDetailNeeded)
            || ( params.toneCurve.hrenabled && params.toneCurve.method != "Color" && imgsrc->IsrgbSourceModified())
            || (!params.toneCurve.hrenabled && params.toneCurve.method == "Color" && imgsrc->IsrgbSourceModified())) {

        if (settings->verbose) {
            if (imgsrc->getSensorType() == ST_BAYER) {
                printf("Demosaic Bayer image using method: %s\n", rp.bayersensor.method.c_str());
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                printf("Demosaic X-Trans image with using method: %s\n", rp.xtranssensor.method.c_str());
            }
        }

        imgsrc->demosaic( rp);//enabled demosaic
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
        LUTf cdcurve (65536, 0);
        LUTf mapcurve (65536, 0);

        imgsrc->retinexPrepareCurves(params.retinex, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, dehacontlutili, mapcontlutili, useHsl, lhist16RETI, histLRETI);
        float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
        imgsrc->retinex( params.icm, params.retinex,  params.toneCurve, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, conversionBuffer, dehacontlutili, mapcontlutili, useHsl, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax, histLRETI);//enabled Retinex

        if(dehaListener) {
            dehaListener->minmaxChanged(maxCD, minCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
        }
    }


    // Updating toneCurve.hrenabled if necessary
    // It has to be done there, because the next 'if' statement will use the value computed here
    if (todo & M_AUTOEXP) {
        if (params.toneCurve.autoexp) {// this enabled HLRecovery
            if (ToneCurveParams::HLReconstructionNecessary(histRedRaw, histGreenRaw, histBlueRaw) && !params.toneCurve.hrenabled) {
                // switching params.toneCurve.hrenabled to true -> shouting in listener's ears!
                params.toneCurve.hrenabled = true;

                // forcing INIT to be done, to reconstruct HL again
                todo |= M_INIT;
            }
        }
    }

    if (todo & (M_INIT | M_LINDENOISE)) {
        MyMutex::MyLock initLock(minit);  // Also used in crop window

        imgsrc->HLRecovery_Global( params.toneCurve); // this handles Color HLRecovery


        if (settings->verbose) {
            printf ("Applying white balance, color correction & sRBG conversion...\n");
        }

        currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

        if (params.wb.method == "Camera") {
            currWB = imgsrc->getWB ();
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

        params.wb.temperature = currWB.getTemp ();
        params.wb.green = currWB.getGreen ();
        if(params.wb.method == "Auto" && awbListener) {
            awbListener->WBChanged(params.wb.temperature, params.wb.green);
        }

        int tr = getCoarseBitMask(params.coarse);

        imgsrc->getFullSize (fw, fh, tr);

        // Will (re)allocate the preview's buffers
        setScale (scale);
        PreviewProps pp (0, 0, fw, fh, scale);
        // Tells to the ImProcFunctions' tools what is the preview scale, which may lead to some simplifications
        ipf.setScale (scale);

        imgsrc->getImage (currWB, tr, orig_prev, pp, params.toneCurve, params.icm, params.raw);
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

        ipf.firstAnalysis (orig_prev, params, vhist16);
    }

    readyphase++;

    progress ("Rotate / Distortion...", 100 * readyphase / numofphases);
    // Remove transformation if unneeded
    bool needstransform = ipf.needsTransform();

    if (!needstransform && !((todo & (M_TRANSFORM | M_RGBCURVE))  && params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled) && orig_prev != oprevi) {
        delete oprevi;
        oprevi = orig_prev;
    }

    if ((needstransform || ((todo & (M_TRANSFORM | M_RGBCURVE))  && params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled)) ) {
        if(!oprevi || oprevi == orig_prev)
            oprevi = new Imagefloat (pW, pH);
        if (needstransform)
            ipf.transform (orig_prev, oprevi, 0, 0, 0, 0, pW, pH, fw, fh, imgsrc->getMetaData()->getFocalLen(),
                           imgsrc->getMetaData()->getFocalLen35mm(), imgsrc->getMetaData()->getFocusDist(), imgsrc->getRotateDegree(), false);
        else
            orig_prev->copyData(oprevi);
    }

    if ((todo & (M_TRANSFORM | M_RGBCURVE))  && params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled) {
        const int W = oprevi->getWidth();
        const int H = oprevi->getHeight();
        LabImage labcbdl(W, H);
        ipf.rgb2lab(*oprevi, labcbdl, params.icm.working);
        ipf.dirpyrequalizer (&labcbdl, scale);
        ipf.lab2rgb(labcbdl, *oprevi, params.icm.working);
    }

    readyphase++;
    progress ("Preparing shadow/highlight map...", 100 * readyphase / numofphases);

    if ((todo & M_BLURMAP) && params.sh.enabled) {
        double radius = sqrt (double(pW * pW + pH * pH)) / 2.0;
        double shradius = params.sh.radius;

        if (!params.sh.hq) {
            shradius *= radius / 1800.0;
        }

        if(!shmap) {
            shmap = new SHMap (pW, pH, true);
        }

        shmap->update (oprevi, shradius, ipf.lumimul, params.sh.hq, scale);
    }



    readyphase++;

    if (todo & M_AUTOEXP) {
        if (params.toneCurve.autoexp) {
            LUTu aehist;
            int aehistcompr;
            imgsrc->getAutoExpHistogram (aehist, aehistcompr);
            ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, params.toneCurve.expcomp,
                            params.toneCurve.brightness, params.toneCurve.contrast, params.toneCurve.black, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh);

            if (aeListener)
                aeListener->autoExpChanged (params.toneCurve.expcomp, params.toneCurve.brightness, params.toneCurve.contrast,
                                            params.toneCurve.black, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, params.toneCurve.hrenabled);
        }
    }

    progress ("Exposure curve & CIELAB conversion...", 100 * readyphase / numofphases);

    if ((todo & M_RGBCURVE) || (todo & M_CROP)) {
//        if (hListener) oprevi->calcCroppedHistogram(params, scale, histCropped);

        //complexCurve also calculated pre-curves histogram depending on crop
        CurveFactory::complexCurve (params.toneCurve.expcomp, params.toneCurve.black / 65535.0,
                                    params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh,
                                    params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast,
                                    params.toneCurve.curveMode, params.toneCurve.curve, params.toneCurve.curveMode2, params.toneCurve.curve2,
                                    vhist16, hltonecurve, shtonecurve, tonecurve, histToneCurve, customToneCurve1, customToneCurve2, 1);

        CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, 1);
        CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, 1);
        CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, 1);


        opautili = false;

        if(params.colorToning.enabled) {
            TMatrix wprof = iccStore->workingSpaceMatrix (params.icm.working);
            double wp[3][3] = {
                {wprof[0][0], wprof[0][1], wprof[0][2]},
                {wprof[1][0], wprof[1][1], wprof[1][2]},
                {wprof[2][0], wprof[2][1], wprof[2][2]}
            };
            TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params.icm.working);
            double wip[3][3] = {
                {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
                {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
                {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
            };
            params.colorToning.getCurves(ctColorCurve, ctOpacityCurve, wp, wip, opautili);
            CurveFactory::curveToning(params.colorToning.clcurve, clToningcurve, scale == 1 ? 1 : 16);
            CurveFactory::curveToning(params.colorToning.cl2curve, cl2Toningcurve, scale == 1 ? 1 : 16);
        }

        if(params.blackwhite.enabled) {
            CurveFactory::curveBW (params.blackwhite.beforeCurve, params.blackwhite.afterCurve, vhist16bw, histToneCurveBW, beforeToneCurveBW, afterToneCurveBW, 1);
        }

        colourToningSatLimit = float(params.colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
        colourToningSatLimitOpacity = 1.f - (float(params.colorToning.saturatedOpacity) / 100.f);

        int satTH = 80;
        int satPR = 30;
        int indi = 0;

        if(params.colorToning.enabled  && params.colorToning.autosat) { //for colortoning evaluation of saturation settings
            float moyS = 0.f;
            float eqty = 0.f;
            ipf.moyeqt (oprevi, moyS, eqty);//return image : mean saturation and standard dev of saturation
            //printf("moy=%f ET=%f\n", moyS,eqty);
            float satp = ((moyS + 1.5f * eqty) - 0.3f) / 0.7f; //1.5 sigma ==> 93% pixels with high saturation -0.3 / 0.7 convert to Hombre scale

            if(satp >= 0.92f) {
                satp = 0.92f;    //avoid values too high (out of gamut)
            }

            if(satp <= 0.15f) {
                satp = 0.15f;    //avoid too low values
            }

            //satTH=(int) 100.f*satp;
            //satPR=(int) 100.f*(moyS-0.85f*eqty);//-0.85 sigma==>20% pixels with low saturation
            colourToningSatLimit = 100.f * satp;
            satTH = (int) 100.f * satp;

            colourToningSatLimitOpacity = 100.f * (moyS - 0.85f * eqty); //-0.85 sigma==>20% pixels with low saturation
            satPR = (int) 100.f * (moyS - 0.85f * eqty);
        }

        if(actListener) {
            //if(params.blackwhite.enabled) {actListener->autoColorTonChanged(0, satTH, satPR);}
            if(params.blackwhite.enabled && params.colorToning.autosat) {
                actListener->autoColorTonChanged(0, satTH, satPR);    //hide sliders only if autosat
                indi = 0;
            } else {
                if(params.colorToning.autosat) {
                    if     (params.colorToning.method == "Lab") {
                        indi = 1;
                    } else if(params.colorToning.method == "RGBCurves") {
                        indi = 1;
                    } else if(params.colorToning.method == "RGBSliders") {
                        indi = 1;
                    } else if(params.colorToning.method == "Splico") {
                        indi = 2;
                    } else if(params.colorToning.method == "Splitlr") {
                        indi = 2;
                    }

                    //actListener->autoColorTonChanged(indi, satTH, satPR);
                }
            }
        }

        // if it's just crop we just need the histogram, no image updates
        if ( todo & M_RGBCURVE ) {
            //initialize rrm bbm ggm different from zero to avoid black screen in some cases
            double rrm = 33.;
            double ggm = 33.;
            double bbm = 33.;

            DCPProfile::ApplyState as;
            DCPProfile *dcpProf = imgsrc->getDCP(params.icm, currWB, as);

            ipf.rgbProc (oprevi, oprevl, nullptr, hltonecurve, shtonecurve, tonecurve, shmap, params.toneCurve.saturation,
                         rCurve, gCurve, bCurve, colourToningSatLimit , colourToningSatLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, beforeToneCurveBW, afterToneCurveBW, rrm, ggm, bbm, bwAutoR, bwAutoG, bwAutoB, params.toneCurve.expcomp, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, dcpProf, as, histToneCurve);

            if(params.blackwhite.enabled && params.blackwhite.autoc && abwListener) {
                if (settings->verbose) {
                    printf("ImProcCoordinator / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", bwAutoR, bwAutoG, bwAutoB);
                }

                abwListener->BWChanged((float) rrm, (float) ggm, (float) bbm);
            }

            if(params.colorToning.autosat && actListener) {
                if (settings->verbose) {
                    printf("ImProcCoordinator / Auto CT:  indi=%d   satH=%d  satPR=%d\n", indi, (int)colourToningSatLimit , (int) colourToningSatLimitOpacity);
                }

                actListener->autoColorTonChanged(indi, (int) colourToningSatLimit, (int)colourToningSatLimitOpacity);//change sliders autosat
            }

            // correct GUI black and white with value
        }

        // compute L channel histogram
        int x1, y1, x2, y2;
        params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
    }

    readyphase++;

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
        CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, lhist16, lumacurve, histLCurve, scale == 1 ? 1 : 16, utili);
    }

    if (todo & M_LUMACURVE) {

        CurveFactory::curveCL(clcutili, params.labCurve.clcurve, clcurve, scale == 1 ? 1 : 16);

        CurveFactory::complexsgnCurve (autili, butili, ccutili, cclutili, params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,
                                       params.labCurve.lccurve, chroma_acurve, chroma_bcurve, satcurve, lhskcurve, scale == 1 ? 1 : 16);
    }

    if (todo & (M_LUMINANCE + M_COLOR) ) {
        nprevl->CopyFrom(oprevl);

        progress ("Applying Color Boost...", 100 * readyphase / numofphases);
        //   ipf.MSR(nprevl, nprevl->W, nprevl->H, 1);
        histCCurve.clear();
        histLCurve.clear();
        ipf.chromiLuminanceCurve (nullptr, pW, nprevl, nprevl, chroma_acurve, chroma_bcurve, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, histCCurve, histLCurve);
        ipf.vibrance(nprevl);

        if((params.colorappearance.enabled && !params.colorappearance.tonecie) ||  (!params.colorappearance.enabled)) {
            ipf.EPDToneMap(nprevl, 5, 1);
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
        if(params.dirpyrequalizer.cbdlMethod == "aft") {
            if(((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) ) {
                progress ("Pyramid wavelet...", 100 * readyphase / numofphases);
                ipf.dirpyrequalizer (nprevl, scale);
                //ipf.Lanczoslab (ip_wavelet(LabImage * lab, LabImage * dst, const procparams::EqualizerParams & eqparams), nprevl, 1.f/scale);
                readyphase++;
            }
        }


        wavcontlutili = false;
        //CurveFactory::curveWavContL ( wavcontlutili,params.wavelet.lcurve, wavclCurve, LUTu & histogramwavcl, LUTu & outBeforeWavCLurveHistogram,int skip);
        CurveFactory::curveWavContL(wavcontlutili, params.wavelet.wavclCurve, wavclCurve, scale == 1 ? 1 : 16);


        if((params.wavelet.enabled)) {
            WaveletParams WaveParams = params.wavelet;
            //      WaveParams.getCurves(wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY);
            WaveParams.getCurves(wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL);

            int kall = 0;
            progress ("Wavelet...", 100 * readyphase / numofphases);
            //  ipf.ip_wavelet(nprevl, nprevl, kall, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, scale);
            ipf.ip_wavelet(nprevl, nprevl, kall, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL, wavclCurve, wavcontlutili, scale);

        }


        if(params.colorappearance.enabled) {
            //L histo  and Chroma histo for ciecam
            // histogram well be for Lab (Lch) values, because very difficult to do with J,Q, M, s, C
            int x1, y1, x2, y2;
            params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
            lhist16CAM.clear();
            lhist16CCAM.clear();

            if(!params.colorappearance.datacie) {
                for (int x = 0; x < pH; x++)
                    for (int y = 0; y < pW; y++) {
                        int pos = CLIP((int)(nprevl->L[x][y]));
                        int posc = CLIP((int)sqrt(nprevl->a[x][y] * nprevl->a[x][y] + nprevl->b[x][y] * nprevl->b[x][y]));
                        lhist16CAM[pos]++;
                        lhist16CCAM[posc]++;
                    }
            }

            CurveFactory::curveLightBrightColor (params.colorappearance.curve, params.colorappearance.curve2, params.colorappearance.curve3,
                                                 lhist16CAM, histLCAM, lhist16CCAM, histCCAM,
                                                 customColCurve1, customColCurve2, customColCurve3, 1);
            float fnum = imgsrc->getMetaData()->getFNumber  ();        // F number
            float fiso = imgsrc->getMetaData()->getISOSpeed () ;       // ISO
            float fspeed = imgsrc->getMetaData()->getShutterSpeed () ; // Speed
            double fcomp = imgsrc->getMetaData()->getExpComp  ();      // Compensation +/-
            double adap;

            if(fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) { //if no exif data or wrong
                adap = 2000.;
            } else {
                double E_V = fcomp + log2 (double((fnum * fnum) / fspeed / (fiso / 100.f)));
                E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                E_V += log2(params.raw.expos);// exposure raw white point ; log2 ==> linear to EV
                adap = powf(2.f, E_V - 3.f); // cd / m2
                // end calculation adaptation scene luminosity
            }

            int begh = 0;
            int endh = pH;
            float d;
            bool execsharp = false;

            if(!ncie) {
                ncie = new CieImage (pW, pH);
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

            ipf.ciecam_02float (ncie, float(adap), begh, endh, pW, 2, nprevl, &params, customColCurve1, customColCurve2, customColCurve3, histLCAM, histCCAM, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, execsharp, d, scale, 1);

            if(params.colorappearance.autodegree && acListener && params.colorappearance.enabled) {
                acListener->autoCamChanged(100.*(double)d);
            }

            if(params.colorappearance.autoadapscen && acListener && params.colorappearance.enabled) {
                acListener->adapCamChanged(adap);    //real value of adapt scene luminosity
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
    if ((todo & M_MONITOR) || (lastOutputProfile!=params.icm.output) || lastOutputIntent!=params.icm.outputIntent || lastOutputBPC!=params.icm.outputBPC) {
        lastOutputProfile = params.icm.output;
        lastOutputIntent = params.icm.outputIntent;
        lastOutputBPC = params.icm.outputBPC;
        ipf.updateColorProfiles(monitorProfile, monitorIntent, softProof, gamutCheck);
    }

    // process crop, if needed
    for (size_t i = 0; i < crops.size(); i++)
        if (crops[i]->hasListener () && cropCall != crops[i] ) {
            crops[i]->update (todo);    // may call ourselves
        }

    progress ("Conversion to RGB...", 100 * readyphase / numofphases);

    if ((todo != CROP && todo != MINUPDATE) || (todo & M_MONITOR)) {
        MyMutex::MyLock prevImgLock(previmg->getMutex());

        try {
            // Computing the preview image, i.e. converting from WCS->Monitor color space (soft-proofing disabled) or WCS->Output profile->Monitor color space (soft-proofing enabled)
            ipf.lab2monitorRgb (nprevl, previmg);

            // Computing the internal image for analysis, i.e. conversion from WCS->Output profile
            delete workimg;
            workimg = ipf.lab2rgb (nprevl, 0, 0, pW, pH, params.icm);
        } catch(char * str) {
            progress ("Error converting file...", 0);
            return;
        }
    }

    if (!resultValid) {
        resultValid = true;

        if (imageListener) {
            imageListener->setImage (previmg, scale, params.crop);
        }
    }

    if (imageListener)
        // TODO: The WB tool should be advertised too in order to get the AutoWB's temp and green values
    {
        imageListener->imageReady (params.crop);
    }

    readyphase++;

    if (hListener) {
        updateLRGBHistograms ();
        hListener->histogramChanged (histRed, histGreen, histBlue, histLuma, histToneCurve, histLCurve, histCCurve, /*histCLurve, histLLCurve,*/ histLCAM, histCCAM, histRedRaw, histGreenRaw, histBlueRaw, histChroma, histLRETI);
    }
}


void ImProcCoordinator::freeAll ()
{

    if (settings->verbose) {
        printf ("freeall starts %d\n", (int)allocated);
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

        if (ncie) {
            delete ncie;
        }

        ncie      = nullptr;

        if (imageListener) {
            imageListener->delImage (previmg);
        } else {
            delete previmg;
        }

        delete workimg;

        if(shmap) {
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
void ImProcCoordinator::setScale (int prevscale)
{

    if (settings->verbose) {
        printf ("setscale before lock\n");
    }

    tr = getCoarseBitMask(params.coarse);

    int nW, nH;
    imgsrc->getFullSize (fw, fh, tr);

    prevscale++;

    do {
        prevscale--;
        PreviewProps pp (0, 0, fw, fh, prevscale);
        imgsrc->getSize (pp, nW, nH);
    } while(nH < 400 && prevscale > 1 && (nW * nH < 1000000) ); // sctually hardcoded values, perhaps a better choice is possible

    if (settings->verbose) {
        printf ("setscale starts (%d, %d)\n", nW, nH);
    }

    if (nW != pW || nH != pH) {

        freeAll ();

        pW = nW;
        pH = nH;

        orig_prev = new Imagefloat (pW, pH);
        oprevi = orig_prev;
        oprevl = new LabImage (pW, pH);
        nprevl = new LabImage (pW, pH);
        //ncie is only used in ImProcCoordinator::updatePreviewImage, it will be allocated on first use and deleted if not used anymore
        previmg = new Image8 (pW, pH);
        workimg = new Image8 (pW, pH);

        if(params.sh.enabled) {
            shmap = new SHMap (pW, pH, true);
        }

        allocated = true;
    }

    scale = prevscale;
    resultValid = false;
    fullw = fw;
    fullh = fh;

    if (settings->verbose) {
        printf ("setscale ends\n");
    }

    if (!sizeListeners.empty())
        for (size_t i = 0; i < sizeListeners.size(); i++) {
            sizeListeners[i]->sizeChanged (fullw, fullh, fw, fh);
        }

    if (settings->verbose) {
        printf ("setscale ends2\n");
    }

}


void ImProcCoordinator::updateLRGBHistograms ()
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
                    histChroma[(int)(sqrtf(SQR(nprevl->a[i][j]) + SQR(nprevl->b[i][j])) / 188.f)]++; //188 = 48000/256
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

void ImProcCoordinator::progress (Glib::ustring str, int pr)
{

    /*  if (plistener) {
        plistener->setProgressStr (str);
        plistener->setProgress ((double)pr / 100.0);
      }*/
}

bool ImProcCoordinator::getAutoWB (double& temp, double& green, double equal, double tempBias)
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

        temp = autoWB.getTemp ();
        green = autoWB.getGreen ();
        return true;
    } else {
        //temp = autoWB.getTemp();
        temp = -1.0;
        green = -1.0;
        return false;
    }
}

void ImProcCoordinator::getCamWB (double& temp, double& green)
{

    if (imgsrc) {
        temp = imgsrc->getWB().getTemp ();
        green = imgsrc->getWB().getGreen ();
    }
}

void ImProcCoordinator::getSpotWB (int x, int y, int rect, double& temp, double& tgreen)
{

    ColorTemp ret;

    {
        MyMutex::MyLock lock(mProcessing);
        std::vector<Coord2D> points, red, green, blue;

        for (int i = y - rect; i <= y + rect; i++)
            for (int j = x - rect; j <= x + rect; j++) {
                points.push_back (Coord2D (j, i));
            }

        ipf.transCoord (fw, fh, points, red, green, blue);

        int tr = getCoarseBitMask(params.coarse);

        ret = imgsrc->getSpotWB (red, green, blue, tr, params.wb.equal);
        currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);
        //double rr,gg,bb;
        //currWB.getMultipliers(rr,gg,bb);

    } // end of mutex lockong

    if (ret.getTemp() > 0) {
        temp = ret.getTemp ();
        tgreen = ret.getGreen ();
    } else {
        temp = currWB.getTemp ();
        tgreen = currWB.getGreen ();
    }
}

void ImProcCoordinator::getAutoCrop (double ratio, int &x, int &y, int &w, int &h)
{

    MyMutex::MyLock lock(mProcessing);

    LCPMapper *pLCPMap = nullptr;

    if (params.lensProf.lcpFile.length() && imgsrc->getMetaData()->getFocalLen() > 0) {
        LCPProfile *pLCPProf = lcpStore->getProfile(params.lensProf.lcpFile);

        if (pLCPProf) pLCPMap = new LCPMapper(pLCPProf, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(), imgsrc->getMetaData()->getFocusDist(),
                                                  0, false, params.lensProf.useDist, fullw, fullh, params.coarse, imgsrc->getRotateDegree());
    }

    double fillscale = ipf.getTransformAutoFill (fullw, fullh, pLCPMap);

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

void ImProcCoordinator::setMonitorProfile (const Glib::ustring& profile, RenderingIntent intent)
{
    monitorProfile = profile;
    monitorIntent = intent;
}

void ImProcCoordinator::getMonitorProfile (Glib::ustring& profile, RenderingIntent& intent) const
{
    profile = monitorProfile;
    intent = monitorIntent;
}

void ImProcCoordinator::setSoftProofing (bool softProof, bool gamutCheck)
{
    this->softProof = softProof;
    this->gamutCheck = gamutCheck;
}

void ImProcCoordinator::getSoftProofing (bool &softProof, bool &gamutCheck)
{
    softProof = this->softProof;
    gamutCheck = this->gamutCheck;
}

void ImProcCoordinator::saveInputICCReference (const Glib::ustring& fname, bool apply_wb)
{

    MyMutex::MyLock lock(mProcessing);

    int fW, fH;

    int tr = getCoarseBitMask(params.coarse);

    imgsrc->getFullSize (fW, fH, tr);
    PreviewProps pp (0, 0, fW, fH, 1);
    ProcParams ppar = params;
    ppar.toneCurve.hrenabled = false;
    ppar.icm.input = "(none)";
    Imagefloat* im = new Imagefloat (fW, fH);
    imgsrc->preprocess( ppar.raw, ppar.lensProf, ppar.coarse );
    imgsrc->demosaic(ppar.raw );
    ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

    if (params.wb.method == "Camera") {
        currWB = imgsrc->getWB ();
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

    imgsrc->getImage (currWB, tr, im, pp, ppar.toneCurve, ppar.icm, ppar.raw);
    ImProcFunctions ipf (&ppar, true);

    if (ipf.needsTransform()) {
        Imagefloat* trImg = new Imagefloat (fW, fH);
        ipf.transform (im, trImg, 0, 0, 0, 0, fW, fH, fW, fH, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(),
                       imgsrc->getMetaData()->getFocusDist(), imgsrc->getRotateDegree(), true);
        delete im;
        im = trImg;
    }

    if (params.crop.enabled) {
        Imagefloat *tmpim = new Imagefloat (params.crop.w, params.crop.h);
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

    for(int i = 0; i < im->getHeight(); i++) {
        for(int j = 0; j < im->getWidth(); j++) {
            im->r(i, j) = CLIP(im->r(i, j));
            im->g(i, j) = CLIP(im->g(i, j));
            im->b(i, j) = CLIP(im->b(i, j));
        }
    }

    Image16* im16 = im->to16();
    delete im;

    int imw, imh;
    double tmpScale = ipf.resizeScale(&params, fW, fH, imw, imh);

    if (tmpScale != 1.0) {
        Image16* tempImage = new Image16 (imw, imh);
        ipf.resize (im16, tempImage, tmpScale);
        delete im16;
        im16 = tempImage;
    }

    im16->saveTIFF (fname, 16, true);
    delete im16;

    if (plistener) {
        plistener->setProgressState (false);
    }

    //im->saveJPEG (fname, 85);
}

void ImProcCoordinator::stopProcessing ()
{

    updaterThreadStart.lock ();

    if (updaterRunning && thread) {
        changeSinceLast = 0;
        thread->join ();
    }

    updaterThreadStart.unlock ();
}

void ImProcCoordinator::startProcessing ()
{

#undef THREAD_PRIORITY_NORMAL

    if (!destroying) {
        if (!updaterRunning) {
            updaterThreadStart.lock ();
            thread = nullptr;
            updaterRunning = true;
            updaterThreadStart.unlock ();

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

    startProcessing ();
}

void ImProcCoordinator::process ()
{
    if (plistener) {
        plistener->setProgressState (true);
    }

    paramsUpdateMutex.lock ();

    while (changeSinceLast) {
        params = nextParams;
        int change = changeSinceLast;
        changeSinceLast = 0;
        paramsUpdateMutex.unlock ();

        // M_VOID means no update, and is a bit higher that the rest
        if (change & (M_VOID - 1)) {
            updatePreviewImage (change);
        }

        paramsUpdateMutex.lock ();
    }

    paramsUpdateMutex.unlock ();
    updaterRunning = false;

    if (plistener) {
        plistener->setProgressState (false);
    }
}

ProcParams* ImProcCoordinator::beginUpdateParams ()
{
    paramsUpdateMutex.lock ();

    return &nextParams;
}

void ImProcCoordinator::endUpdateParams (ProcEvent change)
{
    endUpdateParams( refreshmap[(int)change] );
}

void ImProcCoordinator::endUpdateParams (int changeFlags)
{
    changeSinceLast |= changeFlags;

    paramsUpdateMutex.unlock ();
    startProcessing ();
}


}
