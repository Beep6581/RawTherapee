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
#include "simpleprocess.h"
#include "../rtgui/ppversion.h"
#include "colortemp.h"

namespace rtengine {

extern const Settings* settings;

ImProcCoordinator::ImProcCoordinator ()
    : orig_prev(NULL), oprevi(NULL), oprevl(NULL), nprevl(NULL), previmg(NULL), workimg(NULL),
      ncie(NULL), imgsrc(NULL), shmap(NULL), lastAwbEqual(0.), ipf(&params, true), scale(10),
      highDetailPreprocessComputed(false), highDetailRawComputed(false), allocated(false),
      bwAutoR(-9000.f), bwAutoG(-9000.f), bwAutoB(-9000.f), CAMMean(0.f),

      hltonecurve(65536,0),
      shtonecurve(65536,2),//clip above
      tonecurve(65536,0),//,1);

      lumacurve(65536,0),
      chroma_acurve(65536,0),
      chroma_bcurve(65536,0),
      satcurve(65536,0),
      lhskcurve(65536,0),
      clcurve(65536,0),

      vhist16(65536),vhist16bw(65536),
      lhist16(65536), lhist16Cropped(65536),
      lhist16CAM(65536), lhist16CroppedCAM(65536),
      lhist16CCAM(65536),
      lhist16CCAMAF(65536), lhist16ClabAF(65536),
      histCropped(65536),
      lhist16Clad(65536),lhist16CLlad(65536),
      lhist16LClad(65536), lhist16LLClad(65536),
      histRed(256), histRedRaw(256),
      histGreen(256), histGreenRaw(256),
      histBlue(256), histBlueRaw(256),
      histLuma(256),
      histToneCurve(256),
      histToneCurveBW(256),
      histLCurve(256),
      histCCurve(256),
      histCLurve(256),
      histLLCurve(256),

      histLCAM(256),
      histCCAM(256),
      histClad(256),
      bcabhist(256),
      histChroma(256),

      CAMBrightCurveJ(), CAMBrightCurveQ(),

      rCurve(),
      gCurve(),
      bCurve(),
      rcurvehist(256), rcurvehistCropped(256), rbeforehist(256),
      gcurvehist(256), gcurvehistCropped(256), gbeforehist(256),
      bcurvehist(256), bcurvehistCropped(256), bbeforehist(256),

      pW(-1), pH(-1),
      plistener(NULL), imageListener(NULL), aeListener(NULL), hListener(NULL),acListener(NULL), abwListener(NULL),
      resultValid(false), changeSinceLast(0), updaterRunning(false), destroying(false),utili(false),autili(false),
	  butili(false),ccutili(false),cclutili(false),clcutili(false)

    {}

void ImProcCoordinator::assign (ImageSource* imgsrc) {
    this->imgsrc = imgsrc;
}

ImProcCoordinator::~ImProcCoordinator () {

    destroying = true;
    updaterThreadStart.lock ();
    if (updaterRunning && thread)
        thread->join ();
    mProcessing.lock(); 
    mProcessing.unlock(); 
    freeAll ();

    std::vector<Crop*> toDel = crops;
    for (size_t i=0; i<toDel.size(); i++)
        delete toDel[i];

    imgsrc->decreaseRef ();
    updaterThreadStart.unlock ();
}

DetailedCrop* ImProcCoordinator::createCrop  () { 

    return new Crop (this); 
}


// todo: bitmask containing desired actions, taken from changesSinceLast
// cropCall: calling crop, used to prevent self-updates
void ImProcCoordinator::updatePreviewImage (int todo, Crop* cropCall) {

    MyMutex::MyLock processingLock(mProcessing);

    int numofphases = 14;
    int readyphase = 0;

    bwAutoR = bwAutoG = bwAutoB = -9000.f;

    if (todo==CROP && ipf.needsPCVignetting())
        todo |= TRANSFORM; // Change about Crop does affect TRANSFORM

    // Tells to the ImProcFunctions' tools what is the preview scale, which may lead to some simplifications
    ipf.setScale (scale);

    // Check if any detail crops need high detail. If not, take a fast path short cut
    bool highDetailNeeded = (todo & M_HIGHQUAL);
    if (!highDetailNeeded) {
        for (size_t i=0; i<crops.size(); i++)
            if (crops[i]->get_skip() == 1 ) {  // skip=1 -> full resolution
                highDetailNeeded=true;
                break;
            }
    }

    RAWParams rp = params.raw;
    if( !highDetailNeeded ){
        // if below 100% magnification, take a fast path
        rp.dmethod = RAWParams::methodstring[RAWParams::fast];
        rp.hotdeadpix_filt = false;
        rp.ccSteps = 0;
        //rp.all_enhance = false;
        }

    progress ("Applying white balance, color correction & sRGB conversion...",100*readyphase/numofphases);
    // raw auto CA is bypassed if no high detail is needed, so we have to compute it when high detail is needed
    if ( (todo & M_PREPROC) || (!highDetailPreprocessComputed && highDetailNeeded)) {
        imgsrc->preprocess( rp, params.lensProf, params.coarse );
        imgsrc->getRAWHistogram( histRedRaw, histGreenRaw, histBlueRaw );
        if (highDetailNeeded)
            highDetailPreprocessComputed = true;
        else
            highDetailPreprocessComputed = false;
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
        || ( params.toneCurve.hrenabled && params.toneCurve.method!="Color" && imgsrc->IsrgbSourceModified())
        || (!params.toneCurve.hrenabled && params.toneCurve.method=="Color" && imgsrc->IsrgbSourceModified()))
    {

        if (settings->verbose) printf("Demosaic %s\n",rp.dmethod.c_str());

        imgsrc->demosaic( rp );

        if (highDetailNeeded) {
            highDetailRawComputed = true;
            if (params.toneCurve.hrenabled && params.toneCurve.method=="Color") {
                todo |= M_INIT;
            }
        }
        else
            highDetailRawComputed = false;
    }

    // Updating toneCurve.hrenabled if necessary
    // It has to be done there, because the next 'if' statement will use the value computed here
    if (todo & M_AUTOEXP) {
        if (params.toneCurve.autoexp) {// this enabled HLRecovery
            if (ToneCurveParams::HLReconstructionNecessary(histRedRaw, histGreenRaw, histBlueRaw) && !params.toneCurve.hrenabled) {
                // switching params.toneCurve.hrenabled to true -> shouting in listener's ears!
                params.toneCurve.hrenabled=true;

                // forcing INIT to be done, to reconstruct HL again
                todo |= M_INIT;
            }
        }
    }

    if (todo & (M_INIT|M_LINDENOISE)) {
        MyMutex::MyLock initLock(minit);  // Also used in crop window

        imgsrc->HLRecovery_Global( params.toneCurve ); // this handles Color HLRecovery
        if (settings->verbose) printf ("Applying white balance, color correction & sRBG conversion...\n");
        currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);
        if (params.wb.method=="Camera")
            currWB = imgsrc->getWB ();
        else if (params.wb.method=="Auto") {
            if (lastAwbEqual != params.wb.equal) {
                double rm, gm, bm;
                imgsrc->getAutoWBMultipliers(rm, gm, bm);
                if (rm != -1.) {
                    autoWB.update(rm, gm, bm, params.wb.equal);
                    lastAwbEqual = params.wb.equal;
                }
                else {
                    lastAwbEqual = -1.;
                    autoWB.useDefaults(params.wb.equal);
                }
                //double rr,gg,bb;
                //autoWB.getMultipliers(rr,gg,bb);
            }
            currWB = autoWB;
        }
        params.wb.temperature = currWB.getTemp ();
        params.wb.green = currWB.getGreen ();

        int tr = TR_NONE;
        if (params.coarse.rotate==90)  tr |= TR_R90;
        else if (params.coarse.rotate==180) tr |= TR_R180;
        else if (params.coarse.rotate==270) tr |= TR_R270;

        if (params.coarse.hflip)       tr |= TR_HFLIP;
        if (params.coarse.vflip)       tr |= TR_VFLIP;

        imgsrc->getFullSize (fw, fh, tr);
        PreviewProps pp (0, 0, fw, fh, scale);

        // Will (re)allocate the preview's buffers
        setScale (scale);

        imgsrc->getImage (currWB, tr, orig_prev, pp, params.toneCurve, params.icm, params.raw);
        //ColorTemp::CAT02 (orig_prev, &params)	;

        //imgsrc->convertColorSpace(orig_prev, params.icm, params.raw);

        if (todo & M_LINDENOISE) {
            //printf("denoising!\n");
            if (scale==1 && params.dirpyrDenoise.enabled)
                ipf.RGB_denoise(orig_prev, orig_prev, imgsrc->isRAW(), params.dirpyrDenoise, params.defringe, imgsrc->getDirPyrDenoiseExpComp());
        }
        imgsrc->convertColorSpace(orig_prev, params.icm, currWB, params.raw);

        ipf.firstAnalysis (orig_prev, &params, vhist16, imgsrc->getGamma());
    }
    readyphase++;

    progress ("Rotate / Distortion...",100*readyphase/numofphases);
    // Remove transformation if unneeded
    bool needstransform = ipf.needsTransform();
    if (!needstransform && orig_prev!=oprevi) {
        delete oprevi;
        oprevi = orig_prev;
    }
    if (needstransform && orig_prev==oprevi)
        oprevi = new Imagefloat (pW, pH);
    if ((todo & M_TRANSFORM) && needstransform)
        ipf.transform (orig_prev, oprevi, 0, 0, 0, 0, pW, pH, fw, fh, imgsrc->getMetaData()->getFocalLen(),
                       imgsrc->getMetaData()->getFocalLen35mm(), imgsrc->getMetaData()->getFocusDist(), imgsrc->getRotateDegree(), false);

    readyphase++;

    progress ("Preparing shadow/highlight map...",100*readyphase/numofphases);
    if ((todo & M_BLURMAP) && params.sh.enabled) {
        double radius = sqrt (double(pW*pW+pH*pH)) / 2.0;
        double shradius = params.sh.radius;
        if (!params.sh.hq) shradius *= radius / 1800.0;
        shmap->update (oprevi, shradius, ipf.lumimul, params.sh.hq, scale);
    }
    readyphase++;

    if (todo & M_AUTOEXP) {
        if (params.toneCurve.autoexp) {
            LUTu aehist; int aehistcompr;
            imgsrc->getAutoExpHistogram (aehist, aehistcompr);
            ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, params.toneCurve.expcomp,
                    params.toneCurve.brightness, params.toneCurve.contrast, params.toneCurve.black, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh);
            if (aeListener)
                aeListener->autoExpChanged (params.toneCurve.expcomp, params.toneCurve.brightness, params.toneCurve.contrast,
                                            params.toneCurve.black, params.toneCurve.hlcompr,params.toneCurve.hlcomprthresh, params.toneCurve.hrenabled);
        }
    }

    progress ("Exposure curve & CIELAB conversion...",100*readyphase/numofphases);
    if ((todo & M_RGBCURVE) || (todo & M_CROP)) {
        if (hListener) oprevi->calcCroppedHistogram(params, scale, histCropped);

        // complexCurve also calculated pre-curves histogram dependend on crop
        ipf.g = imgsrc->getGamma();
        ipf.iGamma = true;
        CurveFactory::complexCurve (params.toneCurve.expcomp, params.toneCurve.black/65535.0,
                                    params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh,
                                    params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast,
                                    ipf.g, !ipf.iGamma, params.toneCurve.curveMode, params.toneCurve.curve, params.toneCurve.curveMode2, params.toneCurve.curve2,
                                    vhist16, histCropped, hltonecurve, shtonecurve, tonecurve, histToneCurve, customToneCurve1, customToneCurve2, scale==1 ? 1 : 1);

        CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, scale==1 ? 1 : 1);
        CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, scale==1 ? 1 : 1);
        CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, scale==1 ? 1 : 1);

        CurveFactory::curveBW (params.blackwhite.beforeCurve,params.blackwhite.afterCurve, vhist16bw, histToneCurveBW, beforeToneCurveBW, afterToneCurveBW,scale==1 ? 1 : 1);

        //initialize rrm bbm ggm different from zero to avoid black screen in some cases
        double rrm=33.;
        double ggm=33.;
        double bbm=33.;

        // if it's just crop we just need the histogram, no image updates
        if ( todo!=MINUPDATE ) {
            ipf.rgbProc (oprevi, oprevl, hltonecurve, shtonecurve, tonecurve, shmap, params.toneCurve.saturation,
                         rCurve, gCurve, bCurve, customToneCurve1, customToneCurve2,beforeToneCurveBW, afterToneCurveBW, rrm, ggm, bbm, bwAutoR, bwAutoG, bwAutoB, params.toneCurve.expcomp, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh);
            if(params.blackwhite.enabled && params.blackwhite.autoc && abwListener) {
                if (settings->verbose)
                    printf("ImProcCoordinator / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", bwAutoR, bwAutoG, bwAutoB);
                abwListener->BWChanged((float) rrm, (float) ggm, (float) bbm);
            }
            // correct GUI black and white with value
        }

        // compute L channel histogram
        int x1, y1, x2, y2, pos, poscc;
        params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2); 
        lhist16.clear(); lhist16Cropped.clear();
        lhist16Clad.clear(); lhist16CLlad.clear();lhist16LLClad.clear();
        lhist16ClabAF.clear();
        for (int x=0; x<pH; x++)
            for (int y=0; y<pW; y++) {
                pos=CLIP((int)(oprevl->L[x][y]));
                lhist16[pos]++;
                if (y>=y1 && y<y2 && x>=x1 && x<x2) {lhist16Cropped[pos]++;}
            }
 
    }
    readyphase++;
//    utili=false;
//    autili=false;
//    butili=false;
//    ccutili=false;
//    cclutili=false;
//    clcutili=false;
    if ((todo & M_LUMACURVE) || (todo & M_CROP)) {
		utili = false;
        CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, lhist16, lhist16Cropped,
                                     lumacurve, histLCurve, scale==1 ? 1 : 16, utili);
    }
    if (todo & M_LUMACURVE) {
		autili=false;
		butili=false;
		ccutili=false;
		cclutili=false;
		clcutili=false;

        CurveFactory::curveCL(clcutili, params.labCurve.clcurve, clcurve, lhist16CLlad, histCLurve, scale==1 ? 1 : 16);

        CurveFactory::complexsgnCurve (autili, butili,ccutili,cclutili, params.labCurve.chromaticity, params.labCurve.rstprotection,
                                       params.labCurve.acurve, params.labCurve.bcurve,params.labCurve.cccurve,params.labCurve.lccurve,chroma_acurve, chroma_bcurve, satcurve,lhskcurve,
                                       lhist16Clad, lhist16LLClad, histCCurve, histLLCurve, scale==1 ? 1 : 16);
    }
    if (todo & (M_LUMINANCE+M_COLOR) ) {
        nprevl->CopyFrom(oprevl);

        progress ("Applying Color Boost...",100*readyphase/numofphases);

        ipf.chromiLuminanceCurve (pW,nprevl, nprevl, chroma_acurve, chroma_bcurve, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili,cclutili,clcutili, histCCurve, histCLurve, histLLCurve, histLCurve);
        ipf.vibrance(nprevl);
        if((params.colorappearance.enabled && !params.colorappearance.tonecie) ||  (!params.colorappearance.enabled)) ipf.EPDToneMap(nprevl,5,1);
        // for all treatments Defringe, Sharpening, Contrast detail , Microcontrast they are activated if "CIECAM" function are disabled
        readyphase++;
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
            if((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)){
            //if(params.colorappearance.enabled && !params.colorappearance.sharpcie){
                progress ("Pyramid equalizer...",100*readyphase/numofphases);
                ipf.dirpyrequalizer (nprevl);
                readyphase++;
            }
        }

        //L histo  and Chroma histo for ciecam
        // histogram well be for Lab (Lch) values, because very difficult to do with J,Q, M, s, C
        int x1, y1, x2, y2, pos, posc;
        params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
        lhist16CAM.clear(); lhist16CroppedCAM.clear();
        lhist16CCAM.clear();
        lhist16CCAMAF.clear();
        for (int x=0; x<pH; x++)
            for (int y=0; y<pW; y++) {
                pos=CLIP((int)(nprevl->L[x][y]));
                posc=CLIP((int)sqrt(nprevl->a[x][y]*nprevl->a[x][y] + nprevl->b[x][y]*nprevl->b[x][y]));
                if(!params.colorappearance.datacie) lhist16CCAM[posc]++;
                if(!params.colorappearance.datacie)lhist16CAM[pos]++;
                if (y>=y1 && y<y2 && x>=x1 && x<x2) {lhist16CroppedCAM[pos]++;}
            }
        LUTu dummy;
        CurveFactory::curveLightBrightColor (
                params.colorappearance.curveMode, params.colorappearance.curve,
                params.colorappearance.curveMode2, params.colorappearance.curve2,
                params.colorappearance.curveMode3, params.colorappearance.curve3,
                lhist16CAM, lhist16CroppedCAM,histLCAM,
                lhist16CCAM, histCCAM,
                customColCurve1,
                customColCurve2,
                customColCurve3,
                scale==1 ? 1 : 1
        );
        if(params.colorappearance.enabled){
            float fnum = imgsrc->getMetaData()->getFNumber  ();        // F number
            float fiso = imgsrc->getMetaData()->getISOSpeed () ;       // ISO
            float fspeed = imgsrc->getMetaData()->getShutterSpeed () ; // Speed
            double fcomp = imgsrc->getMetaData()->getExpComp  ();      // Compensation +/-
            double adap;
            if(fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) { //if no exif data or wrong
                adap=2000.;
            }
            else {
                double E_V = fcomp + log2 (double((fnum*fnum) / fspeed / (fiso/100.f)));
                E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                E_V += log2(params.raw.expos);// exposure raw white point ; log2 ==> linear to EV
                adap= powf(2.f, E_V-3.f);// cd / m2
                // end calculation adaptation scene luminosity
            }
            int begh=0;
            int endh=pH;
            float d;
            double dd;
            float **buffer = new float*[pH];
            for (int i=0; i<pH; i++)
                buffer[i] = new float[pW];
            bool execsharp=false;
            if(scale==1) execsharp=true;

            if(!ncie)
                ncie = new CieImage (pW, pH);

            if (!CAMBrightCurveJ && (params.colorappearance.algo=="JC" || params.colorappearance.algo=="JS" || params.colorappearance.algo=="ALL"))
                CAMBrightCurveJ(65536,0);
            if (!CAMBrightCurveQ && (params.colorappearance.algo=="QM" || params.colorappearance.algo=="ALL"))
                CAMBrightCurveQ(65536,0);
            if(settings->ciecamfloat){
                ipf.ciecam_02float (ncie, float(adap), begh, endh, pW, 2, nprevl, &params, customColCurve1,customColCurve2,customColCurve3, histLCAM, histCCAM, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, (float**)buffer, execsharp, d);
                if(params.colorappearance.autodegree && acListener && params.colorappearance.enabled) acListener->autoCamChanged(100.*(double)d);
                if(params.colorappearance.autoadapscen && acListener && params.colorappearance.enabled) acListener->adapCamChanged(adap);//real value of adapt scene luminosity
            }
            else {
                ipf.ciecam_02 (ncie, adap, begh, endh, pW, 2, nprevl, &params, customColCurve1,customColCurve2,customColCurve3, histLCAM, histCCAM, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, (float**)buffer, execsharp, dd);
                if(params.colorappearance.autodegree && acListener && params.colorappearance.enabled) acListener->autoCamChanged(100.*dd);
                if(params.colorappearance.autoadapscen && acListener && params.colorappearance.enabled) acListener->adapCamChanged(adap);
            }

            for (int i=0; i<pH; i++)
                delete [] buffer[i];
            delete [] buffer;
            readyphase++;
        }
        else {
            // CIECAM is disabled, we free up its image buffer to save some space
            if (ncie)
                delete ncie; ncie=NULL;

            if (CAMBrightCurveJ) CAMBrightCurveJ.reset();
            if (CAMBrightCurveQ) CAMBrightCurveQ.reset();
        }
    }
    // process crop, if needed
    for (size_t i=0; i<crops.size(); i++)
        if (crops[i]->hasListener () && cropCall != crops[i] )
            crops[i]->update (todo);  // may call ourselves

    // Flagging some LUT as dirty now, whether they have been freed up or not
    CAMBrightCurveJ.dirty = true;
    CAMBrightCurveQ.dirty = true;


    progress ("Conversion to RGB...",100*readyphase/numofphases);
    if (todo!=CROP && todo!=MINUPDATE) {
        MyMutex::MyLock prevImgLock(previmg->getMutex());
        try
        {
            ipf.lab2monitorRgb (nprevl, previmg);
            delete workimg;
            Glib::ustring outProfile=params.icm.output;
            if (params.icm.output=="" || params.icm.output==ColorManagementParams::NoICMString) outProfile="sRGB";
            workimg = ipf.lab2rgb (nprevl, 0,0,pW,pH, outProfile, true);
        }
        catch(char * str)
        {
            progress ("Error converting file...",0);
            return;
        }
    }
    if (!resultValid) {
        resultValid = true;
        if (imageListener)
            imageListener->setImage (previmg, scale, params.crop);
    }
    if (imageListener)
        // TODO: The WB tool should be advertised too in order to get the AutoWB's temp and green values
        imageListener->imageReady (params.crop);

    readyphase++;

    if (hListener) {
        updateLRGBHistograms ();
        hListener->histogramChanged (histRed, histGreen, histBlue, histLuma, histToneCurve, histLCurve,histCCurve, histCLurve, histLLCurve, histLCAM, histCCAM, histRedRaw, histGreenRaw, histBlueRaw, histChroma);
    }
}


void ImProcCoordinator::freeAll () {

    if (settings->verbose) printf ("freeall starts %d\n", (int)allocated);

    if (allocated) {
        if (orig_prev!=oprevi)
        delete oprevi;    oprevi    = NULL;
        delete orig_prev; orig_prev = NULL;
        delete oprevl;    oprevl    = NULL;
        delete nprevl;    nprevl    = NULL;
        if (ncie)
            delete ncie;  ncie      = NULL;

        if (imageListener) {
            imageListener->delImage (previmg);
        }
        else
            delete previmg;

        delete workimg;
        delete shmap;

    }
    allocated = false;
}

/** @brief Handles image buffer (re)allocation and trigger sizeChanged of SizeListener[s]
 * If the scale change, this method will free all buffers and reallocate ones of the new size.
 * It will then tell to the SizeListener that size has changed (sizeChanged)
 *
 * @param prevscale New Preview's scale.
 */
void ImProcCoordinator::setScale (int prevscale) {

if (settings->verbose) printf ("setscale before lock\n");

    tr = TR_NONE;
    if (params.coarse.rotate==90)  tr |= TR_R90;
    if (params.coarse.rotate==180) tr |= TR_R180;
    if (params.coarse.rotate==270) tr |= TR_R270;
    if (params.coarse.hflip)       tr |= TR_HFLIP;
    if (params.coarse.vflip)       tr |= TR_VFLIP;

    int nW, nH;
    imgsrc->getFullSize (fw, fh, tr);

    PreviewProps pp (0, 0, fw, fh, prevscale);
    imgsrc->getSize (tr, pp, nW, nH);

    if (settings->verbose) printf ("setscale starts (%d, %d)\n", nW, nH);

    if (nW!=pW || nH!=pH) {

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
        shmap = new SHMap (pW, pH, true);

        allocated = true;
    }

    scale = prevscale;
    resultValid = false;
    fullw = fw;
    fullh = fh;
    if (settings->verbose) printf ("setscale ends\n");
    if (!sizeListeners.empty())
        for (size_t i=0; i<sizeListeners.size(); i++)
            sizeListeners[i]->sizeChanged (fullw, fullh, fw, fh);
    if (settings->verbose) printf ("setscale ends2\n");

}


void ImProcCoordinator::updateLRGBHistograms () {

    int x1, y1, x2, y2;
    params.crop.mapToResized(pW, pH, scale, x1, x2, y1, y2); 

    histRed.clear();
    histGreen.clear();
    histBlue.clear();

    for (int i=y1; i<y2; i++) {
        int ofs = (i*pW + x1)*3;
        for (int j=x1; j<x2; j++) {
            int r=workimg->data[ofs++];
            int g=workimg->data[ofs++];
            int b=workimg->data[ofs++];

            histRed[r]++;
            histGreen[g]++;
            histBlue[b]++;
        }
    }
    histLuma.clear();
    histChroma.clear();
    for (int i=y1; i<y2; i++)
        for (int j=x1; j<x2; j++) {
            histChroma[(int)(sqrt(nprevl->a[i][j]*nprevl->a[i][j] + nprevl->b[i][j]*nprevl->b[i][j]))/188]++;//188 = 48000/256
            histLuma[(int)(nprevl->L[i][j]/128)]++;
        }

    /*for (int i=0; i<256; i++) {
        Lhist[i] = (int)(256*sqrt(Lhist[i]));
        rhist[i] = (int)(256*sqrt(rhist[i]));
        ghist[i] = (int)(256*sqrt(ghist[i]));
        bhist[i] = (int)(256*sqrt(bhist[i]));
        bcrgbhist[i] = (int)(256*sqrt(bcrgbhist[i]));
        bcLhist[i] = (int)(256*sqrt(bcLhist[i]));
    }*/
}

void ImProcCoordinator::progress (Glib::ustring str, int pr) {

/*  if (plistener) {
    plistener->setProgressStr (str);
    plistener->setProgress ((double)pr / 100.0);
  }*/
}

bool ImProcCoordinator::getAutoWB (double& temp, double& green, double equal) {

    if (imgsrc) {
        if (lastAwbEqual != equal) {
            MyMutex::MyLock lock(minit);  // Also used in crop window
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers(rm, gm, bm);
            if (rm != -1) {
                autoWB.update(rm, gm, bm, equal);
                lastAwbEqual = equal;
            }
            else {
                lastAwbEqual = -1.;
                autoWB.useDefaults(equal);
            }
        }
        temp = autoWB.getTemp ();
        green = autoWB.getGreen ();
        return true;
    }
    else {
        //temp = autoWB.getTemp();
        temp = -1.0;
        green = -1.0;
        return false;
    }
}

void ImProcCoordinator::getCamWB (double& temp, double& green) {

    if (imgsrc) {
        temp = imgsrc->getWB().getTemp ();
        green = imgsrc->getWB().getGreen ();
    }
}

void ImProcCoordinator::getSpotWB (int x, int y, int rect, double& temp, double& tgreen) {

    ColorTemp ret;

    {
    MyMutex::MyLock lock(mProcessing);
    std::vector<Coord2D> points, red, green, blue;
    for (int i=y-rect; i<=y+rect; i++)
        for (int j=x-rect; j<=x+rect; j++) 
            points.push_back (Coord2D (j, i));

    ipf.transCoord (fw, fh, points, red, green, blue);
    int tr = TR_NONE;
    if (params.coarse.rotate==90)  tr |= TR_R90;
    if (params.coarse.rotate==180) tr |= TR_R180;
    if (params.coarse.rotate==270) tr |= TR_R270;
    if (params.coarse.hflip)       tr |= TR_HFLIP;
    if (params.coarse.vflip)       tr |= TR_VFLIP;
    ret = imgsrc->getSpotWB (red, green, blue, tr, params.wb.equal);
    currWB = ColorTemp (params.wb.temperature, params.wb.green,params.wb.equal, params.wb.method);
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

void ImProcCoordinator::getAutoCrop (double ratio, int &x, int &y, int &w, int &h) {

    MyMutex::MyLock lock(mProcessing);

    LCPMapper *pLCPMap=NULL;
    if (params.lensProf.lcpFile.length() && imgsrc->getMetaData()->getFocalLen()>0) {
        LCPProfile *pLCPProf=lcpStore->getProfile(params.lensProf.lcpFile);
        if (pLCPProf) pLCPMap=new LCPMapper(pLCPProf, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(), imgsrc->getMetaData()->getFocusDist(),
            0, false, params.lensProf.useDist, fullw, fullh, params.coarse, imgsrc->getRotateDegree());
    }

    double fillscale = ipf.getTransformAutoFill (fullw, fullh, pLCPMap);
    if (ratio>0) {
        w = fullw * fillscale;
        h = w / ratio;
        if (h > fullh * fillscale) {
            h = fullh * fillscale;
            w = h * ratio;
        }
    }
    else {
        w = fullw * fillscale;
        h = fullh * fillscale;
    }
    x = (fullw - w) / 2;
    y = (fullh - h) / 2;
}


void ImProcCoordinator::saveInputICCReference (const Glib::ustring& fname) {
	
	MyMutex::MyLock lock(mProcessing);
	
	int fW, fH;
	imgsrc->getFullSize (fW, fH, 0);
	PreviewProps pp (0, 0, fW, fH, 1);
	ProcParams ppar = params;
	ppar.toneCurve.hrenabled = false;
	ppar.icm.input = "(none)";
	Imagefloat* im = new Imagefloat (fW, fH);
	imgsrc->preprocess( ppar.raw, ppar.lensProf, ppar.coarse );
	imgsrc->demosaic(ppar.raw );
	//imgsrc->getImage (imgsrc->getWB(), 0, im, pp, ppar.toneCurve, ppar.icm, ppar.raw);
	ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);
	if (params.wb.method=="Camera")
		currWB = imgsrc->getWB ();
	else if (params.wb.method=="Auto") {
		if (lastAwbEqual != params.wb.equal) {
			double rm, gm, bm;
			imgsrc->getAutoWBMultipliers(rm, gm, bm);
			if (rm != -1.) {
				autoWB.update(rm, gm, bm, params.wb.equal);
				lastAwbEqual = params.wb.equal;
			}
			else {
				lastAwbEqual = -1.;
				autoWB.useDefaults(params.wb.equal);
			}
		}
		currWB = autoWB;
	}
	params.wb.temperature = currWB.getTemp ();
	params.wb.green = currWB.getGreen ();
	imgsrc->getImage (currWB, 0, im, pp, ppar.toneCurve, ppar.icm, ppar.raw);
	imgsrc->convertColorSpace(im, ppar.icm, currWB, params.raw);
	Image16* im16 = im->to16();
	delete im;
	im16->saveTIFF (fname,16,true);
	delete im16;
	//im->saveJPEG (fname, 85);
}

void ImProcCoordinator::stopProcessing () {

    updaterThreadStart.lock ();
    if (updaterRunning && thread) {
        changeSinceLast = 0;
        thread->join ();
    }
    updaterThreadStart.unlock ();
}

void ImProcCoordinator::startProcessing () {

    #undef THREAD_PRIORITY_NORMAL

    if (!destroying) {
        updaterThreadStart.lock ();
        if (!updaterRunning) {
            thread = NULL;
            updaterRunning = true;
            updaterThreadStart.unlock ();

            //batchThread->yield(); //the running batch should wait other threads to avoid conflict
            
            thread = Glib::Thread::create(sigc::mem_fun(*this, &ImProcCoordinator::process), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);

        }
        else
            updaterThreadStart.unlock ();
    }
}

void ImProcCoordinator::startProcessing(int changeCode) {
    paramsUpdateMutex.lock();
    changeSinceLast |= changeCode;
    paramsUpdateMutex.unlock();

    startProcessing ();
}

void ImProcCoordinator::process () {

    if (plistener)
        plistener->setProgressState (true);

    paramsUpdateMutex.lock ();
    while (changeSinceLast) {
        params = nextParams;
        int change = changeSinceLast;
        changeSinceLast = 0;
        paramsUpdateMutex.unlock ();

        // M_VOID means no update, and is a bit higher that the rest
        if (change & (M_VOID-1)) updatePreviewImage (change);
        paramsUpdateMutex.lock ();
    }    
    paramsUpdateMutex.unlock ();
    updaterRunning = false;

    if (plistener)
        plistener->setProgressState (false);
}

ProcParams* ImProcCoordinator::beginUpdateParams () {
    paramsUpdateMutex.lock ();

    return &nextParams;
}

void ImProcCoordinator::endUpdateParams (ProcEvent change) {
    endUpdateParams( refreshmap[(int)change] );
}

void ImProcCoordinator::endUpdateParams (int changeFlags) {
    changeSinceLast |= changeFlags;

    paramsUpdateMutex.unlock ();
    startProcessing ();
}


}
