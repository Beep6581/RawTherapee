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
#include "dcrop.h"
#include "curves.h"
#include "mytime.h"
#include "refreshmap.h"
#include "rt_math.h"
#include "colortemp.h"

#define SKIPS(a,b) ((a) / (b) + ((a) % (b) > 0))

namespace rtengine {

extern const Settings* settings;

Crop::Crop (ImProcCoordinator* parent)
    : origCrop(NULL), transCrop(NULL), laboCrop(NULL), labnCrop(NULL),
      cropImg(NULL), cieCrop(NULL), cbuf_real(NULL), cshmap(NULL),
      cbuffer(NULL), updating(false), newUpdatePending(false),
      skip(10),
      cropx(0), cropy(0), cropw(-1), croph(-1),
      trafx(0), trafy(0), trafw(-1), trafh(-1),
      rqcropx(0), rqcropy(0), rqcropw(-1), rqcroph(-1),
      borderRequested(32), upperBorder(0), leftBorder(0),
      cropAllocated(false),
      cropImageListener(NULL), parent(parent)
{
    parent->crops.push_back (this);
}

Crop::~Crop () {

    MyMutex::MyLock cropLock(cropMutex);
    MyMutex::MyLock processingLock(parent->mProcessing);

    std::vector<Crop*>::iterator i = std::find (parent->crops.begin(), parent->crops.end(), this);
    if (i!=parent->crops.end ())
        parent->crops.erase (i);
    freeAll ();
}

void Crop::destroy () {
    MyMutex::MyLock lock(cropMutex);
    MyMutex::MyLock processingLock(parent->mProcessing);   /////////         RETESTER MAINTENANT QUE CE VERROU EST AJOUTE!!!
    freeAll();
}

void Crop::setListener (DetailedCropListener* il) { 
    // We can make reads in the IF, because the mProcessing lock is only needed for change
    if (cropImageListener!=il) {
        MyMutex::MyLock lock(cropMutex);
        cropImageListener = il;
    }
}

void Crop::update (int todo) {
    MyMutex::MyLock lock(cropMutex);

    ProcParams& params = parent->params;

    // No need to update todo here, since it has already been changed in ImprocCoordinator::updatePreviewImage,
    // and Crop::update ask to do ALL anyway

    // give possibility to the listener to modify crop window (as the full image dimensions are already known at this point)
    int wx, wy, ww, wh, ws;
    bool overrideWindow = false;
    if (cropImageListener)
        overrideWindow = cropImageListener->getWindow (wx, wy, ww, wh, ws);

    // re-allocate sub-images and arrays if their dimensions changed
    bool needsinitupdate = false;
    if (!overrideWindow)
        needsinitupdate = setCropSizes (rqcropx, rqcropy, rqcropw, rqcroph, skip, true);
    else
        needsinitupdate = setCropSizes (wx, wy, ww, wh, ws, true); // this set skip=ws
    // it something has been reallocated, all processing steps have to be performed
    if (needsinitupdate || (todo & M_HIGHQUAL))
        todo = ALL;

    // Tells to the ImProcFunctions' tool what is the preview scale, which may lead to some simplifications
    parent->ipf.setScale (skip);

    Imagefloat* baseCrop = origCrop;

    bool needstransform  = parent->ipf.needsTransform();

    if (todo & (M_INIT|M_LINDENOISE)) {
        MyMutex::MyLock lock(parent->minit);  // Also used in improccoord

        int tr = TR_NONE;
        if (params.coarse.rotate==90)  tr |= TR_R90;
        else if (params.coarse.rotate==180) tr |= TR_R180;
        else if (params.coarse.rotate==270) tr |= TR_R270;

        if (params.coarse.hflip)       tr |= TR_HFLIP;
        if (params.coarse.vflip)       tr |= TR_VFLIP;

        if (!needsinitupdate)
            setCropSizes (rqcropx, rqcropy, rqcropw, rqcroph, skip, true);
        PreviewProps pp (trafx, trafy, trafw*skip, trafh*skip, skip);
        parent->imgsrc->getImage (parent->currWB, tr, origCrop, pp, params.hlrecovery, params.icm, params.raw );
        //ColorTemp::CAT02 (origCrop, &params)	;

        //parent->imgsrc->convertColorSpace(origCrop, params.icm);

        if (todo & M_LINDENOISE) {
            if (skip==1 && params.dirpyrDenoise.enabled)
                parent->ipf.RGB_denoise(origCrop, origCrop, parent->imgsrc->isRAW(), /*Roffset,*/ params.dirpyrDenoise, params.defringe, parent->imgsrc->getDirPyrDenoiseExpComp());
        }
        parent->imgsrc->convertColorSpace(origCrop, params.icm, parent->currWB, params.raw);
    }

    // transform
    if (needstransform) {
        if (!transCrop)
            transCrop = new Imagefloat (cropw, croph);

        if ((todo & M_TRANSFORM) && needstransform)
            parent->ipf.transform (baseCrop, transCrop, cropx/skip, cropy/skip, trafx/skip, trafy/skip, SKIPS(parent->fw,skip), SKIPS(parent->fh,skip), parent->getFullWidth(), parent->getFullHeight(),
                                   parent->imgsrc->getMetaData()->getFocalLen(), parent->imgsrc->getMetaData()->getFocalLen35mm(),
                                   parent->imgsrc->getMetaData()->getFocusDist(), parent->imgsrc->getRotateDegree(), false);
        if (transCrop)
            baseCrop = transCrop;
    }
    else {
        if (transCrop) delete transCrop;
        transCrop = NULL;
    }

    // blurmap for shadow & highlights
    if ((todo & M_BLURMAP) && params.sh.enabled) {
        double radius = sqrt (double(SKIPS(parent->fw,skip)*SKIPS(parent->fw,skip)+SKIPS(parent->fh,skip)*SKIPS(parent->fh,skip))) / 2.0;
        double shradius = params.sh.radius;
        if (!params.sh.hq) shradius *= radius / 1800.0;
        cshmap->update (baseCrop, shradius, parent->ipf.lumimul, params.sh.hq, skip);
        cshmap->forceStat (parent->shmap->max_f, parent->shmap->min_f, parent->shmap->avg);
    }

    // shadows & highlights & tone curve & convert to cielab
    /*int xref,yref;
    xref=000;yref=000;
    if (colortest && cropw>115 && croph>115)
        for(int j=1;j<5;j++){
            xref+=j*30;yref+=j*30;
            if (settings->verbose) printf("before rgbProc RGB Xr%i Yr%i Skip=%d  R=%f  G=%f  B=%f gamma=%f  \n",xref,yref,skip,
                   baseCrop->r[(int)(xref/skip)][(int)(yref/skip)]/256,
                   baseCrop->g[(int)(xref/skip)][(int)(yref/skip)]/256,
                   baseCrop->b[(int)(xref/skip)][(int)(yref/skip)]/256,
                   parent->imgsrc->getGamma());
        }*/
    double rrm, ggm, bbm;

    if (todo & M_RGBCURVE)
        parent->ipf.rgbProc (baseCrop, laboCrop, parent->hltonecurve, parent->shtonecurve, parent->tonecurve, cshmap,
                             params.toneCurve.saturation, parent->rCurve, parent->gCurve, parent->bCurve, parent->customToneCurve1,
                             parent->customToneCurve2, parent->beforeToneCurveBW, parent->afterToneCurveBW,rrm, ggm, bbm,
                             parent->bwAutoR, parent->bwAutoG, parent->bwAutoB);

    /*xref=000;yref=000;
    if (colortest && cropw>115 && croph>115)
    for(int j=1;j<5;j++){
        xref+=j*30;yref+=j*30;
        if (settings->verbose) {
            printf("after rgbProc RGB Xr%i Yr%i Skip=%d  R=%f  G=%f  B=%f  \n",xref,yref,skip,
                   baseCrop->r[(int)(xref/skip)][(int)(yref/skip)]/256,
                   baseCrop->g[(int)(xref/skip)][(int)(yref/skip)]/256,
                   baseCrop->b[(int)(xref/skip)][(int)(yref/skip)]/256);
            printf("after rgbProc Lab Xr%i Yr%i Skip=%d  l=%f  a=%f  b=%f  \n",xref,yref,skip,
                   laboCrop->L[(int)(xref/skip)][(int)(yref/skip)]/327,
                   laboCrop->a[(int)(xref/skip)][(int)(yref/skip)]/327,
                   laboCrop->b[(int)(xref/skip)][(int)(yref/skip)]/327);
        }
    }*/

    // apply luminance operations
    if (todo & (M_LUMINANCE+M_COLOR)) {
        //I made a little change here. Rather than have luminanceCurve (and others) use in/out lab images, we can do more if we copy right here.
        labnCrop->CopyFrom(laboCrop);


        //parent->ipf.luminanceCurve (labnCrop, labnCrop, parent->lumacurve);
        bool utili=true;
        bool autili=true;
        bool butili=true;
        bool ccutili=true;
        bool clcutili=true;
        bool cclutili=true;

 
        LUTu dummy;
        parent->ipf.chromiLuminanceCurve (1,labnCrop, labnCrop, parent->chroma_acurve, parent->chroma_bcurve, parent->satcurve, parent->lhskcurve,  parent->clcurve, parent->lumacurve, utili, autili, butili, ccutili,cclutili, clcutili, dummy, dummy, dummy, dummy);
        parent->ipf.vibrance (labnCrop);
        if((params.colorappearance.enabled && !params.colorappearance.tonecie) ||  (!params.colorappearance.enabled)) parent->ipf.EPDToneMap(labnCrop,5,1);
        //parent->ipf.EPDToneMap(labnCrop, 5, 1);    //Go with much fewer than normal iterates for fast redisplay.
        // for all treatments Defringe, Sharpening, Contrast detail , Microcontrast they are activated if "CIECAM" function are disabled
        if (skip==1) {
            if((params.colorappearance.enabled && !settings->autocielab)  || (!params.colorappearance.enabled)) {
            parent->ipf.impulsedenoise (labnCrop);}
            if((params.colorappearance.enabled && !settings->autocielab) ||(!params.colorappearance.enabled) ) {parent->ipf.defringe (labnCrop);}
            parent->ipf.MLsharpen (labnCrop);
            if((params.colorappearance.enabled && !settings->autocielab)  || (!params.colorappearance.enabled)) {
                    parent->ipf.MLmicrocontrast (labnCrop);
                    parent->ipf.sharpening (labnCrop, (float**)cbuffer);
                    parent->ipf.dirpyrequalizer (labnCrop);
                    }
        }

        if(params.colorappearance.enabled){
            float fnum = parent->imgsrc->getMetaData()->getFNumber  ();        // F number
            float fiso = parent->imgsrc->getMetaData()->getISOSpeed () ;       // ISO
            float fspeed = parent->imgsrc->getMetaData()->getShutterSpeed () ; // Speed
            double fcomp = parent->imgsrc->getMetaData()->getExpComp  ();      // Compensation +/-
            double adap; // Scene's luminosity adaptation factor
            if(fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) { //if no exif data or wrong
                adap=2000.;
            }
            else {
                double E_V = fcomp + log2 (double((fnum*fnum) / fspeed / (fiso/100.f)));
                E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                E_V += log2(params.raw.expos);// exposure raw white point ; log2 ==> linear to EV
                adap= pow(2., E_V-3.);// cd / m2
                // end calculation adaptation scene luminosity
            }

            int begh = 0, endh = labnCrop->H;
            bool execsharp=false;
            if(skip==1) execsharp=true;

            if (!cieCrop)
                cieCrop = new CieImage (cropw, croph);

            if(settings->ciecamfloat) {
                float d; // not used after this block
                parent->ipf.ciecam_02float (cieCrop, float(adap), begh, endh, 1, 2, labnCrop, &params, parent->customColCurve1, parent->customColCurve2, parent->customColCurve3,
                                            dummy, dummy, parent->CAMBrightCurveJ, parent->CAMBrightCurveQ, parent->CAMMean, 5, 1,(float**)cbuffer, execsharp, d);
            }
            else {
                double dd; // not used after this block
                parent->ipf.ciecam_02 (cieCrop,adap, begh, endh, 1, 2, labnCrop, &params, parent->customColCurve1, parent->customColCurve2, parent->customColCurve3,
                                       dummy, dummy, parent->CAMBrightCurveJ, parent->CAMBrightCurveQ, parent->CAMMean, 5, 1,(float**)cbuffer, execsharp, dd);
            }
        }
        else {
            // CIECAM is disbaled, we free up its image buffer to save some space
            if (cieCrop) delete cieCrop; cieCrop=NULL;
        }
    }
    // switch back to rgb
    parent->ipf.lab2monitorRgb (labnCrop, cropImg);

    //parent->ipf.lab2monitorRgb (laboCrop, cropImg);

    //cropImg = baseCrop->to8();
    /*
    //     int xref,yref;
    xref=000;yref=000;
    if (colortest && cropw>115 && croph>115)
    for(int j=1;j<5;j++){
        xref+=j*30;yref+=j*30;
        int rlin = (CurveFactory::igamma2((float)cropImg->data[3*((int)(xref/skip)*cropImg->width+(int)(yref/skip))]/255.0) * 255.0);
        int glin = (CurveFactory::igamma2((float)cropImg->data[3*((int)(xref/skip)*cropImg->width+(int)(yref/skip))+1]/255.0) * 255.0);
        int blin = (CurveFactory::igamma2((float)cropImg->data[3*((int)(xref/skip)*cropImg->width+(int)(yref/skip))+2]/255.0) * 255.0);

        printf("after lab2rgb RGB lab2 Xr%i Yr%i Skip=%d  R=%d  G=%d  B=%d  \n",xref,yref,skip,
               rlin,glin,blin);
               //cropImg->data[3*((int)(xref/skip)*cropImg->width+(int)(yref/skip))],
               //cropImg->data[(3*((int)(xref/skip)*cropImg->width+(int)(yref/skip))+1)],
               //cropImg->data[(3*((int)(xref/skip)*cropImg->width+(int)(yref/skip))+2)]);
        //printf("after lab2rgb Lab lab2 Xr%i Yr%i Skip=%d  l=%f  a=%f  b=%f  \n",xref,yref,skip, labnCrop->L[(int)(xref/skip)][(int)(yref/skip)]/327,labnCrop->a[(int)(xref/skip)][(int)(yref/skip)]/327,labnCrop->b[(int)(xref/skip)][(int)(yref/skip)]/327);
        printf("after lab2rgb Lab Xr%i Yr%i Skip=%d  l=%f  a=%f  b=%f  \n",xref,yref,skip,
               labnCrop->L[(int)(xref/skip)][(int)(yref/skip)]/327,
               labnCrop->a[(int)(xref/skip)][(int)(yref/skip)]/327,
               labnCrop->b[(int)(xref/skip)][(int)(yref/skip)]/327)q;
    }
    */
    /*
    if (colortest && cropImg->height>115 && cropImg->width>115) {//for testing
        xref=000;yref=000;
        printf("dcrop final R= %d  G= %d  B= %d  \n",
               cropImg->data[3*xref/(skip)*(cropImg->width+1)],
               cropImg->data[3*xref/(skip)*(cropImg->width+1)+1],
               cropImg->data[3*xref/(skip)*(cropImg->width+1)+2]);
    }
    */
    if (cropImageListener) {
        // this in output space held in parallel to allow analysis like shadow/highlight
        Glib::ustring outProfile=params.icm.output;
        if (params.icm.output=="" || params.icm.output==ColorManagementParams::NoICMString) outProfile="sRGB";
        Image8 *cropImgtrue = parent->ipf.lab2rgb (labnCrop, 0,0,cropw,croph, outProfile, false);

        int finalW = rqcropw;
        if (cropImg->getWidth()-leftBorder < finalW)
            finalW = cropImg->getWidth()-leftBorder;
        int finalH = rqcroph;
        if (cropImg->getHeight()-upperBorder < finalH)
            finalH = cropImg->getHeight()-upperBorder;

        Image8* final = new Image8 (finalW, finalH);
        Image8* finaltrue = new Image8 (finalW, finalH);
        for (int i=0; i<finalH; i++) {
            memcpy (final->data + 3*i*finalW, cropImg->data + 3*(i+upperBorder)*cropw + 3*leftBorder, 3*finalW);
            memcpy (finaltrue->data + 3*i*finalW, cropImgtrue->data + 3*(i+upperBorder)*cropw + 3*leftBorder, 3*finalW);
        }
        cropImageListener->setDetailedCrop (final, finaltrue, params.icm, params.crop, rqcropx, rqcropy, rqcropw, rqcroph, skip);
        delete final;
        delete finaltrue;
        delete cropImgtrue;
    }
}

void Crop::freeAll () {

    if (settings->verbose) printf ("freeallcrop starts %d\n", (int)cropAllocated);

    if (cropAllocated) {
        if (origCrop ) { delete    origCrop;  origCrop=NULL; }
        if (transCrop) { delete    transCrop; transCrop=NULL; }
        if (laboCrop ) { delete    laboCrop;  laboCrop=NULL; }
        if (labnCrop ) { delete    labnCrop;  labnCrop=NULL; }
        if (cropImg  ) { delete    cropImg;   cropImg=NULL; }
        if (cieCrop  ) { delete    cieCrop;   cieCrop=NULL; }
        if (cbuf_real) { delete [] cbuf_real; cbuf_real=NULL; }
        if (cbuffer  ) { delete [] cbuffer;   cbuffer=NULL; }
        if (cshmap   ) { delete    cshmap;    cshmap=NULL; }
    }
    cropAllocated = false;
}

/** @brief Handles crop's image buffer reallocation and trigger sizeChanged of SizeListener[s]
 * If the scale changes, this method will free all buffers and reallocate ones of the new size.
 * It will then tell to the SizeListener that size has changed (sizeChanged)
 */
bool Crop::setCropSizes (int rcx, int rcy, int rcw, int rch, int skip, bool internal) {

if (settings->verbose) printf ("setcropsizes before lock\n");

    if (!internal)
        cropMutex.lock ();

    bool changed = false;

    rqcropx = rcx;
    rqcropy = rcy;
    rqcropw = rcw;
    rqcroph = rch;

    // store and set requested crop size
    int rqx1 = LIM(rqcropx,0,parent->fullw-1);
    int rqy1 = LIM(rqcropy,0,parent->fullh-1);
    int rqx2 = rqx1 + rqcropw - 1;
    int rqy2 = rqy1 + rqcroph - 1;
    rqx2 = LIM(rqx2,0,parent->fullw-1);
    rqy2 = LIM(rqy2,0,parent->fullh-1);

    this->skip = skip;
    
    // add border, if possible
    int bx1 = rqx1 - skip*borderRequested;
    int by1 = rqy1 - skip*borderRequested;
    int bx2 = rqx2 + skip*borderRequested;
    int by2 = rqy2 + skip*borderRequested;
    // clip it to fit into image area
    bx1 = LIM(bx1,0,parent->fullw-1);
    by1 = LIM(by1,0,parent->fullh-1);
    bx2 = LIM(bx2,0,parent->fullw-1);
    by2 = LIM(by2,0,parent->fullh-1);
    int bw = bx2 - bx1 + 1;
    int bh = by2 - by1 + 1;
    
    // determine which part of the source image is required to compute the crop rectangle
    int orx, ory, orw, orh;
    ProcParams& params = parent->params;
    parent->ipf.transCoord (parent->fw, parent->fh, bx1, by1, bw, bh, orx, ory, orw, orh);
        
    int tr = TR_NONE;
    if (params.coarse.rotate==90)  tr |= TR_R90;
    if (params.coarse.rotate==180) tr |= TR_R180;
    if (params.coarse.rotate==270) tr |= TR_R270;
    if (params.coarse.hflip)       tr |= TR_HFLIP;
    if (params.coarse.vflip)       tr |= TR_VFLIP;

    PreviewProps cp (orx, ory, orw, orh, skip);
    int orW, orH;
    parent->imgsrc->getSize (tr, cp, orW, orH);

    int cw = SKIPS(bw,skip);
    int ch = SKIPS(bh,skip);

    leftBorder  = SKIPS(rqx1-bx1,skip);
    upperBorder = SKIPS(rqy1-by1,skip);

    if (settings->verbose)
        printf ("setsizes starts (%d, %d, %d, %d, %d, %d)\n", orW, orH, trafw, trafh,cw,ch);

    if (cw!=cropw || ch!=croph || orW!=trafw || orH!=trafh) {

        freeAll ();
    
        cropw = cw;
        croph = ch;
        trafw = orW;
        trafh = orH;

        origCrop = new Imagefloat (trafw, trafh);
        //transCrop will be allocated later, if necessary
        laboCrop = new LabImage (cropw, croph);
        labnCrop = new LabImage (cropw, croph);
        cropImg = new Image8 (cropw, croph);
        //cieCrop is only used in Crop::update, it will be allocated on first use and deleted if not used anymore
        cbuffer = new float*[croph];
        cbuf_real= new float[(croph+2)*cropw];
        for (int i=0; i<croph; i++)
            cbuffer[i] = cbuf_real+cropw*i+cropw;
        cshmap = new SHMap (cropw, croph, true);

        cropAllocated = true;

        changed = true;
    }

    cropx = bx1;
    cropy = by1;
    trafx = orx;
    trafy = ory;

    if (settings->verbose) printf ("setsizes ends\n");

    if (!internal)
        cropMutex.unlock ();

    return changed;
}

/** @brief Look out if a new thread has to be started to process the update
  *
  * @return If true, a new updating thread has to be created. If false, the current updating thread will be used
  */
bool Crop::tryUpdate() {
    bool needsNewThread = true;

    if (updating) {
        // tells to the updater thread that a new update is pending
        newUpdatePending = true;
        // no need for a new thread, the current one will do the job
        needsNewThread = false;
    }
    else
        // the crop is now being updated ...well, when fullUpdate will be called
        updating = true;

    return needsNewThread;
}

/* @brief Handles Crop updating in its own thread
 *
 * This method will cycle updates as long as Crop::newUpdatePending will be true. During the processing,
 * intermediary update will be automatically flushed by Crop::tryUpdate.
 *
 * This method is called when the visible part of the crop has changed (resize, zoom, etc..), so it needs a full update
 */
void Crop::fullUpdate () {

    parent->updaterThreadStart.lock ();
    if (parent->updaterRunning && parent->thread) {
        // Do NOT reset changes here, since in a long chain of events it will lead to chroma_scale not being updated,
        // causing Color::lab2rgb to return a black image on some opens
        //parent->changeSinceLast = 0;
        parent->thread->join ();
    }
    parent->updaterThreadStart.unlock ();

    if (parent->plistener)
        parent->plistener->setProgressState (true);

    // If there are more update request, the following WHILE will collect it
    newUpdatePending = true;
    while (newUpdatePending) {
        newUpdatePending = false;
        update (ALL); 
    }
    updating = false;  // end of crop update

    if (parent->plistener)
        parent->plistener->setProgressState (false);
}

}

