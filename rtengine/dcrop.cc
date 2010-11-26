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
#include <dcrop.h>
#include <curves.h>
#include <mytime.h>
#include <refreshmap.h>
#define CLIPTO(a,b,c) ((a)>b?((a)<c?(a):c):b)
#define CLIP(a) ((a)<65535 ? (a) : (65535));
#define SKIPS(a,b) ((a) / (b) + ((a) % (b) > 0))

namespace rtengine {

extern Settings* settings;

Crop::Crop (ImProcCoordinator* parent)
    : resizeCrop(NULL), transCrop(NULL), updating(false),
    cropw(-1), croph(-1), trafw(-1), trafh(-1),
    borderRequested(32), cropAllocated(false),
    cropImageListener(NULL), parent(parent),skip(10)
{
    parent->crops.push_back (this);
}

Crop::~Crop () {

    cropMutex.lock ();
    parent->mProcessing.lock (); 
    std::vector<Crop*>::iterator i = std::find (parent->crops.begin(), parent->crops.end(), this);
    if (i!=parent->crops.end ())
        parent->crops.erase (i);
        
    freeAll ();
    parent->mProcessing.unlock (); 
    cropMutex.unlock ();
}

void Crop::setListener (DetailedCropListener* il) { 

    parent->mProcessing.lock(); 
    cropImageListener = il; 
    parent->mProcessing.unlock(); 
}       

void Crop::update (int todo, bool internal) {

    if (!internal)
        parent->mProcessing.lock ();

    ProcParams& params = parent->params;
    cropMutex.lock ();

    if (!params.resize.enabled)
        params.resize.scale = 1.0;
    else if (params.resize.dataspec==1)
        params.resize.scale = (double)params.resize.width / (params.coarse.rotate==90 || params.coarse.rotate==270 ? parent->fh : parent->fw);
    else if (params.resize.dataspec==2)
        params.resize.scale = (double)params.resize.height / (params.coarse.rotate==90 || params.coarse.rotate==270 ? parent->fw : parent->fh);

    parent->ipf.setScale (skip);

    // give possibility to the listener to modify crop window (as the full image dimensions are already known at this point)
    int wx, wy, ww, wh, ws;
    bool overrideWindow = false;
    if (cropImageListener)
        overrideWindow = cropImageListener->getWindow (wx, wy, ww, wh, ws);

    bool regenHighDetail=false;
    if( ws==1 && skip>1 && !parent->fineDetailsProcessed ){
    	regenHighDetail=true;
    }

    // re-allocate sub-images and arrays if their dimensions changed
    bool needsinitupdate = false;
    if (!overrideWindow)
        needsinitupdate = setCropSizes (rqcropx, rqcropy, rqcropw, rqcroph, skip, true);
    else
        needsinitupdate = setCropSizes (wx, wy, ww, wh, ws, true); // this set skip=ws
    // it something has been reallocated, all processing steps have to be performed
    if (needsinitupdate)
        todo = ALL;

    if( regenHighDetail )
    	parent->updatePreviewImage (ALL,this); // We have just set skip to 1

    if (resizeCrop)
        baseCrop = resizeCrop;
    else
        baseCrop = origCrop;

    bool needstransform  = parent->ipf.needsTransform();

    if (todo & M_INIT) {
        parent->minit.lock ();
        int tr = TR_NONE;
        if (params.coarse.rotate==90)  tr |= TR_R90;
        if (params.coarse.rotate==180) tr |= TR_R180;
        if (params.coarse.rotate==270) tr |= TR_R270;
        if (params.coarse.hflip)       tr |= TR_HFLIP;
        if (params.coarse.vflip)       tr |= TR_VFLIP;

        if (!needsinitupdate)
            setCropSizes (rqcropx, rqcropy, rqcropw, rqcroph, skip, true);
        PreviewProps pp (trafx, trafy, trafw*skip, trafh*skip, skip);
        parent->imgsrc->getImage (parent->currWB, tr, origCrop, pp, params.hlrecovery, params.icm, params.raw );

        if (fabs(params.resize.scale-1.0)<1e-7) {
            if (resizeCrop) {
                delete resizeCrop;
                resizeCrop = NULL;
            }
            baseCrop = origCrop;
        }
        else {
            int rcw = trafw*params.resize.scale;
            int rch = trafh*params.resize.scale;
            if (!needstransform) {
                rcw = cropw;
                rch = croph;
            }            
            if (resizeCrop && (resizeCrop->width!=rcw || resizeCrop->height!=rch)) {
                delete resizeCrop;
                resizeCrop = NULL;
            }
            if (!resizeCrop)
                resizeCrop = new Image16 (rcw, rch);
            parent->ipf.resize (origCrop, resizeCrop);
            baseCrop = resizeCrop;
        }
        parent->minit.unlock ();
    }

    // transform
    if ((!needstransform && transCrop) || (transCrop && (transCrop->width!=cropw || transCrop->height!=croph))) {
        delete transCrop;
        transCrop = NULL;
    }
    if (needstransform && !transCrop)
        transCrop = new Image16 (cropw, croph);
    if ((todo & M_TRANSFORM) && needstransform)
    	parent->ipf.transform (baseCrop, transCrop, cropx/skip, cropy/skip, trafx*params.resize.scale/skip, trafy*params.resize.scale/skip, SKIPS(parent->fw,skip), SKIPS(parent->fh,skip));
    if (transCrop)
        baseCrop = transCrop;

    // blurmap for shadow & highlights
    if ((todo & M_BLURMAP) && params.sh.enabled) {
        double radius = sqrt (double(SKIPS(parent->fw,skip)*SKIPS(parent->fw,skip)+SKIPS(parent->fh,skip)*SKIPS(parent->fh,skip))) / 2.0;
        double shradius = radius / 1800.0 * params.sh.radius;
        cshmap->update (baseCrop, (unsigned short**)cbuffer, shradius, parent->ipf.lumimul, params.sh.hq);
        cshmap->forceStat (parent->shmap->max, parent->shmap->min, parent->shmap->avg);
    }

    // shadows & highlights & tone curve & convert to cielab
    if (todo & M_RGBCURVE)
        parent->ipf.rgbProc (baseCrop, laboCrop, parent->hltonecurve, parent->shtonecurve, parent->tonecurve, cshmap, params.toneCurve.saturation);

	
	// apply luminance operations
    if (todo & (M_LUMINANCE+M_COLOR)) {
        parent->ipf.luminanceCurve (laboCrop, labnCrop, parent->lumacurve, 0, croph);
		parent->ipf.chrominanceCurve (laboCrop, labnCrop, 0, parent->chroma_acurve, 0, croph);
		parent->ipf.chrominanceCurve (laboCrop, labnCrop, 1, parent->chroma_bcurve, 0, croph);

		parent->ipf.colorCurve (labnCrop, labnCrop);

        if (skip==1) {
			parent->ipf.impulsedenoise (labnCrop);
			parent->ipf.defringe (labnCrop);
            parent->ipf.lumadenoise (labnCrop, cbuffer);
            parent->ipf.colordenoise (labnCrop, cbuffer);
			parent->ipf.dirpyrdenoise (labnCrop);
			parent->ipf.sharpening (labnCrop, (unsigned short**)cbuffer);
			parent->ipf.dirpyrequalizer (labnCrop);
            parent->ipf.waveletEqualizer(labnCrop, true, true);
        }

    }
    

    // switch back to rgb
    parent->ipf.lab2rgb (labnCrop, cropImg);

    if (cropImageListener) {
        int finalW = rqcropw;
        if (cropImg->getWidth()-leftBorder < finalW)
            finalW = cropImg->getWidth()-leftBorder;
        int finalH = rqcroph;
        if (cropImg->getHeight()-upperBorder < finalH)
            finalH = cropImg->getHeight()-upperBorder;
            
        Image8* final = new Image8 (finalW, finalH);
        for (int i=0; i<finalH; i++)
            memcpy (final->data + 3*i*finalW, cropImg->data + 3*(i+upperBorder)*cropw + 3*leftBorder, 3*finalW);
        cropImageListener->setDetailedCrop (final, params.crop, rqcropx, rqcropy, rqcropw, rqcroph, skip);
        delete final;
    }

    cropMutex.unlock ();

    if (!internal)
        parent->mProcessing.unlock ();
}

void Crop::freeAll () {

    if (settings->verbose) printf ("freeallcrop starts %d\n", (int)cropAllocated);

    if (cropAllocated) {
        delete origCrop;
        if (transCrop)
            delete transCrop;
        transCrop = NULL;        
        if (resizeCrop)
            delete resizeCrop;
        resizeCrop = NULL;
        delete laboCrop;
        delete labnCrop;
        delete cropImg;
        delete cshmap;
        for (int i=0; i<croph; i++)
            delete [] cbuffer[i];
        delete [] cbuffer;
    }
    cropAllocated = false;
}

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
    int rqx1 = CLIPTO(rqcropx,0,parent->fullw-1);
    int rqy1 = CLIPTO(rqcropy,0,parent->fullh-1);
    int rqx2 = rqx1 + rqcropw - 1;
    int rqy2 = rqy1 + rqcroph - 1;
    rqx2 = CLIPTO(rqx2,0,parent->fullw-1);
    rqy2 = CLIPTO(rqy2,0,parent->fullh-1);

    this->skip = skip;
    
    // add border, if possible
    int bx1 = rqx1 - skip*borderRequested;
    int by1 = rqy1 - skip*borderRequested;
    int bx2 = rqx2 + skip*borderRequested;
    int by2 = rqy2 + skip*borderRequested;
    // clip it to fit into image area
    bx1 = CLIPTO(bx1,0,parent->fullw-1);
    by1 = CLIPTO(by1,0,parent->fullh-1);
    bx2 = CLIPTO(bx2,0,parent->fullw-1);
    by2 = CLIPTO(by2,0,parent->fullh-1);
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

    if (settings->verbose) printf ("setsizes starts (%d, %d, %d, %d)\n", orW, orH, trafw, trafh);

    if (cw!=cropw || ch!=croph || orW!=trafw || orH!=trafh) {

        freeAll ();
    
        cropw = cw;
        croph = ch;
        trafw = orW;
        trafh = orH;

        origCrop = new Image16 (trafw, trafh);
        laboCrop = new LabImage (cropw, croph);    
        labnCrop = new LabImage (cropw, croph);    
        cropImg = new Image8 (cropw, croph);
        cshmap = new SHMap (cropw, croph, true);
        
        cbuffer = new int*[croph];
        for (int i=0; i<croph; i++)
            cbuffer[i] = new int[cropw];

        resizeCrop = NULL;
        transCrop = NULL;
        
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

void Crop::fullUpdate () { 

    if (updating) {
        needsNext = true;
        return;
    }

    updating = true;

    parent->updaterThreadStart.lock ();
    if (parent->updaterRunning && parent->thread) {
		// Do NOT reset changes here, since in a long chain of events it will lead to chroma_scale not being updated,
		// causing ImProcFunctions::lab2rgb to return a black image on some opens
        //parent->changeSinceLast = 0;
        parent->thread->join ();
    }

    if (parent->plistener)
        parent->plistener->setProgressState (1);

    needsNext = true;
    while (needsNext) {
        needsNext = false;
        update (ALL, true); 
    }
    updating = false;

    if (parent->plistener)
        parent->plistener->setProgressState (0);

    parent->updaterThreadStart.unlock ();
}

}

