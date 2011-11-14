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
#include <improccoordinator.h>
#include <curves.h>
#include <mytime.h>
#include <refreshmap.h>
#include <simpleprocess.h>
#define CLIPTO(a,b,c) ((a)>b?((a)<c?(a):c):b)
#define CLIP(a) ((a)>0?((a)<65535?(a):65535):0)

namespace rtengine {

extern const Settings* settings;

ImProcCoordinator::ImProcCoordinator ()
    : awbComputed(false), ipf(&params, true), scale(10), allocated(false),
    pW(-1), pH(-1), plistener(NULL), lastHighDetail(false),
    imageListener(NULL), aeListener(NULL), hListener(NULL), resultValid(false),
    changeSinceLast(0), updaterRunning(false), destroying(false), workimg(NULL) {

    hltonecurve(65536,0);
    shtonecurve(65536,2);//clip above
    tonecurve(65536,0);//,1);

    lumacurve(65536,0);
    chroma_acurve(65536,0);
    chroma_bcurve(65536,0);
	satcurve(65536,0);

    vhist16(65536);
    lhist16(65536); lhist16Cropped(65536);
    histCropped(65536);

    histRed(256); histRedRaw(256);
    histGreen(256); histGreenRaw(256);
    histBlue(256); histBlueRaw(256);
    histLuma(256);
    histToneCurve(256);
    histLCurve(256);
    bcabhist(256);

}

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
    for (int i=0; i<toDel.size(); i++)
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

    mProcessing.lock ();

    int numofphases = 14;
    int readyphase = 0;

    ipf.setScale (scale);

	// Check if any detail crops need high detail. If not, take a fast path short cut
    bool highDetailNeeded = (todo & M_HIGHQUAL);
    if (!highDetailNeeded) {
	for (int i=0; i<crops.size(); i++)
		    if (crops[i]->get_skip() == 1 ) {  // skip=1 -> full resolution
			highDetailNeeded=true;
			break;
		}
    }

	RAWParams rp = params.raw;
	if( !highDetailNeeded ){
        // if below 100% magnification, take a fast path
		rp.dmethod = RAWParams::methodstring[RAWParams::fast];
		rp.ca_autocorrect  = false;
		rp.hotdeadpix_filt = false;
		rp.ccSteps = 0;
		rp.all_enhance = false;
	}

    progress ("Applying white balance, color correction & sRGB conversion...",100*readyphase/numofphases);
    if ( todo & M_PREPROC) {
    	imgsrc->preprocess( rp );
        imgsrc->getRAWHistogram( histRedRaw, histGreenRaw, histBlueRaw );
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
    if ((todo & M_RAW)
    	|| (!lastHighDetail && highDetailNeeded)
    	|| (params.hlrecovery.enabled && params.hlrecovery.method!="Color" && imgsrc->IsrgbSourceModified())
    	|| (!params.hlrecovery.enabled && params.hlrecovery.method=="Color" && imgsrc->IsrgbSourceModified())){

    	if (settings->verbose) printf("Demosaic %s\n",rp.dmethod.c_str());
    	imgsrc->demosaic( rp );
    }
    lastHighDetail=highDetailNeeded;


    if (todo & M_INIT) {
        Glib::Mutex::Lock lock(minit);  // Also used in crop window

        imgsrc->HLRecovery_Global( params.hlrecovery ); // this handles Color HLRecovery

        if (settings->verbose) printf ("Applying white balance, color correction & sRBG conversion...\n");
        currWB = ColorTemp (params.wb.temperature, params.wb.green);
        if (params.wb.method=="Camera")
            currWB = imgsrc->getWB ();
        else if (params.wb.method=="Auto") {
            if (!awbComputed) {
                autoWB = imgsrc->getAutoWB ();
                awbComputed = true;
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
        setScale (scale);
        imgsrc->getImage (currWB, tr, orig_prev, pp, params.hlrecovery, params.icm, params.raw);
        ipf.firstAnalysis (orig_prev, &params, vhist16, imgsrc->getGamma());
    }
    readyphase++;

    progress ("Rotate / Distortion...",100*readyphase/numofphases);
    bool needstransform = ipf.needsTransform();
    // Remove transformation if unneeded
    if (!needstransform && orig_prev!=oprevi) {
        delete oprevi;
        oprevi = orig_prev;
    }
    if (needstransform && orig_prev==oprevi)
        oprevi = new Imagefloat (pW, pH);
    if ((todo & M_TRANSFORM) && needstransform)
    	ipf.transform (orig_prev, oprevi, 0, 0, 0, 0, pW, pH);

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
            ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, params.toneCurve.expcomp, params.toneCurve.black);
            if (aeListener)
                aeListener->autoExpChanged (params.toneCurve.expcomp, params.toneCurve.black);
        }
    }

    progress ("Exposure curve & CIELAB conversion...",100*readyphase/numofphases);
    if ((todo & M_RGBCURVE) || todo==CROP) {
        if (hListener) oprevi->calcCroppedHistogram(params, scale, histCropped);

        // complexCurve also calculated pre-curves histogram dependend on crop
        CurveFactory::complexCurve (params.toneCurve.expcomp, params.toneCurve.black/65535.0, \
									params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, \
									params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast, \
									imgsrc->getGamma(), true, params.toneCurve.curve, \
									vhist16, histCropped, hltonecurve, shtonecurve, tonecurve, histToneCurve, scale==1 ? 1 : 1);
        
        // if it's just crop we just need the histogram, no image updates
        if ( todo!=CROP ) {
            ipf.rgbProc (oprevi, oprevl, hltonecurve, shtonecurve, tonecurve, shmap, params.toneCurve.saturation);
        }

        // compute L channel histogram
        int x1, y1, x2, y2, pos;
        params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2); 

        lhist16.clear(); lhist16Cropped.clear();
        for (int x=0; x<pH; x++)
            for (int y=0; y<pW; y++) {
                pos=CLIP((int)(oprevl->L[x][y]));
                lhist16[pos]++;

                if (y>=y1 && y<y2 && x>=x1 && x<x2) lhist16Cropped[pos]++;
            }
 
    }
    readyphase++;

    if ((todo & M_LUMACURVE) || todo==CROP) {

        CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, lhist16, lhist16Cropped,
                                     lumacurve, histLCurve, scale==1 ? 1 : 16);
    }

    if (todo & M_LUMACURVE) {
		CurveFactory::complexsgnCurve (params.labCurve.saturation, params.labCurve.enable_saturationlimiter, params.labCurve.saturationlimit, \
									   params.labCurve.acurve, params.labCurve.bcurve, chroma_acurve, chroma_bcurve, satcurve, scale==1 ? 1 : 16);
	}
	
    if (todo & (M_LUMINANCE+M_COLOR) ) {
        progress ("Applying Luminance Curve...",100*readyphase/numofphases);

        ipf.luminanceCurve (oprevl, nprevl, lumacurve);

        readyphase++;
		progress ("Applying Color Boost...",100*readyphase/numofphases);
		ipf.chrominanceCurve (oprevl, nprevl, chroma_acurve, chroma_bcurve, satcurve/*, params.labCurve.saturation*/);
        //ipf.colorCurve (nprevl, nprevl);
		ipf.vibrance(nprevl);
        readyphase++;
		if (scale==1) {
            progress ("Denoising luminance impulse...",100*readyphase/numofphases);
            ipf.impulsedenoise (nprevl);
            readyphase++;
			progress ("Defringing...",100*readyphase/numofphases);
            ipf.defringe (nprevl);
            readyphase++;
            progress ("Denoising luma/chroma...",100*readyphase/numofphases);
            ipf.dirpyrdenoise (nprevl);
            readyphase++;
			if (params.sharpenEdge.enabled) {
                progress ("Edge sharpening...",100*readyphase/numofphases);
				ipf.MLsharpen (nprevl);
		        readyphase++;
			}
			if (params.sharpenMicro.enabled) {
                progress ("Microcontrast...",100*readyphase/numofphases);			
				ipf.MLmicrocontrast (nprevl);
		        readyphase++;
			}
            if (params.sharpening.enabled) {
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

            progress ("Pyramid equalizer...",100*readyphase/numofphases);
            ipf.dirpyrequalizer (nprevl);
            readyphase++;
        }
    }

    // process crop, if needed
    for (int i=0; i<crops.size(); i++)
        if (crops[i]->hasListener () && cropCall != crops[i] )
            crops[i]->update (todo);  // may call outselves

    progress ("Conversion to RGB...",100*readyphase/numofphases);
    if (todo!=CROP) {
        previmg->getMutex().lock();
        try
        {
            ipf.lab2rgb (nprevl, previmg);
            delete workimg;
			workimg = ipf.lab2rgb (nprevl, 0,0,pW,pH, params.icm.working);        
        }
        catch(char * str)
        {
            progress ("Error converting file...",0);
            previmg->getMutex().unlock();
            mProcessing.unlock ();
            return;
        }
        previmg->getMutex().unlock();
    }   
    if (!resultValid) {
        resultValid = true;
        if (imageListener)
            imageListener->setImage (previmg, scale, params.crop);
    }
    if (imageListener)
        imageListener->imageReady (params.crop);

    readyphase++;

    if (hListener) {
        updateLRGBHistograms ();
        hListener->histogramChanged (histRed, histGreen, histBlue, histLuma, histToneCurve, histLCurve, histRedRaw, histGreenRaw, histBlueRaw);
    }

    mProcessing.unlock ();
}


void ImProcCoordinator::freeAll () {

    if (settings->verbose) printf ("freeall starts %d\n", (int)allocated);

    if (allocated) {
        if (orig_prev!=oprevi)
            delete oprevi;
        delete orig_prev;
        delete oprevl;
        delete nprevl;
        
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
    if (sizeListeners.size()>0)
        for (int i=0; i<sizeListeners.size(); i++)
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
    for (int i=y1; i<y2; i++)
        for (int j=x1; j<x2; j++) {
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

void ImProcCoordinator::getAutoWB (double& temp, double& green) {

    if (imgsrc) {
        if (!awbComputed) {
            minit.lock ();
            autoWB = imgsrc->getAutoWB ();
            minit.unlock ();
            awbComputed = true;
        }
        temp = autoWB.getTemp ();
        green = autoWB.getGreen ();
    }
}

void ImProcCoordinator::getCamWB (double& temp, double& green) {

    if (imgsrc) {
        temp = imgsrc->getWB().getTemp ();
        green = imgsrc->getWB().getGreen ();
    }
}

void ImProcCoordinator::getSpotWB (int x, int y, int rect, double& temp, double& tgreen) {

    mProcessing.lock ();
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
    
    ColorTemp ret = imgsrc->getSpotWB (red, green, blue, tr);
	currWB = ColorTemp (params.wb.temperature, params.wb.green);
    mProcessing.unlock ();

	if (ret.getTemp() > 0) {
		temp = ret.getTemp ();
		tgreen = ret.getGreen ();
	} else {
		temp = currWB.getTemp ();
		tgreen = currWB.getGreen ();
	}
}

void ImProcCoordinator::getAutoCrop (double ratio, int &x, int &y, int &w, int &h) {

    mProcessing.lock ();

    double fillscale = ipf.getTransformAutoFill (fullw, fullh);
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

    mProcessing.unlock ();
}


void ImProcCoordinator::saveInputICCReference (const Glib::ustring& fname) {
	
	mProcessing.lock ();
	
	int fW, fH;
	imgsrc->getFullSize (fW, fH, 0);
	PreviewProps pp (0, 0, fW, fH, 1);
	ProcParams ppar = params;
	ppar.hlrecovery.enabled = false;
	ppar.icm.input = "(none)";
	Imagefloat* im = new Imagefloat (fW, fH);
	Image16* im16 = new Image16 (fW, fH);
	imgsrc->preprocess( ppar.raw );
	imgsrc->demosaic(ppar.raw );
	//imgsrc->getImage (imgsrc->getWB(), 0, im, pp, ppar.hlrecovery, ppar.icm, ppar.raw);
	ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green);
	if (params.wb.method=="Camera")
		currWB = imgsrc->getWB ();
	else if (params.wb.method=="Auto") {
		if (!awbComputed) {
			autoWB = imgsrc->getAutoWB ();
			awbComputed = true;
		}
		currWB = autoWB;
	}
	params.wb.temperature = currWB.getTemp ();
	params.wb.green = currWB.getGreen ();
	imgsrc->getImage (currWB, 0, im, pp, ppar.hlrecovery, ppar.icm, ppar.raw);
	im16 = im->to16();
	im16->saveTIFF (fname,16,true);
	//im->saveJPEG (fname, 85);
	mProcessing.unlock ();
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
