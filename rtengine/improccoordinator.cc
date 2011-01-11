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
#define CLIP(a) ((a)<65535 ? (a) : (65535));

namespace rtengine {

extern Settings* settings;

ImProcCoordinator::ImProcCoordinator ()
    : awbComputed(false), ipf(&params, true), scale(10), allocated(false),
    pW(-1), pH(-1), plistener(NULL),fineDetailsProcessed(false),
    imageListener(NULL), aeListener(NULL), hListener(NULL), resultValid(false),
    changeSinceLast(0), updaterRunning(false), destroying(false) {

	dummy1        = new float[65536];
    dummy2        = new float[65536];
    hltonecurve   = new float[65536];
    shtonecurve   = new float[65536];
    tonecurve     = new int[65536];

    lumacurve     = new int[65536];
    chroma_acurve = new int[65536];
    chroma_bcurve = new int[65536];

    vhist16 = new unsigned int[65536];
    lhist16 = new unsigned int[65536];

    rhist     = new unsigned int[256];
    ghist     = new unsigned int[256];
    bhist     = new unsigned int[256];
    Lhist     = new unsigned int[256];
    bcrgbhist = new unsigned int[256];
    bcLhist   = new unsigned int[256];
    //bcabhist  = new unsigned int[256];

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

    delete [] dummy1;
    delete [] dummy2;
    delete [] hltonecurve;
    delete [] shtonecurve;
    delete [] tonecurve;

    delete [] lumacurve;
    delete [] chroma_acurve;
    delete [] chroma_bcurve;

    delete [] vhist16;
    delete [] lhist16;

    delete [] rhist;
    delete [] ghist;
    delete [] bhist;
    delete [] Lhist;
    delete [] bcrgbhist;
    delete [] bcLhist;
    //delete [] bcabhist;

    imgsrc->decreaseRef ();
    updaterThreadStart.unlock ();
}

DetailedCrop* ImProcCoordinator::createCrop  () { 

    return new Crop (this); 
}

void ImProcCoordinator::updatePreviewImage (int todo, Crop* cropCall) {

    mProcessing.lock ();

    int numofphases = 10;
    int readyphase = 0;

    ipf.setScale (scale);

    bool highDetailNeeded=false;
	for (int i=0; i<crops.size(); i++)
		if (crops[i]->get_skip() == 1 ){
			highDetailNeeded=true;
			break;
		}
	RAWParams rp = params.raw;
	if( !highDetailNeeded ){
		rp.dmethod = RAWParams::methodstring[RAWParams::fast];
		rp.ca_autocorrect  = false;
		rp.hotdeadpix_filt = false;
		rp.ccSteps = 0;
	}
    progress ("Applying white balance, color correction & sRBG conversion...",100*readyphase/numofphases);
    if ( todo & M_PREPROC)
    	imgsrc->preprocess( rp );
    if( todo & M_RAW){
    	if( !highDetailNeeded ){
    		fineDetailsProcessed = false;
    	}else
    		fineDetailsProcessed = true;
    	imgsrc->demosaic( rp );
    }
    if (todo & M_INIT) {
        minit.lock ();
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
        if (params.coarse.rotate==180) tr |= TR_R180;
        if (params.coarse.rotate==270) tr |= TR_R270;
        if (params.coarse.hflip)       tr |= TR_HFLIP;
        if (params.coarse.vflip)       tr |= TR_VFLIP;

        imgsrc->getFullSize (fw, fh, tr);
        PreviewProps pp (0, 0, fw, fh, scale);
        setScale (scale, true);
        imgsrc->getImage (currWB, tr, orig_prev, pp, params.hlrecovery, params.icm, params.raw);
        ipf.firstAnalysis (orig_prev, &params, vhist16, imgsrc->getGamma());
        minit.unlock ();
    }
    readyphase++;

    progress ("Rotate / Distortion...",100*readyphase/numofphases);
    bool needstransform = ipf.needsTransform();
    if (!needstransform && orig_prev!=oprevi) {
        delete oprevi;
        oprevi = orig_prev;
    }
    if (needstransform && orig_prev==oprevi)
        oprevi = new Image16 (pW, pH);
    if ((todo & M_TRANSFORM) && needstransform)
    	ipf.transform (orig_prev, oprevi, 0, 0, 0, 0, pW, pH);

    readyphase++;

    progress ("Preparing shadow/highlight map...",100*readyphase/numofphases);
    if ((todo & M_BLURMAP) && params.sh.enabled) {
        double radius = sqrt (double(pW*pW+pH*pH)) / 2.0;
        double shradius = radius / 1800.0 * params.sh.radius;
        shmap->update (oprevi, (unsigned short**)buffer, shradius, ipf.lumimul, params.sh.hq);
    }
    readyphase++;

    if (todo & M_AUTOEXP) {
        if (params.toneCurve.autoexp) {
            unsigned int *aehist = new unsigned int[65536]; int aehistcompr;
            imgsrc->getAEHistogram (aehist, aehistcompr);
            ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, params.toneCurve.expcomp, params.toneCurve.black);
            if (aeListener)
                aeListener->autoExpChanged (params.toneCurve.expcomp, params.toneCurve.black);
            delete [] aehist;
        }
    }

    progress ("Exposure curve & CIELAB conversion...",100*readyphase/numofphases);
    if (todo & M_RGBCURVE) {
        CurveFactory::complexCurve (params.toneCurve.expcomp, params.toneCurve.black/65535.0, \
									params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, \
									params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast, \
									imgsrc->getDefGain(), imgsrc->getGamma(), true, params.toneCurve.curve, vhist16, \
									hltonecurve, shtonecurve, tonecurve, bcrgbhist, scale==1 ? 1 : 1);
        ipf.rgbProc (oprevi, oprevl, hltonecurve, shtonecurve, tonecurve, shmap, params.toneCurve.saturation);

        // recompute luminance histogram
        memset (lhist16, 0, 65536*sizeof(int));
        for (int i=0; i<pH; i++)
            for (int j=0; j<pW; j++)
                lhist16[oprevl->L[i][j]]++;
    }
    readyphase++;

    if (todo & M_LUMACURVE) {
        CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, 0.0, params.labCurve.brightness, params.labCurve.contrast, 0.0, 0.0, false, \
									params.labCurve.lcurve, lhist16, dummy1, dummy2, lumacurve, bcLhist, scale==1 ? 1 : 16);
		CurveFactory::complexsgnCurve (0.0, 100.0, params.labCurve.saturation, 1.0, params.labCurve.acurve, chroma_acurve, scale==1 ? 1 : 16);
		CurveFactory::complexsgnCurve (0.0, 100.0, params.labCurve.saturation, 1.0, params.labCurve.bcurve, chroma_bcurve, scale==1 ? 1 : 16);
	}
	
	
    if (todo & (M_LUMINANCE+M_COLOR) ) {
        progress ("Applying Luminance Curve...",100*readyphase/numofphases);
        ipf.luminanceCurve (oprevl, nprevl, lumacurve, 0, pH);

        readyphase++;
		progress ("Applying Color Boost...",100*readyphase/numofphases);
		ipf.chrominanceCurve (oprevl, nprevl, 0, chroma_acurve, 0, pH);
        ipf.chrominanceCurve (oprevl, nprevl, 1, chroma_bcurve, 0, pH);
        ipf.colorCurve (nprevl, nprevl);

        readyphase++;
		if (scale==1) {
            progress ("Denoising luminance impulse...",100*readyphase/numofphases);
            ipf.impulsedenoise (nprevl);
        }
		if (scale==1) {
			progress ("Defringing...",100*readyphase/numofphases);
            ipf.defringe (nprevl);
		}
        if (scale==1) {
            progress ("Denoising luminance...",100*readyphase/numofphases);
            ipf.lumadenoise (nprevl, buffer);
        }
        readyphase++;
		if (scale==1) {
            progress ("Denoising color...",100*readyphase/numofphases);
            ipf.colordenoise (nprevl, buffer);
        }
		if (scale==1) {
            progress ("Denoising luma/chroma...",100*readyphase/numofphases);
            ipf.dirpyrdenoise (nprevl);
        }
		if (scale==1) {
            progress ("Sharpening...",100*readyphase/numofphases);
            ipf.sharpening (nprevl, (unsigned short**)buffer);
        }
        readyphase++;
		if (scale==1) {
            progress ("Pyramid equalizer...",100*readyphase/numofphases);
            ipf.dirpyrequalizer (nprevl);
        }
		if (scale==1) {
            progress ("Wavelet...",100*readyphase/numofphases);
            ipf.waveletEqualizer (nprevl, true, true);
        }
		

    }

    // process crop, if needed
    for (int i=0; i<crops.size(); i++)
        if (crops[i]->hasListener () && cropCall != crops[i] )
            crops[i]->update (todo, true);

    progress ("Conversion to RGB...",100*readyphase/numofphases);
    if (todo!=CROP) {
        previmg->getMutex().lock();
        try
        {
            ipf.lab2rgb (nprevl, previmg);
        }
        catch(char * str)
        {
           progress ("Error converting file...",0);
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
        int hx1 = 0, hx2 = pW, hy1 = 0, hy2 = pH;
        if (params.crop.enabled) {
            hx1 = MIN(pW-1,MAX(0,params.crop.x / scale));
            hy1 = MIN(pH-1,MAX(0,params.crop.y / scale));   
            hx2 = MIN(pW,MAX(0,(params.crop.x+params.crop.w) / scale)); 
            hy2 = MIN(pH,MAX(0,(params.crop.y+params.crop.h) / scale));
        }
        updateHistograms (hx1, hy1, hx2, hy2);
        hListener->histogramChanged (rhist, ghist, bhist, Lhist, bcrgbhist, bcLhist);
    }

    progress ("Ready",100*readyphase/numofphases);
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
        delete shmap;
        for (int i=0; i<pH; i++)
            delete [] buffer[i];
        delete [] buffer;
    }
    allocated = false;
}

void ImProcCoordinator::setScale (int prevscale, bool internal) {

if (settings->verbose) printf ("setscale before lock\n");

    if (!internal)
        mProcessing.lock ();

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
        
        orig_prev = new Image16 (pW, pH);
        oprevi = orig_prev;
        oprevl = new LabImage (pW, pH);    
        nprevl = new LabImage (pW, pH);    
        previmg = new Image8 (pW, pH);
        shmap = new SHMap (pW, pH, true);
        
        buffer = new int*[pH];
        for (int i=0; i<pH; i++)
            buffer[i] = new int[pW];
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

    if (!internal)
        mProcessing.unlock ();
}


void ImProcCoordinator::updateHistograms (int x1, int y1, int x2, int y2) {

    memset (rhist, 0, 256*sizeof(int));
    memset (ghist, 0, 256*sizeof(int));
    memset (bhist, 0, 256*sizeof(int));
	//memset (bcrgbhist, 0, 256*sizeof(int));

    for (int i=y1; i<y2; i++) {
        int ofs = (i*pW + x1)*3;
        for (int j=x1; j<x2; j++) {
			int r=previmg->data[ofs++];
			int g=previmg->data[ofs++];
			int b=previmg->data[ofs++];
			
			//bcrgbhist[(int)(0.299*r + 0.587*g + 0.114*b)]++;
            rhist[r]++;
            ghist[g]++;
            bhist[b]++;
        }
    }

    memset (Lhist, 0, 256*sizeof(int));
    for (int i=y1; i<y2; i++)
        for (int j=x1; j<x2; j++) {
            Lhist[nprevl->L[i][j]/256]++;
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

void ImProcCoordinator::fullUpdatePreviewImage () {

    if (destroying)
        return;

    updaterThreadStart.lock ();
    if (updaterRunning && thread) {
        changeSinceLast = 0;
        thread->join ();
    }

    if (plistener)
        plistener->setProgressState (1);

    updatePreviewImage (ALL); 

    if (plistener)
        plistener->setProgressState (0);

    updaterThreadStart.unlock ();
}

void ImProcCoordinator::fullUpdateDetailedCrops () { 

    if (destroying)
        return;

    updaterThreadStart.lock ();
    if (updaterRunning && thread) {
        changeSinceLast = 0;
        thread->join ();
    }

    if (plistener)
        plistener->setProgressState (1);

    for (int i=0; i<crops.size(); i++)
        crops[i]->update (ALL, true); 

    if (plistener)
        plistener->setProgressState (0);

    updaterThreadStart.unlock ();
}


void ImProcCoordinator::saveInputICCReference (const Glib::ustring& fname) {

    mProcessing.lock ();

    int fW, fH;
    imgsrc->getFullSize (fW, fH, 0);
    PreviewProps pp (0, 0, fW, fH, 1);
    ProcParams ppar = params;
    ppar.hlrecovery.enabled = false;
    ppar.icm.input = "(none)";
    Image16* im = new Image16 (fW, fH);
    imgsrc->preprocess( ppar.raw );
    imgsrc->demosaic(ppar.raw );
    imgsrc->getImage (imgsrc->getWB(), 0, im, pp, ppar.hlrecovery, ppar.icm, ppar.raw);
    im->saveJPEG (fname, 85);
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

            batchThread->yield(); //the running batch should wait other threads to avoid conflict
            
            thread = Glib::Thread::create(sigc::mem_fun(*this, &ImProcCoordinator::process), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);

        }
        else
            updaterThreadStart.unlock ();
    }
}

void ImProcCoordinator::process () {

    if (plistener)
        plistener->setProgressState (1);

    paramsUpdateMutex.lock ();
    while (changeSinceLast) {
        params = nextParams;
        int ch = changeSinceLast;
        changeSinceLast = 0;
        paramsUpdateMutex.unlock ();
        if (ch&32767)
            updatePreviewImage (ch);
        paramsUpdateMutex.lock ();
    }    
    paramsUpdateMutex.unlock ();
    updaterRunning = false;

    if (plistener)
        plistener->setProgressState (0);
}

ProcParams* ImProcCoordinator::getParamsForUpdate (ProcEvent change) {

    paramsUpdateMutex.lock ();
    changeSinceLast |= refreshmap[(int)change];
    return &nextParams;
}

void ImProcCoordinator::paramsUpdateReady () {

    paramsUpdateMutex.unlock ();
    startProcessing ();
}


}
