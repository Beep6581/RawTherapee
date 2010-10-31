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
#include <rtengine.h>
#include <colortemp.h>
#include <imagesource.h>
#include <improcfun.h>
#include <curves.h>
#include <iccstore.h>
#include <processingjob.h>
#include <glibmm.h>

#include <iostream>

#undef THREAD_PRIORITY_NORMAL

namespace rtengine {

IImage16* processImage (ProcessingJob* pjob, int& errorCode, ProgressListener* pl) {

    errorCode = 0;

    ProcessingJobImpl* job = (ProcessingJobImpl*) pjob;

    if (pl) {
        pl->setProgressStr ("Processing...");
        pl->setProgress (0.0);
    }
    
    InitialImage* ii = job->initialImage;
    if (!ii) {
        ii = InitialImage::load (job->fname, job->isRaw, &errorCode);
        if (errorCode) {
            ii->decreaseRef ();
            delete job;
            return NULL;
        }
    }
    procparams::ProcParams& params = job->pparams;
    
    // aquire image from imagesource
    ImageSource* imgsrc = ii->getImageSource ();
    ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green);
    if (params.wb.method=="Camera")
        currWB = imgsrc->getWB ();
    else if (params.wb.method=="Auto")
        currWB = imgsrc->getAutoWB ();

    int tr = TR_NONE;
    if (params.coarse.rotate==90)  tr |= TR_R90;
    if (params.coarse.rotate==180) tr |= TR_R180;
    if (params.coarse.rotate==270) tr |= TR_R270;
    if (params.coarse.hflip)       tr |= TR_HFLIP;
    if (params.coarse.vflip)       tr |= TR_VFLIP;

    int fw, fh;
    imgsrc->getFullSize (fw, fh, tr);

    ImProcFunctions ipf (&params, true);
    
    Image16* baseImg;
    PreviewProps pp (0, 0, fw, fh, 1);
    if (fabs(params.resize.scale-1.0)<1e-5) {
        baseImg = new Image16 (fw, fh);
        imgsrc->getImage (currWB, tr, baseImg, pp, params.hlrecovery, params.icm);
    }
    else {
        Image16* oorig = new Image16 (fw, fh);
        imgsrc->getImage (currWB, tr, oorig, pp, params.hlrecovery, params.icm);
        fw *= params.resize.scale;
        fh *= params.resize.scale;
        baseImg = new Image16 (fw, fh);
        ipf.resize (oorig, baseImg);
        delete oorig;
    }
    if (pl) 
        pl->setProgress (0.25);

    // perform first analysis
    unsigned int* hist16 = new unsigned int[65536];
    ipf.firstAnalysis (baseImg, &params, hist16, imgsrc->getGamma());

    // perform transform
    if (ipf.needsTransform()) {
        Image16* trImg = new Image16 (fw, fh);
        ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh);
        delete baseImg;
        baseImg = trImg;
    }

    // update blurmap
    int** buffer = new int*[fh];
    for (int i=0; i<fh; i++)
        buffer[i] = new int[fw];

    SHMap* shmap = NULL;
    if (params.sh.enabled) {
        shmap = new SHMap (fw, fh, true);
        double radius = sqrt (double(fw*fw+fh*fh)) / 2.0;
        double shradius = radius / 1800.0 * params.sh.radius;
        shmap->update (baseImg, (unsigned short**)buffer, shradius, ipf.lumimul, params.sh.hq);
    }
    // RGB processing
//!!!// auto exposure!!!
    double br = params.toneCurve.expcomp;
    int    bl = params.toneCurve.black;

    if (params.toneCurve.autoexp) {
        unsigned int aehist[65536]; int aehistcompr;
        imgsrc->getAEHistogram (aehist, aehistcompr);
        ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, br, bl);
    }

    int* curve1 = new int [65536];
    int* curve2 = new int [65536];

    CurveFactory::complexCurve (br, bl/65535.0, params.toneCurve.hlcompr, params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast, imgsrc->getDefGain(), imgsrc->getGamma(), true, params.toneCurve.curve, hist16, curve1, curve2, NULL);

    LabImage* labView = new LabImage (baseImg);
    ipf.rgbProc (baseImg, labView, curve1, curve2, shmap);

    if (shmap)
        delete shmap;

    if (pl) 
        pl->setProgress (0.5);


    // luminance histogram update
    memset (hist16, 0, 65536*sizeof(int));
    for (int i=0; i<fh; i++)
        for (int j=0; j<fw; j++)
            hist16[labView->L[i][j]]++;

    // luminance processing
    CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, params.labCurve.brightness, params.labCurve.contrast, 0.0, 0.0, false, params.labCurve.lcurve, hist16, curve1, curve2, NULL);
    ipf.luminanceCurve (labView, labView, curve2, 0, fh);
	CurveFactory::complexsgnCurve (0.0, 100.0, params.labCurve.saturation, 1.0, params.labCurve.acurve, curve1, 1);
    ipf.chrominanceCurve (labView, labView, 0, curve1, 0, fh);    
	CurveFactory::complexsgnCurve (0.0, 100.0, params.labCurve.saturation, 1.0, params.labCurve.bcurve, curve1, 1);
    ipf.chrominanceCurve (labView, labView, 1, curve1, 0, fh);
	
  	ipf.impulsedenoise (labView);
	ipf.lumadenoise (labView, buffer);
    ipf.sharpening (labView, (unsigned short**)buffer);

    delete [] curve1;
	delete [] curve2;
    delete [] hist16;

    // color processing
    ipf.colorCurve (labView, labView);
    ipf.colordenoise (labView, buffer);
	ipf.dirpyrdenoise (labView);

    // wavelet equalizer
    ipf.waveletEqualizer (labView, true, true);
	
	// directional pyramid equalizer
    ipf.dirpyrequalizer (labView);

    for (int i=0; i<fh; i++)
        delete [] buffer[i];
    delete [] buffer;

    if (pl) 
        pl->setProgress (0.75);

    // obtain final image
    Image16* readyImg;
    int cx = 0, cy = 0, cw = labView->W, ch = labView->H;
    if (params.crop.enabled) {
        cx = params.crop.x;
        cy = params.crop.y;
        cw = params.crop.w;
        ch = params.crop.h;
    }
    readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm.output);

    if (pl) 
        pl->setProgress (1.0);

    readyImg->setMetadata (ii->getMetaData()->getExifData (), params.exif, params.iptc);

    ProfileContent pc;
    if (params.icm.output.compare (0, 6, "No ICM") && params.icm.output!="")  
        pc = iccStore->getContent (params.icm.output);

    readyImg->setOutputProfile (pc.data, pc.length);

    delete baseImg;
    
    if (!job->initialImage)
        ii->decreaseRef ();
    
    delete job;

    if (pl) {
        pl->setProgress (1.0);
        pl->setProgressStr ("Ready.");
    }

    return readyImg;
}

void batchProcessingThread (ProcessingJob* job, BatchProcessingListener* bpl) {

    ProcessingJob* currentJob = job;
    
    while (currentJob) {
        int errorCode;
        IImage16* img = processImage (currentJob, errorCode, bpl);
        if (errorCode) 
            bpl->error ("Can not load input image.");
        currentJob = bpl->imageReady (img);
    }
}

void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl) {

    if (bpl)
        Glib::Thread::create(sigc::bind(sigc::ptr_fun(batchProcessingThread), job, bpl), 0, true, true, Glib::THREAD_PRIORITY_LOW);
    
}

}
