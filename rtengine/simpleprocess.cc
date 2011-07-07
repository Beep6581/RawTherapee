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
#include <options.h>
#include <iostream>
#include <rawimagesource.h>
#undef THREAD_PRIORITY_NORMAL

#define CLIP(a) ((a)>0?((a)<65535?(a):65535):0)


namespace rtengine {
IImage16* processImage (ProcessingJob* pjob, int& errorCode, ProgressListener* pl, bool tunnelMetaData) {

    errorCode = 0;

    ProcessingJobImpl* job = (ProcessingJobImpl*) pjob;

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_PROCESSING");
        pl->setProgress (0.0);
    }

    InitialImage* ii = job->initialImage;
    if (!ii) {
        ii = InitialImage::load (job->fname, job->isRaw, &errorCode);
        if (errorCode) {
            delete job;
            return NULL;
        }
    }
    procparams::ProcParams& params = job->pparams;

    // aquire image from imagesource
    ImageSource* imgsrc = ii->getImageSource ();

    int tr = TR_NONE;
    if (params.coarse.rotate==90)  tr |= TR_R90;
    if (params.coarse.rotate==180) tr |= TR_R180;
    if (params.coarse.rotate==270) tr |= TR_R270;
    if (params.coarse.hflip)       tr |= TR_HFLIP;
    if (params.coarse.vflip)       tr |= TR_VFLIP;

    int fw, fh;
    imgsrc->getFullSize (fw, fh, tr);

    // check the crop params
	if (params.crop.x > fw || params.crop.y > fh) {
		// the crop is completely out of the image, so we disable the crop
		params.crop.enabled = false;
		// and we set the values to the defaults
		params.crop.x = 0;
		params.crop.y = 0;
		params.crop.w = fw;
		params.crop.h = fh;
	}
	else {
		if ((params.crop.x + params.crop.w) > fw) {
			// crop overflow in the width dimension ; we trim it
			params.crop.w = fw-params.crop.x;
		}
		if ((params.crop.y + params.crop.h) > fh) {
			// crop overflow in the height dimension ; we trim it
			params.crop.h = fh-params.crop.y;
		}
	}

    ImProcFunctions ipf (&params, true);

	// set the color temperature
    ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green);
    if (params.wb.method=="Camera")
        currWB = imgsrc->getWB ();
    else if (params.wb.method=="Auto")
        currWB = imgsrc->getAutoWB ();

    PreviewProps pp (0, 0, fw, fh, 1);
    imgsrc->preprocess( params.raw, params.hlrecovery );
	if (pl) pl->setProgress (0.20);
    imgsrc->demosaic( params.raw, params.hlrecovery );
	if (pl) pl->setProgress (0.40);
    Imagefloat* baseImg = new Imagefloat (fw, fh);
    imgsrc->getImage (currWB, tr, baseImg, pp, params.hlrecovery, params.icm, params.raw);
    if (pl) pl->setProgress (0.45);


    // perform first analysis
    LUTu hist16 (65536);
    ipf.firstAnalysis (baseImg, &params, hist16, imgsrc->getGamma());

    // perform transform (excepted resizing)
    if (ipf.needsTransform()) {
        Imagefloat* trImg = new Imagefloat (fw, fh);
        ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh);
        delete baseImg;
        baseImg = trImg;
    }

    // update blurmap


    SHMap* shmap = NULL;
    if (params.sh.enabled) {
        shmap = new SHMap (fw, fh, true);
        double radius = sqrt (double(fw*fw+fh*fh)) / 2.0;
		double shradius = params.sh.radius;
		if (!params.sh.hq) shradius *= radius / 1800.0;
		shmap->update (baseImg, shradius, ipf.lumimul, params.sh.hq, 1);
    }
    // RGB processing
//!!!// auto exposure!!!
    double br = params.toneCurve.expcomp;
    int    bl = params.toneCurve.black;

    if (params.toneCurve.autoexp) {
        LUTu aehist; int aehistcompr;
        imgsrc->getAutoExpHistogram (aehist, aehistcompr);
        ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, br, bl);
    }

    LUTf curve1 (65536,0);
    LUTf curve2 (65536,0);
	LUTf curve (65536,0);
	LUTf satcurve (65536,0);
	LUTu dummy;

    CurveFactory::complexCurve (br, bl/65535.0, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast, imgsrc->getGamma(), true, params.toneCurve.curve, 
        hist16, dummy, curve1, curve2, curve, dummy);

	LabImage* labView = new LabImage (fw,fh);

    ipf.rgbProc (baseImg, labView, curve1, curve2, curve, shmap, params.toneCurve.saturation);

    if (shmap)
    delete shmap;

    if (pl) 
        pl->setProgress (0.5);


    // luminance histogram update
    hist16.clear();
    for (int i=0; i<fh; i++)
        for (int j=0; j<fw; j++)
            hist16[CLIP((int)((labView->L[i][j])))]++;

    // luminance processing
	CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, hist16, hist16, curve, dummy, 1);
	ipf.luminanceCurve (labView, labView, curve);
	CurveFactory::complexsgnCurve (params.labCurve.saturation, params.labCurve.enable_saturationlimiter, params.labCurve.saturationlimit, \
								   params.labCurve.acurve, params.labCurve.bcurve, curve1, curve2, satcurve, 1);
	ipf.chrominanceCurve (labView, labView, curve1, curve2, satcurve);

  	ipf.impulsedenoise (labView);
	ipf.defringe (labView);
	//ipf.lumadenoise (labView, buffer);
 	ipf.dirpyrdenoise (labView);
	if (params.clarity.enabled) {
		 ipf.MLsharpen(labView);
		 }
	if (params.clarity.enabledtwo) {		 
		 ipf.MLmicrocontrast (labView);
		 }
    if (params.sharpening.enabled) {
        float** buffer = new float*[fh];
        for (int i=0; i<fh; i++)
            buffer[i] = new float[fw];

        ipf.sharpening (labView, (float**)buffer);

        for (int i=0; i<fh; i++)
            delete [] buffer[i];
        delete [] buffer; buffer=NULL;
    }

    // wavelet equalizer
    //ipf.waveletEqualizer (labView, true, true);
	
	// directional pyramid equalizer
    ipf.dirpyrequalizer (labView);//TODO: this is the luminance tonecurve, not the RGB one


    if (pl) pl->setProgress (0.60);

    // crop and convert to rgb16
    int cx = 0, cy = 0, cw = labView->W, ch = labView->H;
    if (params.crop.enabled) {
        cx = params.crop.x;
        cy = params.crop.y;
        cw = params.crop.w;
        ch = params.crop.h;
    }
  if(params.icm.gamma != "default" || params.icm.freegamma)	
    {	// if select gamma output between BT709, sRGB, linear, low, high, 2.2 , 1.8
	//or selected Free gamma
	Image16* readyImg = ipf.lab2rgb16b (labView, cx, cy, cw, ch, params.icm.output, params.icm.working, params.icm.gamma, params.icm.freegamma, params.icm.gampos, params.icm.slpos);
	
    delete labView;
    if (pl) pl->setProgress (0.70);

    // get the resize parameters
	int refw, refh;
	double tmpScale;

	if (params.resize.enabled) {

		if (params.crop.enabled && params.resize.appliesTo == "Cropped area") {
			// the resize values applies to the crop dimensions
			refw = cw;
			refh = ch;
		}
		else {
			// the resize values applies to the image dimensions
			// if a crop exists, it will be resized to the calculated scale
			refw = fw;
			refh = fh;
		}

		switch(params.resize.dataspec) {
		case (1):
			// Width
			tmpScale = (double)params.resize.width/(double)refw;
			break;
		case (2):
			// Height
			tmpScale = (double)params.resize.height/(double)refh;
			break;
		case (3):
			// FitBox
			if ((double)refw/(double)refh > (double)params.resize.width/(double)params.resize.height) {
				tmpScale = (double)params.resize.width/(double)refw;
			}
			else {
				tmpScale = (double)params.resize.height/(double)refh;
			}
			break;
		default:
			// Scale
			tmpScale = params.resize.scale;
			break;
		}

	    // resize image
	    if (fabs(tmpScale-1.0)>1e-5)
	    {
	    	int imw, imh;
			if (params.crop.enabled && params.resize.appliesTo == "Full image") {
				imw = cw;
				imh = ch;
			}
			else {
				imw = refw;
				imh = refh;
			}
			imw = (int)( (double)imw * tmpScale + 0.5 );
			imh = (int)( (double)imh * tmpScale + 0.5 );
	        Image16* tempImage = new Image16 (imw, imh);
	        ipf.resize (readyImg, tempImage, tmpScale);
	        delete readyImg;
	        readyImg = tempImage;
	    }
	}


    if (tunnelMetaData)
        readyImg->setMetadata (ii->getMetaData()->getExifData ());
    else
        readyImg->setMetadata (ii->getMetaData()->getExifData (), params.exif, params.iptc);
	

    ProfileContent pc;
	Glib::ustring chpro;
	int present_space[7]={0,0,0,0,0,0,0};
	 std::vector<std::string> opnames = rtengine::getOutputProfiles ();
	 //test if files are in system
    for (int j=0; j<7;j++) {
	//one can modify "option" [Color Management] to adapt name of profile ih there are different for windows, MacOS, Linux ?? 
	if(j==0) chpro=options.rtSettings.prophoto;
	else if(j==1) chpro=options.rtSettings.adobe;
	else if(j==2) chpro=options.rtSettings.widegamut;	
	else if(j==3) chpro=options.rtSettings.beta;	
	else if(j==4) chpro=options.rtSettings.best;	
	else if(j==5) chpro=options.rtSettings.bruce;	
	else if(j==6) chpro=options.rtSettings.srgb;	
	for (int i=0; i<opnames.size(); i++)
       if(chpro.compare(opnames[i]) ==0) present_space[j]=1; 
	      if (present_space[j]==0) { if (pl) pl->setProgressStr ("Missing file..");pl->setProgress (0.0);}// display file not present: not very good display information...!!
        }
		//choose output profile
		if(params.icm.working=="ProPhoto" && present_space[0]==1)  params.icm.output=options.rtSettings.prophoto;
		else if(params.icm.working=="Adobe RGB" && present_space[1]==1)  params.icm.output=options.rtSettings.adobe;
		else if(params.icm.working=="WideGamut" && present_space[2]==1)  params.icm.output=options.rtSettings.widegamut;
		else if(params.icm.working=="Beta RGB" && present_space[3]==1)  params.icm.output=options.rtSettings.beta;
		else if(params.icm.working=="BestRGB" && present_space[4]==1)  params.icm.output=options.rtSettings.best;
		else if(params.icm.working=="BruceRGB" && present_space[5]==1)  params.icm.output=options.rtSettings.bruce;
		else params.icm.output=options.rtSettings.srgb; //if not found or choice=srgb

	
    if (params.icm.output.compare (0, 6, "No ICM") && params.icm.output!="")  
        pc = iccStore->getContent (params.icm.output);
	
    readyImg->setOutputProfile (pc.data, pc.length);
	
    delete baseImg;
    
    if (!job->initialImage)
        ii->decreaseRef ();

    delete job;
    if (pl)
        pl->setProgress (0.75);
	
    return readyImg;
	}
	else
	{//if default mode : profil = selection by choice in list (Prophoto.icm, sRGB.icm, etc., etc.) , gamma = gamma of profile or  not selected Free gamma
	
	Image16* readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm.output);
    delete labView;
    if (pl) pl->setProgress (0.70);

    // get the resize parameters
	int refw, refh;
	double tmpScale;

	if (params.resize.enabled) {

		if (params.crop.enabled && params.resize.appliesTo == "Cropped area") {
			// the resize values applies to the crop dimensions
			refw = cw;
			refh = ch;
		}
		else {
			// the resize values applies to the image dimensions
			// if a crop exists, it will be resized to the calculated scale
			refw = fw;
			refh = fh;
		}

		switch(params.resize.dataspec) {
		case (1):
			// Width
			tmpScale = (double)params.resize.width/(double)refw;
			break;
		case (2):
			// Height
			tmpScale = (double)params.resize.height/(double)refh;
			break;
		case (3):
			// FitBox
			if ((double)refw/(double)refh > (double)params.resize.width/(double)params.resize.height) {
				tmpScale = (double)params.resize.width/(double)refw;
			}
			else {
				tmpScale = (double)params.resize.height/(double)refh;
			}
			break;
		default:
			// Scale
			tmpScale = params.resize.scale;
			break;
		}

	    // resize image
	    if (fabs(tmpScale-1.0)>1e-5)
	    {
	    	int imw, imh;
			if (params.crop.enabled && params.resize.appliesTo == "Full image") {
				imw = cw;
				imh = ch;
			}
			else {
				imw = refw;
				imh = refh;
			}
			imw = (int)( (double)imw * tmpScale + 0.5 );
			imh = (int)( (double)imh * tmpScale + 0.5 );
	        Image16* tempImage = new Image16 (imw, imh);
	        ipf.resize (readyImg, tempImage, tmpScale);
	        delete readyImg;
	        readyImg = tempImage;
	    }
	}


    if (tunnelMetaData)
        readyImg->setMetadata (ii->getMetaData()->getExifData ());
    else
        readyImg->setMetadata (ii->getMetaData()->getExifData (), params.exif, params.iptc);
	

    ProfileContent pc;
    if (params.icm.output.compare (0, 6, "No ICM") && params.icm.output!="")  
        pc = iccStore->getContent (params.icm.output);
	
    readyImg->setOutputProfile (pc.data, pc.length);
	
    delete baseImg;
    
    if (!job->initialImage)
        ii->decreaseRef ();

    delete job;
    if (pl)
        pl->setProgress (0.75);
	
    return readyImg;
	
	}
}

void batchProcessingThread (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData) {

    ProcessingJob* currentJob = job;
    
    while (currentJob) {
        int errorCode;
        IImage16* img = processImage (currentJob, errorCode, bpl, tunnelMetaData);
        if (errorCode) 
            bpl->error ("Can not load input image.");
        currentJob = bpl->imageReady (img);
    }
}

void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData) {

    if (bpl)
        Glib::Thread::create(sigc::bind(sigc::ptr_fun(batchProcessingThread), job, bpl, tunnelMetaData), 0, true, true, Glib::THREAD_PRIORITY_LOW);
    
}

}
