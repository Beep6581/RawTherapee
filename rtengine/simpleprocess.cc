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
#include "ppversion.h"
#undef THREAD_PRIORITY_NORMAL

#define CLIP(a) ((a)>0?((a)<65535?(a):65535):0)


namespace rtengine {
extern const Settings* settings;

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
    imgsrc->preprocess( params.raw);
	if (pl) pl->setProgress (0.20);
    imgsrc->demosaic( params.raw);
    if (pl) pl->setProgress (0.30);
    imgsrc->HLRecovery_Global( params.hlrecovery );
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
    double expcomp = params.toneCurve.expcomp;
	int    bright = params.toneCurve.brightness;
	int	   contr = params.toneCurve.contrast;
	int    black = params.toneCurve.black;
	int	   hlcompr = params.toneCurve.hlcompr;
	int    hlcomprthresh = params.toneCurve.hlcomprthresh;

    if (params.toneCurve.autoexp) {
		printf("silpleprocess calling autoexp\n");
		LUTu aehist; int aehistcompr;
		imgsrc->getAutoExpHistogram (aehist, aehistcompr);
		ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, expcomp, bright, contr, black, hlcompr,hlcomprthresh);
    }

    // at this stage, we can flush the raw data to free up quite an important amount of memory
    // commented out because it make the application crash when batch processing...
    // TODO: find a better place to flush rawData and rawRGB
    //imgsrc->flushRawData();
    //imgsrc->flushRGB();

    LUTf curve1 (65536,0);
    LUTf curve2 (65536,0);
	LUTf curve (65536,0);
	LUTf satcurve (65536,0);
	LUTf rCurve (65536,0);
	LUTf gCurve (65536,0);
	LUTf bCurve (65536,0);
	LUTu dummy;

    CurveFactory::complexCurve (expcomp, black/65535.0, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh, params.toneCurve.shcompr, bright, params.toneCurve.contrast, imgsrc->getGamma(), true, params.toneCurve.curve, 
        hist16, dummy, curve1, curve2, curve, dummy);
	
	CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, 1);
	CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, 1);
	CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, 1);

	LabImage* labView = new LabImage (fw,fh);

    ipf.rgbProc (baseImg, labView, curve1, curve2, curve, shmap, params.toneCurve.saturation, rCurve, gCurve, bCurve);

    // Freeing baseImg because not used anymore
    delete baseImg;
    baseImg = NULL;

    if (shmap)
    delete shmap;
    shmap = NULL;

    if (pl) 
        pl->setProgress (0.5);


    // luminance histogram update
    hist16.clear();
    for (int i=0; i<fh; i++)
        for (int j=0; j<fw; j++)
            hist16[CLIP((int)((labView->L[i][j])))]++;

    // luminance processing

	ipf.EPDToneMap(labView);

	CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, hist16, hist16, curve, dummy, 1);

	CurveFactory::complexsgnCurve (params.labCurve.saturation, params.labCurve.enable_saturationlimiter, params.labCurve.saturationlimit, \
								   params.labCurve.acurve, params.labCurve.bcurve, curve1, curve2, satcurve, 1);
	ipf.luminanceCurve (labView, labView, curve);
	ipf.chrominanceCurve (labView, labView, curve1, curve2, satcurve);
	ipf.vibrance(labView);

	ipf.impulsedenoise (labView);
	ipf.defringe (labView);
	ipf.dirpyrdenoise (labView);
	if (params.sharpenEdge.enabled) {
		 ipf.MLsharpen(labView);
	}
	if (params.sharpenMicro.enabled) {
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

    Image16* readyImg = NULL;
    cmsHPROFILE jprof = NULL;
    bool customGamma = false;
    bool useLCMS;

    if(params.icm.gamma != "default" || params.icm.freegamma) { // if select gamma output between BT709, sRGB, linear, low, high, 2.2 , 1.8
        double ga0,ga1,ga2,ga3,ga4,ga5,ga6;
        readyImg = ipf.lab2rgb16b (labView, cx, cy, cw, ch, params.icm.output, params.icm.working, params.icm.gamma, params.icm.freegamma, params.icm.gampos, params.icm.slpos, ga0,ga1,ga2,ga3,ga4,ga5,ga6 );
        customGamma = true;

        //or selected Free gamma
        useLCMS=false;
        bool pro=false;
        Glib::ustring chpro, outProfile;
        bool present_space[9]={false,false,false,false,false,false,false,false,false};
        std::vector<std::string> opnames = iccStore->getOutputProfiles ();
        //test if files are in system
        for (int j=0; j<9; j++) {
            // one can modify "option" [Color Management] to adapt the profile's name if they are different for windows, MacOS, Linux ??
            // some of them are actually provided by RT, thanks to Jacques Desmis
            if     (j==0) chpro=options.rtSettings.prophoto;
            else if(j==1) chpro=options.rtSettings.adobe;
            else if(j==2) chpro=options.rtSettings.widegamut;
            else if(j==3) chpro=options.rtSettings.beta;
            else if(j==4) chpro=options.rtSettings.best;
            else if(j==5) chpro=options.rtSettings.bruce;
            else if(j==6) chpro=options.rtSettings.srgb;
            else if(j==7) chpro=options.rtSettings.srgb10;//gamma 1.0
            else if(j==8) chpro=options.rtSettings.prophoto10;//gamma 1.0

            for (unsigned int i=0; i<opnames.size(); i++) {
                if(chpro.compare(opnames[i]) ==0) present_space[j]=true;
            }
            if (!present_space[j] && settings->verbose) printf("Missing file: %s\n", chpro.c_str());
        }
        if (params.icm.freegamma && params.icm.gampos < 1.35) pro=true; //select profil with gammaTRC modified :
        else if (params.icm.gamma=="linear_g1.0" || (params.icm.gamma=="High_g1.3_s3.35")) pro=true;//pro=0  RT_sRGB || Prophoto

        // Check that output profiles exist, otherwise use LCMS2
        // Use the icc/icm profiles associated to possible working profiles, set in "options"
        if      (params.icm.working=="ProPhoto"  && present_space[0] && !pro)  outProfile=options.rtSettings.prophoto;
        else if (params.icm.working=="Adobe RGB" && present_space[1]        )  outProfile=options.rtSettings.adobe;
        else if (params.icm.working=="WideGamut" && present_space[2]        )  outProfile=options.rtSettings.widegamut;
        else if (params.icm.working=="Beta RGB"  && present_space[3]        )  outProfile=options.rtSettings.beta;
        else if (params.icm.working=="BestRGB"   && present_space[4]        )  outProfile=options.rtSettings.best;
        else if (params.icm.working=="BruceRGB"  && present_space[5]        )  outProfile=options.rtSettings.bruce;
        else if (params.icm.working=="sRGB"      && present_space[6] && !pro)  outProfile=options.rtSettings.srgb;
        else if (params.icm.working=="sRGB"      && present_space[7] &&  pro)  outProfile=options.rtSettings.srgb10;
        else if (params.icm.working=="ProPhoto"  && present_space[8] &&  pro)  outProfile=options.rtSettings.prophoto10;
        else {
            // Should not occurs
            if (settings->verbose) printf("\"%s\": unknown working profile! - use LCMS2 substitution\n", params.icm.working.c_str() );
            useLCMS=true;
        }

        //begin adaptation rTRC gTRC bTRC
        //"jprof" profile has the same characteristics than RGB values, but TRC are adapted... for applying profile
        if (!useLCMS) {
            if (settings->verbose) printf("Output Gamma - profile: \"%s\"\n", outProfile.c_str()  );  //c_str()
            jprof = iccStore->getProfile(outProfile); //get output profile
            if (jprof == NULL) {
                useLCMS = true;
                if (settings->verbose) printf("\"%s\" ICC output profile not found!\n", outProfile.c_str());
            }
            else {
                cmsToneCurve* GammaTRC[3];
                cmsFloat64Number Parameters[7];
                Parameters[0] = ga0;
                Parameters[1] = ga1;
                Parameters[2] = ga2;
                Parameters[3] = ga3;
                Parameters[4] = ga4;
                Parameters[5] = ga5;
                Parameters[6] = ga6;
                // 7 parameters for smoother curves
                // Calculate output profile's rTRC bTRC gTRC
                GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(NULL, 5, Parameters);
                cmsWriteTag(jprof, cmsSigGreenTRCTag, (void*)GammaTRC[1] );
                cmsWriteTag(jprof, cmsSigRedTRCTag,   (void*)GammaTRC[0] );
                cmsWriteTag(jprof, cmsSigBlueTRCTag,  (void*)GammaTRC[2] );
            }
        }
    }
    else {
        // if Default gamma mode: we use the profile selected in the "Output profile" combobox;
        // gamma come from the selected profile, otherwise it comes from "Free gamma" tool
        readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm.output);
        if (settings->verbose) printf("Output profile: \"%s\"\n", params.icm.output.c_str());
    }

    delete labView;
    labView = NULL;

    if (pl) pl->setProgress (0.70);

    if (params.resize.enabled) {

        // get the resize parameters
        int refw, refh;
        double tmpScale;
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
            if ((double)refw/(double)refh > (double)params.resize.width/(double)params.resize.height)
                tmpScale = (double)params.resize.width/(double)refw;
            else
                tmpScale = (double)params.resize.height/(double)refh;
        break;
        default:
            // Scale
            tmpScale = params.resize.scale;
            break;
        }

        // resize image
        if (fabs(tmpScale-1.0)>1e-5) {
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


    // Setting the output curve to readyImg
    if (customGamma) {
        if (!useLCMS) {
            // use corrected sRGB profile in order to apply a good TRC if present, otherwise use LCMS2 profile generate by lab2rgb16b
            ProfileContent pc(jprof);
            readyImg->setOutputProfile (pc.data, pc.length);
        }
    }
    else {
        // use RT_sRGB.icm profile if present, otherwise use LCMS2 profile generate by lab2rgb16b
        Glib::ustring outputProfile;
        if (params.icm.output!="" && params.icm.output!=ColorManagementParams::NoICMString) {
            outputProfile = params.icm.output;
        }
        else {
            // use RT_sRGB.icm profile if present, otherwise use LCMS2 profile generate by lab2rgb16b
            if (settings->verbose) printf("No output profiles set ; looking for the default sRGB profile (\"%s\")...\n", options.rtSettings.srgb.c_str());
            outputProfile = options.rtSettings.srgb;
        }

        // if iccStore->getProfile send back an object, then iccStore->getContent will do too
        cmsHPROFILE jprof = iccStore->getProfile(outputProfile); //get outProfile
        if (jprof == NULL) {
            if (settings->verbose) printf("\"%s\" ICC output profile not found!\n - use LCMS2 substitution\n", outputProfile.c_str());
        }
        else {
            if (settings->verbose) printf("Using \"%s\" output profile\n", outputProfile.c_str());
            ProfileContent pc = iccStore->getContent (outputProfile);
            readyImg->setOutputProfile (pc.data, pc.length);
        }
    }

    if (!job->initialImage)
        ii->decreaseRef ();

    delete job;
    if (pl)
        pl->setProgress (0.75);

    return readyImg;
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
