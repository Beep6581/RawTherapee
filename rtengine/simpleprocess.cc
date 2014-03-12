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
#include "rtengine.h"
#include "colortemp.h"
#include "imagesource.h"
#include "improcfun.h"
#include "curves.h"
#include "iccstore.h"
#include "processingjob.h"
#include <glibmm.h>
#include "../rtgui/options.h"
#include <iostream>
#include "rawimagesource.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/multilangmgr.h"
//#include "mytime.h"

#undef THREAD_PRIORITY_NORMAL
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {
extern const Settings* settings;

IImage16* processImage (ProcessingJob* pjob, int& errorCode, ProgressListener* pl, bool tunnelMetaData) {

    errorCode = 0;

    ProcessingJobImpl* job = static_cast<ProcessingJobImpl*>(pjob);

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

    // acquire image from imagesource
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
//    MyTime t1,t2;
//    t1.set();

    ImProcFunctions ipf (&params, true);

    PreviewProps pp (0, 0, fw, fh, 1);
    imgsrc->preprocess( params.raw, params.lensProf, params.coarse);

    if (params.toneCurve.autoexp) {// this enabled HLRecovery
        LUTu histRedRaw(256), histGreenRaw(256), histBlueRaw(256);
        imgsrc->getRAWHistogram(histRedRaw, histGreenRaw, histBlueRaw);
        if (ToneCurveParams::HLReconstructionNecessary(histRedRaw, histGreenRaw, histBlueRaw) && !params.toneCurve.hrenabled) {
            params.toneCurve.hrenabled=true;
            // WARNING: Highlight Reconstruction is being forced 'on', should we force a method here too?
        }
    }

    if (pl) pl->setProgress (0.20);
    imgsrc->demosaic( params.raw);
    if (pl) pl->setProgress (0.30);
    imgsrc->HLRecovery_Global( params.toneCurve );
    if (pl) pl->setProgress (0.40);
	// set the color temperature
    ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);
    if (params.wb.method=="Camera")
        currWB = imgsrc->getWB ();
    else if (params.wb.method=="Auto") {
        double rm, gm, bm;
        imgsrc->getAutoWBMultipliers(rm, gm, bm);
        currWB.update(rm, gm, bm, params.wb.equal);
    }
    Imagefloat* baseImg = new Imagefloat (fw, fh);
    imgsrc->getImage (currWB, tr, baseImg, pp, params.toneCurve, params.icm, params.raw);
    if (pl) pl->setProgress (0.45);


    // perform luma/chroma denoise
//	CieImage *cieView;    

    if (params.dirpyrDenoise.enabled) {
		ipf.RGB_denoise(baseImg, baseImg, imgsrc->isRAW(), params.dirpyrDenoise, params.defringe, imgsrc->getDirPyrDenoiseExpComp());
    }
    imgsrc->convertColorSpace(baseImg, params.icm, currWB, params.raw);

    // perform first analysis
    LUTu hist16 (65536);
    LUTu hist16C (65536);
	
    ipf.firstAnalysis (baseImg, &params, hist16, imgsrc->getGamma());

    // perform transform (excepted resizing)
    if (ipf.needsTransform()) {
        Imagefloat* trImg = new Imagefloat (fw, fh);
        ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh, fw, fh, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(),
                       imgsrc->getMetaData()->getFocusDist(), imgsrc->getRotateDegree(), true);
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
		LUTu aehist; int aehistcompr;
		imgsrc->getAutoExpHistogram (aehist, aehistcompr);
		ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, expcomp, bright, contr, black, hlcompr,hlcomprthresh);
    }

    // at this stage, we can flush the raw data to free up quite an important amount of memory
    // commented out because it makes the application crash when batch processing...
    // TODO: find a better place to flush rawData and rawRGB
    //imgsrc->flushRawData();
    //imgsrc->flushRGB();

    LUTf curve1 (65536,0);
    LUTf curve2 (65536,0);
	LUTf curve (65536,0);
	LUTf satcurve (65536,0);
	LUTf lhskcurve (65536,0);
	LUTf lumacurve(65536,0);
	LUTf clcurve (65536,0);

	LUTf rCurve (65536,0);
	LUTf gCurve (65536,0);
	LUTf bCurve (65536,0);
	LUTu dummy;

	ToneCurve customToneCurve1, customToneCurve2;
	ColorAppearance customColCurve1, customColCurve2,customColCurve3 ;
	ToneCurve customToneCurvebw1;
	ToneCurve customToneCurvebw2;
	//if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;

	ipf.g = imgsrc->getGamma();
	ipf.iGamma = true;
	CurveFactory::complexCurve (expcomp, black/65535.0, hlcompr, hlcomprthresh, params.toneCurve.shcompr, bright, contr, ipf.g, !ipf.iGamma,
	                            params.toneCurve.curveMode, params.toneCurve.curve, params.toneCurve.curveMode2, params.toneCurve.curve2,
	                            hist16, dummy, curve1, curve2, curve, dummy, customToneCurve1, customToneCurve2 );

	CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, 1);
	CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, 1);
	CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, 1);

    LabImage* labView = new LabImage (fw,fh);


	CurveFactory::curveBW (params.blackwhite.beforeCurve, params.blackwhite.afterCurve, hist16, dummy, customToneCurvebw1, customToneCurvebw2, 1);
	double rrm, ggm, bbm;
    float autor, autog, autob;
    autor = -9000.f; // This will ask to compute the "auto" values for the B&W tool (have to be inferior to -5000)
    ipf.rgbProc (baseImg, labView, curve1, curve2, curve, shmap, params.toneCurve.saturation, rCurve, gCurve, bCurve, customToneCurve1, customToneCurve2,customToneCurvebw1, customToneCurvebw2, rrm, ggm, bbm, autor, autog, autob, expcomp, hlcompr, hlcomprthresh);
    if (settings->verbose)
        printf("Output image / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", autor, autog, autob);

    // Freeing baseImg because not used anymore
    delete baseImg;
    baseImg = NULL;

    if (shmap)
    delete shmap;
    shmap = NULL;

    if (pl) 
        pl->setProgress (0.5);

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// start tile processing...???

    hist16.clear();  hist16C.clear();
	if(params.labCurve.contrast !=0) {//only use hist16 for contrast
	
#ifdef _OPENMP
#pragma omp parallel shared(hist16,labView, fh, fw)
#endif
{
    LUTu hist16thr (65536);  // one temporary lookup table per thread
    hist16thr.clear();
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
    for (int i=0; i<fh; i++)
        for (int j=0; j<fw; j++){
            hist16thr[CLIP((int)((labView->L[i][j])))]++;
        }

#pragma omp critical
{
    for(int i=0;i<65536;i++)
        hist16[i] += hist16thr[i];
}
}
	
	
	}
	bool utili=false;
	bool autili=false;
	bool butili=false;
	bool ccutili=false;
	bool cclutili=false;	
	bool clcutili=false;	

	CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve,hist16, hist16, lumacurve, dummy, 1, utili);
	CurveFactory::curveCL(clcutili, params.labCurve.clcurve, clcurve, hist16C, dummy, 1);

	CurveFactory::complexsgnCurve (autili, butili, ccutili, cclutili, params.labCurve.chromaticity, params.labCurve.rstprotection,
								   params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,params.labCurve.lccurve,curve1, curve2, satcurve,lhskcurve, 
								   hist16C, hist16C, dummy,dummy,
								   1);

	ipf.chromiLuminanceCurve (1,labView, labView, curve1, curve2, satcurve,lhskcurve,clcurve, lumacurve, utili, autili, butili, ccutili,cclutili, clcutili, dummy, dummy, dummy, dummy);
	
 	if((params.colorappearance.enabled && !params.colorappearance.tonecie) || (!params.colorappearance.enabled))ipf.EPDToneMap(labView,5,1);
	

	ipf.vibrance(labView);

	if((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) ipf.impulsedenoise (labView);
	// for all treatments Defringe, Sharpening, Contrast detail ,Microcontrast they are activated if "CIECAM" function are disabled

	if((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) ipf.defringe (labView);
	
	if (params.sharpenEdge.enabled) {
		 ipf.MLsharpen(labView);
	}
	if (params.sharpenMicro.enabled) {
		if((params.colorappearance.enabled && !settings->autocielab) ||  (!params.colorappearance.enabled)) ipf.MLmicrocontrast (labView);//!params.colorappearance.sharpcie
	}
	
	if(((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) && params.sharpening.enabled) {			
                    
        float **buffer = new float*[fh];
            for (int i=0; i<fh; i++)
                buffer[i] = new float[fw];

        ipf.sharpening (labView, (float**)buffer);

            for (int i=0; i<fh; i++)
                delete [] buffer[i];
                delete [] buffer;
    }
	// directional pyramid equalizer
	if((params.colorappearance.enabled && !settings->autocielab)  || !params.colorappearance.enabled) ipf.dirpyrequalizer (labView, 1);//TODO: this is the luminance tonecurve, not the RGB one
	
	//Colorappearance and tone-mapping associated
	
	int f_w=1,f_h=1;
	int begh = 0, endh = fh;
	if(params.colorappearance.tonecie || params.colorappearance.enabled){f_w=fw;f_h=fh;}
	CieImage *cieView = new CieImage (f_w,(f_h));
	begh=0;
	endh=fh;
	CurveFactory::curveLightBrightColor (
					params.colorappearance.curveMode, params.colorappearance.curve,
					params.colorappearance.curveMode2, params.colorappearance.curve2,
					params.colorappearance.curveMode3, params.colorappearance.curve3,
					hist16, hist16,dummy,
					hist16C, dummy,
					customColCurve1,
					customColCurve2,
					customColCurve3,
					1);
	float adap2,adap;
	double ada, ada2;
	if(params.colorappearance.enabled){
		float fnum = imgsrc->getMetaData()->getFNumber  ();// F number
		float fiso = imgsrc->getMetaData()->getISOSpeed () ;// ISO
		float fspeed = imgsrc->getMetaData()->getShutterSpeed () ;//speed
		float fcomp = imgsrc->getMetaData()->getExpComp  ();//compensation + -
		if(fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) {adap=adap2=2000.f;ada=ada2=2000.;}//if no exif data or wrong
		else {
			float E_V = fcomp + log2 ((fnum*fnum) / fspeed / (fiso/100.f));
			float expo2= params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
			E_V += expo2;
			float expo1;//exposure raw white point
			expo1=log2(params.raw.expos);//log2 ==>linear to EV
			E_V += expo1;
			adap2 = adap= powf(2.f, E_V-3.f);//cd / m2
			ada=ada2=(double) adap;
		}
		LUTf CAMBrightCurveJ;
		LUTf CAMBrightCurveQ;
		float CAMMean;
		if (params.sharpening.enabled) {
			float d;
			double dd;

			float** buffer = new float*[fh];
			for (int i=0; i<fh; i++)
			buffer[i] = new float[fw];
			int sk=1;
			if(settings->ciecamfloat) ipf.ciecam_02float (cieView, adap, begh, endh,1,2, labView, &params,customColCurve1,customColCurve2,customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, (float**)buffer, true, d, sk, 1);
			else ipf.ciecam_02 (cieView, ada, begh, endh,1,2, labView, &params,customColCurve1,customColCurve2,customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, (float**)buffer, true, dd, 1, 1);
			for (int i=0; i<fh; i++)
			delete [] buffer[i];
			delete [] buffer;
		}
		else {
			int f_h=2,f_w=2;
			float d;

			double dd;
			float** buffer = new float*[f_h];
			for (int i=0; i<f_h; i++)
			buffer[i] = new float[f_w];
			int sk=1;
			if(settings->ciecamfloat) ipf.ciecam_02float (cieView, adap, begh, endh,1,2, labView, &params,customColCurve1,customColCurve2,customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, (float**)buffer, true, d, sk, 1);
			else ipf.ciecam_02 (cieView, adap, begh, endh,1, 2, labView, &params,customColCurve1,customColCurve2,customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, (float**)buffer, true, dd, 1, 1);

			for (int i=0; i<f_h; i++)
			delete [] buffer[i];
			delete [] buffer;
		}
	}	
    delete cieView;
    cieView = NULL;
	
	
	

	// end tile processing...???
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        cmsMLU *DescriptionMLU, *CopyrightMLU, *DmndMLU, *DmddMLU;// for modification TAG

        cmsToneCurve* GammaTRC[3] = { NULL, NULL, NULL };
        cmsFloat64Number Parameters[7];
        double ga0,ga1,ga2,ga3,ga4,ga5,ga6;
       // wchar_t string[80] ;
        int ns;//numero of stri[]
        if      (params.icm.working=="ProPhoto")   ns=0;
        else if (params.icm.working=="Adobe RGB")  ns=1;
        else if (params.icm.working=="sRGB")  	   ns=2;
        else if (params.icm.working=="WideGamut")  ns=3;
        else if (params.icm.working=="Beta RGB")   ns=4;
        else if (params.icm.working=="BestRGB")    ns=5;
        else if (params.icm.working=="BruceRGB")   ns=6;
	//	if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;
        readyImg = ipf.lab2rgb16b (labView, cx, cy, cw, ch, params.icm.output, params.icm.working, params.icm.gamma, params.icm.freegamma, params.icm.gampos, params.icm.slpos, ga0,ga1,ga2,ga3,ga4,ga5,ga6, params.blackwhite.enabled );
        customGamma = true;

        //or selected Free gamma
        useLCMS=false;
        bool pro=false;
        Glib::ustring chpro, outProfile;
        bool present_space[9]={false,false,false,false,false,false,false,false,false};
        std::vector<Glib::ustring> opnames = iccStore->getOutputProfiles ();
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
                Parameters[0] = ga0;
                Parameters[1] = ga1;
                Parameters[2] = ga2;
                Parameters[3] = ga3;
                Parameters[4] = ga4;
                Parameters[5] = ga5;
                Parameters[6] = ga6;
                // 7 parameters for smoother curves
				//change desc Tag , to "free gamma", or "BT709", etc.
				cmsContext ContextID = cmsGetProfileContextID(jprof);//modification TAG
				DescriptionMLU  = cmsMLUalloc(ContextID, 1);
				CopyrightMLU    = cmsMLUalloc(ContextID, 1);//for ICC
				DmndMLU=cmsMLUalloc(ContextID, 1);//for ICC
				DmddMLU=cmsMLUalloc(ContextID, 1);// for ICC


				// instruction with //ICC are used for generate icc profile
				if (DescriptionMLU == NULL) printf("Description error\n");
				cmsMLUsetWide(CopyrightMLU, "en", "US", L"General Public License - AdobeRGB compatible")	;//adapt to profil
				cmsMLUsetWide(DmndMLU,      "en", "US", L"RawTherapee")	;
				cmsMLUsetWide(DmddMLU,      "en", "US", L"RTMedium")	;	//adapt to profil
				//display Tag desc with : selection of gamma and Primaries
				if (!params.icm.freegamma) {
					std::wstring gammaStr;
					if(params.icm.gamma=="High_g1.3_s3.35") {
						gammaStr = std::wstring(L"GammaTRC: High g=1.3 s=3.35");
					}
					else if (params.icm.gamma=="Low_g2.6_s6.9") {
						gammaStr = std::wstring(L"GammaTRC: Low g=2.6 s=6.9");
					}
					else if (params.icm.gamma=="sRGB_g2.4_s12.92") {
						gammaStr = std::wstring(L"GammaTRC: sRGB g=2.4 s=12.92");
					}
					else if (params.icm.gamma== "BT709_g2.2_s4.5") {
						gammaStr = std::wstring(L"GammaTRC: BT709 g=2.2 s=4.5");
					}
					else if (params.icm.gamma== "linear_g1.0") {
						gammaStr = std::wstring(L"GammaTRC: Linear g=1.0");
					}
					else if (params.icm.gamma== "standard_g2.2") {
						gammaStr = std::wstring(L"GammaTRC: g=2.2");
					}
					else if (params.icm.gamma== "standard_g1.8") {
						gammaStr = std::wstring(L"GammaTRC: g=1.8");
					}
					cmsMLUsetWide(DescriptionMLU,  "en", "US", gammaStr.c_str());

					//for elaboration ICC profiles
				//	else if (params.icm.gamma==	"sRGB_g2.4_s12.92" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Medium gamma sRGB(AdobeRGB compatible)");
				//	else if (params.icm.gamma==	"BT709_g2.2_s4.5" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_sRGB gamma BT709(IEC61966 equivalent)");
				//	else if (params.icm.gamma==	"sRGB_g2.4_s12.92" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_sRGB gamma sRGB(IEC61966 equivalent)");
				//	else if (params.icm.gamma==	"linear_g1.0" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_sRGB gamma Linear1.0(IEC61966 equivalent)");
				//else if (params.icm.gamma==	"BT709_g2.2_s4.5" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Large gamma BT709(Prophoto compatible)");
				//	else if (params.icm.gamma==	"sRGB_g2.4_s12.92" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Large gamma sRGB(Prophoto compatible)");
				//	else if (params.icm.gamma==	"linear_g1.0" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Large gamma Linear1.0(Prophoto compatible)");
				}
				else {
					// create description with gamma + slope + primaries
					std::wostringstream gammaWs;
					gammaWs.precision(2);
					gammaWs<<"Manual GammaTRC: g="<<(float)params.icm.gampos<<" s="<<(float)params.icm.slpos;
					cmsMLUsetWide(DescriptionMLU,  "en", "US", gammaWs.str().c_str());
				}

				cmsWriteTag(jprof, cmsSigProfileDescriptionTag,  DescriptionMLU);//desc changed
			//	cmsWriteTag(jprof, cmsSigCopyrightTag,           CopyrightMLU);    
			//	cmsWriteTag(jprof, cmsSigDeviceMfgDescTag, DmndMLU);
			//	cmsWriteTag(jprof, cmsSigDeviceModelDescTag, DmddMLU);
				
                // Calculate output profile's rTRC bTRC gTRC
                GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(NULL, 5, Parameters);
                cmsWriteTag(jprof, cmsSigGreenTRCTag, (void*)GammaTRC[1] );
                cmsWriteTag(jprof, cmsSigRedTRCTag,   (void*)GammaTRC[0] );
                cmsWriteTag(jprof, cmsSigBlueTRCTag,  (void*)GammaTRC[2] );
				//for generation ICC profiles : here Prophoto ==> Large
			//	if(params.icm.gamma==	"BT709_g2.2_s4.5") cmsSaveProfileToFile(jprof, "RT_sRGB_gBT709.icm");
			//	else if (params.icm.gamma==	"sRGB_g2.4_s12.92") cmsSaveProfileToFile(jprof, "RT_Medium_gsRGB.icc");
			//	else if (params.icm.gamma==	"linear_g1.0") cmsSaveProfileToFile(jprof, "RT_Large_g10.icc");
			
				
            }
        }
        if (GammaTRC[0]) cmsFreeToneCurve(GammaTRC[0]);
    }
    else {
        // if Default gamma mode: we use the profile selected in the "Output profile" combobox;
        // gamma come from the selected profile, otherwise it comes from "Free gamma" tool
		
        readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm.output, params.blackwhite.enabled);
        if (settings->verbose) printf("Output profile: \"%s\"\n", params.icm.output.c_str());
    }

    delete labView;
    labView = NULL;
	
	if(params.blackwhite.enabled) {//force BW r=g=b
		for (int ccw=0;ccw<cw;ccw++) {
			for (int cch=0;cch<ch;cch++) {
			readyImg->r(cch,ccw)=readyImg->g(cch,ccw);
			readyImg->b(cch,ccw)=readyImg->g(cch,ccw);
			}
		}
	}
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
            // use corrected sRGB profile in order to apply a good TRC if present, otherwise use LCMS2 profile generated by lab2rgb16b
            ProfileContent pc(jprof);
            readyImg->setOutputProfile (pc.data, pc.length);
        }
    }
    else {
        // use RT_sRGB.icm profile if present, otherwise use LCMS2 profile generate by lab2rgb16b
        Glib::ustring outputProfile;
        if (params.icm.output!="" && params.icm.output!=ColorManagementParams::NoICMString) {
            outputProfile = params.icm.output;

            /*  if we'd wanted the RT_sRGB profile we would have selected it
        else {
            // use RT_sRGB.icm profile if present, otherwise use LCMS2 profile generate by lab2rgb16b
            if (settings->verbose) printf("No output profiles set ; looking for the default sRGB profile (\"%s\")...\n", options.rtSettings.srgb.c_str());
            outputProfile = options.rtSettings.srgb;
            }*/

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
        } else {
            // No ICM
             readyImg->setOutputProfile (NULL,0);
        }
    }
//    t2.set();
//    if( settings->verbose )
//           printf("Total:- %d usec\n", t2.etime(t1));

    if (!job->initialImage)
        ii->decreaseRef ();

    delete job;
    if (pl)
        pl->setProgress (0.75);
/*	curve1.reset();curve2.reset();
	curve.reset();
	satcurve.reset();
	lhskcurve.reset();
	
	rCurve.reset();
	gCurve.reset();
	bCurve.reset();
	hist16.reset();
	hist16C.reset();
*/
    return readyImg;
}

void batchProcessingThread (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData) {

    ProcessingJob* currentJob = job;
    
    while (currentJob) {
        int errorCode;
        IImage16* img = processImage (currentJob, errorCode, bpl, tunnelMetaData);
        if (errorCode) {
            bpl->error (M("MAIN_MSG_CANNOTLOAD"));
            currentJob = NULL;
        } else {
            try {
                currentJob = bpl->imageReady (img);
            } catch (Glib::Exception& ex) {
                bpl->error (ex.what());
                currentJob = NULL;
            }
        }
    }
}

void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData) {

    if (bpl)
        Glib::Thread::create(sigc::bind(sigc::ptr_fun(batchProcessingThread), job, bpl, tunnelMetaData), 0, true, true, Glib::THREAD_PRIORITY_LOW);
    
}

}
