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
#include <cmath>
#include <glib.h>
#include <glibmm.h>

#include "rtengine.h"
#include "improcfun.h"
#include "curves.h"
#include "colorclip.h"
#include "gauss.h"
#include "bilateral2.h"
#include "mytime.h"
#include "iccstore.h"
#include "impulse_denoise.h"
#include "imagesource.h"
#include "rtthumbnail.h"
#include "utils.h"
#include "iccmatrices.h"
#include "color.h"
#include "calc_distort.h"
#include "rt_math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {
	
	using namespace procparams;
	
#undef ABS
#undef CLIPS
#undef CLIPC

#define ABS(a) ((a)<0?-(a):(a))
#define CLIPS(a) ((a)>-32768?((a)<32767?(a):32767):-32768)
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define CLIP2(a) ((a)<MAXVAL ? a : MAXVAL )
#define FCLIP(a) ((a)>0.0?((a)<65535.5?(a):65535.5):0.0)
	
	
extern const Settings* settings;

LUTf ImProcFunctions::cachef ;
LUTf ImProcFunctions::gamma2curve = 0;

void ImProcFunctions::initCache () {

    int maxindex = 65536;
	cachef(maxindex,0/*LUT_CLIP_BELOW*/);

	gamma2curve(maxindex,0);

    for (int i=0; i<maxindex; i++) {
        if (i>Color::eps_max) {
			cachef[i] = 327.68*( exp(1.0/3.0 * log((double)i / MAXVAL) ));
        }
        else {
			cachef[i] = 327.68*((Color::kappa*i/MAXVAL+16.0)/116.0);
        }
	}

	for (int i=0; i<maxindex; i++) {
		gamma2curve[i] = (CurveFactory::gamma2(i/65535.0) * 65535.0);
	}
}

void ImProcFunctions::cleanupCache () {

}

ImProcFunctions::~ImProcFunctions () {

	if (monitorTransform!=NULL)
		cmsDeleteTransform (monitorTransform);
}

void ImProcFunctions::setScale (double iscale) {
	scale = iscale;
}

// Called from several threads
void ImProcFunctions::firstAnalysisThread (Imagefloat* original, Glib::ustring wprofile, unsigned int* histogram, int row_from, int row_to) {

	TMatrix wprof = iccStore->workingSpaceMatrix (wprofile);

	lumimul[0] = wprof[1][0];
	lumimul[1] = wprof[1][1];
	lumimul[2] = wprof[1][2];
	
    int W = original->width;
    for (int i=row_from; i<row_to; i++) {
        for (int j=0; j<W; j++) {
      
            int r = original->r[i][j];
            int g = original->g[i][j];
            int b = original->b[i][j];

            int y = CLIP((int)(lumimul[0] * r + lumimul[1] * g + lumimul[2] * b)) ;

            if (histogram) { 
                histogram[y]++;
            }
        }
    }
}

void ImProcFunctions::firstAnalysis (Imagefloat* original, const ProcParams* params, LUTu & histogram, double gamma) {

	// set up monitor transform
	Glib::ustring wprofile = params->icm.working;
	if (monitorTransform)
		cmsDeleteTransform (monitorTransform);
	monitorTransform = NULL;

    Glib::ustring monitorProfile=settings->monitorProfile;
    if (settings->autoMonitorProfile) monitorProfile=iccStore->defaultMonitorProfile;

	cmsHPROFILE monitor = iccStore->getProfile ("file:"+monitorProfile);
	if (monitor) {
        cmsHPROFILE iprof = iccStore->getXYZProfile ();       
        lcmsMutex->lock ();
		monitorTransform = cmsCreateTransform (iprof, TYPE_RGB_16, monitor, TYPE_RGB_8, settings->colorimetricIntent,
            cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is for thread safety, NOOPTIMIZE for precision
        lcmsMutex->unlock ();
	}
	
	//chroma_scale = 1;

	// calculate histogram of the y channel needed for contrast curve calculation in exposure adjustments

#ifdef _OPENMP
    int T = omp_get_max_threads();
#else
    int T = 1;
#endif

    unsigned int** hist = new unsigned int* [T];
    for (int i=0; i<T; i++) {
		hist[i] = new unsigned int[65536];
		memset (hist[i], 0, 65536*sizeof(int));
    }

#ifdef _OPENMP
	#pragma omp parallel if (multiThread)
    {
	        int H = original->height;
		int tid = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int blk = H/nthreads;

		if (tid<nthreads-1)
			firstAnalysisThread (original, wprofile, hist[tid], tid*blk, (tid+1)*blk);
		else
			firstAnalysisThread (original, wprofile, hist[tid], tid*blk, H);
    }
#else
    firstAnalysisThread (original, wprofile, hist[0], 0, original->height);
#endif
 
	histogram.clear();
    for (int i=0; i<65536; i++)
    	for (int j=0; j<T; j++)
    		histogram[i] += hist[j][i];

    for (int i=0; i<T; i++)
    	delete [] hist[i];
    delete [] hist;

}

void ImProcFunctions::rgbProc (Imagefloat* working, LabImage* lab, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                               SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve) {
    rgbProc (working, lab, hltonecurve, shtonecurve, tonecurve, shmap, sat, rCurve, gCurve, bCurve, params->toneCurve.expcomp, params->toneCurve.hlcompr, params->toneCurve.hlcomprthresh);
}

// Process RGB image and convert to LAB space
void ImProcFunctions::rgbProc (Imagefloat* working, LabImage* lab, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                               SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve,
                               double expcomp, int hlcompr, int hlcomprthresh) {

    int h_th, s_th;
    if (shmap) {
        h_th = shmap->max_f - params->sh.htonalwidth * (shmap->max_f - shmap->avg) / 100;
        s_th = params->sh.stonalwidth * (shmap->avg - shmap->min_f) / 100;
    }

    bool processSH  = params->sh.enabled && shmap!=NULL && (params->sh.highlights>0 || params->sh.shadows>0);
    bool processLCE = params->sh.enabled && shmap!=NULL && params->sh.localcontrast>0;
    double lceamount = params->sh.localcontrast / 200.0;

    TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);

    double toxyz[3][3] = {
        {
        	( wprof[0][0] / Color::D50x),
        	( wprof[0][1] / Color::D50x),
        	( wprof[0][2] / Color::D50x)
        },{
			( wprof[1][0]		),
			( wprof[1][1]		),
			( wprof[1][2]		)
        },{
			( wprof[2][0] / Color::D50z),
			( wprof[2][1] / Color::D50z),
			( wprof[2][2] / Color::D50z)
        }
    };


    bool mixchannels = (params->chmixer.red[0]!=100	|| params->chmixer.red[1]!=0     || params->chmixer.red[2]!=0   ||
						params->chmixer.green[0]!=0 || params->chmixer.green[1]!=100 || params->chmixer.green[2]!=0 ||
						params->chmixer.blue[0]!=0	|| params->chmixer.blue[1]!=0    || params->chmixer.blue[2]!=100);

    int tW = working->width;
    int tH = working->height;
	double pi = M_PI;
	FlatCurve* hCurve;
	FlatCurve* sCurve;
	FlatCurve* vCurve;
	
	
	float* cossq = new float [8192];
	for (int i=0; i<8192; i++) 
		cossq[i] = SQR(cos(pi*(float)i/16384.0));
	
	FlatCurveType hCurveType = (FlatCurveType)params->hsvequalizer.hcurve.at(0);
	FlatCurveType sCurveType = (FlatCurveType)params->hsvequalizer.scurve.at(0);
	FlatCurveType vCurveType = (FlatCurveType)params->hsvequalizer.vcurve.at(0);
	bool hCurveEnabled = hCurveType > FCT_Linear;
	bool sCurveEnabled = sCurveType > FCT_Linear;
	bool vCurveEnabled = vCurveType > FCT_Linear;

	// TODO: We should create a 'skip' value like for CurveFactory::complexsgnCurve (rtengine/curves.cc)
	if (hCurveEnabled) hCurve = new FlatCurve(params->hsvequalizer.hcurve);
	if (sCurveEnabled) sCurve = new FlatCurve(params->hsvequalizer.scurve);
	if (vCurveEnabled) vCurve = new FlatCurve(params->hsvequalizer.vcurve);
	
	const float exp_scale = pow (2.0, expcomp);
	const float comp = (max(0.0, expcomp) + 1.0)*hlcompr/100.0;
	const float shoulder = ((65536.0/max(1.0f,exp_scale))*(hlcomprthresh/200.0))+0.1;
	const float hlrange = 65536.0-shoulder;
	
	
#pragma omp parallel for if (multiThread)
    for (int i=0; i<tH; i++) {

        for (int j=0; j<tW; j++) {

            float r = working->r[i][j];
            float g = working->g[i][j];
            float b = working->b[i][j];
			
			//if (i==100 & j==100) printf("rgbProc input R= %f  G= %f  B= %f  \n",r,g,b);

            if (mixchannels) {
                float rmix = (r*params->chmixer.red[0]   + g*params->chmixer.red[1]   + b*params->chmixer.red[2]) / 100;
                float gmix = (r*params->chmixer.green[0] + g*params->chmixer.green[1] + b*params->chmixer.green[2]) / 100;
                float bmix = (r*params->chmixer.blue[0]  + g*params->chmixer.blue[1]  + b*params->chmixer.blue[2]) / 100;
				
				r = rmix;
				g = gmix;
				b = bmix;
            }

            if (processSH || processLCE) {
                double mapval = 1.0 + shmap->map[i][j];
                double factor = 1.0;
                
                if (processSH) {
                    if (mapval > h_th) 
                        factor = (h_th + (100.0 - params->sh.highlights) * (mapval - h_th) / 100.0) / mapval; 
                    else if (mapval < s_th) 
                        factor = (s_th - (100.0 - params->sh.shadows) * (s_th - mapval) / 100.0) / mapval; 
                }
                if (processLCE) {
                    double sub = lceamount*(mapval-factor*(r*lumimul[0] + g*lumimul[1] + b*lumimul[2]));
                    r = factor*r-sub;
                    g = factor*g-sub;
                    b = factor*b-sub;
                }
                else {
                    r = factor*r;
                    g = factor*g;
                    b = factor*b;
                }
            }

			//TODO: proper treatment of out-of-gamut colors
			//float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
			float tonefactor=((r<MAXVAL ? hltonecurve[r] : CurveFactory::hlcurve (exp_scale, comp, hlrange, r) ) +
							  (g<MAXVAL ? hltonecurve[g] : CurveFactory::hlcurve (exp_scale, comp, hlrange, g) ) +
							  (b<MAXVAL ? hltonecurve[b] : CurveFactory::hlcurve (exp_scale, comp, hlrange, b) ) )/3.0;
			
			r = (r*tonefactor);
			g = (g*tonefactor);
			b = (b*tonefactor);
			
			//shadow tone curve
			float Y = (0.299*r + 0.587*g + 0.114*b);
			tonefactor = shtonecurve[Y];
			r *= tonefactor;
			g *= tonefactor;
			b *= tonefactor;
			
			//brightness/contrast and user tone curve
			r = rCurve[tonecurve[r]];
			g = gCurve[tonecurve[g]];
			b = bCurve[tonecurve[b]];
			
			//if (r<0 || g<0 || b<0) {
			//	printf("negative values row=%d col=%d  r=%f  g=%f  b=%f  \n", i,j,r,g,b);
			//}

			if (sat!=0 || hCurveEnabled || sCurveEnabled || vCurveEnabled) {
				float h,s,v;
				Color::rgb2hsv(r,g,b,h,s,v);
				if (sat > 0.5) {
					s = (1-(float)sat/100)*s+(float)sat/100*(1-SQR(SQR(1-min(s,1.0f))));
					if (s<0) s=0;
				} else {
					if (sat < -0.5)
						s *= 1+(float)sat/100;	
				}
				//HSV equalizer
				if (hCurveEnabled) {
					h = (hCurve->getVal((double)h) - 0.5) * 2 + h;
					if (h > 1.0)
						h -= 1.0;
					else if (h < 0.0)
						h += 1.0;
				}
				if (sCurveEnabled) {
					//shift saturation
					float satparam = (sCurve->getVal((double)h)-0.5) * 2;
					if (satparam > 0.00001) {
						s = (1-satparam)*s+satparam*(1-SQR(1-min(s,1.0f)));
						if (s<0) s=0;
					} else {
						if (satparam < -0.00001)
							s *= 1+satparam;
					}

				}
				if (vCurveEnabled) {
                    if (v<0) v=0;  // important

					//shift value
					float valparam = vCurve->getVal((double)h)-0.5;
					valparam *= (1-SQR(SQR(1-min(s,1.0f))));
					if (valparam > 0.00001) {
						v = (1-valparam)*v+valparam*(1-SQR(1-min(v,1.0f)));
						if (v<0) v=0;
					} else {
						if (valparam < -0.00001)
							v *= (1+valparam);
					}
					
				}
				Color::hsv2rgb(h,s,v,r,g,b);
			}
			 
            float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
            float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
            float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;
			
			float fx,fy,fz;
			
			//if (x>0) {
				fx = (x<65535.0 ? cachef[x] : (327.68*exp(log(x/MAXVAL)/3.0 )));
			//} else {
			//	fx = (x>-65535.0 ? -cachef[-x] : (-327.68*exp(log(-x/MAXVAL)/3.0 )));
			//}
			//if (y>0) {
				fy = (y<65535.0 ? cachef[y] : (327.68*exp(log(y/MAXVAL)/3.0 )));
			//} else {
			//	fy = (y>-65535.0 ? -cachef[-y] : (-327.68*exp(log(-y/MAXVAL)/3.0 )));
			//}
			//if (z>0) {
				fz = (z<65535.0 ? cachef[z] : (327.68*exp(log(z/MAXVAL)/3.0 )));
			//} else {
			//	fz = (z>-65535.0 ? -cachef[-z] : (-327.68*exp(log(-z/MAXVAL)/3.0 )));
			//}

			lab->L[i][j] = (116.0 * fy - 5242.88); //5242.88=16.0*327.68;
            lab->a[i][j] = (500.0 * (fx - fy) );
            lab->b[i][j] = (200.0 * (fy - fz) );
			

			
			//test for color accuracy
			/*float fy = (0.00862069 * lab->L[i][j])/327.68 + 0.137932; // (L+16)/116
			float fx = (0.002 * lab->a[i][j])/327.68 + fy;
			float fz = fy - (0.005 * lab->b[i][j])/327.68;
			
			float x_ = 65535*Lab2xyz(fx)*Color::D50x;
			float y_ = 65535*Lab2xyz(fy);
			float z_ = 65535*Lab2xyz(fz)*Color::D50z;
			
			int R,G,B;
			xyz2srgb(x_,y_,z_,R,G,B);
			r=(float)R; g=(float)G; b=(float)B;
			float xxx=1;*/

        }
    }
	
	if (hCurveEnabled) delete hCurve;
	if (sCurveEnabled) delete sCurve;
	if (vCurveEnabled) delete vCurve;
	delete [] cossq;
 }

void ImProcFunctions::luminanceCurve (LabImage* lold, LabImage* lnew, LUTf & curve) {

    int W = lold->W;
    int H = lold->H;

#pragma omp parallel for if (multiThread)
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
			float Lin=lold->L[i][j];
			//if (Lin>0 && Lin<65535)
				lnew->L[i][j] = curve[Lin];
		}
}

void ImProcFunctions::chromiLuminanceCurve (LabImage* lold, LabImage* lnew, LUTf & acurve, LUTf & bcurve, LUTf & satcurve/*,LUTf & satbgcurve*/, LUTf & curve, bool utili, bool autili, bool butili, bool ccutili) {
	
	int W = lold->W;
	int H = lold->H;
	//init Flatcurve for C=f(H)
	FlatCurve* chCurve = NULL;
	bool chutili = false;
	if (!params->labCurve.bwtoning) {
		chCurve = new FlatCurve(params->labCurve.chcurve);
		if (chCurve && !chCurve->isIdentity()) {
			chutili=true;
		}//do not use "Munsell" if Chcurve not used
	}

#ifdef _DEBUG
	MyTime t1e,t2e, t3e, t4e;
	t1e.set();
	// init variables to display Munsell corrections
	MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif

	unsigned int N = W*H;
	float *L = lold->L[0];
	float *a=  lold->a[0];
	float *b=  lold->b[0];

	float* Lold = new float [lold->W*lold->H];//to save L before any used
	float* Cold = new float [lold->W*lold->H];//to save C before any used
	float adjustr=1.0f, adjustbg=1.0f;

//	if(params->labCurve.avoidclip ){
	for (unsigned int j=0; j!=N; j++){
		Lold[j]=L[j]/327.68f;
		Cold[j]=sqrt(SQR(a[j]/327.68f)+SQR(b[j]/327.68f));
//		Hr=atan2(b[j],a[j]);
//		if(Hr >-0.15f && Hr < 1.5f && Cold[j]>maxCr)
//			maxCr=Cold[j];	// I do not take into account  "acurve" and "bcurve" to adjust max
//		 else if (Cold[j]>maxCbg)
//			maxCbg=Cold[j];
	}
	// parameter to adapt curve C=f(C) to gamut
	
	if      (params->icm.working=="ProPhoto")   {adjustr =       adjustbg = 1.2f;}// 1.2 instead 1.0 because it's very rare to have C>170..
	else if (params->icm.working=="Adobe RGB")  {adjustr = 1.8f; adjustbg = 1.4f;}
	else if (params->icm.working=="sRGB")  	    {adjustr = 2.0f; adjustbg = 1.7f;}
	else if (params->icm.working=="WideGamut")  {adjustr =       adjustbg = 1.2f;}
	else if (params->icm.working=="Beta RGB")   {adjustr =       adjustbg = 1.4f;}
	else if (params->icm.working=="BestRGB")    {adjustr =       adjustbg = 1.4f;}
	else if (params->icm.working=="BruceRGB")   {adjustr = 1.8f; adjustbg = 1.5f;}


	// reference to the params structure has to be done outside of the parallelization to avoid CPU cache problem
	bool highlight = params->hlrecovery.enabled; //Get the value if "highlight reconstruction" is activated
	int chromaticity = params->labCurve.chromaticity;
	bool bwToning = params->labCurve.bwtoning;
	double rstprotection = 100.-params->labCurve.rstprotection; // Red and Skin Tones Protection
	// avoid color shift is disabled when bwToning is activated
	bool avoidColorShift = params->labCurve.avoidcolorshift && !bwToning;
	int protectRed = settings->protectred;
	double protectRedH = settings->protectredh;
	bool gamutLch = settings->gamutLch;

	// only if user activate Lab adjustements
	if (avoidColorShift) {
		if(autili || butili || ccutili || chutili || utili || chromaticity)
			Color::LabGamutMunsell(lold, Lold, Cold, /*corMunsell*/true, /*lumaMuns*/false, params->hlrecovery.enabled, /*gamut*/true, params->icm.working, multiThread);
	}


#ifdef _DEBUG
#pragma omp parallel default(shared) firstprivate(highlight, chromaticity, bwToning, rstprotection, avoidColorShift, protectRed, protectRedH, gamutLch, lold, lnew, MunsDebugInfo) if (multiThread)
#else
#pragma omp parallel default(shared) firstprivate(highlight, chromaticity, bwToning, rstprotection, avoidColorShift, protectRed, protectRedH, gamutLch, lold, lnew) if (multiThread)
#endif
{

	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);
	//if(utili) curve.dump("Lcurve");

	double wp[3][3] = {
		{wprof[0][0],wprof[0][1],wprof[0][2]},
		{wprof[1][0],wprof[1][1],wprof[1][2]},
		{wprof[2][0],wprof[2][1],wprof[2][2]}};

	//float maxlp=-100.0, minlp=200.0;

#pragma omp for schedule(dynamic, 10)
	for (int i=0; i<H; i++)
		for (int j=0; j<W; j++) {
			float LL=lold->L[i][j]/327.68f;
			float CC=sqrt(SQR(lold->a[i][j]/327.68f) + SQR(lold->b[i][j]/327.68f));
			float HH=atan2(lold->b[i][j],lold->a[i][j]);
			float Chprov=CC;
			float Chprov1=CC;
			float memChprov=Chprov;
			float Lprov2=LL;
			float Lin=lold->L[i][j];
			lnew->L[i][j] = curve[Lin];
			float Lprov1=(lnew->L[i][j])/327.68f;
			float chromaChfactor=1.0f;
			float atmp = acurve[lold->a[i][j]+32768.0f]-32768.0f;// curves Lab a
			float btmp = bcurve[lold->b[i][j]+32768.0f]-32768.0f;// curves Lab b
			//float chromaCredfactor=1.0f;
			//float chromaCbgfactor=1.0f;
//			chromaCfactor=(satcurve[chroma*adjustr])/(chroma*adjustr);//apply C=f(C)
//			chromaCbgfactor=(satbgcurve[chroma*adjustbg])/(chroma*adjustbg);
//			calculate C=f(H)
			if (chutili) {
				double hr;
				//hr=translate Hue Lab value  (-Pi +Pi) in approximative hr (hsv values) (0 1) [red 1/6 yellow 1/6 green 1/6 cyan 1/6 blue 1/6 magenta 1/6 ]
				// with multi linear correspondances (I expect there is no error !!)
				if      (HH<-2.7f) hr=0.020380804*double(HH)+0.970281708; //Lab green                   =>hr # 0.33 ==> 0.33  0.42
				else if (HH<-2.1f) hr=0.266666667*double(HH)+1.14;        //Lab cyan                    =>hr # 0.50 ==> 0.42  0.58
				else if (HH<-0.9f) hr=0.141666   *double(HH)+0.8775;      //Lab blue                    =>hr # 0.67 ==> 0.58  0.75
				else if (HH<-0.1f) hr=0.2125     *double(HH)+0.94125;     //Lab magenta (purple)        =>hr # 0.83 ==> 0.75  0.92
				else if (HH< 1.3f) hr=0.12142857 *double(HH)+0.932142857; //Lab red and skin            =>hr # 0    ==> 0.92  1.09
				else if (HH< 2.2f) hr=0.1666667  *double(HH)-0.1266667;   //Lab yellow and green yellow =>hr # 0.16 ==> 0.09  0.24
				else               hr=0.0955828  *double(HH)+0.02971784;  //Lab green                   =>hr # 0.33 ==> 0.24  0.33

				//allways put h between 0 and 1
				if     (hr<0.0) hr += 1.0;
				else if(hr>1.0) hr -= 1.0;
				float chparam = float((chCurve->getVal(hr)-0.5f) * 2.0f);//get C=f(H)
				chromaChfactor=1.0f+chparam;
			}
			atmp *= chromaChfactor;//apply C=f(H)
			btmp *= chromaChfactor;
//			if (params->labCurve.chromaticity) {// if user use sliders
			if(chromaticity!=0 && !bwToning){
				// approximation in Lab mode to protect skin tones and avoid too big gamut clip for red
				float scale = 100.0f/100.1f;//reduction in normal zone
				float scaleext=1.0f;//reduction in transition zone
				float protect_red,protect_redh;
				float deltaHH;//HH value transition
				float dred=55.0f;//C red value limit
				protect_red=float(protectRed);//default=60  chroma: one can put more or less if necessary...in 'option'  40...160
				if(protect_red < 20.0f) protect_red=20.0; // avoid too low value
				if(protect_red > 180.0f) protect_red=180.0; // avoid too high value
				protect_redh=float(protectRedH);//default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0
				if(protect_redh<0.1f) protect_redh=0.1f;//avoid divide by 0 and negatives values
				if(protect_redh>1.0f) protect_redh=1.0f;//avoid too big values

				deltaHH=protect_redh;//transition hue

				//simulate very approximative gamut f(L) : with pyramid transition
				if     (Lprov1<25.0f)   dred = 40.0f;
				else if(Lprov1<30.0f)   dred = 3.0f*Lprov1 -35.0f;
				else if(Lprov1<70.0f)   dred = 55.0f;
				else if(Lprov1<75.0f)   dred = -3.0f*Lprov1 +265.0f;
				else                    dred = 40.0f;

				if(rstprotection<99.9999) {
					if(chromaticity>0)
						scale = rstprotection/100.1f;
					if((HH< (1.3f+deltaHH) && HH >=1.3f))
						scaleext=HH*(1.0f-scale)/deltaHH + 1.0f - (1.3f+deltaHH)*(1.0f-scale)/deltaHH;    //transition for Hue (red - yellow)
					else if((HH< 0.15f && HH >(0.15f-deltaHH)))
						scaleext=HH*(scale-1.0f)/deltaHH + 1.0f - (0.15f-deltaHH)*(scale-1.0f)/deltaHH;   //transition for hue (red purple)
				}

				//transition for red , near skin tones
				float factorskin, factorsat, factor, factorskinext;
				factorskin=1.0f+(chromaticity*scale)/100.0f;
				factorskinext=1.0f+(chromaticity*scaleext)/100.0f;
				factorsat=1.0f+(chromaticity)/100.0f;/*if(factorsat==1.0f) factorsat=1.1f;*/

				factor = factorsat;
				// Test if chroma is in the normal range first
				if(HH>=0.15f && HH<1.3f) {
					if (Chprov1<dred)
						factor = factorskin;
					else if(Chprov1<(dred+protect_red))
						factor = (factorsat-factorskin)/protect_red*Chprov1+factorsat-(dred+protect_red)*(factorsat-factorskin)/protect_red;
				}
				// then test if chroma is in the extanded range
				else if ( HH>(0.15f-deltaHH) || HH<(1.3f+deltaHH) ) {
					if (Chprov1 < dred)
						factor = factorskinext;// C=dred=55 => real max of skin tones
					else if (Chprov1 < (dred+protect_red))// transition
						factor = (factorsat-factorskinext)/protect_red*Chprov1+factorsat-(dred+protect_red)*(factorsat-factorskinext)/protect_red;
				}

				atmp *= factor;
				btmp *= factor;
				// end approximation
 			}

			// I have placed C=f(C) after all C treatments to assure maximum amplitude of "C"
			if (!bwToning) {
				float chroma = sqrt(SQR(atmp)+SQR(btmp)+0.001f);
				float chromaCfactor = (satcurve[chroma*adjustr])/(chroma*adjustr);//apply C=f(C)
				atmp *= chromaCfactor;
				btmp *= chromaCfactor;
			}
			// end chroma C=f(C)

			Chprov1 = sqrt(SQR(atmp/327.68f)+SQR(btmp/327.68f));

/*
			// modulation of a and b curves with saturation
			if (params->labCurve.chromaticity!=0 && !params->labCurve.bwtoning) {
				float chroma = sqrt(SQR(atmp)+SQR(btmp)+0.001);
				float satfactor = (satcurve[chroma+32768.0f]-32768.0f)/chroma;
				atmp *= satfactor;
				btmp *= satfactor;
			}
*/

			// labCurve.bwtoning option allows to decouple modulation of a & b curves by saturation
			// with bwtoning enabled the net effect of a & b curves is visible
			if (bwToning) {
				atmp -= lold->a[i][j];
				btmp -= lold->b[i][j];
			}

			if (avoidColorShift) {
				//gamutmap Lch ==> preserve Hue,but a little slower than gamutbdy for high values...and little faster for low values
				if(gamutLch) {
					float R,G,B;

#ifdef _DEBUG
					bool neg=false;
					bool more_rgb=false;
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.4f, 0.95f, neg, more_rgb);
#else
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.4f, 0.95f);
#endif
					Lprov2 = Lprov1;

					lnew->L[i][j]=Lprov1*327.68f;
					lnew->a[i][j]=327.68f*Chprov1*cos(HH);
					lnew->b[i][j]=327.68f*Chprov1*sin(HH);
				}
				else {
					//use gamutbdy
					//Luv limiter
					float Y,u,v;
					Color::Lab2Yuv(lnew->L[i][j],atmp,btmp,Y,u,v);
					//Yuv2Lab includes gamut restriction map
					Color::Yuv2Lab(Y,u,v,lnew->L[i][j],lnew->a[i][j],lnew->b[i][j], wp);
				}

				if (utili || autili || butili || ccutili || chutili || chromaticity) {
					float correctionHue=0.0f; // Munsell's correction
					float correctlum=0.0f;

					Lprov1=lnew->L[i][j]/327.68f;
					Chprov=sqrt(SQR(lnew->a[i][j]/327.68f)+ SQR(lnew->b[i][j]/327.68f));

#ifdef _DEBUG
					Color::AllMunsellLch(/*lumaMuns*/true, Lprov1,Lprov2,HH,Chprov,memChprov,correctionHue,correctlum, MunsDebugInfo);
#else
					Color::AllMunsellLch(/*lumaMuns*/true, Lprov1,Lprov2,HH,Chprov,memChprov,correctionHue,correctlum);
#endif

					if(fabs(correctionHue) < 0.015f) HH+=correctlum;	// correct only if correct Munsell chroma very little.

					lnew->a[i][j]=327.68f*Chprov*cos(HH+correctionHue);// apply Munsell
					lnew->b[i][j]=327.68f*Chprov*sin(HH+correctionHue);
				}
			}
			else {
//				if(Lprov1 > maxlp) maxlp=Lprov1;
//				if(Lprov1 < minlp) minlp=Lprov1;
				lnew->L[i][j]=Lprov1*327.68f;

				//Luv limiter only
				lnew->a[i][j] = atmp;
				lnew->b[i][j] = btmp;
			}
		}
} // end of parallelization

#ifdef _DEBUG
	if (settings->verbose) {
		t3e.set();
		printf("Color::AllMunsellLch (correction performed in %d usec):\n", t3e.etime(t1e));
		printf("   Munsell chrominance: MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
		printf("   Munsell luminance  : MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhuelum[0], MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
	}
	delete MunsDebugInfo;
#endif
	delete [] Lold;
	delete [] Cold;

	if (chutili) delete chCurve;
}


//#include "cubic.cc"

void ImProcFunctions::colorCurve (LabImage* lold, LabImage* lnew) {

/*    LUT<double> cmultiplier(181021);

    double boost_a = ((float)params->colorBoost.amount + 100.0) / 100.0;
    double boost_b = ((float)params->colorBoost.amount + 100.0) / 100.0;

    double c, amul = 1.0, bmul = 1.0;
    if (boost_a > boost_b) {
        c = boost_a;
        if (boost_a > 0)
            bmul = boost_b / boost_a;
    }
    else {
        c = boost_b;
        if (boost_b > 0)
            amul = boost_a / boost_b;
    }

    if (params->colorBoost.enable_saturationlimiter && c>1.0) {
        // re-generate color multiplier lookup table
        double d = params->colorBoost.saturationlimit / 3.0;
        double alpha = 0.5;
        double threshold1 = alpha * d;
        double threshold2 = c*d*(alpha+1.0) - d;
        for (int i=0; i<=181020; i++) { // lookup table stores multipliers with a 0.25 chrominance resolution
            double chrominance = (double)i/4.0;
            if (chrominance < threshold1)
                cmultiplier[i] = c;
            else if (chrominance < d)
                cmultiplier[i] = (c / (2.0*d*(alpha-1.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
            else if (chrominance < threshold2) 
                cmultiplier[i] = (1.0 / (2.0*d*(c*(alpha+1.0)-2.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
            else
                cmultiplier[i] = 1.0;
        }
    }
    
	float eps = 0.001;
    double shift_a = params->colorShift.a + eps, shift_b = params->colorShift.b + eps;

    float** oa = lold->a;
    float** ob = lold->b;

	#pragma omp parallel for if (multiThread)
    for (int i=0; i<lold->H; i++)
        for (int j=0; j<lold->W; j++) {

            double wanted_c = c;
            if (params->colorBoost.enable_saturationlimiter && c>1) {
                float chroma = (float)(4.0 * sqrt((oa[i][j]+shift_a)*(oa[i][j]+shift_a) + (ob[i][j]+shift_b)*(ob[i][j]+shift_b)));
                wanted_c = cmultiplier [chroma];
            }

            double real_c = wanted_c;
            if (wanted_c >= 1.0 && params->colorBoost.avoidclip) {
                double cclip = 100000.0;
                double cr = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, 3.079935, -1.5371515, -0.54278342);
                double cg = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, -0.92123418, 1.87599, 0.04524418);
                double cb = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, 0.052889682, -0.20404134, 1.15115166);
                if (cr>1.0 && cr<cclip) cclip = cr;
                if (cg>1.0 && cg<cclip) cclip = cg;
                if (cb>1.0 && cb<cclip) cclip = cb;
                if (cclip<100000.0) {
                    real_c = -cclip + 2.0*cclip / (1.0+exp(-2.0*wanted_c/cclip));
                    if (real_c<1.0)
                        real_c = 1.0;
                }
            }
            
            float nna = ((oa[i][j]+shift_a) * real_c * amul);
            float nnb = ((ob[i][j]+shift_b) * real_c * bmul);
            lnew->a[i][j] = LIM(nna,-32000.0f,32000.0f);
            lnew->b[i][j] = LIM(nnb,-32000.0f,32000.0f);
        }
*/
    //delete [] cmultiplier;
}
	
	void ImProcFunctions::impulsedenoise (LabImage* lab) {
		
		if (params->impulseDenoise.enabled && lab->W>=8 && lab->H>=8)
			
			impulse_nr (lab, (float)params->impulseDenoise.thresh/20.0 );
	}
	
	void ImProcFunctions::defringe (LabImage* lab) {
		
		if (params->defringe.enabled && lab->W>=8 && lab->H>=8)
			
			PF_correct_RT(lab, lab, params->defringe.radius, params->defringe.threshold);
	}
	
	void ImProcFunctions::dirpyrdenoise (LabImage* lab) {
		
		if (params->dirpyrDenoise.enabled && lab->W>=8 && lab->H>=8)
			
			dirpyrLab_denoise(lab, lab, params->dirpyrDenoise );
	}
	
	void ImProcFunctions::dirpyrequalizer (LabImage* lab) {
		
		if (params->dirpyrequalizer.enabled && lab->W>=8 && lab->H>=8) {
			
			//dirpyrLab_equalizer(lab, lab, params->dirpyrequalizer.mult);
			dirpyr_equalizer(lab->L, lab->L, lab->W, lab->H, params->dirpyrequalizer.mult);
			
		}
	}

//Map tones by way of edge preserving decomposition. Is this the right way to include source?
#include "EdgePreservingDecomposition.cc"
void ImProcFunctions::EPDToneMap(LabImage *lab, unsigned int Iterates, int skip){
	//Hasten access to the parameters.
	EPDParams *p = (EPDParams *)(&params->edgePreservingDecompositionUI);

	//Enabled? Leave now if not.
	if(!p->enabled) return;

	//Pointers to whole data and size of it.
	float *L = lab->L[0];
	float *a = lab->a[0];
	float *b = lab->b[0];
	unsigned int i, N = lab->W*lab->H;

	EdgePreservingDecomposition epd = EdgePreservingDecomposition(lab->W, lab->H);

	//Due to the taking of logarithms, L must be nonnegative. Further, scale to 0 to 1 using nominal range of L, 0 to 15 bit.
	float minL = FLT_MAX;
	for(i = 0; i != N; i++)
		if(L[i] < minL) minL = L[i];
	if(minL > 0.0f) minL = 0.0f;		//Disable the shift if there are no negative numbers. I wish there were just no negative numbers to begin with.

	for(i = 0; i != N; i++)
		L[i] = (L[i] - minL)/32767.0f;

	//Some interpretations.
	float Compression = expf(-p->Strength);		//This modification turns numbers symmetric around 0 into exponents.
	float DetailBoost = p->Strength;
	if(p->Strength < 0.0f) DetailBoost = 0.0f;	//Go with effect of exponent only if uncompressing.

	//Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
	if(Iterates == 0) Iterates = (unsigned int)(p->EdgeStopping*15.0);

/* Debuggery. Saves L for toying with outside of RT.
char nm[64];
sprintf(nm, "%ux%ufloat.bin", lab->W, lab->H);
FILE *f = fopen(nm, "wb");
fwrite(L, N, sizeof(float), f);
fclose(f);*/

	epd.CompressDynamicRange(L, (float)p->Scale/skip, (float)p->EdgeStopping, Compression, DetailBoost, Iterates, p->ReweightingIterates, L);

	//Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
	float s = (1.0f + 38.7889f)*powf(Compression, 1.5856f)/(1.0f + 38.7889f*powf(Compression, 1.5856f));
	for(i = 0; i != N; i++)
		a[i] *= s,
		b[i] *= s,
		L[i] = L[i]*32767.0f + minL;
}

	
	void ImProcFunctions::getAutoExp  (LUTu & histogram, int histcompr, double defgain, double clip,
									   double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh) {
		
		float scale = 65536.0;
		float midgray=0.15;//0.18445f;//middle gray in linear gamma = 0.18445*65535
		
		int imax=65536>>histcompr;
		
		float sum = 0, hisum=0, losum=0;
		float ave = 0, hidev=0, lodev=0;
		//find average luminance
		for (int i=0; i<imax; i++) {
			sum += histogram[i];
			ave += histogram[i] * i;
		}
		ave /= (sum);
		
		//find median of luminance
		int median=0, count=histogram[0];
		while (count<sum/2) {
			median++;
			count += histogram[median];
		}
		if (median==0 || ave<1) {//probably the image is a blackframe
			expcomp=0;
			black=0;
			bright=0;
			contr=0;
			hlcompr=0;
			hlcomprthresh=0;
			return;
		}
		
		// compute std dev on the high and low side of median 
		// and octiles of histogram
		float octile[8]={0,0,0,0,0,0,0,0},ospread=0;
		count=0;
		for (int i=0; i<imax; i++) {
			if (count<8) {
				octile[count] += histogram[i];
				if (octile[count]>sum/8 || (count==7 && octile[count]>sum/16)) {
					octile[count]=log(1+i)/log(2);
					count++;// = min(count+1,7);
				}
			}
			if (i<ave) {
				//lodev += SQR(ave-i)*histogram[i];
				lodev += (log(ave+1)-log(i+1))*histogram[i];
				losum += histogram[i];
			} else {
				//hidev += SQR(i-ave)*histogram[i];
				hidev += (log(i+1)-log(ave+1))*histogram[i];
				hisum += histogram[i];
			}
			
		}
		if (losum==0 || hisum==0) {//probably the image is a blackframe
			expcomp=0;
			black=0;
			bright=0;
			contr=0;
			hlcompr=0;
			hlcomprthresh=0;
			return;
		}
		lodev = (lodev/(log(2)*losum));
		hidev = (hidev/(log(2)*hisum));
		if (octile[7]>log(imax+1)/log2(2)) {
			octile[7]=1.5*octile[6]-0.5*octile[5];
		}
		// compute weighted average separation of octiles
		// for future use in contrast setting
		for (int i=1; i<6; i++) {
			ospread += (octile[i+1]-octile[i])/max(0.5f,(i>2 ? (octile[i+1]-octile[3]) : (octile[3]-octile[i])));
		}
		ospread /= 5;
		if (ospread<=0) {//probably the image is a blackframe
			expcomp=0;
			black=0;
			bright=0;
			contr=0;
			hlcompr=0;
			hlcomprthresh=0;
			return;
		}
		
		
		// compute clipping points based on the original histograms (linear, without exp comp.)
		int clipped = 0;
		int rawmax = (imax) - 1;
		while (rawmax>1 && histogram[rawmax]+clipped <= 0) {
			clipped += histogram[rawmax];
			rawmax--;
		}
		
		//compute clipped white point
		int clippable = (int)(sum * clip );
		clipped = 0;
		int whiteclip = (imax) - 1;
		while (whiteclip>1 && histogram[whiteclip]+clipped <= clippable) {
			clipped += histogram[whiteclip];
			whiteclip--;
		}
		
		//compute clipped black point
		clipped = 0;
		int shc = 0;
		while (shc<whiteclip-1 && histogram[shc]+clipped <= clippable) {
			clipped += histogram[shc];
			shc++;
		}
		
		//rescale to 65535 max
		rawmax <<= histcompr;
		whiteclip <<= histcompr;
		ave = ave*(1<<histcompr);
		median <<= histcompr;
		shc <<= histcompr;
		
		//prevent division by 0
		if (lodev==0) lodev=1;
		
		//compute exposure compensation as geometric mean of the amount that
		//sets the mean or median at middle gray, and the amount that sets the estimated top 
		//of the histogram at or near clipping.  
		
		float expcomp1 = log(/*(median/ave)*//*(hidev/lodev)*/midgray*scale/(ave-shc+midgray*shc))/log(2);
		float expcomp2 = 0.5*( (15.5f-histcompr-(2*octile[7]-octile[6])) + log(scale/rawmax)/log(2) );

		/*expcomp = (expcomp1*fabs(expcomp2)+expcomp2*fabs(expcomp1))/(fabs(expcomp1)+fabs(expcomp2));
		if (expcomp<0) {
			min(0.0f,max(expcomp1,expcomp2));
		}*/
		expcomp = 0.5 * (expcomp1 + expcomp2);
		float gain = exp(expcomp*log(2));
		
		
		gain = /*(median/ave)*/sqrt(gain*scale/rawmax);
		black = shc*gain;
		expcomp = log(gain)/log(2.0);//convert to stops
		
		black = shc*gain;
		
		//now tune hlcompr to bring back rawmax to 65535
		hlcomprthresh = 33;
		//this is a series approximation of the actual formula for comp,
		//which is a transcendental equation
		float comp = (gain*((float)whiteclip)/scale - 1)*2;//*(1-shoulder/scale);
		hlcompr=(int)(100*comp/(max(0.0,expcomp) + 1.0));
		hlcompr = max(0,min(100,hlcompr));
		
		//now find brightness if gain didn't bring ave to midgray using
		//the envelope of the actual 'control cage' brightness curve for simplicity
		float midtmp = gain*sqrt(median*ave)/scale;
		if (midtmp<0.1) {
			bright = (midgray-midtmp)*15.0/(midtmp);
		} else {
			bright = (midgray-midtmp)*15.0/(0.10833-0.0833*midtmp);
		}
		
		bright = 0.25*/*(median/ave)*(hidev/lodev)*/max(0,bright);
		
		//compute contrast that spreads the average spacing of octiles
		contr = 50.0*(1.1-ospread);
		contr = max(0,min(100,contr));
		
		//diagnostics
		//printf ("**************** AUTO LEVELS ****************\n");
		//printf ("gain1= %f   gain2= %f		gain= %f\n",expcomp1,expcomp2,gain);
		//printf ("median: %i  average: %f    median/average: %f\n",median,ave, median/ave);
		//printf ("average: %f\n",ave);
		//printf ("median/average: %f\n",median/ave);
		//printf ("lodev: %f   hidev: %f		hidev/lodev: %f\n",lodev,hidev,hidev/lodev);
		//printf ("lodev: %f\n",lodev);
		//printf ("hidev: %f\n",hidev);
		//printf ("rawmax= %d  whiteclip= %d  gain= %f\n",rawmax,whiteclip,gain);
		
		//printf ("octiles: %f %f %f %f %f %f %f %f\n",octile[0],octile[1],octile[2],octile[3],octile[4],octile[5],octile[6],octile[7]);    
		//printf ("ospread= %f\n",ospread);
		
		
		/*
		 // %%%%%%%%%% LEGACY AUTOEXPOSURE CODE %%%%%%%%%%%%%
		 // black point selection is based on the linear result (yielding better visual results)
		 black = (int)(shc * corr);
		 // compute the white point of the exp. compensated gamma corrected image
		 double whiteclipg = (int)(CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);
		 
		 // compute average intensity of the exp compensated, gamma corrected image
		 double gavg = 0;
		 for (int i=0; i<65536>>histcompr; i++) 
		 gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;
		 
		 if (black < gavg) {
		 int maxwhiteclip = (gavg - black) * 4 / 3 + black; // dont let whiteclip be such large that the histogram average goes above 3/4
		 //double mavg = 65536.0 / (whiteclipg-black) * (gavg - black);
		 if (whiteclipg < maxwhiteclip)
		 whiteclipg = maxwhiteclip;
		 }
		 
		 whiteclipg = CurveFactory::igamma2 ((float)(whiteclipg/65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter
		 
		 black = (int)((65535*black)/whiteclipg);
		 expcomp = log(65535.0 / (whiteclipg)) / log(2.0);
		 
		 if (expcomp<0.0)	expcomp = 0.0;*/
		
		if (expcomp<-5.0)	expcomp = -5.0;
		if (expcomp>10.0)	expcomp = 10.0;
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	double ImProcFunctions::getAutoDistor  (const Glib::ustring &fname, int thumb_size) {
		if (fname != "") {
			rtengine::RawMetaDataLocation ri;
			int w_raw=-1, h_raw=thumb_size;
			int w_thumb=-1, h_thumb=thumb_size;
			
			Thumbnail* thumb = rtengine::Thumbnail::loadQuickFromRaw (fname, ri, w_thumb, h_thumb, 1, FALSE);
			if (thumb == NULL)
				return 0.0;
			
			Thumbnail* raw =   rtengine::Thumbnail::loadFromRaw      (fname, ri, w_raw, h_raw, 1, FALSE);
			if (raw == NULL) {
				delete thumb;
				return 0.0;
			}
			
			if (h_thumb != h_raw) {
				delete thumb;
				delete raw;
				return 0.0;
			}
			
			int width;
			
			if (w_thumb > w_raw)
				width = w_raw;
			else
				width = w_thumb;
			
			unsigned char* thumbGray;
			unsigned char* rawGray;
			thumbGray = thumb->getGrayscaleHistEQ (width);
			rawGray = raw->getGrayscaleHistEQ (width);
			
			if (thumbGray == NULL || rawGray == NULL) {
				if (thumbGray) delete thumbGray;
				if (rawGray) delete rawGray;
				delete thumb;
				delete raw;
				return 0.0;
			}
			
			double dist_amount = calcDistortion (thumbGray, rawGray, width, h_thumb);
			delete thumbGray;
			delete rawGray;
			delete thumb;
			delete raw;
			return dist_amount;
		}
		else
			return 0.0;
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

}
