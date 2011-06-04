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
#include <improcfun.h>
#include <curves.h>
#include <math.h>
#include <colorclip.h>
#include <gauss.h>
#include <bilateral2.h>
#include <minmax.h>
#include <mytime.h>
#include <glib.h>
#include <glibmm.h>
#include <iccstore.h>
#include <impulse_denoise.h>
#include <imagesource.h>
#include <rtthumbnail.h>
#include <utils.h>
#include <iccmatrices.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {
	
	using namespace procparams;
	
#undef MAXVAL
#undef CMAXVAL
#undef MAXL
#undef MAX
#undef MIN
#undef ABS
#undef CLIP
#undef CLIPS
#undef CLIPC
#undef CLIPTO
	
#define MAXVAL  0xffff
#define CMAXVAL 0xffff
#define MAXL 	0xffff
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define ABS(a) ((a)<0?-(a):(a))
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPS(a) ((a)>-32768?((a)<32767?(a):32767):-32768)
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIP2(a) ((a)<MAXVAL ? a : MAXVAL )
#define FCLIP(a) ((a)>0.0?((a)<65535.5?(a):65535.5):0.0)
	
#define D50x 0.96422
#define D50z 0.82521
#define u0 4.0*D50x/(D50x+15+3*D50z)
#define v0 9.0/(D50x+15+3*D50z)
	
#define eps_max 580.40756 //(MAXVAL* 216.0f/24389.0);
#define kappa	903.29630 //24389.0/27.0;
	
	
extern const Settings* settings;

LUTf ImProcFunctions::cachef ;
LUTf ImProcFunctions::gamma2curve = 0;

void ImProcFunctions::initCache () {

    int maxindex = 65536;
	cachef(maxindex,0/*LUT_CLIP_BELOW*/);

	gamma2curve(maxindex,0);

    for (int i=0; i<maxindex; i++) {
        if (i>eps_max) {
			cachef[i] = 327.68*( exp(1.0/3.0 * log((double)i / MAXVAL) ));
        }
        else {
			cachef[i] = 327.68*((kappa*i/MAXVAL+16.0)/116.0);
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
		monitorTransform = cmsCreateTransform (iprof, TYPE_RGB_FLT, monitor, TYPE_RGB_8, settings->colorimetricIntent,
            settings->LCMSSafeMode ? cmsFLAGS_NOOPTIMIZE : cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
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

    int H = original->height;
#ifdef _OPENMP
	#pragma omp parallel if (multiThread)
    {
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

// Process RGB image and convert to LAB space
void ImProcFunctions::rgbProc (Imagefloat* working, LabImage* lab, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve, SHMap* shmap, int sat) {

    int h_th, s_th;
    if (shmap) {
        h_th = shmap->max - params->sh.htonalwidth * (shmap->max - shmap->avg) / 100;
        s_th = params->sh.stonalwidth * (shmap->avg - shmap->min) / 100;
    }

    bool processSH  = params->sh.enabled && shmap!=NULL && (params->sh.highlights>0 || params->sh.shadows>0);
    bool processLCE = params->sh.enabled && shmap!=NULL && params->sh.localcontrast>0;
    double lceamount = params->sh.localcontrast / 200.0;

    TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);

	float toxyz[3][3] = {
        {
        	( wprof[0][0] / D50x),
        	( wprof[0][1] / D50x),
        	( wprof[0][2] / D50x)
        },{
			( wprof[1][0]		),
			( wprof[1][1]		),
			( wprof[1][2]		)
        },{
			( wprof[2][0] / D50z),
			( wprof[2][1] / D50z),
			( wprof[2][2] / D50z)
        }
    };


    bool mixchannels = (params->chmixer.red[0]!=100	|| params->chmixer.red[1]!=0     || params->chmixer.red[2]!=0   || \
						params->chmixer.green[0]!=0 || params->chmixer.green[1]!=100 || params->chmixer.green[2]!=0 || \
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
	
	const float exp_scale = pow (2.0, params->toneCurve.expcomp);
	const float comp = (params->toneCurve.expcomp + 1.0)*params->toneCurve.hlcompr/100.0;
	const float shoulder = ((65536.0/exp_scale)*(params->toneCurve.hlcomprthresh/200.0))+0.1;
	const float hlrange = 65536.0-shoulder;
	
	
#pragma omp parallel for if (multiThread)
    for (int i=0; i<tH; i++) {

        for (int j=0; j<tW; j++) {

            float r = working->r[i][j];
            float g = working->g[i][j];
            float b = working->b[i][j];
			
			//if (i==100 & j==100) printf("rgbProc input R= %f  G= %f  B= %f  \n",r,g,b);

            if (mixchannels) {
                r = (r*params->chmixer.red[0]   + g*params->chmixer.red[1]   + b*params->chmixer.red[2]) / 100;
                g = (r*params->chmixer.green[0] + g*params->chmixer.green[1] + b*params->chmixer.green[2]) / 100;
                b = (r*params->chmixer.blue[0]  + g*params->chmixer.blue[1]  + b*params->chmixer.blue[2]) / 100;

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
                    r = CLIP(factor*r-sub);
                    g = CLIP(factor*g-sub);
                    b = CLIP(factor*b-sub);
                }
                else {
                    r = CLIP(factor*r);
                    g = CLIP(factor*g);
                    b = CLIP(factor*b);
                }
            }

			//TODO: proper treatment of out-of-gamut colors
			//float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
			float tonefactor=((r<MAXVAL ? hltonecurve[r] : CurveFactory::hlcurve (exp_scale, comp, hlrange, r) ) + \
							  (g<MAXVAL ? hltonecurve[g] : CurveFactory::hlcurve (exp_scale, comp, hlrange, g) ) + \
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
			r = tonecurve[r];
			g = tonecurve[g];
			b = tonecurve[b];
			
			//if (r<0 || g<0 || b<0) {
			//	printf("negative values row=%d col=%d  r=%f  g=%f  b=%f  \n", i,j,r,g,b);
			//}

			if (abs(sat)>0.5 || hCurveEnabled || sCurveEnabled || vCurveEnabled) {
				float h,s,v;
				rgb2hsv(r,g,b,h,s,v);
				if (sat > 0.5) {
					s = (1-(float)sat/100)*s+(float)sat/100*(1-SQR(SQR(1-s)));
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
						s = (1-satparam)*s+satparam*(1-SQR(1-s));
					} else {
						if (satparam < -0.00001)
							s *= 1+satparam;
					}

				}
				if (vCurveEnabled) {
					//shift value
					float valparam = vCurve->getVal((double)h)-0.5;
					valparam *= (1-SQR(SQR(1-s)));
					if (valparam > 0.00001) {
						v = (1-valparam)*v+valparam*(1-SQR(1-v));
					} else {
						if (valparam < -0.00001)
							v *= (1+valparam);
					}
					
				}
				hsv2rgb(h,s,v,r,g,b);
			}
			 
			//r=FCLIP(r);
			//g=FCLIP(g);
			//b=FCLIP(b);
			
            float x = (toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b) ;
            float y = (toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b) ;
            float z = (toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b) ;
			
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
			
			float x_ = 65535*Lab2xyz(fx)*D50x;
			float y_ = 65535*Lab2xyz(fy);
			float z_ = 65535*Lab2xyz(fz)*D50z;
			
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
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
			float Lin=lold->L[i][j];
			//if (Lin>0 && Lin<65535)
				lnew->L[i][j] = curve[Lin];
		}
}
		
	
void ImProcFunctions::chrominanceCurve (LabImage* lold, LabImage* lnew, LUTf & acurve, LUTf & bcurve, LUTf & satcurve) {
	
	int W = lold->W;
	int H = lold->H;

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);
	
	double wp[3][3] = {
		{wprof[0][0],wprof[0][1],wprof[0][2]},
		{wprof[1][0],wprof[1][1],wprof[1][2]},
		{wprof[2][0],wprof[2][1],wprof[2][2]}};
	
#pragma omp parallel for if (multiThread)
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
			
			float atmp = acurve[lold->a[i][j]+32768.0f]-32768.0f;
			float btmp = bcurve[lold->b[i][j]+32768.0f]-32768.0f;
			
			if (params->labCurve.saturation) {
				float chroma = sqrt(SQR(atmp)+SQR(btmp)+0.001);
				float satfactor = (satcurve[chroma+32768.0f]-32768.0f)/chroma;
				atmp *= satfactor;
				btmp *= satfactor;
			}
			
            //double real_c = 1.0;
            if (params->labCurve.avoidclip) {
				//Luv limiter
				float Y,u,v;
				Lab2Yuv(lnew->L[i][j],atmp,btmp,Y,u,v);
				//Yuv2Lab includes gamut restriction map
				Yuv2Lab(Y,u,v,lnew->L[i][j],lnew->a[i][j],lnew->b[i][j], wp);
				
            } else {
				//Luv limiter only
				lnew->a[i][j] = atmp;
				lnew->b[i][j] = btmp;
			}

        }
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	}
	

#include "cubic.cc"

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
            lnew->a[i][j] = CLIPTO(nna,-32000.0f,32000.0f);
            lnew->b[i][j] = CLIPTO(nnb,-32000.0f,32000.0f);
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
			
			PF_correct_RT(lab, lab, params->defringe.radius, params->defringe.threshold, false /*edges only*/ );
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
	
	void ImProcFunctions::lumadenoise (LabImage* lab, int** b2) {
		
		if (params->lumaDenoise.enabled && lab->W>=8 && lab->H>=8)
#ifdef _OPENMP
#pragma omp parallel
#endif
			bilateral<float, float> (lab->L, lab->L, (float**)b2, lab->W, lab->H, params->lumaDenoise.radius / scale, params->lumaDenoise.edgetolerance, multiThread);
	}
	
	void ImProcFunctions::colordenoise (LabImage* lab, int** b2) {
		
		if (params->colorDenoise.enabled && lab->W>=8 && lab->H>=8) {
#ifdef _OPENMP
#pragma omp parallel
#endif
			{
				AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(lab->W,lab->H));
				gaussHorizontal<float> (lab->a, lab->a, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);
				gaussHorizontal<float> (lab->b, lab->b, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);
				gaussVertical<float>   (lab->a, lab->a, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);
				gaussVertical<float>   (lab->b, lab->b, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);
				
				delete buffer;
			}
		}
	}
	
	void ImProcFunctions::getAutoExp  (LUTu & histogram, int histcompr, double expcomp, double clip, double& br, int& bl) {
		
		double sum = 0;
		for (int i=0; i<65536>>histcompr; i++)
			sum += histogram[i];
		
		// compute clipping points based on the original histograms (linear, without exp comp.)
		int clippable = (int)(sum * clip);
		int clipped = 0;
		int aw = (65536>>histcompr) - 1;
		while (aw>1 && histogram[aw]+clipped <= clippable) {
			clipped += histogram[aw];
			aw--;
		}
		
		clipped = 0;
		int shc = 0;
		while (shc<aw-1 && histogram[shc]+clipped <= clippable) {
			clipped += histogram[shc];
			shc++;
		}
		
		aw <<= histcompr;
		shc <<= histcompr;
		
		double corr = pow(2.0, expcomp);
		
		// black point selection is based on the linear result (yielding better visual results)
		bl = (int)(shc * corr);
		// compute the white point of the exp. compensated gamma corrected image
		double awg = (int)(CurveFactory::gamma2 (aw * corr / 65536.0) * 65536.0);
		
		// compute average intensity of the exp compensated, gamma corrected image
		double gavg = 0;
		for (int i=0; i<65536>>histcompr; i++) 
			gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;
		
		if (bl < gavg) {
			int maxaw = (gavg - bl) * 4 / 3 + bl; // dont let aw be such large that the histogram average goes above 3/4
			//double mavg = 65536.0 / (awg-bl) * (gavg - bl);
			if (awg < maxaw)
				awg = maxaw;
		}
		
		awg = CurveFactory::igamma2 ((float)(awg/65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter
		
		bl = (int)((65535*bl)/awg);
		br = log(65535.0 / (awg)) / log(2.0);
		if (br<0.0)	br = 0.0;
		if (br>10.0)	br=10.0;
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	 #include "calc_distort.h"
	 
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
	
	void ImProcFunctions::rgb2hsv (float r, float g, float b, float &h, float &s, float &v) {
		
		double var_R = r / 65535.0;
		double var_G = g / 65535.0;
		double var_B = b / 65535.0;
		
		double var_Min = MIN(MIN(var_R,var_G),var_B);
		double var_Max = MAX(MAX(var_R,var_G),var_B);
		double del_Max = var_Max - var_Min;
		v = var_Max;
		if (fabs(del_Max)<0.00001) {
			h = 0;
			s = 0;
		}
		else {
			s = del_Max/var_Max;
			
			if      ( var_R == var_Max ) h = (var_G - var_B)/del_Max; 
			else if ( var_G == var_Max ) h = 2.0 + (var_B - var_R)/del_Max; 
			else if ( var_B == var_Max ) h = 4.0 + (var_R - var_G)/del_Max; 
			h /= 6.0;
			
			if ( h < 0 )  h += 1;
			if ( h > 1 )  h -= 1;
		}
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::hsv2rgb (float h, float s, float v, float &r, float &g, float &b) {
		
		float h1 = h*6; // sector 0 to 5
		int i = floor( h1 );
		float f = h1 - i; // fractional part of h
		
		float p = v * ( 1 - s );
		float q = v * ( 1 - s * f );
		float t = v * ( 1 - s * ( 1 - f ) );
		
		float r1,g1,b1;
		
		if (i==0) {r1 = v;  g1 = t;  b1 = p;}
		if (i==1) {r1 = q;  g1 = v;  b1 = p;}
		if (i==2) {r1 = p;  g1 = v;  b1 = t;}
		if (i==3) {r1 = p;  g1 = q;  b1 = v;}
		if (i==4) {r1 = t;  g1 = p;  b1 = v;}
		if (i==5) {r1 = v;  g1 = p;  b1 = q;}
		
		r = ((r1)*65535.0);
		g = ((g1)*65535.0);
		b = ((b1)*65535.0);
	}
	
	void ImProcFunctions::xyz2srgb (float x, float y, float z, float &r, float &g, float &b) {
		
		//Transform to output color.  Standard sRGB is D65, but internal representation is D50
		//Note that it is only at this point that we should have need of clipping color data
		
		/*float x65 = d65_d50[0][0]*x + d65_d50[0][1]*y + d65_d50[0][2]*z ;
		 float y65 = d65_d50[1][0]*x + d65_d50[1][1]*y + d65_d50[1][2]*z ;
		 float z65 = d65_d50[2][0]*x + d65_d50[2][1]*y + d65_d50[2][2]*z ;
		 
		 r = sRGB_xyz[0][0]*x65 + sRGB_xyz[0][1]*y65 + sRGB_xyz[0][2]*z65;
		 g = sRGB_xyz[1][0]*x65 + sRGB_xyz[1][1]*y65 + sRGB_xyz[1][2]*z65;
		 b = sRGB_xyz[2][0]*x65 + sRGB_xyz[2][1]*y65 + sRGB_xyz[2][2]*z65;*/
		
		/*r = sRGBd65_xyz[0][0]*x + sRGBd65_xyz[0][1]*y + sRGBd65_xyz[0][2]*z ;
		 g = sRGBd65_xyz[1][0]*x + sRGBd65_xyz[1][1]*y + sRGBd65_xyz[1][2]*z ;
		 b = sRGBd65_xyz[2][0]*x + sRGBd65_xyz[2][1]*y + sRGBd65_xyz[2][2]*z ;*/
		
		r = ((sRGB_xyz[0][0]*x + sRGB_xyz[0][1]*y + sRGB_xyz[0][2]*z)) ;
		g = ((sRGB_xyz[1][0]*x + sRGB_xyz[1][1]*y + sRGB_xyz[1][2]*z)) ;
		b = ((sRGB_xyz[2][0]*x + sRGB_xyz[2][1]*y + sRGB_xyz[2][2]*z)) ;
		
	}
	
	
	void ImProcFunctions::xyz2rgb (float x, float y, float z, float &r, float &g, float &b, float rgb_xyz[3][3]) {
		
		//Transform to output color.  Standard sRGB is D65, but internal representation is D50
		//Note that it is only at this point that we should have need of clipping color data
		
		/*float x65 = d65_d50[0][0]*x + d65_d50[0][1]*y + d65_d50[0][2]*z ;
		 float y65 = d65_d50[1][0]*x + d65_d50[1][1]*y + d65_d50[1][2]*z ;
		 float z65 = d65_d50[2][0]*x + d65_d50[2][1]*y + d65_d50[2][2]*z ;
		 
		 r = sRGB_xyz[0][0]*x65 + sRGB_xyz[0][1]*y65 + sRGB_xyz[0][2]*z65;
		 g = sRGB_xyz[1][0]*x65 + sRGB_xyz[1][1]*y65 + sRGB_xyz[1][2]*z65;
		 b = sRGB_xyz[2][0]*x65 + sRGB_xyz[2][1]*y65 + sRGB_xyz[2][2]*z65;*/
		
		/*r = sRGBd65_xyz[0][0]*x + sRGBd65_xyz[0][1]*y + sRGBd65_xyz[0][2]*z ;
		 g = sRGBd65_xyz[1][0]*x + sRGBd65_xyz[1][1]*y + sRGBd65_xyz[1][2]*z ;
		 b = sRGBd65_xyz[2][0]*x + sRGBd65_xyz[2][1]*y + sRGBd65_xyz[2][2]*z ;*/
		
		
		
		r = ((rgb_xyz[0][0]*x + rgb_xyz[0][1]*y + rgb_xyz[0][2]*z)) ;
		g = ((rgb_xyz[1][0]*x + rgb_xyz[1][1]*y + rgb_xyz[1][2]*z)) ;
		b = ((rgb_xyz[2][0]*x + rgb_xyz[2][1]*y + rgb_xyz[2][2]*z)) ;
		
	}
 	
void ImProcFunctions::calcGamma (double pwr, double ts, int mode, int imax, double &gamma0, double &gamma1, double &gamma2, double &gamma3, double &gamma4, double &gamma5) {
{//from Dcraw (D.Coffin)
  int i;
  double g[6], bnd[2]={0,0}, r;

  g[0] = pwr;
  g[1] = ts;
  g[2] = g[3] = g[4] = 0;
  bnd[g[1] >= 1] = 1;
  if (g[1] && (g[1]-1)*(g[0]-1) <= 0) {
    for (i=0; i < 48; i++) {
      g[2] = (bnd[0] + bnd[1])/2;
      if (g[0]) bnd[(pow(g[2]/g[1],-g[0]) - 1)/g[0] - 1/g[2] > -1] = g[2];
      else	bnd[g[2]/exp(1-1/g[2]) < g[1]] = g[2];
    }
    g[3] = g[2] / g[1];
    if (g[0]) g[4] = g[2] * (1/g[0] - 1);
  }
  if (g[0]) g[5] = 1 / (g[1]*SQR(g[3])/2 - g[4]*(1 - g[3]) +
		(1 - pow(g[3],1+g[0]))*(1 + g[4])/(1 + g[0])) - 1;
  else      g[5] = 1 / (g[1]*SQR(g[3])/2 + 1
		- g[2] - g[3] -	g[2]*g[3]*(log(g[3]) - 1)) - 1;
  if (!mode--) {
   gamma0=g[0];gamma1=g[1];gamma2=g[2];gamma3=g[3];gamma4=g[4];gamma5=g[5];
    return;
  }
}
}
	
	void ImProcFunctions::Lab2XYZ(float L, float a, float b, float &x, float &y, float &z) {
		float fy = (0.00862069 * L) + 0.137932; // (L+16)/116
		float fx = (0.002 * a) + fy;
		float fz = fy - (0.005 * b);
		
		x = 65535*f2xyz(fx)*D50x;
		y = 65535*f2xyz(fy);
		z = 65535*f2xyz(fz)*D50z;
	}
	
	void ImProcFunctions::XYZ2Lab(float X, float Y, float Z, float &L, float &a, float &b) {
		
		float fx = (X<65535.0 ? cachef[X] : (327.68*exp(log(X/MAXVAL)/3.0 )));
		float fy = (Y<65535.0 ? cachef[Y] : (327.68*exp(log(Y/MAXVAL)/3.0 )));
		float fz = (Z<65535.0 ? cachef[Z] : (327.68*exp(log(Z/MAXVAL)/3.0 )));
		
		L = (116.0 * fy - 5242.88); //5242.88=16.0*327.68;
		a = (500.0 * (fx - fy) );
		b = (200.0 * (fy - fz) );
		
	}
	
	void ImProcFunctions::Lab2Yuv(float L, float a, float b, float &Y, float &u, float &v) {
		float fy = (0.00862069 * L/327.68) + 0.137932; // (L+16)/116
		float fx = (0.002 * a/327.68) + fy;
		float fz = fy - (0.005 * b/327.68);
		
		float X = 65535.0*f2xyz(fx)*D50x;
		Y = 65535.0*f2xyz(fy);
		float Z = 65535.0*f2xyz(fz)*D50z;
		
		u = 4.0*X/(X+15*Y+3*Z)-u0;
		v = 9.0*Y/(X+15*Y+3*Z)-v0;
		
		/*float u0 = 4*D50x/(D50x+15+3*D50z);
		 float v0 = 9/(D50x+15+3*D50z);
		 u -= u0;
		 v -= v0;*/
		
	}
	
	void ImProcFunctions::Yuv2Lab(float Yin, float u, float v, float &L, float &a, float &b, double wp[3][3]) {
		
		float u1 = u + u0;
		float v1 = v + v0;
		
		float Y = Yin;
		float X = (9*u1*Y)/(4*v1*D50x); 
		float Z = (12 - 3*u1 - 20*v1)*Y/(4*v1*D50z);
		
		gamutmap(X,Y,Z,wp);
		
		float fx = (X<65535.0 ? cachef[X] : (327.68*exp(log(X/MAXVAL)/3.0 )));
		float fy = (Y<65535.0 ? cachef[Y] : (327.68*exp(log(Y/MAXVAL)/3.0 )));
		float fz = (Z<65535.0 ? cachef[Z] : (327.68*exp(log(Z/MAXVAL)/3.0 )));
		
		L = (116.0 * fy - 5242.88); //5242.88=16.0*327.68;
		a = (500.0 * (fx - fy) );
		b = (200.0 * (fy - fz) );
		
	}
	
#include "gamutbdy.cc"
}

#undef eps_max
#undef kappa
