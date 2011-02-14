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
#include <glibmm.h>
#include <iccstore.h>
#include <iccmatrices.h>
#include <impulse_denoise.h>

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
	
#define eps_max 580.40756 //(MAXVAL* 216.0f/24389.0);
#define kappa	903.29630 //24389.0/27.0;
	
	
//%%%%%%%%%%%%
#define epsilon 0.00885645 //216/24389
#define kappainv 0.00110706 //inverse of kappa
#define kapeps 8 // kappa*epsilon
#define Lab2xyz(f) (( (g=f*f*f) > epsilon) ? g : (116*f-16)*kappainv)
//%%%%%%%%%%%%


extern const Settings* settings;

LUTf ImProcFunctions::cachef = 0;
/*float* ImProcFunctions::cacheL;
float* ImProcFunctions::cachea;
float* ImProcFunctions::cacheb;
float* ImProcFunctions::xcache;
float* ImProcFunctions::ycache;
float* ImProcFunctions::zcache;*/
//unsigned short ImProcFunctions::gamma2curve[65536];
LUTf ImProcFunctions::gamma2curve = 0;

void ImProcFunctions::initCache () {

    int maxindex = 65536;
	cachef(maxindex,0);
    //cacheL = new float[maxindex];
    //cachea = new float[maxindex];
    //cacheb = new float[maxindex];
	gamma2curve(maxindex);

    for (int i=0; i<maxindex; i++) {
        if (i>eps_max) {
			cachef[i] = 327.68*( exp(1.0/3.0 * log((double)i / MAXVAL) ));
        }
        else {
			cachef[i] = 327.68*((kappa*i/MAXVAL+16.0)/116.0);
        }
	}

	for (int i=0; i<65536; i++) {
		gamma2curve[i] = (CurveFactory::gamma2(i/65535.0) * 65535.0);
	}
}

void ImProcFunctions::cleanupCache () {

	//delete [] cachef;
	//delete [] gamma2curve;
}

ImProcFunctions::~ImProcFunctions () {

	if (monitorTransform!=NULL)
		cmsDeleteTransform (monitorTransform);
}

void ImProcFunctions::setScale (double iscale) {
	scale = iscale;
}

void ImProcFunctions::firstAnalysis_ (Imagefloat* original, Glib::ustring wprofile, unsigned int* histogram, int row_from, int row_to) {

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
	cmsHPROFILE monitor = iccStore->getProfile ("file:"+settings->monitorProfile);
	if (monitor) {
        cmsHPROFILE iprof = iccStore->getXYZProfile ();       
        lcmsMutex->lock ();
		monitorTransform = cmsCreateTransform (iprof, TYPE_RGB_FLT, monitor, TYPE_RGB_8, settings->colorimetricIntent, cmsFLAGS_NOOPTIMIZE);
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
			firstAnalysis_ (original, wprofile, hist[tid], tid*blk, (tid+1)*blk);
		else
			firstAnalysis_ (original, wprofile, hist[tid], tid*blk, H);
    }
#else
    firstAnalysis_ (original, wprofile, hist[0], 0, original->height);
#endif
 
	histogram.clear();
    for (int i=0; i<65536; i++)
    	for (int j=0; j<T; j++)
    		histogram[i] += hist[j][i];

    for (int i=0; i<T; i++)
    	delete [] hist[i];
    delete [] hist;

}

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

//    int mapval;
//    double factor;
    int tW = working->width;
    int tH = working->height;
//    int r, g, b;
//	float h, s, v;
//	float satparam,valparam;
//	int hue, hueband, hueres, nbrband;
	double pi = M_PI;
	
	float* cossq = new float [8192];
	for (int i=0; i<8192; i++) 
		cossq[i] = SQR(cos(pi*(float)i/16384));
	
	
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
                double mapval = shmap->map[i][j];
                double factor = 1.0;
                
                if (processSH) {
                    if (mapval > h_th) 
                        factor = (h_th + (100.0 - params->sh.highlights) * (mapval - h_th) / 100.0) / mapval; 
                    else if (mapval < s_th) 
                        factor = (s_th - (100.0 - params->sh.shadows) * (s_th - mapval) / 100.0) / mapval; 
                }
                if (processLCE) {
                    double sub = lceamount*(mapval-factor*(r*lumimul[0] + g*lumimul[1] + b*lumimul[2]));
                    r = (factor*r-sub);
                    g = (factor*g-sub);
                    b = (factor*b-sub);
                }
                else {
                    r = (factor*r);
                    g = (factor*g);
                    b = (factor*b);
                }
            }

			//TODO: proper treatment of out-of-gamut colors
			//float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
			float tonefactor=(hltonecurve[r] + hltonecurve[g] + hltonecurve[b])/3;

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

			if (abs(sat)>0.5 || params->hsvequalizer.enabled) {
				float h,s,v;
				rgb2hsv(r,g,b,h,s,v);
				if (sat > 0.5) {
					s = (1-(float)sat/100)*s+(float)sat/100*(1-SQR(SQR(1-s)));
				} else {
					if (sat < -0.5)
						s *= 1+(float)sat/100;	
				}
				//HSV equalizer
				if (params->hsvequalizer.enabled) {
					int hue = (int)(65535*h);
					int hueres = hue & 8191;//location of hue within a band
					int hueband = (hue-hueres) >> 13;//divides hue range into 8 bands
					int nbrband = (hueband+1)&7;
					
					//shift hue
					h = fmod(h + 0.0025*(params->hsvequalizer.hue[hueband] * cossq[hueres] + params->hsvequalizer.hue[nbrband] * (1-cossq[hueres])),1);
					if (h<0) h +=1;
					hue = (int)round(65535*h);
					hueres = hue & 8191;//location of hue within a band
					hueband = (hue-hueres) >> 13;//divides hue range into 8 bands
					nbrband = (hueband+1)&7;

					//change saturation
					float satparam = 0.01*(params->hsvequalizer.sat[hueband] * cossq[hueres] + params->hsvequalizer.sat[nbrband] * (1-cossq[hueres]));
					if (satparam > 0.00001) {
						s = (1-satparam)*s+satparam*(1-SQR(1-s));
					} else {
						if (satparam < -0.00001)
							s *= 1+satparam;	
					}
					
					//change value
					float valparam = 0.005*(params->hsvequalizer.val[hueband] * cossq[hueres] + params->hsvequalizer.val[nbrband] * (1-cossq[hueres]));
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
			
			/*lab->L[i][j] = CLIP2(116.0 * (CurveFactory::flinterp(cachef,y)) - 5242.88); //5242.88=16.0*327.68;
            lab->a[i][j] = CLIPS(500.0 * (((CurveFactory::flinterp(cachef,x) - CurveFactory::flinterp(cachef,y)) ) ));
            lab->b[i][j] = CLIPS(200.0 * (((CurveFactory::flinterp(cachef,y) - CurveFactory::flinterp(cachef,z)) ) ));*/
			
			x = (x<65535.0 ? cachef[x] : (327.68*exp(log(x/MAXVAL)/3.0 )));
			y = (y<65535.0 ? cachef[y] : (327.68*exp(log(y/MAXVAL)/3.0 )));
			z = (z<65535.0 ? cachef[z] : (327.68*exp(log(z/MAXVAL)/3.0 )));

			lab->L[i][j] = (116.0 * y - 5242.88); //5242.88=16.0*327.68;
            lab->a[i][j] = (500.0 * (x - y) );
            lab->b[i][j] = (200.0 * (y - z) );

			//float L1 = lab->L[i][j];//for testing
			//float a1 = lab->a[i][j];
			//float b1 = lab->b[i][j];
			//float xxx=1;
			
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
		
	
void ImProcFunctions::chrominanceCurve (LabImage* lold, LabImage* lnew, LUTf & acurve, LUTf & bcurve) {
	
	int W = lold->W;
	int H = lold->H;
	for (int i=0; i<H; i++)
		for (int j=0; j<W; j++) {
			float ain=lold->a[i][j];
			float bin=lold->b[i][j];
			//if (fabs(ain)<32767 && fabs(bin)<32767)
			lnew->a[i][j] = acurve[ain+32768.0f]-32768.0;
			lnew->b[i][j] = bcurve[bin+32768.0f]-32768.0;
		}
}
	
#include "cubic.cc"

void ImProcFunctions::colorCurve (LabImage* lold, LabImage* lnew) {

    LUT<double> cmultiplier(181021);

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

    //delete [] cmultiplier;
}
	
	void ImProcFunctions::impulsedenoise (LabImage* lab) {
		
		if (params->impulseDenoise.enabled && lab->W>=8 && lab->H>=8)
			
			impulse_nr (lab->L, lab->W, lab->H, (float)params->impulseDenoise.thresh/20.0 );
	}
	
	void ImProcFunctions::defringe (LabImage* lab) {
		
		if (params->defringe.enabled && lab->W>=8 && lab->H>=8)
			
			PF_correct_RT(lab, lab, params->defringe.radius, params->defringe.threshold, false /*edges only*/ );
	}
	
	void ImProcFunctions::dirpyrdenoise (LabImage* lab) {
		
		if (params->dirpyrDenoise.enabled && lab->W>=8 && lab->H>=8)
			
			dirpyrLab_denoise(lab, lab, params->dirpyrDenoise.luma, params->dirpyrDenoise.chroma, params->dirpyrDenoise.gamma/3.0 );
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
	
}

#undef eps_max
#undef kappa
