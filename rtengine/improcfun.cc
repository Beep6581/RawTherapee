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

namespace rtengine {

using namespace procparams;

#undef MAX
#undef MIN
#undef MAXVAL
#undef CLIP
#undef CLIPS
#undef CLIPC
#undef CLIPTO
#undef CLIPTOC
#undef THREAD_PRIORITY_NORMAL

#define MAXVAL  0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)
#define CLIPS(a) ((a)>-32768?((a)<32767?(a):32767):-32768)
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIPTOC(a,b,c,d) ((a)>=(b)?((a)<=(c)?(a):((c),d=true)):((b),d=true))

extern const Settings* settings;

//
// STRUCTURES FOR THE SHADOW MAP AND LAB SPACE IMAGE
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LabImage::LabImage (int w, int h) : W(w), H(h), fromImage(false) {

    L = new unsigned short*[H];
    for (int i=0; i<H; i++)
        L[i] = new unsigned short[W];
        
    a = new short*[H];
    for (int i=0; i<H; i++)
        a[i] = new short[W];
        
    b = new short*[H];
    for (int i=0; i<H; i++)
        b[i] = new short[W];
}

LabImage::LabImage (Image16* im) {

    W = im->width;
    H = im->height;
    L = im->r;
    a = (short**) im->g;
    b = (short**) im->b;
    fromImage = true;
}

LabImage::~LabImage () {

    if (!fromImage) {
        for (int i=0; i<H; i++) {
            delete [] L[i];
            delete [] a[i];
            delete [] b[i];
        }
        delete [] L;
        delete [] a;
        delete [] b;
    }
}

//
// IMAGE PROCESSING FUNCTIONS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~

#undef MAX
#undef MIN
#undef CLIP
#undef CMAXVAL
#undef CLIPV
#undef ABS

#define CMAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPV(a,b,c) ((a)<(b)?(b):((a)>(c)?(c):(a)))

#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))

#define MAXL 65535
#define ABS(a) ((a)<0?-(a):(a))


int* ImProcFunctions::cacheL;
int* ImProcFunctions::cachea;
int* ImProcFunctions::cacheb;
int* ImProcFunctions::xcache;
int* ImProcFunctions::ycache;
int* ImProcFunctions::zcache;
unsigned short ImProcFunctions::gamma2curve[65536];

/*const int c00 = (int) (32768.0 * 0.412453 / 0.950456);
const int c01 = (int) (32768.0 * 0.357580 / 0.950456);
const int c02 = (int) (32768.0 * 0.180423 / 0.950456);
const int c10 = (int) (32768.0 * 0.212671);
const int c11 = (int) (32768.0 * 0.715160);
const int c12 = (int) (32768.0 * 0.072169);
const int c20 = (int) (32768.0 * 0.019334 / 1.088754);
const int c21 = (int) (32768.0 * 0.119193 / 1.088754);
const int c22 = (int) (32768.0 * 0.950227 / 1.088754);
*/

void ImProcFunctions::initCache () {

    int maxindex = 2*65536;
    cacheL = new int[maxindex];
    cachea = new int[maxindex];
    cacheb = new int[maxindex];

    int threshold = (int)(0.008856*CMAXVAL);
    for (int i=0; i<maxindex; i++)
        if (i>threshold) {
            cacheL[i] = (int)round(655.35 * (116.0 * exp(1.0/3.0 * log((double)i / CMAXVAL)) - 16.0));
            cachea[i] = (int)round(32768.0 * 500.0 * exp(1.0/3.0 * log((double)i / CMAXVAL)));
            cacheb[i] = (int)round(32768.0 * 200.0 * exp(1.0/3.0 * log((double)i / CMAXVAL)));
        }
        else {
            cacheL[i] = (int)round(9033.0 * (double)i / 1000.0); // assuming CMAXVAL = 65535
            cachea[i] = (int)round(32768.0 * 500.0 * (7.787*i/CMAXVAL+16.0/116.0));
            cacheb[i] = (int)round(32768.0 * 200.0 * (7.787*i/CMAXVAL+16.0/116.0));
        }

    double fY;
    ycache = new int[0x10000];
    for (int i=0; i<0x10000; i++)
        ycache[i] = (int)round(65536.0 * ((fY=((double)i/655.35+16)/116) > 2.0689655172413793e-1 ? fY*fY*fY : 1.107056459879453852e-3*(double)i/655.35));
    for (int i=0; i<0x10000; i++)
        ycache[i] = CLIP(ycache[i]);
    xcache = new int[369621];
    for (int i=-141556; i<228064; i++)
        xcache[i+141556] = (int)round(65536.0 * (i > 15728 ? ((double)i/76021)*((double)i/76021)*((double)i/76021)*0.96422 : (1.2841854934601665e-1*(double)i/76021-1.7712903358071262e-2)*0.96422));
    for (int i=0; i<369620; i++)
        xcache[i] = CLIP(xcache[i]);
    zcache = new int[825747];
    for (int i=-369619; i<456127; i++)
        zcache[i+369619] = (int)round(65536.0 * (i > 15728 ? ((double)i/76021)*((double)i/76021)*((double)i/76021)*0.82521 : (1.2841854934601665e-1*(double)i/76021-1.7712903358071262e-2)*0.82521));
    for (int i=0; i<825747; i++)
        zcache[i] = CLIP(zcache[i]);

	for (int i=0; i<65536; i++) {
		int g = (int)(CurveFactory::gamma2(i/65535.0) * 65535.0);
		gamma2curve[i] = CLIP(g);
	}
}

void ImProcFunctions::release () {

	if (monitorTransform!=NULL)
		cmsDeleteTransform (monitorTransform);
	monitorTransform = NULL;
}

void ImProcFunctions::firstAnalysis_ (Image16* original, Glib::ustring wprofile, int* histogram, int* chroma_radius, int row_from, int row_to) {

	TMatrix wprof = iccStore.workingSpaceMatrix (wprofile);
    int toxyz[3][3];
    toxyz[0][0] = round(32768.0 * wprof[0][0] / 0.96422); 
    toxyz[1][0] = round(32768.0 * wprof[1][0] / 0.96422); 
    toxyz[2][0] = round(32768.0 * wprof[2][0] / 0.96422); 
    toxyz[0][1] = round(32768.0 * wprof[0][1]); 
    toxyz[1][1] = round(32768.0 * wprof[1][1]); 
    toxyz[2][1] = round(32768.0 * wprof[2][1]); 
    toxyz[0][2] = round(32768.0 * wprof[0][2] / 0.82521); 
    toxyz[1][2] = round(32768.0 * wprof[1][2] / 0.82521); 
    toxyz[2][2] = round(32768.0 * wprof[2][2] / 0.82521); 

	lumimul[0] = wprof[0][1];
	lumimul[1] = wprof[1][1];
	lumimul[2] = wprof[2][1];
	
    int W = original->width;
    int cradius = 1;
    for (int i=row_from; i<row_to; i++) {
        for (int j=0; j<W; j++) {
      
            int r = original->r[i][j];
            int g = original->g[i][j];
            int b = original->b[i][j];

            int x = (toxyz[0][0] * r + toxyz[1][0] * g + toxyz[2][0] * b) >> 15;
            int y = (toxyz[0][1] * r + toxyz[1][1] * g + toxyz[2][1] * b) >> 15;
            int z = (toxyz[0][2] * r + toxyz[1][2] * g + toxyz[2][2] * b) >> 15;

            x = CLIPTO(x,0,2*65536-1);
            y = CLIPTO(y,0,2*65536-1);
            z = CLIPTO(z,0,2*65536-1);

            int oa = cachea[x] - cachea[y];
            int ob = cacheb[y] - cacheb[z];

            if (oa<0) oa = -oa;
            if (ob<0) ob = -ob;

            if (oa > cradius)
                cradius = oa;
            if (ob > cradius)
                cradius = ob;

            if (histogram) {
                int hval = CLIP(y); //(306 * original->r[i][j] + 601 * original->g[i][j] + 117 * original->b[i][j]) >> 10;
                histogram[hval]++;
            }
        }
    }
    *chroma_radius = cradius;
}

void ImProcFunctions::firstAnalysis (Image16* original, const ProcParams* params, int* histogram, double gamma) {

    int cr1, cr2;
    int* hist1 = new int[65536]; memset (hist1, 0, 65536*sizeof(int));
    int* hist2 = new int[65536]; memset (hist2, 0, 65536*sizeof(int));
    
    int H = original->height;

    Glib::ustring wprofile = params->icm.working;
	if (monitorTransform)
		cmsDeleteTransform (monitorTransform);
	monitorTransform = NULL;
	cmsHPROFILE monitor = iccStore.getProfile ("file:"+settings->monitorProfile);
	if (monitor) {
        cmsHPROFILE iprof = iccStore.getXYZProfile ();       
		cmsHPROFILE oprof = iccStore.getProfile (params->icm.output);
		if (!oprof)
			oprof = iccStore.getsRGBProfile ();
        lcmsMutex->lock ();
		monitorTransform = cmsCreateTransform (iprof, TYPE_RGB_16, monitor, TYPE_RGB_8, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
	}
    
    if (settings->dualThreadEnabled) {
        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::firstAnalysis_), original, wprofile, hist1, &cr1, 0, H/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::firstAnalysis_), original, wprofile, hist2, &cr2, H/2, H), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else {
      firstAnalysis_ (original, wprofile, hist1, &cr1, 0, H/2);
      firstAnalysis_ (original, wprofile, hist2, &cr2, H/2, H);
    }
    
    if (cr1<cr2) 
        chroma_radius = cr2;
    else
        chroma_radius = cr1;
 
    for (int i=0; i<65536; i++)
        histogram[i] = hist1[i]+hist2[i];

    chroma_scale = 32768*32768 / (3*chroma_radius);
    delete [] hist1;
    delete [] hist2;
}

void ImProcFunctions::rgbProc (Image16* working, LabImage* lab, const ProcParams* params, int* tonecurve, SHMap* shmap) {

    if (settings->dualThreadEnabled) {

        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::rgbProc_), working, lab, params, tonecurve, shmap, 0, working->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::rgbProc_), working, lab, params, tonecurve, shmap, working->height/2, working->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else
        rgbProc_ (working, lab, params, tonecurve, shmap, 0, working->height);
}

void ImProcFunctions::rgbProc_ (Image16* working, LabImage* lab, const ProcParams* params, int* tonecurve, SHMap* shmap, int row_from, int row_to) {

    int r, g, b;

    int h_th, s_th;
    if (shmap) {
        h_th = shmap->max - params->sh.htonalwidth * (shmap->max - shmap->avg) / 100;
        s_th = params->sh.stonalwidth * (shmap->avg - shmap->min) / 100;
    }

    bool processSH  = params->sh.enabled && shmap!=NULL && (params->sh.highlights>0 || params->sh.shadows>0);
    bool processLCE = params->sh.enabled && shmap!=NULL && params->sh.localcontrast>0;
    double lceamount = params->sh.localcontrast / 200.0;

    TMatrix wprof = iccStore.workingSpaceMatrix (params->icm.working);
    int toxyz[3][3] = {
        floor(32768.0 * wprof[0][0] / 0.96422), 
        floor(32768.0 * wprof[0][1]),
        floor(32768.0 * wprof[0][2] / 0.82521),
        floor(32768.0 * wprof[1][0] / 0.96422),
        floor(32768.0 * wprof[1][1]),
        floor(32768.0 * wprof[1][2] / 0.82521),
        floor(32768.0 * wprof[2][0] / 0.96422),
        floor(32768.0 * wprof[2][1]),
        floor(32768.0 * wprof[2][2] / 0.82521)};

    bool mixchannels = params->chmixer.red[0]!=100 || params->chmixer.red[1]!=0 || params->chmixer.red[2]!=0 || params->chmixer.green[0]!=0 || params->chmixer.green[1]!=100 || params->chmixer.green[2]!=0 || params->chmixer.blue[0]!=0 || params->chmixer.blue[1]!=0 || params->chmixer.blue[2]!=100;

    int mapval;
    double factor;
    int tW = working->width;
    for (int i=row_from; i<row_to; i++) {
        for (int j=0; j<tW; j++) {

            r = working->r[i][j];
            g = working->g[i][j];
            b = working->b[i][j];

            if (mixchannels) {
                int newr = (r*params->chmixer.red[0]   + g*params->chmixer.red[1]   + b*params->chmixer.red[2]) / 100;
                int newg = (r*params->chmixer.green[0] + g*params->chmixer.green[1] + b*params->chmixer.green[2]) / 100;
                int newb = (r*params->chmixer.blue[0]  + g*params->chmixer.blue[1]  + b*params->chmixer.blue[2]) / 100;
                r = CLIP(newr);
                g = CLIP(newg);
                b = CLIP(newb);
            }

            if (processSH || processLCE) {
                mapval = shmap->map[i][j];
                factor = 1.0;
                
                if (processSH) {
                    if (mapval > h_th) 
                        factor = (h_th + (100.0 - params->sh.highlights) * (mapval - h_th) / 100.0) / mapval; 
                    else if (mapval < s_th) 
                        factor = (s_th - (100.0 - params->sh.shadows) * (s_th - mapval) / 100.0) / mapval; 
                }
                if (processLCE) {
                    double sub = lceamount*(mapval-factor*(r*lumimul[0] + g*lumimul[1] + b*lumimul[2]));
                    r = CLIP((int)(factor*r-sub));
                    g = CLIP((int)(factor*g-sub));
                    b = CLIP((int)(factor*b-sub));
                }
                else {
                if (i==100 && j==3500)
                    printf ("r=%d, %d, fact=%g, mapval=%d, %d\n", r, (int)(factor*r), factor, mapval, shmap->map[i][j]);
                    r = CLIP((int)(factor*r));
                    g = CLIP((int)(factor*g));
                    b = CLIP((int)(factor*b));
                }
            }
            r = tonecurve[r];
            g = tonecurve[g];
            b = tonecurve[b];

//            int x = (14219 * r + 12328 * g +  6220 * b) >> 15;
//            int y = ( 6968 * r + 23434 * g +  2365 * b) >> 15;
//            int z = (  582 * r +  3587 * g + 28598 * b) >> 15;
            int x = (toxyz[0][0] * r + toxyz[1][0] * g + toxyz[2][0] * b) >> 15;
            int y = (toxyz[0][1] * r + toxyz[1][1] * g + toxyz[2][1] * b) >> 15;
            int z = (toxyz[0][2] * r + toxyz[1][2] * g + toxyz[2][2] * b) >> 15;

            x = CLIPTO(x,0,2*65536-1);
            y = CLIPTO(y,0,2*65536-1);
            z = CLIPTO(z,0,2*65536-1);

            int L = cacheL[y];
            lab->L[i][j] = L;
            lab->a[i][j] = CLIPC(((cachea[x] - cachea[y]) * chroma_scale) >> 15);
            lab->b[i][j] = CLIPC(((cacheb[y] - cacheb[z]) * chroma_scale) >> 15);
        }
    }
 }

void ImProcFunctions::luminanceCurve (LabImage* lold, LabImage* lnew, int* curve, int row_from, int row_to) {

    int W = lold->W;
    int H = lold->H;
    for (int i=row_from; i<row_to; i++)
        for (int j=0; j<W; j++)
            lnew->L[i][j] = curve[lold->L[i][j]];
}

#include "cubic.cc"

void ImProcFunctions::colorCurve (LabImage* lold, LabImage* lnew, const ProcParams* params) {

    double* cmultiplier = new double [181021];

    double boost_a = (params->colorBoost.amount + 100.0) / 100.0;
    double boost_b = (params->colorBoost.amount + 100.0) / 100.0;

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

    if (params->colorBoost.enable_saturationlimiter && c>1) {
        // re-generate color multiplier lookup table
        double d = params->colorBoost.saturationlimit * chroma_scale  / 3.0;
        double alpha = 0.5;
        double threshold1 = alpha * d;
        double threshold2 = c*d*(alpha+1.0) - d;
        for (int i=0; i<=181020; i++) { // lookup table stores multipliers with a 0.25 chrominance resolution
            double chrominance = (double)i/4;
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
    
    if (settings->dualThreadEnabled) {

        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::colorCurve_), lold, lnew, params, 0, lnew->H/2, cmultiplier), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::colorCurve_), lold, lnew, params, lnew->H/2, lnew->H, cmultiplier), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else
        colorCurve_ (lold, lnew, params, 0, lnew->H, cmultiplier);
        
    delete [] cmultiplier;
}

void ImProcFunctions::colorCurve_ (LabImage* lold, LabImage* lnew, const ProcParams* params, int row_from, int row_to, double* cmultiplier) {

    double boost_a = (params->colorBoost.amount + 100.0) / 100.0;
    double boost_b = (params->colorBoost.amount + 100.0) / 100.0;

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

    int nna, nnb;
    double shift_a = params->colorShift.a * chroma_scale, shift_b = params->colorShift.b * chroma_scale;

    short** oa = lold->a;
    short** ob = lold->b;

    for (int i=row_from; i<row_to; i++)
        for (int j=0; j<lold->W; j++) {

            double wanted_c = c;
            if (params->colorBoost.enable_saturationlimiter && c>1) {
                int chroma = (int)(4.0 * sqrt((oa[i][j]+shift_a)*(oa[i][j]+shift_a) + (ob[i][j]+shift_b)*(ob[i][j]+shift_b)));
                wanted_c = cmultiplier [MIN(chroma,181020)];
            }

            double real_c = wanted_c;
            if (wanted_c >= 1.0 && params->colorBoost.avoidclip) {
                double cclip = 100000;
                double cr = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)/chroma_scale*amul, (double)(ob[i][j]+shift_b)/chroma_scale*bmul, 3.079935, -1.5371515, -0.54278342);
                double cg = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)/chroma_scale*amul, (double)(ob[i][j]+shift_b)/chroma_scale*bmul, -0.92123418, 1.87599, 0.04524418);
                double cb = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)/chroma_scale*amul, (double)(ob[i][j]+shift_b)/chroma_scale*bmul, 0.052889682, -0.20404134, 1.15115166);
                if (cr>1.0 && cr<cclip) cclip = cr;
                if (cg>1.0 && cg<cclip) cclip = cg;
                if (cb>1.0 && cb<cclip) cclip = cb;
                if (cclip<100000) {
                    real_c = -cclip + 2.0*cclip / (1.0+exp(-2.0*wanted_c/cclip));
                    if (real_c<1.0)
                        real_c = 1.0;
                }
            }
            
            nna = (int)((oa[i][j]+shift_a) * real_c * amul);
            nnb = (int)((ob[i][j]+shift_b) * real_c * bmul);
            lnew->a[i][j] = CLIPV(nna,-32000,32000);
            lnew->b[i][j] = CLIPV(nnb,-32000,32000);
        }
}

void blur (float** src, float** dst, int W, int H, int r) {

    float** tmpI = new float*[H];
    for (int i=0; i<H; i++) 
        tmpI[i] = new float[W];

    for (int i=0; i<H; i++) 
        for (int j=0; j<W; j++) {
            if (i<r || i>=H-r)
                tmpI[i][j] = src[i][j];
            else {
                int num = 0;
                float sum = 0.0;
                for (int x=-r; x<=r; x++)
                    for (int y=-r; y<=r; y++)
                        if (x*x+y*y<=r*r) {
                            sum += src[i+x][j+y];
                            num++;
                        }
                tmpI[i][j] = sum / num;
            }
        }
    for (int i=0; i<H; i++) 
        for (int j=0; j<W; j++) 
            dst[i][j] = tmpI[i][j];

    for (int i=0; i<H; i++)
        delete [] tmpI[i];
    delete [] tmpI;
}

void ImProcFunctions::damping_ (float** aI, unsigned short** aO, float damping, int W, int rowfrom, int rowto) {

    for (int i=rowfrom; i<rowto; i++)
        for (int j=0; j<W; j++) {
            float I = aI[i][j];
            float O = (float)aO[i][j];
            if (O==0.0 || I==0.0) {
                aI[i][j] = 0.0;
                continue;
            }
            float U = -(O * log(I/O) - I + O) * 2.0 / (damping*damping);
            U = MIN(U,1.0);
            U = U*U*U*U*(5.0-U*4.0);
            aI[i][j] = (O - I) / I * U + 1.0;
        }
}

void ImProcFunctions::deconvsharpening (LabImage* lab, const ProcParams* params, double scale, unsigned short** b2) {

    if (params->sharpening.enabled==false || params->sharpening.deconvamount<1)
        return;

    int W = lab->W, H = lab->H;    

    float** tmpI = new float*[H];
    for (int i=0; i<H; i++) {
        tmpI[i] = new float[W];
        for (int j=0; j<W; j++)
            tmpI[i][j] = (float)lab->L[i][j];
    }

    float** tmp = (float**)b2;
    
    AlignedBuffer<double>* buffer1 = new AlignedBuffer<double> (MAX(W,H)*5);
    AlignedBuffer<double>* buffer2 = new AlignedBuffer<double> (MAX(W,H)*5);
    
    float damping = params->sharpening.deconvdamping / 5.0;
    bool needdamp = params->sharpening.deconvdamping > 0;
    for (int k=0; k<params->sharpening.deconviter; k++) {
    
//        gaussHorizontal<float> (tmpI, tmp, buffer1, W, 0, H, params->sharpening.deconvradius / scale);
//        gaussVertical<float>   (tmp, tmp, buffer1, H, 0, W, params->sharpening.deconvradius / scale);
        // apply blur function (gaussian blur)
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_float), tmpI, tmp, buffer1, W, 0, H/2, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_float), tmpI, tmp, buffer2, W, H/2, H, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            gaussHorizontal_float (tmpI, tmp, buffer1, W, 0, H, params->sharpening.deconvradius / scale);
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_float), tmp, tmp, buffer1, H, 0, W/2, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_float), tmp, tmp, buffer2, H, W/2, W, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            gaussVertical_float (tmp, tmp, buffer1, H, 0, W, params->sharpening.deconvradius / scale);

//        blur (tmpI, tmp, W, H, params->sharpening.radius / scale);
        if (!needdamp) {
            for (int i=0; i<H; i++)
                for (int j=0; j<W; j++) 
                    if (tmp[i][j]>0) 
                        tmp[i][j] = (float)lab->L[i][j] / tmp[i][j];
        }
        else {
            if (settings->dualThreadEnabled) {
                Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::damping_), tmp, lab->L, damping, W, 0, H/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
                Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::damping_), tmp, lab->L, damping, W, H/2, H), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
                thread1->join ();
                thread2->join ();
            }
            else 
                damping_ (tmp, lab->L, damping, W, 0, H);
        }
//        gaussHorizontal<float> (tmp, tmp, buffer1, W, 0, H, params->sharpening.deconvradius / scale);
//        gaussVertical<float>   (tmp, tmp, buffer1, H, 0, W, params->sharpening.deconvradius / scale);
        // apply blur function (gaussian blur)
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_float), tmp, tmp, buffer1, W, 0, H/2, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_float), tmp, tmp, buffer2, W, H/2, H, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            gaussHorizontal_float (tmp, tmp, buffer1, W, 0, H, params->sharpening.deconvradius / scale);
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_float), tmp, tmp, buffer1, H, 0, W/2, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_float), tmp, tmp, buffer2, H, W/2, W, params->sharpening.deconvradius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            gaussVertical_float (tmp, tmp, buffer1, H, 0, W, params->sharpening.deconvradius / scale);

//        blur (tmp, tmp, W, H, params->sharpening.radius / scale);

        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++)
                tmpI[i][j] = tmpI[i][j] * tmp[i][j];
    }
   
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++)
            lab->L[i][j] = lab->L[i][j]*(100-params->sharpening.deconvamount) / 100 + (int)CLIP(tmpI[i][j])*params->sharpening.deconvamount / 100;
   
    for (int i=0; i<H; i++)
        delete [] tmpI[i];
    delete [] tmpI;
}

void ImProcFunctions::sharpening (LabImage* lab, const ProcParams* params, double scale, unsigned short** b2) {

    if (params->sharpening.method=="rld") {
        deconvsharpening (lab, params, scale, b2);
        return;
    }

    if (params->sharpening.enabled==false || params->sharpening.amount<1 || lab->W<8 || lab->H<8)
        return;

    int W = lab->W, H = lab->H;    
    unsigned short** b3;
    if (params->sharpening.edgesonly==false) {
  
        AlignedBuffer<double>* buffer1 = new AlignedBuffer<double> (MAX(W,H)*5);
        AlignedBuffer<double>* buffer2 = new AlignedBuffer<double> (MAX(W,H)*5);
  
	MyTime t1, t2, t3;
	t1.set ();
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_unsigned), lab->L, b2, buffer1, W, 0, H/2, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_unsigned), lab->L, b2, buffer2, W, H/2, H, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            gaussHorizontal_unsigned (lab->L, b2, buffer1, W, 0, H, params->sharpening.radius / scale);

	t2.set ();
//        printf ("Horizontal: %d\n", t2.etime (t1));

        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_unsigned), b2, b2, buffer1, H, 0, W/2, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_unsigned), b2, b2, buffer2, H, W/2, W, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            gaussVertical_unsigned (b2, b2, buffer1, H, 0, W, params->sharpening.radius / scale);

	t3.set ();
//        printf ("Vertical: %d\n", t3.etime (t2));

        delete buffer1;
        delete buffer2;
    }
    else {
        b3 = new unsigned short*[H];
        for (int i=0; i<H; i++)
            b3[i] = new unsigned short[W];

        Dim dim1 (W, H, 0, H/2);
        Dim dim2 (W, H, H/2, H);

        AlignedBuffer<double>* buffer1 = new AlignedBuffer<double> (MAX(W,H)*5);
        AlignedBuffer<double>* buffer2 = new AlignedBuffer<double> (MAX(W,H)*5);

        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(bilateral_unsigned), lab->L, b3, b2, dim1, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(bilateral_unsigned), lab->L, b3, b2, dim2, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
            thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_unsigned), b3, b2, buffer1, W, 0, H/2, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_unsigned), b3, b2, buffer2, W, H/2, H, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
            thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_unsigned), b2, b2, buffer1, H, 0, W/2, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_unsigned), b2, b2, buffer2, H, W/2, W, params->sharpening.radius / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else {
            bilateral_unsigned (lab->L, (unsigned short**)b2, b3, dim1, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance);
            bilateral_unsigned (lab->L, (unsigned short**)b2, b3, dim2, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance);
            gaussHorizontal_unsigned (b2, b2, buffer1, W, 0, H, params->sharpening.radius / scale);
            gaussVertical_unsigned   (b2, b2, buffer1, H, 0, W, params->sharpening.radius / scale);
        }

        delete buffer1;
        delete buffer2;
    }
    unsigned short** base = lab->L;
    if (params->sharpening.edgesonly) 
        base = b3;
        
    if (params->sharpening.halocontrol==false) {
        for (int i=0; i<H; i++) 
            for (int j=0; j<W; j++) {
                int diff = base[i][j] - b2[i][j];
                if (ABS(diff)>params->sharpening.threshold) {
                    int val = lab->L[i][j] + params->sharpening.amount * diff / 100;
                    lab->L[i][j] = CLIP(val);
                }
            }
    }
    else {
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::sharpenHaloCtrl), lab, params, b2, base, W, 2, H/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::sharpenHaloCtrl), lab, params, b2, base, W, H/2, H-2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else
            sharpenHaloCtrl (lab, params, b2, base, W, 2, H-2);
    }

    if (params->sharpening.edgesonly) {
        for (int i=0; i<H; i++)
            delete [] b3[i];
        delete [] b3;
    }
}

void ImProcFunctions::sharpenHaloCtrl (LabImage* lab, const ProcParams* params, unsigned short** blurmap, unsigned short** base, int W, int row_from, int row_to) {

    int scale = 100 * (100-params->sharpening.halocontrol_amount);
    unsigned short** nL = base;
    for (int i=row_from; i<row_to; i++) {
        int max1 = 0, max2 = 0, min1 = 0, min2 = 0, maxn, minn, np1, np2, np3, min, max;
        for (int j=2; j<W-2; j++) {
            int diff = base[i][j] - blurmap[i][j];
            if (ABS(diff) > params->sharpening.threshold) {
                // compute maximum/minimum in a delta environment
                np1 = 2*(nL[i-2][j] + nL[i-2][j+1] + nL[i-2][j+2] + nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i][j] + nL[i][j+1] + nL[i][j+2]) / 27 + nL[i-1][j+1] / 3;
                np2 = 2*(nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i][j] + nL[i][j+1] + nL[i][j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2]) / 27 + nL[i][j+1] / 3;
                np3 = 2*(nL[i][j] + nL[i][j+1] + nL[i][j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2] + nL[i+2][j] + nL[i+2][j+1] + nL[i+2][j+2]) / 27 + nL[i+1][j+1] / 3;
                MINMAX3(np1,np2,np3,maxn,minn);
                MAX3(max1,max2,maxn,max);
                MIN3(min1,min2,minn,min);
                max1 = max2; max2 = maxn;
                min1 = min2; min2 = minn;
                if (max < lab->L[i][j])
                    max = lab->L[i][j];
                if (min > lab->L[i][j])
                    min = lab->L[i][j];
                int val = lab->L[i][j] + params->sharpening.amount * diff / 100;
                int newL = CLIP(val);
                // applying halo control
                if (newL > max)
                    newL = max + (newL-max) * scale / 10000;
                else if (newL<min)
                    newL = min - (min-newL) * scale / 10000;
                lab->L[i][j] = newL;
            }
        }  
    }
}

void ImProcFunctions::lumadenoise (LabImage* lab, const ProcParams* params, double scale, int** b2) {

//    MyTime t1, t2;
//    t1.set ();
        
    if (params->lumaDenoise.enabled && lab->W>=8 && lab->H>=8) {
            
        Dim dim1 (lab->W, lab->H, 0, lab->H/2);
        Dim dim2 (lab->W, lab->H, lab->H/2, lab->H);

        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(bilateral_unsigned), lab->L, lab->L, (unsigned short**)b2, dim1, params->lumaDenoise.radius / scale, params->lumaDenoise.edgetolerance), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(bilateral_unsigned), lab->L, lab->L, (unsigned short**)b2, dim2, params->lumaDenoise.radius / scale, params->lumaDenoise.edgetolerance), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else {
            bilateral_unsigned (lab->L, lab->L, (unsigned short**)b2, dim1, params->lumaDenoise.radius / scale, params->lumaDenoise.edgetolerance);
            bilateral_unsigned (lab->L, lab->L, (unsigned short**)b2, dim2, params->lumaDenoise.radius / scale, params->lumaDenoise.edgetolerance);
        }
    }
//    t2.set ();
//    printf ("Luminance denoising time = %d\n", t2.etime (t1));
}

void ImProcFunctions::colordenoise (LabImage* lab, const ProcParams* params, double scale, int** b2) {

  if (params->colorDenoise.enabled && lab->W>=8 && lab->H>=8) {

/*    if (params->colorDenoise.edgesensitive) {

        short** buffer1 = (short**)b2;
        short** buffer2 = new short*[lab->H];
        for (int i=0; i<lab->H; i++)
            buffer2[i] = buffer1[i]+lab->W;
        Dim dim (lab->W, lab->H, 0, lab->H);
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(bilateral_signed), lab->a, lab->a, buffer1, dim, params->colorDenoise.radius / scale, params->colorDenoise.edgetolerance), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(bilateral_signed), lab->b, lab->b, buffer2, dim, params->colorDenoise.radius / scale, params->colorDenoise.edgetolerance), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else {
            bilateral_signed (lab->a, lab->a, buffer1, dim, params->colorDenoise.radius / scale, params->colorDenoise.edgetolerance);
            bilateral_signed (lab->b, lab->b, buffer1, dim, params->colorDenoise.radius / scale, params->colorDenoise.edgetolerance);
        }
        delete [] buffer2;
    }
    else {
*/
        AlignedBuffer<double>* buffer1 = new AlignedBuffer<double> (MAX(lab->W,lab->H)*5);
        AlignedBuffer<double>* buffer2 = new AlignedBuffer<double> (MAX(lab->W,lab->H)*5);
  
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_signed), lab->a, lab->a, buffer1, lab->W, 0, lab->H, params->colorDenoise.amount / 10.0 / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussHorizontal_signed), lab->b, lab->b, buffer2, lab->W, 0, lab->H, params->colorDenoise.amount / 10.0 / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
            thread1 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_signed), lab->a, lab->a, buffer1, lab->H, 0, lab->W, params->colorDenoise.amount / 10.0 / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread2 = Glib::Thread::create(sigc::bind(sigc::ptr_fun(gaussVertical_signed), lab->b, lab->b, buffer2, lab->H, 0, lab->W, params->colorDenoise.amount / 10.0 / scale), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else {
            gaussHorizontal_signed (lab->a, lab->a, buffer1, lab->W, 0, lab->H, params->colorDenoise.amount / 10.0 / scale);
            gaussHorizontal_signed (lab->b, lab->b, buffer1, lab->W, 0, lab->H, params->colorDenoise.amount / 10.0 / scale);
            gaussVertical_signed (lab->a, lab->a, buffer1, lab->H, 0, lab->W, params->colorDenoise.amount / 10.0 / scale);
            gaussVertical_signed (lab->b, lab->b, buffer1, lab->H, 0, lab->W, params->colorDenoise.amount / 10.0 / scale);
        }
        delete buffer1;
        delete buffer2;
//    }
  }
}

void ImProcFunctions::vignetting_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;

  double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
  double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;
  
  double mul = (1.0-v) / tanh(b);
  
  int val;
  for (int y=row_from; y<row_to; y++) {
      double y_d = (double) (y + cy) - h2 ;
      for (int x=0; x<transformed->width; x++) {
          double x_d = (double) (x + cx) - w2 ;
          double r = sqrt(x_d*x_d + y_d*y_d);
          double vign = v + mul * tanh (b*(maxRadius-r) / maxRadius);
          val = original->r[y][x] / vign;
          transformed->r[y][x] = CLIP(val);
          val =  original->g[y][x] / vign;
          transformed->g[y][x] = CLIP(val);
          val = original->b[y][x] / vign;
          transformed->b[y][x] = CLIP(val);
      }    
  } 
}

void ImProcFunctions::vignetting (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int oW, int oH) {

    STemp sizes;
    sizes.cx = cx;
    sizes.cy = cy;
    sizes.oW = oW;
    sizes.oH = oH;

    if (settings->dualThreadEnabled) {
        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::vignetting_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::vignetting_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else 
        vignetting_ (original, transformed, params, sizes, 0, transformed->height);
}

#include "cubint.cc"
void ImProcFunctions::transform_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;
  int sx = sizes.sx;
  int sy = sizes.sy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate.degree * 3.14/180.0);
  double sint = sin(params->rotate.degree * 3.14/180.0);

  double  max_x = (double) (sx + original->width - 1); 
  double  max_y = (double) (sy + original->height - 1);  
  double  min_x = (double) sx; 
  double  min_y = (double) sy; 
    
  const int n2 = 2;
  const int n = 4;

  int mix  = original->width - 1; // maximum x-index src
  int miy  = original->height - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (oW>oH) ? (double)oW / 2.0 : (double)oH / 2.0 ;
  double radius = sqrt( (double)( oW*oW + oH*oH ) );
  radius /= (oW<oH) ? oW : oH;

  double a = params->distortion.amount;
  
  double d = 1.0 - a;
  
    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan((double)MIN(oH,oW)/MAX(oW,oH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    // 1. check upper and lower border
    double d1 = rotmagn - a*h2/scale;
    double d2 = rotmagn - a*w2/scale;
    double d3 = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;
    d = MIN(d,MIN(d1,MIN(d2,d3)));

    // auxilary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2 / scale;

    double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;
  
    double mul = (1.0-v) / tanh(b);

    // main cycle
    double eps = 1e-10;
    for (int y=row_from; y<row_to; y++) {
        double y_d = (double) (y + cy) - h2 ;
        for (int x=0; x<transformed->width; x++) {
            double x_d = (double) (x + cx) - w2 ;

            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)	
	            s = a * r + d;

            double Dx = s*(x_d * cost - y_d * sint) + w2;
            double Dy = s*(x_d * sint + y_d * cost) + h2;

            if (fabs(Dx)<eps) Dx = 0;
            if (fabs(Dy)<eps) Dy = 0;
            if (fabs(Dx-max_x)<eps) Dx = nextafter(max_x,0);
            if (fabs(Dy-max_y)<eps) Dy = nextafter(max_y,0);

            bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

            // Convert only valid pixels
            if (valid) {
                // Extract integer and fractions of source screen coordinates
                int xc  =  (int) (Dx); Dx -= (double)xc;
                int yc  =  (int) (Dy); Dy -= (double)yc;
                int ys = yc +1 - n2 - sy; // smallest y-index used for interpolation
                int xs = xc +1 - n2 - sx; // smallest x-index used for interpolation

                double vignmul = 1.0 / (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2)   // all interpolation pixels inside image
                    cubint (original, xs, ys, Dx, Dy, &(transformed->r[y][x]), &(transformed->g[y][x]), &(transformed->b[y][x]), vignmul); 
                else { // edge pixels
                    int y1 = (yc>0) ? yc : 0;
                    if (y1>miy) y1 = miy;
                    int y2 = (yc<miy) ? yc+1 : miy;
                    if (y2<0) y2 = 0;
                    int x1 = (xc>0) ? xc : 0;
                    if (x1>mix) x1 = mix;
                    int x2 = (xc<mix) ? xc+1 : mix;
                    if (x2<0) x2 = 0;
                    int r = vignmul*(original->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->r[y1][x2]*Dx*(1.0-Dy) + original->r[y2][x1]*(1.0-Dx)*Dy + original->r[y2][x2]*Dx*Dy);
                    int g = vignmul*(original->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->g[y1][x2]*Dx*(1.0-Dy) + original->g[y2][x1]*(1.0-Dx)*Dy + original->g[y2][x2]*Dx*Dy);
                    int b = vignmul*(original->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->b[y1][x2]*Dx*(1.0-Dy) + original->b[y2][x1]*(1.0-Dx)*Dy + original->b[y2][x2]*Dx*Dy);
                    transformed->r[y][x] = CLIP(r);
                    transformed->g[y][x] = CLIP(g);
                    transformed->b[y][x] = CLIP(b);
                }
            }      
            else {
                // not valid (source pixel x,y not inside source image, etc.)
                transformed->r[y][x] = 0;
                transformed->g[y][x] = 0;
                transformed->b[y][x] = 0;
            }
        }
    } 
}

void ImProcFunctions::simpltransform_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;
  int sx = sizes.sx;
  int sy = sizes.sy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate.degree * 3.14/180.0);
  double sint = sin(params->rotate.degree * 3.14/180.0);

  double  max_x = (double) (sx + original->width - 1); 
  double  max_y = (double) (sy + original->height - 1);  
  double  min_x = (double) sx; 
  double  min_y = (double) sy; 
    
  const int n2 = 2;
  const int n = 2;

  int mix  = original->width - 1; // maximum x-index src
  int miy  = original->height - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (oW>oH) ? (double)oW / 2.0 : (double)oH / 2.0 ;
  double radius = sqrt( (double)( oW*oW + oH*oH ) );
  radius /= (oW<oH) ? oW : oH;

  double a = params->distortion.amount;
  
  double d = 1.0 - a;
  
    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan((double)MIN(oH,oW)/MAX(oW,oH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    // 1. check upper and lower border
    double d1r = rotmagn - a*h2/scale - params->cacorrection.red;
    double d2r = rotmagn - a*w2/scale - params->cacorrection.red;
    double d3r = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.red;   
    double dr = MIN(d,MIN(d1r,MIN(d2r,d3r)));
    double d1b = rotmagn - a*h2/scale - params->cacorrection.blue;
    double d2b = rotmagn - a*w2/scale - params->cacorrection.blue;
    double d3b = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.blue;   
    double db = MIN(d,MIN(d1b,MIN(d2b,d3b)));
    double d1g = rotmagn - a*h2/scale;
    double d2g = rotmagn - a*w2/scale;
    double d3g = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;   
    double dg = MIN(d,MIN(d1g,MIN(d2g,d3g)));

    d = MIN(dg,MIN(dr,db));
        
    // auxilary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2 / scale;

    double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;
  
    double mul = (1.0-v) / tanh(b);

    // main cycle
    double eps = 1e-10;
    for (int y=row_from; y<row_to; y++) {
        double y_d = (double) (y + cy) - h2 ;
        for (int x=0; x<transformed->width; x++) {
            double x_d = (double) (x + cx) - w2 ;

            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)	
	            s = a * r + d;

            double Dx = s*(x_d * cost - y_d * sint) + w2;
            double Dy = s*(x_d * sint + y_d * cost) + h2;

            if (fabs(Dx)<eps) Dx = 0;
            if (fabs(Dy)<eps) Dy = 0;
            if (fabs(Dx-max_x)<eps) Dx = nextafter(max_x,0);
            if (fabs(Dy-max_y)<eps) Dy = nextafter(max_y,0);

            bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

            // Convert only valid pixels
            if (valid) {
                // Extract integer and fractions of source screen coordinates
                int xc  =  (int) (Dx); Dx -= (double)xc;
                int yc  =  (int) (Dy); Dy -= (double)yc;
                int ys = yc +1 - n2 - sy; // smallest y-index used for interpolation
                int xs = xc +1 - n2 - sx; // smallest x-index used for interpolation

                double vignmul = 1.0 / (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2) {   // all interpolation pixels inside image

                    int r = vignmul*(original->r[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->r[yc][xc+1]*Dx*(1.0-Dy) + original->r[yc+1][xc]*(1.0-Dx)*Dy + original->r[yc+1][xc+1]*Dx*Dy);
                    int g = vignmul*(original->g[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->g[yc][xc+1]*Dx*(1.0-Dy) + original->g[yc+1][xc]*(1.0-Dx)*Dy + original->g[yc+1][xc+1]*Dx*Dy);
                    int b = vignmul*(original->b[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->b[yc][xc+1]*Dx*(1.0-Dy) + original->b[yc+1][xc]*(1.0-Dx)*Dy + original->b[yc+1][xc+1]*Dx*Dy);
                    transformed->r[y][x] = CLIP(r);
                    transformed->g[y][x] = CLIP(g);
                    transformed->b[y][x] = CLIP(b);
                }
                else { // edge pixels
                    int y1 = (yc>0) ? yc : 0;
                    if (y1>miy) y1 = miy;
                    int y2 = (yc<miy) ? yc+1 : miy;
                    if (y2<0) y2 = 0;
                    int x1 = (xc>0) ? xc : 0;
                    if (x1>mix) x1 = mix;
                    int x2 = (xc<mix) ? xc+1 : mix;
                    if (x2<0) x2 = 0;
                    int r = vignmul*(original->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->r[y1][x2]*Dx*(1.0-Dy) + original->r[y2][x1]*(1.0-Dx)*Dy + original->r[y2][x2]*Dx*Dy);
                    int g = vignmul*(original->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->g[y1][x2]*Dx*(1.0-Dy) + original->g[y2][x1]*(1.0-Dx)*Dy + original->g[y2][x2]*Dx*Dy);
                    int b = vignmul*(original->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->b[y1][x2]*Dx*(1.0-Dy) + original->b[y2][x1]*(1.0-Dx)*Dy + original->b[y2][x2]*Dx*Dy);
                    transformed->r[y][x] = CLIP(r);
                    transformed->g[y][x] = CLIP(g);
                    transformed->b[y][x] = CLIP(b);
                }
            }      
            else {
                // not valid (source pixel x,y not inside source image, etc.)
                transformed->r[y][x] = 0;
                transformed->g[y][x] = 0;
                transformed->b[y][x] = 0;
            }
        }
    } 
}


#include "cubintch.cc"
void ImProcFunctions::transform_sep_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;
  int sx = sizes.sx;
  int sy = sizes.sy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate.degree * 3.14/180.0);
  double sint = sin(params->rotate.degree * 3.14/180.0);

  double  max_x = (double) (sx + original->width - 1); 
  double  max_y = (double) (sy + original->height - 1);  
  double  min_x = (double) sx; 
  double  min_y = (double) sy; 
    
  const int n2 = 2;
  const int n = 4;

  int mix  = original->width - 1; // maximum x-index src
  int miy  = original->height - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (oW>oH) ? (double)oW / 2.0 : (double)oH / 2.0 ;
  double radius = sqrt( (double)( oW*oW + oH*oH ) );
  radius /= (oW<oH) ? oW : oH;

  double a = params->distortion.amount;
    double d = 1.0 - a;

    double cdist[3];
    cdist[0] = params->cacorrection.red;
    cdist[1] = 0.0;
    cdist[2] = params->cacorrection.blue;

    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan((double)MIN(oH,oW)/MAX(oW,oH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    // 1. check upper and lower border   
    double d1r = rotmagn - a*h2/scale - params->cacorrection.red;
    double d2r = rotmagn - a*w2/scale - params->cacorrection.red;
    double d3r = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.red;   
    double dr = MIN(d,MIN(d1r,MIN(d2r,d3r)));
    double d1b = rotmagn - a*h2/scale - params->cacorrection.blue;
    double d2b = rotmagn - a*w2/scale - params->cacorrection.blue;
    double d3b = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.blue;   
    double db = MIN(d,MIN(d1b,MIN(d2b,d3b)));
    double d1g = rotmagn - a*h2/scale;
    double d2g = rotmagn - a*w2/scale;
    double d3g = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;   
    double dg = MIN(d,MIN(d1g,MIN(d2g,d3g)));

    d = MIN(dg,MIN(dr,db));
    
    unsigned short** chorig[3];
    chorig[0] = original->r;
    chorig[1] = original->g;
    chorig[2] = original->b;
    
    unsigned short** chtrans[3];
    chtrans[0] = transformed->r;
    chtrans[1] = transformed->g;
    chtrans[2] = transformed->b;


    // auxilary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2 / scale;

    double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;
  
    double mul = (1.0-v) / tanh(b);

    // main cycle
    double eps = 1e-10;
    for (int y=row_from; y<row_to; y++) {
        double y_d = (double) (y + cy) - h2 ;
        for (int x=0; x<transformed->width; x++) {
            double x_d = (double) (x + cx) - w2 ;

            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)	
	            s = a * r + d;

            double vignmul = 1.0 / (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

            for (int c=0; c<3; c++) {

                double Dx = (s + cdist[c]) * (x_d * cost - y_d * sint) + w2;
                double Dy = (s + cdist[c]) * (x_d * sint + y_d * cost) + h2;

                if (fabs(Dx)<eps) Dx = 0;
                if (fabs(Dy)<eps) Dy = 0;
                if (fabs(Dx-max_x)<eps) Dx = nextafter(max_x,0);
                if (fabs(Dy-max_y)<eps) Dy = nextafter(max_y,0);

                bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

                // Convert only valid pixels
                if (valid) {
                    // Extract integer and fractions of source screen coordinates
                    int xc  =  (int) (Dx); Dx -= (double)xc;
                    int yc  =  (int) (Dy); Dy -= (double)yc;
                    int ys = yc +1 - n2 - sy; // smallest y-index used for interpolation
                    int xs = xc +1 - n2 - sx; // smallest x-index used for interpolation

                    if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2)  // all interpolation pixels inside image
                        cubintch (chorig[c], xs, ys, Dx, Dy, &(chtrans[c][y][x]), vignmul); 
                    else {// edge pixels, linear interpolation
                        int y1 = (yc>0) ? yc : 0;
                        if (y1>miy) y1 = miy;
                        int y2 = (yc<miy) ? yc+1 : miy;
                        if (y2<0) y2 = 0;
                        int x1 = (xc>0) ? xc : 0;
                        if (x1>mix) x1 = mix;
                        int x2 = (xc<mix) ? xc+1 : mix;
                        if (x2<0) x2 = 0;
                        int val = vignmul*(chorig[c][y1][x1]*(1.0-Dx)*(1.0-Dy) + chorig[c][y1][x2]*Dx*(1.0-Dy) + chorig[c][y2][x1]*(1.0-Dx)*Dy + chorig[c][y2][x2]*Dx*Dy);
                        chtrans[c][y][x] = CLIP(val);
                    }
                }      
                else // not valid (source pixel x,y not inside source image, etc.)
                    chtrans[c][y][x] = 0;
            }
        }
    } 
}

bool ImProcFunctions::transCoord (const ProcParams* params, int W, int H, std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue) {

    bool clipresize = true;
    bool clipped = false;
    
    red.clear ();
    green.clear ();
    blue.clear ();
    bool needstransform  = fabs(params->rotate.degree)>1e-15 || fabs(params->distortion.amount)>1e-15 || fabs(params->cacorrection.red)>1e-15 || fabs(params->cacorrection.blue)>1e-15;
    if (!needstransform) {
        if (clipresize) {
            // Apply resizing
            if (fabs(params->resize.scale-1.0)>=1e-7) {
                for (int i=0; i<src.size(); i++) {
                    red.push_back   (Coord2D (src[i].x / params->resize.scale, src[i].y / params->resize.scale));
                    green.push_back (Coord2D (src[i].x / params->resize.scale, src[i].y / params->resize.scale));
                    blue.push_back  (Coord2D (src[i].x / params->resize.scale, src[i].y / params->resize.scale));
                }
                for (int i=0; i<src.size(); i++) {
                    red[i].x = CLIPTOC(red[i].x,0,W-1,clipped);
                    red[i].y = CLIPTOC(red[i].y,0,H-1,clipped);
                    green[i].x = CLIPTOC(green[i].x,0,W-1,clipped);
                    green[i].y = CLIPTOC(green[i].y,0,H-1,clipped);
                    blue[i].x = CLIPTOC(blue[i].x,0,W-1,clipped);
                    blue[i].y = CLIPTOC(blue[i].y,0,H-1,clipped);
                }
            }
            else 
                for (int i=0; i<src.size(); i++) {
                    red.push_back   (Coord2D (src[i].x, src[i].y));
                    green.push_back (Coord2D (src[i].x, src[i].y));
                    blue.push_back  (Coord2D (src[i].x, src[i].y));
                }
        }
        return clipped;
    }
    double rW = W*params->resize.scale;
    double rH = H*params->resize.scale;
    double  w2 = (double) rW  / 2.0 - 0.5;
    double  h2 = (double) rH  / 2.0 - 0.5;
    double cost = cos(params->rotate.degree * 3.14/180.0);
    double sint = sin(params->rotate.degree * 3.14/180.0);

    double scale = (rW>rH) ? rW / 2.0 : rH / 2.0 ;
    double radius = sqrt ((double)(rW*rW + rH*rH ));
    radius /= (rW<rH) ? rW : rH;
    double a = params->distortion.amount; 
    double d = 1.0 - a;
  
    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan(MIN(rH,rW)/MAX(rW,rH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    if (params->cacorrection.red==0 && params->cacorrection.blue==0) {
        // 1. check upper and lower border
        double d1 = rotmagn - a*h2/scale;
        double d2 = rotmagn - a*w2/scale;
        double d3 = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;
        d = MIN(d,MIN(d1,MIN(d2,d3)));

        for (int i=0; i<src.size(); i++) {
            double y_d = src[i].y - h2 ;
            double x_d = src[i].x - w2 ;
            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)	
                s = a * r + d;
            red.push_back (Coord2D(s*(x_d * cost - y_d * sint) + w2, s*(x_d * sint + y_d * cost) + h2));
            green.push_back (Coord2D(s*(x_d * cost - y_d * sint) + w2, s*(x_d * sint + y_d * cost) + h2));
            blue.push_back (Coord2D(s*(x_d * cost - y_d * sint) + w2, s*(x_d * sint + y_d * cost) + h2));
        }   
    }
    else {
        double cdist[3];
        cdist[0] = params->cacorrection.red;
        cdist[1] = 0.0;
        cdist[2] = params->cacorrection.blue;

        // 1. check upper and lower border   
        double d1r = rotmagn - a*h2/scale - params->cacorrection.red;
        double d2r = rotmagn - a*w2/scale - params->cacorrection.red;
        double d3r = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.red;   
        double dr = MIN(d,MIN(d1r,MIN(d2r,d3r)));
        double d1b = rotmagn - a*h2/scale - params->cacorrection.blue;
        double d2b = rotmagn - a*w2/scale - params->cacorrection.blue;
        double d3b = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.blue;   
        double db = MIN(d,MIN(d1b,MIN(d2b,d3b)));
        double d1g = rotmagn - a*h2/scale;
        double d2g = rotmagn - a*w2/scale;
        double d3g = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;   
        double dg = MIN(d,MIN(d1g,MIN(d2g,d3g)));

        d = MIN(dg,MIN(dr,db));

        for (int i=0; i<src.size(); i++) {
            double y_d = src[i].y - h2 ;
            double x_d = src[i].x - w2 ;
            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)	
                s = a * r + d;
            src[i].x = s*(x_d * cost - y_d * sint) + w2;
            src[i].y  = s*(x_d * sint + y_d * cost) + h2;

            red.push_back (Coord2D((s+cdist[0])*(x_d * cost - y_d * sint) + w2, (s+cdist[0])*(x_d * sint + y_d * cost) + h2));
            green.push_back (Coord2D((s+cdist[1])*(x_d * cost - y_d * sint) + w2, (s+cdist[1])*(x_d * sint + y_d * cost) + h2));
            blue.push_back (Coord2D((s+cdist[2])*(x_d * cost - y_d * sint) + w2, (s+cdist[2])*(x_d * sint + y_d * cost) + h2));
        }
    }
    
    if (clipresize) {
        if (fabs(params->resize.scale-1.0)>=1e-7) {
            for (int i=0; i<src.size(); i++) {
                red[i].x /= params->resize.scale;
                red[i].y /= params->resize.scale;
                green[i].x /= params->resize.scale;
                green[i].y /= params->resize.scale;
                blue[i].x /= params->resize.scale;
                blue[i].y /= params->resize.scale;
            }
        }
        for (int i=0; i<src.size(); i++) {
            red[i].x = CLIPTOC(red[i].x,0,W-1,clipped);
            red[i].y = CLIPTOC(red[i].y,0,H-1,clipped);
            green[i].x = CLIPTOC(green[i].x,0,W-1,clipped);
            green[i].y = CLIPTOC(green[i].y,0,H-1,clipped);
            blue[i].x = CLIPTOC(blue[i].x,0,W-1,clipped);
            blue[i].y = CLIPTOC(blue[i].y,0,H-1,clipped);
        }
    }
    return clipped;
}

bool ImProcFunctions::transCoord (const ProcParams* params, int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv) {

    int x1 = x, y1 = y;
    int x2 = x1 + w - 1;
    int y2 = y1 + h - 1;

    std::vector<Coord2D> corners (8);
    corners[0].set (x1, y1);
    corners[1].set (x1, y2);
    corners[2].set (x2, y2);
    corners[3].set (x2, y1);
    corners[4].set ((x1+x2)/2, y1);
    corners[5].set ((x1+x2)/2, y2);
    corners[6].set (x1, (y1+y2)/2);
    corners[7].set (x2, (y1+y2)/2);

    std::vector<Coord2D> r, g, b;

    bool result = transCoord (params, W, H, corners, r, g, b);

    std::vector<Coord2D> transCorners;
    transCorners.insert (transCorners.end(), r.begin(), r.end());
    transCorners.insert (transCorners.end(), g.begin(), g.end());
    transCorners.insert (transCorners.end(), b.begin(), b.end());
        
    double x1d = transCorners[0].x;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].x<x1d)
            x1d = transCorners[i].x;
   int x1v = (int)(x1d);
            
    double y1d = transCorners[0].y;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].y<y1d)
            y1d = transCorners[i].y;
    int y1v = (int)(y1d);
            
    double x2d = transCorners[0].x;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].x>x2d)
            x2d = transCorners[i].x;
    int x2v = (int)ceil(x2d);

    double y2d = transCorners[0].y;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].y>y2d)
            y2d = transCorners[i].y;
    int y2v = (int)ceil(y2d);

    xv = x1v;
    yv = y1v;
    wv = x2v - x1v + 1;
    hv = y2v - y1v + 1;

    return result;
}

void ImProcFunctions::transform (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int sx, int sy, int oW, int oH) {

    STemp sizes;
    sizes.cx = cx;
    sizes.cy = cy;
    sizes.oW = oW;
    sizes.oH = oH;
    sizes.sx = sx;
    sizes.sy = sy;
    
    if (params->cacorrection.red==0 && params->cacorrection.blue==0) {
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            transform_ (original, transformed, params, sizes, 0, transformed->height);
    }
    else {
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_sep_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_sep_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else 
            transform_sep_ (original, transformed, params, sizes, 0, transformed->height);
    }
}

void ImProcFunctions::simpltransform (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int sx, int sy, int oW, int oH) {

    STemp sizes;
    sizes.cx = cx;
    sizes.cy = cy;
    sizes.oW = oW;
    sizes.oH = oH;
    sizes.sx = sx;
    sizes.sy = sy;
    
    if (settings->dualThreadEnabled) {
        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::simpltransform_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::simpltransform_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else 
        simpltransform_ (original, transformed, params, sizes, 0, transformed->height);
}
/*void ImProcFunctions::transform (Image16* original, Image16* transformed, const ProcParams* params, int ox, int oy) {

  if (!transformed)
    return;
    
  int oW = W, oH = H, tW = W, tH = H;

  double  w2 = (double) tW / 2.0 - 0.5;
  double  h2 = (double) tH / 2.0 - 0.5;
  double  sw2 = (double) oW  / 2.0 - 0.5;
  double  sh2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate_fine * 3.14/180.0);
  double sint = sin(params->rotate_fine * 3.14/180.0);

  double  max_x = (double) oW; 
  double  max_y = (double) oH; 
  double  min_x =  0.0; 
  double  min_y =  0.0; 

  const int n2 = 2;
  const int n = 4;

  int mix  = oW - 1; // maximum x-index src
  int miy  = oH - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (tW>tH) ? (double)tW / 2.0 : (double)tH / 2.0 ;
  double radius = sqrt( (double)( tW*tW + tH*tH ) );
  radius /= (tW<tH) ? tW : tH;

  double a = params->lens_distortion;

  for (int y=0; y<transformed->height; y++) {
    double y_d = (double) y + oy - h2 ;

    for (int x=0; x<transformed->width; x++) {
      double x_d = (double) x + ox - w2 ;

      double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
      double s = 10000.0;
      if (r<radius)	
	        s = a * r + 1.0 - a;

      double Dx = s*(x_d * cost - y_d * sint) + sw2;
      double Dy = s*(x_d * sint + y_d * cost) + sh2;

      bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

      // Convert only valid pixels
      if (valid) {
        // Extract integer and fractions of source screen coordinates
        int xc  =  (int) floor (Dx) ; Dx -= (double)xc;
        int yc  =  (int) floor (Dy) ; Dy -= (double)yc;
        int ys = yc +1 - n2 ; // smallest y-index used for interpolation
        int xs = xc +1 - n2 ; // smallest x-index used for interpolation

        unsigned short sr[2][2], sg[2][2], sb[2][2];

        if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2)   // all interpolation pixels inside image
          cubint (original, xs, ys, Dx, Dy, &(transformed->r[y][x]), &(transformed->g[y][x]), &(transformed->b[y][x])); 
        else { // edge pixels
          transformed->r[y][x] = 0;
          transformed->g[y][x] = 0;
          transformed->b[y][x] = 0;
        }
      }      
      else {
        // not valid (source pixel x,y not inside source image, etc.)
        transformed->r[y][x] = 0;
        transformed->g[y][x] = 0;
        transformed->b[y][x] = 0;
      }
    }
  } 
}*/

void ImProcFunctions::lab2rgb (LabImage* lab, Image8* image) {

    if (settings->dualThreadEnabled) {

        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::lab2rgb_), lab, image, 0, lab->H/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::lab2rgb_), lab, image, lab->H/2, lab->H), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else
        lab2rgb_ (lab, image, 0, lab->H);
}

void ImProcFunctions::lab2rgb_ (LabImage* lab, Image8* image, int row_from, int row_to) {

    int X, Y, Z;
    unsigned short** nL = lab->L;
    short** na = lab->a;
    short** nb = lab->b;
    int tW = lab->W;
    int ix = row_from*tW*3;

	if (monitorTransform) {
        short* buffer = new short [3*tW];
		for (int i=row_from; i<row_to; i++) {
			unsigned short* rL = nL[i];
			short* ra = na[i];
			short* rb = nb[i];
			int iy = 0;
			for (register int j=0; j<tW; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);
                    
                Y = ycache[y_];
				X = xcache[x_];
				Z = zcache[z_];

                buffer[iy++] = CLIP(X);
                buffer[iy++] = CLIP(Y);
                buffer[iy++] = CLIP(Z);
			}
            cmsDoTransform (monitorTransform, buffer, image->data + ix, tW);
            ix += 3*tW;
		}
        delete [] buffer;
	}
	else {
		for (int i=row_from; i<row_to; i++) {
			unsigned short* rL = nL[i];
			short* ra = na[i];
			short* rb = nb[i];
			for (register int j=0; j<tW; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);
                    
                Y = ycache[y_];
				X = xcache[x_];
				Z = zcache[z_];

				/* XYZ-D50 to RGB */
				int R = (25689*X-13261*Y-4022*Z) >> 13;
				int G = (-8017*X+15697*Y+274*Z) >> 13;
				int B = (590*X-1877*Y+11517*Z) >> 13;

				/* copy RGB */
				image->data[ix++] = gamma2curve[CLIP(R)] >> 8;
				image->data[ix++] = gamma2curve[CLIP(G)] >> 8;
				image->data[ix++] = gamma2curve[CLIP(B)] >> 8;
			}
		}
	}
}

Image8* ImProcFunctions::lab2rgb (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {

    int tW = lab->W;
    int tH = lab->H;

    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>tW) cw = tW-cx;
    if (cy+ch>tH) ch = tH-cy;
    
    Image8* image = new Image8 (cw, ch);

    int X, Y, Z;
    int ix = 0;
    unsigned short** nL = lab->L;
    short** na = lab->a;
    short** nb = lab->b;
    
    cmsHPROFILE oprof = iccStore.getProfile (profile);
    
    if (oprof) {
        cmsHPROFILE iprof = iccStore.getXYZProfile ();       
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprof, TYPE_RGB_8, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
        short* buffer = new short [3*cw];
        for (int i=cy; i<cy+ch; i++) {
            unsigned short* rL = nL[i];
            short* ra = na[i];
            short* rb = nb[i];
            int iy = 0;
            for (register int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);
                    
                Y = ycache[y_];
				X = xcache[x_];
				Z = zcache[z_];

                buffer[iy++] = CLIP(X);
                buffer[iy++] = CLIP(Y);
                buffer[iy++] = CLIP(Z);
            }
            cmsDoTransform (hTransform, buffer, image->data + ix, cw);
            ix += 3*cw;
        }
        delete [] buffer;
        cmsDeleteTransform(hTransform);
    }
    else {   
        for (int i=cy; i<cy+ch; i++) {
            unsigned short* rL = nL[i];
            short* ra = na[i];
            short* rb = nb[i];
            for (register int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);
                    
                Y = ycache[y_];
				X = xcache[x_];
				Z = zcache[z_];

                int R = (25689*X-13261*Y-4022*Z) >> 13;
                int G = (-8017*X+15697*Y+274*Z) >> 13;
                int B = (590*X-1877*Y+11517*Z) >> 13;

                image->data[ix++] = gamma2curve[CLIP(R)] >> 8;
                image->data[ix++] = gamma2curve[CLIP(G)] >> 8;
                image->data[ix++] = gamma2curve[CLIP(B)] >> 8;
            }
        }
    }
    return image;
}

Image16* ImProcFunctions::lab2rgb16 (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {

    int tW = lab->W;
    int tH = lab->H;

    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>tW) cw = tW-cx;
    if (cy+ch>tH) ch = tH-cy;
    
    Image16* image = new Image16 (cw, ch);

    int X, Y, Z;
    int ix = 0;
    unsigned short** nL = lab->L;
    short** na = lab->a;
    short** nb = lab->b;
    
    cmsHPROFILE oprof = iccStore.getProfile (profile);
    
    if (oprof) {
		for (int i=cy; i<cy+ch; i++) {
			unsigned short* rL = nL[i];
			short* ra = na[i];
			short* rb = nb[i];
			short* xa = (short*)image->r[i-cy];
			short* ya = (short*)image->g[i-cy];
			short* za = (short*)image->b[i-cy];
			for (register int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);
                    
                Y = ycache[y_];
				X = xcache[x_];
				Z = zcache[z_];

				xa[j-cx] = CLIP(X);
				ya[j-cx] = CLIP(Y);
				za[j-cx] = CLIP(Z);
			}
		}
        cmsHPROFILE iprof = iccStore.getXYZProfile ();       
        lcmsMutex->lock ();
		cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprof, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
		cmsDoTransform (hTransform, image->data, image->data, image->planestride/2);
		cmsDeleteTransform(hTransform);	
	}
	else {
		for (int i=cy; i<cy+ch; i++) {
			unsigned short* rL = nL[i];
			short* ra = na[i];
			short* rb = nb[i];
			for (register int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);
                    
                Y = ycache[y_];
				X = xcache[x_];
				Z = zcache[z_];

				int R = (25689*X-13261*Y-4022*Z) >> 13;
				int G = (-8017*X+15697*Y+274*Z) >> 13;
				int B = (590*X-1877*Y+11517*Z) >> 13;

				image->r[i-cy][j-cx] = gamma2curve[CLIP(R)];
				image->g[i-cy][j-cx] = gamma2curve[CLIP(G)];
				image->b[i-cy][j-cx] = gamma2curve[CLIP(B)];
			}
		}
	}
    return image;
}

void ImProcFunctions::resize (Image16* src, Image16* dst, ResizeParams params) {

    if (settings->dualThreadEnabled) {

        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::resize_), src, dst, params, 0, dst->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::resize_), src, dst, params, dst->height/2, dst->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else
        resize_ (src, dst, params, 0, dst->height);
}

void ImProcFunctions::resize_ (Image16* src, Image16* dst, ResizeParams params, int row_from, int row_to) {

    if (params.method.substr(0,7)=="Bicubic") {
        double Av = -0.5;
        if (params.method=="Bicubic (Sharper)")
            Av = -0.75;
        else if (params.method=="Bicubic (Softer)")
            Av = -0.25;
        double wx[4], wy[4];
        for (int i=row_from; i<row_to; i++) {
            double Dy = i / params.scale;
            int yc  =  (int) Dy; Dy -= (double)yc;
            int ys = yc - 1; // smallest y-index used for interpolation
            // compute vertical weights
            double t1y = -Av*(Dy-1.0)*Dy;
            double t2y = (3.0-2.0*Dy)*Dy*Dy;
            wy[3] = t1y*Dy;
            wy[2] = t1y*(Dy-1.0) + t2y;
            wy[1] = -t1y*Dy + 1.0 - t2y;
            wy[0] = -t1y*(Dy-1.0);
            for (int j=0; j<dst->width; j++) {
                double Dx = j / params.scale;
                int xc  =  (int) Dx; Dx -= (double)xc;
                int xs = xc - 1; // smallest x-index used for interpolation
                if (ys >= 0 && ys <src->height-3 && xs >= 0 && xs <= src->width-3) {
                    // compute horizontal weights
                    double t1 = -Av*(Dx-1.0)*Dx;
                    double t2 = (3.0-2.0*Dx)*Dx*Dx;
                    wx[3] = t1*Dx;
                    wx[2] = t1*(Dx-1.0) + t2;
                    wx[1] = -t1*Dx + 1.0 - t2;
                    wx[0] = -t1*(Dx-1.0);
                    // compute weighted sum
                    int r = 0;
                    int g = 0;
                    int b = 0;
/*                    r = wx[0]*wy[0]*src->r[ys+0][xs+0] + wx[0]*wy[1]*src->r[ys+1][xs+0] + wx[0]*wy[2]*src->r[ys+2][xs+0] + wx[0]*wy[3]*src->r[ys+3][xs+0] +
                        wx[1]*wy[0]*src->r[ys+0][xs+1] + wx[1]*wy[1]*src->r[ys+1][xs+1] + wx[1]*wy[2]*src->r[ys+2][xs+1] + wx[1]*wy[3]*src->r[ys+3][xs+1] +
                        wx[2]*wy[0]*src->r[ys+0][xs+2] + wx[2]*wy[1]*src->r[ys+1][xs+1] + wx[2]*wy[2]*src->r[ys+2][xs+2] + wx[2]*wy[3]*src->r[ys+3][xs+2] +
                        wx[3]*wy[0]*src->r[ys+0][xs+3] + wx[3]*wy[1]*src->r[ys+1][xs+1] + wx[3]*wy[2]*src->r[ys+2][xs+3] + wx[3]*wy[3]*src->r[ys+3][xs+3];
                    g = wx[0]*wy[0]*src->g[ys+0][xs+0] + wx[0]*wy[1]*src->g[ys+1][xs+0] + wx[0]*wy[2]*src->g[ys+2][xs+0] + wx[0]*wy[3]*src->g[ys+3][xs+0] +
                        wx[1]*wy[0]*src->g[ys+0][xs+1] + wx[1]*wy[1]*src->g[ys+1][xs+1] + wx[1]*wy[2]*src->g[ys+2][xs+1] + wx[1]*wy[3]*src->g[ys+3][xs+1] +
                        wx[2]*wy[0]*src->g[ys+0][xs+2] + wx[2]*wy[1]*src->g[ys+1][xs+1] + wx[2]*wy[2]*src->g[ys+2][xs+2] + wx[2]*wy[3]*src->g[ys+3][xs+2] +
                        wx[3]*wy[0]*src->g[ys+0][xs+3] + wx[3]*wy[1]*src->g[ys+1][xs+1] + wx[3]*wy[2]*src->g[ys+2][xs+3] + wx[3]*wy[3]*src->g[ys+3][xs+3];
                    b = wx[0]*wy[0]*src->b[ys+0][xs+0] + wx[0]*wy[1]*src->b[ys+1][xs+0] + wx[0]*wy[2]*src->b[ys+2][xs+0] + wx[0]*wy[3]*src->b[ys+3][xs+0] +
                        wx[1]*wy[0]*src->b[ys+0][xs+1] + wx[1]*wy[1]*src->b[ys+1][xs+1] + wx[1]*wy[2]*src->b[ys+2][xs+1] + wx[1]*wy[3]*src->b[ys+3][xs+1] +
                        wx[2]*wy[0]*src->b[ys+0][xs+2] + wx[2]*wy[1]*src->b[ys+1][xs+1] + wx[2]*wy[2]*src->b[ys+2][xs+2] + wx[2]*wy[3]*src->b[ys+3][xs+2] +
                        wx[3]*wy[0]*src->b[ys+0][xs+3] + wx[3]*wy[1]*src->b[ys+1][xs+1] + wx[3]*wy[2]*src->b[ys+2][xs+3] + wx[3]*wy[3]*src->b[ys+3][xs+3];*/
                    for (int x=0; x<4; x++)
                        for (int y=0; y<4; y++) {
                            double w = wx[x]*wy[y];
                            r += w*src->r[ys+y][xs+x];
                            g += w*src->g[ys+y][xs+x];
                            b += w*src->b[ys+y][xs+x];
                        }
                    dst->r[i][j] = CLIP(r);
                    dst->g[i][j] = CLIP(g);
                    dst->b[i][j] = CLIP(b);
                }
                else {
                    xc = CLIPTO(xc, 0, src->width-1);
                    yc = CLIPTO(yc, 0, src->height-1);
                    int nx = xc + 1;
                    if (nx>=src->width)
                        nx = xc;
                    int ny = yc + 1;
                    if (ny>=src->height)
                        ny = yc;
                    dst->r[i][j] = (1-Dx)*(1-Dy)*src->r[yc][xc] + (1-Dx)*Dy*src->r[ny][xc] + Dx*(1-Dy)*src->r[yc][nx] + Dx*Dy*src->r[ny][nx];
                    dst->g[i][j] = (1-Dx)*(1-Dy)*src->g[yc][xc] + (1-Dx)*Dy*src->g[ny][xc] + Dx*(1-Dy)*src->g[yc][nx] + Dx*Dy*src->g[ny][nx];
                    dst->b[i][j] = (1-Dx)*(1-Dy)*src->b[yc][xc] + (1-Dx)*Dy*src->b[ny][xc] + Dx*(1-Dy)*src->b[yc][nx] + Dx*Dy*src->b[ny][nx];
                }
            }
        }
    }
    else if (params.method=="Bilinear") {
        for (int i=row_from; i<row_to; i++) {
            int sy = i/params.scale;
            sy = CLIPTO(sy, 0, src->height-1);
            double dy = i/params.scale - sy;
            int ny = sy+1;
            if (ny>=src->height)
                ny = sy;
            for (int j=0; j<dst->width; j++) {
                int sx = j/params.scale;
                sx = CLIPTO(sx, 0, src->width-1);
                double dx = j/params.scale - sx;
                int nx = sx+1;
                if (nx>=src->width)
                    nx = sx;
                dst->r[i][j] = (1-dx)*(1-dy)*src->r[sy][sx] + (1-dx)*dy*src->r[ny][sx] + dx*(1-dy)*src->r[sy][nx] + dx*dy*src->r[ny][nx];
                dst->g[i][j] = (1-dx)*(1-dy)*src->g[sy][sx] + (1-dx)*dy*src->g[ny][sx] + dx*(1-dy)*src->g[sy][nx] + dx*dy*src->g[ny][nx];
                dst->b[i][j] = (1-dx)*(1-dy)*src->b[sy][sx] + (1-dx)*dy*src->b[ny][sx] + dx*(1-dy)*src->b[sy][nx] + dx*dy*src->b[ny][nx];
            }
        }
    }
    else {
        for (int i=row_from; i<row_to; i++) {
            int sy = i/params.scale;
            sy = CLIPTO(sy, 0, src->height-1);
            for (int j=0; j<dst->width; j++) {
                int sx = j/params.scale;
                sx = CLIPTO(sx, 0, src->width-1);
                dst->r[i][j] = src->r[sy][sx];
                dst->g[i][j] = src->g[sy][sx];
                dst->b[i][j] = src->b[sy][sx];
            }
        }
    }
}

void ImProcFunctions::getAutoExp  (int* histogram, int histcompr, double expcomp, double clip, double& br, int& bl) {

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
        double mavg = 65536.0 / (awg-bl) * (gavg - bl);
        if (awg < maxaw)
            awg = maxaw;
    }

    br = log(65535.0 / (awg-bl)) / log(2.0);   
    if (br<0)
        br = 0;

//printf ("br=%g, bl=%d, %g\n", br, bl, expcomp);
        
/*
    if (shc<avg) {
        int maxaw = (avg-shc) * 4 / 3 + shc; // dont let aw be such large that the histogram average goes above 3/4
        double mavg = 65536.0 / (aw-shc) * (avg - shc);
        if (aw < maxaw)
            aw = maxaw;
    }
    br = log(65535.0 / (aw-shc)) / log(2.0);   
*/
}
}

