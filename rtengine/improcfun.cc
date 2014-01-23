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
#include "cplx_wavelet_dec.h"
#include "boxblur.h"
#include "rt_math.h"
#include "EdgePreservingDecomposition.h"
#include "improccoordinator.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#undef CLIPD
#define CLIPD(a) ((a)>0.0f?((a)<1.0f?(a):1.0f):0.0f)
#define PIX_SORT(a,b) { if ((a)>(b)) {temp=(a);(a)=(b);(b)=temp;} }

#define med3(a0,a1,a2,a3,a4,a5,a6,a7,a8,median) { \
pp[0]=a0; pp[1]=a1; pp[2]=a2; pp[3]=a3; pp[4]=a4; pp[5]=a5; pp[6]=a6; pp[7]=a7; pp[8]=a8; \
PIX_SORT(pp[1],pp[2]); PIX_SORT(pp[4],pp[5]); PIX_SORT(pp[7],pp[8]); \
PIX_SORT(pp[0],pp[1]); PIX_SORT(pp[3],pp[4]); PIX_SORT(pp[6],pp[7]); \
PIX_SORT(pp[1],pp[2]); PIX_SORT(pp[4],pp[5]); PIX_SORT(pp[7],pp[8]); \
PIX_SORT(pp[0],pp[3]); PIX_SORT(pp[5],pp[8]); PIX_SORT(pp[4],pp[7]); \
PIX_SORT(pp[3],pp[6]); PIX_SORT(pp[1],pp[4]); PIX_SORT(pp[2],pp[5]); \
PIX_SORT(pp[4],pp[7]); PIX_SORT(pp[4],pp[2]); PIX_SORT(pp[6],pp[4]); \
PIX_SORT(pp[4],pp[2]); median=pp[4];} //pp4 = median


	

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
LUTf ImProcFunctions::cachef;
LUTf ImProcFunctions::gamma2curve;
void ImProcFunctions::initCache () {

    const int maxindex = 65536;
	cachef(maxindex,0/*LUT_CLIP_BELOW*/);

	gamma2curve(maxindex,0);

    for (int i=0; i<maxindex; i++) {
        if (i>Color::eps_max) {
            cachef[i] = 327.68*( exp(1.0/3.0 * log((double)i / MAXVALD) ));
        }
        else {
            cachef[i] = 327.68*((Color::kappa*i/MAXVALD+16.0)/116.0);
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
      
            int r = original->r(i,j);
            int g = original->g(i,j);
            int b = original->b(i,j);

            int y = CLIP((int)(lumimul[0] * r + lumimul[1] * g + lumimul[2] * b)) ;

            if (histogram) { 
                histogram[y]++;
            }
        }
    }
}
/*
void ImProcFunctions::CAT02 (Imagefloat* baseImg, const ProcParams* params)
{
	const double toxyz[3][3] = 		{{0.7976749,  0.1351917,  0.0313534},
									{0.2880402,  0.7118741,  0.0000857},
									{0.0000000,  0.0000000,  0.8252100}};

	const double xyzto[3][3] = 		{{1.3459433, -0.2556075, -0.0511118},
									{-0.5445989,  1.5081673,  0.0205351},
									{0.0000000,  0.0000000,  1.2118128}};
  int fw = baseImg->width;
  int fh = baseImg->height;

  double CAM02BB00,CAM02BB01,CAM02BB02,CAM02BB10,CAM02BB11,CAM02BB12,CAM02BB20,CAM02BB21,CAM02BB22;
  double Xxx,Yyy,Zzz;
 // Xxx=1.09844;
 // Yyy=1.0;
 // Zzz=0.355961;
 //params.wb.temperature, params.wb.green, params.wb.method
 double Xxyz, Zxyz;
// ColorTemp::temp2mulxyz (params->wb.temperature, params->wb.green, params->wb.method, Xxyz, Zxyz);
  ColorTemp::temp2mulxyz (5000.0, 1.0, "Camera", Xxyz, Zxyz);

 ColorTemp::cieCAT02(Xxx, Yyy, Zzz, CAM02BB00,CAM02BB01,CAM02BB02,CAM02BB10,CAM02BB11,CAM02BB12,CAM02BB20,CAM02BB21,CAM02BB22);
  printf("00=%f 01=%f 11=%f 20=%f 22=%f\n", CAM02BB00,CAM02BB01,CAM02BB11,CAM02BB20,CAM02BB22);

 
 	for (int i=0; i<fh; i++) {
		for (int j=0; j<fw; j++) {
			float r = baseImg->r(i,j);
			float g = baseImg->g(i,j);
			float b = baseImg->b(i,j);
			
			float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
			float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
			float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;
			float Xcam=CAM02BB00* x +CAM02BB01* y + CAM02BB02* z ;
			float Ycam=CAM02BB10* x +CAM02BB11* y + CAM02BB12* z ;
			float Zcam=CAM02BB20* x +CAM02BB21* y + CAM02BB22* z ;	
			baseImg->r(i,j) = xyzto[0][0] * Xcam + xyzto[0][1] * Ycam + xyzto[0][2] * Zcam;
			baseImg->g(i,j) = xyzto[1][0] * Xcam + xyzto[1][1] * Ycam + xyzto[1][2] * Zcam;
			baseImg->b(i,j) = xyzto[2][0] * Xcam + xyzto[2][1] * Ycam + xyzto[2][2] * Zcam;
		}
	}
}
*/
void ImProcFunctions::firstAnalysis (Imagefloat* original, const ProcParams* params, LUTu & histogram, double gamma) {

	// set up monitor transform
	Glib::ustring wprofile = params->icm.working;
	if (monitorTransform)
		cmsDeleteTransform (monitorTransform);
	monitorTransform = NULL;

#if !defined(__APPLE__) // No support for monitor profiles on OS X, all data is sRGB
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
#endif
	// calculate histogram of the y channel needed for contrast curve calculation in exposure adjustments

	int T = 1;
#ifdef _OPENMP
	if(multiThread)
		T = omp_get_max_threads();
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
	for (int j=0; j<T; j++)
		for (int i=0; i<65536; i++)
			histogram[i] += hist[j][i];

    for (int i=0; i<T; i++)
    	delete [] hist[i];
    delete [] hist;
}

// Copyright (c) 2012 Jacques Desmis <jdesmis@gmail.com>
void ImProcFunctions::ciecam_02 (CieImage* ncie, double adap, int begh, int endh, int pW, int pwb, LabImage* lab, const ProcParams* params ,
								const ColorAppearance & customColCurve1, const ColorAppearance & customColCurve2,const ColorAppearance & customColCurve3,
								LUTu & histLCAM, LUTu & histCCAM, LUTf & CAMBrightCurveJ, LUTf & CAMBrightCurveQ, float &mean, int Iterates, int scale, float** buffer, bool execsharp, double &d)
{
if(params->colorappearance.enabled) {

#ifdef _DEBUG
	MyTime t1e,t2e;
	t1e.set();
#endif
	LUTf dLcurve;
	LUTu hist16JCAM;
	bool jp=false;
	float val;
	//preparate for histograms CIECAM
	if(pW!=1){//only with improccoordinator
		dLcurve(65536,0);
		dLcurve.clear();
		hist16JCAM(65536,0);
		hist16JCAM.clear();
		for (int i=0; i<32768; i++) {  //# 32768*1.414  approximation maxi for chroma
			val = (double)i / 32767.0;
			dLcurve[i] = CLIPD(val);
		}
	}
	LUTf dCcurve;
	LUTu hist16_CCAM;
	bool chropC=false;
	float valc;

	if(pW!=1){//only with improccoordinator
		dCcurve(65536,0);
		hist16_CCAM(65536);
		hist16_CCAM.clear();
		for (int i=0; i<48000; i++) {  //# 32768*1.414  approximation maxi for chroma
			valc = (double)i / 47999.0;
			dCcurve[i] = CLIPD(valc);
		}
	}
	//end preparate histogram
	int width = lab->W, height = lab->H;
	float minQ=10000.f;
	float minM=10000.f;
	float maxQ= -1000.f;
	float maxM= -1000.f;
	float w_h;
	float a_w;
	float c_;
	float f_l;
	double Yw;
	Yw=1.0;
	double Xw, Zw;
	double f,c,nc,yb,la,xw,yw,zw,f2,c2,nc2,yb2,la2;
	double fl,n,nbb,ncb,aw;
	double xwd,ywd,zwd;
	int alg=0;
	bool algepd=false;	
	float sum=0.f;

	bool ciedata=params->colorappearance.datacie;

	ColorTemp::temp2mulxyz (params->wb.temperature, params->wb.green, params->wb.method, Xw, Zw); //compute white Xw Yw Zw  : white current WB
	//viewing condition for surround
	if(params->colorappearance.surround=="Average") { f  = 1.00; c  = 0.69; nc = 1.00;f2=1.0,c2=0.69,nc2=1.0;}
	else if(params->colorappearance.surround=="Dim"){ f2  = 0.9; c2  = 0.59; nc2 = 0.9;f=1.0,c=0.69,nc=1.0;}
	else if(params->colorappearance.surround=="Dark"){f2  = 0.8; c2  = 0.525;nc2 = 0.8;f=1.0,c=0.69,nc=1.0;}
	else if(params->colorappearance.surround=="ExtremelyDark"){f2  = 0.8; c2  = 0.41;nc2 = 0.8;f=1.0,c=0.69,nc=1.0;}

	//scene condition for surround
	if(params->colorappearance.surrsource==true)  {f  = 0.85; c  = 0.55; nc = 0.85;}// if user => source image has surround very dark
	//with which algorithme
	if     (params->colorappearance.algo=="JC")  alg=0;
	else if(params->colorappearance.algo=="JS")  alg=1;
	else if(params->colorappearance.algo=="QM")  {alg=2;algepd=true;}
	else if(params->colorappearance.algo=="ALL")  {alg=3;algepd=true;}

	bool needJ = (alg==0 || alg==1 || alg==3);
	bool needQ = (alg==2 || alg==3);
	//settings white point of output device - or illuminant viewing
	if(settings->viewingdevice==0) {xwd=96.42;ywd=100.0;zwd=82.52;}//5000K
	else if(settings->viewingdevice==1) {xwd=95.68;ywd=100.0;zwd=92.15;}//5500
	else if(settings->viewingdevice==2) {xwd=95.24;ywd=100.0;zwd=100.81;}//6000
	else if(settings->viewingdevice==3)  {xwd=95.04;ywd=100.0;zwd=108.88;}//6500
	else if(settings->viewingdevice==4)  {xwd=109.85;ywd=100.0;zwd=35.58;}//tungsten
	else if(settings->viewingdevice==5)  {xwd=99.18;ywd=100.0;zwd=67.39;}//fluo F2
	else if(settings->viewingdevice==6)  {xwd=95.04;ywd=100.0;zwd=108.75;}//fluo F7
	else if(settings->viewingdevice==7)  {xwd=100.96;ywd=100.0;zwd=64.35;}//fluo F11


	//settings mean Luminance Y of output device or viewing
	if(settings->viewingdevicegrey==0) {yb2=5.0;}
	else if(settings->viewingdevicegrey==1) {yb2=10.0;}
	else if(settings->viewingdevicegrey==2) {yb2=15.0;}
	else if(settings->viewingdevicegrey==3) {yb2=18.0;}
	else if(settings->viewingdevicegrey==4) {yb2=23.0;}
	else if(settings->viewingdevicegrey==5)  {yb2=30.0;}
	else if(settings->viewingdevicegrey==6)  {yb2=40.0;}

	//La and la2 = ambiant luminosity scene and viewing
	la=double(params->colorappearance.adapscen);
	if(pwb==2){
	if(params->colorappearance.autoadapscen) la=adap;
	}
	
	la2=double(params->colorappearance.adaplum);

	// level of adaptation
	double deg=(params->colorappearance.degree)/100.0;
	double pilot=params->colorappearance.autodegree ? 2.0 : deg;

	//algoritm's params
	float jli=params->colorappearance.jlight;
	float chr=params->colorappearance.chroma;
	float contra=params->colorappearance.contrast;
	float qbri=params->colorappearance.qbright;
	float schr=params->colorappearance.schroma;
	float mchr=params->colorappearance.mchroma;
	float qcontra=params->colorappearance.qcontrast;
	float hue=params->colorappearance.colorh;
	double rstprotection = 100.-params->colorappearance.rstprotection;
	if(schr>0.0) schr=schr/2.0f;//divide sensibility for saturation

    // extracting datas from 'params' to avoid cache flush (to be confirmed)
    ColorAppearanceParams::eTCModeId curveMode = params->colorappearance.curveMode;
    ColorAppearanceParams::eTCModeId curveMode2 = params->colorappearance.curveMode2;
    bool hasColCurve1 = bool(customColCurve1);
    bool hasColCurve2 = bool(customColCurve2);
    ColorAppearanceParams::eCTCModeId curveMode3 = params->colorappearance.curveMode3;
    bool hasColCurve3 = bool(customColCurve3);


	if(CAMBrightCurveJ.dirty || CAMBrightCurveQ.dirty){
		LUTu hist16J;
		LUTu hist16Q;
		if (needJ) {
			hist16J (65536);
			hist16J.clear();
		}
		if (needQ) {
			hist16Q (65536);
			hist16Q.clear();
		}
		float koef=1.0f;//rough correspondence between L and J
		for (int i=0; i<height; i++)
	 //   for (int i=begh; i<endh; i++)
			for (int j=0; j<width; j++) {//rough correspondence between L and J
				float currL = lab->L[i][j]/327.68f;
				if     (currL>95.) koef=1.f;
				else if(currL>85.) koef=0.97f;
				else if(currL>80.) koef=0.93f;
				else if(currL>70.) koef=0.87f;
				else if(currL>60.) koef=0.85f;
				else if(currL>50.) koef=0.8f;
				else if(currL>40.) koef=0.75f;
				else if(currL>30.) koef=0.7f;
				else if(currL>20.) koef=0.7f;
				else if(currL>10.) koef=0.9f;
				else if(currL>0.) koef=1.0f;

				if (needJ)
					hist16J[CLIP((int)((koef*lab->L[i][j])))]++;//evaluate histogram luminance L # J
				if (needQ)
					hist16Q[CLIP((int) (32768.f*sqrt((koef*(lab->L[i][j]))/32768.f)))]++;	//for brightness Q : approximation for Q=wh*sqrt(J/100)  J not equal L
				sum+=koef*lab->L[i][j];//evaluate mean J to calcualte Yb
		}
		//mean=(sum/((endh-begh)*width))/327.68f;//for Yb  for all image...if one day "pipette" we can adapt Yb for each zone
		mean=(sum/((height)*width))/327.68f;//for Yb  for all image...if one day "pipette" we can adapt Yb for each zone

		//evaluate lightness, contrast
		if (needJ) {
			if (!CAMBrightCurveJ) {
				CAMBrightCurveJ(65536,0);
				CAMBrightCurveJ.dirty = false;
			}
			ColorTemp::curveJ (jli, contra, 1, CAMBrightCurveJ, hist16J);//lightness and contrast J
		}
		if (needQ) {
			if (!CAMBrightCurveQ) {
				CAMBrightCurveQ(65536,0);
				CAMBrightCurveQ.dirty = false;
			}
			ColorTemp::curveJ (qbri, qcontra, 1, CAMBrightCurveQ, hist16Q);//brightness and contrast Q
		}
	}
	if     (mean<15.f) yb=3.0;
	else if(mean<30.f) yb=5.0;
	else if(mean<40.f) yb=10.0;
	else if(mean<45.f) yb=15.0;
	else if(mean<50.f) yb=18.0;
	else if(mean<55.f) yb=23.0;
	else if(mean<60.f) yb=30.0;
	else if(mean<70.f) yb=40.0;
	else if(mean<80.f) yb=60.0;
	else if(mean<90.f) yb=80.0;
	else               yb=90.0;

	int gamu=0;
	bool highlight = params->toneCurve.hrenabled; //Get the value if "highlight reconstruction" is activated

	if(params->colorappearance.gamut==true) gamu=1;//enabled gamut control
	xw=100.0*Xw;
	yw=100.0*Yw;
	zw=100.0*Zw;
	double xw1,yw1,zw1,xw2,yw2,zw2;
	// settings of WB: scene and viewing
    if(params->colorappearance.wbmodel=="RawT") {xw1=96.46;yw1=100.0;zw1=82.445;xw2=xwd;yw2=ywd;zw2=zwd;}	//use RT WB; CAT 02 is used for output device (see prefreneces)
    else if(params->colorappearance.wbmodel=="RawTCAT02") {xw1=xw;yw1=yw;zw1=zw;xw2=xwd;yw2=ywd;zw2=zwd;}	// Settings RT WB are used for CAT02 => mix , CAT02 is use for output device (screen: D50 D65, projector: lamp, LED) see preferences
	double cz,wh, pfl;
	ColorTemp::initcam1(gamu, yb, pilot, f, la,xw, yw, zw, n, d, nbb, ncb,cz, aw, wh, pfl, fl, c);
	double nj,dj,nbbj,ncbj,czj,awj,flj;
	ColorTemp::initcam2(gamu,yb2, f2,  la2,  xw2,  yw2,  zw2, nj, dj, nbbj, ncbj,czj, awj, flj);




#ifndef _DEBUG
#pragma omp parallel default(shared) firstprivate(lab,xw1,xw2,yw1,yw2,zw1,zw2,pilot,jli,chr,yb,la,yb2,la2,fl,nc,f,c, height,width,begh, endh,nc2,f2,c2, alg,algepd, gamu, highlight, rstprotection, pW)
#endif
{	//matrix for current working space
	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};

#ifndef _DEBUG
#pragma omp for schedule(dynamic, 10)
#endif
	for (int i=0; i<height; i++)
//	for (int i=begh; i<endh; i++)
		for (int j=0; j<width; j++) {

			float L=lab->L[i][j];
			float a=lab->a[i][j];
			float b=lab->b[i][j];
			float x1,y1,z1;
			double x,y,z;
			double epsil=0.0001;
			//convert Lab => XYZ
			Color::Lab2XYZ(L, a, b, x1, y1, z1);
		//	double J, C, h, Q, M, s, aw, fl, wh;
			double J, C, h, Q, M, s;

			double Jpro,Cpro, hpro, Qpro, Mpro, spro;
			bool t1L=false;
			bool t1B=false;
			bool t2B=false;
			int c1s=0;
			int c1co=0;
			//double n,nbb,ncb,pfl,cz,d;
			x=(double)x1/655.35;
			y=(double)y1/655.35;
			z=(double)z1/655.35;
			//process source==> normal
			ColorTemp::xyz2jchqms_ciecam02( J, C,  h,
                           Q,  M,  s, aw, fl, wh,
                           x,  y,  z,
                           xw1, yw1,  zw1,
                           yb,  la,
                           f, c,  nc,  pilot, gamu , n, nbb, ncb, pfl, cz, d );
			Jpro=J;
			Cpro=C;
			hpro=h;
			Qpro=Q;
			Mpro=M;
			spro=s;
			w_h=wh+epsil;
			a_w=aw;
			c_=c;
			f_l=fl;
			// we cannot have all algoritms with all chroma curves
			if(alg==1) 	{
				// Lightness saturation
				if(Jpro > 99.9f)
					Jpro = 99.9f;
				Jpro=(CAMBrightCurveJ[(float)(Jpro*327.68)])/327.68;//ligthness CIECAM02 + contrast
				double sres;
				double Sp=spro/100.0;
				double parsat=1.5;
				parsat=1.5;//parsat=1.5 =>saturation  ; 1.8 => chroma ; 2.5 => colorfullness (personal evaluation)
				if(schr==-100.0) schr=-99.8;
				ColorTemp::curvecolor(schr, Sp , sres, parsat);
				double coe=pow(fl,0.25);
				float dred=100.f;// in C mode
				float protect_red=80.0f; // in C mode
				dred = 100.0 * sqrt((dred*coe)/Qpro);
				protect_red=100.0 * sqrt((protect_red*coe)/Qpro);
				int sk=0;
				float ko=100.f;
				Color::skinred(Jpro, hpro, sres, Sp, dred, protect_red,sk,rstprotection,ko, spro);
				Qpro= ( 4.0 / c ) * sqrt( Jpro / 100.0 ) * ( aw + 4.0 ) ;
				Cpro=(spro*spro*Qpro)/(10000.0);
				}
			else if(alg==3 || alg==0  || alg==2) {
				double coef=32760./wh;
				if(alg==3 || alg==2) {
					if(Qpro*coef > 32767.0f)
						Qpro=(CAMBrightCurveQ[(float)32767.0f])/coef;//brightness and contrast
					else
						Qpro=(CAMBrightCurveQ[(float)(Qpro*coef)])/coef;//brightness and contrast
				}
				double Mp, sres;
				double coe=pow(fl,0.25);
				Mp=Mpro/100.0;
				double parsat=2.5;
				if(mchr==-100.0) mchr=-99.8 ;
				if(mchr==100.0) mchr=99.9;
				if(alg==3 || alg==2) ColorTemp::curvecolor(mchr, Mp , sres, parsat); else ColorTemp::curvecolor(0.0, Mp , sres, parsat);//colorfullness
				float dred=100.f;//in C mode
				float protect_red=80.0f;// in C mode
				dred *=coe;//in M mode
				protect_red	*=coe;//M mode
				int sk=0;
				float ko=100.f;
				Color::skinred(Jpro, hpro, sres, Mp, dred, protect_red,sk,rstprotection,ko, Mpro);
				Jpro=(100.0* Qpro*Qpro) /(wh*wh);
				Cpro= Mpro/coe;
				spro = 100.0 * sqrt( Mpro / Qpro );
				if(alg!=2) {
					if(Jpro > 99.9f)
						Jpro = 99.9f;
					Jpro=(CAMBrightCurveJ[(float)(Jpro*327.68f)])/327.68f;//ligthness CIECAM02 + contrast
				}
				double Cp;
				double Sp=spro/100.0;
				parsat=1.5;
				if(schr==-100.0) schr=-99.;
				if(schr==100.0) schr=98.;
				if(alg==3) ColorTemp::curvecolor(schr, Sp , sres, parsat);	else ColorTemp::curvecolor(0.0, Sp , sres, parsat);	//saturation
				dred=100.f;// in C mode
				protect_red=80.0f; // in C mode
				dred = 100.0 * sqrt((dred*coe)/Q);
				protect_red=100.0 * sqrt((protect_red*coe)/Q);
				sk=0;
				Color::skinred(Jpro, hpro, sres, Sp, dred, protect_red,sk,rstprotection,ko, spro);
				//double Q1;
				Qpro= ( 4.0 / c ) * sqrt( Jpro / 100.0 ) * ( aw + 4.0 ) ;
				Cpro=(spro*spro*Qpro)/(10000.0);
				Cp=Cpro/100.0;
				parsat=1.8;//parsat=1.5 =>saturation  ; 1.8 => chroma ; 2.5 => colorfullness (personal evaluation : for not)
				if(chr==-100.0) chr=-99.8;
				if(alg!=2) ColorTemp::curvecolor(chr, Cp , sres, parsat);else ColorTemp::curvecolor(0.0, Cp , sres, parsat);	//chroma
				dred=55.f;
				protect_red=30.0f;
				sk=1;
				Color::skinred(Jpro, hpro, sres, Cp, dred, protect_red,sk,rstprotection, ko, Cpro);

				hpro=hpro+hue;if( hpro < 0.0 ) hpro += 360.0;//hue
			}

	 if (hasColCurve1) {//curve 1 with Lightness and Brightness
		if (curveMode==ColorAppearanceParams::TC_MODE_LIGHT){
		  /*  float Jj=(float) Jpro*327.68;
			float Jold=Jj;
			const Lightcurve& userColCurve = static_cast<const Lightcurve&>(customColCurve1);
				userColCurve.Apply(Jj);
				Jj=0.7f*(Jj-Jold)+Jold;//divide sensibility
			*/	
		    float Jj=(float) Jpro*327.68f;
			float Jold=Jj;
			float Jold100=(float) Jpro;
			float redu=25.f;
			float reduc=1.f;		
			const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&>(customColCurve1);
				userColCurveJ1.Apply(Jj);
				if(Jj>Jold)	{
					if(Jj<65535.f)	{
							if(Jold < 327.68f*redu) Jj=0.3f*(Jj-Jold)+Jold;//divide sensibility
							else 		{
										reduc=LIM((100.f-Jold100)/(100.f-redu),0.f,1.f);
										Jj=0.3f*reduc*(Jj-Jold)+Jold;//reduct sensibility in highlights
										}
									}
							}
				else if(Jj>10.f) Jj=0.8f*(Jj-Jold)+Jold;
				else if (Jj>=0.f) Jj=0.90f*(Jj-Jold)+Jold;// not zero ==>artifacts	
				
				
			Jpro=(double)(Jj/327.68f);
			t1L=true;
		}
	else if (curveMode==ColorAppearanceParams::TC_MODE_BRIGHT){
			//attention! Brightness curves are open - unlike Lightness or Lab or RGB==> rendering  and algoritms will be different
			float coef=((aw+4.f)*(4.f/c))/100.f;
			float Qq=(float) Qpro*327.68f*(1.f/coef);
			float Qold100=(float) Qpro/coef;
			
			float Qold=Qq;
			float redu=20.f;
			float reduc=1.f;		
			
			const Brightcurve& userColCurveB1 = static_cast<const Brightcurve&>(customColCurve1);
					userColCurveB1.Apply(Qq);
				if(Qq>Qold)	{
					if(Qq<65535.f)	{
						if(Qold < 327.68f*redu) Qq=0.25f*(Qq-Qold)+Qold;//divide sensibility
						else 			{
										reduc=LIM((100.f-Qold100)/(100.f-redu),0.f,1.f);
										Qq=0.25f*reduc*(Qq-Qold)+Qold;//reduct sensibility in highlights
										}
									}
							}
				else if(Qq>10.f) Qq=0.5f*(Qq-Qold)+Qold;
				else if (Qq>=0.f) Qq=0.7f*(Qq-Qold)+Qold;// not zero ==>artifacts
			Qpro=(double)(Qq*(coef)/327.68f);
			Jpro=100.*(Qpro*Qpro)/((4.0/c)*(4.0/c)*(aw+4.0)*(aw+4.0));
			t1B=true;		
		}
	}

	if (hasColCurve2) {//curve 2 with Lightness and Brightness
		if (curveMode2==ColorAppearanceParams::TC_MODE_LIGHT){
			float Jj=(float) Jpro*327.68;
			float Jold=Jj;
			/*
			const Lightcurve& userColCurve = static_cast<const Lightcurve&>(customColCurve2);
					userColCurve.Apply(Jj);
					Jj=0.7f*(Jj-Jold)+Jold;//divide sensibility
					*/
			float Jold100=(float) Jpro;
			float redu=25.f;
			float reduc=1.f;					
			const Lightcurve& userColCurveJ2 = static_cast<const Lightcurve&>(customColCurve2);
					userColCurveJ2.Apply(Jj);
				if(Jj>Jold)	{
					if(Jj<65535.f)	{
							if(Jold < 327.68f*redu) Jj=0.3f*(Jj-Jold)+Jold;//divide sensibility
							else 		{
										reduc=LIM((100.f-Jold100)/(100.f-redu),0.f,1.f);
										Jj=0.3f*reduc*(Jj-Jold)+Jold;//reduct sensibility in highlights
										}
									}
							}
				else if(Jj>10.f) {if(!t1L)Jj=0.8f*(Jj-Jold)+Jold;else Jj=0.4f*(Jj-Jold)+Jold;}
				else if (Jj>=0.f){if(!t1L)Jj=0.90f*(Jj-Jold)+Jold;else Jj=0.5f*(Jj-Jold)+Jold;}// not zero ==>artifacts	
					
			Jpro=(double)(Jj/327.68f);
		}
	else if (curveMode2==ColorAppearanceParams::TC_MODE_BRIGHT){ //
			float coef=((aw+4.f)*(4.f/c))/100.f;
			float Qq=(float) Qpro*327.68f*(1.f/coef);
			float Qold100=(float) Qpro/coef;
			
			float Qold=Qq;
			float redu=20.f;
			float reduc=1.f;		
			
			const Brightcurve& userColCurveB2 = static_cast<const Brightcurve&>(customColCurve2);
					userColCurveB2.Apply(Qq);
				if(Qq>Qold)	{
					if(Qq<65535.f)	{
						if(Qold < 327.68f*redu) Qq=0.25f*(Qq-Qold)+Qold;//divide sensibility
						else 			{
										reduc=LIM((100.f-Qold100)/(100.f-redu),0.f,1.f);
										Qq=0.25f*reduc*(Qq-Qold)+Qold;//reduct sensibility in highlights
										}
									}
							}
				else if(Qq>10.f) Qq=0.5f*(Qq-Qold)+Qold;
				else if (Qq>=0.f) Qq=0.7f*(Qq-Qold)+Qold;// not zero ==>artifacts
			Qpro=(double)(Qq*(coef)/327.68f);
			Jpro=100.*(Qpro*Qpro)/((4.0/c)*(4.0/c)*(aw+4.0)*(aw+4.0));
			t2B=true;
			
			if(t1L){//to workaround the problem if we modify curve1-lightnees after curve2 brightness(the cat that bites its own tail!) in fact it's another type of curve only for this case
			coef=2.f;//adapt Q to J approximation
			Qq=(float) Qpro*coef;
			Qold=Qq;
			const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&>(customColCurve1);
				userColCurveJ1.Apply(Qq);
				Qq=0.1f*(Qq-Qold)+Qold;//approximative adaptation	
			Qpro=(double)(Qq/coef);
			Jpro=100.*(Qpro*Qpro)/((4.0/c)*(4.0/c)*(aw+4.0)*(aw+4.0));			
			} 
			
			}
	}

	if (hasColCurve3) {//curve 3 with chroma saturation colorfullness
		if (curveMode3==ColorAppearanceParams::TC_MODE_CHROMA){
		    double parsat=0.8;//0.68;
			double coef=327.68/parsat;
			float Cc=(float) Cpro*coef;
			float Ccold=Cc;
			const Chromacurve& userColCurve = static_cast<const Chromacurve&>(customColCurve3);
				userColCurve.Apply(Cc);
				float dred=55.f;
				float protect_red=30.0f;
				float sk=1;
				float ko=1.f/coef;
				Color::skinred(Jpro, hpro, Cc, Ccold, dred, protect_red,sk,rstprotection,ko, Cpro);
			//	Cpro=Cc/coef;
		}
	else if (curveMode3==ColorAppearanceParams::TC_MODE_SATUR){ //
				double parsat=0.8;//0.6
				double coef=327.68/parsat;
				float Ss=(float) spro*coef;
				float Sold=Ss;
				const Saturcurve& userColCurve = static_cast<const Saturcurve&>(customColCurve3);
					userColCurve.Apply(Ss);
					Ss=0.6f*(Ss-Sold)+Sold;//divide sensibility saturation
				double coe=pow(fl,0.25);
				float dred=100.f;// in C mode
				float protect_red=80.0f; // in C mode
				dred = 100.0 * sqrt((dred*coe)/Qpro);
				protect_red=100.0 * sqrt((protect_red*coe)/Qpro);
				int sk=0;
				float ko=1.f/coef;
				Color::skinred(Jpro, hpro, Ss, Sold, dred, protect_red,sk,rstprotection,ko, spro);
				Qpro= ( 4.0 / c ) * sqrt( Jpro / 100.0 ) * ( aw + 4.0 ) ;
				Cpro=(spro*spro*Qpro)/(10000.0);
				c1s=1;

			}
	else if (curveMode3==ColorAppearanceParams::TC_MODE_COLORF){ //
				double parsat=0.8;//0.68;
				double coef=327.68/parsat;
				float Mm=(float) Mpro*coef;
				float Mold=Mm;
				const Colorfcurve& userColCurve = static_cast<const Colorfcurve&>(customColCurve3);
					userColCurve.Apply(Mm);
				double coe=pow(fl,0.25);
				float dred=100.f;//in C mode
				float protect_red=80.0f;// in C mode
				dred *=coe;//in M mode
				protect_red	*=coe;
				int sk=0;
				float ko=1.f/coef;
				Color::skinred(Jpro, hpro, Mm, Mold, dred, protect_red,sk,rstprotection,ko, Mpro);
				Cpro= Mpro/coe;
				c1co=1;
			}
	}
			//to retrieve the correct values of variables
			if(t2B && t1B) Jpro=(100.0* Qpro*Qpro) /(wh*wh);// for brightness curve
			if(c1s==1) {
				Qpro= ( 4.0 / c ) * sqrt( Jpro / 100.0 ) * ( aw + 4.0 ) ;//for saturation curve
				Cpro=(spro*spro*Qpro)/(10000.0);
				}
			if(c1co==1) {	double coe=pow(fl,0.25);Cpro= Mpro/coe;}	// for colorfullness curve
			//retrieve values C,J...s
			C=Cpro;
			J=Jpro;
			Q=Qpro;
			M=Mpro;
			h=hpro;
			s=spro;

		if(params->colorappearance.tonecie  || settings->autocielab){//use pointer for tonemapping with CIECAM and also sharpening , defringe, contrast detail
	//	if(params->colorappearance.tonecie  || params->colorappearance.sharpcie){//use pointer for tonemapping with CIECAM and also sharpening , defringe, contrast detail
			float Qred= ( 4.0 / c)  * ( aw + 4.0 );//estimate Q max if J=100.0
			ncie->Q_p[i][j]=(float)Q+epsil;//epsil to avoid Q=0
			ncie->M_p[i][j]=(float)M+epsil;
			ncie->J_p[i][j]=(float)J+epsil;
			ncie->h_p[i][j]=(float)h;
			ncie->C_p[i][j]=(float)C+epsil;
	//		ncie->s_p[i][j]=(float)s;
			ncie->sh_p[i][j]=(float) 32768.*(( 4.0 / c )*sqrt( J / 100.0 ) * ( aw + 4.0 ))/Qred ;
	//		ncie->ch_p[i][j]=(float) 327.68*C;
			if(ncie->Q_p[i][j]<minQ) minQ=ncie->Q_p[i][j];//minima
			if(ncie->Q_p[i][j]>maxQ) maxQ=ncie->Q_p[i][j];//maxima
			}
			if(!params->colorappearance.tonecie  || !settings->autocielab  || !params->edgePreservingDecompositionUI.enabled ){						
			
//			if(!params->edgePreservingDecompositionUI.enabled || !params->colorappearance.tonecie  || !settings->autocielab){
		//	if(!params->edgePreservingDecompositionUI.enabled || !params->colorappearance.tonecie  || !params->colorappearance.sharpcie){
			int posl, posc;
			double brli=327.;
			double chsacol=327.;
			int libr=0;
			int colch=0;
			if(curveMode==ColorAppearanceParams::TC_MODE_BRIGHT) {brli=70.0; libr=1;}
			else if(curveMode==ColorAppearanceParams::TC_MODE_LIGHT) {brli=327.;libr=0;}
			if (curveMode3==ColorAppearanceParams::TC_MODE_CHROMA) {chsacol=327.;colch=0;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_SATUR) {chsacol=450.0;colch=1;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_COLORF) {chsacol=327.0;colch=2;}
			if(ciedata) {
			// Data for J Q M s and C histograms
				//update histogram
				jp=true;
                if(pW!=1){//only with improccoordinator
				if(libr==1) posl=CLIP((int)(Q*brli));//40.0 to 100.0 approximative factor for Q  - 327 for J
				else if(libr==0) posl=CLIP((int)(J*brli));//327 for J
				hist16JCAM[posl]++;
				}
				chropC=true;
                if(pW!=1){//only with improccoordinator
				if(colch==0) posc=CLIP((int)(C*chsacol));//450.0 approximative factor for s    320 for M
				else if(colch==1) posc=CLIP((int)(s*chsacol));
				else if(colch==2) posc=CLIP((int)(M*chsacol));
				hist16_CCAM[posc]++;
				}
			}
			double xx,yy,zz;
			//double nj, nbbj, ncbj, flj, czj, dj, awj;
			//process normal==> viewing
			ColorTemp::jch2xyz_ciecam02( xx, yy, zz,
			                             J,  C, h,
			                             xw2, yw2,  zw2,
			                             yb2, la2,
			                             f2,  c2, nc2, gamu, nj, nbbj, ncbj, flj, czj, dj, awj);
			x=(float)xx*655.35;
			y=(float)yy*655.35;
			z=(float)zz*655.35;
			float Ll,aa,bb;
			//convert xyz=>lab
			Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
			lab->L[i][j]=Ll;
			lab->a[i][j]=aa;
			lab->b[i][j]=bb;
		// gamut control in Lab mode; I must study how to do with cIECAM only
		if(gamu==1) {
					float R,G,B;
					float HH, Lprov1, Chprov1;
					Lprov1=lab->L[i][j]/327.68f;
					Chprov1=sqrt(SQR(lab->a[i][j]/327.68f) + SQR(lab->b[i][j]/327.68f));
					HH=atan2(lab->b[i][j],lab->a[i][j]);

#ifdef _DEBUG
					bool neg=false;
					bool more_rgb=false;
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif

					lab->L[i][j]=Lprov1*327.68f;
					lab->a[i][j]=327.68f*Chprov1*cos(HH);
					lab->b[i][j]=327.68f*Chprov1*sin(HH);

		}
		}
		}
	}
	// End of parallelization
if(!params->edgePreservingDecompositionUI.enabled || !params->colorappearance.tonecie   || !settings->autocielab){//normal
//if(!params->edgePreservingDecompositionUI.enabled || !params->colorappearance.tonecie   || !params->colorappearance.sharpcie){//normal

	if(ciedata) {
    //update histogram J
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=32768; i++) {//
			if (jp) {
				float hval = dLcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				histLCAM[hi] += hist16JCAM[i] ;
			}
		}
	}
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=48000; i++) {//
			if (chropC) {
				float hvalc = dCcurve[i];
				int hic = (int)(255.0*CLIPD(hvalc)); //
				histCCAM[hic] += hist16_CCAM[i] ;
			}
		}
	}
	}
}
#ifdef _DEBUG
	if (settings->verbose) {
		t2e.set();
		printf("CIECAM02 performed in %d usec:\n", t2e.etime(t1e));
		//	printf("minc=%f maxc=%f minj=%f maxj=%f\n",minc,maxc,minj,maxj);
	}
#endif

if(settings->autocielab) {
//if(params->colorappearance.sharpcie) {

//all this treatments reduce artifacts, but can lead to slightly  different results
if(params->defringe.enabled) if(execsharp) ImProcFunctions::defringecam (ncie);// 

if(params->colorappearance.badpixsl > 0) if(execsharp){ int mode=params->colorappearance.badpixsl;
											ImProcFunctions::badpixcam (ncie, 3.4, 5, mode);//for bad pixels
										}	


if (params->sharpenMicro.enabled)if(execsharp) ImProcFunctions::MLmicrocontrastcam(ncie);

if(params->sharpening.enabled)  if(execsharp) {ImProcFunctions::sharpeningcam (ncie, (float**)buffer);}			//sharpening adapted to CIECAM 

if(params->dirpyrequalizer.enabled) if(execsharp) dirpyr_equalizercam(ncie, ncie->sh_p, ncie->sh_p, ncie->W, ncie->H, params->dirpyrequalizer.mult, params->dirpyrequalizer.threshold, true);//contrast by detail adapted to CIECAM

		   float Qredi= ( 4.0 / c_)  * ( a_w + 4.0 );
		   float co_e=(pow(f_l,0.25f));

#ifndef _DEBUG	
#pragma omp parallel default(shared) firstprivate(height,width, Qredi,a_w,c_)
#endif
{	
#ifndef _DEBUG
		#pragma omp for schedule(dynamic, 10)
#endif
		   		for (int i=0; i<height; i++) // update CieImages with new values after sharpening, defringe, contrast by detail level
					for (int j=0; j<width; j++) {
						float interm=Qredi*ncie->sh_p[i][j]/(32768.f);
						ncie->J_p[i][j]=100.0* interm*interm/((a_w+4.)*(a_w+4.)*(4./c_)*(4./c_));
						ncie->Q_p[i][j]=( 4.0 / c_)  * ( a_w + 4.0 ) *  sqrt(ncie->J_p[i][j]/100.f);
						ncie->M_p[i][j]=ncie->C_p[i][j]*co_e;
					}
			}
}	

if((params->colorappearance.tonecie || (params->colorappearance.tonecie && params->edgePreservingDecompositionUI.enabled)) || (params->sharpening.enabled && settings->autocielab)
		|| (params->dirpyrequalizer.enabled && settings->autocielab) ||(params->defringe.enabled && settings->autocielab)  || (params->sharpenMicro.enabled && settings->autocielab)
		||  (params->colorappearance.badpixsl > 0 && settings->autocielab)) {
		
		if(params->edgePreservingDecompositionUI.enabled  && params->colorappearance.tonecie  && algepd) ImProcFunctions::EPDToneMapCIE(ncie, a_w, c_, w_h, width, height, begh, endh, minQ, maxQ, Iterates, scale );
			//EPDToneMapCIE adapted to CIECAM

	
#ifndef _DEBUG	
#pragma omp parallel default(shared) firstprivate(lab,xw2,yw2,zw2,chr,yb,la2,yb2, height,width,begh, endh, nc2,f2,c2, gamu, highlight,pW)
#endif
{	
	TMatrix wiprofa = iccStore->workingSpaceInverseMatrix (params->icm.working);
	double wipa[3][3] = {
		{wiprofa[0][0],wiprofa[0][1],wiprofa[0][2]},
		{wiprofa[1][0],wiprofa[1][1],wiprofa[1][2]},
		{wiprofa[2][0],wiprofa[2][1],wiprofa[2][2]}
	};
	
	
#ifndef _DEBUG
		#pragma omp for schedule(dynamic, 10)
#endif
		for (int i=0; i<height; i++) // update CIECAM with new values after tone-mapping
	//	for (int i=begh; i<endh; i++) 
			for (int j=0; j<width; j++) {
			double xx,yy,zz;
			float x,y,z;
			const float eps=0.0001;
			float co_e=(pow(f_l,0.25f))+eps;
	//		if(params->edgePreservingDecompositionUI.enabled) ncie->J_p[i][j]=(100.0* ncie->Q_p[i][j]*ncie->Q_p[i][j])/(w_h*w_h);
			if(params->edgePreservingDecompositionUI.enabled) ncie->J_p[i][j]=(100.0* ncie->Q_p[i][j]*ncie->Q_p[i][j])/SQR((4./c)*(aw+4.));
			
			ncie->C_p[i][j]	=(ncie->M_p[i][j])/co_e;
			//show histogram in CIECAM mode (Q,J, M,s,C)
			int posl, posc;
			double brli=327.;
			double chsacol=327.;
			int libr=0;
			int colch=0;
			float sa_t;
			if(curveMode==ColorAppearanceParams::TC_MODE_BRIGHT) {brli=70.0; libr=1;}
			else if(curveMode==ColorAppearanceParams::TC_MODE_LIGHT) {brli=327.;libr=0;}
			if (curveMode3==ColorAppearanceParams::TC_MODE_CHROMA) {chsacol=327.;colch=0;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_SATUR) {chsacol=450.0;colch=1;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_COLORF) {chsacol=327.0;colch=2;}
			if(ciedata) {
			// Data for J Q M s and C histograms
				//update histogram
				jp=true;
                if(pW!=1){//only with improccoordinator
				if(libr==1) posl=CLIP((int)(ncie->Q_p[i][j]*brli));//40.0 to 100.0 approximative factor for Q  - 327 for J
				else if(libr==0) posl=CLIP((int)(ncie->J_p[i][j]*brli));//327 for J
				hist16JCAM[posl]++;
				}
				chropC=true;
                if(pW!=1){//only with improccoordinator
			if(colch==0) posc=CLIP((int)(ncie->C_p[i][j]*chsacol));//450.0 approximative factor for s    320 for M
				else if(colch==1) {sa_t=100.f*sqrt(ncie->C_p[i][j]/ncie->Q_p[i][j]); posc=CLIP((int)(sa_t*chsacol));}//Q_p always > 0
				else if(colch==2) posc=CLIP((int)(ncie->M_p[i][j]*chsacol));
				hist16_CCAM[posc]++;
				}
			}
			//end histograms
		//	double nd, nbbd, ncbd, fld, czd, dd, awd;
			ColorTemp::jch2xyz_ciecam02( xx, yy, zz,
			                             ncie->J_p[i][j],  ncie->C_p[i][j], ncie->h_p[i][j],
			                             xw2, yw2,  zw2,
			                             yb2, la2,
			                             f2,  c2, nc2, gamu, nj, nbbj, ncbj, flj, czj, dj, awj);
			x=(float)xx*655.35;
			y=(float)yy*655.35;
			z=(float)zz*655.35;
			float Ll,aa,bb;
			//convert xyz=>lab
			Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
			lab->L[i][j]=Ll;
			lab->a[i][j]=aa;
			lab->b[i][j]=bb;
			if(gamu==1) {
					float R,G,B;
					float HH, Lprov1, Chprov1;
					Lprov1=lab->L[i][j]/327.68f;
					Chprov1=sqrt(SQR(lab->a[i][j]/327.68f) + SQR(lab->b[i][j]/327.68f));
					HH=atan2(lab->b[i][j],lab->a[i][j]);

#ifdef _DEBUG
					bool neg=false;
					bool more_rgb=false;
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wipa, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wipa, highlight, 0.15f, 0.96f);
#endif

					lab->L[i][j]=Lprov1*327.68f;
					lab->a[i][j]=327.68f*Chprov1*cos(HH);
					lab->b[i][j]=327.68f*Chprov1*sin(HH);
				}
			}


		}
		//end parallelization
	//show CIECAM histograms
	if(ciedata) {
    //update histogram J and Q
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=32768; i++) {//
			if (jp) {
				float hval = dLcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				histLCAM[hi] += hist16JCAM[i] ;
			}
		}
	}
	//update color histogram M,s,C
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=48000; i++) {//
			if (chropC) {
				float hvalc = dCcurve[i];
				int hic = (int)(255.0*CLIPD(hvalc)); //
				histCCAM[hic] += hist16_CCAM[i] ;
			}
		}
	}
	}

	}

}
}
//end CIECAM


// Copyright (c) 2012 Jacques Desmis <jdesmis@gmail.com>
void ImProcFunctions::ciecam_02float (CieImage* ncie, float adap, int begh, int endh, int pW, int pwb, LabImage* lab, const ProcParams* params,
								const ColorAppearance & customColCurve1, const ColorAppearance & customColCurve2,const ColorAppearance & customColCurve3,
								LUTu & histLCAM, LUTu & histCCAM, LUTf & CAMBrightCurveJ, LUTf & CAMBrightCurveQ, float &mean, int Iterates, int scale, float** buffer, bool execsharp, float &d)
{
if(params->colorappearance.enabled) {
//printf("ciecam float\n");
#ifdef _DEBUG
	MyTime t1e,t2e;
	t1e.set();
#endif
	LUTf dLcurve;
	LUTu hist16JCAM;
	bool jp=false;
	float val;
	//preparate for histograms CIECAM
	if(pW!=1){//only with improccoordinator
		dLcurve(65536, 0);
		dLcurve.clear();
		hist16JCAM(65536);
		hist16JCAM.clear();
		for (int i=0; i<32768; i++) {  //# 32768*1.414  approximation maxi for chroma
			val = (double)i / 32767.0;
			dLcurve[i] = CLIPD(val);
		}
	}
	LUTf dCcurve(65536,0);
	LUTu hist16_CCAM(65536);
	bool chropC=false;
	float valc;
	
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<48000; i++) {  //# 32768*1.414  approximation maxi for chroma
			valc = (double)i / 47999.0;
			dCcurve[i] = CLIPD(valc);
		}
		hist16_CCAM.clear();
	}
	//end preparate histogram
	int width = lab->W, height = lab->H;
	float minQ=10000.f;
	float minM=10000.f;
	float maxQ= -1000.f;
	float maxM= -1000.f;
	float w_h;
	float a_w;
	float c_;
	float f_l;
	float Yw;
	Yw=1.0;
	double Xw, Zw;
	float f,nc,yb,la,c,xw,yw,zw,f2,c2,nc2,yb2,la2;
	float fl,n,nbb,ncb,aw;//d
	float xwd,ywd,zwd;
	int alg=0;
	bool algepd=false;
	float sum=0.f;

	bool ciedata=params->colorappearance.datacie;

	ColorTemp::temp2mulxyz (params->wb.temperature, params->wb.green, params->wb.method, Xw, Zw); //compute white Xw Yw Zw  : white current WB
	//viewing condition for surround
	if(params->colorappearance.surround=="Average") { f  = 1.00f; c  = 0.69f; nc = 1.00f;f2=1.0f,c2=0.69f,nc2=1.0f;}
	else if(params->colorappearance.surround=="Dim"){ f2  = 0.9f; c2  = 0.59f; nc2 = 0.9f;f=1.0f,c=0.69f,nc=1.0f;}
	else if(params->colorappearance.surround=="Dark"){f2  = 0.8f; c2  = 0.525f;nc2 = 0.8f;f=1.0f,c=0.69f,nc=1.0f;}
	else if(params->colorappearance.surround=="ExtremelyDark"){f2  = 0.8f; c2  = 0.41f;nc2 = 0.8f;f=1.0f,c=0.69f,nc=1.0f;}

	//scene condition for surround
	if(params->colorappearance.surrsource==true)  {f  = 0.85f; c  = 0.55f; nc = 0.85f;}// if user => source image has surround very dark
	//with which algorithm
	if     (params->colorappearance.algo=="JC")  alg=0;
	else if(params->colorappearance.algo=="JS")  alg=1;
	else if(params->colorappearance.algo=="QM")  {alg=2;algepd=true;}
	else if(params->colorappearance.algo=="ALL") {alg=3;algepd=true;}

	bool needJ = (alg==0 || alg==1 || alg==3);
	bool needQ = (alg==2 || alg==3);
	//settings white point of output device - or illuminant viewing
	if(settings->viewingdevice==0) {xwd=96.42f;ywd=100.0f;zwd=82.52f;}//5000K
	else if(settings->viewingdevice==1) {xwd=95.68f;ywd=100.0f;zwd=92.15f;}//5500
	else if(settings->viewingdevice==2) {xwd=95.24f;ywd=100.0f;zwd=100.81f;}//6000
	else if(settings->viewingdevice==3)  {xwd=95.04f;ywd=100.0f;zwd=108.88f;}//6500
	else if(settings->viewingdevice==4)  {xwd=109.85f;ywd=100.0f;zwd=35.58f;}//tungsten
	else if(settings->viewingdevice==5)  {xwd=99.18f;ywd=100.0f;zwd=67.39f;}//fluo F2
	else if(settings->viewingdevice==6)  {xwd=95.04f;ywd=100.0f;zwd=108.75f;}//fluo F7
	else if(settings->viewingdevice==7)  {xwd=100.96f;ywd=100.0f;zwd=64.35f;}//fluo F11


	//settings mean Luminance Y of output device or viewing
	if(settings->viewingdevicegrey==0) {yb2=5.0f;}
	else if(settings->viewingdevicegrey==1) {yb2=10.0f;}
	else if(settings->viewingdevicegrey==2) {yb2=15.0f;}
	else if(settings->viewingdevicegrey==3) {yb2=18.0f;}
	else if(settings->viewingdevicegrey==4) {yb2=23.0f;}
	else if(settings->viewingdevicegrey==5)  {yb2=30.0f;}
	else if(settings->viewingdevicegrey==6)  {yb2=40.0f;}

	//La and la2 = ambiant luminosity scene and viewing
	la=float(params->colorappearance.adapscen);
	if(pwb==2){
	if(params->colorappearance.autoadapscen) la=adap;
	}
	
	la2=float(params->colorappearance.adaplum);

	// level of adaptation
	float deg=(params->colorappearance.degree)/100.0f;
	float pilot=params->colorappearance.autodegree ? 2.0f : deg;

	//algoritm's params
	float jli=params->colorappearance.jlight;
	float chr=params->colorappearance.chroma;
	float contra=params->colorappearance.contrast;
	float qbri=params->colorappearance.qbright;
	float schr=params->colorappearance.schroma;
	float mchr=params->colorappearance.mchroma;
	float qcontra=params->colorappearance.qcontrast;
	float hue=params->colorappearance.colorh;
	float rstprotection = 100.-params->colorappearance.rstprotection;
	if(schr>0.0) schr=schr/2.0f;//divide sensibility for saturation

    // extracting datas from 'params' to avoid cache flush (to be confirmed)
    ColorAppearanceParams::eTCModeId curveMode = params->colorappearance.curveMode;
    ColorAppearanceParams::eTCModeId curveMode2 = params->colorappearance.curveMode2;
    bool hasColCurve1 = bool(customColCurve1);
    bool hasColCurve2 = bool(customColCurve2);
    ColorAppearanceParams::eCTCModeId curveMode3 = params->colorappearance.curveMode3;
    bool hasColCurve3 = bool(customColCurve3);


	if(CAMBrightCurveJ.dirty || CAMBrightCurveQ.dirty){
		LUTu hist16J;
		LUTu hist16Q;
		if (needJ) {
			hist16J (65536);
			hist16J.clear();
		}
		if (needQ) {
			hist16Q (65536);
			hist16Q.clear();
		}
		float koef=1.0f;//rough correspondence between L and J
		for (int i=0; i<height; i++)
	//	for (int i=begh; i<endh; i++)
			for (int j=0; j<width; j++) {//rough correspondence between L and J
				float currL = lab->L[i][j]/327.68f;
				if     (currL>95.) koef=1.f;
				else if(currL>85.) koef=0.97f;
				else if(currL>80.) koef=0.93f;
				else if(currL>70.) koef=0.87f;
				else if(currL>60.) koef=0.85f;
				else if(currL>50.) koef=0.8f;
				else if(currL>40.) koef=0.75f;
				else if(currL>30.) koef=0.7f;
				else if(currL>20.) koef=0.7f;
				else if(currL>10.) koef=0.9f;
				else if(currL>0.) koef=1.0f;

				if (needJ)
					hist16J[CLIP((int)((koef*lab->L[i][j])))]++;//evaluate histogram luminance L # J
				if (needQ)
					hist16Q[CLIP((int) (32768.f*sqrt((koef*(lab->L[i][j]))/32768.f)))]++;	//for brightness Q : approximation for Q=wh*sqrt(J/100)  J not equal L
				sum+=koef*lab->L[i][j];//evaluate mean J to calculate Yb
		}
		//mean=(sum/((endh-begh)*width))/327.68f;//for Yb  for all image...if one day "pipette" we can adapt Yb for each zone
		mean=(sum/((height)*width))/327.68f;//for Yb  for all image...if one day "pipette" we can adapt Yb for each zone

		//evaluate lightness, contrast
		if (needJ) {
			if (!CAMBrightCurveJ) {
				CAMBrightCurveJ(65536,0);
				CAMBrightCurveJ.dirty = false;
			}
			ColorTemp::curveJfloat (jli, contra, 1, CAMBrightCurveJ, hist16J);//lightness and contrast J
		}
		if (needQ) {
			if (!CAMBrightCurveQ) {
				CAMBrightCurveQ(65536,0);
				CAMBrightCurveQ.dirty = false;
			}
			ColorTemp::curveJfloat (qbri, qcontra, 1, CAMBrightCurveQ, hist16Q);//brightness and contrast Q
		}
	}
	if     (mean<15.f) yb=3.0f;
	else if(mean<30.f) yb=5.0f;
	else if(mean<40.f) yb=10.0f;
	else if(mean<45.f) yb=15.0f;
	else if(mean<50.f) yb=18.0f;
	else if(mean<55.f) yb=23.0f;
	else if(mean<60.f) yb=30.0f;
	else if(mean<70.f) yb=40.0f;
	else if(mean<80.f) yb=60.0f;
	else if(mean<90.f) yb=80.0f;
	else               yb=90.0f;

	int gamu=0;
	bool highlight = params->toneCurve.hrenabled; //Get the value if "highlight reconstruction" is activated

	if(params->colorappearance.gamut==true) gamu=1;//enabled gamut control
	xw=100.0f*Xw;
	yw=100.0f*Yw;
	zw=100.0f*Zw;
	float xw1,yw1,zw1,xw2,yw2,zw2;
	// settings of WB: scene and viewing
    if(params->colorappearance.wbmodel=="RawT") {xw1=96.46f;yw1=100.0f;zw1=82.445f;xw2=xwd;yw2=ywd;zw2=zwd;}	//use RT WB; CAT 02 is used for output device (see prefreneces)
    else if(params->colorappearance.wbmodel=="RawTCAT02") {xw1=xw;yw1=yw;zw1=zw;xw2=xwd;yw2=ywd;zw2=zwd;}	// Settings RT WB are used for CAT02 => mix , CAT02 is use for output device (screen: D50 D65, projector: lamp, LED) see preferences
	float cz,wh, pfl;
	ColorTemp::initcam1float(gamu, yb, pilot, f, la,xw, yw, zw, n, d, nbb, ncb,cz, aw, wh, pfl, fl, c);
	float nj,dj,nbbj,ncbj,czj,awj,flj;
	ColorTemp::initcam2float(gamu,yb2, f2,  la2,  xw2,  yw2,  zw2, nj, dj, nbbj, ncbj,czj, awj, flj);

	
#ifndef _DEBUG	
#pragma omp parallel default(shared) firstprivate(lab,xw1,xw2,yw1,yw2,zw1,zw2,pilot,jli,chr,yb,la,yb2,la2,fl,nc,f,c, height,width,begh, endh,nc2,f2,c2, alg, algepd, gamu, highlight, rstprotection, pW,nj, nbbj, ncbj, flj, czj, dj, awj, n, nbb, ncb, pfl, cz)
#endif
{	//matrix for current working space
	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};

#ifndef _DEBUG
#pragma omp for schedule(dynamic, 10)
#endif
	for (int i=0; i<height; i++)
//	for (int i=begh; i<endh; i++)
		for (int j=0; j<width; j++) {

			float L=lab->L[i][j];
			float a=lab->a[i][j];
			float b=lab->b[i][j];
			float x1,y1,z1;
			float x,y,z;
			float epsil=0.0001f;
			//convert Lab => XYZ
			Color::Lab2XYZ(L, a, b, x1, y1, z1);
			float J, C, h, Q, M, s;
			float Jp;
			float Jpro,Cpro, hpro, Qpro, Mpro, spro;
			bool t1L=false;
			bool t1B=false;
			bool t2B=false;
			int c1s=0;
			int c1co=0;

			x=(float)x1/655.35f;
			y=(float)y1/655.35f;
			z=(float)z1/655.35f;
			//process source==> normal
			ColorTemp::xyz2jchqms_ciecam02float( J, C,  h,
                           Q,  M,  s, aw, fl, wh,
                           x,  y,  z,
                           xw1, yw1,  zw1,
                           yb,  la,
                           f, c,  nc,  pilot, gamu, n, nbb, ncb, pfl, cz, d);
			Jpro=J;
			Cpro=C;
			hpro=h;
			Qpro=Q;
			Mpro=M;
			spro=s;
			w_h=wh+epsil;
			a_w=aw;
			c_=c;
			f_l=fl;
			// we cannot have all algoritms with all chroma curves
			if(alg==1) 	{
				// Lightness saturation
				if(Jpro > 99.9f)
					Jpro = 99.9f;
				Jpro=(CAMBrightCurveJ[(float)(Jpro*327.68f)])/327.68f;//ligthness CIECAM02 + contrast
				float sres;
				float Sp=spro/100.0f;
				float parsat=1.5f;
				parsat=1.5f;//parsat=1.5 =>saturation  ; 1.8 => chroma ; 2.5 => colorfullness (personal evaluation)
				if(schr==-100.0f) schr=-99.8f;
				ColorTemp::curvecolorfloat(schr, Sp , sres, parsat);
				float coe=pow_F(fl,0.25f);
				float dred=100.f;// in C mode
				float protect_red=80.0f; // in C mode
				dred = 100.0f * sqrt((dred*coe)/Qpro);
				protect_red=100.0f * sqrt((protect_red*coe)/Qpro);
				int sk=0;
				float ko=100.f;
				Color::skinredfloat(Jpro, hpro, sres, Sp, dred, protect_red,sk,rstprotection,ko, spro);
				Qpro= ( 4.0f / c ) * sqrt( Jpro / 100.0f ) * ( aw + 4.0f ) ;
				Cpro=(spro*spro*Qpro)/(10000.0f);
				}
			else if(alg==3 || alg==0  || alg==2) {
				float coef=32760.f/wh;
				if(alg==3 || alg==2) {
					if(Qpro*coef > 32767.0f)
						Qpro=(CAMBrightCurveQ[(float)32767.0f])/coef;//brightness and contrast
					else
						Qpro=(CAMBrightCurveQ[(float)(Qpro*coef)])/coef;//brightness and contrast
				}
				float Mp, sres;
				float coe=pow_F(fl,0.25f);
				Mp=Mpro/100.0f;
				float parsat=2.5f;
				if(mchr==-100.0f) mchr=-99.8f ;
				if(mchr==100.0f) mchr=99.9f;
				if(alg==3 || alg==2) ColorTemp::curvecolorfloat(mchr, Mp , sres, parsat); else ColorTemp::curvecolorfloat(0.0, Mp , sres, parsat);//colorfullness
				float dred=100.f;//in C mode
				float protect_red=80.0f;// in C mode
				dred *=coe;//in M mode
				protect_red	*=coe;//M mode
				int sk=0;
				float ko=100.f;
				Color::skinredfloat(Jpro, hpro, sres, Mp, dred, protect_red,sk,rstprotection,ko, Mpro);
				Jpro=(100.0f* Qpro*Qpro) /(wh*wh);
				Cpro= Mpro/coe;
				spro = 100.0f * sqrt( Mpro / Qpro );
				if(alg!=2) {
					if(Jpro > 99.9f)
						Jpro = 99.9f;
					Jpro=(CAMBrightCurveJ[(float)(Jpro*327.68f)])/327.68f;//ligthness CIECAM02 + contrast
				}
				float Cp;
				float Sp=spro/100.0f;
				parsat=1.5f;
				if(schr==-100.0f) schr=-99.f;
				if(schr==100.0f) schr=98.f;
				if(alg==3) ColorTemp::curvecolorfloat(schr, Sp , sres, parsat);	else ColorTemp::curvecolorfloat(0.0f, Sp , sres, parsat);	//saturation
				dred=100.f;// in C mode
				protect_red=80.0f; // in C mode
				dred = 100.0f * sqrt((dred*coe)/Q);
				protect_red=100.0f * sqrt((protect_red*coe)/Q);
				sk=0;
				Color::skinredfloat(Jpro, hpro, sres, Sp, dred, protect_red,sk,rstprotection,ko, spro);
				Qpro= ( 4.0f / c ) * sqrt( Jpro / 100.0f ) * ( aw + 4.0f ) ;
				Cpro=(spro*spro*Qpro)/(10000.0f);
				Cp=Cpro/100.0f;
				parsat=1.8f;//parsat=1.5 =>saturation  ; 1.8 => chroma ; 2.5 => colorfullness (personal evaluation : for not)
				if(chr==-100.0f) chr=-99.8f;
				if(alg!=2) ColorTemp::curvecolorfloat(chr, Cp , sres, parsat);else ColorTemp::curvecolorfloat(0.0f, Cp , sres, parsat);	//chroma
				dred=55.f;
				protect_red=30.0f;
				sk=1;
				Color::skinredfloat(Jpro, hpro, sres, Cp, dred, protect_red,sk,rstprotection, ko, Cpro);

				hpro=hpro+hue;if( hpro < 0.0f ) hpro += 360.0f;//hue
			}

	 if (hasColCurve1) {//curve 1 with Lightness and Brightness
		if (curveMode==ColorAppearanceParams::TC_MODE_LIGHT){
		    float Jj=(float) Jpro*327.68f;
			float Jold=Jj;
			float Jold100=(float) Jpro;
			float redu=25.f;
			float reduc=1.f;		
			const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&>(customColCurve1);
				userColCurveJ1.Apply(Jj);
				if(Jj>Jold)	{
					if(Jj<65535.f)	{
							if(Jold < 327.68f*redu) Jj=0.3f*(Jj-Jold)+Jold;//divide sensibility
							else 		{
										reduc=LIM((100.f-Jold100)/(100.f-redu),0.f,1.f);
										Jj=0.3f*reduc*(Jj-Jold)+Jold;//reduct sensibility in highlights
										}
									}
							}
				else if(Jj>10.f) Jj=0.8f*(Jj-Jold)+Jold;
				else if (Jj>=0.f) Jj=0.90f*(Jj-Jold)+Jold;// not zero ==>artifacts	
			Jpro=(float)(Jj/327.68f);			
			t1L=true;
		}
	else if (curveMode==ColorAppearanceParams::TC_MODE_BRIGHT){
			//attention! Brightness curves are open - unlike Lightness or Lab or RGB==> rendering  and algoritms will be different
			float coef=((aw+4.f)*(4.f/c))/100.f;
			float Qq=(float) Qpro*327.68f*(1.f/coef);
			float Qold100=(float) Qpro/coef;
			
			float Qold=Qq;
			float redu=20.f;
			float reduc=1.f;		
			
			const Brightcurve& userColCurveB1 = static_cast<const Brightcurve&>(customColCurve1);
					userColCurveB1.Apply(Qq);
				if(Qq>Qold)	{
					if(Qq<65535.f)	{
						if(Qold < 327.68f*redu) Qq=0.25f*(Qq-Qold)+Qold;//divide sensibility
						else 			{
										reduc=LIM((100.f-Qold100)/(100.f-redu),0.f,1.f);
										Qq=0.25f*reduc*(Qq-Qold)+Qold;//reduct sensibility in highlights
										}
									}
							}
				else if(Qq>10.f) Qq=0.5f*(Qq-Qold)+Qold;
				else if (Qq>=0.f) Qq=0.7f*(Qq-Qold)+Qold;// not zero ==>artifacts
			Qpro=(float)(Qq*(coef)/327.68f);
			Jpro=100.f*(Qpro*Qpro)/((4.0f/c)*(4.0f/c)*(aw+4.0f)*(aw+4.0f));
			t1B=true;		
		}
	}

	if (hasColCurve2) {//curve 2 with Lightness and Brightness
		if (curveMode2==ColorAppearanceParams::TC_MODE_LIGHT){
			float Jj=(float) Jpro*327.68f;
			float Jold=Jj;
			float Jold100=(float) Jpro;
			float redu=25.f;
			float reduc=1.f;					
			const Lightcurve& userColCurveJ2 = static_cast<const Lightcurve&>(customColCurve2);
					userColCurveJ2.Apply(Jj);
				if(Jj>Jold)	{
					if(Jj<65535.f)	{
							if(Jold < 327.68f*redu) Jj=0.3f*(Jj-Jold)+Jold;//divide sensibility
							else 		{
										reduc=LIM((100.f-Jold100)/(100.f-redu),0.f,1.f);
										Jj=0.3f*reduc*(Jj-Jold)+Jold;//reduct sensibility in highlights
										}
									}
							}
				else if(Jj>10.f) {if(!t1L)Jj=0.8f*(Jj-Jold)+Jold;else Jj=0.4f*(Jj-Jold)+Jold;}
				else if (Jj>=0.f){if(!t1L)Jj=0.90f*(Jj-Jold)+Jold;else Jj=0.5f*(Jj-Jold)+Jold;}// not zero ==>artifacts	
			Jpro=(float)(Jj/327.68f);
		}
	else if (curveMode2==ColorAppearanceParams::TC_MODE_BRIGHT){ //
			float coef=((aw+4.f)*(4.f/c))/100.f;
			float Qq=(float) Qpro*327.68f*(1.f/coef);
			float Qold100=(float) Qpro/coef;
			
			float Qold=Qq;
			float redu=20.f;
			float reduc=1.f;		
			
			const Brightcurve& userColCurveB2 = static_cast<const Brightcurve&>(customColCurve2);
					userColCurveB2.Apply(Qq);
				if(Qq>Qold)	{
					if(Qq<65535.f)	{
						if(Qold < 327.68f*redu) Qq=0.25f*(Qq-Qold)+Qold;//divide sensibility
						else 			{
										reduc=LIM((100.f-Qold100)/(100.f-redu),0.f,1.f);
										Qq=0.25f*reduc*(Qq-Qold)+Qold;//reduct sensibility in highlights
										}
									}
							}
				else if(Qq>10.f) Qq=0.5f*(Qq-Qold)+Qold;
				else if (Qq>=0.f) Qq=0.7f*(Qq-Qold)+Qold;// not zero ==>artifacts
			Qpro=(float)(Qq*(coef)/327.68f);
			Jpro=100.f*(Qpro*Qpro)/((4.0f/c)*(4.0f/c)*(aw+4.0f)*(aw+4.0f));
			t2B=true;
			
			if(t1L){//to workaround the problem if we modify curve1-lightnees after curve2 brightness(the cat that bites its own tail!) in fact it's another type of curve only for this case
			coef=2.f;//adapt Q to J approximation
			Qq=(float) Qpro*coef;
			Qold=Qq;
			const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&>(customColCurve1);
				userColCurveJ1.Apply(Qq);
				Qq=0.05f*(Qq-Qold)+Qold;//approximative adaptation	
			Qpro=(float)(Qq/coef);
			Jpro=100.f*(Qpro*Qpro)/((4.0f/c)*(4.0f/c)*(aw+4.0f)*(aw+4.0f));			
			} 
			
			}
	}

	if (hasColCurve3) {//curve 3 with chroma saturation colorfullness
		if (curveMode3==ColorAppearanceParams::TC_MODE_CHROMA){
		    float parsat=0.8f;//0.68;
			float coef=327.68f/parsat;
			float Cc=(float) Cpro*coef;
			float Ccold=Cc;
			const Chromacurve& userColCurve = static_cast<const Chromacurve&>(customColCurve3);
				userColCurve.Apply(Cc);
				float dred=55.f;
				float protect_red=30.0f;
				float sk=1;
				float ko=1.f/coef;
				Color::skinredfloat(Jpro, hpro, Cc, Ccold, dred, protect_red,sk,rstprotection,ko, Cpro);
		}
	else if (curveMode3==ColorAppearanceParams::TC_MODE_SATUR){ //
				float parsat=0.8f;//0.6
				float coef=327.68f/parsat;
				float Ss=(float) spro*coef;
				float Sold=Ss;
				const Saturcurve& userColCurve = static_cast<const Saturcurve&>(customColCurve3);
					userColCurve.Apply(Ss);
					Ss=0.6f*(Ss-Sold)+Sold;//divide sensibility saturation
				float coe=pow_F(fl,0.25f);
				float dred=100.f;// in C mode
				float protect_red=80.0f; // in C mode
				dred = 100.0f * sqrt((dred*coe)/Qpro);
				protect_red=100.0f * sqrt((protect_red*coe)/Qpro);
				int sk=0;
				float ko=1.f/coef;
				Color::skinredfloat(Jpro, hpro, Ss, Sold, dred, protect_red,sk,rstprotection,ko, spro);
				Qpro= ( 4.0f / c ) * sqrt( Jpro / 100.0f ) * ( aw + 4.0f ) ;
				Cpro=(spro*spro*Qpro)/(10000.0f);
				c1s=1;

			}
	else if (curveMode3==ColorAppearanceParams::TC_MODE_COLORF){ //
				float parsat=0.8f;//0.68;
				float coef=327.68f/parsat;
				float Mm=(float) Mpro*coef;
				float Mold=Mm;
				const Colorfcurve& userColCurve = static_cast<const Colorfcurve&>(customColCurve3);
					userColCurve.Apply(Mm);
				float coe=pow_F(fl,0.25f);
				float dred=100.f;//in C mode
				float protect_red=80.0f;// in C mode
				dred *=coe;//in M mode
				protect_red	*=coe;
				int sk=0;
				float ko=1.f/coef;
				Color::skinredfloat(Jpro, hpro, Mm, Mold, dred, protect_red,sk,rstprotection,ko, Mpro);
				Cpro= Mpro/coe;
				c1co=1;
			}
	}
			//to retrieve the correct values of variables
			
			if(c1s==1) {
				Qpro= ( 4.0f / c ) * sqrt( Jpro / 100.0f ) * ( aw + 4.0f ) ;//for saturation curve
				Cpro=(spro*spro*Qpro)/(10000.0f);
				}
			if(c1co==1) {	float coe=pow_F(fl,0.25f);Cpro= Mpro/coe;}	// for colorfullness curve
			//retrieve values C,J...s
			C=Cpro;
			J=Jpro;
			Q=Qpro;
			M=Mpro;
			h=hpro;
			s=spro;

		if(params->colorappearance.tonecie  || settings->autocielab){//use pointer for tonemapping with CIECAM and also sharpening , defringe, contrast detail
			float Qred= ( 4.0f / c)  * ( aw + 4.0f );//estimate Q max if J=100.0
			ncie->Q_p[i][j]=(float)Q+epsil;//epsil to avoid Q=0
			ncie->M_p[i][j]=(float)M+epsil;
			ncie->J_p[i][j]=(float)J+epsil;
			ncie->h_p[i][j]=(float)h;
			ncie->C_p[i][j]=(float)C+epsil;
			ncie->sh_p[i][j]=(float) 32768.f*(( 4.0f / c )*sqrt( J / 100.0f ) * ( aw + 4.0f ))/Qred ;
			if(ncie->Q_p[i][j]<minQ) minQ=ncie->Q_p[i][j];//minima
			if(ncie->Q_p[i][j]>maxQ) maxQ=ncie->Q_p[i][j];//maxima
			}
			if(!params->colorappearance.tonecie  || !settings->autocielab || !params->edgePreservingDecompositionUI.enabled){
			int posl, posc;
			float brli=327.f;
			float chsacol=327.f;
			int libr=0;
			int colch=0;
			if(curveMode==ColorAppearanceParams::TC_MODE_BRIGHT) {brli=70.0f; libr=1;}
			else if(curveMode==ColorAppearanceParams::TC_MODE_LIGHT) {brli=327.f;libr=0;}
			if (curveMode3==ColorAppearanceParams::TC_MODE_CHROMA) {chsacol=327.f;colch=0;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_SATUR) {chsacol=450.0f;colch=1;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_COLORF) {chsacol=327.0f;colch=2;}
			
			if(ciedata) {
			// Data for J Q M s and C histograms
				//update histogram
				jp=true;
                if(pW!=1){//only with improccoordinator
				if(libr==1) posl=CLIP((int)(Q*brli));//40.0 to 100.0 approximative factor for Q  - 327 for J
				else if(libr==0) posl=CLIP((int)(J*brli));//327 for J
				hist16JCAM[posl]++;
				}
				chropC=true;
                if(pW!=1){//only with improccoordinator
				if(colch==0) posc=CLIP((int)(C*chsacol));//450.0 approximative factor for s    320 for M
				else if(colch==1) posc=CLIP((int)(s*chsacol));
				else if(colch==2) posc=CLIP((int)(M*chsacol));
				hist16_CCAM[posc]++;
				}
			}
			float xx,yy,zz;
			//process normal==> viewing

			ColorTemp::jch2xyz_ciecam02float( xx, yy, zz,
			                             J,  C, h,
			                             xw2, yw2,  zw2,
			                             yb2, la2,
			                             f2,  c2, nc2, gamu, nj, nbbj, ncbj, flj, czj, dj, awj);
			x=(float)xx*655.35f;
			y=(float)yy*655.35f;
			z=(float)zz*655.35f;
			float Ll,aa,bb;
			//convert xyz=>lab
			Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
			lab->L[i][j]=Ll;
			lab->a[i][j]=aa;
			lab->b[i][j]=bb;
		// gamut control in Lab mode; I must study how to do with cIECAM only
		if(gamu==1) {
					float R,G,B;
					float HH, Lprov1, Chprov1;
					Lprov1=lab->L[i][j]/327.68f;
					Chprov1=sqrt(SQR(lab->a[i][j]/327.68f) + SQR(lab->b[i][j]/327.68f));
					HH=xatan2f(lab->b[i][j],lab->a[i][j]);

#ifdef _DEBUG
					bool neg=false;
					bool more_rgb=false;
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif

					lab->L[i][j]=Lprov1*327.68f;
					float2 sincosval = xsincosf(HH);
					
					lab->a[i][j]=327.68f*Chprov1*sincosval.y;
					lab->b[i][j]=327.68f*Chprov1*sincosval.x;

		}
		}
		}
	}
	// End of parallelization
if(!params->colorappearance.tonecie   || !settings->autocielab){//normal

	if(ciedata) {
    //update histogram J
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=32768; i++) {//
			if (jp) {
				float hval = dLcurve[i];
				int hi = (int)(255.0f*CLIPD(hval)); //
				histLCAM[hi] += hist16JCAM[i] ;
			}
		}
	}
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=48000; i++) {//
			if (chropC) {
				float hvalc = dCcurve[i];
				int hic = (int)(255.0f*CLIPD(hvalc)); //
				histCCAM[hic] += hist16_CCAM[i] ;
			}
		}
	}
	}
}
#ifdef _DEBUG
	if (settings->verbose) {
		t2e.set();
		printf("CIECAM02 performed in %d usec:\n", t2e.etime(t1e));
		//	printf("minc=%f maxc=%f minj=%f maxj=%f\n",minc,maxc,minj,maxj);
	}
#endif

if(settings->autocielab) {



//all this treatments reduce artefacts, but can leed to slighty  different results

if(params->defringe.enabled) if(execsharp) ImProcFunctions::defringecam (ncie);//defringe adapted to CIECAM

if(params->colorappearance.badpixsl > 0) if(execsharp){ int mode=params->colorappearance.badpixsl;
											ImProcFunctions::badpixcam (ncie, 3.0, 10, mode);//for bad pixels
										}	

if(params->impulseDenoise.enabled) if(execsharp) ImProcFunctions::impulsedenoisecam (ncie);//impulse adapted to CIECAM

if (params->sharpenMicro.enabled)if(execsharp) ImProcFunctions::MLmicrocontrastcam(ncie);

if(params->sharpening.enabled)  if(execsharp) {ImProcFunctions::sharpeningcam (ncie, (float**)buffer);}			//sharpening adapted to CIECAM

if(params->dirpyrequalizer.enabled) if(execsharp) dirpyr_equalizercam(ncie, ncie->sh_p, ncie->sh_p, ncie->W, ncie->H, params->dirpyrequalizer.mult, params->dirpyrequalizer.threshold, true);//contrast by detail adapted to CIECAM

		   float Qredi= ( 4.0f / c_)  * ( a_w + 4.0f );
		   float co_e=(pow_F(f_l,0.25f));

#ifndef _DEBUG
#pragma omp parallel default(shared) firstprivate(height,width, Qredi,a_w,c_)
#endif
{
#ifndef _DEBUG
		#pragma omp for schedule(dynamic, 10)
#endif
		   		for (int i=0; i<height; i++) // update CieImages with new values after sharpening, defringe, contrast by detail level
					for (int j=0; j<width; j++) {
						float interm=Qredi*ncie->sh_p[i][j]/(32768.f);
						ncie->J_p[i][j]=100.0f* interm*interm/((a_w+4.f)*(a_w+4.f)*(4.f/c_)*(4.f/c_));
						ncie->Q_p[i][j]=( 4.0f / c_)  * ( a_w + 4.0f ) *  sqrt(ncie->J_p[i][j]/100.f);
						ncie->M_p[i][j]=ncie->C_p[i][j]*co_e;
					}
			}
}
if((params->colorappearance.tonecie && (params->edgePreservingDecompositionUI.enabled)) || (params->sharpening.enabled && settings->autocielab) 
		|| (params->dirpyrequalizer.enabled && settings->autocielab) ||(params->defringe.enabled && settings->autocielab)  || (params->sharpenMicro.enabled && settings->autocielab)
		|| (params->impulseDenoise.enabled && settings->autocielab) ||  (params->colorappearance.badpixsl >0 && settings->autocielab)){
		
		if(params->edgePreservingDecompositionUI.enabled  && params->colorappearance.tonecie && algepd) ImProcFunctions::EPDToneMapCIE(ncie, a_w, c_, w_h, width, height, begh, endh, minQ, maxQ, Iterates, scale );
			//EPDToneMapCIE adated to CIECAM

	
#ifndef _DEBUG	
#pragma omp parallel default(shared) firstprivate(lab,xw2,yw2,zw2,chr,yb,la2,yb2, height,width,begh, endh, nc2,f2,c2, gamu, highlight,pW,nj, nbbj, ncbj, flj, czj, dj, awj)
#endif
{	
	TMatrix wiprofa = iccStore->workingSpaceInverseMatrix (params->icm.working);
	double wipa[3][3] = {
		{wiprofa[0][0],wiprofa[0][1],wiprofa[0][2]},
		{wiprofa[1][0],wiprofa[1][1],wiprofa[1][2]},
		{wiprofa[2][0],wiprofa[2][1],wiprofa[2][2]}
	};
	
	
#ifndef _DEBUG
		#pragma omp for schedule(dynamic, 10)
#endif
		for (int i=0; i<height; i++) // update CIECAM with new values after tone-mapping
			for (int j=0; j<width; j++) {
			float xx,yy,zz;
			float x,y,z;
			float eps=0.0001f;
			float co_e=(pow_F(f_l,0.25f))+eps;
		//	if(params->edgePreservingDecompositionUI.enabled) ncie->J_p[i][j]=(100.0f* ncie->Q_p[i][j]*ncie->Q_p[i][j])/(w_h*w_h);
			if(params->edgePreservingDecompositionUI.enabled) ncie->J_p[i][j]=(100.0f* ncie->Q_p[i][j]*ncie->Q_p[i][j])/SQR((4.f/c)*(aw+4.f));
			
			ncie->C_p[i][j]	=(ncie->M_p[i][j])/co_e;
			//show histogram in CIECAM mode (Q,J, M,s,C)
			int posl, posc;
			float brli=327.f;
			float chsacol=327.f;
			int libr=0;
			int colch=0;
			float sa_t;
			if(curveMode==ColorAppearanceParams::TC_MODE_BRIGHT) {brli=70.0f; libr=1;}
			else if(curveMode==ColorAppearanceParams::TC_MODE_LIGHT) {brli=327.f;libr=0;}
			if (curveMode3==ColorAppearanceParams::TC_MODE_CHROMA) {chsacol=327.f;colch=0;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_SATUR) {chsacol=450.0f;colch=1;}
			else if(curveMode3==ColorAppearanceParams::TC_MODE_COLORF) {chsacol=327.0f;colch=2;}
			if(ciedata) {
			// Data for J Q M s and C histograms
				//update histogram
				jp=true;
                if(pW!=1){//only with improccoordinator
				if(libr==1) posl=CLIP((int)(ncie->Q_p[i][j]*brli));//40.0 to 100.0 approximative factor for Q  - 327 for J
				else if(libr==0) posl=CLIP((int)(ncie->J_p[i][j]*brli));//327 for J
				hist16JCAM[posl]++;
				}
				chropC=true;
                if(pW!=1){//only with improccoordinator
			if(colch==0) posc=CLIP((int)(ncie->C_p[i][j]*chsacol));//450.0 approximative factor for s    320 for M
				else if(colch==1) {sa_t=100.f*sqrt(ncie->C_p[i][j]/ncie->Q_p[i][j]); posc=CLIP((int)(sa_t*chsacol));}//Q_p always > 0
				else if(colch==2) posc=CLIP((int)(ncie->M_p[i][j]*chsacol));
				hist16_CCAM[posc]++;
				}
			}
			//end histograms

			ColorTemp::jch2xyz_ciecam02float( xx, yy, zz,
			                             ncie->J_p[i][j],  ncie->C_p[i][j], ncie->h_p[i][j],
			                             xw2, yw2,  zw2,
			                             yb2, la2,
			                             f2,  c2, nc2, gamu, nj, nbbj, ncbj, flj, czj, dj, awj);
			x=(float)xx*655.35f;
			y=(float)yy*655.35f;
			z=(float)zz*655.35f;
			float Ll,aa,bb;
			//convert xyz=>lab
			Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
			lab->L[i][j]=Ll;
			lab->a[i][j]=aa;
			lab->b[i][j]=bb;
			if(gamu==1) {
					float R,G,B;
					float HH, Lprov1, Chprov1;
					Lprov1=lab->L[i][j]/327.68f;
					Chprov1=sqrt(SQR(lab->a[i][j]/327.68f) + SQR(lab->b[i][j]/327.68f));
					HH=xatan2f(lab->b[i][j],lab->a[i][j]);

#ifdef _DEBUG
					bool neg=false;
					bool more_rgb=false;
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wipa, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wipa, highlight, 0.15f, 0.96f);
#endif
					float2 sincosval = xsincosf(HH);

					lab->L[i][j]=Lprov1*327.68f;
					lab->a[i][j]=327.68f*Chprov1*sincosval.y;
					lab->b[i][j]=327.68f*Chprov1*sincosval.x;
						}
										}


		}
		//end parallelization
	//show CIECAM histograms
	if(ciedata) {
    //update histogram J and Q
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=32768; i++) {//
			if (jp) {
				float hval = dLcurve[i];
				int hi = (int)(255.0f*CLIPD(hval)); //
				histLCAM[hi] += hist16JCAM[i] ;
			}
		}
	}
	//update color histogram M,s,C
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=48000; i++) {//
			if (chropC) {
				float hvalc = dCcurve[i];
				int hic = (int)(255.0f*CLIPD(hvalc)); //
				histCCAM[hic] += hist16_CCAM[i] ;
			}
		}
	}
	}

	}
}
}
//end CIECAM
void ImProcFunctions::rgbProc (Imagefloat* working, LabImage* lab, EditBuffer *editBuffer, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                               SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve,
                               const ToneCurve & customToneCurve1,const ToneCurve & customToneCurve2, const ToneCurve & customToneCurvebw1,const ToneCurve & customToneCurvebw2, double &rrm, double &ggm, double &bbm, float &autor, float &autog, float &autob ) {
    rgbProc (working, lab, editBuffer, hltonecurve, shtonecurve, tonecurve, shmap, sat, rCurve, gCurve, bCurve, customToneCurve1, customToneCurve2,  customToneCurvebw1, customToneCurvebw2,rrm, ggm, bbm, autor, autog, autob, params->toneCurve.expcomp, params->toneCurve.hlcompr, params->toneCurve.hlcomprthresh);
}

// Process RGB image and convert to LAB space
void ImProcFunctions::rgbProc (Imagefloat* working, LabImage* lab, EditBuffer *editBuffer, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                               SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve, const ToneCurve & customToneCurve1,
                               const ToneCurve & customToneCurve2,  const ToneCurve & customToneCurvebw1,const ToneCurve & customToneCurvebw2,double &rrm, double &ggm, double &bbm, float &autor, float &autog, float &autob, double expcomp, int hlcompr, int hlcomprthresh) {

    LUTf iGammaLUTf;
    Imagefloat *tmpImage;

    // NOTE: We're getting all 3 pointers here, but this function may not need them all, so one could optimize this
    Imagefloat* editImgFloat = NULL;
    LabImage* editLab = NULL;
    PlanarWhateverData<float>* editWhatever = NULL;
    EditUniqueID editID = editBuffer ? editBuffer->getEditID() : EUID_None;
    if (editID != EUID_None) {
        switch  (editBuffer->getDataProvider()->getCurrSubscriber()->getEditBufferType()) {
        case (BT_IMAGEFLOAT):
            editImgFloat = editBuffer->getImgFloatBuffer();
            break;
        case (BT_LABIMAGE):
            editLab = editBuffer->getLabBuffer();
            break;
        case (BT_SINGLEPLANE_FLOAT):
            editWhatever = editBuffer->getSinglePlaneBuffer();
            break;
        }
    }

    int h_th, s_th;
    if (shmap) {
        h_th = shmap->max_f - params->sh.htonalwidth * (shmap->max_f - shmap->avg) / 100;
        s_th = params->sh.stonalwidth * (shmap->avg - shmap->min_f) / 100;
    }

    bool processSH  = params->sh.enabled && shmap!=NULL && (params->sh.highlights>0 || params->sh.shadows>0);
    bool processLCE = params->sh.enabled && shmap!=NULL && params->sh.localcontrast>0;
    double lceamount = params->sh.localcontrast / 200.0;

    TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);
    TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);

    double toxyz[3][3] = {
        {
            ( wprof[0][0] / Color::D50x),
            ( wprof[0][1] / Color::D50x),
            ( wprof[0][2] / Color::D50x)
        },{
            ( wprof[1][0]),
            ( wprof[1][1]),
            ( wprof[1][2])
        },{
            ( wprof[2][0] / Color::D50z),
            ( wprof[2][1] / Color::D50z),
            ( wprof[2][2] / Color::D50z)
        }
    };

	//inverse matrix user select
	double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};

	double wp[3][3] = {
		{wprof[0][0],wprof[0][1],wprof[0][2]},
		{wprof[1][0],wprof[1][1],wprof[1][2]},
		{wprof[2][0],wprof[2][1],wprof[2][2]}};


    bool mixchannels = (params->chmixer.red[0]!=100	|| params->chmixer.red[1]!=0     || params->chmixer.red[2]!=0   ||
						params->chmixer.green[0]!=0 || params->chmixer.green[1]!=100 || params->chmixer.green[2]!=0 ||
						params->chmixer.blue[0]!=0	|| params->chmixer.blue[1]!=0    || params->chmixer.blue[2]!=100);

	double pi = M_PI;
	FlatCurve* hCurve;
	FlatCurve* sCurve;
	FlatCurve* vCurve;
	FlatCurve* bwlCurve;

	FlatCurveType hCurveType = (FlatCurveType)params->hsvequalizer.hcurve.at(0);
	FlatCurveType sCurveType = (FlatCurveType)params->hsvequalizer.scurve.at(0);
	FlatCurveType vCurveType = (FlatCurveType)params->hsvequalizer.vcurve.at(0);
	FlatCurveType bwlCurveType = (FlatCurveType)params->blackwhite.luminanceCurve.at(0);
	bool hCurveEnabled = hCurveType > FCT_Linear;
	bool sCurveEnabled = sCurveType > FCT_Linear;
	bool vCurveEnabled = vCurveType > FCT_Linear;
	bool bwlCurveEnabled = bwlCurveType > FCT_Linear;

	// TODO: We should create a 'skip' value like for CurveFactory::complexsgnCurve (rtengine/curves.cc)
	if (hCurveEnabled) hCurve = new FlatCurve(params->hsvequalizer.hcurve);
	if (sCurveEnabled) sCurve = new FlatCurve(params->hsvequalizer.scurve);
	if (vCurveEnabled) vCurve = new FlatCurve(params->hsvequalizer.vcurve);
	if (bwlCurveEnabled) {
		bwlCurve = new FlatCurve(params->blackwhite.luminanceCurve);
		if (bwlCurve->isIdentity()) {
			delete bwlCurve;
			bwlCurve = NULL;
			bwlCurveEnabled = false;
		}
	}

	const float exp_scale = pow (2.0, expcomp);
	const float comp = (max(0.0, expcomp) + 1.0)*hlcompr/100.0;
	const float shoulder = ((65536.0/max(1.0f,exp_scale))*(hlcomprthresh/200.0))+0.1;
	const float hlrange = 65536.0-shoulder;

    // extracting datas from 'params' to avoid cache flush (to be confirmed)
    ToneCurveParams::eTCModeId curveMode = params->toneCurve.curveMode;
    ToneCurveParams::eTCModeId curveMode2 = params->toneCurve.curveMode2;
	bool highlight = params->toneCurve.hrenabled;//Get the value if "highlight reconstruction" is activated
    bool hasToneCurve1 = bool(customToneCurve1);
    bool hasToneCurve2 = bool(customToneCurve2);
    BlackWhiteParams::eTCModeId beforeCurveMode = params->blackwhite.beforeCurveMode;
    BlackWhiteParams::eTCModeId afterCurveMode = params->blackwhite.afterCurveMode;

    bool hasToneCurvebw1 = bool(customToneCurvebw1);
    bool hasToneCurvebw2 = bool(customToneCurvebw2);

    float chMixRR = float(params->chmixer.red[0]);
    float chMixRG = float(params->chmixer.red[1]);
    float chMixRB = float(params->chmixer.red[2]);
    float chMixGR = float(params->chmixer.green[0]);
    float chMixGG = float(params->chmixer.green[1]);
    float chMixGB = float(params->chmixer.green[2]);
    float chMixBR = float(params->chmixer.blue[0]);
    float chMixBG = float(params->chmixer.blue[1]);
    float chMixBB = float(params->chmixer.blue[2]);
    int shHighlights = params->sh.highlights;
    int shShadows = params->sh.shadows;
    bool blackwhite = params->blackwhite.enabled;
    bool complem = params->blackwhite.enabledcc;
    float bwr = float(params->blackwhite.mixerRed);
    float bwg = float(params->blackwhite.mixerGreen);
    float bwb = float(params->blackwhite.mixerBlue);
    float bwrgam = float(params->blackwhite.gammaRed);
    float bwggam = float(params->blackwhite.gammaGreen);
    float bwbgam = float(params->blackwhite.gammaBlue);
    float mixerOrange = float(params->blackwhite.mixerOrange);
    float mixerYellow = float(params->blackwhite.mixerYellow);
    float mixerCyan = float(params->blackwhite.mixerCyan);
    float mixerMagenta = float(params->blackwhite.mixerMagenta);
    float mixerPurple = float(params->blackwhite.mixerPurple);
	int algm=0;
	if     (params->blackwhite.method=="Desaturation")  algm=0;
	else if(params->blackwhite.method=="LumEqualizer")  algm=1;
	else if(params->blackwhite.method=="ChannelMixer")  algm=2;
	float kcorec=1.f;
	//gamma correction of each channel
	float gamvalr=125.f;
	float gamvalg=125.f;
	float gamvalb=125.f;
	double nr=0;
	double ng=0;
	double nb=0;
	bool computeMixerAuto = params->blackwhite.autoc && (autor < -5000.f);
	if(bwrgam < 0) gamvalr=100.f;
	if(bwggam < 0) gamvalg=100.f;
	if(bwbgam < 0) gamvalb=100.f;
	float gammabwr=1.f;
	float gammabwg=1.f;
	float gammabwb=1.f;
	//if     (params->blackwhite.setting=="Ma" || params->blackwhite.setting=="Mr" || params->blackwhite.setting=="Fr" || params->blackwhite.setting=="Fa")  {
	{
		gammabwr=1.f -bwrgam/gamvalr;
		gammabwg=1.f -bwggam/gamvalg;
		gammabwb=1.f -bwbgam/gamvalb;
	}
	bool hasgammabw = gammabwr!=1.f || gammabwg!=1.f || gammabwb!=1.f;

	//normalize gamma to sRGB
	double start = exp(g*log( -0.055 / ((1.0/g-1.0)*1.055 )));
	double slope = 1.055 * pow (start, 1.0/g-1) - 0.055/start;
	double mul = 1.055;
	double add = 0.055;

	if (iGamma && g > 1.) {
		iGammaLUTf(65535);
#pragma omp parallel for
		for (int i=0; i<65536; i++) {
			iGammaLUTf[i] = float(CurveFactory::igamma (double(i)/65535., g, start, slope, mul, add)*65535.);
		}
	}

	if(blackwhite)
		tmpImage = new Imagefloat(working->width,working->height);

#define TS 112

#ifdef _OPENMP
#pragma omp parallel if (multiThread)
#endif
{
	char *buffer;
	char *editIFloatBuffer = NULL;
	char *editWhateverBuffer = NULL;

	buffer = (char *) malloc(3*sizeof(float)*TS*TS + 20*64 + 63);
	char *data;
	data = (char*)( ( uintptr_t(buffer) + uintptr_t(63)) / 64 * 64);

	float *rtemp = (float(*))data;
	float *gtemp = (float (*))         ((char*)rtemp + sizeof(float)*TS*TS + 4*64);
	float *btemp = (float (*))         ((char*)gtemp + sizeof(float)*TS*TS + 8*64);
	int istart;
	int jstart;
	int tW;
	int tH;

	// Allocating buffer for the EditBuffer
	float *editIFloatTmpR, *editIFloatTmpG, *editIFloatTmpB, *editWhateverTmp;
	char *editIFBuffer = NULL;
	if (editImgFloat) {
		editIFloatBuffer = (char *) malloc(3*sizeof(float)*TS*TS + 20*64 + 63);
		data = (char*)( ( uintptr_t(editIFloatBuffer) + uintptr_t(63)) / 64 * 64);

		editIFloatTmpR = (float(*))data;
		editIFloatTmpG = (float (*))         ((char*)editIFloatTmpR + sizeof(float)*TS*TS + 4*64);
		editIFloatTmpB = (float (*))         ((char*)editIFloatTmpG + sizeof(float)*TS*TS + 8*64);
	}
	if (editWhatever) {
		editWhateverBuffer = (char *) malloc(sizeof(float)*TS*TS + 20*64 + 63);
		data = (char*)( ( uintptr_t(editWhateverBuffer) + uintptr_t(63)) / 64 * 64);

		editWhateverTmp = (float(*))data;
	}

#pragma omp for schedule(dynamic) collapse(2) nowait
	for(int ii=0;ii<working->height;ii+=TS)
		for(int jj=0;jj<working->width;jj+=TS) {
			istart = ii;
			jstart = jj;
			tH = min(ii+TS,working->height);
			tW = min(jj+TS,working->width);


			for (int i=istart,ti=0; i<tH; i++,ti++) {
				for (int j=jstart,tj=0; j<tW; j++,tj++) {
					rtemp[ti*TS+tj] = working->r(i,j);
					gtemp[ti*TS+tj] = working->g(i,j);
					btemp[ti*TS+tj] = working->b(i,j);
				}
			}

			if (mixchannels) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						float r = rtemp[ti*TS+tj];
						float g = gtemp[ti*TS+tj];
						float b = btemp[ti*TS+tj];

						//if (i==100 & j==100) printf("rgbProc input R= %f  G= %f  B= %f  \n",r,g,b);
						float rmix = (r*chMixRR + g*chMixRG + b*chMixRB) / 100.f;
						float gmix = (r*chMixGR + g*chMixGG + b*chMixGB) / 100.f;
						float bmix = (r*chMixBR + g*chMixBG + b*chMixBB) / 100.f;

						rtemp[ti*TS+tj] = rmix;
						gtemp[ti*TS+tj] = gmix;
						btemp[ti*TS+tj] = bmix;
					}
				}
			}
	
			if (processSH || processLCE) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {

						float r = rtemp[ti*TS+tj];
						float g = gtemp[ti*TS+tj];
						float b = btemp[ti*TS+tj];

						double mapval = 1.0 + shmap->map[i][j];
						double factor = 1.0;

						if (processSH) {
							if (mapval > h_th)
								factor = (h_th + (100.0 - shHighlights) * (mapval - h_th) / 100.0) / mapval;
							else if (mapval < s_th)
								factor = (s_th - (100.0 - shShadows) * (s_th - mapval) / 100.0) / mapval;
						}
						if (processLCE) {
							double sub = lceamount*(mapval-factor*(r*lumimul[0] + g*lumimul[1] + b*lumimul[2]));
							rtemp[ti*TS+tj] = factor*r-sub;
							gtemp[ti*TS+tj] = factor*g-sub;
							btemp[ti*TS+tj] = factor*b-sub;
						}
						else {
							rtemp[ti*TS+tj] = factor*r;
							gtemp[ti*TS+tj] = factor*g;
							btemp[ti*TS+tj] = factor*b;
						}
					}
				}
			}

			for (int i=istart,ti=0; i<tH; i++,ti++) {
				for (int j=jstart,tj=0; j<tW; j++,tj++) {

					float r = rtemp[ti*TS+tj];
					float g = gtemp[ti*TS+tj];
					float b = btemp[ti*TS+tj];

				//TODO: proper treatment of out-of-gamut colors
					//float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
					float tonefactor=((r<MAXVALF ? hltonecurve[r] : CurveFactory::hlcurve (exp_scale, comp, hlrange, r) ) +
									  (g<MAXVALF ? hltonecurve[g] : CurveFactory::hlcurve (exp_scale, comp, hlrange, g) ) +
									  (b<MAXVALF ? hltonecurve[b] : CurveFactory::hlcurve (exp_scale, comp, hlrange, b) ) )/3.0;

					rtemp[ti*TS+tj] = r*tonefactor;
					gtemp[ti*TS+tj] = g*tonefactor;
					btemp[ti*TS+tj] = b*tonefactor;
				}
			}

			for (int i=istart,ti=0; i<tH; i++,ti++) {
				for (int j=jstart,tj=0; j<tW; j++,tj++) {

					float r = rtemp[ti*TS+tj];
					float g = gtemp[ti*TS+tj];
					float b = btemp[ti*TS+tj];

					//shadow tone curve
					float Y = (0.299f*r + 0.587f*g + 0.114f*b);
					float tonefactor = shtonecurve[Y];
					rtemp[ti*TS+tj] = rtemp[ti*TS+tj]*tonefactor;
					gtemp[ti*TS+tj] = gtemp[ti*TS+tj]*tonefactor;
					btemp[ti*TS+tj] = btemp[ti*TS+tj]*tonefactor;
				}
			}

			for (int i=istart,ti=0; i<tH; i++,ti++) {
				for (int j=jstart,tj=0; j<tW; j++,tj++) {

					//brightness/contrast
					rtemp[ti*TS+tj] = tonecurve[ rtemp[ti*TS+tj] ];
					gtemp[ti*TS+tj] = tonecurve[ gtemp[ti*TS+tj] ];
					btemp[ti*TS+tj] = tonecurve[ btemp[ti*TS+tj] ];
				}
			}

			if (editID == EUID_ToneCurve1) {  // filling the pipette buffer
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						editIFloatTmpR[ti*TS+tj] = CLIP(rtemp[ti*TS+tj]/65535.f);
						editIFloatTmpG[ti*TS+tj] = CLIP(gtemp[ti*TS+tj]/65535.f);
						editIFloatTmpB[ti*TS+tj] = CLIP(btemp[ti*TS+tj]/65535.f);
					}
				}
			}

			if (hasToneCurve1) {
				if (curveMode==ToneCurveParams::TC_MODE_STD){ // Standard
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							const StandardToneCurve& userToneCurve = static_cast<const StandardToneCurve&>(customToneCurve1);
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
				else if (curveMode==ToneCurveParams::TC_MODE_FILMLIKE){ // Adobe like
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							const AdobeToneCurve& userToneCurve = static_cast<const AdobeToneCurve&>(customToneCurve1);
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
				else if (curveMode==ToneCurveParams::TC_MODE_SATANDVALBLENDING){ // apply the curve on the saturation and value channels
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							const SatAndValueBlendingToneCurve& userToneCurve = static_cast<const SatAndValueBlendingToneCurve&>(customToneCurve1);
							rtemp[ti*TS+tj] = CLIP<float>(rtemp[ti*TS+tj]);
							gtemp[ti*TS+tj] = CLIP<float>(gtemp[ti*TS+tj]);
							btemp[ti*TS+tj] = CLIP<float>(btemp[ti*TS+tj]);
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
				else if (curveMode==ToneCurveParams::TC_MODE_WEIGHTEDSTD){ // apply the curve to the rgb channels, weighted
					const WeightedStdToneCurve& userToneCurve = static_cast<const WeightedStdToneCurve&>(customToneCurve1);
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							rtemp[ti*TS+tj] = CLIP<float>(rtemp[ti*TS+tj]);
							gtemp[ti*TS+tj] = CLIP<float>(gtemp[ti*TS+tj]);
							btemp[ti*TS+tj] = CLIP<float>(btemp[ti*TS+tj]);
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
			}

			if (editID == EUID_ToneCurve2) {  // filling the pipette buffer
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						editIFloatTmpR[ti*TS+tj] = CLIP(rtemp[ti*TS+tj]/65535.f);
						editIFloatTmpG[ti*TS+tj] = CLIP(gtemp[ti*TS+tj]/65535.f);
						editIFloatTmpB[ti*TS+tj] = CLIP(btemp[ti*TS+tj]/65535.f);
					}
				}
			}

			if (hasToneCurve2) {
				if (curveMode2==ToneCurveParams::TC_MODE_STD){ // Standard
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							const StandardToneCurve& userToneCurve = static_cast<const StandardToneCurve&>(customToneCurve2);
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
				else if (curveMode2==ToneCurveParams::TC_MODE_FILMLIKE){ // Adobe like
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							const AdobeToneCurve& userToneCurve = static_cast<const AdobeToneCurve&>(customToneCurve2);
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
				else if (curveMode2==ToneCurveParams::TC_MODE_SATANDVALBLENDING){ // apply the curve on the saturation and value channels
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							const SatAndValueBlendingToneCurve& userToneCurve = static_cast<const SatAndValueBlendingToneCurve&>(customToneCurve2);
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
				else if (curveMode2==ToneCurveParams::TC_MODE_WEIGHTEDSTD){ // apply the curve to the rgb channels, weighted
					const WeightedStdToneCurve& userToneCurve = static_cast<const WeightedStdToneCurve&>(customToneCurve2);
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							userToneCurve.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
						}
					}
				}
			}

			if (iGammaLUTf) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						// apply inverse gamma
						rtemp[ti*TS+tj] = iGammaLUTf[ rtemp[ti*TS+tj] ];
						gtemp[ti*TS+tj] = iGammaLUTf[ gtemp[ti*TS+tj] ];
						btemp[ti*TS+tj] = iGammaLUTf[ btemp[ti*TS+tj] ];
					}
				}
			}

			if (editID == EUID_RGB_R) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						editWhateverTmp[ti*TS+tj] = rtemp[ti*TS+tj]/65536.f;
					}
				}
			}
			else if (editID == EUID_RGB_G) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						editWhateverTmp[ti*TS+tj] = gtemp[ti*TS+tj]/65536.f;
					}
				}
			}
			else if (editID == EUID_RGB_B) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						editWhateverTmp[ti*TS+tj] = btemp[ti*TS+tj]/65536.f;
					}
				}
			}

			if (rCurve || gCurve || bCurve) { // if any of the RGB curves is engaged
				if (!params->rgbCurves.lumamode){ // normal RGB mode
				
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
								// individual R tone curve
								if (rCurve) rtemp[ti*TS+tj] = rCurve[ rtemp[ti*TS+tj] ];
								// individual G tone curve
								if (gCurve) gtemp[ti*TS+tj] = gCurve[ gtemp[ti*TS+tj] ];
								// individual B tone curve
								if (bCurve) btemp[ti*TS+tj] = bCurve[ btemp[ti*TS+tj] ];
						}
					}
				} else { //params->rgbCurves.lumamode==true (Luminosity mode)
				// rCurve.dump("r_curve");//debug

					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							float r1,g1,b1, r2,g2,b2, L_1,L_2, Lfactor,a_1,b_1,x_,y_,z_,R,G,B ;
							float y,fy, yy,fyy,x,z,fx,fz;
							
							// rgb values before RGB curves
							r1 = rtemp[ti*TS+tj] ;
							g1 = gtemp[ti*TS+tj] ;
							b1 = btemp[ti*TS+tj] ;
							//convert to Lab to get a&b before RGB curves
							x = toxyz[0][0] * r1 + toxyz[0][1] * g1 + toxyz[0][2] * b1;
							y = toxyz[1][0] * r1 + toxyz[1][1] * g1 + toxyz[1][2] * b1;
							z = toxyz[2][0] * r1 + toxyz[2][1] * g1 + toxyz[2][2] * b1;
					
							fx = (x<65535.0f ? cachef[std::max(x,0.f)] : (327.68f*float(exp(log(x/MAXVALF)/3.0f ))));
							fy = (y<65535.0f ? cachef[std::max(y,0.f)] : (327.68f*float(exp(log(y/MAXVALF)/3.0f ))));
							fz = (z<65535.0f ? cachef[std::max(z,0.f)] : (327.68f*float(exp(log(z/MAXVALF)/3.0f ))));

							L_1 = (116.0f *  fy - 5242.88f); //5242.88=16.0*327.68;
							a_1 = (500.0f * (fx - fy) );
							b_1 = (200.0f * (fy - fz) );
							
							// rgb values after RGB curves
							if (rCurve) r2 = rCurve[ rtemp[ti*TS+tj]]; else r2=r1;
							if (gCurve) g2 = gCurve[ gtemp[ti*TS+tj]]; else g2=g1;
							if (bCurve) b2 = bCurve[ btemp[ti*TS+tj]]; else b2=b1;
							
							// Luminosity after
							// only Luminance in Lab
							yy = toxyz[1][0] * r2 + toxyz[1][1] * g2 + toxyz[1][2] * b2;
							fyy = (yy<65535.0f ? cachef[std::max(yy,0.f)] : (327.68f*float(exp(log(yy/MAXVALF)/3.0f ))));
							L_2 = (116.0f *  fyy - 5242.88f);

							//gamut control
							if(settings->rgbcurveslumamode_gamut) {
								float RR,GG,BB;
								float HH, Lpro, Chpro;
								Lpro=L_2/327.68f;					
								Chpro=sqrt(SQR(a_1) + SQR(b_1))/327.68f;
								HH=xatan2f(b_1,a_1);
								// According to mathematical laws we can get the sin and cos of HH by simple operations
								float2  sincosval;
								if(Chpro==0.0f) {
									sincosval.y = 1.0f;
									sincosval.x = 0.0f;
								} else {
									sincosval.y = a_1/(Chpro*327.68f);
									sincosval.x = b_1/(Chpro*327.68f);
								}

#ifdef _DEBUG
								bool neg=false;
								bool more_rgb=false;
								//gamut control : Lab values are in gamut
								Color::gamutLchonly(HH,Lpro,Chpro, RR, GG, BB, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
								//gamut control : Lab values are in gamut
								Color::gamutLchonly(HH,Lpro,Chpro, RR, GG, BB, wip, highlight, 0.15f, 0.96f);
#endif
								
								//Color::gamutLchonly(HH,Lpro,Chpro, RR, GG, BB, wip, highlight, 0.15f, 0.96f);



//								float2 sincosval = xsincosf(HH);

								L_2=Lpro*327.68f;
								a_1=327.68f*Chpro*sincosval.y;
								b_1=327.68f*Chpro*sincosval.x;
							} //end of gamut control
							
							//calculate RGB with L_2 and old value of a and b
							Color::Lab2XYZ(L_2, a_1, b_1, x_, y_, z_) ;
							Color::xyz2rgb(x_,y_,z_,R,G,B,wip);
							
							rtemp[ti*TS+tj] =R;
							gtemp[ti*TS+tj] =G;
							btemp[ti*TS+tj] =B;
						}
					}
				}
			}

			if (editID == EUID_HSV_H || editID == EUID_HSV_S || editID == EUID_HSV_V) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						float h,s,v;
						Color::rgb2hsv(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj], h, s, v);
						editWhateverTmp[ti*TS+tj] = h;
					}
				}
			}

			if (sat!=0 || hCurveEnabled || sCurveEnabled || vCurveEnabled) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {

						const float satby100 = sat/100.f;
						float r = rtemp[ti*TS+tj];
						float g = gtemp[ti*TS+tj];
						float b = btemp[ti*TS+tj];
						float h,s,v;
						Color::rgb2hsv(r,g,b,h,s,v);
						if (sat > 0) {
							s = (1.f-satby100)*s+satby100*(1.f-SQR(SQR(1.f-min(s,1.0f))));
							if (s<0.f) s=0.f;
						} else /*if (sat < 0)*/
							s *= 1.f+satby100;

						//HSV equalizer
						if (hCurveEnabled) {
							h = (hCurve->getVal(double(h)) - 0.5) * 2.f + h;
							if (h > 1.0f)
								h -= 1.0f;
							else if (h < 0.0f)
								h += 1.0f;
						}
						if (sCurveEnabled) {
							//shift saturation
							float satparam = (sCurve->getVal(double(h))-0.5) * 2;
							if (satparam > 0.00001f) {
								s = (1.f-satparam)*s+satparam*(1.f-SQR(1.f-min(s,1.0f)));
								if (s<0.f) s=0.f;
							} else if (satparam < -0.00001f)
								s *= 1.f+satparam;

						}
						if (vCurveEnabled) {
							if (v<0) v=0;  // important

							//shift value
							float valparam = vCurve->getVal((double)h)-0.5f;
							valparam *= (1.f-SQR(SQR(1.f-min(s,1.0f))));
							if (valparam > 0.00001f) {
								v = (1.f-valparam)*v+ valparam*(1.f-SQR(1.f-min(v,1.0f)));// SQR (SQR  to increase action and avoid artefacts
								if (v<0) v=0;
							} else {
								if (valparam < -0.00001f)
									v *= (1.f+ valparam);//1.99 to increase action
							}

						}
						Color::hsv2rgb(h, s, v, rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);	
					}
				}
			}

			// filling the pipette buffer
			if (editID == EUID_BlackWhiteBeforeCurve) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						editIFloatTmpR[ti*TS+tj] = CLIP(rtemp[ti*TS+tj]/65535.f);
						editIFloatTmpG[ti*TS+tj] = CLIP(gtemp[ti*TS+tj]/65535.f);
						editIFloatTmpB[ti*TS+tj] = CLIP(btemp[ti*TS+tj]/65535.f);
					}
				}
			}
			else if (editID==EUID_BlackWhiteLuminance) {
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						float X,Y,Z,L,aa,bb;
						//rgb=>lab
						Color::rgbxyz(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj], X, Y, Z, wp);
						//convert Lab
						Color::XYZ2Lab(X, Y, Z, L, aa, bb);
						//end rgb=>lab
						float HH=xatan2f(bb,aa);// HH hue in -3.141  +3.141

						editWhateverTmp[ti*TS+tj] = float(Color::huelab_to_huehsv2(HH));
					}
				}
			}

			//black and white
			if(blackwhite){
				if (hasToneCurvebw1) {
					if (beforeCurveMode==BlackWhiteParams::TC_MODE_STD_BW){ // Standard
						for (int i=istart,ti=0; i<tH; i++,ti++) {
							for (int j=jstart,tj=0; j<tW; j++,tj++) {
								const StandardToneCurvebw& userToneCurvebw = static_cast<const StandardToneCurvebw&>(customToneCurvebw1);
								userToneCurvebw.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
							}
						}
					}
					else if (beforeCurveMode==BlackWhiteParams::TC_MODE_FILMLIKE_BW){ // Adobe like
						for (int i=istart,ti=0; i<tH; i++,ti++) {
							for (int j=jstart,tj=0; j<tW; j++,tj++) {
								const AdobeToneCurvebw& userToneCurvebw = static_cast<const AdobeToneCurvebw&>(customToneCurvebw1);
								userToneCurvebw.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
							}
						}
					}
					else if (beforeCurveMode==BlackWhiteParams::TC_MODE_SATANDVALBLENDING_BW){ // apply the curve on the saturation and value channels
						for (int i=istart,ti=0; i<tH; i++,ti++) {
							for (int j=jstart,tj=0; j<tW; j++,tj++) {
								const SatAndValueBlendingToneCurvebw& userToneCurvebw = static_cast<const SatAndValueBlendingToneCurvebw&>(customToneCurvebw1);
								rtemp[ti*TS+tj] = CLIP<float>(rtemp[ti*TS+tj]);
								gtemp[ti*TS+tj] = CLIP<float>(gtemp[ti*TS+tj]);
								btemp[ti*TS+tj] = CLIP<float>(btemp[ti*TS+tj]);
								userToneCurvebw.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
							}
						}
					}
					else if (beforeCurveMode==BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW){ // apply the curve to the rgb channels, weighted
						for (int i=istart,ti=0; i<tH; i++,ti++) {
							for (int j=jstart,tj=0; j<tW; j++,tj++) {
								const WeightedStdToneCurvebw& userToneCurvebw = static_cast<const WeightedStdToneCurvebw&>(customToneCurvebw1);
								rtemp[ti*TS+tj] = CLIP<float>(rtemp[ti*TS+tj]);
								gtemp[ti*TS+tj] = CLIP<float>(gtemp[ti*TS+tj]);
								btemp[ti*TS+tj] = CLIP<float>(btemp[ti*TS+tj]);

								userToneCurvebw.Apply(rtemp[ti*TS+tj], gtemp[ti*TS+tj], btemp[ti*TS+tj]);
							}
						}
					}
				}

				if (algm==0){//lightness
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {

							float r = rtemp[ti*TS+tj];
							float g = gtemp[ti*TS+tj];
							float b = btemp[ti*TS+tj];

							// --------------------------------------------------

							// Method 1: Luminosity (code taken from Gimp)
							/*
							float maxi = max(r, g, b);
							float mini = min(r, g, b);
							r = g = b = (maxi+mini)/2;
							*/

							// Method 2: Luminance (former RT code)
							r = g = b = (0.299f*r + 0.587f*g + 0.114f*b);

							// --------------------------------------------------

							//gamma correction: pseudo TRC curve
							if (hasgammabw) Color::trcGammaBW (r, g, b, gammabwr, gammabwg, gammabwb);

							rtemp[ti*TS+tj] = r;
							gtemp[ti*TS+tj] = g;
							btemp[ti*TS+tj] = b;
						}
					}
				}
				else if (algm==1) {//Luminance mixer in Lab mode to avoid artifacts
					for (int i=istart,ti=0; i<tH; i++,ti++) {
						for (int j=jstart,tj=0; j<tW; j++,tj++) {
							//rgb=>lab
							float r = rtemp[ti*TS+tj];
							float g = gtemp[ti*TS+tj];
							float b = btemp[ti*TS+tj];
							float X,Y,Z;
							float L,aa,bb;
							Color::rgbxyz(r,g,b,X,Y,Z,wp);
							//convert Lab
							Color::XYZ2Lab(X, Y, Z, L, aa, bb);
							//end rgb=>lab
							//lab ==> Ch
							float CC=sqrt(SQR(aa/327.68f) + SQR(bb/327.68f));//CC chromaticity in 0..180 or more
							float HH=xatan2f(bb,aa);// HH hue in -3.141  +3.141
							float l_r;//Luminance Lab in 0..1
							l_r = L/32768.f;

							if (bwlCurveEnabled) {
								double hr = Color::huelab_to_huehsv2(HH);
								float valparam = float((bwlCurve->getVal(hr)-0.5f) * 2.0f);//get l_r=f(H)
								float kcc=(CC/70.f);//take Chroma into account...70 "middle" of chromaticity (arbitrary and simple), one can imagine other algorithme
								//reduct action for low chroma and increase action for high chroma
								valparam *= kcc; 
								if(valparam > 0.f) { l_r = (1.f-valparam)*l_r+ valparam*(1.f-SQR(SQR(SQR(SQR(1.f-min(l_r,1.0f))))));}// SQR (SQR((SQR)  to increase action in low light	
								else l_r *= (1.f + valparam);//for negative
								}
							L=l_r*32768.f;
							float RR,GG,BB;
							float Lr;
							Lr=L/327.68f;//for gamutlch
#ifdef _DEBUG
								bool neg=false;
								bool more_rgb=false;
								//gamut control : Lab values are in gamut
								Color::gamutLchonly(HH,Lr,CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
								//gamut control : Lab values are in gamut
								Color::gamutLchonly(HH,Lr,CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f);
#endif			
							//convert CH ==> ab
							L=Lr*327.68f;
							float a_,b_;
							a_=0.f;//grey
							b_=0.f;//grey
							//convert lab=>rgb
							Color::Lab2XYZ(L, a_, b_, X, Y, Z);
							float rr_,gg_,bb_;
							Color::xyz2rgb(X,Y,Z,rr_,gg_,bb_,wip);

							//gamma correction: pseudo TRC curve
							if (hasgammabw)
								Color::trcGammaBW (rr_, gg_, bb_, gammabwr, gammabwg, gammabwb);
							rtemp[ti*TS+tj] = rr_;
							gtemp[ti*TS+tj] = gg_;
							btemp[ti*TS+tj] = bb_;
						}
					}
				}
			}
			if(!blackwhite){
				// ready, fill lab
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {

						// filling the pipette buffer by the content of the temp pipette buffers
						if (editImgFloat) {
							editImgFloat->r(i,j) = editIFloatTmpR[ti*TS+tj];
							editImgFloat->g(i,j) = editIFloatTmpG[ti*TS+tj];
							editImgFloat->b(i,j) = editIFloatTmpB[ti*TS+tj];
						}
						else if (editWhatever) {
							editWhatever->v(i,j) = editWhateverTmp[ti*TS+tj];
						}

						float r = rtemp[ti*TS+tj];
						float g = gtemp[ti*TS+tj];
						float b = btemp[ti*TS+tj];

						float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
						float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
						float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;
						
						float fx,fy,fz;
						
						fx = (x<65535.0f ? cachef[std::max(x,0.f)] : (327.68f*float(exp(log(x/MAXVALF)/3.0f ))));
						fy = (y<65535.0f ? cachef[std::max(y,0.f)] : (327.68f*float(exp(log(y/MAXVALF)/3.0f ))));
						fz = (z<65535.0f ? cachef[std::max(z,0.f)] : (327.68f*float(exp(log(z/MAXVALF)/3.0f ))));

						lab->L[i][j] = (116.0f *  fy - 5242.88f); //5242.88=16.0*327.68;
						lab->a[i][j] = (500.0f * (fx - fy) );
						lab->b[i][j] = (200.0f * (fy - fz) );

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
			} else { // black & white
				// Auto channel mixer needs whole image, so we now copy to tmpImage and close the tiled processing
				for (int i=istart,ti=0; i<tH; i++,ti++) {
					for (int j=jstart,tj=0; j<tW; j++,tj++) {
						// filling the pipette buffer by the content of the temp pipette buffers
						if (editImgFloat) {
							editImgFloat->r(i,j) = editIFloatTmpR[ti*TS+tj];
							editImgFloat->g(i,j) = editIFloatTmpG[ti*TS+tj];
							editImgFloat->b(i,j) = editIFloatTmpB[ti*TS+tj];
						}
						else if (editWhatever) {
							editWhatever->v(i,j) = editWhateverTmp[ti*TS+tj];
						}

						tmpImage->r(i,j) = rtemp[ti*TS+tj];
						tmpImage->g(i,j) = gtemp[ti*TS+tj];
						tmpImage->b(i,j) = btemp[ti*TS+tj];
					}
				}
			}
		}
	free(buffer);
	if (editIFloatBuffer) free (editIFloatBuffer);
	if (editWhateverBuffer) free (editWhateverBuffer);
}

	//black and white
	if(blackwhite){
		int tW = working->width;
		int tH = working->height;

/*
	else if (algm==1) {//Luminance mixer
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 5)
#endif
		for (int i=0; i<tH; i++) {
			for (int j=0; j<tW; j++) {
				float r = tmpImage->r(i,j);
				float g = tmpImage->g(i,j);
				float b = tmpImage->b(i,j);
				if (bwlCurveEnabled) {
					float h,s,v;
					Color::rgb2hsv(r,g,b,h,s,v);
						if (v<0) v=0;
						float valparam = bwlCurve->getVal((double)h)-0.5f;
						valparam *= (1.f-(SQR(SQR(1.f-min(s,1.0f)))));
						//float valcor=1.f/(0.501f-valparam);

						if (valparam < -0.00001f || valparam > 0.00001f)
							v *= (1.f + 4.f*valparam);
						v = CLIPD(v);
					Color::hsv2rgb(h, s, v, r, g, b);
				}

				// get luminance
				r = g = b = (0.2126f*r + 0.7152f*g + 0.0722f*b); // (constant taken from Gimp, see https://git.gnome.org/browse/gimp/tree/libgimpcolor/gimprgb.h)
				//r = g = b = (0.299f*r + 0.587f*g + 0.114f*b);  // (obsolete constant)

				//gamma correction: pseudo TRC curve
				if (hasgammabw) Color::trcGammaBW (r, g, b, gammabwr, gammabwg, gammabwb);

				tmpImage->r(i,j) = r;
				tmpImage->g(i,j) = g;
				tmpImage->b(i,j) = b;
			}
		}
	}
*/

		if (algm==2) {//channel-mixer
			//end auto chmix
			float mix[3][3];

			if (computeMixerAuto) {
				// auto channel-mixer

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 5) reduction(+:nr,ng,nb)
#endif
				for (int i=0; i<tH; i++) {
					for (int j=0; j<tW; j++) {
						nr += tmpImage->r(i,j);
						ng += tmpImage->g(i,j);
						nb += tmpImage->b(i,j);
					}
				}

				double srgb = nr+ng+nb;
				double knr = srgb/nr;
				double kng = srgb/ng;
				double knb = srgb/nb;
				double sk = knr+kng+knb;
				autor=(float)(100.0*knr/sk);
				autog=(float)(100.0*kng/sk);
				autob=(float)(100.0*knb/sk);

			}

			if (params->blackwhite.autoc) {
				// auto channel-mixer
				bwr = autor;
				bwg = autog;
				bwb = autob;
				mixerOrange  = 33.f;
				mixerYellow  = 33.f;
				mixerMagenta = 33.f;
				mixerPurple  = 33.f;
				mixerCyan    = 33.f;
			}

			Color::computeBWMixerConstants(params->blackwhite.setting, params->blackwhite.filter,
					bwr, bwg, bwb, mixerOrange, mixerYellow, mixerCyan, mixerPurple, mixerMagenta,
					params->blackwhite.autoc, complem, kcorec, rrm, ggm, bbm);

			mix[0][0] = bwr;
			mix[1][0] = bwr;
			mix[2][0] = bwr;
			mix[0][1] = bwg;
			mix[1][1] = bwg;
			mix[2][1] = bwg;
			mix[0][2] = bwb;
			mix[1][2] = bwb;
			mix[2][2] = bwb;

			float in[3], val[3];

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 5)
#endif
			for (int i=0; i<tH; i++) {
				for (int j=0; j<tW; j++) {
					in[0] = tmpImage->r(i,j);
					in[1] = tmpImage->g(i,j);
					in[2] = tmpImage->b(i,j);
					//mix channel
					for (int end=0; end < 3 ; end++){
						val[end]=0.f;
						for (int beg=0; beg < 3 ; beg++) {
							val[end] += mix[end][beg] *in[beg];
						}
					}
					tmpImage->r(i,j) = tmpImage->g(i,j) = tmpImage->b(i,j) = CLIP(val[0]*kcorec);

					//gamma correction: pseudo TRC curve
					if (hasgammabw) Color::trcGammaBW (tmpImage->r(i,j), tmpImage->g(i,j), tmpImage->b(i,j), gammabwr, gammabwg, gammabwb);
				}
			} 
		}

		if (editID == EUID_BlackWhiteAfterCurve) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 5)
#endif
			for (int i=0; i<tH; i++) {
				for (int j=0; j<tW; j++) {
					editWhatever->v(i,j) = CLIP(tmpImage->r(i,j)/65535.f);  // assuming that r=g=b
				}
			}
		}

		if (hasToneCurvebw2) {

			if (afterCurveMode==BlackWhiteParams::TC_MODE_STD_BW){ // Standard
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 5)
#endif
				for (int i=0; i<tH; i++) {
					for (int j=0; j<tW; j++) {
						const StandardToneCurvebw& userToneCurve = static_cast<const StandardToneCurvebw&>(customToneCurvebw2);
						userToneCurve.Apply(tmpImage->r(i,j), tmpImage->g(i,j), tmpImage->b(i,j));
					}
				}
			}
			else if (afterCurveMode==BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW){ // apply the curve to the rgb channels, weighted
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 5)
#endif
				for (int i=0; i<tH; i++) {//for ulterior usage if bw data modified
					for (int j=0; j<tW; j++) {
						const WeightedStdToneCurvebw& userToneCurve = static_cast<const WeightedStdToneCurvebw&>(customToneCurvebw2);
						
						tmpImage->r(i,j) = CLIP<float>(tmpImage->r(i,j));
						tmpImage->g(i,j) = CLIP<float>(tmpImage->g(i,j));
						tmpImage->b(i,j) = CLIP<float>(tmpImage->b(i,j));

						userToneCurve.Apply(tmpImage->r(i,j), tmpImage->g(i,j), tmpImage->b(i,j));
					}
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 5)
#endif
		for (int i=0; i<tH; i++) {
			for (int j=0; j<tW; j++) {
				float r = tmpImage->r(i,j);
				float g = tmpImage->g(i,j);
				float b = tmpImage->b(i,j);

				float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
				float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
				float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;
				
				float fx,fy,fz;
				
				fx = (x<65535.0f ? cachef[std::max(x,0.f)] : (327.68f*float(exp(log(x/MAXVALF)/3.0f ))));
				fy = (y<65535.0f ? cachef[std::max(y,0.f)] : (327.68f*float(exp(log(y/MAXVALF)/3.0f ))));
				fz = (z<65535.0f ? cachef[std::max(z,0.f)] : (327.68f*float(exp(log(z/MAXVALF)/3.0f ))));

				lab->L[i][j] = (116.0f *  fy - 5242.88f); //5242.88=16.0*327.68;
				lab->a[i][j] = (500.0f * (fx - fy) );
				lab->b[i][j] = (200.0f * (fy - fz) );



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
		if(tmpImage)
			delete tmpImage;
	}


//	delete tmpImage;

	if (hCurveEnabled) delete hCurve;
	if (sCurveEnabled) delete sCurve;
	if (vCurveEnabled) delete vCurve;
	
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



void ImProcFunctions::chromiLuminanceCurve (EditBuffer *editBuffer, int pW, LabImage* lold, LabImage* lnew, LUTf & acurve, LUTf & bcurve, LUTf & satcurve,LUTf & lhskcurve, LUTf & clcurve, LUTf & curve, bool utili, bool autili, bool butili, bool ccutili, bool cclutili, bool clcutili, LUTu &histCCurve, LUTu &histCLurve, LUTu &histLLCurve, LUTu &histLCurve) {
	int W = lold->W;
	int H = lold->H;
   // lhskcurve.dump("lh_curve");
	//init Flatcurve for C=f(H)

	// NOTE: We're getting all 3 pointers here, but this function may not need them all, so one could optimize this
	Imagefloat* editImgFloat = NULL;
	LabImage* editLab = NULL;
	PlanarWhateverData<float>* editWhatever = NULL;
	EditUniqueID editID = EUID_None;
	if (editBuffer) {
		editID = editBuffer->getEditID();
		if (editID != EUID_None) {
			switch  (editBuffer->getDataProvider()->getCurrSubscriber()->getEditBufferType()) {
			case (BT_IMAGEFLOAT):
				editImgFloat = editBuffer->getImgFloatBuffer();
				break;
			case (BT_LABIMAGE):
				editLab = editBuffer->getLabBuffer();
				break;
			case (BT_SINGLEPLANE_FLOAT):
				editWhatever = editBuffer->getSinglePlaneBuffer();
				break;
			}
		}
	}

	FlatCurve* chCurve = NULL;// curve C=f(H)
	bool chutili = false;
	if (params->labCurve.chromaticity > -100) {
		chCurve = new FlatCurve(params->labCurve.chcurve);
		if (!chCurve || chCurve->isIdentity()) {
			if (chCurve) {
				delete chCurve;
				chCurve = NULL;
			}
		}//do not use "Munsell" if Chcurve not used
		else
			chutili=true;
	}
	FlatCurve* lhCurve = NULL;//curve L=f(H)
	bool lhutili = false;
	if (params->labCurve.chromaticity > -100) {
		lhCurve = new FlatCurve(params->labCurve.lhcurve);
		if (!lhCurve || lhCurve->isIdentity()) {
			if (lhCurve) {
				delete lhCurve;
				lhCurve = NULL;
			}
		}//do not use "Munsell" if Chcurve not used
		else
			lhutili=true;
	}

	FlatCurve* hhCurve = NULL;//curve H=f(H)
	bool hhutili = false;
	if (params->labCurve.chromaticity > -100) {
		hhCurve = new FlatCurve(params->labCurve.hhcurve);
		if (!hhCurve || hhCurve->isIdentity()) {
			if (hhCurve) {
				delete hhCurve;
				hhCurve = NULL;
			}
		}//do not use "Munsell" if Chcurve not used
		else
			hhutili = true;
	}

	LUTf dCcurve;
	LUTf dLcurve;

	LUTu hist16Clad;
	LUTu hist16CLlad;
	LUTu hist16LLClad;

	bool chrop=false;
	float val;
	//preparate for histograms CIECAM
	if(pW!=1){//only with improccoordinator
		dCcurve(65536,0);
		dLcurve(65536,0);
		hist16Clad(65536);
		hist16CLlad(65536);
		hist16LLClad(65536);
		chrop = true;
		for (int i=0; i<48000; i++) {  //# 32768*1.414  approximation maxi for chroma
			val = (double)i / 47999.0;
			dCcurve[i] = CLIPD(val);
		}
		for (int i=0; i<65535; i++) {  //  a
			val = (double)i / 65534.0;
			dLcurve[i] = CLIPD(val);
		}

		hist16Clad.clear();
		hist16CLlad.clear();
		hist16LLClad.clear();

	}
#ifdef _DEBUG
	MyTime t1e,t2e;
	t1e.set();
	// init variables to display Munsell corrections
	MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif


	float adjustr=1.0f, adjustbg=1.0f;

//	if(params->labCurve.avoidclip ){
	// parameter to adapt curve C=f(C) to gamut

	if      (params->icm.working=="ProPhoto")   {adjustr =       adjustbg = 1.2f;}// 1.2 instead 1.0 because it's very rare to have C>170..
	else if (params->icm.working=="Adobe RGB")  {adjustr = 1.8f; adjustbg = 1.4f;}
	else if (params->icm.working=="sRGB")  	    {adjustr = 2.0f; adjustbg = 1.7f;}
	else if (params->icm.working=="WideGamut")  {adjustr =       adjustbg = 1.2f;}
	else if (params->icm.working=="Beta RGB")   {adjustr =       adjustbg = 1.4f;}
	else if (params->icm.working=="BestRGB")    {adjustr =       adjustbg = 1.4f;}
	else if (params->icm.working=="BruceRGB")   {adjustr = 1.8f; adjustbg = 1.5f;}


	// reference to the params structure has to be done outside of the parallelization to avoid CPU cache problem
	bool highlight = params->toneCurve.hrenabled; //Get the value if "highlight reconstruction" is activated
	int chromaticity = params->labCurve.chromaticity;
	bool bwToning = params->labCurve.chromaticity==-100  /*|| params->blackwhite.method=="Ch" */ || params->blackwhite.enabled;
	//if(chromaticity==-100) chromaticity==-99;
	bool LCredsk = params->labCurve.lcredsk;
	bool ccut = ccutili;
	bool clut = clcutili;
	double rstprotection = 100.-params->labCurve.rstprotection; // Red and Skin Tones Protection
	// avoid color shift is disabled when bwToning is activated and enabled if gamut is true in colorappearanace
	bool avoidColorShift = (params->labCurve.avoidcolorshift || (params->colorappearance.gamut && params->colorappearance.enabled)) && !bwToning ;
	int protectRed = settings->protectred;
	double protectRedH = settings->protectredh;
	bool gamutLch = settings->gamutLch;
	float amountchroma = (float) settings->amchroma;
	// only if user activate Lab adjustements
	if (avoidColorShift) {
		if(autili || butili || ccutili ||  cclutili || chutili || lhutili || hhutili || clcutili || utili || chromaticity) {
			unsigned int N = W*H;
			float *L = lold->L[0];
			float *a=  lold->a[0];
			float *b=  lold->b[0];
			float* Lold = new float [N];//to save L before any used
			float* Cold = new float [N];//to save C before any used
#ifdef _OPENMP
#pragma omp parallel for if (multiThread)
#endif // _OPENMP
			for (unsigned int j=0; j<N; j++){
				Lold[j]=L[j]/327.68f;
				Cold[j]=sqrt(SQR(a[j]/327.68f)+SQR(b[j]/327.68f));
			}
			Color::LabGamutMunsell(lold, Lold, Cold, /*corMunsell*/true, /*lumaMuns*/false, params->toneCurve.hrenabled, /*gamut*/true, params->icm.working, multiThread);
			delete [] Lold;
			delete [] Cold;
		}
	}

#ifdef _DEBUG
#pragma omp parallel default(shared) firstprivate(highlight, ccut, clut, chromaticity, bwToning, rstprotection, avoidColorShift, LCredsk, protectRed, protectRedH, gamutLch, lold, lnew, MunsDebugInfo, pW) if (multiThread)
#else
#pragma omp parallel default(shared) firstprivate(highlight, ccut, clut, chromaticity, bwToning, rstprotection, avoidColorShift, LCredsk, protectRed, protectRedH, gamutLch, lold, lnew, pW) if (multiThread)
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

	double wp[3][3] = {
		{wprof[0][0],wprof[0][1],wprof[0][2]},
		{wprof[1][0],wprof[1][1],wprof[1][2]},
		{wprof[2][0],wprof[2][1],wprof[2][2]}};

#pragma omp for schedule(dynamic, 16)
	for (int i=0; i<H; i++)
		for (int j=0; j<W; j++) {
			float LL=lold->L[i][j]/327.68f;
			float CC=sqrt(SQR(lold->a[i][j]) + SQR(lold->b[i][j]))/327.68f;
			float HH=xatan2f(lold->b[i][j],lold->a[i][j]);
			// According to mathematical laws we can get the sin and cos of HH by simple operations
			float2  sincosval;
			if(CC==0.0f) {
				sincosval.y = 1.0f;
				sincosval.x = 0.0f;
			} else {
				sincosval.y = lold->a[i][j]/(CC*327.68f);
				sincosval.x = lold->b[i][j]/(CC*327.68f);
			}

			float Chprov=CC;
			float Chprov1=CC;
			float memChprov=Chprov;
			
			float Lin=lold->L[i][j];
			float Lprov2=Lin/327.68f;
			lnew->L[i][j] = curve[Lin];
			float Lprov1=(lnew->L[i][j])/327.68f;
			float chromaChfactor=1.0f;
			float atmp = acurve[lold->a[i][j]+32768.0f]-32768.0f;// curves Lab a
			float btmp = bcurve[lold->b[i][j]+32768.0f]-32768.0f;// curves Lab b
			float Chprov2=Chprov1;
			int poscc,posp,posl;
			bool inRGB;
			const float ClipLevel = 65535.0f;

			if (lhutili) {  // L=f(H)
				float l_r;//Luminance Lab in 0..1
				l_r = Lprov1/100.f;	
				{
				double lr;
				float khue=1.9f;//in reserve in case of!
				float valparam = float((lhCurve->getVal(lr=Color::huelab_to_huehsv2(HH))-0.5f));//get l_r=f(H)
				float valparamneg;
				valparamneg=valparam;		
				float kcc=(CC/amountchroma);//take Chroma into account...40 "middle low" of chromaticity (arbitrary and simple), one can imagine other algorithme
				//reduce action for low chroma and increase action for high chroma 
				valparam *= 2.f*kcc; 
				valparamneg*= kcc;//slightly different for negative
				if(valparam > 0.f)
					l_r = (1.f-valparam)*l_r+ valparam*(1.f-SQR(((SQR(1.f-min(l_r,1.0f))))));
				else
					//for negative
					l_r *= (1.f+khue*valparamneg);
				}

				Lprov1=l_r*100.f;

				Chprov2 = sqrt(SQR(atmp/327.68f)+SQR(btmp/327.68f));
				//Gamut control especialy fot negative values slightly different of gamutlchonly
				do {
					inRGB=true;	
					float aprov1=Chprov2*sincosval.y;
					float bprov1=Chprov2*sincosval.x;

					float fy = (0.00862069f *Lprov1 )+ 0.137932f;
					float fx = (0.002f * aprov1) + fy;
					float fz = fy - (0.005f * bprov1);

					float x_ = 65535.0f * Color::f2xyz(fx)*Color::D50x;
					float z_ = 65535.0f * Color::f2xyz(fz)*Color::D50z;
					float y_=(Lprov1>Color::epskap) ? 65535.0*fy*fy*fy : 65535.0*Lprov1/Color::kappa;
					float R,G,B;
					Color::xyz2rgb(x_,y_,z_,R,G,B,wip);
					if (R<0.0f || G<0.0f || B<0.0f) {
						if(Lprov1 < 0.1f) Lprov1=0.1f;
						Chprov2*=0.95f;
						inRGB=false;
					}
					else if (!highlight && (R>ClipLevel || G>ClipLevel || B>ClipLevel)) {
						if (Lprov1 > 99.999f) Lprov1 = 99.98f;
						Chprov2 *= 0.95f;
						inRGB = false;
					}	
				}
				while (!inRGB)	;
					
//				float2  sincosval = xsincosf(HH);
				atmp=327.68f*Chprov2*sincosval.y;
				btmp=327.68f*Chprov2*sincosval.x;
			}

			if (editID == EUID_Lab_LCurve)
				editWhatever->v(i,j) = LIM01<float>(Lprov1/100.f);

//			calculate C=f(H)
			if (chutili) {
				double hr = Color::huelab_to_huehsv2(HH);
				float chparam = float((chCurve->getVal(hr)-0.5f) * 2.0f);//get C=f(H)
				chromaChfactor=1.0f+chparam;
			}
			
			atmp *= chromaChfactor;//apply C=f(H)
			btmp *= chromaChfactor;
			
			if (hhutili) {  // H=f(H)
				//hue Lab in -PI +PI
				double hr;
				float valparam = float((hhCurve->getVal(hr=Color::huelab_to_huehsv2(HH))-0.5f) * 1.7f) +HH;//get H=f(H)  1.7 optimisation !
				HH = valparam;
			}
			
			//simulate very approximative gamut f(L) : with pyramid transition
			float dred=55.0f;//C red value limit
			if     (Lprov1<25.0f)   dred = 40.0f;
			else if(Lprov1<30.0f)   dred = 3.0f*Lprov1 -35.0f;
			else if(Lprov1<70.0f)   dred = 55.0f;
			else if(Lprov1<75.0f)   dred = -3.0f*Lprov1 +265.0f;
			else                    dred = 40.0f;
			// end pyramid
			if(params->dirpyrDenoise.enabled && chromaticity ==0) chromaticity = 0.5f;
				
			if(!bwToning){
				float factorskin, factorsat, factor, factorskinext, interm;
				float scale = 100.0f/100.1f;//reduction in normal zone
				float scaleext=1.0f;//reduction in transition zone
				float protect_red,protect_redh;
				float deltaHH;//HH value transition
				protect_red=float(protectRed);//default=60  chroma: one can put more or less if necessary...in 'option'  40...160
				if(protect_red < 20.0f) protect_red=20.0; // avoid too low value
				if(protect_red > 180.0f) protect_red=180.0; // avoid too high value
				protect_redh=float(protectRedH);//default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0
				if(protect_redh<0.1f) protect_redh=0.1f;//avoid divide by 0 and negatives values
				if(protect_redh>1.0f) protect_redh=1.0f;//avoid too big values

				deltaHH=protect_redh;//transition hue
				float chromapro	= (chromaticity	+ 100.0f)/100.0f;
				if(chromapro>0.0) Color::scalered ( rstprotection, chromapro, 0.0, HH, deltaHH, scale, scaleext);//1.0
				if(chromapro>1.0) {
					interm=(chromapro-1.0f)*100.0f;
					factorskin= 1.0f+(interm*scale)/100.0f;
					factorskinext=1.0f+(interm*scaleext)/100.0f;
				}
				else {
					//factorskin= chromapro*scale;
					//factorskinext= chromapro*scaleext;
					factorskin= chromapro ; // +(chromapro)*scale;
					factorskinext= chromapro ;// +(chromapro)*scaleext;

					}
				factorsat=chromapro;
				//increase saturation after denoise : ...approximation
				float factnoise=1.f;
				if(params->dirpyrDenoise.enabled) {
					factnoise=(1.f+params->dirpyrDenoise.chroma/500.f);//levels=5
					
					
			//		if(yyyy) factnoise=(1.f+params->dirpyrDenoise.chroma/100.f);//levels=7
				}
				factorsat*=factnoise;
				
				factor=factorsat;
				// Test if chroma is in the normal range first
				Color::transitred ( HH, Chprov1, dred, factorskin, protect_red, factorskinext, deltaHH, factorsat, factor);
				atmp *= factor;
				btmp *= factor;

				if (clut) { // begin C=f(L)
					float factorskin,factorsat,factor,factorskinext,interm;
					float chroma = sqrt(SQR(atmp)+SQR(btmp)+0.001f);
				//	float chromaCfactor=(clcurve[Lprov1*327.68f])/(Lprov1*327.68f);//apply C=f(L)
					float chromaCfactor=(clcurve[Lprov2*327.68f])/(Lprov2*327.68f);//apply C=f(L)
					float curf=0.7f;//empirical coeff because curve is more progressive
					float scale = 100.0f/100.1f;//reduction in normal zone for curve C
					float scaleext=1.0f;//reduction in transition zone for curve C
					float protect_redcur,protect_redhcur;//perhaps the same value than protect_red and protect_redh
					float deltaHH;//HH value transition for C curve
					protect_redcur=curf*float(protectRed);//default=60  chroma: one can put more or less if necessary...in 'option'  40...160==> curf =because curve is more progressive
					if(protect_redcur < 20.0f) protect_redcur=20.0; // avoid too low value
					if(protect_redcur > 180.0f) protect_redcur=180.0; // avoid too high value
					protect_redhcur=curf*float(protectRedH);//default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0 ==> curf =because curve is more progressive
					if(protect_redhcur<0.1f) protect_redhcur=0.1f;//avoid divide by 0 and negatives values
					if(protect_redhcur>1.0f) protect_redhcur=1.0f;//avoid too big values

					deltaHH=protect_redhcur;//transition hue
					if(chromaCfactor>0.0) Color::scalered ( rstprotection, chromaCfactor, 0.0, HH, deltaHH, scale, scaleext);//1.0
					if(chromaCfactor>1.0) {
						interm=(chromaCfactor-1.0f)*100.0f;
						factorskin= 1.0f+(interm*scale)/100.0f;
						factorskinext=1.0f+(interm*scaleext)/100.0f;
						}
					else {
						factorskin= chromaCfactor; // +(1.0f-chromaCfactor)*scale;
						factorskinext= chromaCfactor ; //+(1.0f-chromaCfactor)*scaleext;
						}

					factorsat=chromaCfactor;
					factor=factorsat;
					Color::transitred ( HH, Chprov1, dred, factorskin, protect_redcur, factorskinext, deltaHH, factorsat, factor);
					atmp *= factor;
					btmp *= factor;
				}
				// end C=f(L)

				// I have placed C=f(C) after all C treatments to assure maximum amplitude of "C"
				if (ccut) {
					float factorskin,factorsat,factor,factorskinext,interm;
					float chroma = sqrt(SQR(atmp)+SQR(btmp)+0.001f);
					float chromaCfactor=(satcurve[chroma*adjustr])/(chroma*adjustr);//apply C=f(C)
					float curf=0.7f;//empirical coeff because curve is more progressive
					float scale = 100.0f/100.1f;//reduction in normal zone for curve CC
					float scaleext=1.0f;//reduction in transition zone for curve CC
					float protect_redcur,protect_redhcur;//perhaps the same value than protect_red and protect_redh
					float deltaHH;//HH value transition for CC curve
					protect_redcur=curf*float(protectRed);//default=60  chroma: one can put more or less if necessary...in 'option'  40...160==> curf =because curve is more progressive
					if(protect_redcur < 20.0f) protect_redcur=20.0; // avoid too low value
					if(protect_redcur > 180.0f) protect_redcur=180.0; // avoid too high value
					protect_redhcur=curf*float(protectRedH);//default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0 ==> curf =because curve is more progressive
					if(protect_redhcur<0.1f) protect_redhcur=0.1f;//avoid divide by 0 and negatives values
					if(protect_redhcur>1.0f) protect_redhcur=1.0f;//avoid too big values

					deltaHH=protect_redhcur;//transition hue
					if(chromaCfactor>0.0) Color::scalered ( rstprotection, chromaCfactor, 0.0, HH, deltaHH, scale, scaleext);//1.0
					if(chromaCfactor>1.0) {
						interm=(chromaCfactor-1.0f)*100.0f;
						factorskin= 1.0f+(interm*scale)/100.0f;
						factorskinext=1.0f+(interm*scaleext)/100.0f;
						}
					else {
						//factorskin= chromaCfactor*scale;
						//factorskinext=chromaCfactor*scaleext;
						factorskin= chromaCfactor; // +(1.0f-chromaCfactor)*scale;
						factorskinext= chromaCfactor ; //+(1.0f-chromaCfactor)*scaleext;

						}

					factorsat=chromaCfactor;
					factor=factorsat;
					Color::transitred ( HH, Chprov1, dred, factorskin, protect_redcur, factorskinext, deltaHH, factorsat, factor);
					atmp *= factor;
					btmp *= factor;
				}
			}
			// end chroma C=f(C)

			//update histogram C
			if(pW!=1){//only with improccoordinator
				posp=CLIP((int)sqrt((atmp*atmp + btmp*btmp)));
				hist16Clad[posp]++;
				hist16CLlad[posp]++;
			}

			if (!bwToning) {	//apply curve L=f(C) for skin and rd...but also for extended color ==> near green and blue (see 'curf')

				const float xx=0.25f;//soft : between 0.2 and 0.4
				float protect_redhcur;
				float curf=1.0f;
				float deltaHH;
				protect_redhcur=curf*float(protectRedH);//default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1
				if(protect_redhcur<0.1f) protect_redhcur=0.1f;//avoid divide by 0 and negatives values:minimal protection for transition
				if(protect_redhcur>3.5f) protect_redhcur=3.5f;//avoid too big values

				deltaHH=protect_redhcur;//transition hue

				float skbeg=-0.05f;//begin hue skin
				float skend=1.60f;//end hue skin
				const float chrmin=50.0f;//to avoid artifact, because L curve is not a real curve for luminance
				float aa,bb;
				float zz=0.0f;
				float yy=0.0f;
				if(Chprov1 < chrmin) yy=(Chprov1/chrmin)*(Chprov1/chrmin)*xx;else yy=xx;//avoid artifact for low C
				if(!LCredsk) {skbeg=-3.1415; skend=3.14159; deltaHH=0.001f;}
					if(HH>skbeg && HH < skend ) zz=yy;
					else if(HH>skbeg-deltaHH && HH<=skbeg) {aa=yy/deltaHH;bb=-aa*(skbeg-deltaHH); zz=aa*HH+bb;}//transition
					else if(HH>=skend && HH < skend+deltaHH) {aa=-yy/deltaHH;bb=-aa*(skend+deltaHH);zz=aa*HH+bb;}//transition

				float chroma=sqrt(SQR(atmp)+SQR(btmp)+0.001f);
				float Lc = (lhskcurve[chroma*adjustr])/(chroma*adjustr);//apply L=f(C)
				Lc=(Lc-1.0f)*zz+1.0f;//reduct action
				Lprov1*=Lc;//adjust luminance
				}
				//update histo L
                if(pW!=1){//only with improccoordinator
					posl=CLIP((int(Lprov1*327.68f)));
					hist16LLClad[posl]++;
				}

			Chprov1 = sqrt(SQR(atmp/327.68f)+SQR(btmp/327.68f));

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
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif
					lnew->L[i][j]=Lprov1*327.68f;
//					float2 sincosval = xsincosf(HH);
					lnew->a[i][j]=327.68f*Chprov1*sincosval.y;
					lnew->b[i][j]=327.68f*Chprov1*sincosval.x;
				}
				else {
					//use gamutbdy
					//Luv limiter
					float Y,u,v;
					Color::Lab2Yuv(lnew->L[i][j],atmp,btmp,Y,u,v);
					//Yuv2Lab includes gamut restriction map
					Color::Yuv2Lab(Y,u,v,lnew->L[i][j],lnew->a[i][j],lnew->b[i][j], wp);
				}

				if (utili || autili || butili || ccut || clut || cclutili || chutili || lhutili || hhutili || clcutili || chromaticity) {
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
			/*		if((HH>0.0f && HH < 1.6f)	&& memChprov < 70.0f) HH+=correctlum;//skin correct
					else if(fabs(correctionHue) < 0.3f) HH+=0.08f*correctlum;
					else if(fabs(correctionHue) < 0.2f) HH+=0.25f*correctlum;
					else if(fabs(correctionHue) < 0.1f) HH+=0.35f*correctlum;
					else if(fabs(correctionHue) < 0.015f) HH+=correctlum;	// correct only if correct Munsell chroma very little.
			*/
					float2 sincosval = xsincosf(HH+correctionHue);
			
					lnew->a[i][j]=327.68f*Chprov*sincosval.y;// apply Munsell
					lnew->b[i][j]=327.68f*Chprov*sincosval.x;
				}
			}
			else {
//				if(Lprov1 > maxlp) maxlp=Lprov1;
//				if(Lprov1 < minlp) minlp=Lprov1;
				if(!bwToning){
					lnew->L[i][j]=Lprov1*327.68f;
//					float2 sincosval = xsincosf(HH);
					lnew->a[i][j]=327.68f*Chprov1*sincosval.y;
					lnew->b[i][j]=327.68f*Chprov1*sincosval.x;
				}
				else {
					//Luv limiter only
					lnew->a[i][j] = atmp;
					lnew->b[i][j] = btmp;
				}
			}
		}

} // end of parallelization

    //update histogram C  with data chromaticity and not with CC curve
	if(pW!=1){//only with improccoordinator
		for (int i=0; i<=48000; i++) {//32768*1.414  + ...
			if (chrop) {
				float hval = dCcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				histCCurve[hi] += hist16Clad[i] ;
				histCLurve[hi] += hist16CLlad[i] ;
			}
		}
		   //update histogram L  with data luminance
		for (int i=0; i<=65535; i++) {
			if (chrop) {
				float hlval = dLcurve[i];
				int hli = (int)(255.0*CLIPD(hlval));
				histLLCurve[hli] += hist16LLClad[i] ;
				histLCurve[hli] += hist16LLClad[i] ;
				
			}
		}
		
	}

#ifdef _DEBUG
	if (settings->verbose) {
		t2e.set();
		printf("Color::AllMunsellLch (correction performed in %d usec):\n", t2e.etime(t1e));
		printf("   Munsell chrominance: MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
		printf("   Munsell luminance  : MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhuelum[0], MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
	}
	delete MunsDebugInfo;
#endif

	if (chCurve) delete chCurve;
	if (lhCurve) delete lhCurve;
	if (hhCurve) delete hhCurve;
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

	void ImProcFunctions::impulsedenoisecam (CieImage* ncie) {

		if (params->impulseDenoise.enabled && ncie->W>=8 && ncie->H>=8)

			impulse_nrcam (ncie, (float)params->impulseDenoise.thresh/20.0 );
	}
	
	void ImProcFunctions::defringe (LabImage* lab) {

		if (params->defringe.enabled && lab->W>=8 && lab->H>=8)

			PF_correct_RT(lab, lab, params->defringe.radius, params->defringe.threshold);
	}

	void ImProcFunctions::defringecam (CieImage* ncie) {
		if (params->defringe.enabled && ncie->W>=8 && ncie->H>=8) PF_correct_RTcam(ncie, ncie, params->defringe.radius, params->defringe.threshold);
		
	}
	
	void ImProcFunctions::badpixcam(CieImage* ncie, double rad, int thr, int mode){
		if(ncie->W>=8 && ncie->H>=8) Badpixelscam(ncie, ncie, rad, thr, mode);
	}

	void ImProcFunctions::dirpyrequalizer (LabImage* lab) {

		if (params->dirpyrequalizer.enabled && lab->W>=8 && lab->H>=8) {

			//dirpyrLab_equalizer(lab, lab, params->dirpyrequalizer.mult);
			dirpyr_equalizer(lab->L, lab->L, lab->W, lab->H, params->dirpyrequalizer.mult, params->dirpyrequalizer.threshold);
		}
	}
void ImProcFunctions::EPDToneMapCIE(CieImage *ncie, float a_w, float c_, float w_h, int Wid, int Hei, int begh, int endh, float minQ, float maxQ, unsigned int Iterates, int skip){

if(!params->edgePreservingDecompositionUI.enabled) return;
		float stren=params->edgePreservingDecompositionUI.Strength;
		float edgest=params->edgePreservingDecompositionUI.EdgeStopping;
		float sca=params->edgePreservingDecompositionUI.Scale;
		float rew=params->edgePreservingDecompositionUI.ReweightingIterates;
		unsigned int i, N = Wid*Hei;
		float Qpro= ( 4.0 / c_)  * ( a_w + 4.0 ) ;//estimate Q max if J=100.0
		float *Qpr=ncie->Q_p[0];
		float eps=0.0001;
		if (settings->verbose) printf("minQ=%f maxQ=%f  Qpro=%f\n",minQ,maxQ, Qpro);
		if(maxQ>Qpro) Qpro=maxQ;
		for (int i=0; i<Hei; i++)
			for (int j=0; j<Wid; j++) { Qpr[i*Wid+j]=ncie->Q_p[i][j];}

		EdgePreservingDecomposition epd = EdgePreservingDecomposition(Wid, Hei);

		for(i = 0; i != N; i++) Qpr[i] = (Qpr[i]+eps)/(Qpro);

		float Compression = expf(-stren);		//This modification turns numbers symmetric around 0 into exponents.
		float DetailBoost = stren;
		if(stren < 0.0f) DetailBoost = 0.0f;	//Go with effect of exponent only if uncompressing.

		//Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
		if(Iterates == 0) Iterates = (unsigned int)(edgest*15.0);
		//Jacques Desmis : always Iterates=5 for compatibility images between preview and output

		epd.CompressDynamicRange(Qpr, (float)sca/skip, (float)edgest, Compression, DetailBoost, Iterates, rew, Qpr);

		//Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
		float s = (1.0f + 38.7889f)*powf(Compression, 1.5856f)/(1.0f + 38.7889f*powf(Compression, 1.5856f));
		#ifndef _DEBUG
		#pragma omp parallel for schedule(dynamic,10)
		#endif
		for (int i=0; i<Hei; i++)
			for (int j=0; j<Wid; j++) {
			ncie->Q_p[i][j]=(Qpr[i*Wid+j]+eps)*Qpro;
			ncie->M_p[i][j]*=s;
		}
/*
	float *Qpr2 = new float[Wid*((heir)+1)];

		for (int i=heir; i<Hei; i++)
			for (int j=0; j<Wid; j++) { Qpr2[(i-heir)*Wid+j]=ncie->Q_p[i][j];}
	if(minQ>0.0) minQ=0.0;//normaly minQ always > 0...
//	EdgePreservingDecomposition epd = EdgePreservingDecomposition(Wid, Hei);
//EdgePreservingDecomposition epd = EdgePreservingDecomposition(Wid, Hei/2);
	for(i = N2; i != N; i++)
//	for(i = begh*Wid; i != N; i++)
		//Qpr[i] = (Qpr[i]-minQ)/(maxQ+1.0);
		Qpr2[i-N2] = (Qpr2[i-N2]-minQ)/(Qpro+1.0);

	float Compression2 = expf(-stren);		//This modification turns numbers symmetric around 0 into exponents.
	float DetailBoost2 = stren;
	if(stren < 0.0f) DetailBoost2 = 0.0f;	//Go with effect of exponent only if uncompressing.

	//Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
	if(Iterates == 0) Iterates = (unsigned int)(edgest*15.0);


	epd.CompressDynamicRange(Qpr2, (float)sca/skip, (float)edgest, Compression2, DetailBoost2, Iterates, rew, Qpr2);

	//Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
	 float s2 = (1.0f + 38.7889f)*powf(Compression, 1.5856f)/(1.0f + 38.7889f*powf(Compression, 1.5856f));
		for (int i=heir; i<Hei; i++)
	//	for (int i=begh; i<endh; i++)
			for (int j=0; j<Wid; j++) {
			ncie->Q_p[i][j]=Qpr2[(i-heir)*Wid+j]*Qpro + minQ;
		//	Qpr[i*Wid+j]=Qpr[i*Wid+j]*maxQ + minQ;
		//	ncie->J_p[i][j]=(100.0* Qpr[i*Wid+j]*Qpr[i*Wid+j]) /(w_h*w_h);

			ncie->M_p[i][j]*=s2;
		}
				delete [] Qpr2;

*/
}


//Map tones by way of edge preserving decomposition. Is this the right way to include source?
//#include "EdgePreservingDecomposition.cc"
void ImProcFunctions::EPDToneMap(LabImage *lab, unsigned int Iterates, int skip){
	//Hasten access to the parameters.
//	EPDParams *p = (EPDParams *)(&params->edgePreservingDecompositionUI);

	//Enabled? Leave now if not.
//	if(!p->enabled) return;
if(!params->edgePreservingDecompositionUI.enabled) return;
float stren=params->edgePreservingDecompositionUI.Strength;
float edgest=params->edgePreservingDecompositionUI.EdgeStopping;
float sca=params->edgePreservingDecompositionUI.Scale;
float rew=params->edgePreservingDecompositionUI.ReweightingIterates;
	//Pointers to whole data and size of it.
	float *L = lab->L[0];
	float *a = lab->a[0];
	float *b = lab->b[0];
	unsigned int i, N = lab->W*lab->H;

	EdgePreservingDecomposition epd = EdgePreservingDecomposition(lab->W, lab->H);

	//Due to the taking of logarithms, L must be nonnegative. Further, scale to 0 to 1 using nominal range of L, 0 to 15 bit.
    float minL = FLT_MAX;
#pragma omp parallel
{
	float lminL = FLT_MAX;
#pragma omp for
	for(i = 0; i < N; i++)
		if(L[i] < lminL) lminL = L[i];
#pragma omp critical
    if(lminL < minL) minL = lminL;
}
	if(minL > 0.0f) minL = 0.0f;		//Disable the shift if there are no negative numbers. I wish there were just no negative numbers to begin with.
#pragma omp parallel for
	for(i = 0; i < N; i++)
		L[i] = (L[i] - minL)/32767.0f;

	//Some interpretations.
	float Compression = expf(-stren);		//This modification turns numbers symmetric around 0 into exponents.
	float DetailBoost = stren;
	if(stren < 0.0f) DetailBoost = 0.0f;	//Go with effect of exponent only if uncompressing.

	//Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
	if(Iterates == 0) Iterates = (unsigned int)(edgest*15.0f);

/* Debuggery. Saves L for toying with outside of RT.
char nm[64];
sprintf(nm, "%ux%ufloat.bin", lab->W, lab->H);
FILE *f = fopen(nm, "wb");
fwrite(L, N, sizeof(float), f);
fclose(f);*/

	epd.CompressDynamicRange(L, sca/float(skip), edgest, Compression, DetailBoost, Iterates, rew, L);

	//Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
	float s = (1.0f + 38.7889f)*powf(Compression, 1.5856f)/(1.0f + 38.7889f*powf(Compression, 1.5856f));
	#ifdef _OPENMP
	#pragma omp parallel for            // removed schedule(dynamic,10)
	#endif
	for(int ii = 0; ii < N; ii++)
		a[ii] *= s,
		b[ii] *= s,
		L[ii] = L[ii]*32767.0f + minL;
}


	void ImProcFunctions::getAutoExp  (LUTu & histogram, int histcompr, double defgain, double clip,
									   double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh) {

		float scale = 65536.0f;
		float midgray=0.1842f;//middle gray in linear gamma =1 50% luminance

		int imax=65536>>histcompr;
		int overex=0;
		float sum = 0.f, hisum=0.f, losum=0.f;
		float ave = 0.f, hidev=0.f, lodev=0.f;
		//find average luminance
		for (int i=0; i<imax; i++) {
			sum += histogram[i];
			ave += histogram[i] *(float)i;
		}
		ave /= (sum);

		//find median of luminance
		int median=0, count=histogram[0];
		while (count<sum/2) {
			median++;
			count += histogram[median];
		}
		if (median==0 || ave<1.f) {//probably the image is a blackframe
			expcomp=0.;
			black=0;
			bright=0;
			contr=0;
			hlcompr=0;
			hlcomprthresh=0;
			return;
		}

		// compute std dev on the high and low side of median
		// and octiles of histogram
		float octile[8]={0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f},ospread=0.f;
		count=0;
		for (int i=0; i<imax; i++) {
			if (count<8) {
				octile[count] += histogram[i];
				if (octile[count]>sum/8.f || (count==7 && octile[count]>sum/16.f)) {
					octile[count]=log(1.+(float)i)/log(2.f);
					count++;// = min(count+1,7);
				}
			}
			if (i<ave) {
				//lodev += SQR(ave-i)*histogram[i];
				lodev += (log(ave+1.f)-log((float)i+1.))*histogram[i];
				losum += histogram[i];
			} else {
				//hidev += SQR(i-ave)*histogram[i];
				hidev += (log((float)i+1.)-log(ave+1.f))*histogram[i];
				hisum += histogram[i];
			}
		}
		if (losum==0 || hisum==0) {//probably the image is a blackframe
			expcomp=0.;
			black=0;
			bright=0;
			contr=0;
			hlcompr=0;
			hlcomprthresh=0;
			return;
		}

		lodev = (lodev/(log(2.f)*losum));
		hidev = (hidev/(log(2.f)*hisum));
		if (octile[6]>log((float)imax+1.f)/log2(2.f)) {//if very overxposed image
			octile[6]=1.5f*octile[5]-0.5f*octile[4];
			overex=2;
		}
		
		if (octile[7]>log((float)imax+1.f)/log2(2.f)) {//if overexposed
			octile[7]=1.5f*octile[6]-0.5f*octile[5];
			overex=1;
		}

		// store values of octile[6] and octile[7] for calculation of exposure compensation
		// if we don't do this and the pixture is underexposed, calculation of exposure compensation assumes
		// that it's overexposed and calculates the wrong direction
		float oct6,oct7;
		oct6 = octile[6];
		oct7 = octile[7];
		

		for(int i=1; i<8; i++) {
			if (octile[i] == 0.0f)
				octile[i] = octile[i-1];
		}
		// compute weighted average separation of octiles
		// for future use in contrast setting
		for (int i=1; i<6; i++) {
			ospread += (octile[i+1]-octile[i])/max(0.5f,(i>2 ? (octile[i+1]-octile[3]) : (octile[3]-octile[i])));
		}
		ospread /= 5.f;
		if (ospread<=0.f) {//probably the image is a blackframe
			expcomp=0.;
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
		int clippable = (int)(sum * clip/100.f );
		int somm=sum;
		clipped = 0;
		int whiteclip = (imax) - 1;
		while (whiteclip>1 && histogram[whiteclip]+clipped <= clippable) {
			clipped += histogram[whiteclip];
			whiteclip--;
		}
		int clipwh=clipped;
		//compute clipped black point
		clipped = 0;
		int shc = 0;
		
		while (shc<whiteclip-1 && histogram[shc]+clipped <= clippable) {
			clipped += histogram[shc];
			shc++;
		}
		int clipbl=clipped;
		
		//rescale to 65535 max
		rawmax <<= histcompr;
		whiteclip <<= histcompr;
		ave = ave*(1<<histcompr);
		median <<= histcompr;
		shc <<= histcompr;

		//prevent division by 0
		if (lodev==0.f)
			lodev=1.f;

		//compute exposure compensation as geometric mean of the amount that
		//sets the mean or median at middle gray, and the amount that sets the estimated top
		//of the histogram at or near clipping.
		float expo=log(midgray*scale/(ave-shc+midgray*shc));
		//float expcomp1 = (log(/*(median/ave)*//*(hidev/lodev)*/midgray*scale/(ave-shc+midgray*shc))+log((hidev/lodev)))/log(2.f);
		float expcomp1 = (log(/*(median/ave)*//*(hidev/lodev)*/midgray*scale/(ave-shc+midgray*shc)))/log(2.f);
		float expcomp2;
		
		if(overex == 0) { // image is not overexposed
			expcomp2 = 0.5f*( (15.5f-histcompr-(2.f*oct7-oct6)) + log(scale/rawmax)/log(2.f) );
		}
		else {
			expcomp2 = 0.5f*( (15.5f-histcompr-(2.f*octile[7]-octile[6])) + log(scale/rawmax)/log(2.f) );
		}
		
		if(fabs(expcomp1)-fabs(expcomp2)> 1.f) {//for great expcomp
			expcomp = (expcomp1*fabs(expcomp2)+expcomp2*fabs(expcomp1))/(fabs(expcomp1)+fabs(expcomp2));
		}
		else {
			expcomp = 0.5 * (double)expcomp1 + 0.5 *(double) expcomp2;//for small expcomp
		}
		float gain = exp((float)expcomp*log(2.f));

		float corr = sqrt(gain*scale/rawmax);
		black = (int) shc*corr;


		//now tune hlcompr to bring back rawmax to 65535
		hlcomprthresh = 33;
		//this is a series approximation of the actual formula for comp,
		//which is a transcendental equation
		float comp = (gain*((float)whiteclip)/scale - 1.f)*2.3f;// 2.3 instead of 2 to increase slightly comp
		hlcompr=(int)(100.*comp/(max(0.0,expcomp) + 1.0));
		hlcompr = max(0,min(100,hlcompr));

		//now find brightness if gain didn't bring ave to midgray using
		//the envelope of the actual 'control cage' brightness curve for simplicity
		float midtmp = gain*sqrt(median*ave)/scale;
		if (midtmp<0.1f) {
			bright = (midgray-midtmp)*15.0/(midtmp);
		} else {
			bright = (midgray-midtmp)*15.0/(0.10833-0.0833*midtmp);
		}

		bright = 0.25*/*(median/ave)*(hidev/lodev)*/max(0,bright);

		//compute contrast that spreads the average spacing of octiles
		contr = (int) 50.0f*(1.1f-ospread);
		contr = max(0,min(100,contr));
		//take gamma into account
		double whiteclipg = (int)(CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);
		
		double gavg = 0.;
		for (int i=0; i<65536>>histcompr; i++)
		gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;
		if (black < gavg) {
			int maxwhiteclip = (gavg - black) * 4 / 3 + black; // dont let whiteclip be such large that the histogram average goes above 3/4
			if (whiteclipg < maxwhiteclip)
				whiteclipg = maxwhiteclip;
		 }
		whiteclipg = CurveFactory::igamma2 ((float)(whiteclipg/65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter
		
		//corection with gamma
		black = (int)((65535*black)/whiteclipg);
		//expcomp = log(65535.0 / (whiteclipg)) / log(2.0);

		//diagnostics
		//printf ("**************** AUTO LEVELS ****************\n");
		
		if (settings->verbose) {
			printf("sum=%i clip=%f clippable=%i  clipWh=%i  clipBl=%i\n",somm, clip, clippable,clipwh, clipbl);
			printf ("expcomp1= %f   expcomp2= %f gain= %f  expcomp=%f\n",expcomp1,expcomp2,gain,expcomp);
			printf ("expo=%f\n",expo);
			printf ("median: %i  average: %f    median/average: %f\n",median,ave, median/ave);
			printf ("average: %f\n",ave);
			printf("comp=%f hlcomp: %i\n",comp, hlcompr);
			printf ("median/average: %f\n",median/ave);
			printf ("lodev: %f   hidev: %f		hidev/lodev: %f\n",lodev,hidev,hidev/lodev);
			printf ("lodev: %f\n",lodev);
			printf ("hidev: %f\n",hidev);
			printf ("imax=%d rawmax= %d  whiteclip= %d  gain= %f\n",imax,rawmax,whiteclip,gain);
			printf ("octiles: %f %f %f %f %f %f %f %f\n",octile[0],octile[1],octile[2],octile[3],octile[4],octile[5],octile[6],octile[7]);
			printf ("ospread= %f\n",ospread);
			printf ("overexp= %i\n",overex);
		}

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
		if (expcomp>12.0)	expcomp = 12.0;
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

			Thumbnail* raw =   rtengine::Thumbnail::loadFromRaw      (fname, ri, w_raw, h_raw, 1, 1.0, FALSE);
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
#undef PIX_SORT
#undef med3x3
