////////////////////////////////////////////////////////////////
//
//			CFA denoise by wavelet transform, FT filtering
//
//	copyright (c) 2008-2012  Emil Martinec <ejmartin@uchicago.edu>
//
//
//  code dated: March 9, 2012
//
//	FTblockDN.cc is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include <math.h>
#include <fftw3.h>
#include "../rtgui/threadutils.h"

#include "rtengine.h"
#include "improcfun.h"
#include "LUT.h"
#include "array2D.h"
#include "iccmatrices.h"
#include "boxblur.h"
#include "rt_math.h"
#include "mytime.h"
#include "sleef.c"
#include "opthelper.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cplx_wavelet_dec.h"

//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
//#define LIM(x,min,max) MAX(min,MIN(x,max))
//#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
//#define CLIP(x) LIM(x,0,65535)

#define TS 64		// Tile size
#define offset 25	// shift between tiles
#define fTS ((TS/2+1))	// second dimension of Fourier tiles
#define blkrad 1	// radius of block averaging

#define epsilon 0.001f/(TS*TS) //tolerance
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

#define med2(a0,a1,a2,a3,a4,median) { \
pp[0]=a0; pp[1]=a1; pp[2]=a2; pp[3]=a3; pp[4]=a4;  \
PIX_SORT(pp[0],pp[1]) ; PIX_SORT(pp[3],pp[4]) ; PIX_SORT(pp[0],pp[3]) ;\
PIX_SORT(pp[1],pp[4]) ; PIX_SORT(pp[1],pp[2]) ; PIX_SORT(pp[2],pp[3]) ;\
PIX_SORT(pp[1],pp[2]) ; median=pp[2] ;}

#define med5(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,median) { \
pp[0]=a0; pp[1]=a1; pp[2]=a2; pp[3]=a3; pp[4]=a4; pp[5]=a5; pp[6]=a6; pp[7]=a7; pp[8]=a8; pp[9]=a9; pp[10]=a10; pp[11]=a11; pp[12]=a12; \
pp[13]=a13; pp[14]=a14; pp[15]=a15; pp[16]=a16; pp[17]=a17; pp[18]=a18; pp[19]=a19; pp[20]=a20; pp[21]=a21; pp[22]=a22; pp[23]=a23; pp[24]=a24; \
PIX_SORT(pp[0], pp[1]) ;   PIX_SORT(pp[3], pp[4]) ;   PIX_SORT(pp[2], pp[4]) ;\
PIX_SORT(pp[2], pp[3]) ;   PIX_SORT(pp[6], pp[7]) ;   PIX_SORT(pp[5], pp[7]) ;\
PIX_SORT(pp[5], pp[6]) ;   PIX_SORT(pp[9], pp[10]) ;  PIX_SORT(pp[8], pp[10]) ;\
PIX_SORT(pp[8], pp[9]) ;   PIX_SORT(pp[12], pp[13]) ; PIX_SORT(pp[11], pp[13]) ;\
PIX_SORT(pp[11], pp[12]) ; PIX_SORT(pp[15], pp[16]) ; PIX_SORT(pp[14], pp[16]) ;\
PIX_SORT(pp[14], pp[15]) ; PIX_SORT(pp[18], pp[19]) ; PIX_SORT(pp[17], pp[19]) ;\
PIX_SORT(pp[17], pp[18]) ; PIX_SORT(pp[21], pp[22]) ; PIX_SORT(pp[20], pp[22]) ;\
PIX_SORT(pp[20], pp[21]) ; PIX_SORT(pp[23], pp[24]) ; PIX_SORT(pp[2], pp[5]) ;\
PIX_SORT(pp[3], pp[6]) ;   PIX_SORT(pp[0], pp[6]) ;   PIX_SORT(pp[0], pp[3]) ;\
PIX_SORT(pp[4], pp[7]) ;   PIX_SORT(pp[1], pp[7]) ;   PIX_SORT(pp[1], pp[4]) ;\
PIX_SORT(pp[11], pp[14]) ; PIX_SORT(pp[8], pp[14]) ;  PIX_SORT(pp[8], pp[11]) ;\
PIX_SORT(pp[12], pp[15]) ; PIX_SORT(pp[9], pp[15]) ;  PIX_SORT(pp[9], pp[12]) ;\
PIX_SORT(pp[13], pp[16]) ; PIX_SORT(pp[10], pp[16]) ; PIX_SORT(pp[10], pp[13]) ;\
PIX_SORT(pp[20], pp[23]) ; PIX_SORT(pp[17], pp[23]) ; PIX_SORT(pp[17], pp[20]) ;\
PIX_SORT(pp[21], pp[24]) ; PIX_SORT(pp[18], pp[24]) ; PIX_SORT(pp[18], pp[21]) ;\
PIX_SORT(pp[19], pp[22]) ; PIX_SORT(pp[8], pp[17]) ;  PIX_SORT(pp[9], pp[18]) ;\
PIX_SORT(pp[0], pp[18]) ;  PIX_SORT(pp[0], pp[9]) ;   PIX_SORT(pp[10], pp[19]) ;\
PIX_SORT(pp[1], pp[19]) ;  PIX_SORT(pp[1], pp[10]) ;  PIX_SORT(pp[11], pp[20]) ;\
PIX_SORT(pp[2], pp[20]) ;  PIX_SORT(pp[2], pp[11]) ;  PIX_SORT(pp[12], pp[21]) ;\
PIX_SORT(pp[3], pp[21]) ;  PIX_SORT(pp[3], pp[12]) ;  PIX_SORT(pp[13], pp[22]) ;\
PIX_SORT(pp[4], pp[22]) ;  PIX_SORT(pp[4], pp[13]) ;  PIX_SORT(pp[14], pp[23]) ;\
PIX_SORT(pp[5], pp[23]) ;  PIX_SORT(pp[5], pp[14]) ;  PIX_SORT(pp[15], pp[24]) ;\
PIX_SORT(pp[6], pp[24]) ;  PIX_SORT(pp[6], pp[15]) ;  PIX_SORT(pp[7], pp[16]) ;\
PIX_SORT(pp[7], pp[19]) ;  PIX_SORT(pp[13], pp[21]) ; PIX_SORT(pp[15], pp[23]) ;\
PIX_SORT(pp[7], pp[13]) ;  PIX_SORT(pp[7], pp[15]) ;  PIX_SORT(pp[1], pp[9]) ;\
PIX_SORT(pp[3], pp[11]) ;  PIX_SORT(pp[5], pp[17]) ;  PIX_SORT(pp[11], pp[17]) ;\
PIX_SORT(pp[9], pp[17]) ;  PIX_SORT(pp[4], pp[10]) ;  PIX_SORT(pp[6], pp[12]) ;\
PIX_SORT(pp[7], pp[14]) ;  PIX_SORT(pp[4], pp[6]) ;   PIX_SORT(pp[4], pp[7]) ;\
PIX_SORT(pp[12], pp[14]) ; PIX_SORT(pp[10], pp[14]) ; PIX_SORT(pp[6], pp[7]) ;\
PIX_SORT(pp[10], pp[12]) ; PIX_SORT(pp[6], pp[10]) ;  PIX_SORT(pp[6], pp[17]) ;\
PIX_SORT(pp[12], pp[17]) ; PIX_SORT(pp[7], pp[17]) ;  PIX_SORT(pp[7], pp[10]) ;\
PIX_SORT(pp[12], pp[18]) ; PIX_SORT(pp[7], pp[12]) ;  PIX_SORT(pp[10], pp[18]) ;\
PIX_SORT(pp[12], pp[20]) ; PIX_SORT(pp[10], pp[20]) ; PIX_SORT(pp[10], pp[12]) ;\
median=pp[12];} 

#define ELEM_FLOAT_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }


namespace rtengine {

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*
	 Structure of the algorithm:

	 1. Compute an initial denoise of the image via undecimated wavelet transform
	 and universal thresholding modulated by user input.
	 2. Decompose the residual image into TSxTS size tiles, shifting by 'offset' each step
	 (so roughly each pixel is in (TS/offset)^2 tiles); Discrete Cosine transform the tiles.
	 3. Filter the DCT data to pick out patterns missed by the wavelet denoise
	 4. Inverse DCT the denoised tile data and combine the tiles into a denoised output image.

	 */

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


extern const Settings* settings;

// Median calculation using quicksort
float fq_sort2(float arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low)
            return arr[median] ;
        if (high == low + 1) {
            if (arr[low] > arr[high])
                ELEM_FLOAT_SWAP(arr[low], arr[high]) ;
           return arr[median] ;
        }

    middle = (low + high) / 2;
    if (arr[middle] > arr[high]) ELEM_FLOAT_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high]) ELEM_FLOAT_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low]) ELEM_FLOAT_SWAP(arr[middle], arr[low]) ;

    ELEM_FLOAT_SWAP(arr[middle], arr[low+1]) ;
    ll = low + 1;
    hh = high;

    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_FLOAT_SWAP(arr[ll], arr[hh]) ;
   }

    ELEM_FLOAT_SWAP(arr[low], arr[hh]) ;

    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
 }

float media(float *elements, int N)
{
  int i,j,min;
  float temp;

  //   Order elements (only half of them)
   for (i = 0; i < (N >> 1) + 1; ++i)
   {
      //   Find position of minimum element
      min = i;
      for (j = i + 1; j < N; ++j)
         if (elements[j] < elements[min])
            min = j;
      //   Put found minimum element in its place
      temp = elements[i];
      elements[i] = elements[min];
      elements[min] = temp;
   }
   //   Get result - the middle element
   return elements[N >> 1];
}

void ImProcFunctions::Tile_calc (int tilesize, int overlap, int kall, int imwidth, int imheight, int &numtiles_W, int &numtiles_H, int &tilewidth, int &tileheight, int &tileWskip, int &tileHskip)

{
		if(kall==2) {

		if (imwidth<tilesize) {
			numtiles_W = 1;
			tileWskip = imwidth;
			tilewidth = imwidth;
		} else {
			numtiles_W = ceil(((float)(imwidth))/(tilesize-overlap));
			tilewidth  = ceil(((float)(imwidth))/(numtiles_W))+overlap;
			tilewidth += (tilewidth&1);
			tileWskip = tilewidth-overlap;
		}
		if (imheight<tilesize) {
			numtiles_H = 1;
			tileHskip = imheight;
			tileheight = imheight;
		} else {
			numtiles_H = ceil(((float)(imheight))/(tilesize-overlap));
			tileheight = ceil(((float)(imheight))/(numtiles_H))+overlap;
			tileheight += (tileheight&1);
			tileHskip = tileheight-overlap;
		}
		}
		if(kall==0) {
			numtiles_W = 1;
			tileWskip = imwidth;
			tilewidth = imwidth;
			numtiles_H = 1;
			tileHskip = imheight;
			tileheight = imheight;	
		}
		
	//	printf("Nw=%d NH=%d tileW=%d tileH=%d\n",numtiles_W,numtiles_H,tileWskip,tileHskip);
}

int denoiseNestedLevels = 1;


SSEFUNCTION void ImProcFunctions::RGB_denoise(int kall, Imagefloat * src, Imagefloat * dst,Imagefloat * calclum, float * ch_M, float *max_r, float *max_b, bool isRAW, const procparams::DirPyrDenoiseParams & dnparams, const double expcomp, const NoiseCurve & noiseLCurve, const NoiseCurve & noiseCCurve, float &chaut, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &nresi, float &highresi)
{
//#ifdef _DEBUG
	MyTime t1e,t2e;
	t1e.set();
//#endif

	if (dnparams.luma==0 && dnparams.chroma==0  && !dnparams.median && !noiseLCurve && !noiseCCurve) {
		//nothing to do; copy src to dst or do nothing in case src == dst
		if(src != dst)
			memcpy(dst->data,src->data,dst->width*dst->height*3*sizeof(float));
		if(calclum) {
			delete calclum;
			calclum = NULL;
		}
			
		return;
	}
	static MyMutex FftwMutex;
	MyMutex::MyLock lock(FftwMutex);
	int hei,wid;
	float** lumcalc;
	float** acalc;
	float** bcalc;
	if(noiseLCurve || noiseCCurve) {
		hei=calclum->height;
		wid=calclum->width;
		TMatrix wprofi = iccStore->workingSpaceMatrix (params->icm.working);

		double wpi[3][3] = {
			{wprofi[0][0],wprofi[0][1],wprofi[0][2]},
			{wprofi[1][0],wprofi[1][1],wprofi[1][2]},
			{wprofi[2][0],wprofi[2][1],wprofi[2][2]}
			};
		lumcalc = new float*[(hei+1)/2];
			for (int i=0; i<(hei+1)/2; i++)
				lumcalc[i] = new float[(wid+1)/2];
		acalc = new float*[(hei+1)/2];
			for (int i=0; i<(hei+1)/2; i++)
				acalc[i] = new float[(wid+1)/2];
		bcalc = new float*[(hei+1)/2];
			for (int i=0; i<(hei+1)/2; i++)
				bcalc[i] = new float[(wid+1)/2];
#ifdef _OPENMP
#pragma omp parallel for
#endif		
		for(int ii=0;ii<hei;ii+=2){
			for(int jj=0;jj<wid;jj+=2){
				float LLum,AAum,BBum;
			
				float RL = calclum->r(ii,jj);
				float GL = calclum->g(ii,jj);
				float BL = calclum->b(ii,jj);
							// determine luminance for noisecurve
				float XL,YL,ZL;
				Color::rgbxyz(RL,GL,BL,XL,YL,ZL,wpi);
				Color::XYZ2Lab(XL, YL, ZL, LLum, AAum, BBum);							
				lumcalc[ii>>1][jj>>1]=LLum;
				acalc[ii>>1][jj>>1]=AAum;
				bcalc[ii>>1][jj>>1]=BBum;
			}
		}
	delete calclum;
	calclum = NULL;
}		

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	/*if (plistener) {
	 plistener->setProgressStr ("Denoise...");
	 plistener->setProgress (0.0);
	 }*/

//		volatile double progress = 0.0;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	const short int imheight=src->height, imwidth=src->width;
	Qhigh=1.0f;
	if(dnparams.smethod=="shalbi")
		Qhigh=1.f/(float) settings->nrhigh;
	if (dnparams.luma!=0 || dnparams.chroma!=0 || dnparams.methodmed=="Lab" || dnparams.methodmed=="Lonly" ) {
		const bool perf = (dnparams.dmethod=="RGB");
		// gamma transform for input data
		float gam = dnparams.gamma;
		float gamthresh = 0.001f;
		if(!isRAW) {//reduce gamma under 1 for Lab mode ==> TIF and JPG
			if(gam <1.9f)
				gam=1.f - (1.9f-gam)/3.f;//minimum gamma 0.7
			else if (gam >= 1.9f && gam <= 3.f)
				gam=(1.4f/1.1f)*gam - 1.41818f;
		}
		float gamslope = exp(log((double)gamthresh)/gam)/gamthresh;

		LUTf gamcurve(65536,LUT_CLIP_BELOW);
		if(perf) {
			for (int i=0; i<65536; i++) {
				gamcurve[i] = (Color::gamma((double)i/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)) * 32768.0f;
			}
		} else {
			for (int i=0; i<65536; i++) {
				gamcurve[i] = (Color::gamman((double)i/65535.0,gam)) * 32768.0f;
			}
		}

		// inverse gamma transform for output data
		float igam = 1.f/gam;
		float igamthresh = gamthresh*gamslope;
		float igamslope = 1.f/gamslope;

		LUTf igamcurve(65536,LUT_CLIP_BELOW);
		if(perf) {
			for (int i=0; i<65536; i++) {
				igamcurve[i] = (Color::gamma((float)i/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
			}
		} else {
			for (int i=0; i<65536; i++) {
				igamcurve[i] = (Color::gamman((float)i/32768.0f,igam) * 65535.0f);
			}
		}

		const float gain = pow (2.0f, float(expcomp));
		float incr=1.f;
		float noisevar_Ldetail = SQR((float)(SQR(100.-dnparams.Ldetail) + 50.*(100.-dnparams.Ldetail)) * TS * 0.5f * incr);
		bool enhance_denoise = dnparams.enhance;
		int gamlab = settings->denoiselabgamma;//gamma lab essentialy for Luminance detail
		if(gamlab > 2)
			gamlab=2;
		if(settings->verbose)
			printf("Denoise Lab=%i\n",gamlab);

		array2D<float> tilemask_in(TS,TS);
		array2D<float> tilemask_out(TS,TS);

		const int border = MAX(2,TS/16);
		

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<TS; i++) {
		float i1 = abs((i>TS/2 ? i-TS+1 : i));
		float vmask = (i1<border ? SQR(sin((M_PI*i1)/(2*border))) : 1.0f);
		float vmask2 = (i1<2*border ? SQR(sin((M_PI*i1)/(2*border))) : 1.0f);
		for (int j=0; j<TS; j++) {
			float j1 = abs((j>TS/2 ? j-TS+1 : j));
			tilemask_in[i][j] = (vmask * (j1<border ? SQR(sin((M_PI*j1)/(2*border))) : 1.0f)) + epsilon;
			tilemask_out[i][j] = (vmask2 * (j1<2*border ? SQR(sin((M_PI*j1)/(2*border))) : 1.0f)) + epsilon;

		}
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// begin tile processing of image

	//output buffer
	Imagefloat * dsttmp = new Imagefloat(imwidth,imheight);
	for (int n=0; n<3*imwidth*imheight; n++) dsttmp->data[n] = 0;

int tilesize;
int overlap;
if(settings->leveldnti ==0) {
	tilesize = 1024;
	overlap = 128;
}
if(settings->leveldnti ==1) {
	tilesize = 768;
	overlap = 96;
}

	int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;
	
	Tile_calc (tilesize, overlap, 2, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
	//now we have tile dimensions, overlaps
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// According to FFTW-Doc 'it is safe to execute the same plan in parallel by multiple threads', so we now create 4 plans
	// outside the parallel region and use them inside the parallel region.

	// calculate max size of numblox_W.
	int max_numblox_W = ceil(((float)(MIN(imwidth,tilewidth)))/(offset))+2*blkrad;
	// calculate min size of numblox_W.
	int min_numblox_W = ceil(((float)((MIN(imwidth,((numtiles_W - 1) * tileWskip) + tilewidth) ) - ((numtiles_W - 1) * tileWskip)))/(offset))+2*blkrad;

	// these are needed only for creation of the plans and will be freed before entering the parallel loop
	float * Lbloxtmp;
	float * fLbloxtmp;
	Lbloxtmp  = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof (float));
	fLbloxtmp = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof (float));

	int nfwd[2]={TS,TS};

	//for DCT:
	const fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
	const fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};

	fftwf_plan plan_forward_blox[2];
	fftwf_plan plan_backward_blox[2];

	// Creating the plans with FFTW_MEASURE instead of FFTW_ESTIMATE speeds up the execute a bit
	plan_forward_blox[0]  = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, Lbloxtmp, NULL, 1, TS*TS, fLbloxtmp, NULL, 1, TS*TS, fwdkind, FFTW_MEASURE );
	plan_backward_blox[0] = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, fLbloxtmp, NULL, 1, TS*TS, Lbloxtmp, NULL, 1, TS*TS, bwdkind, FFTW_MEASURE );
	plan_forward_blox[1]  = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, Lbloxtmp, NULL, 1, TS*TS, fLbloxtmp, NULL, 1, TS*TS, fwdkind, FFTW_MEASURE );
	plan_backward_blox[1] = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, fLbloxtmp, NULL, 1, TS*TS, Lbloxtmp, NULL, 1, TS*TS, bwdkind, FFTW_MEASURE );
	fftwf_free ( Lbloxtmp );
	fftwf_free ( fLbloxtmp );

	int numthreads = 1;

#ifdef _OPENMP
	// Calculate number of tiles. If less than omp_get_max_threads(), then limit num_threads to number of tiles
	int numtiles = numtiles_W * numtiles_H;
	numthreads = MIN(numtiles,omp_get_max_threads());
	if(options.rgbDenoiseThreadLimit > 0)
		numthreads = MIN(numthreads,options.rgbDenoiseThreadLimit);
	denoiseNestedLevels = omp_get_max_threads() / numthreads;
	bool oldNested = omp_get_nested();
	if(denoiseNestedLevels < 2)
		denoiseNestedLevels = 1;
	else
		omp_set_nested(true);
	if(options.rgbDenoiseThreadLimit > 0)
		while(denoiseNestedLevels*numthreads > options.rgbDenoiseThreadLimit)
			denoiseNestedLevels--;
	if(settings->verbose)
		printf("RGB_denoise uses %d main thread(s) and up to %d nested thread(s) for each main thread\n",numthreads,denoiseNestedLevels);
#endif

	float *LbloxArray[denoiseNestedLevels*numthreads];
	float *fLbloxArray[denoiseNestedLevels*numthreads];

	for(int i=0;i<denoiseNestedLevels*numthreads;i++) {
		LbloxArray[i]  = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof(float));
		fLbloxArray[i] = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof(float));
	}

	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	//inverse matrix user select
	const float wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};

	TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);

	const float wp[3][3] = {
		{wprof[0][0],wprof[0][1],wprof[0][2]},
		{wprof[1][0],wprof[1][1],wprof[1][2]},
		{wprof[2][0],wprof[2][1],wprof[2][2]}
	};


#ifdef _OPENMP
#pragma omp parallel num_threads(numthreads)
#endif
	{
	float resid=0.f;
	float nbresid=0;
	float maxredresid=0.f;
	float maxblueresid=0.f;
	float residred=0.f;
	float residblue=0.f;
	
	int pos;
	float** noisevarlum = new float*[(tileheight+1)/2];
	for (int i=0; i<(tileheight+1)/2; i++)
		noisevarlum[i] = new float[(tilewidth+1)/2];
	float** noisevarchrom = new float*[(tileheight+1)/2];
	for (int i=0; i<(tileheight+1)/2; i++)
		noisevarchrom[i] = new float[(tilewidth+1)/2];


#ifdef _OPENMP
#pragma omp for schedule(dynamic) collapse(2)
#endif

	for (int tiletop=0; tiletop<imheight; tiletop+=tileHskip) {
		for (int tileleft=0; tileleft<imwidth ; tileleft+=tileWskip) {
			//printf("titop=%d tileft=%d\n",tiletop/tileHskip, tileleft/tileWskip);
			pos= (tiletop/tileHskip)*numtiles_W + tileleft/tileWskip ;
	//	    printf("posw=%d\n", pos);
	/*		if(kall==0) {//with Dcrop
				//int trafx, int trafy, int trafw, int trafh, // trafx trafy = begin window trafw trafh = width and heigh window
				//int widIm, int heiIm : width and height full image
				int centerwCrop=trafx+trafw/2;//center crop window
				int centerhCrop=trafy+trafh/2;
				int poscalculate; 
		//	for(int wcr=0;wcr< ;wcr++) {
		//	for(int hcr=0;hcr<numtiles_H;hcr++) {
		//	    int posca= hcr*numtiles_W + wcr;
		
			}
			*/
			int tileright = MIN(imwidth,tileleft+tilewidth);
			int tilebottom = MIN(imheight,tiletop+tileheight);
			int width  = tileright-tileleft;
			int height = tilebottom-tiletop;
			bool ponder=false;
			float ponderCC=1.f;
			if(settings->leveldnautsimpl==1 && params->dirpyrDenoise.Cmethod=="PON") {ponder=true;ponderCC=0.5f;}
			if(settings->leveldnautsimpl==1 && params->dirpyrDenoise.Cmethod=="PRE") {ponderCC=0.5f;}
			if(settings->leveldnautsimpl==0 && params->dirpyrDenoise.Cmethod=="PREV") {ponderCC=0.5f;}
			
			float realred2, realred, realblue, realblue2;
			float interm_med =(float) dnparams.chroma/10.0;
			float intermred, intermblue;
			if(dnparams.redchro > 0.) intermred=(dnparams.redchro/10.); else intermred= (float) dnparams.redchro/7.0;//increase slower than linear for more sensit
			if(dnparams.bluechro > 0.) intermblue=(dnparams.bluechro/10.); else intermblue= (float) dnparams.bluechro/7.0;//increase slower than linear for more sensit
			if(ponder && kall==2){interm_med=ch_M[pos]/10.f; intermred=max_r[pos]/10.f;intermblue=max_b[pos]/10.f;}							
			if(ponder && kall==0){interm_med=0.01f; intermred=0.f;intermblue=0.f;}							
			realred = interm_med + intermred; if (realred < 0.f) realred=0.001f;
			realblue = interm_med + intermblue; if (realblue < 0.f) realblue=0.001f;
			const float noisevarab_r = SQR(realred);
			const float noisevarab_b = SQR(realblue);
		//	printf("Ch=%f red=%f blu=%f\n",interm_med, intermred,intermblue );

			//input L channel
			array2D<float> Lin(width,height);
			//wavelet denoised image
			LabImage * labdn = new LabImage(width,height);
			float* mad_LL = new float [height*width];

			float* mad_aa = new float [height*width];
			float* mad_bb = new float [height*width];

			//residual between input and denoised L channel
			array2D<float> Ldetail(width,height,ARRAY2D_CLEAR_DATA);
			//pixel weight
			array2D<float> totwt(width,height,ARRAY2D_CLEAR_DATA);//weight for combining DCT blocks
			// init luma noisevarL
			float noiseluma=(float) dnparams.luma;
			if(perf && noiseLCurve)
				noiseluma += 1.f;
			float noisevarL	 = (float) (SQR((noiseluma/125.0)*(1.+ noiseluma/25.0)));
		//	printf("nova=%f\n",noisevarL);
			//TODO: implement using AlignedBufferMP
			//fill tile from image; convert RGB to "luma/chroma"
			if (isRAW) {//image is raw; use channel differences for chroma channels

				if(!perf){//lab mode
						//modification Jacques feb 2013 and july 2014
#ifdef _OPENMP
#pragma omp parallel for num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
				for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
					int i1 = i - tiletop;
					for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
						int j1 = j - tileleft;
						float R_ = gain*src->r(i,j);
						float G_ = gain*src->g(i,j);
						float B_ = gain*src->b(i,j);

						//modify arbitrary data for Lab..I have test : nothing, gamma 2.6 11 - gamma 4 5 - gamma 5.5 10
						//we can put other as gamma g=2.6 slope=11, etc.
						// but noting to do with real gamma !!!: it's only for data Lab # data RGB
						//finally I opted fot gamma55 and with options we can change
						if (gamlab == 0) {// options 12/2013
							R_ = Color::igammatab_26_11[R_];
							G_ = Color::igammatab_26_11[G_];
							B_ = Color::igammatab_26_11[B_];
						}
						else if (gamlab == 1) {
						//other new gamma 4 5 
							R_ = Color::igammatab_4[R_];
							G_ = Color::igammatab_4[G_];
							B_ = Color::igammatab_4[B_];
						}
						else if (gamlab == 2) {							
						//new gamma 5.5 10 better for detail luminance..it is a compromise...which depends on the image (distribution BL, ML, HL ...)
							R_ = Color::igammatab_55[R_];
							G_ = Color::igammatab_55[G_];
							B_ = Color::igammatab_55[B_];
						}
						//apply gamma noise	standard (slider)
						R_ = R_<65535.0f ? gamcurve[R_] : (Color::gammanf(R_/65535.f, gam)*32768.0f);
						G_ = G_<65535.0f ? gamcurve[G_] : (Color::gammanf(G_/65535.f, gam)*32768.0f);
						B_ = B_<65535.0f ? gamcurve[B_] : (Color::gammanf(B_/65535.f, gam)*32768.0f);
						
						//true conversion xyz=>Lab
						float L,a,b;
						float X,Y,Z;
						Color::rgbxyz(R_,G_,B_,X,Y,Z,wp);

						//convert to Lab
						Color::XYZ2Lab(X, Y, Z, L, a, b);

						labdn->L[i1][j1] = L;
						labdn->a[i1][j1] = a;
						labdn->b[i1][j1] = b;
					
						Lin[i1][j1] = L;

						if(((i1|j1)&1) == 0) {
							if(noiseLCurve) {
								float kN = lumcalc[i>>1][j>>1]; //with no gamma and take into account working profile
								float epsi = 0.01f;
								if(kN<2.f)
									kN = 2.f;//avoid divided by zero
								if(kN>32768.f)
									kN = 32768.f;	// not strictly necessary							
								float kinterm = epsi + noiseLCurve[xdivf(kN,15)*500.f];
								float ki = kinterm*100.f;
								ki += noiseluma;
								noisevarlum[i1>>1][j1>>1] = SQR((ki/125.f)*(1.f+ki/25.f));
							}
							if(noiseCCurve) {
								float aN = acalc[i>>1][j>>1];
								float bN = bcalc[i>>1][j>>1];
								float cN = sqrtf(SQR(aN)+SQR(bN));
								if(cN < 100.f)
									cN = 100.f;//avoid divided by zero
								float Cinterm = 1.f + ponderCC*4.f*noiseCCurve[cN/60.f];//C=f(C)
								noisevarchrom[i1>>1][j1>>1] = max(noisevarab_b,noisevarab_r)*SQR(Cinterm);
							}
						}
						//end chroma

					}
				}
				}
				else {//RGB mode
#ifdef _OPENMP
#pragma omp parallel for num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
				for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
					int i1 = i - tiletop;
					for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
						int j1 = j - tileleft;

						float X = gain*src->r(i,j);
						float Y = gain*src->g(i,j);
						float Z = gain*src->b(i,j);
						//conversion colorspace to determine luminance with no gamma
						X = X<65535.0f ? gamcurve[X] : (Color::gamma((double)X/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						Y = Y<65535.0f ? gamcurve[Y] : (Color::gamma((double)Y/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						Z = Z<65535.0f ? gamcurve[Z] : (Color::gamma((double)Z/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						if(((i1|j1)&1) == 0) {
							if(noiseLCurve) {
							//	float noiseluma=(float) dnparams.luma;	
								float kN = lumcalc[i>>1][j>>1];
								float epsi = 0.01f;
								if (kN<2.f) 
									kN = 2.f;
								if (kN>32768.f)
									kN = 32768.f;
								float kinterm = epsi + noiseLCurve[xdivf(kN,15)*500.f];
								float ki = kinterm*100.f;
								ki += noiseluma;
								noisevarlum[i1>>1][j1>>1] = SQR((ki/125.f)*(1.f+ki/25.f));
							}
							if(noiseCCurve) {
								float aN = acalc[i>>1][j>>1];
								float bN = bcalc[i>>1][j>>1];
								float cN = sqrtf(SQR(aN)+SQR(bN));
								if(cN < 100.f)
									cN = 100.f;//avoid divided by zero
								float Cinterm = 1.f + ponderCC*4.f*noiseCCurve[cN/60.f];
								noisevarchrom[i1>>1][j1>>1] = max(noisevarab_b,noisevarab_r)*SQR(Cinterm);
							}
						}
						//end chroma					
						labdn->L[i1][j1] = Y;
						labdn->a[i1][j1] = (X-Y);
						labdn->b[i1][j1] = (Y-Z);

						Lin[i1][j1] = Y;
					}
				}

				}
				
			} else {//image is not raw; use Lab parametrization
#ifdef _OPENMP
#pragma omp parallel for num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
				for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
					int i1 = i - tiletop;
					for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
						int j1 = j - tileleft;
						float L,a,b;
						float rLum=src->r(i,j) ;//for luminance denoise curve
						float gLum=src->g(i,j) ;
						float bLum=src->b(i,j) ;
						
						//use gamma sRGB, not good if TIF (JPG) Output profil not with gamma sRGB  (eg : gamma =1.0, or 1.8...)
						//very difficult to solve !
						// solution ==> save TIF with gamma sRGB and re open
						float rtmp = Color::igammatab_srgb[ src->r(i,j) ];
						float gtmp = Color::igammatab_srgb[ src->g(i,j) ];
						float btmp = Color::igammatab_srgb[ src->b(i,j) ];
						//modification Jacques feb 2013
						// gamma slider different from raw
						rtmp = rtmp<65535.0f ? gamcurve[rtmp] : (Color::gamman((double)rtmp/65535.0, gam)*32768.0f);
						gtmp = gtmp<65535.0f ? gamcurve[gtmp] : (Color::gamman((double)gtmp/65535.0, gam)*32768.0f);
						btmp = btmp<65535.0f ? gamcurve[btmp] : (Color::gamman((double)btmp/65535.0, gam)*32768.0f);

						float X,Y,Z;
						Color::rgbxyz(rtmp,gtmp,btmp,X,Y,Z,wp);

						//convert Lab
						Color::XYZ2Lab(X, Y, Z, L, a, b);
						float Llum,alum,blum;
						if(((i1|j1)&1) == 0) {
							if(noiseLCurve || noiseCCurve) {
								float XL,YL,ZL;
								Color::rgbxyz(rLum,gLum,bLum,XL,YL,ZL,wp);
								Color::XYZ2Lab(XL, YL, ZL, Llum, alum, blum);
							}
					
							if(noiseLCurve) {
								float kN = Llum;
								float epsi=0.01f;
								if(kN<2.f) kN=2.f;
								if(kN>32768.f) kN=32768.f;								
								float kinterm=epsi + noiseLCurve[xdivf(kN,15)*500.f];
								float ki=kinterm*100.f;
								ki+=noiseluma;
								noiseluma += 1.f;
								noisevarlum[i1>>1][j1>>1]=SQR((ki/125.f)*(1.f+ki/25.f));
							}	
							if(noiseCCurve) {
								float aN=alum;
								float bN=blum;
								float cN=sqrtf(SQR(aN)+SQR(bN));
								if(cN < 100.f)
									cN=100.f;//avoid divided by zero ???
								float Cinterm=1.f + ponderCC*4.f*noiseCCurve[cN/60.f];
								noisevarchrom[i1>>1][j1>>1]=max(noisevarab_b,noisevarab_r)*SQR(Cinterm);
							}
						}
						labdn->L[i1][j1] = L;
						labdn->a[i1][j1] = a;
						labdn->b[i1][j1] = b;

						Lin[i1][j1] = L;
					}
				}
			}

			//initial impulse denoise, removed in Issue 2557
//				if (dnparams.luma>0.01) {
//					impulse_nr (labdn, float(MIN(50.0,dnparams.luma))/20.0f);
//				}

			int datalen = labdn->W * labdn->H;

			//now perform basic wavelet denoise
			//last two arguments of wavelet decomposition are max number of wavelet decomposition levels;
			//and whether to subsample the image after wavelet filtering.  Subsampling is coded as
			//binary 1 or 0 for each level, eg subsampling = 0 means no subsampling, 1 means subsample
			//the first level only, 7 means subsample the first three levels, etc.
			float interm_medT = (float) dnparams.chroma/10.0;
			bool execwavelet = true;
			bool autoch = false;
			if(noisevarL < 0.000007f && interm_medT < 0.05f && dnparams.median && (dnparams.methodmed=="Lab" || dnparams.methodmed=="Lonly"))
				execwavelet=false;//do not exec wavelet if sliders luminance and chroma are very small and median need
			//we considered user don't want wavelet
			if(settings->leveldnautsimpl==1 && dnparams.Cmethod!="MAN")
				execwavelet=true;
			if(settings->leveldnautsimpl==0 && dnparams.C2method!="MANU")
				execwavelet=true;
			if(settings->leveldnautsimpl==1 && (dnparams.Cmethod=="AUT" || dnparams.Cmethod=="PRE"))
				autoch=true;
			if(settings->leveldnautsimpl==0 && dnparams.C2method=="AUTO" || dnparams.C2method=="PREV")
				autoch=true;
			
			if(execwavelet) {//gain time if user choose only median  sliders L <=1  slider chrom master < 1
			{ // enclosing this code in a block frees about 120 MB before allocating 20 MB after this block (measured with D700 NEF)
			wavelet_decomposition* Ldecomp;
			wavelet_decomposition* adecomp;
			wavelet_decomposition* bdecomp;
			int schoice=0;//shrink method
			if(dnparams.smethod=="shal") schoice=0;
			else if(dnparams.smethod=="shalbi") schoice=2;

			int levwav=5;
			float maxreal = max(realred, realblue);
			//increase the level of wavelet if user increase much or very much sliders
			if( maxreal < 8.f) levwav=5;
			else if( maxreal < 10.f)levwav=6;
			else if( maxreal < 15.f)levwav=7;
			else levwav=8;//maximum ==> I have increase Maxlevel in cplx_wavelet_dec.h from 8 to 9
			if(schoice==2) levwav+=settings->nrwavlevel;//increase level for enhanced mode
			if(levwav>8) levwav=8;
			
			//	if (settings->verbose) printf("levwavelet=%i  noisevarA=%f noisevarB=%f \n",levwav, noisevarab_r, noisevarab_b );
#ifdef _OPENMP
#pragma omp parallel sections num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
{
#ifdef _OPENMP
#pragma omp section
#endif
{
			Ldecomp = new wavelet_decomposition (labdn->data, labdn->W, labdn->H, levwav/*maxlevels*/, 1/*subsampling*/ );
}
#ifdef _OPENMP
#pragma omp section
#endif
{
			adecomp = new wavelet_decomposition (labdn->data+datalen, labdn->W, labdn->H,levwav, 1 );
}
#ifdef _OPENMP
#pragma omp section
#endif
{
			bdecomp = new wavelet_decomposition (labdn->data+2*datalen, labdn->W, labdn->H, levwav, 1 );
}
}
			if(schoice==0) WaveletDenoiseAll(*Ldecomp, *adecomp, *bdecomp, noisevarL, noisevarlum, noisevarchrom, width, height, mad_LL, mad_aa, mad_bb, noisevarab_r, noisevarab_b,labdn, noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, schoice, autoch, perf );//enhance mode
			if(schoice==2) {
				WaveletDenoiseAll_BiShrink(*Ldecomp, *adecomp, *bdecomp, noisevarL, noisevarlum, noisevarchrom, width, height, mad_LL,  mad_aa, mad_bb, noisevarab_r, noisevarab_b,labdn, noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, schoice, autoch, perf );//enhance mode
				WaveletDenoiseAll(*Ldecomp, *adecomp, *bdecomp, noisevarL, noisevarlum, noisevarchrom, width, height, mad_LL,  mad_aa, mad_bb, noisevarab_r, noisevarab_b,labdn, noiseLCurve, noiseCCurve, chaut ,redaut, blueaut, maxredaut, maxblueaut, schoice, autoch, perf );
				}
			float chresid,chmaxredresid,chmaxblueresid,chresidred, chresidblue;	
			//kall=0 call by Dcrop
			resid=0.f;residred=0.f;residblue=0.f;
			if(kall==0) Noise_residual(*Ldecomp, *adecomp, *bdecomp,  width, height, chresid, chmaxredresid,chmaxblueresid, chresidred, chresidblue, resid, residblue, residred, maxredresid, maxblueresid, nbresid);
			//printf("NoiRESID=%3.1f maxR=%3.1f maxB=%3.1f red=%3.1f  blue=%3.1f\n",chresid, chmaxredresid,chmaxblueresid, chresidred, chresidblue);	
			nresi=chresid;
			highresi=chresid + 0.66f*(max(chmaxredresid,chmaxblueresid) - chresid);//evaluate sigma
#ifdef _OPENMP
#pragma omp parallel sections num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
{
#ifdef _OPENMP
#pragma omp section
#endif
{
			Ldecomp->reconstruct(labdn->data);
			delete Ldecomp;
}
#ifdef _OPENMP
#pragma omp section
#endif
{
			adecomp->reconstruct(labdn->data+datalen);
			delete adecomp;
}
#ifdef _OPENMP
#pragma omp section
#endif
{
			bdecomp->reconstruct(labdn->data+2*datalen);
			delete bdecomp;
}
}
			}
			}
			//TODO: at this point wavelet coefficients storage can be freed
			//Issue 1680: Done now

			//second impulse denoise, removed in Issue 2557
//				if (dnparams.luma>0.01) {
//					impulse_nr (labdn, MIN(50.0f,(float)dnparams.luma)/20.0f);
//				}
			//PF_correct_RT(dst, dst, defringe.radius, defringe.threshold);

			
			int metchoice=0;
			if(dnparams.methodmed=="Lonly") metchoice=1;
			else if(dnparams.methodmed=="Lab") metchoice=2;
			else if(dnparams.methodmed=="ab") metchoice=3;
			else if(dnparams.methodmed=="Lpab") metchoice=4;
			
			//median on Luminance Lab only
		if( (metchoice==1 || metchoice==2 || metchoice==3 || metchoice==4)   && dnparams.median) {
		//printf("Lab et Lonly \n");
			float** tmL;
			int wid=labdn->W;
			int hei=labdn->H;
			tmL = new float*[hei];
			for (int i=0; i<hei; i++)
				tmL[i] = new float[wid];
			int methmedL=0;
			int methmedAB=0;
			int borderL = 1;
			if(dnparams.medmethod=="soft")
				{if(metchoice!=4) methmedL=methmedAB=0;
				else {methmedL=0;methmedAB=0;}
				}
			else if(dnparams.medmethod=="33")
				{if(metchoice!=4) methmedL=methmedAB=1;
				else {methmedL=0;methmedAB=1;borderL = 1;}
				}

			//	methmedL=methmedAB=1;
			else if(dnparams.medmethod=="55") {
				{if(metchoice!=4) methmedL=methmedAB=3;
				else {methmedL=1;methmedAB=3;}
				}

			//	methmedL =methmedAB= 3;
				borderL = 2;
			}
			else if(dnparams.medmethod=="55soft") {
				{if(metchoice!=4) methmedL=methmedAB=2;
				else {methmedL=0;methmedAB=2;}
				}

			//	methmedL = methmedAB= 2;
				borderL = 2;
			}
			else if(dnparams.medmethod=="77") {
				{if(metchoice!=4) methmedL=methmedAB=4;
				else {methmedL=1;methmedAB=4;}
				}

			//	methmedL = methmedAB=4;
				borderL = 3;
			}
			else if(dnparams.medmethod=="99") {
				{if(metchoice!=4) methmedL=methmedAB=5;
				else {methmedL=2;methmedAB=5;}
				}

			//	methmedL = methmedAB = 5;
				borderL = 4;
			}
		for(int iteration=1;iteration<=dnparams.passes;iteration++){
		  //printf("pas=%i\n",iteration);



			if (metchoice==1 || metchoice==2 || metchoice==4)
			{ /*printf("LONLY methmedL=%d\n", methmedL);*/

				if(methmedL < 2) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
					for (int i=1; i<hei-1; i++) {
						float pp[9],results[5],temp;
						if(methmedL == 0) 
							for (int j=1; j<wid-1; j++) {
								med2(labdn->L[i][j] ,labdn->L[i-1][j], labdn->L[i+1][j] ,labdn->L[i][j+1],labdn->L[i][j-1], tmL[i][j]);//3x3 soft
							}
						else
							for (int j=1; j<wid-1; j++) {//hybrid median filter ==> impluse denoise !I let here to get trace
								/*	pp[0] = labdn->L[i][j-1] ;
									pp[1] = labdn->L[i-1][j] ;
									pp[2] = labdn->L[i][j] ;
									pp[3] = labdn->L[i+1][j] ;
									pp[4] = labdn->L[i][j+1] ;
								//   Get median
									results[0] = media(pp, 5);
								//   Pick up x-window elements
									pp[0] = labdn->L[i-1][j-1] ;
									pp[1] = labdn->L[i+1][j-1] ;
									pp[2] = labdn->L[i][j] ;
									pp[3] = labdn->L[i-1][j+1] ;
									pp[4] = labdn->L[i+1][j+1] ;
									  //   Get median
									results[1] = media(pp, 5);
									//   Pick up leading element
									results[2] = labdn->L[i][j] ;;
									//   Get result
									tmL[i][j] = media(results, 3);
							*/
							med3(labdn->L[i][j] ,labdn->L[i-1][j], labdn->L[i+1][j] ,labdn->L[i][j+1],labdn->L[i][j-1], labdn->L[i-1][j-1],labdn->L[i-1][j+1],labdn->L[i+1][j-1],labdn->L[i+1][j+1],tmL[i][j]);//3x3 soft
							}
					}
				}
				else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
					for (int i=borderL; i<hei-borderL; i++) {
						float pp[81],temp;
					if(methmedL == 4)						
						for (int j=borderL; j<wid-borderL; j++) {
								int kk=0;
								for (int ii=-3;ii<=3;ii++) {
									for (int jj=-3;jj<=3;jj++) {
										kk++;
										pp[kk]=labdn->L[i+ii][j+jj];
									}
								}
							fq_sort2(pp,49);
							tmL[i][j]=pp[24];//7x7
						}
					else if(methmedL == 5)						
						for (int j=borderL; j<wid-borderL; j++) {
								int kk=0;
								for (int ii=-4;ii<=4;ii++) {
									for (int jj=-4;jj<=4;jj++) {
										kk++;
										pp[kk]=labdn->L[i+ii][j+jj];
									}
								}
							fq_sort2(pp,81);
							tmL[i][j]=pp[40];//9x9
						}
						
					else if(methmedL == 3)	
						for (int j=2; j<wid-2; j++) {
							med5(labdn->L[i][j],labdn->L[i-1][j],labdn->L[i+1][j],labdn->L[i][j+1],labdn->L[i][j-1],labdn->L[i-1][j-1],labdn->L[i-1][j+1], labdn->L[i+1][j-1],labdn->L[i+1][j+1],
							labdn->L[i-2][j],labdn->L[i+2][j],labdn->L[i][j+2],labdn->L[i][j-2],labdn->L[i-2][j-2],labdn->L[i-2][j+2],labdn->L[i+2][j-2],labdn->L[i+2][j+2],	
							labdn->L[i-2][j+1],labdn->L[i+2][j+1],labdn->L[i-1][j+2],labdn->L[i-1][j-2],labdn->L[i-2][j-1],labdn->L[i+2][j-1],labdn->L[i+1][j+2],labdn->L[i+1][j-2],	
							tmL[i][j]);//5x5
						}
					else
						for (int j=2; j<wid-2; j++) {
							pp[0]=labdn->L[i][j];pp[1]=labdn->L[i-1][j]; pp[2]=labdn->L[i+1][j];pp[3]=labdn->L[i][j+1];pp[4]=labdn->L[i][j-1];pp[5]=labdn->L[i-1][j-1];pp[6]=labdn->L[i-1][j+1];
							pp[7]=labdn->L[i+1][j-1];pp[8]=labdn->L[i+1][j+1];pp[9]=labdn->L[i+2][j];pp[10]=labdn->L[i-2][j];pp[11]=labdn->L[i][j+2];pp[12]=labdn->L[i][j-2];
							fq_sort2(pp,13);
							tmL[i][j]=pp[6];//5x5 soft
						}
					}
				}

				
	
				for(int i = borderL; i < hei-borderL; i++ ) {
					for(int j = borderL; j < wid-borderL; j++) {
						labdn->L[i][j] = tmL[i][j];
					}
				}
		 }    
			if(metchoice==2   || metchoice==3 || metchoice==4) {/*printf(" AB methmedL=%d\n", methmedL);*/
				if(methmedAB < 2) {
					for (int i=1; i<hei-1; i++) {
						float pp[9],temp;
						if(methmedAB == 0) 
							for (int j=1; j<wid-1; j++) {
								med2(labdn->a[i][j] ,labdn->a[i-1][j], labdn->a[i+1][j] ,labdn->a[i][j+1],labdn->a[i][j-1], tmL[i][j]);//3x3 soft
							}
						else
							for (int j=1; j<wid-1; j++) {
							med3(labdn->a[i][j] ,labdn->a[i-1][j], labdn->a[i+1][j] ,labdn->a[i][j+1],labdn->a[i][j-1], labdn->a[i-1][j-1],labdn->a[i-1][j+1],labdn->a[i+1][j-1],labdn->a[i+1][j+1],tmL[i][j]);//3x3 soft
							}
					}
				}
				else {
					for (int i=borderL; i<hei-borderL; i++) {
						float pp[81],temp;
					if(methmedAB == 4)						
						for (int j=borderL; j<wid-borderL; j++) {
								int kk=0;
								for (int ii=-3;ii<=3;ii++) {
									for (int jj=-3;jj<=3;jj++) {
										kk++;
										pp[kk]=labdn->a[i+ii][j+jj];
									}
								}
							fq_sort2(pp,49);
							tmL[i][j]=pp[24];//7x7
						}
					else if(methmedAB == 5)						
						for (int j=borderL; j<wid-borderL; j++) {
								int kk=0;
								for (int ii=-4;ii<=4;ii++) {
									for (int jj=-4;jj<=4;jj++) {
										kk++;
										pp[kk]=labdn->a[i+ii][j+jj];
									}
								}
							fq_sort2(pp,81);
							tmL[i][j]=pp[40];//9
						}
						
					else if(methmedAB == 3)
						for (int j=2; j<wid-2; j++) {						
							med5(labdn->a[i][j],labdn->a[i-1][j],labdn->a[i+1][j],labdn->a[i][j+1],labdn->a[i][j-1],labdn->a[i-1][j-1],labdn->a[i-1][j+1], labdn->a[i+1][j-1],labdn->a[i+1][j+1],
							labdn->a[i-2][j],labdn->a[i+2][j],labdn->a[i][j+2],labdn->a[i][j-2],labdn->a[i-2][j-2],labdn->a[i-2][j+2],labdn->a[i+2][j-2],labdn->a[i+2][j+2],	
							labdn->a[i-2][j+1],labdn->a[i+2][j+1],labdn->a[i-1][j+2],labdn->a[i-1][j-2],labdn->a[i-2][j-1],labdn->a[i+2][j-1],labdn->a[i+1][j+2],labdn->a[i+1][j-2],	
							tmL[i][j]);//5x5
							}
					else
						for (int j=2; j<wid-2; j++) {
							pp[0]=labdn->a[i][j];pp[1]=labdn->a[i-1][j]; pp[2]=labdn->a[i+1][j];pp[3]=labdn->a[i][j+1];pp[4]=labdn->a[i][j-1];pp[5]=labdn->a[i-1][j-1];pp[6]=labdn->a[i-1][j+1];
							pp[7]=labdn->a[i+1][j-1];pp[8]=labdn->a[i+1][j+1];pp[9]=labdn->a[i+2][j];pp[10]=labdn->a[i-2][j];pp[11]=labdn->a[i][j+2];pp[12]=labdn->a[i][j-2];
							fq_sort2(pp,13);
							tmL[i][j]=pp[6];//5x5 soft
						}
					/*	for (int j=3; j<wid-3; j++) {
								int kk=0;
								for (int ii=-2;ii<=2;ii++) {
									for (int jj=-2;jj<=2;jj++){kk++;
									 pp[kk]=labdn->a[i+ii][j+jj];
									}
									pp[kk+1]=labdn->a[i-3][j-3];pp[kk+2]=labdn->a[i-3][j+3];pp[kk+3]=labdn->a[i+3][j-3];pp[kk+4]=labdn->a[i+3][j+3];
									pp[kk+5]=labdn->a[i-3][j];pp[kk+6]=labdn->a[i+3][j];pp[kk+7]=labdn->a[i][j-3];pp[kk+8]=labdn->a[i][j+3];
									
								}	
							fq_sort2(pp,33);
							tmL[i][j]=pp[16];//7x7
						}
						*/
					}
				}

	
				for(int i = borderL; i < hei-borderL; i++ ) {
					for(int j = borderL; j < wid-borderL; j++) {
						labdn->a[i][j] = tmL[i][j];
					}
				}
				
				
//b
				if(methmedAB < 2) {
					for (int i=1; i<hei-1; i++) {
						float pp[9],temp;
						if(methmedAB == 0) 
							for (int j=1; j<wid-1; j++) {
								med2(labdn->b[i][j] ,labdn->b[i-1][j], labdn->b[i+1][j] ,labdn->b[i][j+1],labdn->b[i][j-1], tmL[i][j]);//3x3 soft
							}
						else
							for (int j=1; j<wid-1; j++) {
							med3(labdn->b[i][j] ,labdn->b[i-1][j], labdn->b[i+1][j] ,labdn->b[i][j+1],labdn->b[i][j-1], labdn->b[i-1][j-1],labdn->b[i-1][j+1],labdn->b[i+1][j-1],labdn->b[i+1][j+1],tmL[i][j]);//3x3 soft
							}
					}
				}
				else {
					for (int i=borderL; i<hei-borderL; i++) {
						float pp[81],temp;
					if(methmedL == 4)			
						for (int j=borderL; j<wid-borderL; j++) {
								int kk=0;
								for (int ii=-3;ii<=3;ii++) {
									for (int jj=-3;jj<=3;jj++) {
										kk++;
										pp[kk]=labdn->b[i+ii][j+jj];
									}
								}
							fq_sort2(pp,49);
							tmL[i][j]=pp[24];//7x7
						
							
						}
					else if(methmedAB == 5)			
						for (int j=borderL; j<wid-borderL; j++) {
								int kk=0;
								for (int ii=-4;ii<=4;ii++) {
									for (int jj=-4;jj<=4;jj++) {
										kk++;
										pp[kk]=labdn->b[i+ii][j+jj];
									}
								}
							fq_sort2(pp,81);
							tmL[i][j]=pp[40];//9
						}
						
					else if(methmedAB == 3)			
						for (int j=2; j<wid-2; j++) {
							med5(labdn->b[i][j],labdn->b[i-1][j],labdn->b[i+1][j],labdn->b[i][j+1],labdn->b[i][j-1],labdn->b[i-1][j-1],labdn->b[i-1][j+1], labdn->b[i+1][j-1],labdn->b[i+1][j+1],
							labdn->b[i-2][j],labdn->b[i+2][j],labdn->b[i][j+2],labdn->b[i][j-2],labdn->b[i-2][j-2],labdn->b[i-2][j+2],labdn->b[i+2][j-2],labdn->b[i+2][j+2],	
							labdn->b[i-2][j+1],labdn->b[i+2][j+1],labdn->b[i-1][j+2],labdn->b[i-1][j-2],labdn->b[i-2][j-1],labdn->b[i+2][j-1],labdn->b[i+1][j+2],labdn->b[i+1][j-2],	
							tmL[i][j]);//5x5
						
						}
					else
						for (int j=2; j<wid-2; j++) {
						/*for (int j=3; j<wid-3; j++) {
								int kk=0;
								for (int ii=-2;ii<=2;ii++) {
									for (int jj=-2;jj<=2;jj++){kk++;
									 pp[kk]=labdn->b[i+ii][j+jj];
									}
									pp[kk+1]=labdn->b[i-3][j-3];pp[kk+2]=labdn->b[i-3][j+3];pp[kk+3]=labdn->b[i+3][j-3];pp[kk+4]=labdn->b[i+3][j+3];
									pp[kk+5]=labdn->b[i-3][j];pp[kk+6]=labdn->b[i+3][j];pp[kk+7]=labdn->b[i][j-3];pp[kk+8]=labdn->b[i][j+3];
									
								}	
							fq_sort2(pp,33);
							tmL[i][j]=pp[16];//7x7
						*/
						
						
							pp[0]=labdn->b[i][j];pp[1]=labdn->b[i-1][j]; pp[2]=labdn->b[i+1][j];pp[3]=labdn->b[i][j+1];pp[4]=labdn->b[i][j-1];pp[5]=labdn->b[i-1][j-1];pp[6]=labdn->b[i-1][j+1];
							pp[7]=labdn->b[i+1][j-1];pp[8]=labdn->b[i+1][j+1];pp[9]=labdn->b[i+2][j];pp[10]=labdn->b[i-2][j];pp[11]=labdn->b[i][j+2];pp[12]=labdn->b[i][j-2];
							fq_sort2(pp,13);
							tmL[i][j]=pp[6];//5x5 soft
						
						}
					}
				}

				for(int i = borderL; i < hei-borderL; i++ ) {
					for(int j = borderL; j < wid-borderL; j++) {
						labdn->b[i][j] = tmL[i][j];
					}
				}
			}			

		}
			for (int i=0; i<hei; i++)
				delete [] tmL[i];
			delete [] tmL;
		}

			//wavelet denoised L channel
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// now do detail recovery using block DCT to detect
			// patterns missed by wavelet denoise
			// blocks are not the same thing as tiles!



			// calculation for detail recovery blocks
			const int numblox_W = ceil(((float)(width))/(offset))+2*blkrad;
			const int numblox_H = ceil(((float)(height))/(offset))+2*blkrad;

			//const int nrtiles = numblox_W*numblox_H;
			// end of tiling calc
			{

			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Main detail recovery algorithm: Block loop
				//DCT block data storage


#ifdef _OPENMP
			int masterThread = omp_get_thread_num();
#endif
#ifdef _OPENMP
#pragma omp parallel num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
{
#ifdef _OPENMP
			int subThread = masterThread * denoiseNestedLevels + omp_get_thread_num();
#else
			int subThread = 0;
#endif
			float blurbuffer[TS*TS] ALIGNED64;
			float * Lblox = LbloxArray[subThread];
			float * fLblox = fLbloxArray[subThread];
			float pBuf[width + TS + 2*blkrad*offset] ALIGNED16;
			float nbrwt[TS*TS] ALIGNED64;
#ifdef _OPENMP
#pragma omp for
#endif
			for (int vblk=0; vblk<numblox_H; vblk++) {

				int top = (vblk-blkrad)*offset;
				float * datarow = pBuf +blkrad*offset;

				for (int i=0/*, row=top*/; i<TS; i++/*, row++*/) {
					int row = top + i;
					int rr = row;
					if (row<0) {
						rr = MIN(-row,height-1);
					} else if (row>=height) {
						rr = MAX(0,2*height-2-row);
					}

					for (int j=0; j<labdn->W; j++) {
						datarow[j] = (Lin[rr][j]-labdn->L[rr][j]);
					}

					for (int j=-blkrad*offset; j<0; j++) {
						datarow[j] = datarow[MIN(-j,width-1)];
					}
					for (int j=width; j<width+TS+blkrad*offset; j++) {
						datarow[j] = datarow[MAX(0,2*width-2-j)];
					}//now we have a padded data row

					//now fill this row of the blocks with Lab high pass data
					for (int hblk=0; hblk<numblox_W; hblk++) {
						int left = (hblk-blkrad)*offset;
						int indx = (hblk)*TS;//index of block in malloc

						for (int j=0; j<TS; j++) {
							Lblox[(indx + i)*TS+j] = tilemask_in[i][j]*datarow[left+j];// luma data
							if (top+i>=0 && top+i<height && left+j>=0 && left+j<width) {
								totwt[top+i][left+j] += tilemask_in[i][j]*tilemask_out[i][j];
							}
						}
					}
				}//end of filling block row

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				//fftwf_print_plan (plan_forward_blox);
				if(numblox_W == max_numblox_W)
					fftwf_execute_r2r(plan_forward_blox[0],Lblox,fLblox);		// DCT an entire row of tiles
				else
					fftwf_execute_r2r(plan_forward_blox[1],Lblox,fLblox);		// DCT an entire row of tiles
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// now process the vblk row of blocks for noise reduction


				for (int hblk=0; hblk<numblox_W; hblk++) {
					RGBtile_denoise (fLblox, hblk, noisevar_Ldetail, nbrwt, blurbuffer );
				}//end of horizontal block loop

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				//now perform inverse FT of an entire row of blocks
				if(numblox_W == max_numblox_W)
					fftwf_execute_r2r(plan_backward_blox[0],fLblox,Lblox);	//for DCT
				else
					fftwf_execute_r2r(plan_backward_blox[1],fLblox,Lblox);	//for DCT

				int topproc = (vblk-blkrad)*offset;

				//add row of blocks to output image tile
				RGBoutput_tile_row (Lblox, Ldetail, tilemask_out, height, width, topproc );

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			}//end of vertical block loop
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

}
			}
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			for (int i=0; i<height; i++) {
				for (int j=0; j<width; j++) {
					//may want to include masking threshold for large hipass data to preserve edges/detail
					float hpdn = Ldetail[i][j]/totwt[i][j];//note that labdn initially stores the denoised hipass data

					labdn->L[i][j] += hpdn;

				}
			}

			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// transform denoised "Lab" to output RGB

			//calculate mask for feathering output tile overlaps
			float Vmask[height+1] ALIGNED16;
			float Hmask[width+1] ALIGNED16;

			for (int i=0; i<height; i++) {
				Vmask[i] = 1;
			}
			float newGain = 1.f;
			if(isRAW)
				newGain = gain;
			for (int j=0; j<width; j++) {
				Hmask[j] = 1.f/newGain;
			}
				
			for (int i=0; i<overlap; i++) {
				float mask = SQR(xsinf((M_PI*i)/(2*overlap)));
				if (tiletop>0) Vmask[i] = mask;
				if (tilebottom<imheight) Vmask[height-i] = mask;
				if (tileleft>0) Hmask[i] = mask/newGain;
				if (tileright<imwidth) Hmask[width-i] = mask/newGain;
			}
			//convert back to RGB and write to destination array
			if (isRAW) {
			if(!perf) {//Lab mode
					realred /= 100.f;
					realblue /= 100.f;
#ifdef _OPENMP
#pragma omp parallel for num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
				for (int i=tiletop; i<tilebottom; i++){
					int i1 = i-tiletop;
					float X,Y,Z,L,a,b;

					for (int j=tileleft; j<tileright; j++) {
						int j1=j-tileleft;
						//modification Jacques feb 2013
						//true conversion Lab==>xyz
						L = labdn->L[i1][j1];
						a = labdn->a[i1][j1];
						b = labdn->b[i1][j1];
						float c_h=SQR(a)+SQR(b);
						if(c_h>9000000.f){
							a *= 1.f + Qhigh*realred;
							b *= 1.f + Qhigh*realblue;
						}
						//convert XYZ
						Color::Lab2XYZ(L, a, b, X, Y, Z);
						//apply inverse gamma noise
						float r_,g_,b_;
						Color::xyz2rgb(X,Y,Z,r_,g_,b_,wip);
						//inverse gamma standard (slider)
						r_ = r_<32768.f ? igamcurve[r_] : (Color::gammanf(r_/32768.f, igam) * 65535.f);
						g_ = g_<32768.f ? igamcurve[g_] : (Color::gammanf(g_/32768.f, igam) * 65535.f);
						b_ = b_<32768.f ? igamcurve[b_] : (Color::gammanf(b_/32768.f, igam) * 65535.f);
						//readapt arbitrary gamma (inverse from beginning)
						if (gamlab == 0) {
							r_ = Color::gammatab_26_11[r_];
							g_ = Color::gammatab_26_11[g_];
							b_ = Color::gammatab_26_11[b_];
						}
						else if (gamlab == 1) {
							r_ = Color::gammatab_4[r_];
							g_ = Color::gammatab_4[g_];
							b_ = Color::gammatab_4[b_];
						}
						else if (gamlab == 2) {
							r_ = Color::gammatab_55[r_];
							g_ = Color::gammatab_55[g_];
							b_ = Color::gammatab_55[b_];
						}
						float factor = Vmask[i1]*Hmask[j1];

						dsttmp->r(i,j) += factor*r_;
						dsttmp->g(i,j) += factor*g_;
						dsttmp->b(i,j) += factor*b_;

					}
				}
				}
				else {//RGB mode
				for (int i=tiletop; i<tilebottom; i++){
					int i1 = i-tiletop;
					float X,Y,Z;
					for (int j=tileleft; j<tileright; j++) {
						int j1=j-tileleft;
						float c_h=sqrt(SQR(labdn->a[i1][j1])+SQR(labdn->b[i1][j1]));
						if(c_h>3000.f){
							labdn->a[i1][j1]*=1.f + Qhigh*realred/100.f;
							labdn->b[i1][j1]*=1.f + Qhigh*realblue/100.f;
						}
						Y = labdn->L[i1][j1];
						X = (labdn->a[i1][j1]) + Y;
						Z = Y - (labdn->b[i1][j1]);
						

						X = X<32768.0f ? igamcurve[X] : (Color::gamma((float)X/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
						Y = Y<32768.0f ? igamcurve[Y] : (Color::gamma((float)Y/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
						Z = Z<32768.0f ? igamcurve[Z] : (Color::gamma((float)Z/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);

						float factor = Vmask[i1]*Hmask[j1];

						dsttmp->r(i,j) += factor*X;
						dsttmp->g(i,j) += factor*Y;
						dsttmp->b(i,j) += factor*Z;

					}
				}

				}
			} else {
				for (int i=tiletop; i<tilebottom; i++){
					int i1 = i-tiletop;
					float X,Y,Z,L,a,b;
					for (int j=tileleft; j<tileright; j++) {
						int j1=j-tileleft;
						//modification Jacques feb 2013
						L = labdn->L[i1][j1];
						a = labdn->a[i1][j1];
						b = labdn->b[i1][j1];
						float c_h=sqrt(SQR(a)+SQR(b));
						if(c_h>3000.f){
							a*=1.f + Qhigh*realred/100.f;
							b*=1.f + Qhigh*realblue/100.f;
						}
						
						Color::Lab2XYZ(L, a, b, X, Y, Z);

						float factor = Vmask[i1]*Hmask[j1];
						float r_,g_,b_;
						Color::xyz2rgb(X,Y,Z,r_,g_,b_,wip);
						//gamma slider is different from Raw
						r_ = r_<32768.0f ? igamcurve[r_] : (Color::gamman((float)r_/32768.0f, igam) * 65535.0f);
						g_ = g_<32768.0f ? igamcurve[g_] : (Color::gamman((float)g_/32768.0f, igam) * 65535.0f);
						b_ = b_<32768.0f ? igamcurve[b_] : (Color::gamman((float)b_/32768.0f, igam) * 65535.0f);

						dsttmp->r(i,j) += factor*r_;
						dsttmp->g(i,j) += factor*g_;
						dsttmp->b(i,j) += factor*b_;

					}
				}
			}

			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			delete labdn;
			delete [] mad_LL;
			delete [] mad_aa;
			delete [] mad_bb;

		}//end of tile row
	}//end of tile loop
	for (int i=0; i<(tileheight+1)/2; i++)
		delete [] noisevarlum[i];
	delete [] noisevarlum;
	for (int i=0; i<(tileheight+1)/2; i++)
		delete [] noisevarchrom[i];
	delete [] noisevarchrom;

	}
	for(int i=0;i<denoiseNestedLevels*numthreads;i++) {
		fftwf_free(LbloxArray[i]);
		fftwf_free(fLbloxArray[i]);
	}

#ifdef _OPENMP
omp_set_nested(oldNested);
#endif
	//copy denoised image to output
	memcpy (dst->data, dsttmp->data, 3*dst->width*dst->height*sizeof(float));
	if (!isRAW) {//restore original image gamma
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<3*dst->width*dst->height; i++) {
			dst->data[i] = Color::gammatab_srgb[ dst->data[i] ];
		}
	}

	delete dsttmp;
// destroy the plans
fftwf_destroy_plan( plan_forward_blox[0] );
fftwf_destroy_plan( plan_backward_blox[0] );
fftwf_destroy_plan( plan_forward_blox[1] );
fftwf_destroy_plan( plan_backward_blox[1] );
fftwf_cleanup();
}


//median 3x3 in complement on RGB
if(dnparams.methodmed=="RGB" && dnparams.median) {
//printf("RGB den\n");
	int wid=dst->width, hei=dst->height;
	float** tm;
	tm = new float*[hei];
	for (int i=0; i<hei; i++)
		tm[i] = new float[wid];

	Imagefloat *source;
	if (dnparams.luma==0 && dnparams.chroma==0)
		source = dst;
	else
		source = src;

	int methmed=0;
	int border = 1;
	if(dnparams.rgbmethod=="soft")
		methmed=0;
	else if(dnparams.rgbmethod=="33")
		methmed=1;
	else if(dnparams.rgbmethod=="55") {
		methmed = 3;
		border = 2;
	}
	else if(dnparams.rgbmethod=="55soft") {
		methmed = 2;
		border = 2;
	}

for(int iteration=1;iteration<=dnparams.passes;iteration++){

#pragma omp parallel
{
	if(methmed < 2) {
#pragma omp for	
		for (int i=1; i<hei-1; i++) {
			float pp[9],temp;
			if(methmed == 0)
				for (int j=1; j<wid-1; j++) {
					med2(source->r(i,j),source->r(i-1,j),source->r(i+1,j),source->r(i,j+1),source->r(i,j-1),tm[i][j]);//3x3 soft
					}
			else
				for (int j=1; j<wid-1; j++) {
					med3(source->r(i,j),source->r(i-1,j),source->r(i+1,j),source->r(i,j+1),source->r(i,j-1),source->r(i-1,j-1),source->r(i-1,j+1),source->r(i+1,j-1),source->r(i+1,j+1),tm[i][j]);//3x3
				}
		}
	} else {
#pragma omp for	
		for (int i=2; i<hei-2; i++) {
			float pp[25],temp;
			if(methmed == 3)			
				for (int j=2; j<wid-2; j++) {
					med5(source->r(i,j),source->r(i-1,j),source->r(i+1,j),source->r(i,j+1),source->r(i,j-1),source->r(i-1,j-1),source->r(i-1,j+1),source->r(i+1,j-1),source->r(i+1,j+1),
					source->r(i-2,j),source->r(i+2,j),source->r(i,j+2),source->r(i,j-2),source->r(i-2,j-2),source->r(i-2,j+2),source->r(i+2,j-2),source->r(i+2,j+2),	
					source->r(i-2,j+1),source->r(i+2,j+1),source->r(i-1,j+2),source->r(i-1,j-2),source->r(i-2,j-1),source->r(i+2,j-1),source->r(i+1,j+2),source->r(i+1,j-2),	
					tm[i][j]);//5x5
				}
			else
				for (int j=2; j<wid-2; j++) {
					pp[0]=source->r(i,j);pp[1]=source->r(i-1,j); pp[2]=source->r(i+1,j);pp[3]=source->r(i,j+1);pp[4]=source->r(i,j-1);pp[5]=source->r(i-1,j-1);pp[6]=source->r(i-1,j+1);
					pp[7]=source->r(i+1,j-1);pp[8]=source->r(i+1,j+1);pp[9]=source->r(i+2,j);pp[10]=source->r(i-2,j);pp[11]=source->r(i,j+2);pp[12]=source->r(i,j-2);
					fq_sort2(pp,13);
					tm[i][j]=pp[6];//5x5 soft
				}
			}
	}
#ifdef _OPENMP
#pragma omp for nowait
#endif
	for(int i = border; i < hei-border; i++ ) {
		for(int j = border; j < wid-border; j++) {
			dst->r(i,j) = tm[i][j];
		}
	}

	if(methmed < 2) {
#pragma omp for		
		for (int i=1; i<hei-1; i++) {
			float pp[9],temp;
			if(methmed == 0)
				for (int j=1; j<wid-1; j++) {
					med2(source->b(i,j),source->b(i-1,j),source->b(i+1,j),source->b(i,j+1),source->b(i,j-1),tm[i][j]);
				}
			else 
				for (int j=1; j<wid-1; j++) {
					med3(source->b(i,j),source->b(i-1,j),source->b(i+1,j),source->b(i,j+1),source->b(i,j-1),source->b(i-1,j-1),source->b(i-1,j+1),source->b(i+1,j-1),source->b(i+1,j+1),tm[i][j]);
				}
		}
	} else {
#pragma omp for	
		for (int i=2; i<hei-2; i++) {
			float pp[25],temp;
			if(methmed == 3)							
				for (int j=2; j<wid-2; j++) {
					med5(source->b(i,j),source->b(i-1,j),source->b(i+1,j),source->b(i,j+1),source->b(i,j-1),source->b(i-1,j-1),source->b(i-1,j+1),source->b(i+1,j-1),source->b(i+1,j+1),
					source->b(i-2,j),source->b(i+2,j),source->b(i,j+2),source->b(i,j-2),source->b(i-2,j-2),source->b(i-2,j+2),source->b(i+2,j-2),source->b(i+2,j+2),	
					source->b(i-2,j+1),source->b(i+2,j+1),source->b(i-1,j+2),source->b(i-1,j-2),source->b(i-2,j-1),source->b(i+2,j-1),source->b(i+1,j+2),source->b(i+1,j-2),	
					tm[i][j]);//5x5
				}
			else
				for (int j=2; j<wid-2; j++) {
					pp[0]=source->b(i,j);pp[1]=source->b(i-1,j); pp[2]=source->b(i+1,j);pp[3]=source->b(i,j+1);pp[4]=source->b(i,j-1);pp[5]=source->b(i-1,j-1);pp[6]=source->b(i-1,j+1);
					pp[7]=source->b(i+1,j-1);pp[8]=source->b(i+1,j+1);pp[9]=source->b(i+2,j);pp[10]=source->b(i-2,j);pp[11]=source->b(i,j+2);pp[12]=source->b(i,j-2);
					fq_sort2(pp,13);
					tm[i][j]=pp[6];//5x5 soft
				}
		}
	}

#ifdef _OPENMP
#pragma omp for nowait
#endif
	for(int i = border; i < hei-border; i++ ) {
		for(int j = border; j < wid-border; j++) {
			dst->b(i,j) = tm[i][j];	
		}
	}


	if(methmed < 2) {
#pragma omp for		
		for (int i=1; i<hei-1; i++) {
			float pp[9],temp;
			if(methmed == 0)
				for (int j=1; j<wid-1; j++) {
					med2(source->g(i,j),source->g(i-1,j),source->g(i+1,j),source->g(i,j+1),source->g(i,j-1),tm[i][j]);
				}
			else
				for (int j=1; j<wid-1; j++) {
					med3(source->g(i,j),source->g(i-1,j),source->g(i+1,j),source->g(i,j+1),source->g(i,j-1),source->g(i-1,j-1),source->g(i-1,j+1),source->g(i+1,j-1),source->g(i+1,j+1),tm[i][j]);
				}
		}
	} else {
#pragma omp for	
		for (int i=2; i<hei-2; i++) {
			float pp[25],temp;
			if(methmed == 3)											
				for (int j=2; j<wid-2; j++) {
					med5(source->g(i,j),source->g(i-1,j),source->g(i+1,j),source->g(i,j+1),source->g(i,j-1),source->g(i-1,j-1),source->g(i-1,j+1),source->g(i+1,j-1),source->g(i+1,j+1),
					source->g(i-2,j),source->g(i+2,j),source->g(i,j+2),source->g(i,j-2),source->g(i-2,j-2),source->g(i-2,j+2),source->g(i+2,j-2),source->g(i+2,j+2),	
					source->g(i-2,j+1),source->g(i+2,j+1),source->g(i-1,j+2),source->g(i-1,j-2),source->g(i-2,j-1),source->g(i+2,j-1),source->g(i+1,j+2),source->g(i+1,j-2),	
					tm[i][j]);//5x5
				}
			else
				for (int j=2; j<wid-2; j++) {
					pp[0]=source->g(i,j);pp[1]=source->g(i-1,j); pp[2]=source->g(i+1,j);pp[3]=source->g(i,j+1);pp[4]=source->g(i,j-1);pp[5]=source->g(i-1,j-1);pp[6]=source->g(i-1,j+1);
					pp[7]=source->g(i+1,j-1);pp[8]=source->g(i+1,j+1);pp[9]=source->g(i+2,j);pp[10]=source->g(i-2,j);pp[11]=source->g(i,j+2);pp[12]=source->g(i,j-2);
					fq_sort2(pp,13);
					tm[i][j]=pp[6];//5x5 soft
				}
		}
	}

#ifdef _OPENMP
#pragma omp for
#endif
	for(int i = border; i < hei-border; i++ ) {
		for(int j = border; j < wid-border; j++) {
			dst->g(i,j) = tm[i][j];
		}
	}
}
}		
			for (int i=0; i<hei; i++)
				delete [] tm[i];
			delete [] tm;
		
}		
	//end median
	if(noiseLCurve || noiseCCurve) {
			for (int i=0; i<(hei+1)/2; i++)
				delete [] lumcalc[i];
			delete [] lumcalc;
			for (int i=0; i<(hei+1)/2; i++)
				delete [] acalc[i];
			delete [] acalc;
			for (int i=0; i<(hei+1)/2; i++)
				delete [] bcalc[i];
			delete [] bcalc;
			
			}

//#ifdef _DEBUG
	if (settings->verbose) {
		t2e.set();
		printf("Denoise performed in %d usec:\n", t2e.etime(t1e));
	}
//#endif

}//end of main RGB_denoise


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SSEFUNCTION	void ImProcFunctions::RGBtile_denoise (float * fLblox, int hblproc, float noisevar_Ldetail, float * nbrwt, float * blurbuffer )	//for DCT
	{
		int blkstart = hblproc*TS*TS;

		boxabsblur(fLblox+blkstart, nbrwt, 3, 3, TS, TS, blurbuffer);//blur neighbor weights for more robust estimation	//for DCT

#ifdef __SSE2__
		__m128	tempv;
		__m128	noisevar_Ldetailv = _mm_set1_ps( noisevar_Ldetail );
		__m128	onev = _mm_set1_ps( 1.0f );
		for (int n=0; n<TS*TS; n+=4) {		//for DCT
			tempv  = onev - xexpf( -SQRV( LVF(nbrwt[n]))/noisevar_Ldetailv);
			_mm_storeu_ps( &fLblox[blkstart+n], LVFU(fLblox[blkstart+n]) * tempv );
		}//output neighbor averaged result
#else
		for (int n=0; n<TS*TS; n++) {		//for DCT
			fLblox[blkstart+n] *= (1-xexpf(-SQR(nbrwt[n])/noisevar_Ldetail));
		}//output neighbor averaged result
#endif

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//printf("vblk=%d  hlk=%d  wsqave=%f   ||   ",vblproc,hblproc,wsqave);

	}//end of function tile_denoise


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	void ImProcFunctions::RGBoutput_tile_row (float *bloxrow_L, float ** Ldetail, float ** tilemask_out, int height, int width, int top )
	{
		const int numblox_W = ceil(((float)(width))/(offset));
		const float DCTnorm = 1.0f/(4*TS*TS); //for DCT

		int imin = MAX(0,-top);
		int bottom = MIN( top+TS,height);
		int imax = bottom - top;

		//add row of tiles to output image
		for (int hblk=0; hblk < numblox_W; hblk++) {
			int left = (hblk-blkrad)*offset;
			int right  = MIN(left+TS, width);
			int jmin = MAX(0,-left);
			int jmax = right - left;
			int indx = hblk*TS;

			for (int i=imin; i<imax; i++)
				for (int j=jmin; j<jmax; j++) {	// this loop gets auto vectorized by gcc
					Ldetail[top+i][left+j] += tilemask_out[i][j]*bloxrow_L[(indx + i)*TS+j]*DCTnorm; //for DCT

				}
		}
	}
/*
#undef TS
#undef fTS
#undef offset
#undef epsilon
*/

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	float ImProcFunctions::MadMax(float * DataList, int & max, int datalen) {

		//computes Median Absolute Deviation and Maximum of DataList
		//DataList values should mostly have abs val < 65535

		int * histo = new int[65536];
		//memset(histo, 0, 65536*sizeof(histo));
		for (int i=0; i<65536; i++) histo[i]=0;

		//calculate histogram of absolute values of HH wavelet coeffs
		for (int i=0; i<datalen; i++) {
			histo[MAX(0,MIN(65535,abs((int)DataList[i])))]++;
		}

		//find median of histogram
		int median=0, count=0;
		while (count<datalen/2) {
			count += histo[median];
			median++;
		}

		//find max of histogram
		max=65535;
		while (histo[max]==0) {
			max--;
		}

		int count_ = count - histo[median-1];

		delete[] histo;

		// interpolate
		return (( (median-1) + (datalen/2-count_)/((float)(count-count_)) )/0.6745);

	}

	float ImProcFunctions::Mad( float * DataList, const int datalen) {

		//computes Median Absolute Deviation
		//DataList values should mostly have abs val < 256 because we are in Lab mode
		int histo[256] ALIGNED64 = {0};
		
		//calculate histogram of absolute values of wavelet coeffs
		for (int i=0; i<datalen; i++) {
			histo[min(255,abs((int)DataList[i]))]++;			
		}

		//find median of histogram
		int median = 0, count = 0;
		while (count<datalen/2) {
			count += histo[median];
			median++;
		}

		int count_ = count - histo[median-1];

		// interpolate
		return (( (median-1) + (datalen/2-count_)/((float)(count-count_)) )/0.6745);
	}

	float ImProcFunctions::MadRgb( float * DataList, const int datalen) {

		//computes Median Absolute Deviation
		//DataList values should mostly have abs val < 65536 because we are in RGB mode
		int * histo = new int[65536];
		for(int i=0;i<65536;i++)
			histo[i] = 0;
		
		//calculate histogram of absolute values of wavelet coeffs
		int i;

		for (i=0; i<datalen; i++) {
			histo[min(65536,abs((int)DataList[i]))]++;			
		}

		//find median of histogram
		int median = 0, count = 0;
		while (count<datalen/2) {
			count += histo[median];
			median++;
		}

		int count_ = count - histo[median-1];

		// interpolate
		delete [] histo;
		return (( (median-1) + (datalen/2-count_)/((float)(count-count_)) )/0.6745);
	}



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ImProcFunctions::Noise_residual(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a, wavelet_decomposition &WaveletCoeffs_b,  int width, int height, float &chresid, float &chmaxredresid,float &chmaxblueresid, float &chresidred, float & chresidblue, float resid, float residblue, float residred, float maxredresid, float maxblueresid, float nbresid)
{
		int maxlvl = WaveletCoeffs_L.maxlevel();
		const float eps = 0.01f;

		float madaC,madbC;
		for (int lvl=0; lvl<maxlvl; lvl++) {
			// compute median absolute deviation (MAD) of detail coefficients as robust noise estimator

			int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
			int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);

			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);
			for (int dir=1; dir<4; dir++) {
				madaC = SQR(Mad(WavCoeffs_a[dir], Wlvl_ab*Hlvl_ab));
				madbC = SQR(Mad(WavCoeffs_b[dir], Wlvl_ab*Hlvl_ab));

				resid +=(madaC+madbC);
				residred+=madaC;
				residblue+=madbC;
				
				if(madaC >maxredresid ) maxredresid=madaC;
				if(madbC > maxblueresid) maxblueresid=madbC;
				nbresid++;	
			}
			chresid=sqrt(resid/(2*nbresid));
			chmaxredresid=sqrt(maxredresid);
			chmaxblueresid=sqrt(maxblueresid);
			chresidred=sqrt(residred/nbresid);
			chresidblue=sqrt(residblue/nbresid);
		}
}

SSEFUNCTION	void ImProcFunctions::WaveletDenoiseAll_BiShrink(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a,
													 wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float **noisevarlum, float **noisevarchrom, int width, int height, float *mad_LL, float *mad_aa, float *mad_bb, float noisevar_abr, float noisevar_abb, LabImage * noi, const NoiseCurve & noiseLCurve, const NoiseCurve & noiseCCurve, float &chaut, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, int schoice, bool autoch, bool perf)
	{
		int maxlvl = WaveletCoeffs_L.maxlevel();
		const float eps = 0.01f;
		if(autoch && noisevar_abr <=0.001f) noisevar_abr=0.02f;
		if(autoch && noisevar_abb <=0.001f) noisevar_abb=0.02f;

		float madL[8][3], mada[8][3], madb[8][3];
		
		int maxWL = 0, maxHL = 0;
		for (int lvl=0; lvl<maxlvl; lvl++) {
			if(WaveletCoeffs_L.level_W(lvl) > maxWL)
				maxWL = WaveletCoeffs_L.level_W(lvl);
			if(WaveletCoeffs_L.level_H(lvl) > maxHL)
				maxHL = WaveletCoeffs_L.level_H(lvl);
		}

#ifdef _OPENMP
#pragma omp parallel num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
{
		float *buffer[5];
		buffer[0] = new float[maxWL*maxHL+32];
		buffer[1] = new float[maxWL*maxHL+64];
		buffer[2] = new float[maxWL*maxHL+96];
		buffer[3] = new float[maxWL*maxHL+128];
		buffer[4] = new float[maxWL*maxHL+160];
		

#ifdef _OPENMP
#pragma omp for
#endif
		for (int lvl=0; lvl<maxlvl; lvl++) {
			// compute median absolute deviation (MAD) of detail coefficients as robust noise estimator

			int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
			int Hlvl_L = WaveletCoeffs_L.level_H(lvl);
			
			int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
			int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);

			float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);

			for (int dir=1; dir<4; dir++) {
				if(!perf) {
					madL[lvl][dir-1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L*Hlvl_L));
					mada[lvl][dir-1] = SQR(Mad(WavCoeffs_a[dir], Wlvl_ab*Hlvl_ab));
					madb[lvl][dir-1] = SQR(Mad(WavCoeffs_b[dir], Wlvl_ab*Hlvl_ab));
				} else {
					madL[lvl][dir-1] = SQR(MadRgb(WavCoeffs_L[dir], Wlvl_L*Hlvl_L));
					mada[lvl][dir-1] = SQR(MadRgb(WavCoeffs_a[dir], Wlvl_ab*Hlvl_ab));
					madb[lvl][dir-1] = SQR(MadRgb(WavCoeffs_b[dir], Wlvl_ab*Hlvl_ab));
				}
			}	
		}
		
#ifdef _OPENMP
#pragma omp for
#endif
		for (int lvl=maxlvl-1; lvl>=0; lvl--) {//for levels less than max, use level diff to make edge mask
			int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
			int Hlvl_L = WaveletCoeffs_L.level_H(lvl);
			
			int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
			int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);

			int skip_L = WaveletCoeffs_L.level_stride(lvl);
			int skip_ab = WaveletCoeffs_a.level_stride(lvl);

			float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);
			int callby=1;
			if (lvl==maxlvl-1) {
				ShrinkAll(WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, buffer, lvl, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,
						  skip_L, skip_ab, noisevar_L, noisevarlum,  noisevarchrom, width, height, mad_LL, mad_aa, mad_bb,noisevar_abr, noisevar_abb, noi, noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, callby, autoch, perf, mada[lvl], madb[lvl], madL[lvl], true);

			} else {
				//simple wavelet shrinkage
				float * sfave = buffer[0]+32;
				float * sfaved = buffer[3]+128;
				float * blurBuffer = buffer[1]+64;
				for (int dir=1; dir<4; dir++) {
					float mad_Lr = madL[lvl][dir-1];
					float mad_ar = noisevar_abr*mada[lvl][dir-1];
					float mad_br = noisevar_abb*madb[lvl][dir-1];

					if (!noiseCCurve || noiseCCurve.getSum() < 5.f ){
						for (int i=0; i<Hlvl_ab; i++) {
							for (int j=0; j<Wlvl_ab; j++) {
								mad_aa[i*Wlvl_ab+j]=mad_ar*(noisevar_abr);
								mad_bb[i*Wlvl_ab+j]=mad_br*(noisevar_abb);
							}//noisevarchrom	
						}
					}	

					if (noiseCCurve  && noiseCCurve.getSum() > 5.f ){
						for (int i=0; i<Hlvl_ab; i++) {
							for (int j=0; j<Wlvl_ab; j++) {
									mad_aa[i*Wlvl_ab+j]=mad_ar*(noisevarchrom[i][j]);
									mad_bb[i*Wlvl_ab+j]=mad_br*(noisevarchrom[i][j]);
							}//noisevarchrom
						}
					}
					if (noisevar_abr>0.001f  || noisevar_abb>0.001f) {
					
#ifdef __SSE2__
						__m128 onev = _mm_set1_ps(1.f);
						__m128 rmad_Lm9v = onev / _mm_set1_ps(mad_Lr * 9.f);
						__m128 mad_av;
						__m128 mad_bv;
						__m128 mag_Lv, mag_av, mag_bv;
						__m128 tempav, tempbv;
						int coeffloc_ab;
						for (coeffloc_ab=0; coeffloc_ab<Hlvl_ab*Wlvl_ab-3; coeffloc_ab+=4) {
							mad_av = LVFU(mad_aa[coeffloc_ab]);
							mad_bv = LVFU(mad_bb[coeffloc_ab]);
						
							tempav = LVFU(WavCoeffs_a[dir][coeffloc_ab]);
							tempbv = LVFU(WavCoeffs_b[dir][coeffloc_ab]);
							mag_Lv = LVFU(WavCoeffs_L[dir][coeffloc_ab]);
							mag_av = SQRV(tempav);
							mag_bv = SQRV(tempbv);
							mag_Lv = SQRV(mag_Lv) * rmad_Lm9v;
							_mm_storeu_ps(&WavCoeffs_a[dir][coeffloc_ab], tempav * SQRV((onev-xexpf(-(mag_av/mad_av)-(mag_Lv)))));
							_mm_storeu_ps(&WavCoeffs_b[dir][coeffloc_ab], tempbv * SQRV((onev-xexpf(-(mag_bv/mad_bv)-(mag_Lv)))));
						}
						// few remaining pixels
						for (; coeffloc_ab<Hlvl_ab*Wlvl_ab; coeffloc_ab++) {
							float mag_L = SQR(WavCoeffs_L[dir][coeffloc_ab ]);
							float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab]);
							float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab]);
							WavCoeffs_a[dir][coeffloc_ab] *= SQR(1.f-xexpf(-(mag_a/mad_aa[coeffloc_ab])-(mag_L/(9.f*mad_Lr)))/*satfactor_a*/);
							WavCoeffs_b[dir][coeffloc_ab] *= SQR(1.f-xexpf(-(mag_b/mad_bb[coeffloc_ab])-(mag_L/(9.f*mad_Lr)))/*satfactor_b*/);
						}//now chrominance coefficients are denoised
#else
						for (int i=0; i<Hlvl_ab; i++) {
							for (int j=0; j<Wlvl_ab; j++) {
								int coeffloc_ab = i*Wlvl_ab+j;

								float mag_L = SQR(WavCoeffs_L[dir][coeffloc_ab ])+eps;
								float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
								float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;

								WavCoeffs_a[dir][coeffloc_ab] *= SQR(1.f-xexpf(-(mag_a/mad_aa[coeffloc_ab])-(mag_L/(9.f*mad_Lr)))/*satfactor_a*/);
								WavCoeffs_b[dir][coeffloc_ab] *= SQR(1.f-xexpf(-(mag_b/mad_bb[coeffloc_ab])-(mag_L/(9.f*mad_Lr)))/*satfactor_b*/);

							}
						}//now chrominance coefficients are denoised
#endif
					}


					float levelFactor = mad_Lr*5.f/(lvl+1);
					if (noisevar_L>0.00001f) {
						if (!noiseLCurve || noiseLCurve.getSum() < 7.f ) {
							for (int i=0; i<Hlvl_L; i++)
								for (int j=0; j<Wlvl_L; j++) {
									int coeffloc_L = i*Wlvl_L+j;
									mad_LL[coeffloc_L] = noisevar_L*levelFactor;//noisevarlum
								}
							}
						else { 	if (noiseLCurve && noiseLCurve.getSum() >= 7.f) {	
							for (int i=0; i<Hlvl_L; i++)
								for (int j=0; j<Wlvl_L; j++) {
									int coeffloc_L = i*Wlvl_L+j;
									mad_LL[coeffloc_L] = noisevarlum[i][j]*levelFactor;//noisevarlum
								}
							}
						}
#ifdef __SSE2__
						__m128 mad_Lv;
						__m128 ninev = _mm_set1_ps( 9.0f );
						__m128 epsv = _mm_set1_ps(eps);
						__m128 mag_Lv;
						int coeffloc_L;
						for (coeffloc_L=0; coeffloc_L<Hlvl_L*Wlvl_L-3; coeffloc_L+=4) {
							mad_Lv = LVFU(mad_LL[coeffloc_L]);
							mag_Lv = SQRV(LVFU(WavCoeffs_L[dir][coeffloc_L]));
							_mm_storeu_ps(&sfave[coeffloc_L], mag_Lv / ( mag_Lv + mad_Lv * xexpf(-mag_Lv/(mad_Lv*ninev) )+ epsv));
						}	
						for (; coeffloc_L<Hlvl_L*Wlvl_L; coeffloc_L++) {
							float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
							sfave[coeffloc_L] = mag_L/(mag_L+mad_LL[coeffloc_L]*xexpf(-mag_L/(9.f*mad_LL[coeffloc_L]))+eps);
						}
#else
						for (int i=0; i<Hlvl_L; i++)
							for (int j=0; j<Wlvl_L; j++) {

								int coeffloc_L = i*Wlvl_L+j;
								float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
								sfave[coeffloc_L] = mag_L/(mag_L+mad_LL[coeffloc_L]*xexpf(-mag_L/(9.f*mad_LL[coeffloc_L]))+eps);
							}
#endif
						boxblur(sfave, sfaved, blurBuffer, lvl+2, lvl+2, Wlvl_L, Hlvl_L);//increase smoothness by locally averaging shrinkage
#ifdef __SSE2__
						__m128 sfavev;
						__m128 sf_Lv;
						for (coeffloc_L=0; coeffloc_L<Hlvl_L*Wlvl_L-3; coeffloc_L+=4) {
							sfavev = LVFU(sfaved[coeffloc_L]);
							sf_Lv = LVFU(sfave[coeffloc_L]);
							_mm_storeu_ps(&WavCoeffs_L[dir][coeffloc_L], LVFU(WavCoeffs_L[dir][coeffloc_L]) * (SQRV(sfavev)+SQRV(sf_Lv))/(sfavev+sf_Lv+epsv));
							//use smoothed shrinkage unless local shrinkage is much less
						}
						// few remaining pixels
						for (; coeffloc_L<Hlvl_L*Wlvl_L; coeffloc_L++) {
							float sf_L = sfave[coeffloc_L];
							//use smoothed shrinkage unless local shrinkage is much less
							WavCoeffs_L[dir][coeffloc_L] *= (SQR(sfaved[coeffloc_L])+SQR(sf_L))/(sfaved[coeffloc_L]+sf_L+eps);
						}//now luminance coeffs are denoised
#else
						for (int i=0; i<Hlvl_L; i++)
							for (int j=0; j<Wlvl_L; j++) {
								int coeffloc_L = i*Wlvl_L+j;
								float sf_L = sfave[coeffloc_L];
								//use smoothed shrinkage unless local shrinkage is much less
								WavCoeffs_L[dir][coeffloc_L] *= (SQR(sfaved[coeffloc_L])+SQR(sf_L))/(sfaved[coeffloc_L]+sf_L+eps);
							}//now luminance coeffs are denoised
#endif
					}
				}
			}
		}
		for(int i=4;i>=0;i--)
			delete [] buffer[i];

}
	}


	void ImProcFunctions::WaveletDenoiseAll(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a,
											wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float **noisevarlum, float **noisevarchrom, int width, int height, float *mad_LL, float *mad_aa, float *mad_bb, float noisevar_abr, float noisevar_abb, LabImage * noi, const NoiseCurve & noiseLCurve, const NoiseCurve & noiseCCurve, float &chaut,float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, int schoice, bool autoch, bool perf)//mod JD

	{

		int maxlvl = WaveletCoeffs_L.maxlevel();
		int maxWL = 0, maxHL = 0;
		for (int lvl=0; lvl<maxlvl; lvl++) {
			if(WaveletCoeffs_L.level_W(lvl) > maxWL)
				maxWL = WaveletCoeffs_L.level_W(lvl);
			if(WaveletCoeffs_L.level_H(lvl) > maxHL)
				maxHL = WaveletCoeffs_L.level_H(lvl);
		}


#ifdef _OPENMP
#pragma omp parallel num_threads(denoiseNestedLevels) if(denoiseNestedLevels>1)
#endif
{
		float *buffer[5];
		buffer[0] = new float[maxWL*maxHL+32];
		buffer[1] = new float[maxWL*maxHL+64];
		buffer[2] = new float[maxWL*maxHL+96];
		buffer[3] = new float[maxWL*maxHL+128];
		buffer[4] = new float[maxWL*maxHL+160];
		
		
#ifdef _OPENMP
#pragma omp for
#endif

		for (int lvl=0; lvl<maxlvl; lvl++) {

			int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
			int Hlvl_L = WaveletCoeffs_L.level_H(lvl);

			int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
			int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);


			int skip_L = WaveletCoeffs_L.level_stride(lvl);
			int skip_ab = WaveletCoeffs_a.level_stride(lvl);

			float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);

			int callby=0;
			if(schoice==0) callby=0;
			if(schoice==2) callby=1;
			
			ShrinkAll(WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, buffer, lvl, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,
					  skip_L, skip_ab, noisevar_L, noisevarlum,  noisevarchrom, width, height, mad_LL, mad_aa, mad_bb, noisevar_abr, noisevar_abb, noi, noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, callby, autoch, perf);

		}
		for(int i=4;i>=0;i--)
			delete [] buffer[i];
}
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SSEFUNCTION	void ImProcFunctions::ShrinkAll(float ** WavCoeffs_L, float ** WavCoeffs_a, float ** WavCoeffs_b, float **buffer, int level, 
									int W_L, int H_L, int W_ab, int H_ab,int skip_L, int skip_ab, float noisevar_L, float **noisevarlum, float **noisevarchrom, int width, int height, float * mad_LL, float * mad_aa, float * mad_bb, float noisevar_abr, float noisevar_abb, LabImage * noi, const NoiseCurve & noiseLCurve, const NoiseCurve & noiseCCurve, float &chaut, float &redaut, float &blueaut,
									float &maxredaut, float &maxblueaut, int callby, bool autoch, bool perf, float * madaa, float * madab, float * madaL, bool madCalculated )

									{
		//simple wavelet shrinkage
		const float eps = 0.01f;
		if(autoch && noisevar_abr <=0.001f)
			noisevar_abr=0.02f;
		if(autoch && noisevar_abb <=0.001f)
			noisevar_abb=0.02f;
		
		float * sfavea = buffer[0]+32;
		float * sfaveb = buffer[1]+64;
		float * sfavead = buffer[3]+128;
		float * sfavebd = buffer[4]+160;
		float * sfave = sfavea; // we can safely reuse sfavea here, because they are not used together
		float * sfaved = sfavead; // we can safely reuse sfavead here, because they are not used together
		float * blurBuffer = buffer[2]+96;

		for (int dir=1; dir<4; dir++) {
			float mada, madb, madL;
			if(madCalculated) {
				mada = madaa[dir-1];
				madb = madab[dir-1];
				madL = madaL[dir-1] ;
			} else {
				if(!perf) {
					madL = SQR(Mad(WavCoeffs_L[dir], W_L*H_L));
					mada = SQR(Mad(WavCoeffs_a[dir], W_ab*H_ab));
					madb = SQR(Mad(WavCoeffs_b[dir], W_ab*H_ab));
				} else {
					madL = SQR(MadRgb(WavCoeffs_L[dir], W_L*H_L));
					mada = SQR(MadRgb(WavCoeffs_a[dir], W_ab*H_ab));
					madb = SQR(MadRgb(WavCoeffs_b[dir], W_ab*H_ab));
				}
			}
		
			if (!noiseCCurve || noiseCCurve.getSum() < 5.f ){
				for (int i=0; i<H_ab; i++) {
					for (int j=0; j<W_ab; j++) {
						mad_aa[i*W_ab+j]=mada*(noisevar_abr);
						mad_bb[i*W_ab+j]=madb*(noisevar_abb);
					}//noisevarchrom	
				}
			}	
			if (noiseCCurve  && noiseCCurve.getSum() > 5.f ){
				for (int i=0; i<H_ab; i++) {
					for (int j=0; j<W_ab; j++) {
						mad_aa[i*W_ab+j]=mada*(noisevarchrom[i][j]);
						mad_bb[i*W_ab+j]=madb*(noisevarchrom[i][j]);
					}//noisevarchrom
				}
			}
			
			if (noisevar_abr>0.001f  ||  noisevar_abb>0.001f ) {
#ifdef __SSE2__
				__m128 onev = _mm_set1_ps(1.f);
				__m128 rmadLm9v = onev / _mm_set1_ps(madL * 9.f);
				__m128 mad_av ;
				__m128 mad_bv;
				__m128 mag_Lv, mag_av, mag_bv;
				int coeffloc_ab;
				for (coeffloc_ab=0; coeffloc_ab<H_ab*W_ab-3; coeffloc_ab+=4) {
					mad_av = LVFU(mad_aa[coeffloc_ab]);
					mad_bv = LVFU(mad_bb[coeffloc_ab]);

					mag_Lv = LVFU(WavCoeffs_L[dir][coeffloc_ab]);
					mag_av = SQRV(LVFU(WavCoeffs_a[dir][coeffloc_ab]));
					mag_bv = SQRV(LVFU(WavCoeffs_b[dir][coeffloc_ab]));
					mag_Lv = (SQRV(mag_Lv)) * rmadLm9v;
					_mm_storeu_ps(&sfavea[coeffloc_ab], (onev-xexpf(-(mag_av/mad_av)-(mag_Lv))));
					_mm_storeu_ps(&sfaveb[coeffloc_ab], (onev-xexpf(-(mag_bv/mad_bv)-(mag_Lv))));
				}
				// few remaining pixels
				for (; coeffloc_ab<H_ab*W_ab; coeffloc_ab++) {
					float mag_L = SQR(WavCoeffs_L[dir][coeffloc_ab]);
					float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab]);
					float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab]);
					sfavea[coeffloc_ab] = (1.f-xexpf(-(mag_a/mad_aa[coeffloc_ab])-(mag_L/(9.f*madL))));
					sfaveb[coeffloc_ab] = (1.f-xexpf(-(mag_b/mad_bb[coeffloc_ab])-(mag_L/(9.f*madL))));
				}//now chrominance coefficients are denoised

#else
				for (int i=0; i<H_ab; i++) {
					for (int j=0; j<W_ab; j++) {
						int coeffloc_ab = i*W_ab+j;
						int coeffloc_L	= ((i*skip_L)/skip_ab)*W_L + ((j*skip_L)/skip_ab);
						float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L ]);
						float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab]);
						float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab]);
						sfavea[coeffloc_ab] = (1.f-xexpf(-(mag_a/mad_aa[coeffloc_ab])-(mag_L/(9.f*madL))));
						sfaveb[coeffloc_ab] = (1.f-xexpf(-(mag_b/mad_bb[coeffloc_ab])-(mag_L/(9.f*madL))));
					}
				}//now chrominance coefficients are denoised
#endif

				boxblur(sfavea, sfavead, blurBuffer, level+2, level+2, W_ab, H_ab);//increase smoothness by locally averaging shrinkage
				boxblur(sfaveb, sfavebd, blurBuffer, level+2, level+2, W_ab, H_ab);//increase smoothness by locally averaging shrinkage

#ifdef __SSE2__
				__m128 epsv = _mm_set1_ps(eps);
				__m128 sfav, sfbv;
				__m128 sfaveav, sfavebv;
				for (coeffloc_ab=0; coeffloc_ab<H_ab*W_ab-3; coeffloc_ab+=4) {
					sfav = LVFU(sfavea[coeffloc_ab]);
					sfbv = LVFU(sfaveb[coeffloc_ab]);
					sfaveav = LVFU(sfavead[coeffloc_ab]);
					sfavebv = LVFU(sfavebd[coeffloc_ab]);

					//use smoothed shrinkage unless local shrinkage is much less
					_mm_storeu_ps( &WavCoeffs_a[dir][coeffloc_ab], LVFU(WavCoeffs_a[dir][coeffloc_ab]) * (SQRV(sfaveav)+SQRV(sfav))/(sfaveav+sfav+epsv));
					_mm_storeu_ps( &WavCoeffs_b[dir][coeffloc_ab], LVFU(WavCoeffs_b[dir][coeffloc_ab]) * (SQRV(sfavebv)+SQRV(sfbv))/(sfavebv+sfbv+epsv));
				}
				// few remaining pixels
				for (; coeffloc_ab<H_ab*W_ab; coeffloc_ab++) {
					//modification Jacques feb 2013
					float sfa = sfavea[coeffloc_ab];
					float sfb = sfaveb[coeffloc_ab];

					//use smoothed shrinkage unless local shrinkage is much less
					WavCoeffs_a[dir][coeffloc_ab] *= (SQR(sfavead[coeffloc_ab])+SQR(sfa))/(sfavead[coeffloc_ab]+sfa+eps);
					WavCoeffs_b[dir][coeffloc_ab] *= (SQR(sfavebd[coeffloc_ab])+SQR(sfb))/(sfavebd[coeffloc_ab]+sfb+eps);
				}//now chrominance coefficients are denoised
#else
				for (int i=0; i<H_ab; i++) {
					for (int j=0; j<W_ab; j++) {
						int coeffloc_ab = i*W_ab+j;
						float sfa = sfavea[coeffloc_ab];
						float sfb = sfaveb[coeffloc_ab];

						//use smoothed shrinkage unless local shrinkage is much less
						WavCoeffs_a[dir][coeffloc_ab] *= (SQR(sfavead[coeffloc_ab])+SQR(sfa))/(sfavead[coeffloc_ab]+sfa+eps);
						WavCoeffs_b[dir][coeffloc_ab] *= (SQR(sfavebd[coeffloc_ab])+SQR(sfb))/(sfavebd[coeffloc_ab]+sfb+eps);
					}//now chrominance coefficients are denoised
				}
#endif
			}
			
			if (noisevar_L>0.00001f) {
				if (!noiseLCurve || noiseLCurve.getSum() < 7.f ){//under 7 quasi no action
					for (int i=0; i<H_L; i++) {
						for (int j=0; j<W_L; j++) {
							mad_LL[i*W_L+j]=madL*(noisevar_L)*5/(level+1);
						}//noisevarlum
					}
				}
				else if (noiseLCurve && noiseLCurve.getSum() >= 7.f) {
					float precalc = madL * 5.f / (float)(level+1);
					for (int i=0; i<H_L; i++) {
						for (int j=0; j<W_L; j++) {
							mad_LL[i*W_L+j]=precalc*(noisevarlum[i][j]);
						}//noisevarlum
					}
				}
#ifdef __SSE2__
				__m128	magv;
				__m128  mad_Lv;
				__m128	ninev = _mm_set1_ps( 9.0f );
				__m128 	epsv = _mm_set1_ps( eps );
				int i;
				for (i=0; i<W_L*H_L-3; i+=4) {
					mad_Lv = LVFU(mad_LL[i]);
					magv = SQRV(LVFU(WavCoeffs_L[dir][i]));
					_mm_storeu_ps( &sfave[i], magv / (magv + mad_Lv*xexpf(-magv/(ninev * mad_Lv)) + epsv));
				}
				// few remaining pixels
				for (; i<W_L*H_L; i++) {
					float mag = SQR(WavCoeffs_L[dir][i]);
					sfave[i] = mag/(mag+mad_LL[i]*xexpf(-mag/(9*mad_LL[i]))+eps);
				}
#else
				for (int i=0; i<W_L*H_L; i++) {

					float mag = SQR(WavCoeffs_L[dir][i]);
					float shrinkfactor = mag/(mag+mad_LL[i]*xexpf(-mag/(9*mad_LL[i]))+eps);
					sfave[i] = shrinkfactor;
				}
#endif
				boxblur(sfave, sfaved, blurBuffer, level+2, level+2, W_L, H_L);//increase smoothness by locally averaging shrinkage

#ifdef __SSE2__ 
				__m128	sfv;
				for (i=0; i<W_L*H_L-3; i+=4) {
					sfv = LVFU(sfave[i]);
					//use smoothed shrinkage unless local shrinkage is much less
					_mm_storeu_ps( &WavCoeffs_L[dir][i], _mm_loadu_ps( &WavCoeffs_L[dir][i]) * (SQRV( LVFU(sfaved[i] )) + SQRV(sfv)) / (LVFU(sfaved[i])+sfv+epsv));
				}
				// few remaining pixels
				for (; i<W_L*H_L; i++) {
					float sf = sfave[i];

					//use smoothed shrinkage unless local shrinkage is much less
					WavCoeffs_L[dir][i] *= (SQR(sfaved[i])+SQR(sf))/(sfaved[i]+sf+eps);
				}//now luminance coefficients are denoised

#else 
				for (int i=0; i<W_L*H_L; i++) {
					float sf = sfave[i];

					//use smoothed shrinkage unless local shrinkage is much less
					WavCoeffs_L[dir][i] *= (SQR(sfaved[i])+SQR(sf))/(sfaved[i]+sf+eps);

				}//now luminance coefficients are denoised
#endif
			}
		}
	}


SSEFUNCTION	void ImProcFunctions::ShrinkAll_info(float ** WavCoeffs_a, float ** WavCoeffs_b, int level,
									int W_ab, int H_ab, int skip_ab, float **noisevarlum, float **noisevarchrom, float **noisevarhue, int width, int height, float * mad_aa, float * mad_bb, float noisevar_abr, float noisevar_abb, LabImage * noi, float &chaut, int &Nb, float &redaut, float &blueaut,
									float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut,  bool autoch, int schoice, int lvl, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel,float &skinc, float &nsknc,
									float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool perf, bool multiThread) {

	//simple wavelet shrinkage
	if(lvl==1){//only one time
		float chro = 0.f;
		float dev = 0.f;
		float devL = 0.f;
		int nc = 0;
		int nL = 0;
		int nry = 0;
		float lume = 0.f;
		float red_yel = 0.f;
		float skin_c = 0.f;
		int nsk = 0;
		float nsk_nc = 0.f;
		for (int i=0; i<H_ab; i++) {
			for (int j=0; j<W_ab; j++) {
				chro += noisevarchrom[i][j];
				nc++;
				dev += SQR(noisevarchrom[i][j]-(chro/nc));
				if(noisevarhue[i][j] > -0.8f && noisevarhue[i][j] < 2.0f && noisevarchrom[i][j] > 10000.f) {//saturated red yellow
					red_yel += noisevarchrom[i][j];
					nry++;
				}	
				if(noisevarhue[i][j] > 0.f && noisevarhue[i][j] < 1.6f && noisevarchrom[i][j] < 10000.f) {//skin
					skin_c += noisevarchrom[i][j];
					nsk++;
				}	
				lume += noisevarlum[i][j];
				nL++;
				devL += SQR(noisevarlum[i][j]-(lume/nL));
			}
		}
		if(nc>0) {
			chromina=chro/nc;
			sigma=sqrt(dev/nc);
			nsknc=(float)nsk/(float)nc;
		} else {
			nsknc=(float)nsk;
		}
		if(nL>0) {
			lumema=lume/nL;
			sigma_L=sqrt(devL/nL);
		}
		if(nry>0)
			redyel=red_yel/nry;
		if(nsk>0)
			skinc=skin_c/nsk;
	}

	const float reduc = (schoice == 2) ? (float) settings->nrhigh : 1.f;
	for (int dir=1; dir<4; dir++) {
		float mada, madb;
		if(!perf)
			mada = SQR(Mad(WavCoeffs_a[dir], W_ab*H_ab));
		else
			mada = SQR(MadRgb(WavCoeffs_a[dir], W_ab*H_ab));
		chred += mada;
		if(mada > maxchred)
			maxchred = mada;
		if(mada < minchred)
			minchred = mada;
		maxredaut = sqrt(reduc*maxchred);
		minredaut = sqrt(reduc*minchred);

		if(!perf)
			madb = SQR(Mad(WavCoeffs_b[dir], W_ab*H_ab));
		else
			madb = SQR(MadRgb(WavCoeffs_b[dir], W_ab*H_ab));
		chblue += madb;
		if(madb > maxchblue)
			maxchblue = madb;
		if(madb < minchblue)
			minchblue = madb;
		maxblueaut = sqrt(reduc*maxchblue);
		minblueaut = sqrt(reduc*minchblue);
				
		chau += (mada+madb);
		nb++;
		//here evaluation of automatic 
		chaut = sqrt(reduc*chau/(nb + nb));
		redaut = sqrt(reduc*chred/nb);
		blueaut = sqrt(reduc*chblue/nb);
		Nb = nb;
	}

}


void ImProcFunctions::WaveletDenoiseAll_info(int levwav, wavelet_decomposition &WaveletCoeffs_a,
											wavelet_decomposition &WaveletCoeffs_b, float **noisevarlum, float **noisevarchrom, float **noisevarhue, int width, int height, float *mad_aa, float *mad_bb, float noisevar_abr, float noisevar_abb, LabImage * noi, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut,int schoice, bool autoch, 
											float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc, float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool perf, bool multiThread ){

	int maxlvl = levwav;
	for (int lvl=0; lvl<maxlvl; lvl++) {

		int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
		int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);

		int skip_ab = WaveletCoeffs_a.level_stride(lvl);

		float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
		float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);

		ShrinkAll_info(WavCoeffs_a, WavCoeffs_b, lvl, Wlvl_ab, Hlvl_ab,
		  skip_ab, noisevarlum,  noisevarchrom, noisevarhue, width, height, mad_aa, mad_bb, noisevar_abr, noisevar_abb, noi, chaut, Nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, 
		 autoch, schoice, lvl, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc, maxchred, maxchblue, minchred, minchblue, nb, chau, chred, chblue, perf, multiThread );

	}
}
	
void ImProcFunctions::RGB_denoise_infoGamCurve(const procparams::DirPyrDenoiseParams & dnparams, bool isRAW, LUTf &gamcurve, float &gam, float &gamthresh, float &gamslope) {
	gam = dnparams.gamma;
	gamthresh = 0.001f;
	if(!isRAW) {//reduce gamma under 1 for Lab mode ==> TIF and JPG
		if(gam <1.9f)
			gam=1.f - (1.9f-gam)/3.f;//minimum gamma 0.7
		else if (gam >= 1.9f && gam <= 3.f)
			gam=(1.4f/1.1f)*gam - 1.41818f;
	}
	gamslope = exp(log((double)gamthresh)/gam)/gamthresh;
	bool perf = (dnparams.dmethod=="RGB");
	if(perf) {
		for (int i=0; i<65536; i++) {
			gamcurve[i] = (Color::gamma((double)i/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)) * 32768.0f;
		}
	}
	else {
		for (int i=0; i<65536; i++) {
			gamcurve[i] = (Color::gamman((double)i/65535.0,gam)) * 32768.0f;
		}
	}
}

void ImProcFunctions::calcautodn_info (float &chaut, float &delta, int Nb, int levaut, float maxmax, float lumema, float chromina, int mode, int lissage, float redyel, float skinc, float nsknc) {

	float reducdelta=1.f;
	if (params->dirpyrDenoise.smethod == "shalbi")
		reducdelta = (float)settings->nrhigh;

	chaut = (chaut*Nb-maxmax)/(Nb-1);//suppress maximum for chaut calcul
	if (redyel > 5000.f || skinc > 1000.f && nsknc < 0.4f  && chromina > 3000.f)
		chaut *= 0.45f;//reduct action in red zone, except skin for high / med chroma
	else if (redyel>12000.f || skinc > 1200.f && nsknc < 0.3f && chromina > 3000.f)
		chaut *= 0.3f;
	
	if (mode==0 || mode==2) {//Preview or Auto multizone
		if (chromina > 10000.f)	chaut *= 0.7f;//decrease action for high chroma  (visible noise)
		else if (chromina > 6000.f)	chaut *= 0.9f;	
		else if (chromina < 3000.f)	chaut *= 1.2f;//increase action in low chroma==> 1.2  /==>2.0 ==> curve CC
		else if (chromina < 2000.f)	chaut *= 1.5f;//increase action in low chroma==> 1.5 / ==>2.7
		
		if (lumema < 2500.f) chaut *= 1.3f;//increase action for low light
		else if (lumema < 5000.f) chaut *= 1.2f;
		else if (lumema > 20000.f) chaut *= 0.9f;//decrease for high light		
	} else if (mode == 1) {//auto ==> less coefficient because interaction
		if (chromina > 10000.f)	chaut *= 0.8f;//decrease action for high chroma  (visible noise)
		else if (chromina > 6000.f)	chaut *= 0.9f;	
		else if (chromina < 3000.f)	chaut *= 1.5f;//increase action in low chroma
		else if (chromina < 2000.f)	chaut *= 2.2f;//increase action in low chroma
		if (lumema < 2500.f) chaut *= 1.2f;//increase action for low light
		else if (lumema < 5000.f) chaut *= 1.1f;
		else if (lumema > 20000.f) chaut *= 0.9f;//decrease for high light		
	}

	if(levaut==0) {//Low denoise
		if(chaut > 300.f)
			chaut = 0.714286f*chaut + 85.71428f;
	}
	
	delta = maxmax-chaut;
	delta *= reducdelta;
	
	if (lissage==1 || lissage==2) {
		if (chaut < 200.f && delta < 200.f ) delta *= 0.95f;
		else if (chaut < 200.f && delta < 400.f ) delta *= 0.5f;
		else if (chaut < 200.f && delta >= 400.f ) delta = 200.f;
		else if (chaut < 400.f && delta < 400.f ) delta *= 0.4f;
		else if (chaut < 400.f && delta >= 400.f ) delta = 120.f;
		else if (chaut < 550.f) delta *= 0.15f;
		else if (chaut < 650.f) delta *= 0.1f;
		else if (chaut >= 650.f) delta *= 0.07f;
		if (mode==0 || mode==2) {//Preview or Auto multizone
			if (chromina < 6000.f) delta *= 1.4f;//increase maxi
			if (lumema < 5000.f) delta *= 1.4f;
		}
		else if (mode==1) {//Auto
			if (chromina < 6000.f) delta *= 1.2f;//increase maxi
			if (lumema < 5000.f) delta *= 1.2f;
		}
	}
	if (lissage==0) {
		if (chaut < 200.f && delta < 200.f ) delta *= 0.95f;
		else if (chaut < 200.f && delta < 400.f ) delta *= 0.7f;
		else if (chaut < 200.f && delta >= 400.f ) delta = 280.f;
		else if (chaut < 400.f && delta < 400.f ) delta *= 0.6f;
		else if (chaut < 400.f && delta >= 400.f ) delta = 200.f;
		else if (chaut < 550.f) delta *= 0.3f;
		else if (chaut < 650.f) delta *= 0.2f;
		else if (chaut >= 650.f) delta *= 0.15f;
		if (mode==0 || mode==2) {//Preview or Auto multizone
			if (chromina < 6000.f) delta *= 1.4f;//increase maxi
			if (lumema < 5000.f) delta *= 1.4f;
		}
		else if (mode==1) {//Auto
			if (chromina < 6000.f) delta *= 1.2f;//increase maxi
			if (lumema < 5000.f) delta *= 1.2f;
		}
	}
	
}
	
SSEFUNCTION void ImProcFunctions::RGB_denoise_info(Imagefloat * src, Imagefloat * provicalc, const bool isRAW, LUTf &gamcurve, float gam, float gamthresh, float gamslope, const procparams::DirPyrDenoiseParams & dnparams, const double expcomp, float &chaut, int &Nb,  float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut, float &nresi, float &highresi, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc, bool multiThread)
{
	if ((settings->leveldnautsimpl==1 && dnparams.Cmethod == "MAN") || (settings->leveldnautsimpl==0 && dnparams.C2method == "MANU")) {
		//nothing to do
		return;
	}

	int hei,wid;
	float** lumcalc;
	float** acalc;
	float** bcalc;
	hei=provicalc->height;
	wid=provicalc->width;
	TMatrix wprofi = iccStore->workingSpaceMatrix (params->icm.working);

	double wpi[3][3] = {
		{wprofi[0][0],wprofi[0][1],wprofi[0][2]},
		{wprofi[1][0],wprofi[1][1],wprofi[1][2]},
		{wprofi[2][0],wprofi[2][1],wprofi[2][2]}
		};
	lumcalc = new float*[hei];
	for (int i=0; i<hei; i++)
		lumcalc[i] = new float[wid];
	acalc = new float*[hei];
	for (int i=0; i<hei; i++)
		acalc[i] = new float[wid];
	bcalc = new float*[hei];
	for (int i=0; i<hei; i++)
		bcalc[i] = new float[wid];

#ifdef _OPENMP
#pragma omp parallel for if(multiThread)
#endif
	for(int ii=0;ii<hei;ii++){
		for(int jj=0;jj<wid;jj++){
			float LLum,AAum,BBum;
			float RL = provicalc->r(ii,jj);
			float GL = provicalc->g(ii,jj);
			float BL = provicalc->b(ii,jj);
			// determine luminance for noisecurve
			float XL,YL,ZL;
			Color::rgbxyz(RL,GL,BL,XL,YL,ZL,wpi);
			Color::XYZ2Lab(XL, YL, ZL, LLum, AAum, BBum);							
			lumcalc[ii][jj]=LLum;
			acalc[ii][jj]=AAum;
			bcalc[ii][jj]=BBum;
		}
	}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	const int imheight=src->height, imwidth=src->width;

	bool perf = (dnparams.dmethod == "RGB");

	const float gain = pow (2.0f, float(expcomp));
	int gamlab = settings->denoiselabgamma;//gamma lab essentialy for Luminance detail
	if(gamlab > 2)
		gamlab = 2;
	
	int tilesize;
	int overlap;
	if(settings->leveldnti == 0) {
		tilesize = 1024;
		overlap = 128;
	}
	if(settings->leveldnti == 1) {
		tilesize = 768;
		overlap = 96;
	}

	int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;

	//always no Tiles
	int kall=0;
	Tile_calc (tilesize, overlap, kall, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	//DCT block data storage
	TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);

	double wp[3][3] = {
		{wprof[0][0],wprof[0][1],wprof[0][2]},
		{wprof[1][0],wprof[1][1],wprof[1][2]},
		{wprof[2][0],wprof[2][1],wprof[2][2]}
	};

	float chau=0.f;
	float chred=0.f;
	float chblue=0.f;
	float maxchred=0.f;
	float maxchblue=0.f;
	float minchred =100000000.f;
	float minchblue=100000000.f;
	int nb=0;
	int comptlevel=0;

	for (int tiletop=0; tiletop<imheight; tiletop+=tileHskip) {
		for (int tileleft=0; tileleft<imwidth; tileleft+=tileWskip) {

			int tileright = MIN(imwidth,tileleft+tilewidth);
			int tilebottom = MIN(imheight,tiletop+tileheight);
			int width  = tileright-tileleft;
			int height = tilebottom-tiletop;
			LabImage * labdn = new LabImage(width,height);
			float** noisevarlum = new float*[height];
			for (int i=0; i<height; i++)
				noisevarlum[i] = new float[width];
						
			float* mad_aa = new float [height*width];
			float* mad_bb = new float [height*width];
			float** noisevarchrom = new float*[height];
			for (int i=0; i<height; i++)
				noisevarchrom[i] = new float[width];
			float** noisevarhue = new float*[height];
			for (int i=0; i<height; i++)
				noisevarhue[i] = new float[width];
			
			// init luma noisevarL
			float noisevarab_b, noisevarab_r;
			
			float realred, realblue;
			float interm_med =(float) dnparams.chroma/10.0;
			float intermred, intermblue;
			if(dnparams.redchro > 0.)
				intermred = (dnparams.redchro/10.);
			else
				intermred = (float) dnparams.redchro/7.0;//increase slower than linear for more sensit
			if(dnparams.bluechro > 0.)
				intermblue =(dnparams.bluechro/10.);
			else
				intermblue = (float) dnparams.bluechro/7.0;//increase slower than linear for more sensit
			realred = interm_med + intermred;
			if (realred < 0.f)
				realred = 0.001f;
			realblue = interm_med + intermblue;
			if (realblue < 0.f)
				realblue = 0.001f;
			//TODO: implement using AlignedBufferMP
			//fill tile from image; convert RGB to "luma/chroma"
			float noiseluma=(float) dnparams.luma;
			noiseluma += 1.f;
			float noisevarL = SQR((noiseluma/125.f)*(1.f+noiseluma/25.f));
			if (isRAW) {//image is raw; use channel differences for chroma channels
#ifdef _OPENMP
#pragma omp parallel for if(multiThread)
#endif
					for (int i=tiletop; i<tilebottom; i++) {
						int i1 = i - tiletop;
#ifdef __SSE2__
						__m128 aNv,bNv;
						__m128 c100v = _mm_set1_ps(100.f);
						int j;
						for (j=tileleft; j<tileright-3; j+=4) {
							int j1 = j - tileleft;
							aNv = LVFU(acalc[i][j]);
							bNv = LVFU(bcalc[i][j]);
							_mm_storeu_ps(&noisevarhue[i1][j1], xatan2f(bNv,aNv));
							_mm_storeu_ps(&noisevarchrom[i1][j1], _mm_max_ps(c100v,_mm_sqrt_ps(SQRV(aNv)+SQRV(bNv))));
						}
						for (; j<tileright; j++) {
							int j1 = j - tileleft;
							float aN = acalc[i][j];
							float bN = bcalc[i][j];
							float cN = sqrtf(SQR(aN)+SQR(bN));
							noisevarhue[i1][j1] = xatan2f(bN,aN);
							if(cN < 100.f)
								cN=100.f;//avoid divided by zero
							noisevarchrom[i1][j1] = cN;
						}
#else
						for (int j=tileleft; j<tileright; j++) {
							int j1 = j - tileleft;
							float aN = acalc[i][j];
							float bN = bcalc[i][j];
							float cN = sqrtf(SQR(aN)+SQR(bN));
							float hN = xatan2f(bN,aN);
							if(cN < 100.f)
								cN = 100.f;//avoid divided by zero
							noisevarchrom[i1][j1] = cN;
							noisevarhue[i1][j1] = hN;
						}
#endif
					}
#ifdef _OPENMP
#pragma omp parallel for if(multiThread)
#endif
					for (int i=tiletop; i<tilebottom; i++) {
						int i1 = i - tiletop;
						for (int j=tileleft; j<tileright; j++) {
							int j1 = j - tileleft;
							float Llum = lumcalc[i][j];
							Llum = Llum < 2.f ? 2.f : Llum; //avoid divided by zero ?
							Llum = Llum > 32768.f ? 32768.f : Llum; // not strictly necessary							
							noisevarlum[i1][j1]= Llum;
						}
					}
				if (!perf){//lab mode, modification Jacques feb 2013 and july 2014

#ifdef _OPENMP
#pragma omp parallel for if(multiThread)
#endif
				for (int i=tiletop; i<tilebottom; i++) {
					int i1 = i - tiletop;
					for (int j=tileleft; j<tileright; j++) {
						int j1 = j - tileleft;
						float R_ = gain*src->r(i,j);
						float G_ = gain*src->g(i,j);
						float B_ = gain*src->b(i,j);

						//finally I opted fot gamma55 and with options we can change
						if (gamlab == 0) {// options 12/2013
							R_ = Color::igammatab_26_11[R_];
							G_ = Color::igammatab_26_11[G_];
							B_ = Color::igammatab_26_11[B_];
						}
						else if (gamlab == 1) {
							//other new gamma 4 5 
							R_ = Color::igammatab_4[R_];
							G_ = Color::igammatab_4[G_];
							B_ = Color::igammatab_4[B_];
						}
						else if (gamlab == 2) {							
							//new gamma 5.5 10 better for detail luminance..it is a compromise...which depends on the image (distribution BL, ML, HL ...)
							R_ = Color::igammatab_55[R_];
							G_ = Color::igammatab_55[G_];
							B_ = Color::igammatab_55[B_];
						}
						//apply gamma noise	standard (slider)
						R_ = R_<65535.0f ? gamcurve[R_] : (Color::gamman((double)R_/65535.0, gam)*32768.0f);
						G_ = G_<65535.0f ? gamcurve[G_] : (Color::gamman((double)G_/65535.0, gam)*32768.0f);
						B_ = B_<65535.0f ? gamcurve[B_] : (Color::gamman((double)B_/65535.0, gam)*32768.0f);
						//true conversion xyz=>Lab
						float X,Y,Z;
						Color::rgbxyz(R_,G_,B_,X,Y,Z,wp);

						//convert to Lab
						float L,a,b;
						Color::XYZ2Lab(X, Y, Z, L, a, b);
						
						labdn->a[i1][j1] = a;
						labdn->b[i1][j1] = b;
					}
				}
				}
				else {//RGB mode
				
				for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
					int i1 = i - tiletop;
					for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
						int j1 = j - tileleft;

						float X = gain*src->r(i,j);
						float Y = gain*src->g(i,j);
						float Z = gain*src->b(i,j);
						
						X = X<65535.0f ? gamcurve[X] : (Color::gamma((double)X/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						Y = Y<65535.0f ? gamcurve[Y] : (Color::gamma((double)Y/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						Z = Z<65535.0f ? gamcurve[Z] : (Color::gamma((double)Z/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);

						labdn->a[i1][j1] = (X-Y);
						labdn->b[i1][j1] = (Y-Z);
					}
				}
				}
				
			} else {//image is not raw; use Lab parametrization
				for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
					int i1 = i - tiletop;
					for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
						int j1 = j - tileleft;
						float L,a,b;
						float rLum=src->r(i,j) ;//for luminance denoise curve
						float gLum=src->g(i,j) ;
						float bLum=src->b(i,j) ;
						
						//use gamma sRGB, not good if TIF (JPG) Output profil not with gamma sRGB  (eg : gamma =1.0, or 1.8...)
						//very difficult to solve !
						// solution ==> save TIF with gamma sRGB and re open
						float rtmp = Color::igammatab_srgb[ src->r(i,j) ];
						float gtmp = Color::igammatab_srgb[ src->g(i,j) ];
						float btmp = Color::igammatab_srgb[ src->b(i,j) ];
						//modification Jacques feb 2013
						// gamma slider different from raw
						rtmp = rtmp<65535.0f ? gamcurve[rtmp] : (Color::gamman((double)rtmp/65535.0, gam)*32768.0f);
						gtmp = gtmp<65535.0f ? gamcurve[gtmp] : (Color::gamman((double)gtmp/65535.0, gam)*32768.0f);
						btmp = btmp<65535.0f ? gamcurve[btmp] : (Color::gamman((double)btmp/65535.0, gam)*32768.0f);
						
						float X,Y,Z;
						Color::rgbxyz(rtmp,gtmp,btmp,X,Y,Z,wp);

						//convert Lab
						Color::XYZ2Lab(X, Y, Z, L, a, b);
						float Llum,alum,blum;
						float XL,YL,ZL;
						Color::rgbxyz(rLum,gLum,bLum,XL,YL,ZL,wp);
						Color::XYZ2Lab(XL, YL, ZL, Llum, alum, blum);
						
						float kN=Llum;
						float epsi=0.01f;

						if(kN<2.f) kN=2.f;
						if(kN>32768.f) kN=32768.f;								
						noisevarlum[i1][j1]=kN;
						float aN=alum;
						float bN=blum;
						float hN=xatan2f(bN,aN);
						float cN=sqrt(SQR(aN)+SQR(bN));
						if(cN < 100.f) cN=100.f;//avoid divided by zero
						noisevarchrom[i1][j1]=cN;
						noisevarhue[i1][j1]=hN;

						labdn->a[i1][j1] = a;
						labdn->b[i1][j1] = b;
					}
				}
			}
			int datalen = labdn->W * labdn->H;

			//now perform basic wavelet denoise
			//last two arguments of wavelet decomposition are max number of wavelet decomposition levels;
			//and whether to subsample the image after wavelet filtering.  Subsampling is coded as
			//binary 1 or 0 for each level, eg subsampling = 0 means no subsampling, 1 means subsample
			//the first level only, 7 means subsample the first three levels, etc.
			noisevarab_r = SQR(realred) + 0.01f;
			noisevarab_b = SQR(realblue) + 0.01f;
			
			wavelet_decomposition* adecomp;
			wavelet_decomposition* bdecomp;

			int schoice = 0;//shrink method
			if (dnparams.smethod == "shalbi")
				schoice=2;

			const int levwav=5;
#ifdef _OPENMP
#pragma omp parallel sections if(multiThread)
#endif
{
#ifdef _OPENMP
#pragma omp section
#endif
{
			adecomp = new wavelet_decomposition (labdn->data+datalen, labdn->W, labdn->H, levwav, 1 );
}
#ifdef _OPENMP
#pragma omp section
#endif
{
			bdecomp = new wavelet_decomposition (labdn->data+2*datalen, labdn->W, labdn->H, levwav, 1 );
}
}
			bool autoch = dnparams.autochroma;
			if(comptlevel==0) WaveletDenoiseAll_info(levwav, *adecomp, *bdecomp, noisevarlum, noisevarchrom, noisevarhue, width, height, mad_aa, mad_bb, noisevarab_r, noisevarab_b, labdn, chaut, Nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, schoice, autoch, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc,maxchred, maxchblue, minchred, minchblue, nb,chau ,chred, chblue, perf, multiThread);//enhance mode
			comptlevel+=1;
			float chresid,chmaxredresid,chmaxblueresid,chresidred, chresidblue;	
			nresi=chresid;
			highresi=chresid + 0.66f*(max(chmaxredresid,chmaxblueresid) - chresid);//evaluate sigma
			delete adecomp;
			delete bdecomp;
			delete labdn;
			for (int i=0; i<height; i++)
				delete [] noisevarlum[i];
			delete [] noisevarlum;
			for (int i=0; i<height; i++)
				delete [] noisevarchrom[i];
			delete [] noisevarchrom;
			for (int i=0; i<height; i++)
				delete [] noisevarhue[i];
			delete [] noisevarhue;
			
			delete [] mad_aa;
			delete [] mad_bb;

		}//end of tile row
	}//end of tile loop

	for (int i=0; i<hei; i++)
		delete [] lumcalc[i];
	delete [] lumcalc;
	for (int i=0; i<hei; i++)
		delete [] acalc[i];
	delete [] acalc;
	for (int i=0; i<hei; i++)
		delete [] bcalc[i];
	delete [] bcalc;

#undef TS
#undef fTS
#undef offset
#undef epsilon

}//end of main RGB_denoise

}
