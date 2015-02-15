////////////////////////////////////////////////////////////////
//
//
//
//
//  code dated: December , 2014
//
//	Ipwaveletcc is free software: you can redistribute it and/or modify
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
// *  2014 Jacques Desmis <jdesmis@gmail.com>
// *  2014 Ingo Weyrich <heckflosse@i-weyrich.de>

//
////////////////////////////////////////////////////////////////



#include <math.h>
#include "../rtgui/threadutils.h"

#include "rtengine.h"
#include "improcfun.h"
#include "LUT.h"
#include "array2D.h"
#include "boxblur.h"
#include "rt_math.h"
#include "mytime.h"
#include "sleef.c"
#include "opthelper.h"
#include "StopWatch.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cplx_wavelet_dec.h"

#define TS 64		// Tile size
#define offset 25	// shift between tiles
#define fTS ((TS/2+1))	// second dimension of Fourier tiles
#define blkrad 1	// radius of block averaging

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

#define epsilon 0.001f/(TS*TS) //tolerance


namespace rtengine {

extern const Settings* settings;

struct cont_params {
  float mul[10];
  int chrom;
  int chro;
  int unif;
  float th;
  float thH;
  float conres;
  float conresH;
  float chrores;
  float sky;
  float b_l,t_l,b_r,t_r;
  float b_ly,t_ly,b_ry,t_ry;
  float b_lsl,t_lsl,b_rsl,t_rsl;
  float b_lhl,t_lhl,b_rhl,t_rhl;
  
  float b_lpast,t_lpast,b_rpast,t_rpast;
  float b_lsat,t_lsat,b_rsat,t_rsat;
  
  int numlevH, numlevS;	
  float mulC[9];
  float mulopaRG[9];
  float mulopaBY[9];
  bool curv;
  bool opaBY;
  bool opaRG;
  int CHmet;
  bool HSmet;
  bool avoi;
};

int wavNestedLevels = 1;

	void ImProcFunctions::ip_wavelet(LabImage * lab, LabImage * dst, int kall, const procparams::WaveletParams & waparams, const WavCurve & wavCLVCcurve, const WavOpacityCurveRG & waOpacityCurveRG, const WavOpacityCurveBY & waOpacityCurveBY, int skip)
	{
	MyTime t1e,t2e;
	t1e.set();
	
#ifdef _DEBUG
	// init variables to display Munsell corrections
	MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif
	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};
	const short int imheight=lab->H, imwidth=lab->W;
	struct cont_params cp;
	cp.avoi=false;
	if(params->wavelet.avoid)	cp.avoi=true;			
	
	int N=imheight*imwidth;
	int maxmul=params->wavelet.thres;
	cp.curv=false;
	cp.opaRG=false;
	cp.opaBY=false;
	cp.CHmet=0;	
	cp.HSmet=false;
	if(params->wavelet.CHmethod=="with")	cp.CHmet=1;			
	if(params->wavelet.CHmethod=="link")	cp.CHmet=2;			
	if(params->wavelet.HSmethod=="with")  cp.HSmet=true;					

	if(wavCLVCcurve) cp.curv=true; 
	
	if(cp.curv) {//convert curve in discret values
		cp.mulC[0]=200.f*(wavCLVCcurve[0]-0.5f);
		cp.mulC[1]=200.f*(wavCLVCcurve[62]-0.5f);
		cp.mulC[2]=200.f*(wavCLVCcurve[125]-0.5f);
		cp.mulC[3]=200.f*(wavCLVCcurve[187]-0.5f);
		cp.mulC[4]=200.f*(wavCLVCcurve[250]-0.5f);
		cp.mulC[5]=200.f*(wavCLVCcurve[312]-0.5f);
		cp.mulC[6]=200.f*(wavCLVCcurve[375]-0.5f);
		cp.mulC[7]=200.f*(wavCLVCcurve[438]-0.5f);
		cp.mulC[8]=200.f*(wavCLVCcurve[500]-0.5f);
	}
	else {
         for(int level=0;level<9;level++)
			cp.mulC[level] = 0.f;
    }
	if(waOpacityCurveRG) cp.opaRG=true;
		
	if(cp.opaRG) {
		cp.mulopaRG[0]=200.f*(waOpacityCurveRG[0]-0.5f);
		cp.mulopaRG[1]=200.f*(waOpacityCurveRG[62]-0.5f);
		cp.mulopaRG[2]=200.f*(waOpacityCurveRG[125]-0.5f);
		cp.mulopaRG[3]=200.f*(waOpacityCurveRG[187]-0.5f);
		cp.mulopaRG[4]=200.f*(waOpacityCurveRG[250]-0.5f);
		cp.mulopaRG[5]=200.f*(waOpacityCurveRG[312]-0.5f);
		cp.mulopaRG[6]=200.f*(waOpacityCurveRG[375]-0.5f);
		cp.mulopaRG[7]=200.f*(waOpacityCurveRG[438]-0.5f);
		cp.mulopaRG[8]=200.f*(waOpacityCurveRG[500]-0.5f);
	}
	else {
         for(int level=0;level<9;level++)
			cp.mulopaRG[level] = 0.f;
    }
	
	if(waOpacityCurveBY) cp.opaBY=true;
	if(cp.opaBY) {
		cp.mulopaBY[0]=200.f*(waOpacityCurveBY[0]-0.5f);
		cp.mulopaBY[1]=200.f*(waOpacityCurveBY[62]-0.5f);
		cp.mulopaBY[2]=200.f*(waOpacityCurveBY[125]-0.5f);
		cp.mulopaBY[3]=200.f*(waOpacityCurveBY[187]-0.5f);
		cp.mulopaBY[4]=200.f*(waOpacityCurveBY[250]-0.5f);
		cp.mulopaBY[5]=200.f*(waOpacityCurveBY[312]-0.5f);
		cp.mulopaBY[6]=200.f*(waOpacityCurveBY[375]-0.5f);
		cp.mulopaBY[7]=200.f*(waOpacityCurveBY[438]-0.5f);
		cp.mulopaBY[8]=200.f*(waOpacityCurveBY[500]-0.5f);
	}
	else {
         for(int level=0;level<9;level++)
			cp.mulopaBY[level] = 0.f;
    }
	
	for(int m=0;m<maxmul;m++)
        cp.mul[m]=waparams.c[m];
        cp.mul[9]=(float) waparams.sup;
	
		cp.chro=waparams.chro;
		cp.chrom=waparams.chroma;
		cp.unif=waparams.unif;
		cp.conres=waparams.rescon;
		cp.conresH=waparams.resconH;
		cp.chrores=waparams.reschro;
		cp.th=float(waparams.thr);
		cp.thH=float(waparams.thrH);
		cp.sky=waparams.sky;
		//skin
		cp.b_l = static_cast<float>(params->wavelet.hueskin.value[0]) / 100.0f;
		cp.t_l = static_cast<float>(params->wavelet.hueskin.value[1]) / 100.0f;
		cp.b_r = static_cast<float>(params->wavelet.hueskin.value[2]) / 100.0f;
		cp.t_r = static_cast<float>(params->wavelet.hueskin.value[3]) / 100.0f;
		
		cp.b_ly = static_cast<float>(params->wavelet.hueskin2.value[0]) / 100.0f;
		cp.t_ly = static_cast<float>(params->wavelet.hueskin2.value[1]) / 100.0f;
		cp.b_ry = static_cast<float>(params->wavelet.hueskin2.value[2]) / 100.0f;
		cp.t_ry = static_cast<float>(params->wavelet.hueskin2.value[3]) / 100.0f;
		cp.numlevH=params->wavelet.threshold;
		
		cp.numlevH=params->wavelet.threshold;
		//shadows
		cp.b_lsl = static_cast<float>(params->wavelet.bllev.value[0]);
		cp.t_lsl = static_cast<float>(params->wavelet.bllev.value[1]);
		cp.b_rsl = static_cast<float>(params->wavelet.bllev.value[2]);
		cp.t_rsl = static_cast<float>(params->wavelet.bllev.value[3]);
		cp.numlevS=params->wavelet.threshold2;
		int maxlevS=9-cp.numlevH;
		cp.numlevS = MIN(cp.numlevS,maxlevS);
		//highlight
		cp.b_lhl = static_cast<float>(params->wavelet.hllev.value[0]);
		cp.t_lhl = static_cast<float>(params->wavelet.hllev.value[1]);
		cp.b_rhl = static_cast<float>(params->wavelet.hllev.value[2]);
		cp.t_rhl = static_cast<float>(params->wavelet.hllev.value[3]);
		//printf("H=%d S=%d\n",cp.numlevH,cp.numlevS);
		//pastel
		cp.b_lpast = static_cast<float>(params->wavelet.pastlev.value[0]);
		cp.t_lpast = static_cast<float>(params->wavelet.pastlev.value[1]);
		cp.b_rpast = static_cast<float>(params->wavelet.pastlev.value[2]);
		cp.t_rpast = static_cast<float>(params->wavelet.pastlev.value[3]);
		//saturated
		cp.b_lsat = static_cast<float>(params->wavelet.satlev.value[0]);
		cp.t_lsat = static_cast<float>(params->wavelet.satlev.value[1]);
		cp.b_rsat = static_cast<float>(params->wavelet.satlev.value[2]);
		cp.t_rsat = static_cast<float>(params->wavelet.satlev.value[3]);
		
		
		int minwin=min(imwidth,imheight);
		int maxlevelcrop=9;
		if(cp.mul[9]!=0)
			maxlevelcrop=10;
		// adap maximum level wavelet to size of crop
		if(minwin*skip < 1024) maxlevelcrop = 9;//sampling wavelet 512
		if(minwin*skip < 512) maxlevelcrop = 8;//sampling wavelet 256
		if(minwin*skip < 256) maxlevelcrop = 7;//sampling 128
		if(minwin*skip < 128) maxlevelcrop = 6;
		if(minwin < 64) maxlevelcrop = 5;
	//	printf("minwin=%d maxcrop=%d\n",minwin, maxlevelcrop);
		
		int levwav=params->wavelet.thres;
		if(levwav==9 && cp.mul[9]!=0) levwav=10;
		levwav=min(maxlevelcrop,levwav);
		// determine number of levels to process.
	//	for(levwav=min(maxlevelcrop,levwav);levwav>0;levwav--)
	//		if(cp.mul[levwav-1]!=0.f  || cp.curv)
		//	if(cp.mul[levwav-1]!=0.f)
	//			break;
	// I suppress this fonctionality ==> crash for level < 3
		if(levwav<1)
			return; // nothing to do  

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// begin tile processing of image

		//output buffer
		int realtile;
		if(params->wavelet.Tilesmethod=="big") realtile=22;
		if(params->wavelet.Tilesmethod=="lit") realtile=12;
		
		int tilesize;
		int overlap;
		tilesize = 1024;
		overlap = 128;
		//tilesize=128*params->wavelet.tiles;
		tilesize=128*realtile;
		//overlap=(int) tilesize*params->wavelet.overl;
		overlap=(int) tilesize*0.125f;
		//	printf("overl=%d\n",overlap);
		int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;
		if(params->wavelet.Tilesmethod=="full") kall=0;
		Tile_calc (tilesize, overlap, kall, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);

		const int numtiles = numtiles_W*numtiles_H;
		LabImage * dsttmp;
		if(numtiles == 1) {
			dsttmp = dst;
		} else {
			dsttmp = new LabImage(imwidth,imheight);
			for (int n=0; n<3*imwidth*imheight; n++) dsttmp->data[n] = 0;
		}
		//now we have tile dimensions, overlaps
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		int minsizetile=min(tilewidth, tileheight);
		int maxlev2=10;
		if(minsizetile < 1024 && levwav==10) maxlev2 = 9;
		if(minsizetile < 512) maxlev2 = 8;
		if(minsizetile < 256) maxlev2 = 7;
		if(minsizetile < 128) maxlev2 = 6;
		levwav=min(maxlev2,levwav);
		
		//printf("levwav = %d\n",levwav);

		int numthreads = 1;
		int maxnumberofthreadsforwavelet =0;
		//reduce memory for big tile size
		if(kall!=0) {
			if(realtile <= 22) maxnumberofthreadsforwavelet=2;
			if(realtile <= 20) maxnumberofthreadsforwavelet=3;
			if(realtile <= 18) maxnumberofthreadsforwavelet=4;
			if(realtile <= 16) maxnumberofthreadsforwavelet=6;
			if(realtile <= 14) maxnumberofthreadsforwavelet=8;
			//printf("maxNRT=%d\n",maxnumberofthreadsforwavelet);
			if((maxnumberofthreadsforwavelet==6 || maxnumberofthreadsforwavelet==8)  && levwav==10) maxnumberofthreadsforwavelet-=2;
			if(levwav <=7 && maxnumberofthreadsforwavelet ==8) maxnumberofthreadsforwavelet=0;
		}
		//printf("maxthre=%d\n",maxnumberofthreadsforwavelet);

#ifdef _OPENMP
	// Calculate number of tiles. If less than omp_get_max_threads(), then limit num_threads to number of tiles
		if( options.rgbDenoiseThreadLimit>0)
			maxnumberofthreadsforwavelet = min(max(options.rgbDenoiseThreadLimit / 2, 1), maxnumberofthreadsforwavelet);	
	
		numthreads = MIN(numtiles,omp_get_max_threads());
		if(maxnumberofthreadsforwavelet > 0)
			numthreads = MIN(numthreads,maxnumberofthreadsforwavelet);
		wavNestedLevels = omp_get_max_threads() / numthreads;
		bool oldNested = omp_get_nested();
		if(wavNestedLevels < 2)
			wavNestedLevels = 1;
		else
			omp_set_nested(true);
		if(maxnumberofthreadsforwavelet > 0)
			while(wavNestedLevels*numthreads > maxnumberofthreadsforwavelet)
				wavNestedLevels--;
		if(settings->verbose)
			printf("Ip Wavelet uses %d main thread(s) and up to %d nested thread(s) for each main thread\n",numthreads,wavNestedLevels);


#endif			

#pragma omp parallel num_threads(numthreads)
	{
	   float *mean = new float [9];
	   float *meanN = new float [9];
	   float *sigma = new float [9];
	   float *sigmaN = new float [9];
		
		float** varhue = new float*[tileheight];
		for (int i=0; i<tileheight; i++)
			varhue[i] = new float[tilewidth];
		float** varchro = new float*[tileheight];
		for (int i=0; i<tileheight; i++)
			varchro[i] = new float[tilewidth];

#ifdef _OPENMP
#pragma omp for schedule(dynamic) collapse(2)
#endif	
		for (int tiletop=0; tiletop<imheight; tiletop+=tileHskip) {
			for (int tileleft=0; tileleft<imwidth ; tileleft+=tileWskip) {
				int tileright = MIN(imwidth,tileleft+tilewidth);
				int tilebottom = MIN(imheight,tiletop+tileheight);
				int width  = tileright-tileleft;
				int height = tilebottom-tiletop;
				LabImage * labco;
				if(numtiles == 1 && !cp.avoi) // untiled processing and no 'Avoid Colour Shift' => we can use output buffer for labco
					labco = dst;
				else
					labco = new LabImage(width,height);
				
#ifdef _OPENMP
#pragma omp parallel for num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
				
					for (int i=tiletop; i<tilebottom; i++) {
						int i1 = i - tiletop;
						for (int j=tileleft; j<tileright; j++) {
							int j1 = j - tileleft;
							float L=lab->L[i][j];
							float a=lab->a[i][j];
							float b=lab->b[i][j];
							labco->L[i1][j1] = L;
							labco->a[i1][j1] = a;
							labco->b[i1][j1] = b;
							varhue[i1][j1]=xatan2f(b,a);
							varchro[i1][j1]=(sqrt(a*a+b*b))/327.68f;	
						}
					}
				//to avoid artifacts in blue sky
				if(params->wavelet.median) {	
				float** tmL;
				int wid=labco->W;
				int hei=labco->H;
				int borderL = 1;
				tmL = new float*[hei];
				for (int i=0; i<hei; i++)
					tmL[i] = new float[wid];
				for(int i = borderL; i < hei-borderL; i++ ) {
						for(int j = borderL; j < wid-borderL; j++) {
							tmL[i][j] = labco->L[i][j];
					}
				}
				
#ifdef _OPENMP
#pragma omp parallel for num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
				for (int i=1; i<hei-1; i++) {
						float pp[9],temp;
							for (int j=1; j<wid-1; j++) {
							if((varhue[i][j] < -1.3f && varhue[i][j] > - 2.6f)  && (varchro[i][j] > 15.f && varchro[i][j] < 55.f) && labco->L[i][j] > 5000.f) //blue sky + med3x3  ==> after for more effect use denoise
							med3(labco->L[i][j] ,labco->L[i-1][j], labco->L[i+1][j] ,labco->L[i][j+1],labco->L[i][j-1], labco->L[i-1][j-1],labco->L[i-1][j+1],labco->L[i+1][j-1],labco->L[i+1][j+1],tmL[i][j]);//3x3
							}
				}
				for(int i = borderL; i < hei-borderL; i++ ) {
					for(int j = borderL; j < wid-borderL; j++) {
						labco->L[i][j] = tmL[i][j];
					}
				}

			for (int i=0; i<hei; i++)
				delete [] tmL[i];
			delete [] tmL;
			// end blue sky	
		}
				if(params->wavelet.skinprotect != 0.0  || (cp.curv && cp.CHmet!=2)) // reduce the arrays to get faster access in following processing
					for (int i=0; i<(tileheight)/2; i++) {
						for (int j=0; j<(tilewidth)/2; j++) {
							varhue[i][j]=varhue[i*2][j*2];
						}
					}
				int datalen = labco->W * labco->H;

                wavelet_decomposition* Ldecomp = new wavelet_decomposition (labco->data, labco->W, labco->H, levwav, 1, skip, max(1,wavNestedLevels) );
                if(!Ldecomp->memoryAllocationFailed) {
					WaveletcontAllL(labco, varhue, varchro, *Ldecomp, cp);
					Ldecomp->reconstruct(labco->data);
                }
				delete Ldecomp;

				wavelet_decomposition* adecomp = new wavelet_decomposition (labco->data+datalen, labco->W, labco->H,levwav, 1, skip, max(1,wavNestedLevels) );
                if(!adecomp->memoryAllocationFailed) {
					WaveletcontAllAB(labco, varhue, varchro, *adecomp, cp, true);
					adecomp->reconstruct(labco->data+datalen);
                }
				delete adecomp;

                wavelet_decomposition* bdecomp = new wavelet_decomposition (labco->data+2*datalen, labco->W, labco->H, levwav, 1, skip, max(1,wavNestedLevels) );
                if(!bdecomp->memoryAllocationFailed) {
					WaveletcontAllAB(labco, varhue, varchro, *bdecomp, cp, false);
					bdecomp->reconstruct(labco->data+2*datalen);
                }
				delete bdecomp;

				if(numtiles > 1 || (numtiles == 1 && cp.avoi)) {
					//calculate mask for feathering output tile overlaps
					float Vmask[height+overlap] ALIGNED16;
					float Hmask[width+overlap] ALIGNED16;

					if(numtiles > 1) {
						for (int i=0; i<height; i++) {
							Vmask[i] = 1;
						}
						for (int j=0; j<width; j++) {
							Hmask[j] = 1;
						}
						for (int i=0; i<overlap; i++) {
							float mask = SQR(sin((M_PI*i)/(2*overlap)));
							if (tiletop>0) Vmask[i] = mask;
							if (tilebottom<imheight) Vmask[height-i] = mask;
							if (tileleft>0) Hmask[i] = mask;
							if (tileright<imwidth) Hmask[width-i] = mask;
						}
					}
				
#ifdef _OPENMP
#pragma omp parallel for num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif				
					for (int i=tiletop; i<tilebottom; i++){
							int i1 = i-tiletop;
							float X,Y,Z,L,a,b;

							for (int j=tileleft; j<tileright; j++) {
								int j1=j-tileleft;
								L = labco->L[i1][j1];
								a = labco->a[i1][j1];
								b = labco->b[i1][j1];
								if(cp.avoi){//Gamut and Munsell
									float HH=xatan2f(b,a);
									float Chprov1=sqrt(SQR(a/327.68f) + SQR(b/327.68f));
									float Lprov1=L/327.68f;
									float Lprov2 = lab->L[i][j]/327.68f;
									float memChprov=varchro[i1][j1];
									bool highlight = params->toneCurve.hrenabled;
									float R,G,B;
	#ifdef _DEBUG
									bool neg=false;
									bool more_rgb=false;
									Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
	#else
									Color::gamutLchonly(HH,Lprov1,Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
	#endif
									L=Lprov1*327.68f;
									float2 sincosv = xsincosf(HH);
									
									a=327.68f*Chprov1*sincosv.y;//gamut
									b=327.68f*Chprov1*sincosv.x;//gamut
									{
									float correctionHue=0.0f; // Munsell's correction
									float correctlum=0.0f;
									Lprov1=L/327.68f;
									float Chprov=sqrt(SQR(a/327.68f)+ SQR(b/327.68f));
	#ifdef _DEBUG
									Color::AllMunsellLch(true, Lprov1,Lprov2,HH,Chprov,memChprov,correctionHue,correctlum, MunsDebugInfo);
	#else
									Color::AllMunsellLch(true, Lprov1,Lprov2,HH,Chprov,memChprov,correctionHue,correctlum);
	#endif

									if(fabs(correctionHue) < 0.015f) HH+=correctlum;	// correct only if correct Munsell chroma very little.
									float2 sincosval = xsincosf(HH+correctionHue);
				
									a=327.68f*Chprov*sincosval.y;// apply Munsell
									b=327.68f*Chprov*sincosval.x;//aply Munsell
									}
								}
								if(numtiles > 1) {
									float factor = Vmask[i1]*Hmask[j1];
									dsttmp->L[i][j]+= factor*L;
									dsttmp->a[i][j]+= factor*a;
									dsttmp->b[i][j]+= factor*b;
								} else {
									dsttmp->L[i][j] = L;
									dsttmp->a[i][j] = a;
									dsttmp->b[i][j] = b;
									
								}
							}
					}
				}
				if(numtiles>1 || cp.avoi)
					delete labco;
			}
		}
		for (int i=0; i<tileheight; i++)
			delete [] varhue[i];
		delete [] varhue;
		for (int i=0; i<tileheight; i++)
		delete [] varchro[i];
		delete [] varchro;
	
		delete [] mean;
		delete [] meanN;
		delete [] sigma;
		delete [] sigmaN;
		
	}
#ifdef _OPENMP
omp_set_nested(oldNested);
#endif	
		if(numtiles > 1) {
			dst->CopyFrom(dsttmp);
			delete dsttmp;
		}
		
	if (settings->verbose) {
		t2e.set();
		printf("Wavelet performed in %d usec:\n", t2e.etime(t1e));
	}
		
}//end o




#undef TS
#undef fTS
#undef offset
#undef epsilon

	void ImProcFunctions::Aver( float *  RESTRICT DataList, int datalen, float &averagePlus, float &averageNeg, float &max, float &min) {
		//find absolute mean
		int averaP=0, averaN=0, count=0, countP=0, countN=0;
		max=0.f;min=0.f;
		averagePlus=0.f;averageNeg=0.f;
		while (count<datalen) {
			if(DataList[count] >= 0.f) {averaP += abs((int)DataList[count]);	
			if(abs((int)DataList[count])> max) max=abs((int)DataList[count]);
			countP++;
			}
			if(DataList[count] < 0.f) {averaN += abs((int)DataList[count]);	
			if(abs((int)DataList[count])> min) min=abs((int)DataList[count]);
			countN++;
			}
			
			count++;
		}	
		averagePlus=averaP/countP;
		averageNeg=averaN/countN;
		
	}
	
	
	void ImProcFunctions::Sigma( float *  RESTRICT DataList, int datalen, float averagePlus, float averageNeg, float &sigmaPlus, float &sigmaNeg) {
	    int count=0, countP=0, countN=0;
		float variP=0.f,variN=0.f;
		while (count<datalen) {
			if(DataList[count] >= 0.f) {variP += SQR(DataList[count] - averagePlus);	
			countP++;
			}
			else if(DataList[count] < 0.f) {variN += SQR(DataList[count] - averageNeg);	
			countN++;
			}
		count++;
		}
		sigmaPlus=sqrt(variP/countP);
		sigmaNeg=sqrt(variN/countN);

	}
	
	void ImProcFunctions::Evaluate(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a,
											wavelet_decomposition &WaveletCoeffs_b, float *av_LL, float *av_aa, float *av_bb,struct cont_params cp, int ind, float *mean, float *meanN, float *sigma, float *sigmaN){
		int maxlvl = WaveletCoeffs_L.maxlevel();
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

			Eval (WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, lvl, cp, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,skip_L, skip_ab, av_LL, av_aa, av_bb, ind, mean, meanN, sigma, sigmaN);		
		}
											
	}										
	void ImProcFunctions::Eval (float ** WavCoeffs_L, float ** WavCoeffs_a, float ** WavCoeffs_b, int level,struct cont_params cp,
									int W_L, int H_L, int W_ab, int H_ab,int skip_L, int skip_ab, float * av_LL, float * av_aa, float * av_bb, int ind, float *mean, float *meanN, float *sigma, float *sigmaN)
	{
		const float eps = 0.01f;
		
		float ava[4], avb[4], avLP[4], avLN[4];
		float maxL[4], minL[4], maxa[4], maxb[4];
		float sigP[4], sigN[4];
		float AvL,AvN,SL,SN;
		
		float thr= params->wavelet.thres;
		for (int dir=1; dir<4; dir++) {
			{
				float averagePlus=0.f,averageNeg=0.f, max, min;
			//	Aver(WavCoeffs_L[dir], W_L*H_L,  averagePlus, averageNeg, max, min);
				Aver(WavCoeffs_b[dir], W_L*H_L,  averagePlus, averageNeg, max, min);
				avLP[dir] = fabs(averagePlus);
				avLN[dir] = -fabs(averageNeg);
				maxL[dir] = max;
				minL[dir] = -min;
				float sigmaPlus, sigmaNeg;
				Sigma(WavCoeffs_b[dir], W_L*H_L, avLP[dir], -avLN[dir], sigmaPlus, sigmaNeg);
				sigP[dir]=sigmaPlus;
				sigN[dir]=sigmaNeg;
			//	printf("dir=%d level=%d avLP=%f max=%f avLN=%f min=%f sigP=%f sigN=%f\n",dir,level,avLP[dir] ,maxL[dir], avLN[dir] ,minL[dir], sigP[dir], sigN[dir]);
			}
		
		}
		AvL=0.f;AvN=0.f;SL=0.f;SN=0.f;
		for (int dir=1; dir<4; dir++) {
			AvL +=avLP[dir];
			AvN +=avLN[dir];
			SL +=sigP[dir];
			SN +=sigN[dir];
		}
		AvL/=3;
		AvN/=3;
		SL/=3;
		SN/=3;
		mean[level]=AvL;
		meanN[level]=AvN;
		sigma[level]=SL;
		sigmaN[level]=SN;
		
		printf("Ind=%d Level=%d AvL=%f AvN=%f SL=%f SN=%f\n",ind, level,mean[level],meanN[level],sigma[level],sigmaN[level]);
		
	}	
	
	void ImProcFunctions::WaveletcontAllL(LabImage * labco, float ** varhue, float **varchrom, wavelet_decomposition &WaveletCoeffs_L, 
											struct cont_params cp){
		
		int maxlvl = WaveletCoeffs_L.maxlevel();
		int W_L = WaveletCoeffs_L.level_W(1);
		int H_L = WaveletCoeffs_L.level_H(1);
		float * WavCoeffs_L0 = WaveletCoeffs_L.coeff0;

		float maxh=2.5f;//amplification contrast above mean
		float maxl=2.5f; //reduction contrast under mean
		float contrast=cp.unif;
		float multL=(float)contrast*(maxl-1.f)/100.f + 1.f;	
		float multH=(float) contrast*(maxh-1.f)/100.f + 1.f;

		double avedbl=0.f; // use double precision for big summations
#ifdef _OPENMP
#pragma omp parallel for reduction(+:avedbl) num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
		for (int i=0; i<W_L*H_L; i++) {
			avedbl += WavCoeffs_L0[i];
		}

		float ave = avedbl / (double)(W_L*H_L);
		float av=ave/327.68f;
		float ah=(multH-1.f)/(av-100.f);//av ==> lumaref
		float bh=1.f-100.f*ah;
		float al=(multL-1.f)/av;
		float bl=1.f;
		float factorx=1.f;

		
#ifdef _OPENMP
#pragma omp parallel num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
{
	
#ifdef _OPENMP
#pragma omp for
#endif						
			for (int i=0; i<W_L*H_L; i++) {//contrast
				if(WavCoeffs_L0[i] < 32768.f) {
					float prov;
					if( WavCoeffs_L0[i]> ave) {
						float kh = ah*(WavCoeffs_L0[i]/327.68f)+bh;
						prov=WavCoeffs_L0[i];
						WavCoeffs_L0[i]=ave+kh*(WavCoeffs_L0[i]-ave);
					} else {
						float kl = al*(WavCoeffs_L0[i]/327.68f)+1.f;
						prov=WavCoeffs_L0[i];
						WavCoeffs_L0[i]=ave-kl*(ave-WavCoeffs_L0[i]);
					}
					float diflc=WavCoeffs_L0[i]-prov;
					diflc*=factorx; 
					WavCoeffs_L0[i] =  prov + diflc; 	
				}					
			}
		
#ifdef _OPENMP
#pragma omp for nowait
#endif				
			for (int i=0; i<W_L*H_L; i++) {
				int ii = i/W_L;
				int jj = i-ii*W_L;
				float LL = labco->L[ii*2][jj*2];						
				float LL100 = LL/327.68f;
				float tran = 5.f;//transition
				//shadow
				float alp=3.f;//increase contrast sahdow in lowlights  between 1 and ??
				if(cp.th > (100.f-tran))
					tran=100.f-cp.th;
				if(LL100 < cp.th){
					float aalp=(1.f-alp)/cp.th;//no changes for LL100 = cp.th
					float kk=aalp*LL100+alp;
					WavCoeffs_L0[i] *= (1.f+kk*cp.conres/200.f);
				}
				else if(LL100 < cp.th + tran) {
					float ath = -cp.conres/tran;
					float bth = cp.conres-ath*cp.th;
					WavCoeffs_L0[i] *= (1.f+(LL100*ath+bth)/200.f);
				}
				//highlight
				tran=5.f;
				if(cp.thH < (tran))
					tran = cp.thH;
				if(LL100 > cp.thH)
					WavCoeffs_L0[i] *= (1.f+cp.conresH/200.f);
				else if(LL100 > (cp.thH - tran)) {
					float athH = cp.conresH/tran;
					float bthH = cp.conresH-athH*cp.thH;
					WavCoeffs_L0[i] *= (1.f+(LL100*athH+bthH)/200.f);
				}
			}

#ifdef _OPENMP
#pragma omp for schedule(dynamic) collapse(2)
#endif		
		for (int dir=1; dir<4; dir++) {
			for (int lvl=0; lvl<maxlvl; lvl++) {

				int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
				int Hlvl_L = WaveletCoeffs_L.level_H(lvl);

				float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);

				ContAllL (labco,  varhue, varchrom, WavCoeffs_L, WavCoeffs_L0, lvl, dir, cp, Wlvl_L, Hlvl_L);		
			}
		}
}		
	}

	void ImProcFunctions::WaveletcontAllAB(LabImage * labco, float ** varhue, float **varchrom, wavelet_decomposition &WaveletCoeffs_ab,
											 struct cont_params cp, const bool useChannelA){
		
		int maxlvl = WaveletCoeffs_ab.maxlevel();
		int W_L = WaveletCoeffs_ab.level_W(1);
		int H_L = WaveletCoeffs_ab.level_H(1);
	
		float * WavCoeffs_ab0 = WaveletCoeffs_ab.coeff0;
		
		
		
		
#ifdef _OPENMP
#pragma omp parallel num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
{
#ifdef _OPENMP
#pragma omp for nowait
#endif				
			for (int i=0; i<W_L*H_L; i++) {
				const float skyprot = cp.sky;
				//chroma
				int ii = i/W_L;
				int jj = i-ii*W_L;
				float modhue = varhue[ii][jj];
				float scale=1.f;
				if(skyprot > 0.f){
					if((modhue < cp.t_ry && modhue > cp.t_ly)) {
						scale=(100.f-cp.sky)/100.1f;
					} else if((modhue >= cp.t_ry && modhue < cp.b_ry)) {
						scale=(100.f-cp.sky)/100.1f;
						float ar=(scale-1.f)/(cp.t_ry- cp.b_ry);
						float br=scale-cp.t_ry*ar;
						scale=ar*modhue+br;
					} else if((modhue > cp.b_ly && modhue < cp.t_ly)) {
						scale=(100.f-cp.sky)/100.1f;
						float al=(scale-1.f)/(-cp.b_ly + cp.t_ly);
						float bl=scale-cp.t_ly*al;
						scale=al*modhue+bl;
					}
				} else if(skyprot < 0.f){
					if((modhue > cp.t_ry || modhue < cp.t_ly)){
						scale=(100.f+cp.sky)/100.1f;
					}
				/*	else if((modhue >= cp.t_ry && modhue < cp.b_ry)) {
						scale=(100.f+cp.sky)/100.1f;
						float ar=(scale-1.f)/(cp.t_ry- cp.b_ry);
						float br=scale-cp.t_ry*ar;
						scale=ar*modhue+br;
					}
					else if((modhue > cp.b_ly && modhue < cp.t_ly)) {
						scale=(100.f+cp.sky)/100.1f;
						float al=(scale-1.f)/(-cp.b_ly + cp.t_ly);
						float bl=scale-cp.t_ly*al;
						scale=al*modhue+bl;
					}
				*/		
				}
				WavCoeffs_ab0[i]*=(1.f+cp.chrores*(scale)/100.f);
			}

#ifdef _OPENMP
#pragma omp for schedule(dynamic) collapse(2)
#endif		
		for (int dir=1; dir<4; dir++) {
			for (int lvl=0; lvl<maxlvl; lvl++) {

				int Wlvl_ab = WaveletCoeffs_ab.level_W(lvl);
				int Hlvl_ab = WaveletCoeffs_ab.level_H(lvl);


				int skip_ab = WaveletCoeffs_ab.level_stride(lvl);
				//printf("lev=%d skipL=%d  skipab=%d\n",lvl, skip_L,skip_ab);
				float ** WavCoeffs_ab = WaveletCoeffs_ab.level_coeffs(lvl);
				ContAllAB (labco,  varhue, varchrom, WavCoeffs_ab, WavCoeffs_ab0, lvl, dir, cp, Wlvl_ab, Hlvl_ab, useChannelA);		
			}
		}
}		
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	void ImProcFunctions::ContAllL (LabImage * labco, float ** varhue, float **varchrom, float ** WavCoeffs_L, float * WavCoeffs_L0, int level, int dir, struct cont_params cp,
									int W_L, int H_L)
	{
		const float eps = 0.01f;
		
		//float thr= params->wavelet.thres;
		float cpMul = cp.mul[level];
		if(cpMul != 0.f) { // cpMul == 0.f means all will be multiplied by 1.f, so we can skip this
		
			const float skinprot = params->wavelet.skinprotect;
			const float skinprotneg = -skinprot;
			const float factorHard = (1.f - skinprotneg/100.f);
		
			//to adjust increase contrast with local contrast
	
			//for each pixel
		//	float k[8]={0.85f, 0.7f, 0.55f, 0.4f, 0.3f, 0.25f, 0.2f, 0.1f};//values to tested with several images
			float k[9]={0.95f, 0.85f, 0.7f, 0.6f, 0.45f, 0.3f, 0.2f, 0.15f, 0.1f};//values to tested with several images
			//float meath[8]={700.f, 1400.f, 1900.f, 2200.f, 2800.f, 3500.f, 4500.f, 6000.f};//values to tested with several images
			float meath[9]={1000.f, 1500.f, 2000.f, 2500.f, 3000.f, 3500.f, 4000.f, 4500.f, 6000.f};//values to tested with several images
			float ampli[9]={1.2f, 1.4f, 1.7f, 2.2f, 2.5f, 3.f, 3.5f, 4.f, 4.5f};
			float mea[9];
			float tr=cp.th;//suppress 2 slider 
			tr=90.f;
		
			for(int j=0;j<9;j++) mea[j]=meath[j]*(1.f+(ampli[j]-1.f)*(tr/100.f));
			//
			//float uni=(float) cp.unif;
			float uni = 95.f;
			float bbet=1.f;
			float abet[9];
			for(int h=0;h<9;h++) abet[h]=((k[h]-1.f)/100.f)*uni+bbet;
			float beta;
		
			bool skinControl = (skinprot != 0.f);
			bool useChromAndHue = (skinprot != 0.f || cp.HSmet);
			float modchro, kLlev;
	
			for (int i=0; i<W_L*H_L; i++) {
				kLlev=1.f;

				if(cpMul<0.f) {
					beta=1.f;// disabled for negatives values "less contrast"
				} else {
					float WavCL = fabsf(WavCoeffs_L[dir][i]);
					if(WavCL < 30.f) beta=0.6f;//preserve very low contrast (sky...)
					else if(WavCL < 100.f) beta=0.8f;
					else if(WavCL < mea[0]) beta=1.f;//no changes
					else if(WavCL < mea[1]) beta=abet[0];//linear regression
					else if(WavCL < mea[2]) beta=abet[1];
					else if(WavCL < mea[3]) beta=abet[2];
					else if(WavCL < mea[4]) beta=abet[3];
					else if(WavCL < mea[5]) beta=abet[4];
					else if(WavCL < mea[6]) beta=abet[5];
					else if(WavCL < mea[7]) beta=abet[6];
					else if(WavCL < mea[8]) beta=abet[7];
					// next condition is automatically true, so skip the if
					else beta=abet[8];
				}
				float scale = 1.f;
				float LL100;
				if(useChromAndHue) {
					int ii=i/W_L;
					int jj=i-ii*W_L;
					float LL = labco->L[ii*2][jj*2];
					LL100=LL/327.68f;
					float modhue = varhue[ii][jj];
					modchro = varchrom[ii*2][jj*2];
					// hue chroma skin with initial lab datas
					scale=1.f;
					if(skinprot > 0.f){
						Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprot, scale, true, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0);	//0 for skin and extand
					} else if(skinprot < 0.f){
						Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprotneg, scale, false, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0);	
						if (scale == 1.f)
							scale=factorHard;
						else
							scale=1.f;
					}
				}
				//linear transition HL
				float alpha = (1024.f + 15.f *(float) cpMul*scale*beta)/1024.f ;	
				if(cp.HSmet){					
					float aaal=(1.f-alpha)/(cp.b_lhl-cp.t_lhl);
					float bbal=1.f-aaal*cp.b_lhl;
					float aaar=(alpha-1.f)/(cp.t_rhl-cp.b_rhl);
					float bbbr=1.f-cp.b_rhl*aaar;
					//linear transition Shadows
					float aaalS=(1.f-alpha)/(cp.b_lsl-cp.t_lsl);
					float bbalS=1.f-aaalS*cp.b_lsl;
					float aaarS=(alpha-1.f)/(cp.t_rsl-cp.b_rsl);
					float bbbrS=1.f-cp.b_rsl*aaarS;
					if(level <=cp.numlevH) {//in function of levels
						if((LL100 > cp.t_lhl && LL100 < cp.t_rhl)) kLlev=alpha; 
						else if((LL100 > cp.b_lhl && LL100 <= cp.t_lhl)) kLlev=aaal*LL100+bbal; 	
						else if((LL100 > cp.t_rhl && LL100 <= cp.b_rhl)) kLlev=aaar*LL100+bbbr; 	
						else 	kLlev=1.f;
					}
					if(level >=(9-cp.numlevS)) {
						if((LL100 > cp.t_lsl && LL100 < cp.t_rsl)) kLlev=alpha; 
						else if((LL100 > cp.b_lsl && LL100 <= cp.t_lsl)) kLlev=aaalS*LL100+bbalS; 	
						else if((LL100 > cp.t_rsl && LL100 <= cp.b_rsl)) kLlev=aaarS*LL100+bbbrS; 	
						else 	kLlev=1.f;
					} else
						kLlev=alpha;
				}
				else kLlev=alpha;	
				WavCoeffs_L[dir][i]*=(kLlev);
			}
		}

		
		
		// to see each level of wavelet ...level from 0 to 7
		int choicelevel=0;
		if(params->wavelet.Lmethod=="0_") choicelevel=0;	
		else if(params->wavelet.Lmethod=="1_") choicelevel=1;	
		else if(params->wavelet.Lmethod=="2_") choicelevel=2;	
		else if(params->wavelet.Lmethod=="3_") choicelevel=3;	
		else if(params->wavelet.Lmethod=="4_") choicelevel=4;	
		else if(params->wavelet.Lmethod=="5_") choicelevel=5;	
		else if(params->wavelet.Lmethod=="6_") choicelevel=6;	
		else if(params->wavelet.Lmethod=="7_") choicelevel=7;	
		else if(params->wavelet.Lmethod=="8_") choicelevel=8;	
		int choiceClevel=0;
		if(params->wavelet.CLmethod=="one") choiceClevel=0;	
		else if(params->wavelet.CLmethod=="inf") choiceClevel=1;	
		else if(params->wavelet.CLmethod=="sup") choiceClevel=2;	
		else if(params->wavelet.CLmethod=="all") choiceClevel=3;	
		int choiceDir=0;
		if(params->wavelet.Dirmethod=="one") choiceDir=1;	
		else if(params->wavelet.Dirmethod=="two") choiceDir=2;	
		else if(params->wavelet.Dirmethod=="thr") choiceDir=3;	
		else if(params->wavelet.Dirmethod=="all") choiceDir=0;	

	//	printf("LUm lev=%d clev=%d dir=%d\n",choicelevel,choiceClevel,choiceDir);
	if(choiceClevel==0){
		if(choiceDir==0){
			if(level != choicelevel){
				for (int dir=1; dir<4; dir++) {	
					for (int i=0; i<W_L*H_L; i++) {
						WavCoeffs_L[dir][i] =0.f;
						WavCoeffs_L0[i] =0.f;
					}
				}
			}
		}
		else {
		if(level != choicelevel){
		
				for (int i=0; i<W_L*H_L; i++) {
					WavCoeffs_L0[i] =0.f;	
					
					if(choiceDir==1)	{WavCoeffs_L[2][i] =0.f;WavCoeffs_L[3][i] =0.f;}
					else if(choiceDir==2)	{WavCoeffs_L[1][i] =0.f;WavCoeffs_L[3][i] =0.f;}
					else if(choiceDir==3)	{WavCoeffs_L[2][i] =0.f;WavCoeffs_L[1][i] =0.f;}
				}
			}
		}
	}
	else if	(choiceClevel==1){
		if(choiceDir==0){
			if(level >= choicelevel){
				for (int dir=1; dir<4; dir++) {	
					for (int i=0; i<W_L*H_L; i++) {
						WavCoeffs_L[dir][i] =0.f;
						WavCoeffs_L0[i] =0.f;
					}
				}
			}
		}
		else {
		if(level >= choicelevel){
		
				for (int i=0; i<W_L*H_L; i++) {
					WavCoeffs_L0[i] =0.f;
					
					if(choiceDir==1)	{WavCoeffs_L[2][i] =0.f;WavCoeffs_L[3][i] =0.f;}
					else if(choiceDir==2)	{WavCoeffs_L[1][i] =0.f;WavCoeffs_L[3][i] =0.f;}
					else if(choiceDir==3)	{WavCoeffs_L[2][i] =0.f;WavCoeffs_L[1][i] =0.f;}
				}
			}
		}
	}
	else if	(choiceClevel==2){
		if(choiceDir==0){
			if(level <= choicelevel){
				for (int dir=1; dir<4; dir++) {	
					for (int i=0; i<W_L*H_L; i++) {
					WavCoeffs_L[dir][i] =0.f;
					}
				}
			}
		}
		else {
			if(level <= choicelevel){
				for (int i=0; i<W_L*H_L; i++) {
					
					if(choiceDir==1)	{WavCoeffs_L[2][i] =0.f;WavCoeffs_L[3][i] =0.f;}
					else if(choiceDir==2)	{WavCoeffs_L[1][i] =0.f;WavCoeffs_L[3][i] =0.f;}
					else if(choiceDir==3)	{WavCoeffs_L[2][i] =0.f;WavCoeffs_L[1][i] =0.f;}
				}
			}
		}
	}
		
	}

	void ImProcFunctions::ContAllAB (LabImage * labco, float ** varhue, float **varchrom, float ** WavCoeffs_ab, float * WavCoeffs_ab0, int level, int dir, struct cont_params cp,
									int W_ab, int H_ab, const bool useChannelA)
	{
		const float eps = 0.01f;
		
		//float thr= params->wavelet.thres;
		float cpMul = cp.mul[level];
		if(cpMul != 0.f && cp.CHmet==2 && cp.chro != 0.f) { // cpMul == 0.f means all will be multiplied by 1.f, so we can skip this
			const float skinprot = params->wavelet.skinprotect;
			const float skinprotneg = -skinprot;
			const float factorHard = (1.f - skinprotneg/100.f);
			const float cpChrom = cp.chro;

			//to adjust increase contrast with local contrast
			bool useChromAndHue = (skinprot != 0.f);
			float modchro;
	
			for (int i=0; i<W_ab*H_ab; i++) {
				float scale = 1.f;
				float scaleSK = 1.f ;
				if(useChromAndHue) {
					int ii=i/W_ab;
					int jj=i-ii*W_ab;
					float LL = labco->L[ii*2][jj*2];
					float LL100=LL/327.68f;
					float modhue = varhue[ii][jj];
					modchro = varchrom[ii*2][jj*2];
					// hue chroma skin with initial lab datas
					scale=1.f;
					if(skinprot > 0.f){
						Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprot, scale, true, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0);	//0 for skin and extand
					} else if(skinprot < 0.f){
						Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprotneg, scale, false, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0);	
						if (scale == 1.f)
							scale=factorHard;
						else
							scale=1.f;
					}
				}
				float alphaC =(1024.f + 15.f *cpMul*cpChrom*scale*scaleSK/50.f)/1024.f ;
				WavCoeffs_ab[dir][i]*=alphaC;
			}
		}
		//Curve chro
		bool useOpacity;
		float mulOpacity;
		if(useChannelA) {
			useOpacity = cp.opaRG;
			mulOpacity = cp.mulopaRG[level];
		}
		else {
			useOpacity = cp.opaBY;
			mulOpacity = cp.mulopaBY[level];
		}
		if(cp.curv && cp.CHmet!=2 && level < 9) {
			float modchro, modhue, kLlev, kClev;
			float cpMulC = cp.mulC[level];
			const float skinprot = params->wavelet.skinprotect;
			const float skinprotneg = -skinprot;
			const float factorHard = (1.f - skinprotneg/100.f);
			bool skinControl = (skinprot != 0.f);
		
			for (int i=0; i<W_ab*H_ab; i++) {
				int ii=i/W_ab;
				int jj=i-ii*W_ab;
				//WL and W_ab are identical
				float LL = labco->L[ii*2][jj*2];
				float LL100=LL/327.68f;
				float scale=1.f;
				modchro = varchrom[ii*2][jj*2];
				float modhue = varhue[ii][jj];
				
				if(skinControl) {
					// hue chroma skin with initial lab datas
					modhue = varhue[ii][jj];
				
					scale=1.f;
					if(skinprot > 0.f){
						Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprot, scale, true, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 1);	//1 for curve
					}
					else if(skinprot < 0.f){
						Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprotneg, scale, false, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 1);	
						if (scale == 1.f) { scale=factorHard;}
						else scale=1.f;
					}
				}
				float scaleSK =1.f;
				float beta = (1024.f + 20.f *(float) cpMulC*scale*scaleSK)/1024.f ;
				if(beta < 0.02f) beta=0.02f;	  
				//linear for saturated
				float aaal=(1.f-beta)/(cp.b_lsat-cp.t_lsat);
				float bbal=1.f-aaal*cp.b_lsat;
				float aaar=(beta-1.f)/(cp.t_rsat-cp.b_rsat);
				float bbbr=1.f-cp.b_rsat*aaar;
				//linear for pastel
				float aaalS=(1.f-beta)/(cp.b_lpast-cp.t_lpast);
				float bbalS=1.f-aaalS*cp.b_lpast;
				float aaarS=(beta-1.f)/(cp.t_rpast-cp.b_rpast);
				float bbbrS=1.f-cp.b_rpast*aaarS;
				kClev=beta;
				if(cp.CHmet==1){					
					if(level < cp.chrom) {
							if((modchro > cp.t_lsat && modchro < cp.t_rsat)) kClev=beta; 
							else if((modchro > cp.b_lsat && modchro <= cp.t_lsat)) kClev=aaal*modchro+bbal; 	
							else if((modchro > cp.t_rsat &&  modchro <= cp.b_rsat)) kClev=aaar*modchro+bbbr; 	
							else 	kClev=1.f;
					}
					if(level >= cp.chrom) {
							if((modchro > cp.t_lpast && modchro < cp.t_rpast)) kClev=beta; 
							else if((modchro > cp.b_lpast && modchro <= cp.t_lpast)) kClev=aaalS*modchro+bbalS; 	
							else if((modchro > cp.t_rpast &&  modchro <= cp.b_rpast)) kClev=aaarS*modchro+bbbrS; 	
							else 	kClev=1.f;
					}
				}
				else if(cp.CHmet==0)kClev=beta;	
				WavCoeffs_ab[dir][i]*=kClev;
			}
		}	

		if(useOpacity && level < 9) { //toning
			float betaRG = (1024.f + 20.f *(float) mulOpacity)/1024.f ;
			for (int i=0; i<W_ab*H_ab; i++)
				WavCoeffs_ab[dir][i]*=betaRG;			
		}
		
		
		// to see each level of wavelet ...level from 0 to 7
		int choicelevel=0;
		if(params->wavelet.Lmethod=="0_") choicelevel=0;	
		else if(params->wavelet.Lmethod=="1_") choicelevel=1;	
		else if(params->wavelet.Lmethod=="2_") choicelevel=2;	
		else if(params->wavelet.Lmethod=="3_") choicelevel=3;	
		else if(params->wavelet.Lmethod=="4_") choicelevel=4;	
		else if(params->wavelet.Lmethod=="5_") choicelevel=5;	
		else if(params->wavelet.Lmethod=="6_") choicelevel=6;	
		else if(params->wavelet.Lmethod=="7_") choicelevel=7;	
		else if(params->wavelet.Lmethod=="8_") choicelevel=8;	
		int choiceClevel=0;
		if(params->wavelet.CLmethod=="one") choiceClevel=0;	
		else if(params->wavelet.CLmethod=="inf") choiceClevel=1;	
		else if(params->wavelet.CLmethod=="sup") choiceClevel=2;	
		else if(params->wavelet.CLmethod=="all") choiceClevel=3;	
		int choiceDir=0;
		if(params->wavelet.Dirmethod=="one") choiceDir=1;	
		else if(params->wavelet.Dirmethod=="two") choiceDir=2;	
		else if(params->wavelet.Dirmethod=="thr") choiceDir=3;	
		else if(params->wavelet.Dirmethod=="all") choiceDir=0;	
	//	printf("CHRO lev=%d clev=%d dir=%d\n",choicelevel,choiceClevel,choiceDir);

		
	if(choiceClevel==0){
		if(choiceDir==0){
			if(level != choicelevel){
				for (int dir=1; dir<4; dir++) {	
					for (int i=0; i<W_ab*H_ab; i++) {
					WavCoeffs_ab0[i] =0.f;
					
					}
				}
			}
		}
		else {
		if(level != choicelevel){
		
				for (int i=0; i<W_ab*H_ab; i++) {
					WavCoeffs_ab0[i] =0.f;
					
				}
			}
		}
	}
	else if	(choiceClevel==1){
		if(choiceDir==0){
			if(level >= choicelevel){
				for (int dir=1; dir<4; dir++) {	
					for (int i=0; i<W_ab*H_ab; i++) {
					WavCoeffs_ab0[i] =0.f;
					
					}
				}
			}
		}
		else {
		if(level >= choicelevel){
		
				for (int i=0; i<W_ab*H_ab; i++) {
					WavCoeffs_ab0[i] =0.f;
				}
			}
		}
	}
	else if	(choiceClevel==2){
		if(choiceDir==0){
			if(level <= choicelevel){
				for (int dir=1; dir<4; dir++) {	
					for (int i=0; i<W_ab*H_ab; i++) {
					WavCoeffs_ab[dir][i] =0.f;
					}
				}
			}
		}
		else {
			if(level <= choicelevel){
				for (int i=0; i<W_ab*H_ab; i++) {
					
					if(choiceDir==1)	{WavCoeffs_ab[2][i] =0.f;WavCoeffs_ab[3][i] =0.f;}
					else if(choiceDir==2)	{WavCoeffs_ab[1][i] =0.f;WavCoeffs_ab[3][i] =0.f;}
					else if(choiceDir==3)	{WavCoeffs_ab[2][i] =0.f;WavCoeffs_ab[1][i] =0.f;}
				}
			}
		}
	}
		
	}

}
