/*
 *  This file is part of RawTherapee.
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
 *
 *  � 2010 Emil Martinec <ejmartin@uchicago.edu>
 *    
 */

//#include "rtengine.h"
#include <cstddef>
#include <cmath>
#include "curves.h"
#include "labimage.h"
#include "improcfun.h"
#include "array2D.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define SQR(x) ((x)*(x))
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define CLIP(a) (CLIPTO(a,0,65535))


#define DIRWT_L(i1,j1,i,j) (  rangefn_L[abs(data_fine->L[i1][j1]-data_fine->L[i][j])] )

#define DIRWT_AB(i1,j1,i,j) (  rangefn_ab[/*abs(data_fine->L[i1][j1]-data_fine->L[i][j])*/0 + \
abs(data_fine->a[i1][j1]-data_fine->a[i][j]) + \
abs(data_fine->b[i1][j1]-data_fine->b[i][j])] )


#define NRWT_AB (nrwt_ab[abs(hipass[1])] * nrwt_ab[abs(hipass[2])])


#define med3(a,b,c) (a<b ? (b<c ? b : (a<c ? c : a)) : (a<c ? a : (b<c ? c : b)))

#define hmf(a11,a12,a13,a21,a22,a23,a31,a32,a33) \
(med3(a22,med3(a22,med3(a12,a22,a32),med3(a21,a22,a23)), \
		med3(a22,med3(a11,a22,a33),med3(a13,a22,a31))) )

#define PIX_SORT(a,b) { if ((a)>(b)) {temp=(a);(a)=(b);(b)=temp;} }

#define med3x3(a0,a1,a2,a3,a4,a5,a6,a7,a8,median) { \
p[0]=a0; p[1]=a1; p[2]=a2; p[3]=a3; p[4]=a4; p[5]=a5; p[6]=a6; p[7]=a7; p[8]=a8; \
PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]); \
PIX_SORT(p[0],p[1]); PIX_SORT(p[3],p[4]); PIX_SORT(p[6],p[7]); \
PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]); \
PIX_SORT(p[0],p[3]); PIX_SORT(p[5],p[8]); PIX_SORT(p[4],p[7]); \
PIX_SORT(p[3],p[6]); PIX_SORT(p[1],p[4]); PIX_SORT(p[2],p[5]); \
PIX_SORT(p[4],p[7]); PIX_SORT(p[4],p[2]); PIX_SORT(p[6],p[4]); \
PIX_SORT(p[4],p[2]); median=p[4];} //a4 is the median


namespace rtengine {
	
	static const int maxlevel = 4;
	
	//sequence of scales
	//static const int scales[8] = {1,2,4,8,16,32,64,128};
	//sequence of pitches
	//static const int pitches[8] = {1,1,1,1,1,1,1,1};
	
	//sequence of scales
	//static const int scales[8] = {1,1,1,1,1,1,1,1};
	//sequence of pitches
	//static const int pitches[8] = {2,2,2,2,2,2,2,2};
	
	//sequence of scales
	//static const int scales[8] = {1,1,2,2,4,4,8,8};
	//sequence of pitches
	//static const int pitches[8] = {2,1,2,1,2,1,2,1};
	
	//sequence of scales
	static const int scales[8] = {1,1,2,4,8,16,32,64};
	//sequence of pitches
	static const int pitches[8] = {2,1,1,1,1,1,1,1};
	
	//pitch is spacing of subsampling
	//scale is spacing of directional averaging weights
	//example 1: no subsampling at any level -- pitch=1, scale=2^n
	//example 2: subsampling by 2 every level -- pitch=2, scale=1 at each level
	//example 3: no subsampling at first level, subsampling by 2 thereafter -- 
	//	pitch =1, scale=1 at first level; pitch=2, scale=2 thereafter
	
	
	
	
	void ImProcFunctions :: dirpyrLab_denoise(LabImage * src, LabImage * dst, const procparams::DirPyrDenoiseParams & dnparams )
	{

		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		LUTf rangefn_L(65536);
		LUTf nrwt_l(1);
		
		LUTf rangefn_ab(65536);
		LUTf nrwt_ab(65536);
		
		//set up NR weight functions
		
		float noise_L = 10.0*dnparams.luma;
		float noisevar_L = SQR(noise_L);
		
		float noise_ab = 10.0*dnparams.chroma;
		float noisevar_ab = SQR(noise_ab);
		
		
		//set up range functions
		for (int i=0; i<65536; i++) 
			rangefn_L[i] = exp(-(double)i / (1.0f+noise_L));// * (1.0+noisevar_L)/((double)(i*i) + noisevar_L+1.0); 
		for (int i=0; i<65536; i++) 
			rangefn_ab[i]= exp(-SQR((double)i) / (1.0f+3*noisevar_ab));// * (1.0+noisevar_ab)/((double)(i*i) + noisevar_ab+1.0); 
		
		
		for (int i=0; i<65536; i++) 
			nrwt_ab[i] = ((1.0+abs(i-32768)/(1.0+8*noise_ab)) * exp(-(double)fabs(i-32768)/ (1.0+8*noise_ab) ) );
		
		
		//for (int i=0; i<65536; i+=100)  printf("%d %d \n",i,gamcurve[i]);
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		int level;
		
		LabImage * dirpyrLablo[maxlevel]; 
		int w = (int)((src->W-1)/pitches[0])+1;
		int h = (int)((src->H-1)/pitches[0])+1;
		dirpyrLablo[0] = new LabImage(w, h);
		for (level=1; level<maxlevel; level++) {
			w = (int)((w-1)/pitches[level])+1;
			h = (int)((h-1)/pitches[level])+1;
			dirpyrLablo[level] = new LabImage(w, h);
		};
		
		
		//////////////////////////////////////////////////////////////////////////////
		
		level = 0;
		
		int scale = scales[level];
		int pitch = pitches[level];
		//int thresh = 10 * c[8];
		//impulse_nr (src, src, m_w1, m_h1, thresh, noisevar);
		
		dirpyr(src, dirpyrLablo[0], 0, rangefn_L, rangefn_ab, pitch, scale, dnparams.luma, dnparams.chroma );
		
		level = 1;
		
		while(level < maxlevel)
		{
			scale = scales[level];
			pitch = pitches[level];
			
			dirpyr(dirpyrLablo[level-1], dirpyrLablo[level], level, rangefn_L, rangefn_ab, pitch, scale, dnparams.luma, dnparams.chroma );
			
			level ++;
		}
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		for(int level = maxlevel - 1; level > 0; level--)
		{
			
			int scale = scales[level];
			int pitch = pitches[level];
			idirpyr(dirpyrLablo[level], dirpyrLablo[level-1], level, rangefn_L, nrwt_l, nrwt_ab, pitch, scale, dnparams.luma, dnparams.chroma/*, Lcurve, abcurve*/ );
		}
		
		
		scale = scales[0];
		pitch = pitches[0];
		// freeing as much memory as possible since the next call to idirpyr will need lots
		for(int i = 1; i < maxlevel; i++) {
			delete dirpyrLablo[i];
		}

		idirpyr(dirpyrLablo[0], dst, 0, rangefn_L, nrwt_l, nrwt_ab, pitch, scale, dnparams.luma, dnparams.chroma/*, Lcurve, abcurve*/ );

		// freeing the last bunch of memory
		delete dirpyrLablo[0];

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
	};
	
	void ImProcFunctions::dirpyr(LabImage* data_fine, LabImage* data_coarse, int level, \
								 LUTf & rangefn_L, LUTf & rangefn_ab, int pitch, int scale, \
								 const int luma, const int chroma )
	{
		
		//pitch is spacing of subsampling
		//scale is spacing of directional averaging weights
		//example 1: no subsampling at any level -- pitch=1, scale=2^n
		//example 2: subsampling by 2 every level -- pitch=2, scale=1 at each level
		//example 3: no subsampling at first level, subsampling by 2 thereafter -- 
		//	pitch =1, scale=1 at first level; pitch=2, scale=2 thereafter
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// calculate weights, compute directionally weighted average
		
		int width = data_fine->W;
		int height = data_fine->H;
		
		//generate domain kernel 
		int halfwin = 3;//MIN(ceil(2*sig),3);
		int scalewin = halfwin*scale;
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
 
		for(int i = 0; i < height; i+=pitch ) { int i1=i/pitch;
			for(int j = 0, j1=0; j < width; j+=pitch, j1++)
			{	
				float dirwt_l, dirwt_ab, norm_l, norm_ab;
				//float Lmed,Lhmf;
				//float lops,aops,bops;
				float Lout, aout, bout;
				norm_l = norm_ab = 0;//if we do want to include the input pixel in the sum
				Lout = 0;
				aout = 0;
				bout = 0;
				
				for(int inbr=MAX(0,(i-scalewin)); inbr<=MIN(height-1,(i+scalewin)); inbr+=scale) {
					for (int jnbr=MAX(0,(j-scalewin)); jnbr<=MIN(width-1,(j+scalewin)); jnbr+=scale) {
				/*for(int inbr=(i-scalewin); inbr<=(i+scalewin); inbr+=scale) {
					if (inbr<0 || inbr>height-1) continue;
					for (int jnbr=(j-scalewin); jnbr<=(j+scalewin); jnbr+=scale) {
						if (jnbr<0 || jnbr>width-1) continue;*/
						dirwt_l = DIRWT_L(inbr, jnbr, i, j);
						dirwt_ab = DIRWT_AB(inbr, jnbr, i, j);
						Lout += dirwt_l*data_fine->L[inbr][jnbr];
						aout += dirwt_l*dirwt_ab*data_fine->a[inbr][jnbr];
						bout += dirwt_l*dirwt_ab*data_fine->b[inbr][jnbr];
						norm_l += dirwt_l;
						norm_ab += dirwt_l*dirwt_ab;
					}
				}
				//lops = Lout/norm;//diagnostic
				//aops = aout/normab;//diagnostic
				//bops = bout/normab;//diagnostic
				
				data_coarse->L[i1][j1]=Lout/norm_l;//low pass filter
				data_coarse->a[i1][j1]=aout/norm_ab;
				data_coarse->b[i1][j1]=bout/norm_ab;
				
			}
		}
		
		
		
		
	};
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::idirpyr(LabImage* data_coarse, LabImage* data_fine, int level, LUTf &rangefn_L, LUTf & nrwt_l, LUTf & nrwt_ab, \
								  int pitch, int scale, const int luma, const int chroma/*, LUTf & Lcurve, LUTf & abcurve*/ )
	{
		
		int width = data_fine->W;
		int height = data_fine->H;
		
		//array2D<float> nrfactorL (width,height);
		
		float noisevar_L = 4*SQR(25.0 * luma);
		float noisevar_ab = 2*SQR(100.0 * chroma);
		float scalefactor = 1.0/pow(2.0,(level+1)*2);//change the last 2 to 1 for longer tail of higher scale NR

		noisevar_L *= scalefactor;
				
		// for coarsest level, take non-subsampled lopass image and subtract from lopass_fine to generate hipass image
		
		// denoise hipass image, add back into lopass_fine to generate denoised image at fine scale
		
		// now iterate:
		// (1) take denoised image at level n, expand and smooth using gradient weights from lopass image at level n-1
		//     the result is the smoothed image at level n-1
		// (2) subtract smoothed image at level n-1 from lopass image at level n-1 to make hipass image at level n-1
		// (3) denoise the hipass image at level n-1
		// (4) add the denoised image at level n-1 to the smoothed image at level n-1 to make the denoised image at level n-1
		
		// note that the coarsest level amounts to skipping step (1) and doing (2,3,4).
		// in other words, skip step one if pitch=1
		
		// step (1)
		
		if (pitch==1) { 
			
			// step (1-2-3-4) 
			
#ifdef _OPENMP
#pragma omp parallel
#endif
{
			
#ifdef _OPENMP
#pragma omp for
#endif
			for(int  i = 0; i < height; i++)
				for(int  j = 0; j < width; j++) {
					double wtdsum[3], norm;
					float hipass[3], hpffluct[3], tonefactor, nrfactora, nrfactorb;
										
					hipass[1] = data_fine->a[i][j]-data_coarse->a[i][j];
					hipass[2] = data_fine->b[i][j]-data_coarse->b[i][j];
					
					//Wiener filter
					//chroma
					hpffluct[1]=SQR(hipass[1])+0.001;
					hpffluct[2]=SQR(hipass[2])+0.001;
					nrfactora = (hpffluct[1]) /((hpffluct[1]) + noisevar_ab * NRWT_AB);
					nrfactorb = (hpffluct[2]) /((hpffluct[2]) + noisevar_ab * NRWT_AB);

					hipass[1] *= nrfactora;
					hipass[2] *= nrfactorb;
					
					data_fine->a[i][j] = hipass[1]+data_coarse->a[i][j];
					data_fine->b[i][j] = hipass[2]+data_coarse->b[i][j];
				}

	
}//end of pitch=1
			
		} else {//pitch>1
			
			LabImage* smooth;
			
			smooth = new LabImage(width, height);
#ifdef _OPENMP
#pragma omp parallel
#endif

{
		
#ifdef _OPENMP
#pragma omp for
#endif
			for(int  i = 0; i < height; i+=pitch)
			{
				int ix=i/pitch;
				for(int  j = 0, jx=0; j < width; j+=pitch, jx++) {
					
					//copy common pixels
					smooth->L[i][j] = data_coarse->L[ix][jx];
					smooth->a[i][j] = data_coarse->a[ix][jx];
					smooth->b[i][j] = data_coarse->b[ix][jx];
				}
			}
			//if (pitch>1) {//pitch=2; step (1) expand coarse image, fill in missing data
#ifdef _OPENMP
#pragma omp for	
#endif
			for(int  i = 0; i < height-1; i+=2)
				for(int j = 0; j < width-1; j+=2) {
					//do midpoint first
					double norm=0.0,wtdsum[3]={0.0,0.0,0.0};
					//wtdsum[0]=wtdsum[1]=wtdsum[2]=0.0;
					for(int ix=i; ix<MIN(height,i+3); ix+=2)
						for (int jx=j; jx<MIN(width,j+3); jx+=2) {
							wtdsum[0] += smooth->L[ix][jx];
							wtdsum[1] += smooth->a[ix][jx];
							wtdsum[2] += smooth->b[ix][jx];
							norm++;
						}
					norm = 1/norm;
					smooth->L[i+1][j+1]=wtdsum[0]*norm;
					smooth->a[i+1][j+1]=wtdsum[1]*norm;
					smooth->b[i+1][j+1]=wtdsum[2]*norm;
				}
#ifdef _OPENMP
#pragma omp for
#endif
			
			for(int i = 0; i < height-1; i+=2)
				for(int j = 0; j < width-1; j+=2) {
					//now right neighbor
					if (j+1==width) continue;
					double norm=0.0,wtdsum[3]={0.0,0.0,0.0};

					for (int jx=j; jx<MIN(width,j+3); jx+=2) {
						wtdsum[0] += smooth->L[i][jx];
						wtdsum[1] += smooth->a[i][jx];
						wtdsum[2] += smooth->b[i][jx];
						norm++;
					}
					for (int ix=MAX(0,i-1); ix<MIN(height,i+2); ix+=2) {
						wtdsum[0] += smooth->L[ix][j+1];
						wtdsum[1] += smooth->a[ix][j+1];
						wtdsum[2] += smooth->b[ix][j+1];
						norm++;
					}
					norm = 1/norm;
					smooth->L[i][j+1]=wtdsum[0]*norm;
					smooth->a[i][j+1]=wtdsum[1]*norm;
					smooth->b[i][j+1]=wtdsum[2]*norm;
					
					//now down neighbor
					if (i+1==height) continue;
					norm=0.0;wtdsum[0]=wtdsum[1]=wtdsum[2]=0.0;
					for (int ix=i; ix<MIN(height,i+3); ix+=2) {
						wtdsum[0] += smooth->L[ix][j];
						wtdsum[1] += smooth->a[ix][j];
						wtdsum[2] += smooth->b[ix][j];
						norm++;
					}
					for (int jx=MAX(0,j-1); jx<MIN(width,j+2); jx+=2) {
						wtdsum[0] += smooth->L[i+1][jx];
						wtdsum[1] += smooth->a[i+1][jx];
						wtdsum[2] += smooth->b[i+1][jx];
						norm++;
					}
					norm=1/norm;
					smooth->L[i+1][j]=wtdsum[0]*norm;
					smooth->a[i+1][j]=wtdsum[1]*norm;
					smooth->b[i+1][j]=wtdsum[2]*norm;
					
				}
			
#ifdef _OPENMP
#pragma omp for	
#endif
		
			// step (2-3-4) 
			for( int i = 0; i < height; i++)
				for(int j = 0; j < width; j++) {
					
					float hipass[3], hpffluct[3], nrfactora, nrfactorb;
					
					hipass[1] = data_fine->a[i][j]-smooth->a[i][j];
					hipass[2] = data_fine->b[i][j]-smooth->b[i][j];
					
					//Wiener filter
					//chroma
					hpffluct[1]=SQR(hipass[1])+0.001;
					hpffluct[2]=SQR(hipass[2])+0.001;
					nrfactora = (hpffluct[1]) /((hpffluct[1]) + noisevar_ab * NRWT_AB /* * abcurve[smooth->L[i][j]]*/);
					nrfactorb = (hpffluct[2]) /((hpffluct[2]) + noisevar_ab * NRWT_AB /* * abcurve[smooth->L[i][j]]*/);

					hipass[1] *= nrfactora;
					hipass[2] *= nrfactorb;
					
					data_fine->a[i][j] = hipass[1]+smooth->a[i][j];
					data_fine->b[i][j] = hipass[2]+smooth->b[i][j];
				}
	
}	// end parallel		
			delete smooth;
		}//end of pitch>1
		
	};
	
	
#undef DIRWT_L
#undef DIRWT_AB
	
//#undef NRWT_L	
#undef NRWT_AB	
	
}

