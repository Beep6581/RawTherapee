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
 *  © 2010 Emil Martinec <ejmartin@uchicago.edu>
 *    
 */

//#include <rtengine.h>
#include <cstddef>
#include <math.h>
#include <curves.h>
#include <labimage.h>
#include <improcfun.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SQR(x) ((x)*(x))
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define CLIP(a) (CLIPTO(a,0,65535))



#define DIRWT_L(i1,j1,i,j) (/*domker[(i1-i)/scale+halfwin][(j1-j)/scale+halfwin] */  rangefn_L[(int)(data_fine->L[i1][j1]-data_fine->L[i][j]+0x10000)] )

#define DIRWT_AB(i1,j1,i,j) ( /*domker[(i1-i)/scale+halfwin][(j1-j)/scale+halfwin]*/ rangefn_ab[(int)(data_fine->a[i1][j1]-data_fine->a[i][j]+0x10000)] *  \
rangefn_ab[(int)(data_fine->L[i1][j1]-data_fine->L[i][j]+0x10000)] * \
rangefn_ab[(int)(data_fine->b[i1][j1]-data_fine->b[i][j]+0x10000)] )

#define NRWT_L(a) (nrwt_l[a] )

#define NRWT_AB (nrwt_ab[(int)((hipass[1]+0x10000))] * nrwt_ab[(int)((hipass[2]+0x10000))])



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
	
	
	
	
	void ImProcFunctions :: dirpyrLab_denoise(LabImage * src, LabImage * dst, const int luma, const int chroma, float gam )
	{
		//float gam = 2.0;//MIN(3.0, 0.1*fabs(c[4])/3.0+0.001);
		float gamthresh = 0.03;
		float gamslope = exp(log((double)gamthresh)/gam)/gamthresh;
		unsigned short gamcurve[65536];
		for (int i=0; i<65536; i++) {
			int g = (int)(CurveFactory::gamma((double)i/65535.0, gam, gamthresh, gamslope, 1.0, 0.0) * 65535.0);
			//if (i<500)  printf("%d %d \n",i,g);
			gamcurve[i] = CLIP(g);
		}
		
		
		//#pragma omp parallel for if (multiThread)
		for (int i=0; i<src->H; i++) {
			for (int j=0; j<src->W; j++) {
				src->L[i][j] = gamcurve[src->L[i][j] ];
			}
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		int * rangefn_L = new int [0x20000];
		float * nrwt_l = new float [0x10000];
		
		
		int * rangefn_ab = new int [0x20000];
		float * nrwt_ab = new float [0x20000];
		
		int intfactor = 1024;//16384;
		
		
		//set up NR weight functions
		
		//gamma correction for chroma in shadows
		float nrwtl_norm = ((CurveFactory::gamma((double)65535.0/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)) - \
							(CurveFactory::gamma((double)75535.0/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)));
		for (int i=0; i<0x10000; i++) {
			nrwt_l[i] = ((CurveFactory::gamma((double)i/65535.0, gam, gamthresh, gamslope, 1.0, 0.0) - \
						  CurveFactory::gamma((double)(i+10000)/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)) )/nrwtl_norm;
			//if (i % 100 ==0) printf("%d %f \n",i,nrwt_l[i]);
		}
		
		float tonefactor = nrwt_l[32768];
		
		float noise_L = 25.0*luma;
		float noisevar_L = 4*SQR(noise_L);
		
		float noise_ab = 25*chroma;
		float noisevar_ab = SQR(noise_ab);
		
		
		//set up range functions
		for (int i=0; i<0x20000; i++) 
			rangefn_L[i] = (int)(( exp(-(double)fabs(i-0x10000) * tonefactor / (1+3*noise_L)) * noisevar_L/((double)(i-0x10000)*(double)(i-0x10000) + noisevar_L))*intfactor); 
		for (int i=0; i<0x20000; i++) 
			rangefn_ab[i] = (int)(( exp(-(double)fabs(i-0x10000) * tonefactor / (1+3*noise_ab)) * noisevar_ab/((double)(i-0x10000)*(double)(i-0x10000) + noisevar_ab))*intfactor); 
		
		
		for (int i=0; i<0x20000; i++) 
			nrwt_ab[i] = ((1+abs(i-0x10000)/(1+8*noise_ab)) * exp(-(double)fabs(i-0x10000)/ (1+8*noise_ab) ) );
		
		
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
		
		
		// c[0] = luma = noise_L
		// c[1] = chroma = noise_ab
		// c[2] decrease of noise var with scale
		// c[3] radius of domain blur at each level
		// c[4] shadow smoothing
		// c[5] edge preservation
		
		level = 0;
		
		int scale = scales[level];
		int pitch = pitches[level];
		//int thresh = 10 * c[8];
		//impulse_nr (src, src, m_w1, m_h1, thresh, noisevar);
		
		dirpyr(src, dirpyrLablo[0], 0, rangefn_L, rangefn_ab, pitch, scale, luma, chroma );
		
		level = 1;
		
		while(level < maxlevel)
		{
			scale = scales[level];
			pitch = pitches[level];
			
			dirpyr(dirpyrLablo[level-1], dirpyrLablo[level], level, rangefn_L, rangefn_ab, pitch, scale, luma, chroma );
			
			level ++;
		}
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		for(int level = maxlevel - 1; level > 0; level--)
		{
			
			int scale = scales[level];
			int pitch = pitches[level];
			idirpyr(dirpyrLablo[level], dirpyrLablo[level-1], level, nrwt_l, nrwt_ab, pitch, scale, luma, chroma );
		}
		
		
		scale = scales[0];
		pitch = pitches[0];
		idirpyr(dirpyrLablo[0], dst, 0, nrwt_l, nrwt_ab, pitch, scale, luma, chroma );
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		float igam = 1/gam;
		float igamthresh = gamthresh*gamslope;
		float igamslope = 1/gamslope;
		for (int i=0; i<65536; i++) {
			int g = (int)(CurveFactory::gamma((float)i/65535.0, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0);
			gamcurve[i] = CLIP(g);
		}
		
		
		for (int i=0; i<dst->H; i++) 
			for (int j=0; j<dst->W; j++) {
				
				dst->L[i][j] = gamcurve[CLIP(dst->L[i][j]) ];
				
			}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		for(int i = 0; i < maxlevel; i++)
		{
			delete dirpyrLablo[i];
		}
		
		delete [] rangefn_L;
		delete [] rangefn_ab;
		delete [] nrwt_l;
		delete [] nrwt_ab;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	};
	
	void ImProcFunctions::dirpyr(LabImage* data_fine, LabImage* data_coarse, int level, int * rangefn_L, int * rangefn_ab, int pitch, int scale, const int luma, const int chroma )
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
		//int intfactor = 16384;
		
		/*float domker[7][7];
		 for (int i=-halfwin; i<=halfwin; i++)
		 for (int j=-halfwin; j<=halfwin; j++) {
		 domker[i+halfwin][j+halfwin] = (int)(exp(-(i*i+j*j)/(2*sig*sig))*intfactor); //or should we use a value that depends on sigma???
		 }*/
		//float domker[5][5] = {{1,1,1,1,1},{1,2,2,2,1},{1,2,4,2,1},{1,2,2,2,1},{1,1,1,1,1}};
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
 
		for(int i = 0; i < height; i+=pitch ) { int i1=i/pitch;
			for(int j = 0, j1=0; j < width; j+=pitch, j1++)
			{				
				//norm = DIRWT(i, j, i, j);
				//Lout = -norm*data_fine->L[i][j];//if we don't want to include the input pixel in the sum
				//aout = -norm*data_fine->a[i][j];
				//bout = -norm*data_fine->b[i][j];
				//or
				float dirwt_l, dirwt_ab, norm_l, norm_ab;
				//float lops,aops,bops;
				float Lout, aout, bout;
				norm_l = norm_ab = 0;//if we do want to include the input pixel in the sum
				Lout = 0;
				aout = 0;
				bout = 0;
				//normab = 0;
				
				for(int inbr=MAX(0,i-scalewin); inbr<=MIN(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=MAX(0,j-scalewin); jnbr<=MIN(width-1,j+scalewin); jnbr+=scale) {
						dirwt_l = DIRWT_L(inbr, jnbr, i, j);
						dirwt_ab = DIRWT_AB(inbr, jnbr, i, j);
						Lout += dirwt_l*data_fine->L[inbr][jnbr];
						aout += dirwt_ab*data_fine->a[inbr][jnbr];
						bout += dirwt_ab*data_fine->b[inbr][jnbr];
						norm_l += dirwt_l;
						norm_ab += dirwt_ab;
					}
				}
				//lops = Lout/norm;//diagnostic
				//aops = aout/normab;//diagnostic
				//bops = bout/normab;//diagnostic
				
				//data_coarse->L[i1][j1]=0.5*(data_fine->L[i][j]+Lout/norm_l);//low pass filter
				//data_coarse->a[i1][j1]=0.5*(data_fine->a[i][j]+aout/norm_ab);
				//data_coarse->b[i1][j1]=0.5*(data_fine->b[i][j]+bout/norm_ab);
				//or
				data_coarse->L[i1][j1]=Lout/norm_l;//low pass filter
				data_coarse->a[i1][j1]=aout/norm_ab;
				data_coarse->b[i1][j1]=bout/norm_ab;
			}
		}
		
		
		
		
	};
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::idirpyr(LabImage* data_coarse, LabImage* data_fine, int level, float * nrwt_l, float * nrwt_ab, int pitch, int scale, const int luma, const int chroma )
	{
		
		int width = data_fine->W;
		int height = data_fine->H;
		
		//float eps = 0.0;
		
		// c[0] noise_L
		// c[1] noise_ab (relative to noise_L)
		// c[2] decrease of noise var with scale
		// c[3] radius of domain blur at each level
		// c[4] shadow smoothing
		
		float noisevar_L = 4*SQR(25.0 * luma);
		float noisevar_ab = 2*SQR(100.0 * chroma);
		float scalefactor = 1.0/pow(2.0,(level+1)*2);//change the last 2 to 1 for longer tail of higher scale NR
		//float recontrast = (1+((float)(c[6])/100.0));
		//float resaturate = 10*(1+((float)(c[7])/100.0));
		noisevar_L *= scalefactor;
		
		//int halfwin = 3;//MIN(ceil(2*sig),3);
		//int intfactor= 16384; 
		//int winwidth=1+2*halfwin;//this belongs in calling function
		/*float domker[7][7];
		 for (int i=-halfwin; i<=halfwin; i++)
		 for (int j=-halfwin; j<=halfwin; j++) {
		 domker[i][j] = (int)(exp(-(i*i+j*j)/(2*sig*sig))*intfactor); //or should we use a value that depends on sigma???
		 }*/
		//float domker[5][5] = {{1,1,1,1,1},{1,2,2,2,1},{1,2,4,2,1},{1,2,2,2,1},{1,1,1,1,1}};
		
		
		
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
#pragma omp parallel for
#endif
			for(int  i = 0; i < height; i++)
				for(int  j = 0; j < width; j++) {
					double wtdsum[3], norm;
					float hipass[3], hpffluct[3], tonefactor, nrfactor;
				
					tonefactor = ((NRWT_L(data_coarse->L[i][j])));
					
					//Wiener filter
					//luma
					if (level<2) {
						hipass[0] = data_fine->L[i][j]-data_coarse->L[i][j];
						hpffluct[0]=SQR(hipass[0])+0.001;
						hipass[0] *= hpffluct[0]/(hpffluct[0]+noisevar_L);
						data_fine->L[i][j] = CLIP(hipass[0]+data_coarse->L[i][j]);
					}
					
					//chroma
					hipass[1] = data_fine->a[i][j]-data_coarse->a[i][j];
					hipass[2] = data_fine->b[i][j]-data_coarse->b[i][j];
					hpffluct[1]=SQR(hipass[1]*tonefactor)+0.001;
					hpffluct[2]=SQR(hipass[2]*tonefactor)+0.001;
					nrfactor = (hpffluct[1]+hpffluct[2]) /((hpffluct[1]+hpffluct[2]) + noisevar_ab * NRWT_AB);

					hipass[1] *= nrfactor;
					hipass[2] *= nrfactor;
					
					data_fine->a[i][j] = hipass[1]+data_coarse->a[i][j];
					data_fine->b[i][j] = hipass[2]+data_coarse->b[i][j];
				}
			
		} else {
			
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
					
					double tonefactor = ((NRWT_L(smooth->L[i][j])));
					//double wtdsum[3], norm;
					float hipass[3], hpffluct[3],  nrfactor;
					//Wiener filter
					//luma
					if (level<2) {
						hipass[0] = data_fine->L[i][j]-smooth->L[i][j];
						hpffluct[0]=SQR(hipass[0])+0.001;
						hipass[0] *= hpffluct[0]/(hpffluct[0]+noisevar_L);
						data_fine->L[i][j] = CLIP(hipass[0]+smooth->L[i][j]);
					}
					
					//chroma
					hipass[1] = data_fine->a[i][j]-smooth->a[i][j];
					hipass[2] = data_fine->b[i][j]-smooth->b[i][j];
					hpffluct[1]=SQR(hipass[1]*tonefactor)+0.001;
					hpffluct[2]=SQR(hipass[2]*tonefactor)+0.001;
					nrfactor = (hpffluct[1]+hpffluct[2]) /((hpffluct[1]+hpffluct[2]) + noisevar_ab * NRWT_AB);

					hipass[1] *= nrfactor;
					hipass[2] *= nrfactor;
					
					data_fine->a[i][j] = hipass[1]+smooth->a[i][j];
					data_fine->b[i][j] = hipass[2]+smooth->b[i][j];
				}
}	// end parallel		
			delete smooth;
		}//end of pitch>1
		
	};
	
	
#undef DIRWT_L
#undef DIRWT_AB
	
#undef NRWT_L	
#undef NRWT_AB	
	
}

