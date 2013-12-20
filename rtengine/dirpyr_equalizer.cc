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

#include <cstddef>
#include <cmath>
#include "curves.h"
#include "labimage.h"
#include "improcfun.h"
#include "rawimagesource.h"
#include "array2D.h"
#include "rt_math.h"
#ifdef __SSE2__
#include "sleefsseavx.c"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#define CLIPI(a) ((a)>0 ?((a)<32768 ?(a):32768):0)

#define RANGEFN(i) ((1000.0f / (i + 1000.0f)))
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define DIRWT(i1,j1,i,j) ( domker[(i1-i)/scale+halfwin][(j1-j)/scale+halfwin] * RANGEFN(fabsf((data_fine[i1][j1]-data_fine[i][j]))) )

namespace rtengine {
	
	static const int maxlevel = 5;
	static const float noise = 2000;
	static const float thresh = 1000;
	
	//sequence of scales
	static const int scales[8] = {1,2,4,8,16,32,64,128};
	
	//sequence of scales
	//static const int scales[8] = {1,2,3,6,15,21,28,36};
	//scale is spacing of directional averaging weights
	
	
	void ImProcFunctions :: dirpyr_equalizer(float ** src, float ** dst, int srcwidth, int srcheight, const double * mult, const double dirpyrThreshold )
	{
		int lastlevel=maxlevel;
		
		while (fabs(mult[lastlevel-1]-1)<0.001 && lastlevel>0) {
			lastlevel--;
			//printf("last level to process %d \n",lastlevel);
		}
		if (lastlevel==0) return;
		
		int level;
		
		multi_array2D<float,maxlevel> dirpyrlo (srcwidth, srcheight);

		level = 0;
		
		int scale = scales[level];
		//int thresh = 100 * mult[5];
				
		dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, 0, scale );
		
		level = 1;
		
		while(level < lastlevel)
		{
			scale = scales[level];
						
			dirpyr_channel(dirpyrlo[level-1], dirpyrlo[level], srcwidth, srcheight, level, scale );
			
			level ++;
		}
		
		// with the current implementation of idirpyr_eq_channel we can safely use the buffer from last level as buffer, saves some memory
		float ** buffer = dirpyrlo[lastlevel-1];
		
		for(int level = lastlevel - 1; level > 0; level--)
		{
			idirpyr_eq_channel(dirpyrlo[level], dirpyrlo[level-1], buffer, srcwidth, srcheight, level, mult, dirpyrThreshold );
		}
		
		
		scale = scales[0];
		
		idirpyr_eq_channel(dirpyrlo[0], dst, buffer, srcwidth, srcheight, 0, mult, dirpyrThreshold );
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int i=0; i<srcheight; i++) 
			for (int j=0; j<srcwidth; j++) {
				dst[i][j] = CLIP(  buffer[i][j] );  // TODO: Really a clip necessary?
								
			}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
				
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}


	
	void ImProcFunctions :: dirpyr_equalizercam (CieImage *ncie, float ** src, float ** dst, int srcwidth, int srcheight, const double * mult, const double dirpyrThreshold, bool execdir )
	{
		int lastlevel=maxlevel;
		
		while (fabs(mult[lastlevel-1]-1)<0.001 && lastlevel>0) {
			lastlevel--;
			//printf("last level to process %d \n",lastlevel);
		}
		if (lastlevel==0) return;
		
		int level;
	
		multi_array2D<float,maxlevel> dirpyrlo (srcwidth, srcheight);

		level = 0;
		
		int scale = scales[level];
		//int thresh = 100 * mult[5];
				
		dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, 0, scale );
		
		level = 1;
		
		while(level < lastlevel)
		{
			scale = scales[level];
						
			dirpyr_channel(dirpyrlo[level-1], dirpyrlo[level], srcwidth, srcheight, level, scale );
			
			level ++;
		}
		
		
		// with the current implementation of idirpyr_eq_channel we can safely use the buffer from last level as buffer, saves some memory
		float ** buffer = dirpyrlo[lastlevel-1];
		
		for(int level = lastlevel - 1; level > 0; level--)
		{
			idirpyr_eq_channel(dirpyrlo[level], dirpyrlo[level-1], buffer, srcwidth, srcheight, level, mult, dirpyrThreshold );
		}
		
		
		scale = scales[0];
		
		idirpyr_eq_channel(dirpyrlo[0], dst, buffer, srcwidth, srcheight, 0, mult, dirpyrThreshold );
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(execdir)
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (int i=0; i<srcheight; i++) 
				for (int j=0; j<srcwidth; j++) {
					if(ncie->J_p[i][j] > 8.f && ncie->J_p[i][j] < 92.f)
						dst[i][j] = CLIP( buffer[i][j] );  // TODO: Really a clip necessary?
					else
						dst[i][j]=src[i][j];
				}
		else
			for (int i=0; i<srcheight; i++) 
				for (int j=0; j<srcwidth; j++) {
					dst[i][j] = CLIP( buffer[i][j] );  // TODO: Really a clip necessary?
				}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}


#if defined( __SSE2__ ) && defined( WIN32 )
__attribute__((force_align_arg_pointer)) void ImProcFunctions::dirpyr_channel(float ** data_fine, float ** data_coarse, int width, int height, int level, int scale  )
#else
void ImProcFunctions::dirpyr_channel(float ** data_fine, float ** data_coarse, int width, int height, int level, int scale )
#endif
{
		//scale is spacing of directional averaging weights
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// calculate weights, compute directionally weighted average
		
	int halfwin;
	int scalewin;
		
	if(level > 1) {
		//generate domain kernel 
		int domker[5][5] = {{1,1,1,1,1},{1,2,2,2,1},{1,2,2,2,1},{1,2,2,2,1},{1,1,1,1,1}};
		halfwin=2;
		scalewin = halfwin*scale;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef __SSE2__
	__m128 thousandv = _mm_set1_ps( 1000.0f );
	__m128 dirwtv, valv, normv;
	float domkerv[5][5][4] = {{{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}},{{1,1,1,1},{2,2,2,2},{2,2,2,2},{2,2,2,2},{1,1,1,1}},{{1,1,1,1},{2,2,2,2},{2,2,2,2},{2,2,2,2},{1,1,1,1}},{{1,1,1,1},{2,2,2,2},{2,2,2,2},{2,2,2,2},{1,1,1,1}},{{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}}};
#endif // __SSE2__
	int j;
#ifdef _OPENMP
#pragma omp for
#endif
		for(int i = 0; i < height; i++) {
			float dirwt;
			for(j = 0; j < scalewin; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=max(0,j-scalewin); jnbr<=j+scalewin; jnbr+=scale) {
						dirwt = DIRWT(inbr, jnbr, i, j);
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
#ifdef __SSE2__
			for(; j < width-scalewin-3; j+=4)
			{
				valv = _mm_setzero_ps();
				normv = _mm_setzero_ps();
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
						dirwtv = _mm_loadu_ps((float*)&domkerv[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin]) * (thousandv / (vabsf(LVFU(data_fine[inbr][jnbr])-(LVFU(data_fine[i][j]))) + thousandv));
						valv += dirwtv*LVFU(data_fine[inbr][jnbr]);
						normv += dirwtv;
					}
				}
				_mm_storeu_ps( &data_coarse[i][j],valv/normv);//low pass filter
			}
			for(; j < width-scalewin; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
						dirwt = DIRWT(inbr, jnbr, i, j);
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
#else
			for(; j < width-scalewin; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
						dirwt = DIRWT(inbr, jnbr, i, j);
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
#endif
			for(; j < width; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=min(width-1,j+scalewin); jnbr+=scale) {
						dirwt = DIRWT(inbr, jnbr, i, j);
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
		}
}
	} else {	// level <=1 means that all values of domker would be 1.0f, so no need for multiplication
		halfwin = 1;
		scalewin = halfwin*scale;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef __SSE2__
	__m128 thousandv = _mm_set1_ps( 1000.0f );
	__m128 dirwtv, valv, normv;
#endif // __SSE2__
	int j;
#ifdef _OPENMP
#pragma omp for
#endif
		for(int i = 0; i < height; i++) {
			float dirwt;
			for(j = 0; j < scalewin; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=max(0,j-scalewin); jnbr<=j+scalewin; jnbr+=scale) {
						dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr]-data_fine[i][j]));
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
#ifdef __SSE2__
			for(; j < width-scalewin-3; j+=4)
			{
				valv = _mm_setzero_ps();
				normv = _mm_setzero_ps();
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
						dirwtv = thousandv / (vabsf(LVFU(data_fine[inbr][jnbr])-(LVFU(data_fine[i][j]))) + thousandv);
						valv += dirwtv*LVFU(data_fine[inbr][jnbr]);
						normv += dirwtv;
					}
				}
				_mm_storeu_ps( &data_coarse[i][j], valv/normv);//low pass filter
			}

			for(; j < width-scalewin; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
						dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr]-data_fine[i][j]));
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
#else
			for(; j < width-scalewin; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
						dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr]-data_fine[i][j]));
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
#endif
			for(; j < width; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=max(0,i-scalewin); inbr<=min(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=j-scalewin; jnbr<=min(width-1,j+scalewin); jnbr+=scale) {
						dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr]-data_fine[i][j]));
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
		}
}
	}
}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::idirpyr_eq_channel(float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, const double * mult, const double dirpyrThreshold )
	{
		float noisehi = 1.33*noise*dirpyrThreshold/expf(level*log(3.0)), noiselo = 0.66*noise*dirpyrThreshold/expf(level*log(3.0));
		LUTf irangefn (0x20000);

		for (int i=0; i<0x20000; i++) {
			if (abs(i-0x10000)>noisehi || mult[level]<1.0) {
				irangefn[i] = mult[level] ;
			} else {
				if (abs(i-0x10000)<noiselo) {
					irangefn[i] = 1.f ;
				} else {
					irangefn[i] = 1.f + (mult[level]-1) * (noisehi-abs(i-0x10000))/(noisehi-noiselo+0.01) ;
				}
			}
		}
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(int i = 0; i < height; i++) {
			for(int j = 0; j < width; j++) {
				float hipass = (data_fine[i][j]-data_coarse[i][j]);
				buffer[i][j] += irangefn[hipass+0x10000] * hipass ;
			}
		}
		
	}
	
	
#undef DIRWT_L
#undef DIRWT_AB
	
#undef NRWT_L	
#undef NRWT_AB	
	
}

