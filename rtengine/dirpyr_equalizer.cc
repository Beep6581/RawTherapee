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
#include <rawimagesource.h>
#include <array2D.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SQR(x) ((x)*(x))
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define CLIP(a) (CLIPTO(a,0,65535))



#define DIRWT(i1,j1,i,j) ( domker[(i1-i)/scale+halfwin][(j1-j)/scale+halfwin] * rangefn[abs((int)data_fine[i1][j1]-data_fine[i][j])] )

namespace rtengine {
	
	static const int maxlevel = 4;
	static const float noise = 2000;
	static const float thresh = 1000;
	
	//sequence of scales
	static const int scales[8] = {1,2,4,8,16,32,64,128};
	
	//sequence of scales
	//static const int scales[8] = {1,2,3,6,15,21,28,36};
	//scale is spacing of directional averaging weights
	
	
	void ImProcFunctions :: dirpyr_equalizer(float ** src, float ** dst, int srcwidth, int srcheight, const double * mult )
	{
		int lastlevel=maxlevel;
		
		while (fabs(mult[lastlevel-1]-1)<0.001 && lastlevel>0) {
			lastlevel--;
			//printf("last level to process %d \n",lastlevel);
		}
		if (lastlevel==0) return;
		

		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		LUTf rangefn(0x10000);
				
		int intfactor = 1024;//16384;
		
		
		//set up range functions
		
		for (int i=0; i<0x10000; i++) {
			rangefn[i] = (int)((thresh/((double)(i) + thresh))*intfactor);
		}
				
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		int level;
		array2D<float> buffer (srcwidth, srcheight);
		
		multi_array2D<float,maxlevel> dirpyrlo (srcwidth, srcheight);

				
		for (int i=0; i<srcheight; i++)
			for (int j=0; j<srcwidth; j++) {
				buffer[i][j]=0;
			}
		
		level = 0;
		
		int scale = scales[level];
		//int thresh = 100 * mult[5];
				
		dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, rangefn, 0, scale, mult );
		
		level = 1;
		
		while(level < lastlevel)
		{
			scale = scales[level];
						
			dirpyr_channel(dirpyrlo[level-1], dirpyrlo[level], srcwidth, srcheight, rangefn, level, scale, mult );
			
			level ++;
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//initiate buffer for final image
		for(int i = 0; i < srcheight; i++)
			for(int j = 0; j < srcwidth; j++) {
				
				//copy pixels
				buffer[i][j] = dirpyrlo[lastlevel-1][i][j];
				
			}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		for(int level = lastlevel - 1; level > 0; level--)
		{
			idirpyr_eq_channel(dirpyrlo[level], dirpyrlo[level-1], buffer, srcwidth, srcheight, level, mult );
		}
		
		
		scale = scales[0];
		
		idirpyr_eq_channel(dirpyrlo[0], dst, buffer, srcwidth, srcheight, 0, mult );
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		for (int i=0; i<srcheight; i++) 
			for (int j=0; j<srcwidth; j++) {
				
				dst[i][j] = CLIP((int)(  buffer[i][j]  ));  // TODO: Really a clip necessary?
								
			}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
				
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}
	
	void ImProcFunctions::dirpyr_channel(float ** data_fine, float ** data_coarse, int width, int height, LUTf & rangefn, int level, int scale, const double * mult  )
	{
		//scale is spacing of directional averaging weights
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// calculate weights, compute directionally weighted average
		
		int halfwin=2;
		int domker[5][5] = {{1,1,1,1,1},{1,2,2,2,1},{1,2,2,2,1},{1,2,2,2,1},{1,1,1,1,1}};
		
		//generate domain kernel 
		if (level<2) {
			halfwin = 1;
			domker[1][1]=domker[1][2]=domker[2][1]=domker[2][2]=1;
		}
		
		
		int scalewin = halfwin*scale;
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(int i = 0; i < height; i++) {
			for(int j = 0; j < width; j++)
			{
				float val=0;
				float norm=0;
				
				for(int inbr=MAX(0,i-scalewin); inbr<=MIN(height-1,i+scalewin); inbr+=scale) {
					for (int jnbr=MAX(0,j-scalewin); jnbr<=MIN(width-1,j+scalewin); jnbr+=scale) {
						float dirwt = DIRWT(inbr, jnbr, i, j);
						val += dirwt*data_fine[inbr][jnbr];
						norm += dirwt;
					}
				}
				data_coarse[i][j]=val/norm;//low pass filter
			}
		}
		
		
		
		
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::idirpyr_eq_channel(float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, const double * mult )
	{
		float noisehi = 1.33*noise*mult[4]/pow(3,level), noiselo = 0.66*noise*mult[4]/pow(3,level);
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
				register float hipass = (data_fine[i][j]-data_coarse[i][j]);
				buffer[i][j] += irangefn[hipass+0x10000] * hipass ;
			}
		}
		
	}
	
	
#undef DIRWT_L
#undef DIRWT_AB
	
#undef NRWT_L	
#undef NRWT_AB	
	
}

