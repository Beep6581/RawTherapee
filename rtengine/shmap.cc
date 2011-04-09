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
#include <shmap.h>
#include <gauss.h>
#include <bilateral2.h>
#include <rtengine.h>

#include <rawimagesource.h>//for dirpyr


#undef MAXVAL
#undef CLIP
#undef MAX
#undef MIN
#undef SQR

#undef THREAD_PRIORITY_NORMAL
#define MAXVAL  0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)
#define SQR(x) ((x)*(x))

#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))

namespace rtengine {

extern const Settings* settings;

SHMap::SHMap (int w, int h, bool multiThread) : W(w), H(h), multiThread(multiThread) {

    map = new float*[H];
    for (int i=0; i<H; i++)
        map[i] = new float[W];
}

SHMap::~SHMap () {

    for (int i=0; i<H; i++)
        delete [] map[i];
    delete [] map;
}

void SHMap::update (Imagefloat* img, float** buffer, double radius, double lumi[3], bool hq, int skip) {

    // fill with luminance
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
			int val = lumi[0]*img->r[i][j] + lumi[1]*img->g[i][j] + lumi[2]*img->b[i][j];
			map[i][j] = CLIP(val);
		}
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    if (!hq) {
    	AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(W,H));
    	gaussHorizontal<float> (map, map, buffer, W, H, radius, multiThread);
		gaussVertical<float>   (map, map, buffer, W, H, radius, multiThread);

        delete buffer;
    }
    else {
/*		
#if 0
// the new OpenMP method does not need thread number specific code.
//    	#ifdef _OPENMP
		#pragma omp parallel if (multiThread)
    	{
    		int tid = omp_get_thread_num();
    		int nthreads = omp_get_num_threads();
    		int blk = H/nthreads;

    		if (tid<nthreads-1)
    			bilateral<float> (map, buffer, W, H, 8000, radius, tid*blk, (tid+1)*blk);
    		else
    			bilateral<float> (map, buffer, W, H, 8000, radius, tid*blk, H);
		}
#else
    	bilateral<float> (map, buffer, W, H, 8000, radius, 0, H);
#endif
*/
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//experimental dirpyr shmap
		
		float thresh = 100*radius;//1000;
		LUTf rangefn(0x10000);
		float ** dirpyrlo[2];

		int intfactor = 1024;//16384;
		
		//set up range functions
		
		for (int i=0; i<0x10000; i++) {
			//rangefn[i] = (int)(((thresh)/((double)(i) + (thresh)))*intfactor);
			rangefn[i] = (int)(exp(-(MIN(10,((float)i*i) / (thresh*thresh))))*intfactor);
			//if (rangefn[i]<0 || rangefn[i]>intfactor) 
				//printf("i=%d rangefn=%d arg=%f \n",i,rangefn[i], float(i*i) / (thresh*thresh));
		}
		dirpyrlo[0] = allocArray<float> (W, H);
		dirpyrlo[1] = allocArray<float> (W, H);

		int scale=1;
		int level=0;
		int indx=0;
		dirpyr_shmap(map, dirpyrlo[indx], W, H, rangefn, 0, scale );
		scale *= 2;
		level += 1;
		indx = 1-indx;
		while (skip*scale<16) {
			dirpyr_shmap(dirpyrlo[1-indx], dirpyrlo[indx], W, H, rangefn, level, scale );
			scale *= 2;
			level += 1;
			indx = 1-indx;
		}
		/*dirpyr_shmap(dirpyrlo[0], dirpyrlo[1], W, H, rangefn, 1, scale );
		scale = 4;
		dirpyr_shmap(dirpyrlo[1], dirpyrlo[0], W, H, rangefn, 2, scale );
		scale = 8;
		dirpyr_shmap(dirpyrlo[0], dirpyrlo[1], W, H, rangefn, 3, scale );
		scale = 16;
		dirpyr_shmap(dirpyrlo[1], map, W, H, rangefn, 3, scale );*/

		dirpyr_shmap(dirpyrlo[1-indx], map, W, H, rangefn, level, scale );
		

		freeArray<float>(dirpyrlo[0], H);
		freeArray<float>(dirpyrlo[1], H);
		

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
        // anti-alias filtering the result
#ifdef _OPENMP
#pragma omp for
#endif
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++)
                if (i>0 && j>0 && i<H-1 && j<W-1)
                    map[i][j] = (buffer[i-1][j-1]+buffer[i-1][j]+buffer[i-1][j+1]+buffer[i][j-1]+buffer[i][j]+buffer[i][j+1]+buffer[i+1][j-1]+buffer[i+1][j]+buffer[i+1][j+1])/9;
                else
                    map[i][j] = buffer[i][j];
*/		
		
    }
    } // end parallel enclosure
    // update average, minimum, maximum
    double _avg = 0;
    int n = 1;
    min = 65535;
    max = 0;
    for (int i=32; i<H-32; i++)
        for (int j=32; j<W-32; j++) {
            int val = map[i][j];
            if (val < min)
                min = val;
            if (val > max)
                max = val;
            _avg = 1.0/n * val + (1.0 - 1.0/n) * _avg;
            n++;
        }
    avg = (int) _avg;
}

void SHMap::forceStat (float max_, float min_, float avg_) {

    max = max_;
    min = min_;
    avg = avg_;
}


void SHMap::dirpyr_shmap(float ** data_fine, float ** data_coarse, int width, int height, LUTf & rangefn, int level, int scale)
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
			
			for(int inbr=(i-scalewin); inbr<=(i+scalewin); inbr+=scale) {
				if (inbr<0 || inbr>height-1) continue;
				for (int jnbr=(j-scalewin); jnbr<=(j+scalewin); jnbr+=scale) {
					if (jnbr<0 || jnbr>width-1) continue;
					float dirwt = ( domker[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin] * rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
					/*if (val<0 || norm<0) {
						printf("val=%f norm=%f \n",val,norm);
						printf("i=%d j=%d inbr=%d jnbr=%d domker=%d val=%d nbrval=%d rangefn=%d \n",i,j,inbr,jnbr, \
							   domker[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin], \
							   data_fine[i][j], data_fine[inbr][jnbr], \
							   rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])]);
					}*/
				}
			}
			data_coarse[i][j]=CLIP((int)(val/norm));//low pass filter
			if (val<=0 || norm<=0)
				printf("val=%f norm=%f \n",val,norm);
		}
	}
	
}
	
	
}//end of SHMap
