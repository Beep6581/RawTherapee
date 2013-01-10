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
#include "shmap.h"
#include "gauss.h"
#include "bilateral2.h"
#include "rtengine.h"
#include "rt_math.h"
#include "rawimagesource.h"

#undef THREAD_PRIORITY_NORMAL

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

void SHMap::update (Imagefloat* img, double radius, double lumi[3], bool hq, int skip) {
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    // fill with luminance
    #pragma omp for
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
            map[i][j] = lumi[0]*std::max(img->r(i,j),0.f) + lumi[1]*std::max(img->g(i,j),0.f) + lumi[2]*std::max(img->b(i,j),0.f);
		}

    if (!hq) {
        AlignedBufferMP<double>* pBuffer = new AlignedBufferMP<double> (max(W,H));
    	gaussHorizontal<float> (map, map, *pBuffer, W, H, radius);
		gaussVertical<float>   (map, map, *pBuffer, W, H, radius);
        delete pBuffer;
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
			rangefn[i] = static_cast<int>(exp(-(min(10.0f,(static_cast<float>(i)*i) / (thresh*thresh))))*intfactor);
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
    min_f = 65535;
    max_f = 0;
    for (int i=32; i<H-32; i++)
        for (int j=32; j<W-32; j++) {
            int val = map[i][j];
            if (val < min_f)
                min_f = val;
            if (val > max_f)
                max_f = val;
            _avg = 1.0/n * val + (1.0 - 1.0/n) * _avg;
            n++;
        }
    avg = (int) _avg;
}

void SHMap::forceStat (float max_, float min_, float avg_) {

    max_f = max_;
    min_f = min_;
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
			data_coarse[i][j] = val/norm; // low pass filter
			/*if (val<=0 || norm<=0)
				printf("val=%f norm=%f \n",val,norm); */
		}
	}
	
}
	
	
}//end of SHMap
