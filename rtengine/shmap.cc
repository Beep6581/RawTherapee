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
#include "rtengine.h"
#include "rt_math.h"
#include "rawimagesource.h"
#include "sleef.c"
#undef THREAD_PRIORITY_NORMAL
#ifdef __SSE2__
#include "sleefsseavx.c"
#endif // __SSE2__

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

    // fill with luminance
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
            map[i][j] = lumi[0]*std::max(img->r(i,j),0.f) + lumi[1]*std::max(img->g(i,j),0.f) + lumi[2]*std::max(img->b(i,j),0.f);
		}

    if (!hq) {
#ifdef _OPENMP
#pragma omp parallel
#endif
{
        AlignedBufferMP<double>* pBuffer = new AlignedBufferMP<double> (max(W,H));
    	gaussHorizontal<float> (map, map, *pBuffer, W, H, radius);
		gaussVertical<float>   (map, map, *pBuffer, W, H, radius);
        delete pBuffer;
}
    }

    else {
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//experimental dirpyr shmap

		float thresh = 100*radius;//1000;
		LUTf rangefn(0x10000);
		float ** dirpyrlo[2];

		int intfactor = 1024;//16384;

		//set up range functions
		for (int i=0; i<0x10000; i++) {
			//rangefn[i] = (int)(((thresh)/((double)(i) + (thresh)))*intfactor);
			rangefn[i] = static_cast<int>(xexpf(-(min(10.0f,(static_cast<float>(i)*i) / (thresh*thresh))))*intfactor);
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
    // update average, minimum, maximum

    float _avg = 0.0f;
    min_f = 65535;
    max_f = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    float _min_f = 65535.0f;
    float _max_f = 0.0f;
    float _val;
#ifdef _OPENMP
#pragma omp for reduction(+:_avg) nowait
#endif
    for (int i=32; i<H-32; i++)
        for (int j=32; j<W-32; j++) {
            _val = map[i][j];
            if (_val < _min_f)
                _min_f = _val;
            if (_val > _max_f)
                _max_f = _val;
            _avg += _val;
        }
#ifdef _OPENMP
#pragma omp critical
#endif
{
	if(_min_f < min_f )
		min_f = _min_f;
	if(_max_f > max_f )
		max_f = _max_f;
}
}
    _avg /= ((H-64)*(W-64));
    avg = _avg;

}

void SHMap::forceStat (float max_, float min_, float avg_) {

    max_f = max_;
    min_f = min_;
    avg = avg_;
}

#if defined( __SSE__ ) && defined( WIN32 )
__attribute__((force_align_arg_pointer)) void SHMap::dirpyr_shmap(float ** data_fine, float ** data_coarse, int width, int height, LUTf & rangefn, int level, int scale)
#else
void SHMap::dirpyr_shmap(float ** data_fine, float ** data_coarse, int width, int height, LUTf & rangefn, int level, int scale)
#endif
{
	//scale is spacing of directional averaging weights

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// calculate weights, compute directionally weighted average

	int scalewin, halfwin;

	if(level < 2) {
		halfwin = 1;
		scalewin = halfwin*scale;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
#if defined( __SSE2__ ) && defined( __x86_64__ )
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

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j%scale; jnbr<=j+scalewin; jnbr+=scale) {
					dirwt = ( rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}
#if defined( __SSE2__ ) && defined( __x86_64__ )
		for(; j < (width-scalewin)-3; j+=4)
		{
			valv= _mm_setzero_ps();
			normv= _mm_setzero_ps();

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
					dirwtv = ( rangefn[_mm_cvttps_epi32(vabsf(LVFU(data_fine[inbr][jnbr])-LVFU(data_fine[i][j])))] );
					valv += dirwtv*LVFU(data_fine[inbr][jnbr]);
					normv += dirwtv;
				}
			}
			_mm_storeu_ps( &data_coarse[i][j], valv/normv);
		}
		for(; j < width-scalewin; j++)
		{
			float val=0;
			float norm=0;

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
					dirwt = ( rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}

#else
		for(; j < width-scalewin; j++)
		{
			float val=0;
			float norm=0;

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
					dirwt = ( rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}
#endif
		for(; j < width; j++)
		{
			float val=0;
			float norm=0;

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<width; jnbr+=scale) {
					dirwt = ( rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}
	}
}
}
else {
	halfwin=2;
	scalewin = halfwin*scale;
	int domker[5][5] = {{1,1,1,1,1},{1,2,2,2,1},{1,2,2,2,1},{1,2,2,2,1},{1,1,1,1,1}};
	//generate domain kernel

#ifdef _OPENMP
#pragma omp parallel
#endif
{
#if defined( __SSE2__ ) && defined( __x86_64__ )
	__m128 dirwtv, valv, normv;
	float domkerv[5][5][4] __attribute__ ((aligned (16))) = {{{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}},{{1,1,1,1},{2,2,2,2},{2,2,2,2},{2,2,2,2},{1,1,1,1}},{{1,1,1,1},{2,2,2,2},{2,2,2,2},{2,2,2,2},{1,1,1,1}},{{1,1,1,1},{2,2,2,2},{2,2,2,2},{2,2,2,2},{1,1,1,1}},{{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}}};

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

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j%scale; jnbr<=j+scalewin; jnbr+=scale) {
					dirwt = ( domker[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin] * rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}
#if defined( __SSE2__ ) && defined( __x86_64__ )
		for(; j < width-scalewin-3; j+=4)
		{
			valv = _mm_setzero_ps();
			normv = _mm_setzero_ps();

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
					dirwtv = ( _mm_load_ps((float*)&domkerv[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin]) * rangefn[_mm_cvttps_epi32(vabsf(LVFU(data_fine[inbr][jnbr])-LVFU(data_fine[i][j])))] );
					valv += dirwtv*LVFU(data_fine[inbr][jnbr]);
					normv += dirwtv;
				}
			}
			_mm_storeu_ps( &data_coarse[i][j], valv/normv);
		}
		for(; j < width-scalewin; j++)
		{
			float val=0;
			float norm=0;

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
					dirwt = ( domker[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin] * rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}

#else
		for(; j < width-scalewin; j++)
		{
			float val=0;
			float norm=0;

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<=j+scalewin; jnbr+=scale) {
					dirwt = ( domker[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin] * rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}
#endif
		for(; j < width; j++)
		{
			float val=0;
			float norm=0;

			for(int inbr=max(i-scalewin,i%scale); inbr<=min(i+scalewin, height-1); inbr+=scale) {
				for (int jnbr=j-scalewin; jnbr<width; jnbr+=scale) {
					dirwt = ( domker[(inbr-i)/scale+halfwin][(jnbr-j)/scale+halfwin] * rangefn[abs(data_fine[inbr][jnbr]-data_fine[i][j])] );
					val += dirwt*data_fine[inbr][jnbr];
					norm += dirwt;
				}
			}
			data_coarse[i][j] = val/norm; // low pass filter
		}
	}
}

}

}

}//end of SHMap
