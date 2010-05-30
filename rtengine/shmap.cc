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

#undef THREAD_PRIORITY_NORMAL
#define MAXVAL  0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)

namespace rtengine {

extern const Settings* settings;

SHMap::SHMap (int w, int h, bool multiThread) : W(w), H(h), multiThread(multiThread) {

    map = new unsigned short*[H];
    for (int i=0; i<H; i++)
        map[i] = new unsigned short[W];
}

SHMap::~SHMap () {

    for (int i=0; i<H; i++)
        delete [] map[i];
    delete [] map;
}

void SHMap::update (Image16* img, unsigned short** buffer, double radius, double lumi[3], bool hq) {

    // fill with luminance
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
			int val = lumi[0]*img->r[i][j] + lumi[1]*img->g[i][j] + lumi[2]*img->b[i][j];
			map[i][j] = CLIP(val);
		}

    if (!hq) {
#ifdef _OPENMP
    	AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(W,H)*omp_get_max_threads());
#else
    	AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(W,H));
#endif
    	gaussHorizontal<unsigned short> (map, map, buffer, W, H, radius, multiThread);
		gaussVertical<unsigned short>   (map, map, buffer, W, H, radius, multiThread);

        delete buffer;
    }
    else {
#ifdef _OPENMP
		#pragma omp parallel if (multiThread)
    	{
    		int tid = omp_get_thread_num();
    		int nthreads = omp_get_num_threads();
    		int blk = H/nthreads;

    		if (tid<nthreads-1)
    			bilateral<unsigned short> (map, buffer, W, H, 8000, radius, tid*blk, (tid+1)*blk);
    		else
    			bilateral<unsigned short> (map, buffer, W, H, 8000, radius, tid*blk, H);
		}
#else
    	bilateral<unsigned short> (map, buffer, W, H, 8000, radius, 0, H);
#endif
        // anti-alias filtering the result
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++)
                if (i>0 && j>0 && i<H-1 && j<W-1)
                    map[i][j] = (buffer[i-1][j-1]+buffer[i-1][j]+buffer[i-1][j+1]+buffer[i][j-1]+buffer[i][j]+buffer[i][j+1]+buffer[i+1][j-1]+buffer[i+1][j]+buffer[i+1][j+1])/9;
                else
                    map[i][j] = buffer[i][j];
    }

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

void SHMap::forceStat (unsigned short max_, unsigned short min_, unsigned short avg_) {

    max = max_;
    min = min_;
    avg = avg_;
}}

