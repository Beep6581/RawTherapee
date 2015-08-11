
////////////////////////////////////////////////////////////////
//
//		//exposure correction before interpolation
//
//  code dated: December 27, 2010
//
//	Expo_before.cc is free software: you can redistribute it and/or modify
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



//		     Jacques Desmis <jdesmis@gmail.com>
//   	     use fast-demo(provisional) from Emil Martinec
//			 inspired from work Guillermo Luijk and Manuel LLorens(Perfectraw)
//			 Ingo Weyrich (2014-07-07)
//			    optimized the highlight protection case
//				needs 2*width*height*sizeof(float) byte less memory than before
//			    needs about 60% less processing time than before
//
// This function uses parameters:
//      exposure (linear): 2^(-8..0..8): currently 0.5 +3
//      preserve (log)   : 0..8 : currently 0.1 1

#include "rtengine.h"
#include "rawimagesource.h"
#include "mytime.h"
#include "rt_math.h"

namespace rtengine {

extern const Settings* settings;

void RawImageSource::processRawWhitepoint(float expos, float preser) {
	MyTime t1e,t2e;
	if (settings->verbose)
		t1e.set();
	
	int width=W, height=H;
    // exposure correction inspired from G.Luijk

	for (int c=0; c<4; c++)
		chmax[c] *= expos;

    if (fabs(preser)<0.001f) {
        // No highlight protection - simple mutiplication

        if (ri->getSensorType()==ST_BAYER || ri->getSensorType()==ST_FUJI_XTRANS)
			#pragma omp parallel for
			for (int row=0;row<height;row++)
				for (int col=0;col<width;col++)
					rawData[row][col] *= expos;
		else
			#pragma omp parallel for
			for (int row=0;row<height;row++)
				for (int col=0;col<width;col++) {
					rawData[row][col*3] *= expos;
					rawData[row][col*3+1] *= expos;
					rawData[row][col*3+2] *= expos;
				}
    } else {
        if (ri->getSensorType()==ST_BAYER || ri->getSensorType()==ST_FUJI_XTRANS) {
			// Demosaic to allow calculation of luminosity.
			if(ri->getSensorType()==ST_BAYER)
				fast_demosaic (0,0,W,H);
			else
				fast_xtrans_interpolate();
        }
		
        // Find maximum to adjust LUTs. New float engines clips only at the very end
        float maxValFloat = 0.f;
#pragma omp parallel
{
		float maxValFloatThr = 0.f;
        if (ri->getSensorType()==ST_BAYER || ri->getSensorType()==ST_FUJI_XTRANS)
#pragma omp for schedule(dynamic,16) nowait
			for(int row=0;row<height;row++)
				for (int col=0;col<width;col++) {
					if (rawData[row][col]>maxValFloatThr)
						maxValFloatThr = rawData[row][col];
				}
		else
#pragma omp for schedule(dynamic,16) nowait
			for(int row=0;row<height;row++)
				for (int col=0;col<width;col++) {
                    for (int c=0;c<3;c++)
						if (rawData[row][col*3+c]>maxValFloatThr)
							maxValFloatThr = rawData[row][col*3+c];
                }

#pragma omp critical
{
		if(maxValFloatThr > maxValFloat)
			maxValFloat = maxValFloatThr;
}
}

		// Exposure correction with highlight preservation
		int maxVal = maxValFloat;
        LUTf lut(maxVal+1);
        float K;
		if(expos>1){
            // Positive exposure
            K = (float) maxVal / expos*exp(-preser*log(2.0));
            for (int j=max(1,(int)K);j<=maxVal;j++) 
                lut[(int)j]=(((float)maxVal-K*expos)/((float)maxVal-K)*(j-maxVal)+(float) maxVal) / j;
        } else {
            // Negative exposure
			float EV=log(expos)/log(2.0);                              // Convert exp. linear to EV
            K = (float)maxVal * exp(-preser * log(2.0));
			
            for (int j=0;j<=maxVal;j++) 
                lut[(int)j] = exp(EV*((float)maxVal-j) / ((float)maxVal-K) * log(2.0));
        }	

		if (ri->getSensorType()==ST_BAYER || ri->getSensorType()==ST_FUJI_XTRANS)
			#pragma omp parallel for schedule(dynamic,16)
			for(int row=0;row<height;row++)
				for(int col=0;col<width;col++){
					float lumi = 0.299f * red[row][col] + 0.587f * green[row][col] + 0.114f * blue[row][col];
					rawData[row][col] *= lumi < K ? expos : lut[lumi];
				}
		else
			#pragma omp parallel for
			for(int row=0;row<height;row++)
				for(int col=0;col<width;col++){
					float lumi = 0.299f * rawData[row][col*3] + 0.587f * rawData[row][col*3+1] + 0.114f * rawData[row][col*3+2];
					float fac = lumi < K ? expos : lut[lumi];
					for (int c=0;c<3;c++)
						rawData[row][col*3+c] *= fac;
				}
		
    }
	if (settings->verbose) {
		t2e.set();
		printf("Exposure before %d usec\n", t2e.etime(t1e));
	}
	
}

} //namespace
