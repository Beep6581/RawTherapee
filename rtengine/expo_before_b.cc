
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
//
// This function uses parameters:
//      exposure (linear): 2^(-8..0..8): currently 0.5 +3
//      preserve (log)   : 0..8 : currently 0.1 1

void RawImageSource::processRawWhitepoint(float expos, float preser) {
	MyTime t1e,t2e;
	t1e.set();
	
	int width=W, height=H;
	
    // exposure correction inspired from G.Luijk
    if (fabs(preser)<0.001) {	
        // No highlight protection - simple mutiplication
		for (int c=0; c<4; c++) {
			chmax[c] *= expos;
		}
#pragma omp parallel for shared(expos)
        for (int row=0;row<height;row++)
            for (int col=0;col<width;col++) {
                rawData[row][col] *= expos;
			}
    } else {
        // save old image as it's overwritten by demosaic
        float** imgd = allocArray< float >(W,H); 
		
        // with memcpy it's faster than for (...)
        for (int i=0; i<H; i++) memcpy (imgd[i], rawData[i], W*sizeof(**imgd));
		
        // Demosaic to calc luminosity
        fast_demosaic (0,0,W,H);
		
		// calculate CIE luminosity
        float* luminosity = (float *) new float[width*height];
		
#pragma omp parallel default(shared)  
		{
            // CIE luminosity
#pragma omp for  
			for(int row=0;row<height;row++)
				for(int col=0;col<width;col++)
                    luminosity[row*width+col] = 
                    0.299f * (float)red[row][col] + 0.587f * (float)green[row][col] + 0.114f * (float)blue[row][col]; 
		}	
		
        // restore image destroyed by demosaic
        for (int i=0; i<H; i++) memcpy (rawData[i], imgd[i], W*sizeof(**imgd));
        freeArray<float>(imgd, H);
		
        // Find maximum to adjust LUTs. New float engines clips only at the very end
        int maxVal=0;
		for(int row=0;row<height;row++)
			for(int col=0;col<width;col++)
                if (rawData[row][col]>maxVal) maxVal = rawData[row][col];
		
		// Exposure correction with highlight preservation
        LUTf lut(maxVal+1);
		if(expos>1){
            float K = (float) maxVal / expos*exp(-preser*log(2.0));
            for (int j=0;j<=maxVal;j++) 
                lut[(int)j]=(((float)maxVal-K*expos)/((float)maxVal-K)*(j-maxVal)+(float) maxVal) / j;
			
			for (int c=0; c<4; c++) {
				chmax[c] *= expos;
			}
			
#pragma omp parallel for  shared(expos)
			for(int row=0;row<height;row++)
				for(int col=0;col<width;col++){
                    if (luminosity[row*width + col] < K) {
                        rawData[row][col] *= expos;
                    } else {
                        rawData[row][col] *= lut[luminosity[row*width+col]];
					}
                }
        } else {
            // Negative exposure
			float EV=log(expos)/log(2.0);                              // Convert exp. linear to EV
            float K = (float)maxVal * exp(-preser * log(2.0));
			
            for (int j=0;j<=maxVal;j++) 
                lut[(int)j] = exp(EV*((float)maxVal-j) / ((float)maxVal-K) * log(2.0));
			
#pragma omp parallel for  shared(expos)	  
			for(int row=0;row<height;row++)
				for(int col=0;col<width;col++){
                    if (luminosity[row*width+col]<K) {
                        rawData[row][col] *= expos;
                    } else {
                        rawData[row][col] *= lut[luminosity[row*width+col]];
					}	
				}
			for (int c=0; c<4; c++) {
				chmax[c] *= expos;
			}
        }	
		
		
        delete[] luminosity;
    }
	t2e.set();
	if( settings->verbose )
		printf("Exposure before  %d usec\n", t2e.etime(t1e));
	
}