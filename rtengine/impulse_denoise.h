/*
 *  This file is part of RawTherapee.
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but widthITheightOUT ANY widthARRANTY; without even the implied warranty of
 *  MERCheightANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  2010 Emil Martinec <ejmartin@uchicago.edu>
 *    
 */

#define SQR(x) ((x)*(x))

#include <cstddef>
#include <algorithm>
#include <labimage.h>
#include <improcfun.h>


namespace rtengine {

void ImProcFunctions::impulse_nr (LabImage* lab, double thresh) {
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// impulse noise removal
	// local variables
	
	int width = lab->W;
	int height = lab->H;
	
	float hpfabs, hfnbrave;
	
	// buffer for the lowpass image
    float ** lpf = new float *[height];
	// buffer for the highpass image
    float ** impish = new float *[height];
    for (int i=0; i<height; i++) {
        lpf[i] = new float [width];
        //memset (lpf[i], 0, width*sizeof(float));
		impish[i] = new float [width];
		//memset (impish[i], 0, width*sizeof(unsigned short));
    }

	
	//The cleaning algorithm starts here
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// modified bilateral filter for lowpass image, omitting input pixel; or Gaussian blur
	
	static float eps = 1.0;
	float wtdsum[3], dirwt, norm;
	int i1, j1;	
	
	//rangeblur<unsigned short, unsigned int> (lab->L, lpf, impish /*used as buffer here*/, width, height, thresh, false);
	
	AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(width,height));

	gaussHorizontal<float> (lab->L, lpf, buffer, width, height, MAX(2.0,thresh-1.0), false /*multiThread*/);
	gaussVertical<float>   (lpf, lpf, buffer, width, height, MAX(2.0,thresh-1.0), false);

	delete buffer;

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	float impthr = MAX(1.0,5.5-thresh);
	
	for (int i=0; i < height; i++)
		for (int j=0; j < width; j++) {
			
			hpfabs = fabs(lab->L[i][j]-lpf[i][j]);
			//block average of high pass data
			for (i1=MAX(0,i-2), hfnbrave=0; i1<=MIN(i+2,height-1); i1++ )
				for (j1=MAX(0,j-2); j1<=MIN(j+2,width-1); j1++ ) {
					hfnbrave += fabs(lab->L[i1][j1]-lpf[i1][j1]);
				}
			hfnbrave = (hfnbrave-hpfabs)/24;
			hpfabs>(hfnbrave*impthr) ? impish[i][j]=1 : impish[i][j]=0;
			
		}//now impulsive values have been identified
	
	for (int i=0; i < height; i++)
		for (int j=0; j < width; j++) {
			if (!impish[i][j]) continue;
			norm=0.0;
			wtdsum[0]=wtdsum[1]=wtdsum[2]=0.0;
			for (i1=MAX(0,i-2), hfnbrave=0; i1<=MIN(i+2,height-1); i1++ )
				for (j1=MAX(0,j-2); j1<=MIN(j+2,width-1); j1++ ) {
					if (i1==i && j1==j) continue;
					if (impish[i1][j1]) continue;
					dirwt = 1/(SQR(lab->L[i1][j1]-lab->L[i][j])+eps);//use more sophisticated rangefn???
					wtdsum[0] += dirwt*lab->L[i1][j1];
					wtdsum[1] += dirwt*lab->a[i1][j1];
					wtdsum[2] += dirwt*lab->b[i1][j1];
					norm += dirwt;
			}
			//wtdsum /= norm;
			if (norm) {
				lab->L[i][j]=wtdsum[0]/norm;//low pass filter
				lab->a[i][j]=wtdsum[1]/norm;//low pass filter
				lab->b[i][j]=wtdsum[2]/norm;//low pass filter
			} 
			
		}//now impulsive values have been corrected
	
    for (int i=0; i<height; i++) {
        delete [] lpf[i];
		delete [] impish[i];
	}
	delete [] lpf;
	delete [] impish;
	
}

}



