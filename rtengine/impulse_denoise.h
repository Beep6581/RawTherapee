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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Gabor's implementation of bilateral filtering, without input pixel

#define NBRWT(a,b) (src[i - a][j - b] * ec[src[i - a][j - b]-src[i][j]+0x10000])
#define NORM(a,b) (1 + ec[src[i - a][j - b]-src[i][j]+0x10000])

#define RB_BEGIN(a,b)   double scale = (a); \
	int* ec = new int [0x20000]; \
	for (int i=0; i<0x20000; i++) \
	ec[i] = (int)(exp(-(double)(i-0x10000)*(double)(i-0x10000) / (2.0*rangewidth*rangewidth))*scale); \
	int rstart = b; \
	int rend = H-b; \
	int cstart = b; \
	int cend = W-b;

#define RB_END(b)       buffer[i][j] = v; }} delete [] ec; \
	for (int i=0; i<H; i++)  \
	for (int j=0; j<W; j++)  \
	if (i<rstart || j<cstart || i>=rend || j>=cend) \
	dst[i][j] = src[i][j]; \
	else \
	dst[i][j] = buffer[i][j];

#define RB_OPER5 for (int i=rstart; i<rend; i++) { \
for (int j=cstart; j<cend; j++) { \
	A v = NBRWT(-2,-2) + NBRWT(-2,-1) + NBRWT(-2,0) + NBRWT(-2,1) + NBRWT(-2,2) + \
		NBRWT(-1,-2) + NBRWT(-1,-1) + NBRWT(-1,0) + NBRWT(-1,1) + NBRWT(-1,2) + \
		NBRWT(0,-2) + NBRWT(0,-1) /*+ NBRWT(0,0)*/ + NBRWT(0,1) + NBRWT(0,2) + \
		NBRWT(1,-2) + NBRWT(1,-1) + NBRWT(1,0) + NBRWT(1,1) + NBRWT(1,2) + \
		NBRWT(2,-2) + NBRWT(2,-1) + NBRWT(2,0) + NBRWT(2,1) + NBRWT(2,2); \
	v /= NORM(-2,-2) + NORM(-2,-1) + NORM(-2,0) + NORM(-2,1) + NORM(-2,2) + \
		NORM(-1,-2) + NORM(-1,-1) + NORM(-1,0) + NORM(-1,1) + NORM(-1,2) + \
		NORM(0,-2) + NORM(0,-1) /*+ NORM(0,0)*/ + NORM(0,1) + NORM(0,2) + \
		NORM(1,-2) + NORM(1,-1) + NORM(1,0) + NORM(1,1) + NORM(1,2) + \
		NORM(2,-2) + NORM(2,-1) + NORM(2,0) + NORM(2,1) + NORM(2,2);


template<class T, class A> void rangeblur (T** src, T** dst, T** buffer, int W, int H, double rangewidth, bool multiThread) {
	
	RB_BEGIN(753,2)
#pragma omp parallel for if (multiThread)
    RB_OPER5
    RB_END(2)
	
}



template<class T> void impulse_nr (T** src, T** dst, int width, int height, double thresh) {
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// impulse noise removal
	// local variables
	
	float hpfabs, hfnbrave;
	
	// buffer for the lowpass image
    unsigned short ** lpf = new unsigned short *[height];
    for (int i=0; i<height; i++) {
        lpf[i] = new unsigned short [width];
        //memset (lpf[i], 0, width*sizeof(float));
    }
	
	// buffer for the highpass image
    unsigned short ** impish = new unsigned short *[height];
	 for (int i=0; i<height; i++) {
	 impish[i] = new unsigned short [width];
	 //memset (impish[i], 0, width*sizeof(unsigned short));
	 }
	
	//The cleaning algorithm starts here
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// modified bilateral filter for lowpass image, omitting input pixel
	
	static float eps = 1.0;
	float wtdsum, dirwt, norm;
	int i1, j1;	
	
	//rangeblur<unsigned short, unsigned int> (src, lpf, impish /*used as buffer here*/, width, height, thresh, false);
	
	AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(width,height));

	gaussHorizontal<unsigned short> (src, lpf, buffer, width, height, 2.0, false /*multiThread*/);
	gaussVertical<unsigned short>   (lpf, lpf, buffer, width, height, 2.0, false);

	delete buffer;

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	for (int i=0; i < height; i++)
		for (int j=0; j < width; j++) {
			
			hpfabs = fabs(src[i][j]-lpf[i][j]);
			//block average of high pass data
			for (i1=MAX(0,i-2), hfnbrave=0; i1<=MIN(i+2,height-1); i1++ )
				for (j1=MAX(0,j-2); j1<=MIN(j+2,width-1); j1++ ) {
					hfnbrave += fabs(src[i1][j1]-lpf[i1][j1]);
				}
			hfnbrave = (hfnbrave-hpfabs)/24;
			hpfabs>(hfnbrave*(5.5-thresh)) ? impish[i][j]=1 : impish[i][j]=0;
			
		}//now impulsive values have been corrected
	
	for (int i=0; i < height; i++)
		for (int j=0; j < width; j++) {
			if (!impish[i][j]) continue;
			norm=0.0;
			wtdsum=0.0;
			for (i1=MAX(0,i-2), hfnbrave=0; i1<=MIN(i+2,height-1); i1++ )
				for (j1=MAX(0,j-2); j1<=MIN(j+2,width-1); j1++ ) {
					if (i1==i && j1==j) continue;
					if (impish[i1][j1]) continue;
					dirwt = 1/(SQR(src[i1][j1]-src[i][j])+eps);//use more sophisticated rangefn???
					wtdsum += dirwt*src[i1][j1];
					norm += dirwt;
			}
			//wtdsum /= norm;
			if (norm) {
				src[i][j]=wtdsum/norm;//low pass filter
			} 
			
		}//now impulsive values have been corrected
	
    for (int i=0; i<height; i++) 
        delete [] lpf[i];
	delete [] lpf;
	
	for (int i=0; i<height; i++) 
        delete [] impish[i];
	delete [] impish;
	
}




