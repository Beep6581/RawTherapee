////////////////////////////////////////////////////////////////
//
//			Green Equilibration via directional average
//
//	copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: February 12, 2011
//
//	green_equil_RT.cc is free software: you can redistribute it and/or modify
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
#define TS 256	 // Tile size

#include <cmath>
#include <cstdlib>
#include <ctime>


#include "rt_math.h"
#include "rawimagesource.h"

namespace rtengine {

//void green_equilibrate()//for dcraw implementation
void RawImageSource::green_equilibrate(float thresh)
{  
	// thresh = threshold for performing green equilibration; max percentage difference of G1 vs G2  
	// G1-G2 differences larger than this will be assumed to be Nyquist texture, and left untouched
	
	int height=H, width=W; 

	// local variables
	float** rawptr = rawData;
	array2D<float> cfa (width,height,rawptr);
	//array2D<int> checker (width,height,ARRAY2D_CLEAR_DATA);
	
	
	//int verbose=1;
	
	static const float eps=1.0;	//tolerance to avoid dividing by zero
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// Fill G interpolated values with border interpolation and input values

	//int vote1, vote2;
	//int counter, vtest;
	
	//The green equilibration algorithm starts here
    /*
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int rr=1; rr < height-1; rr++)
		for (int cc=3-(FC(rr,2)&1); cc < width-2; cc+=2) {
			
			float pcorr = (cfa[rr+1][cc+1]-cfa[rr][cc])*(cfa[rr-1][cc-1]-cfa[rr][cc]);
			float mcorr = (cfa[rr-1][cc+1]-cfa[rr][cc])*(cfa[rr+1][cc-1]-cfa[rr][cc]);
			
			if (pcorr>0 && mcorr>0) {checker[rr][cc]=1;} else {checker[rr][cc]=0;}
			
			checker[rr][cc]=1;//test what happens if we always interpolate
		}
	
	counter=vtest=0;
	*/
	//now smooth the cfa data
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int rr=4; rr < height-4; rr++)
		for (int cc=5-(FC(rr,2)&1); cc < width-6; cc+=2) {
			//if (checker[rr][cc]) {
				//%%%%%%%%%%%%%%%%%%%%%%
				//neighbor checking code from Manuel Llorens Garcia
				float o1_1=cfa[(rr-1)][cc-1];
				float o1_2=cfa[(rr-1)][cc+1];
				float o1_3=cfa[(rr+1)][cc-1];
				float o1_4=cfa[(rr+1)][cc+1];
				float o2_1=cfa[(rr-2)][cc];
				float o2_2=cfa[(rr+2)][cc];
				float o2_3=cfa[(rr)][cc-2];
				float o2_4=cfa[(rr)][cc+2];
				
				float d1=(o1_1+o1_2+o1_3+o1_4)*0.25f;
				float d2=(o2_1+o2_2+o2_3+o2_4)*0.25f;
				
				float c1=(fabs(o1_1-o1_2)+fabs(o1_1-o1_3)+fabs(o1_1-o1_4)+fabs(o1_2-o1_3)+fabs(o1_3-o1_4)+fabs(o1_2-o1_4))/6.0;
				float c2=(fabs(o2_1-o2_2)+fabs(o2_1-o2_3)+fabs(o2_1-o2_4)+fabs(o2_2-o2_3)+fabs(o2_3-o2_4)+fabs(o2_2-o2_4))/6.0;
				//%%%%%%%%%%%%%%%%%%%%%%
				
				//vote1=(checker[rr-2][cc]+checker[rr][cc-2]+checker[rr][cc+2]+checker[rr+2][cc]);
				//vote2=(checker[rr+1][cc-1]+checker[rr+1][cc+1]+checker[rr-1][cc-1]+checker[rr-1][cc+1]);
				//if ((vote1==0 || vote2==0) && (c1+c2)<2*thresh*fabs(d1-d2)) vtest++;
				//if (vote1>0 && vote2>0 && (c1+c2)<4*thresh*fabs(d1-d2)) {
                if ((c1+c2)<4*thresh*fabs(d1-d2)) {
					//pixel interpolation
					float gin=cfa[rr][cc];
					
					float gse=(cfa[rr+1][cc+1])+0.5*(cfa[rr][cc]-cfa[rr+2][cc+2]);
					float gnw=(cfa[rr-1][cc-1])+0.5*(cfa[rr][cc]-cfa[rr-2][cc-2]);
					float gne=(cfa[rr-1][cc+1])+0.5*(cfa[rr][cc]-cfa[rr-2][cc+2]);
					float gsw=(cfa[rr+1][cc-1])+0.5*(cfa[rr][cc]-cfa[rr+2][cc-2]);
					
					
					
					float wtse=1.0f/(eps+SQR(cfa[rr+2][cc+2]-cfa[rr][cc])+SQR(cfa[rr+3][cc+3]-cfa[rr+1][cc+1]));
					float wtnw=1.0f/(eps+SQR(cfa[rr-2][cc-2]-cfa[rr][cc])+SQR(cfa[rr-3][cc-3]-cfa[rr-1][cc-1]));
					float wtne=1.0f/(eps+SQR(cfa[rr-2][cc+2]-cfa[rr][cc])+SQR(cfa[rr-3][cc+3]-cfa[rr-1][cc+1]));
					float wtsw=1.0f/(eps+SQR(cfa[rr+2][cc-2]-cfa[rr][cc])+SQR(cfa[rr+3][cc-3]-cfa[rr+1][cc-1]));
					
					float ginterp=(gse*wtse+gnw*wtnw+gne*wtne+gsw*wtsw)/(wtse+wtnw+wtne+wtsw);
					
					if ( ((ginterp-gin) < thresh*(ginterp+gin)) ) {
						rawData[rr][cc]=0.5f*(ginterp+gin);
						//counter++;
					}
					
				}
		//	}
		}
	//printf("pixfix count= %d; vtest= %d \n",counter,vtest);
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	// done
	/*t2 = clock();
	 dt = ((double)(t2-t1)) / CLOCKS_PER_SEC;
	 if (verbose) {
	 fprintf(stderr,_("elapsed time = %5.3fs\n"),dt);
	 }*/
	
	
}
}

#undef TS
