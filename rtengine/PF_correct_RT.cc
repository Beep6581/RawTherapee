////////////////////////////////////////////////////////////////
//
//		Chromatic Aberration Auto-correction
//
//		copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: November 24, 2010
//
//	PF_correct_RT.cc is free software: you can redistribute it and/or modify
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
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "gauss.h"
#include "improcfun.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "rt_math.h"

using namespace std;

namespace rtengine {

void ImProcFunctions::PF_correct_RT(LabImage * src, LabImage * dst, double radius, int thresh) { 
	int halfwin = ceil(2*radius)+1;
	
#include "rt_math.h"
		
	// local variables
	int width=src->W, height=src->H;
	//temporary array to store chromaticity
	int (*fringe);
	fringe = (int (*)) calloc ((height)*(width), sizeof *fringe);
	
	LabImage * tmp1; 
	tmp1 = new LabImage(width, height);
	
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		AlignedBuffer<double>* buffer = new AlignedBuffer<double> (max(src->W,src->H));
		gaussHorizontal<float> (src->a, tmp1->a, buffer, src->W, src->H, radius, multiThread);
		gaussHorizontal<float> (src->b, tmp1->b, buffer, src->W, src->H, radius, multiThread);
		gaussVertical<float>   (tmp1->a, tmp1->a, buffer, src->W, src->H, radius, multiThread);
		gaussVertical<float>   (tmp1->b, tmp1->b, buffer, src->W, src->H, radius, multiThread);

		//gaussHorizontal<float> (src->L, tmp1->L, buffer, src->W, src->H, radius, multiThread);
		//gaussVertical<float>   (tmp1->L, tmp1->L, buffer, src->W, src->H, radius, multiThread);

		delete buffer;
	}
	
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
	float chromave=0;
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			float chroma = SQR(src->a[i][j]-tmp1->a[i][j])+SQR(src->b[i][j]-tmp1->b[i][j]);
			chromave += chroma;
			fringe[i*width+j]=chroma;
		}
	}
	chromave /= (height*width);
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
	
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			tmp1->a[i][j] = src->a[i][j];
			tmp1->b[i][j] = src->b[i][j];
			//test for pixel darker than some fraction of neighborhood ave, near an edge, more saturated than average
			/*if (100*tmp1->L[i][j]>50*src->L[i][j] && \*/
				/*1000*abs(tmp1->L[i][j]-src->L[i][j])>thresh*(tmp1->L[i][j]+src->L[i][j]) && \*/
			if (33*fringe[i*width+j]>thresh*chromave) {
				float atot=0;
				float btot=0;
				float norm=0;
				float wt;
				for (int i1=max(0,i-halfwin+1); i1<min(height,i+halfwin); i1++)
					for (int j1=max(0,j-halfwin+1); j1<min(width,j+halfwin); j1++) {
						//neighborhood average of pixels weighted by chrominance
						wt = 1/(fringe[i1*width+j1]+chromave);
						atot += wt*src->a[i1][j1];
						btot += wt*src->b[i1][j1];
						norm += wt;
					}
				tmp1->a[i][j] = (int)(atot/norm);
				tmp1->b[i][j] = (int)(btot/norm);
			}//end of ab channel averaging
		}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			dst->L[i][j] = src->L[i][j];
			dst->a[i][j] = tmp1->a[i][j];
			dst->b[i][j] = tmp1->b[i][j];
		}
	}
	
	delete tmp1;
	free(fringe);
}

}

