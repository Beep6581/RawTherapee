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
#include "sleef.c"
#include "mytime.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "rt_math.h"
#define PIX_SORT(a,b) { if ((a)>(b)) {temp=(a);(a)=(b);(b)=temp;} }


#define med3(a0,a1,a2,a3,a4,a5,a6,a7,a8,median) { \
pp[0]=a0; pp[1]=a1; pp[2]=a2; pp[3]=a3; pp[4]=a4; pp[5]=a5; pp[6]=a6; pp[7]=a7; pp[8]=a8; \
PIX_SORT(pp[1],pp[2]); PIX_SORT(pp[4],pp[5]); PIX_SORT(pp[7],pp[8]); \
PIX_SORT(pp[0],pp[1]); PIX_SORT(pp[3],pp[4]); PIX_SORT(pp[6],pp[7]); \
PIX_SORT(pp[1],pp[2]); PIX_SORT(pp[4],pp[5]); PIX_SORT(pp[7],pp[8]); \
PIX_SORT(pp[0],pp[3]); PIX_SORT(pp[5],pp[8]); PIX_SORT(pp[4],pp[7]); \
PIX_SORT(pp[3],pp[6]); PIX_SORT(pp[1],pp[4]); PIX_SORT(pp[2],pp[5]); \
PIX_SORT(pp[4],pp[7]); PIX_SORT(pp[4],pp[2]); PIX_SORT(pp[6],pp[4]); \
PIX_SORT(pp[4],pp[2]); median=pp[4];} //pp4 = median

using namespace std;

namespace rtengine {
extern const Settings* settings;

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
		AlignedBufferMP<double> buffer(max(src->W,src->H));

		gaussHorizontal<float> (src->a, tmp1->a, buffer, src->W, src->H, radius);
		gaussHorizontal<float> (src->b, tmp1->b, buffer, src->W, src->H, radius);
		gaussVertical<float>   (tmp1->a, tmp1->a, buffer, src->W, src->H, radius);
		gaussVertical<float>   (tmp1->b, tmp1->b, buffer, src->W, src->H, radius);

//		gaussHorizontal<float> (src->L, tmp1->L, buffer, src->W, src->H, radius);
//		gaussVertical<float>   (tmp1->L, tmp1->L, buffer, src->W, src->H, radius);
	}

float chromave=0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:chromave)
#endif
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			float chroma = SQR(src->a[i][j]-tmp1->a[i][j])+SQR(src->b[i][j]-tmp1->b[i][j]);
			chromave += chroma;
			fringe[i*width+j]=chroma;
		}
	}
	chromave /= (height*width);
	float threshfactor = (thresh*chromave)/33.f;      // Calculated once to eliminate mult inside the next loop
//	printf("Chro %f \n",chromave);

// Issue 1674:
// often, CA isn't evenly distributed, e.g. a lot in contrasty regions and none in the sky.
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,16)
#endif
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			tmp1->a[i][j] = src->a[i][j];
			tmp1->b[i][j] = src->b[i][j];
			//test for pixel darker than some fraction of neighborhood ave, near an edge, more saturated than average
			/*if (100*tmp1->L[i][j]>50*src->L[i][j] && \*/
				/*1000*abs(tmp1->L[i][j]-src->L[i][j])>thresh*(tmp1->L[i][j]+src->L[i][j]) && \*/
			if (fringe[i*width+j]>threshfactor) {
				float atot=0.f;
				float btot=0.f;
				float norm=0.f;
				float wt;
				for (int i1=max(0,i-halfwin+1); i1<min(height,i+halfwin); i1++)
					for (int j1=max(0,j-halfwin+1); j1<min(width,j+halfwin); j1++) {
						//neighborhood average of pixels weighted by chrominance
						wt = 1.f/(fringe[i1*width+j1]+chromave);
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
void ImProcFunctions::PF_correct_RTcam(CieImage * src, CieImage * dst, double radius, int thresh) {
	int halfwin = ceil(2*radius)+1;

#include "rt_math.h"

	// local variables
	int width=src->W, height=src->H;
	float piid=3.14159265f/180.f;
	static float eps2=0.01f;
	//temporary array to store chromaticity
	int (*fringe);
	fringe = (int (*)) calloc ((height)*(width), sizeof *fringe);

	float** sraa;
		sraa = new float*[height];
		for (int i=0; i<height; i++)
			sraa[i] = new float[width];
/*
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++) {
				sraa[i][j]=src->C_p[i][j]*cos(piid*src->h_p[i][j]);
			}
			*/
	float** tmaa;
		tmaa = new float*[height];
		for (int i=0; i<height; i++)
			tmaa[i] = new float[width];

	float** srbb;
	srbb = new float*[height];
		for (int i=0; i<height; i++)
			srbb[i] = new float[width];
					

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++) {
				float2 sincosval = xsincosf(piid*src->h_p[i][j]);
				sraa[i][j]=src->C_p[i][j]*sincosval.y;			
				srbb[i][j]=src->C_p[i][j]*sincosval.x;
			}
	float** tmbb;
		tmbb = new float*[height];
	for (int i=0; i<height; i++)
			tmbb[i] = new float[width];

	/*float** tmL;
		tmL = new float*[height];
	for (int i=0; i<height; i++)
			tmL[i] = new float[width];
*/

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		AlignedBufferMP<double> buffer(max(src->W,src->H));
		gaussHorizontal<float> (sraa, tmaa, buffer, src->W, src->H, radius);
		gaussHorizontal<float> (srbb, tmbb, buffer, src->W, src->H, radius);
		gaussVertical<float>   (tmaa, tmaa, buffer, src->W, src->H, radius);
		gaussVertical<float>   (tmbb, tmbb, buffer, src->W, src->H, radius);
	//	gaussHorizontal<float> (src->sh_p, tmL, buffer, src->W, src->H, radius);
	//	gaussVertical<float>   (tmL, tmL, buffer, src->W, src->H, radius);

	}

float chromave=0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:chromave)
#endif
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			float chroma =SQR(sraa[i][j]-tmaa[i][j])+SQR(srbb[i][j]-tmbb[i][j]);
			chromave += chroma;
			fringe[i*width+j]=chroma;
		}
	}
	chromave /= (height*width);
	float threshfactor = (thresh*chromave)/33.f;      // Calculated once to eliminate mult inside the next loop
	//	printf("Chromave CAM %f \n",chromave);

// Issue 1674:
// often, CA isn't evenly distributed, e.g. a lot in contrasty regions and none in the sky.
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,16)
#endif
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			tmaa[i][j] = sraa[i][j];
			tmbb[i][j] = srbb[i][j];

			//test for pixel darker than some fraction of neighborhood ave, near an edge, more saturated than average
			/*if (100*tmp1->L[i][j]>50*src->L[i][j] && \*/
				/*1000*abs(tmp1->L[i][j]-src->L[i][j])>thresh*(tmp1->L[i][j]+src->L[i][j]) && \*/
			if (fringe[i*width+j]>threshfactor) {
				float atot=0.f;
				float btot=0.f;
				float norm=0.f;
				float wt;
				for (int i1=max(0,i-halfwin+1); i1<min(height,i+halfwin); i1++)
					for (int j1=max(0,j-halfwin+1); j1<min(width,j+halfwin); j1++) {
						//neighborhood average of pixels weighted by chrominance
						wt = 1.f/(fringe[i1*width+j1]+chromave+eps2);
						atot += wt*sraa[i1][j1];
						btot += wt*srbb[i1][j1];
						norm += wt;
					}
				if(norm > 0.f){	
				tmaa[i][j] = (atot/norm);
				tmbb[i][j] = (btot/norm);
				}
			}//end of ab channel averaging
		}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			dst->sh_p[i][j] = src->sh_p[i][j];
			float intera = tmaa[i][j];
			float interb = tmbb[i][j];
			dst->h_p[i][j]=(xatan2f(interb,intera))/piid;
			dst->C_p[i][j]=sqrt(SQR(interb)+SQR(intera));
		}
	}
                for (int i=0; i<height; i++)
                    delete [] sraa[i];
                delete [] sraa;
                for (int i=0; i<height; i++)
                    delete [] srbb[i];
                delete [] srbb;
                for (int i=0; i<height; i++)
                    delete [] tmaa[i];
                delete [] tmaa;
                for (int i=0; i<height; i++)
                    delete [] tmbb[i];
                delete [] tmbb;
           /*     for (int i=0; i<height; i++)
                    delete [] tmL[i];
                delete [] tmL;
			*/
	free(fringe);
}
void ImProcFunctions::Badpixelscam(CieImage * src, CieImage * dst, double radius, int thresh, int mode) {
	#include "rt_math.h"

	int halfwin = ceil(2*radius)+1;
    MyTime t1,t2;
    t1.set();
	//bool algogauss = settings->ciebadpixgauss;
	int width=src->W, height=src->H;
	float piid=3.14159265f/180.f;
	float shfabs, shmed;
	int i1, j1, tot;
	static float eps = 1.0f;
	static float eps2 =0.01f;
	float shsum, dirsh, norm, sum;	

	float** sraa;
		sraa = new float*[height];
		for (int i=0; i<height; i++)
			sraa[i] = new float[width];

	float** tmaa;
		tmaa = new float*[height];
		for (int i=0; i<height; i++)
			tmaa[i] = new float[width];

	float** srbb;
	srbb = new float*[height];
		for (int i=0; i<height; i++)
			srbb[i] = new float[width];

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++) {
				float2 sincosval = xsincosf(piid*src->h_p[i][j]);
				sraa[i][j]=src->C_p[i][j]*sincosval.y;			
				srbb[i][j]=src->C_p[i][j]*sincosval.x;
			}		
	float** tmbb;
		tmbb = new float*[height];
	for (int i=0; i<height; i++)
			tmbb[i] = new float[width];
    float ** badpix = new float *[height];
	float** tmL;
		tmL = new float*[height];
	for (int i=0; i<height; i++) {
			tmL[i] = new float[width];
			badpix[i] = new float [width];

		}
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		AlignedBufferMP<double> buffer(max(src->W,src->H));
		//chroma a and b
		if(mode==2) {//choice of gaussian blur 
		gaussHorizontal<float> (sraa, tmaa, buffer, src->W, src->H, radius);
		gaussHorizontal<float> (srbb, tmbb, buffer, src->W, src->H, radius);
		gaussVertical<float>   (tmaa, tmaa, buffer, src->W, src->H, radius);
		gaussVertical<float>   (tmbb, tmbb, buffer, src->W, src->H, radius);
		}
		//luma sh_p
		gaussHorizontal<float> (src->sh_p, tmL, buffer, src->W, src->H, 2.0);//low value to avoid artifacts
		gaussVertical<float>   (tmL, tmL, buffer, src->W, src->H, 2.0);

	}
if(mode==1){	//choice of median
#pragma omp parallel
	{
#pragma omp for	
	for (int i=0; i<height; i++) {
		int ip,in,jp,jn;
		float pp[9],temp;
		if (i<2) {ip=i+2;} else {ip=i-2;}
		if (i>height-3) {in=i-2;} else {in=i+2;}
		for (int j=0; j<width; j++) {
			if (j<2) {jp=j+2;} else {jp=j-2;}
			if (j>width-3) {jn=j-2;} else {jn=j+2;}
			med3(sraa[ip][jp],sraa[ip][j],sraa[ip][jn],sraa[i][jp],sraa[i][j],sraa[i][jn],sraa[in][jp],sraa[in][j],sraa[in][jn],tmaa[i][j]);
		}
	}
#pragma omp for		
	for (int i=0; i<height; i++) {
		int ip,in,jp,jn;
		float pp[9],temp;
		if (i<2) {ip=i+2;} else {ip=i-2;}
		if (i>height-3) {in=i-2;} else {in=i+2;}
		for (int j=0; j<width; j++) {
			if (j<2) {jp=j+2;} else {jp=j-2;}
			if (j>width-3) {jn=j-2;} else {jn=j+2;}
			med3(srbb[ip][jp],srbb[ip][j],srbb[ip][jn],srbb[i][jp],srbb[i][j],srbb[i][jn],srbb[in][jp],srbb[in][j],srbb[in][jn],tmbb[i][j]);
		}
	}
	}
}	
	
//luma badpixels	
float sh_thr = 4.5f;//low value for luma sh_p to avoid artifacts
float shthr = sh_thr / 24.0f;	      

#ifdef _OPENMP
  #pragma omp parallel for private(shfabs, shmed,i1,j1)
#endif
	for (int i=0; i < height; i++)
		for (int j=0; j < width; j++) {
			shfabs = fabs(src->sh_p[i][j]-tmL[i][j]);
			for (i1=max(0,i-2), shmed=0; i1<=min(i+2,height-1); i1++ )
				for (j1=max(0,j-2); j1<=min(j+2,width-1); j1++ ) {
					shmed += fabs(src->sh_p[i1][j1]-tmL[i1][j1]);
				}
			badpix[i][j] = (shfabs>((shmed-shfabs)*shthr));
		}
#ifdef _OPENMP
  #pragma omp parallel for private(shsum,norm,dirsh,sum,i1,j1) schedule(dynamic,16)
#endif
	for (int i=0; i < height; i++)
		for (int j=0; j < width; j++) {
			if (!badpix[i][j]) continue;
			norm=0.0f;
			shsum=0.0f;
			sum=0.0f;
			tot=0;
			for (i1=max(0,i-2), shmed=0; i1<=min(i+2,height-1); i1++ )
				for (j1=max(0,j-2); j1<=min(j+2,width-1); j1++ ) {
					if (i1==i && j1==j) continue;
					if (badpix[i1][j1]) continue;
					sum += src->sh_p[i1][j1]; tot++;
					dirsh = 1.f/(SQR(src->sh_p[i1][j1]-src->sh_p[i][j])+eps);
					shsum += dirsh*src->sh_p[i1][j1];
					norm += dirsh;
			}
			if (norm > 0.f) {
				src->sh_p[i][j]=shsum/norm;
			}
			else {
				if(tot > 0) src->sh_p[i][j]=sum / tot;
				}
		}
// end luma badpixels		

// begin chroma badpixels
float chrommed=0.f;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:chrommed)
#endif
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			float chroma =SQR(sraa[i][j]-tmaa[i][j])+SQR(srbb[i][j]-tmbb[i][j]);
			chrommed += chroma;
			badpix[i][j]=chroma;
		}
	}
	chrommed /= (height*width);
	float threshfactor = (thresh*chrommed)/33.f;   

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,16)
#endif
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			tmaa[i][j] = sraa[i][j];
			tmbb[i][j] = srbb[i][j];

			if (badpix[i][j]>threshfactor) {
				float atot=0.f;
				float btot=0.f;
				float norm=0.f;
				float wt;
				for (int i1=max(0,i-halfwin+1); i1<min(height,i+halfwin); i1++)
					for (int j1=max(0,j-halfwin+1); j1<min(width,j+halfwin); j1++) {
						wt = 1.f/(badpix[i1][j1]+chrommed+eps2);
						atot += wt*sraa[i1][j1];
						btot += wt*srbb[i1][j1];
						norm += wt;
					}
				if(norm > 0.f){	
				tmaa[i][j] = (atot/norm);
				tmbb[i][j] = (btot/norm);
				}
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i = 0; i < height; i++ ) {
		for(int j = 0; j < width; j++) {
			dst->sh_p[i][j] = src->sh_p[i][j];
			float intera = tmaa[i][j];
			float interb = tmbb[i][j];
			dst->h_p[i][j]=(xatan2f(interb,intera))/piid;
			dst->C_p[i][j]=sqrt(SQR(interb)+SQR(intera));
		}
	}
                for (int i=0; i<height; i++)
                    delete [] sraa[i];
                delete [] sraa;
                for (int i=0; i<height; i++)
                    delete [] srbb[i];
                delete [] srbb;
                for (int i=0; i<height; i++)
                    delete [] tmaa[i];
                delete [] tmaa;
                for (int i=0; i<height; i++)
                    delete [] tmbb[i];
                delete [] tmbb;
                for (int i=0; i<height; i++){
                    delete [] tmL[i];
					delete [] badpix[i];
				}
                delete [] tmL;
                delete [] badpix;
			
    t2.set();
    if( settings->verbose )
           printf("Ciecam badpixels:- %d usec\n", t2.etime(t1));
	
	
}
}
#undef PIX_SORT
#undef med3

