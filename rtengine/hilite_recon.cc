////////////////////////////////////////////////////////////////
//
//		Highlight reconstruction
//
//		copyright (c) 2008-2011  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: June 16, 2011
//
//	hilite_recon.cc is free software: you can redistribute it and/or modify
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

//#include "rtengine.h"
#include <cstddef>
#include <cmath>
#include "curves.h"
#include "array2D.h"
#include "improcfun.h"
#include "rawimagesource.h"
//#include "stack1.h"


#ifdef _OPENMP
#include <omp.h>
#endif

#include "rt_math.h"
#include "rawimagesource.h"
#ifdef __SSE2__
#include "sleefsseavx.c"
#endif // __SSE2__




#define FOREACHCOLOR for (int c=0; c < ColorCount; c++)				

//#include "RGBdefringe.cc"

namespace rtengine {


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if defined( __SSE2__ ) && defined( WIN32 )
__attribute__((force_align_arg_pointer)) void RawImageSource::boxblur2(float** src, float** dst, int H, int W, int box ) 
#else
void RawImageSource::boxblur2(float** src, float** dst, int H, int W, int box ) 
#endif
{
	
	array2D<float> temp(W,H);
	
	//box blur image channel; box size = 2*box+1
	//horizontal blur
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int row = 0; row < H; row++) {
		int len = box + 1;
		temp[row][0] = src[row][0]/len;
		for (int j=1; j<=box; j++) {
			temp[row][0] += src[row][j]/len;
		}
		for (int col=1; col<=box; col++) {
			temp[row][col] = (temp[row][col-1]*len + src[row][col+box])/(len+1);
			len ++;
		}
		for (int col = box+1; col < W-box; col++) {
			temp[row][col] = temp[row][col-1] + (src[row][col+box] - src[row][col-box-1])/len;
		}
		for (int col=W-box; col<W; col++) {
			temp[row][col] = (temp[row][col-1]*len - src[row][col-box-1])/(len-1);
			len --;
		}
	}
	
#ifdef __SSE2__
	//vertical blur
#ifdef _OPENMP
#pragma omp parallel
#endif
{
	float len = box + 1;
	__m128 lenv = _mm_set1_ps( len );
	__m128 lenp1v = _mm_set1_ps( len + 1.0f );
	__m128 onev = _mm_set1_ps( 1.0f );
	__m128 tempv;
#ifdef _OPENMP
#pragma omp for
#endif
	for (int col = 0; col < W-3; col+=4) {
		tempv = LVFU(temp[0][col]) / lenv;
		for (int i=1; i<=box; i++) {
			tempv = tempv + LVFU(temp[i][col]) / lenv;
		}
		_mm_storeu_ps( &dst[0][col], tempv);
		for (int row=1; row<=box; row++) {
			_mm_storeu_ps( &dst[row][col], (LVFU(dst[(row-1)][col])*lenv + LVFU(temp[(row+box)][col]))/lenp1v);
			lenv = lenp1v;
			lenp1v = lenp1v + onev;
		}
		for (int row = box+1; row < H-box; row++) {
			_mm_storeu_ps( &dst[row][col], LVFU(dst[(row-1)][col]) + (LVFU(temp[(row+box)][col]) - LVFU(temp[(row-box-1)][col]))/lenv );
		}
		for (int row=H-box; row<H; row++) {
			lenp1v = lenv;
			lenv = lenv - onev;
			_mm_storeu_ps( &dst[row][col], (LVFU(dst[(row-1)][col])*lenp1v - LVFU(temp[(row-box-1)][col])) / lenv );
		}
	}
}
	for (int col = W-(W%4); col < W; col++) {
		int len = box + 1;
		dst[0][col] = temp[0][col]/len;
		for (int i=1; i<=box; i++) {
			dst[0][col] += temp[i][col]/len;
		}
		for (int row=1; row<=box; row++) {
			dst[row][col] = (dst[(row-1)][col]*len + temp[(row+box)][col])/(len+1);
			len ++;
		}
		for (int row = box+1; row < H-box; row++) {
			dst[row][col] = dst[(row-1)][col] + (temp[(row+box)][col] - temp[(row-box-1)][col])/len;
		}
		for (int row=H-box; row<H; row++) {
			dst[row][col] = (dst[(row-1)][col]*len - temp[(row-box-1)][col])/(len-1);
			len --;
		}
	}

#else
	//vertical blur
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int col = 0; col < W; col++) {
		int len = box + 1;
		dst[0][col] = temp[0][col]/len;
		for (int i=1; i<=box; i++) {
			dst[0][col] += temp[i][col]/len;
		}
		for (int row=1; row<=box; row++) {
			dst[row][col] = (dst[(row-1)][col]*len + temp[(row+box)][col])/(len+1);
			len ++;
		}
		for (int row = box+1; row < H-box; row++) {
			dst[row][col] = dst[(row-1)][col] + (temp[(row+box)][col] - temp[(row-box-1)][col])/len;
		}
		for (int row=H-box; row<H; row++) {
			dst[row][col] = (dst[(row-1)][col]*len - temp[(row-box-1)][col])/(len-1);
			len --;
		}
	}
#endif
	
}

void RawImageSource::boxblur_resamp(float **src, float **dst, float & max_f, int H, int W, int box, int samp )
{
	
	array2D<float> temp((W/samp)+ ((W%samp)==0 ? 0 : 1),H);
	
	float maxtmp=0.0f;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
	float maxtmpthr = 0;
	float tempval;
#ifdef _OPENMP
#pragma omp for
#endif
	//box blur image channel; box size = 2*box+1
	//horizontal blur
	for (int row = 0; row < H; row++) {
		int len = box + 1;
		tempval = src[row][0]/len;
		maxtmpthr = max(maxtmpthr,src[row][0]);
		for (int j=1; j<=box; j++) {
			tempval += src[row][j]/len;
			maxtmpthr = max(maxtmpthr,src[row][j]);
		}
		temp[row][0] = tempval;
		for (int col=1; col<=box; col++) {
			tempval = (tempval*len + src[row][col+box])/(len+1);
			if(col%samp == 0)
				temp[row][col/samp] = tempval;
			maxtmpthr = max(maxtmpthr,src[row][col]);
			len ++;
		}
		for (int col = box+1; col < W-box; col++) {
			tempval = tempval + (src[row][col+box] - src[row][col-box-1])/len;
			if(col%samp == 0)
				temp[row][col/samp] = tempval;
			maxtmpthr = max(maxtmpthr,src[row][col]);
		}
		for (int col=W-box; col<W; col++) {
			tempval = (tempval*len - src[row][col-box-1])/(len-1);
			if(col%samp == 0)
				temp[row][col/samp] = tempval;
			maxtmpthr = max(maxtmpthr,src[row][col]);
			len --;
		}
	}
#ifdef _OPENMP
#pragma omp critical
#endif
{
	maxtmp = max(maxtmp,maxtmpthr);
}
}
	
	max_f = maxtmp;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
	float tempval;
#ifdef _OPENMP
#pragma omp for
#endif
	//vertical blur
	for (int col = 0; col < W/samp; col++) {
		int len = box + 1;
		tempval = temp[0][col]/len;
		for (int i=1; i<=box; i++) {
			tempval += temp[i][col]/len;
		}
		dst[0][col] = tempval;
		for (int row=1; row<=box; row++) {
			tempval = (tempval*len + temp[(row+box)][col])/(len+1);
			if(row%samp == 0)
				dst[row/samp][col] = tempval;
			len ++;
		}
		for (int row = box+1; row < H-box; row++) {
			tempval = tempval + (temp[(row+box)][col] - temp[(row-box-1)][col])/len;
			if(row%samp == 0)
				dst[row/samp][col] = tempval;
		}
		for (int row=H-box; row<H; row++) {
			tempval = (tempval*len - temp[(row-box-1)][col])/(len-1);
			if(row%samp == 0)
				dst[row/samp][col] = tempval;
			len --;
		}
	}
}

	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
void RawImageSource :: HLRecovery_inpaint (float** red, float** green, float** blue)
{
	double progress = 0.0;

	if (plistener) {
		plistener->setProgressStr ("HL reconstruction...");
		plistener->setProgress (progress);
	}
	
	int height = H;
	int width = W;
	
	const int range = 2;
	const int pitch = 4;

	
	int hfh = (height-(height%pitch))/pitch;
	int hfw = (width-(width%pitch))/pitch;
	
	static const int numdirs = 4;				
	
	static const float threshpct = 0.25;
	static const float fixthreshpct = 0.7;
	static const float maxpct = 0.95;

	//%%%%%%%%%%%%%%%%%%%%
	//for blend algorithm:
	static const float blendthresh=1.0;
	const int ColorCount=3;
	// Transform matrixes rgb>lab and back
	static const float trans[2][ColorCount][ColorCount] =
	{ { { 1,1,1 }, { 1.7320508,-1.7320508,0 }, { -1,-1,2 } },
		{ { 1,1,1 }, { 1,-1,1 }, { 1,1,-1 } } };
	static const float itrans[2][ColorCount][ColorCount] =
	{ { { 1,0.8660254,-0.5 }, { 1,-0.8660254,-0.5 }, { 1,0,1 } },
		{ { 1,1,1 }, { 1,-1,1 }, { 1,1,-1 } } };
	//%%%%%%%%%%%%%%%%%%%%

	
	float max_f[3], thresh[3], fixthresh[3], norm[3];
		
	//float red1, green1, blue1;//diagnostic
	float chmaxalt[4]={0,0,0,0};//diagnostic
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//halfsize demosaic
	
	
	multi_array2D<float,3> hfsize (hfw+1,hfh+1,ARRAY2D_CLEAR_DATA);
	
	boxblur_resamp(red,hfsize[0],chmaxalt[0],height,width,range,pitch);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	boxblur_resamp(green,hfsize[1],chmaxalt[1],height,width,range,pitch);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	boxblur_resamp(blue,hfsize[2],chmaxalt[2],height,width,range,pitch);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	
	//blur image
	//for (int m=0; m<3; m++)
	//	boxblur2(hfsize[m],hfsizeblur[m],hfh,hfw,3);
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	for (int c=0; c<3; c++) {
		thresh[c] = chmax[c]*threshpct;
		fixthresh[c] = chmax[c]*fixthreshpct;
		max_f[c] = chmax[c]*maxpct;//min(chmax[0],chmax[1],chmax[2])*maxpct;
		norm[c] = 1.0/(max_f[c]-fixthresh[c]);
	}
	float whitept = max(max_f[0],max_f[1],max_f[2]);
	float clippt  = min(max_f[0],max_f[1],max_f[2]);
	float medpt   = max_f[0]+max_f[1]+max_f[2]-whitept-clippt;
	float blendpt = blendthresh*clippt;
	
	float camwb[4];
	for (int c=0; c<4; c++) camwb[c]=ri->get_cam_mul(c);

	multi_array2D<float,3> channelblur(width,height,ARRAY2D_CLEAR_DATA);
	multi_array2D<float,5> hilite_full(width,height,ARRAY2D_CLEAR_DATA);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
		
	// blur RGB channels
	boxblur2(red  ,channelblur[0],height,width,4);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	boxblur2(green,channelblur[1],height,width,4);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	boxblur2(blue ,channelblur[2],height,width,4);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	
	float hipass_sum=0, hipass_norm=0.00;
	
	// set up which pixels are clipped or near clipping
#ifdef _OPENMP
#pragma omp parallel for reduction(+:hipass_sum,hipass_norm)
#endif
	for (int i=0; i<height; i++) {
		for (int j=0; j<width; j++) {
			
			//if one or more channels is highlight but none are blown, add to highlight accumulator
			
			if ((red[i][j]>thresh[0] || green[i][j]>thresh[1] || blue[i][j]>thresh[2]) &&
				(red[i][j]<max_f[0] && green[i][j]<max_f[1] && blue[i][j]<max_f[2])) {
				
				hipass_sum += fabs(channelblur[0][i][j]-red[i][j]) + fabs(channelblur[1][i][j]-green[i][j]) + fabs(channelblur[2][i][j]-blue[i][j]);
				hipass_norm++;
				
				hilite_full[0][i][j] = red[i][j];
				hilite_full[1][i][j] = green[i][j];
				hilite_full[2][i][j] = blue[i][j];
				hilite_full[3][i][j] = 1;
				hilite_full[4][i][j] = 1;
				
			}
			//if (i%100==0 && j%100==0)
			//	printf("row=%d  col=%d  r=%f  g=%f  b=%f hilite=%f  \n",i,j,hilite_full[0][i][j],hilite_full[1][i][j],hilite_full[2][i][j],hilite_full[3][i][j]);
		}
	}//end of filling highlight array

	hipass_norm += 0.01;

	float hipass_ave = (hipass_sum/hipass_norm);
	
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//blur highlight data
	boxblur2(hilite_full[4],hilite_full[4],height,width,1);
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<height; i++) {
		for (int j=0; j<width; j++) {
			
			float hipass = fabs(channelblur[0][i][j]-red[i][j]) + fabs(channelblur[1][i][j]-green[i][j]) + fabs(channelblur[2][i][j]-blue[i][j]);

			if (hipass > 2*hipass_ave) {
				//too much variation
				hilite_full[0][i][j] = hilite_full[1][i][j] = hilite_full[2][i][j] = hilite_full[3][i][j] = 0;
				continue;
			}
						
			if (hilite_full[4][i][j]>0.00001 && hilite_full[4][i][j]<0.95) {
				//too near an edge, could risk using CA affected pixels, therefore omit
				hilite_full[0][i][j] = hilite_full[1][i][j] = hilite_full[2][i][j] = hilite_full[3][i][j] = 0;
			}
		}
	}
	
	for (int c=0; c<3; c++) channelblur[c](1,1);//free up some memory
	
	multi_array2D<float,4> hilite(hfw+1,hfh+1,ARRAY2D_CLEAR_DATA);
		
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// blur and resample highlight data; range=size of blur, pitch=sample spacing
	for (int m=0; m<4; m++) {
		boxblur_resamp(hilite_full[m],hilite[m],chmaxalt[m],height,width,range,pitch);
		if(plistener){
			progress += 0.05;
			plistener->setProgress(progress);
		}
	}
	for (int c=0; c<5; c++) hilite_full[c](1,1);//free up some memory
	
	multi_array2D<float,4*numdirs> hilite_dir(hfw,hfh,ARRAY2D_CLEAR_DATA);
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//blur highlights
	//for (int m=0; m<4; m++)
	//	boxblur2(hilite[m],hilite[m],hfh,hfw,4);
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	
	LUTf invfn(0x10000);
	
	for (int i=0; i<0x10000; i++)
		invfn[i] = 1.0/(1.0+i);
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	//fill gaps in highlight map by directional extension
	//raster scan from four corners
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j=1; j<hfw-1; j++) 
		for (int i=2; i<hfh-2; i++) {
			//from left 
			if (hilite[3][i][j]>0.01) {
				for (int c=0; c<4; c++) {
					hilite_dir[c][i][j] = hilite[c][i][j]/hilite[3][i][j];
				}
			} else {
				for (int c=0; c<4; c++) {
					hilite_dir[c][i][j] = 0.1*((hilite_dir[0+c][i-2][j-1]+hilite_dir[0+c][i-1][j-1]+hilite_dir[0+c][i][j-1]+hilite_dir[0+c][i+1][j-1]+hilite_dir[0+c][i+2][j-1])/
											  (hilite_dir[0+3][i-2][j-1]+hilite_dir[0+3][i-1][j-1]+hilite_dir[0+3][i][j-1]+hilite_dir[0+3][i+1][j-1]+hilite_dir[0+3][i+2][j-1]+0.00001));
					hilite_dir[4+c][i][j+1]  += hilite_dir[c][i][j];
					hilite_dir[8+c][i-2][j]  += hilite_dir[c][i][j];
					hilite_dir[12+c][i+2][j] += hilite_dir[c][i][j];
				}
			}
		}
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j=hfw-2; j>0; j--)  
		for (int i=2; i<hfh-2; i++) {
			//from right 
			if (hilite[3][i][j]>0.01) {
				for (int c=0; c<4; c++) {
					hilite_dir[4+c][i][j] = hilite[c][i][j]/hilite[3][i][j];
				}
			} else {
				for (int c=0; c<4; c++) {
					hilite_dir[4+c][i][j] = 0.1*((hilite_dir[4+c][(i-2)][(j+1)]+hilite_dir[4+c][(i-1)][(j+1)]+hilite_dir[4+c][(i)][(j+1)]+hilite_dir[4+c][(i+1)][(j+1)]+hilite_dir[4+c][(i+2)][(j+1)])/
														  (hilite_dir[4+3][(i-2)][(j+1)]+hilite_dir[4+3][(i-1)][(j+1)]+hilite_dir[4+3][(i)][(j+1)]+hilite_dir[4+3][(i+1)][(j+1)]+hilite_dir[4+3][(i+2)][(j+1)]+0.00001));
					hilite_dir[8+c][i-2][j] += hilite_dir[4+c][i][j];
					hilite_dir[12+c][i+2][j] += hilite_dir[4+c][i][j];
				}
			}
		}
		
	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=1; i<hfh-1; i++) 
		for (int j=2; j<hfw-2; j++) {
			//if (i%100==0 && j%100==0)
			//	printf("row=%d  col=%d  r=%f  g=%f  b=%f hilite=%f  \n",i,j,hilite[0][i][j],hilite[1][i][j],hilite[2][i][j],hilite[3][i][j]);
			
			//from top 
			if (hilite[3][i][j]>0.01) {
				for (int c=0; c<4; c++) {
					hilite_dir[8+c][i][j] = hilite[c][i][j]/hilite[3][i][j];
				}
			} else {
				for (int c=0; c<4; c++) {
					hilite_dir[8+c][i][j] = 0.1*((hilite_dir[8+c][i-1][j-2]+hilite_dir[8+c][i-1][j-1]+hilite_dir[8+c][i-1][j]+hilite_dir[8+c][i-1][j+1]+hilite_dir[8+c][i-1][j+2])/
											  (hilite_dir[8+3][i-1][j-2]+hilite_dir[8+3][i-1][j-1]+hilite_dir[8+3][i-1][j]+hilite_dir[8+3][i-1][j+1]+hilite_dir[8+3][i-1][j+2]+0.00001));
					hilite_dir[12+c][i+1][j] += hilite_dir[8+c][i][j];
				}
			}
		}

	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=hfh-2; i>0; i--) 
		for (int j=2; j<hfw-2; j++) {
			//from bottom 
			if (hilite[3][i][j]>0.01) {
				for (int c=0; c<4; c++) {
					hilite_dir[12+c][i][j] = hilite[c][i][j]/hilite[3][i][j];
				}
			} else {
				for (int c=0; c<4; c++) {
					hilite_dir[12+c][i][j] = 0.1*((hilite_dir[12+c][(i+1)][(j-2)]+hilite_dir[12+c][(i+1)][(j-1)]+hilite_dir[12+c][(i+1)][(j)]+hilite_dir[12+c][(i+1)][(j+1)]+hilite_dir[12+c][(i+1)][(j+2)])/
														  (hilite_dir[12+3][(i+1)][(j-2)]+hilite_dir[12+3][(i+1)][(j-1)]+hilite_dir[12+3][(i+1)][(j)]+hilite_dir[12+3][(i+1)][(j+1)]+hilite_dir[12+3][(i+1)][(j+2)]+0.00001));
				}
			}
	
		}

	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
	//fill in edges
	for (int dir=0; dir<numdirs; dir++) {
		for (int i=1; i<hfh-1; i++) 
			for (int c=0; c<4; c++) {
				hilite_dir[dir*4+c][i][0] = hilite_dir[dir*4+c][i][1];
				hilite_dir[dir*4+c][i][hfw-1] = hilite_dir[dir*4+c][i][hfw-2];
			}
		
		for (int j=1; j<hfw-1; j++) 
			for (int c=0; c<4; c++) {
				hilite_dir[dir*4+c][0][j] = hilite_dir[dir*4+c][1][j];
				hilite_dir[dir*4+c][hfh-1][j] = hilite_dir[dir*4+c][hfh-2][j];
			}
		
		for (int c=0; c<4; c++) {
			hilite_dir[dir*4+c][0][0] = hilite_dir[dir*4+c][1][0] = hilite_dir[dir*4+c][0][1] = hilite_dir[dir*4+c][1][1] = hilite_dir[dir*4+c][2][2];
			hilite_dir[dir*4+c][0][hfw-1] = hilite_dir[dir*4+c][1][hfw-1] = hilite_dir[dir*4+c][0][hfw-2] = hilite_dir[dir*4+c][1][hfw-2] = hilite_dir[dir*4+c][2][hfw-3];
			hilite_dir[dir*4+c][hfh-1][0] = hilite_dir[dir*4+c][hfh-2][0] = hilite_dir[dir*4+c][hfh-1][1] = hilite_dir[dir*4+c][hfh-2][1] = hilite_dir[dir*4+c][hfh-3][2];
			hilite_dir[dir*4+c][hfh-1][hfw-1] = hilite_dir[dir*4+c][hfh-2][hfw-1] = hilite_dir[dir*4+c][hfh-1][hfw-2] = hilite_dir[dir*4+c][hfh-2][hfw-2] = hilite_dir[dir*4+c][hfh-3][hfw-3];
		}
	}

	if(plistener){
		progress += 0.05;
		plistener->setProgress(progress);
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	/*for (int m=0; m<4*numdirs; m++) {
	 boxblur2(hilite_dir[m],hilite_dir[m],hfh,hfw,4);
	 }*/
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//now we have highlight data extended in various directions
	//next step is to build back local data by averaging
	

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// now reconstruct clipped channels using color ratios
	
	//const float Yclip = 0.3333*(max[0] + max[1] + max[2]);
	//float sumwt=0, counts=0;
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<height; i++) {
		int i1 = min((i-(i%pitch))/pitch,hfh-1);
		for (int j=0; j<width; j++) {
			int j1 = min((j-(j%pitch))/pitch,hfw-1);
			
			float pixel[3]={red[i][j],green[i][j],blue[i][j]};
			if (pixel[0]<max_f[0] && pixel[1]<max_f[1] && pixel[2]<max_f[2]) continue;//pixel not clipped
			//if (pixel[0]<fixthresh[0] && pixel[1]<fixthresh[1] && pixel[2]<fixthresh[2]) continue;//pixel not clipped
			
			//%%%%%%%%%%%%%%%%%%%%%%%
			//estimate recovered values using modified HLRecovery_blend algorithm
			
			float rgb[ColorCount], rgb_blend[ColorCount] = {}, cam[2][ColorCount], lab[2][ColorCount], sum[2], chratio;
			
			// Copy input pixel to rgb so it's easier to access in loops
			rgb[0] = pixel[0]; rgb[1] = pixel[1]; rgb[2] = pixel[2];
			
			// Initialize cam with raw input [0] and potentially clipped input [1]
			FOREACHCOLOR {
				cam[0][c] = rgb[c];
				cam[1][c] = min(cam[0][c],clippt);
			}
			
			// Calculate the lightness correction ratio (chratio)
			for (int i2=0; i2<2; i2++) {
				FOREACHCOLOR {
					lab[i2][c]=0;
					for (int j=0; j < ColorCount; j++)
						lab[i2][c] += trans[ColorCount-3][c][j] * cam[i2][j];
				}
				
				sum[i2]=0;
				for (int c=1; c < ColorCount; c++)
					sum[i2] += SQR(lab[i2][c]);
			}
			if(sum[0] == 0.f)		// avoid division by zero
				sum[0] = 0.0001f;
			chratio = sqrt(sum[1]/sum[0]);
			
			
			// Apply ratio to lightness in lab space
			for (int c=1; c < ColorCount; c++) 
				lab[0][c] *= chratio;
			
			// Transform back from lab to RGB
			FOREACHCOLOR {
				cam[0][c]=0;
				for (int j=0; j < ColorCount; j++) {
					cam[0][c] += itrans[ColorCount-3][c][j] * lab[0][j];
				}
			}
			FOREACHCOLOR rgb[c] = cam[0][c] / ColorCount;
			
			// Copy converted pixel back
			float rfrac = min(1.0f,max(1.0f,max_f[0]/medpt)*(pixel[0]-blendpt)/(hlmax[0]-blendpt));
			float gfrac = min(1.0f,max(1.0f,max_f[1]/medpt)*(pixel[1]-blendpt)/(hlmax[1]-blendpt));
			float bfrac = min(1.0f,max(1.0f,max_f[2]/medpt)*(pixel[2]-blendpt)/(hlmax[2]-blendpt));
			if (pixel[0] > blendpt) rgb_blend[0]= rfrac*rgb[0]+(1-rfrac)*pixel[0]; 
			if (pixel[1] > blendpt) rgb_blend[1]= gfrac*rgb[1]+(1-gfrac)*pixel[1]; 
			if (pixel[2] > blendpt) rgb_blend[2]= bfrac*rgb[2]+(1-bfrac)*pixel[2];
			
			//end of HLRecovery_blend estimation
			//%%%%%%%%%%%%%%%%%%%%%%%
			
			//float pixref[3]={min(Yclip,hfsize[0][i1][j1]),min(Yclip,hfsize[1][i1][j1]),min(Yclip,hfsize[2][i1][j1])};
			
			//there are clipped highlights
			//first, determine weighted average of unclipped extensions (weighting is by 'hue' proximity)
			float dirwt, factor, Y;
			float totwt=0;//0.0003;
			float clipfix[3]={0,0,0};//={totwt*rgb_blend[0],totwt*rgb_blend[1],totwt*rgb_blend[2]};
			for (int dir=0; dir<numdirs; dir++) {
				float Yhi = 0.001+(hilite_dir[dir*4+0][i1][j1] + hilite_dir[dir*4+1][i1][j1] + hilite_dir[dir*4+2][i1][j1]);
				float Y   = 0.001+(rgb_blend[0]+rgb_blend[1]+rgb_blend[2]);
				if (hilite_dir[dir*4+0][i1][j1]+hilite_dir[dir*4+1][i1][j1]+hilite_dir[dir*4+2][i1][j1]>0.5) {
					dirwt = invfn[65535*(SQR(rgb_blend[0]/Y-hilite_dir[dir*4+0][i1][j1]/Yhi) +
										 SQR(rgb_blend[1]/Y-hilite_dir[dir*4+1][i1][j1]/Yhi) +
										 SQR(rgb_blend[2]/Y-hilite_dir[dir*4+2][i1][j1]/Yhi))];
					totwt += dirwt;
					clipfix[0] += dirwt*hilite_dir[dir*4+0][i1][j1]/(hilite_dir[dir*4+3][i1][j1]+0.00001);
					clipfix[1] += dirwt*hilite_dir[dir*4+1][i1][j1]/(hilite_dir[dir*4+3][i1][j1]+0.00001);
					clipfix[2] += dirwt*hilite_dir[dir*4+2][i1][j1]/(hilite_dir[dir*4+3][i1][j1]+0.00001);
				}
			}

			if(totwt == 0.f)
				continue;
			
			clipfix[0] /= totwt;
			clipfix[1] /= totwt;
			clipfix[2] /= totwt;
			//sumwt += totwt;
			//counts ++;
					
			//now correct clipped channels
			if (pixel[0]>max_f[0] && pixel[1]>max_f[1] && pixel[2]>max_f[2]) {
				//all channels clipped
				float Y = (0.299*clipfix[0] + 0.587*clipfix[1] + 0.114*clipfix[2]);
				//float Y = (clipfix[0] + clipfix[1] + clipfix[2]);
				factor = whitept/Y;
				red[i][j]   = clipfix[0]*factor;
				green[i][j] = clipfix[1]*factor;
				blue[i][j]  = clipfix[2]*factor;
			} else {//some channels clipped
				int notclipped[3] = {pixel[0]<max_f[0] ? 1 : 0, pixel[1]<max_f[1] ? 1 : 0, pixel[2]<max_f[2] ? 1 : 0};
				
				if (notclipped[0]==0) {//red clipped
					red[i][j]  = max(red[i][j],(clipfix[0]*((notclipped[1]*pixel[1] + notclipped[2]*pixel[2])/
																(notclipped[1]*clipfix[1] + notclipped[2]*clipfix[2]))));
				}
				if (notclipped[1]==0) {//green clipped
					green[i][j] = max(green[i][j],(clipfix[1]*((notclipped[2]*pixel[2] + notclipped[0]*pixel[0])/
																   (notclipped[2]*clipfix[2] + notclipped[0]*clipfix[0]))));
				}
				if (notclipped[2]==0) {//blue clipped
					blue[i][j]  = max(blue[i][j],(clipfix[2]*((notclipped[0]*pixel[0] + notclipped[1]*pixel[1])/
																 (notclipped[0]*clipfix[0] + notclipped[1]*clipfix[1]))));
				}
			}
			
			/*if (hilite[3][i1][j1]>0.01) {
				red[i][j]	= (red[i][j]  + hilite[0][i1][j1])/(1+hilite[3][i1][j1]);
				green[i][j] = (green[i][j]+ hilite[1][i1][j1])/(1+hilite[3][i1][j1]);
				blue[i][j]	= (blue[i][j] + hilite[2][i1][j1])/(1+hilite[3][i1][j1]);
			}*/
			
			Y = (0.299 * red[i][j] + 0.587 * green[i][j] + 0.114 * blue[i][j]);
			
			if (Y>whitept) {
				factor = whitept/Y;
				
				/*I = (0.596 * red[i][j] - 0.275 * green[i][j] - 0.321 * blue[i][j]);
				Q = (0.212 * red[i][j] - 0.523 * green[i][j] + 0.311 * blue[i][j]);
				
				Y *= factor;
				I *= factor;//max(0,min(1,(whitept-Y)/(whitept-clippt)));
				Q *= factor;//max(0,min(1,(whitept-Y)/(whitept-clippt)));
				
				red[i][j]   = Y + 0.956*I + 0.621*Q;
				green[i][j] = Y - 0.272*I - 0.647*Q;
				blue[i][j]  = Y - 1.105*I + 1.702*Q;*/
				
				red[i][j]   *= factor;
				green[i][j] *= factor;
				blue[i][j]  *= factor;
			}
		}
	}
	
	if(plistener) plistener->setProgress(1.00);
	//printf("ave wt=%f\n",sumwt/counts);

	
	// diagnostic output
	/*for (int i=0; i<height; i++) {
		int i1 = min(hfh-1,(i-(i%pitch))/pitch);
		for (int j=0; j<width; j++) {
			int j1 = min(hfw-1,(j-(j%pitch))/pitch);
			
			//red[i][j]  =hfsize[0][i1][j1];
			//green[i][j]=hfsize[1][i1][j1];
			//blue[i][j] =hfsize[2][i1][j1];
	 
			//red[i][j]  =clippt/2+red[i][j]-channelblur[0][i][j];
			//green[i][j]=clippt/2+green[i][j]-channelblur[1][i][j];
			//blue[i][j] =clippt/2+blue[i][j]-channelblur[2][i][j];
			
			red[i][j]=  hilite[0][i1][j1]/(hilite[3][i1][j1]+0.001);
			green[i][j]=hilite[1][i1][j1]/(hilite[3][i1][j1]+0.001);
			blue[i][j]= hilite[2][i1][j1]/(hilite[3][i1][j1]+0.001);
			
			//red[i][j]=  hilite_dir[0+0][i1][j1]/hilite_dir[0+3][i1][j1];
			//green[i][j]=hilite_dir[0+1][i1][j1]/hilite_dir[0+3][i1][j1];
			//blue[i][j]= hilite_dir[0+2][i1][j1]/hilite_dir[0+3][i1][j1];
		
			//red[i][j]=  clipfix[0][i1][j1];
			//green[i][j]=clipfix[1][i1][j1];
			//blue[i][j]= clipfix[2][i1][j1];
			
		}
	}*/
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
}// end of HLReconstruction
	
/*
	void RawImageSource::halfsize()
	{
		int ex,ey;
		//determine GRBG coset; (ey,ex) is the offset of the R subarray
		if (FC(0,0)==1) {//first pixel is G
			if (FC(0,1)==0) {ey=0; ex=1;} else {ey=1; ex=0;}
		} else {//first pixel is R or B
			if (FC(0,0)==0) {ey=0; ex=0;} else {ey=1; ex=1;}
		}
		
		for (int i=0; i<(H-(H&1)); i++) {
			for (int j=0; j<(W-(W&1)); j++){
				red[i][j] = rawData[i+ey-(i&1)][j+ex-(j&1)];
				green[i][j] = (rawData[i+(1-ey)-(i&1)][j+ex-(j&1)]+rawData[i+ey-(i&1)][j+(1-ex)-(j&1)])/2;
				blue[i][j] = rawData[i+(1-ey)-(i&1)][j+(1-ex)-(j&1)];
			}
			if (W&1) {
				red[i][W-1]=red[i][W-2];
				green[i][W-1]=green[i][W-2];
				blue[i][W-1]=blue[i][W-2];
			}
		}
		if (H&1) {
			for (int j=0; j<W; j++){
				red[H-1][j] = red[H-2][j];
				green[H-1][j] = green[H-2][j];
				blue[H-1][j] = blue[H-2][j];
			}
		}
	}
*/
	
}

