/*
 *  This file is part of RawTherapee.
 *
 *  Copyright © 2010 Emil Martinec <ejmartin@uchicago.edu>
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
#ifndef _BOXBLUR_H_
#define _BOXBLUR_H_

#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include "alignedbuffer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define SQR(x) ((x)*(x))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// classical filtering if the support window is small:

template<class T, class A> void boxblur (T** src, A** dst, int radx, int rady, int W, int H) {
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//box blur image; box range = (radx,rady)
	
	AlignedBuffer<float>* buffer = new AlignedBuffer<float> (W*H);
	float* temp = buffer->data;
	
	if (radx==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				temp[row*H+col] = (float)src[row][col];
			}
	} else {
		//horizontal blur
//OpenMP here		
		for (int row = 0; row < H; row++) {
			int len = radx + 1;
			temp[row*W+0] = (float)src[row][0]/len;
			for (int j=1; j<=radx; j++) {
				temp[row*W+0] += (float)src[row][j]/len;
			}
			for (int col=1; col<=radx; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len + (float)src[row][col+radx])/(len+1);
				len ++;
			}
			for (int col = radx+1; col < W-radx; col++) {
				temp[row*W+col] = temp[row*W+col-1] + ((float)(src[row][col+radx] - src[row][col-radx-1]))/len;
			}
			for (int col=W-radx; col<W; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len - src[row][col-radx-1])/(len-1);
				len --;
			}
		}
	}
	
	if (rady==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				dst[row][col] = temp[row*W+col];
			}
	} else {
		//vertical blur
//OpenMP here		
		for (int col = 0; col < W; col++) {
			int len = rady + 1;
			dst[0][col] = temp[0*W+col]/len;
			for (int i=1; i<=rady; i++) {
				dst[0][col] += temp[i*W+col]/len;
			}
			for (int row=1; row<=rady; row++) {
				dst[row][col] = (dst[(row-1)][col]*len + temp[(row+rady)*W+col])/(len+1);
				len ++;
			}
			for (int row = rady+1; row < H-rady; row++) {
				dst[row][col] = dst[(row-1)][col] + (temp[(row+rady)*W+col] - temp[(row-rady-1)*W+col])/len;
			}
			for (int row=H-rady; row<H; row++) {
				dst[row][col] = (dst[(row-1)][col]*len - temp[(row-rady-1)*W+col])/(len-1);
				len --;
			}
		}
	}
	
	delete buffer;
	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<class T, class A> void boxblur (T* src, A* dst, int radx, int rady, int W, int H) {
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//box blur image; box range = (radx,rady) i.e. box size is (2*radx+1)x(2*rady+1)
	
	AlignedBuffer<float>* buffer = new AlignedBuffer<float> (W*H);
	float* temp = buffer->data;
	
	if (radx==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				temp[row*H+col] = src[row*W+col];
			}
	} else {
		//horizontal blur
//OpenMP here		
		for (int row = 0; row < H; row++) {
			int len = radx + 1;
			temp[row*W+0] = (float)src[row*W+0]/len;
			for (int j=1; j<=radx; j++) {
				temp[row*W+0] += (float)src[row*W+j]/len;
			}
			for (int col=1; col<=radx; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len + src[row*W+col+radx])/(len+1);
				len ++;
			}
			for (int col = radx+1; col < W-radx; col++) {
				temp[row*W+col] = temp[row*W+col-1] + ((float)(src[row*W+col+radx] - src[row*W+col-radx-1]))/len;
			}
			for (int col=W-radx; col<W; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len - src[row*W+col-radx-1])/(len-1);
				len --;
			}
		}
	}
	
	if (rady==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				dst[row*W+col] = temp[row*W+col];
			}
	} else {
		//vertical blur
//OpenMP here		
		for (int col = 0; col < W; col++) {
			int len = rady + 1;
			dst[0*W+col] = temp[0*W+col]/len;
			for (int i=1; i<=rady; i++) {
				dst[0*W+col] += temp[i*W+col]/len;
			}
			for (int row=1; row<=rady; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len + temp[(row+rady)*W+col])/(len+1);
				len ++;
			}
			for (int row = rady+1; row < H-rady; row++) {
				dst[row*W+col] = dst[(row-1)*W+col] + (temp[(row+rady)*W+col] - temp[(row-rady-1)*W+col])/len;
			}
			for (int row=H-rady; row<H; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len - temp[(row-rady-1)*W+col])/(len-1);
				len --;
			}
		}
	}
	
	delete buffer;
	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<typename T> void boxvar (T* src, T* dst, int radx, int rady, int W, int H) {
	
	AlignedBuffer<float>* buffer1 = new AlignedBuffer<float> (W*H);
	AlignedBuffer<float>* buffer2 = new AlignedBuffer<float> (W*H);
	float* tempave = buffer1->data;
	float* tempsqave = buffer2->data;
	
	//float *tempave2 = new float[H];
	AlignedBuffer<float>* buffer3 = new AlignedBuffer<float> (H);
	float* tempave2 = buffer3->data;

	//float image_ave = 0;
	
	//box blur image channel; box size = 2*box+1
	//horizontal blur
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int row = 0; row < H; row++) {
		int len = radx + 1;
		tempave[row*W+0] = src[row*W+0]/len;
		tempsqave[row*W+0] = SQR(src[row*W+0])/len;
		for (int j=1; j<=radx; j++) {
			tempave[row*W+0] += src[row*W+j]/len;
			tempsqave[row*W+0] += SQR(src[row*W+j])/len;
		}
		for (int col=1; col<=radx; col++) {
			tempave[row*W+col] = (tempave[row*W+col-1]*len + src[row*W+col+radx])/(len+1);
			tempsqave[row*W+col] = (tempsqave[row*W+col-1]*len + SQR(src[row*W+col+radx]))/(len+1);
			len ++;
		}
		for (int col = radx+1; col < W-radx; col++) {
			tempave[row*W+col] = tempave[row*W+col-1] + (src[row*W+col+radx] - src[row*W+col-radx-1])/len;
			tempsqave[row*W+col] = tempsqave[row*W+col-1] + (SQR(src[row*W+col+radx]) - SQR(src[row*W+col-radx-1]))/len;
		}
		for (int col=W-radx; col<W; col++) {
			tempave[row*W+col] = (tempave[row*W+col-1]*len - src[row*W+col-radx-1])/(len-1);
			tempsqave[row*W+col] = (tempsqave[row*W+col-1]*len - SQR(src[row*W+col-radx-1]))/(len-1);
			len --;
		}
	}
	
	//vertical blur
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int col = 0; col < W; col++) {
		int len = rady + 1;
		tempave2[0] = tempave[0*W+col]/len;
		dst[0*W+col] = tempsqave[0*W+col]/len;
		for (int i=1; i<=rady; i++) {
			tempave2[0] += tempave[i*W+col]/len;
			dst[0*W+col] += tempsqave[i*W+col]/len;
		}
		for (int row=1; row<=rady; row++) {
			tempave2[row] = (tempave2[(row-1)]*len + tempave[(row+rady)*W+col])/(len+1);
			dst[row*W+col] = (dst[(row-1)*W+col]*len + tempsqave[(row+rady)*W+col])/(len+1);
			len ++;
		}
		for (int row = rady+1; row < H-rady; row++) {
			tempave2[row] = tempave2[(row-1)] + (tempave[(row+rady)*W+col] - tempave[(row-rady-1)*W+col])/len;
			dst[row*W+col] = dst[(row-1)*W+col] + (tempsqave[(row+rady)*W+col] - tempsqave[(row-rady-1)*W+col])/len;
		}
		for (int row=H-rady; row<H; row++) {
			tempave2[row] = (tempave2[(row-1)]*len - tempave[(row-rady-1)*W+col])/(len-1);
			dst[row*W+col] = (dst[(row-1)*W+col]*len - tempsqave[(row-rady-1)*W+col])/(len-1);
			len --;
		}
		//now finish off
		for (int row=0; row<H; row++) {
			dst[row*W+col] = fabs(dst[row*W+col] - SQR(tempave2[row]));
			//image_ave += src[row*W+col];
		}
	}
	
	//image_ave /= (W*H);
	
	delete buffer1;
	delete buffer2;
	delete buffer3;
	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<typename T> void boxdev (T* src, T* dst, int radx, int rady, int W, int H) {
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//box blur image; box range = (radx,rady) i.e. box size is (2*radx+1)x(2*rady+1)
	
	AlignedBuffer<float>* buffer1 = new AlignedBuffer<float> (W*H);
	float* temp = buffer1->data;
	
	AlignedBuffer<float>* buffer2 = new AlignedBuffer<float> (W*H);
	float* tempave = buffer2->data;
	
	if (radx==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				temp[row*H+col] = src[row*W+col];
			}
	} else {
		//horizontal blur
//OpenMP here		
		for (int row = 0; row < H; row++) {
			int len = radx + 1;
			temp[row*W+0] = (float)src[row*W+0]/len;
			for (int j=1; j<=radx; j++) {
				temp[row*W+0] += (float)src[row*W+j]/len;
			}
			for (int col=1; col<=radx; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len + src[row*W+col+radx])/(len+1);
				len ++;
			}
			for (int col = radx+1; col < W-radx; col++) {
				temp[row*W+col] = temp[row*W+col-1] + ((float)(src[row*W+col+radx] - src[row*W+col-radx-1]))/len;
			}
			for (int col=W-radx; col<W; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len - src[row*W+col-radx-1])/(len-1);
				len --;
			}
		}
	}
	
	if (rady==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				tempave[row*W+col] = temp[row*W+col];
			}
	} else {
		//vertical blur
//OpenMP here		
		for (int col = 0; col < W; col++) {
			int len = rady + 1;
			tempave[0*W+col] = temp[0*W+col]/len;
			for (int i=1; i<=rady; i++) {
				tempave[0*W+col] += temp[i*W+col]/len;
			}
			for (int row=1; row<=rady; row++) {
				tempave[row*W+col] = (tempave[(row-1)*W+col]*len + temp[(row+rady)*W+col])/(len+1);
				len ++;
			}
			for (int row = rady+1; row < H-rady; row++) {
				tempave[row*W+col] = tempave[(row-1)*W+col] + (temp[(row+rady)*W+col] - temp[(row-rady-1)*W+col])/len;
			}
			for (int row=H-rady; row<H; row++) {
				tempave[row*W+col] = (tempave[(row-1)*W+col]*len - temp[(row-rady-1)*W+col])/(len-1);
				len --;
			}
		}
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//box blur absolute deviation

	
	if (radx==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				temp[row*H+col] = fabs(src[row*W+col]-tempave[row*W+col]);
			}
	} else {
		//horizontal blur
//OpenMP here		
		for (int row = 0; row < H; row++) {
			int len = radx + 1;
			temp[row*W+0] = fabs(src[row*W+0]-tempave[row*W+0])/len;
			for (int j=1; j<=radx; j++) {
				temp[row*W+0] += fabs(src[row*W+j]-tempave[row*W+j])/len;
			}
			for (int col=1; col<=radx; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len + fabs(src[row*W+col+radx]-tempave[row*W+col+radx]))/(len+1);
				len ++;
			}
			for (int col = radx+1; col < W-radx; col++) {
				temp[row*W+col] = temp[row*W+col-1] + (fabs(src[row*W+col+radx]-tempave[row*W+col+radx]) - \
													   fabs(src[row*W+col-radx-1]-tempave[row*W+col-radx-1]))/len;
			}
			for (int col=W-radx; col<W; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len - fabs(src[row*W+col-radx-1]-tempave[row*W+col-radx-1]))/(len-1);
				len --;
			}
		}
	}
	
	if (rady==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				dst[row*W+col] = temp[row*W+col];
			}
	} else {
		//vertical blur
//OpenMP here		
		for (int col = 0; col < W; col++) {
			int len = rady + 1;
			dst[0*W+col] = temp[0*W+col]/len;
			for (int i=1; i<=rady; i++) {
				dst[0*W+col] += temp[i*W+col]/len;
			}
			for (int row=1; row<=rady; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len + temp[(row+rady)*W+col])/(len+1);
				len ++;
			}
			for (int row = rady+1; row < H-rady; row++) {
				dst[row*W+col] = dst[(row-1)*W+col] + (temp[(row+rady)*W+col] - temp[(row-rady-1)*W+col])/len;
			}
			for (int row=H-rady; row<H; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len - temp[(row-rady-1)*W+col])/(len-1);
				len --;
			}
		}
	}
	
	delete buffer1;
	delete buffer2;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<class T, class A> void boxsqblur (T* src, A* dst, int radx, int rady, int W, int H) {
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//box blur image; box range = (radx,rady) i.e. box size is (2*radx+1)x(2*rady+1)
	
	AlignedBuffer<float>* buffer = new AlignedBuffer<float> (W*H);
	float* temp = buffer->data;
	
	if (radx==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				temp[row*H+col] = SQR(src[row*W+col]);
			}
	} else {
		//horizontal blur
//OpenMP here		
		for (int row = 0; row < H; row++) {
			int len = radx + 1;
			temp[row*W+0] = SQR((float)src[row*W+0])/len;
			for (int j=1; j<=radx; j++) {
				temp[row*W+0] += SQR((float)src[row*W+j])/len;
			}
			for (int col=1; col<=radx; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len + SQR(src[row*W+col+radx]))/(len+1);
				len ++;
			}
			for (int col = radx+1; col < W-radx; col++) {
				temp[row*W+col] = temp[row*W+col-1] + ((float)(SQR(src[row*W+col+radx]) - SQR(src[row*W+col-radx-1])))/len;
			}
			for (int col=W-radx; col<W; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len - SQR(src[row*W+col-radx-1]))/(len-1);
				len --;
			}
		}
	}
	
	if (rady==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				dst[row*W+col] = temp[row*W+col];
			}
	} else {
		//vertical blur
//OpenMP here		
		for (int col = 0; col < W; col++) {
			int len = rady + 1;
			dst[0*W+col] = temp[0*W+col]/len;
			for (int i=1; i<=rady; i++) {
				dst[0*W+col] += temp[i*W+col]/len;
			}
			for (int row=1; row<=rady; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len + temp[(row+rady)*W+col])/(len+1);
				len ++;
			}
			for (int row = rady+1; row < H-rady; row++) {
				dst[row*W+col] = dst[(row-1)*W+col] + (temp[(row+rady)*W+col] - temp[(row-rady-1)*W+col])/len;
			}
			for (int row=H-rady; row<H; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len - temp[(row-rady-1)*W+col])/(len-1);
				len --;
			}
		}
	}
	
	delete buffer;
	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<class T, class A> void boxcorrelate (T* src, A* dst, int dx, int dy, int radx, int rady, int W, int H) {
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//box blur image; box range = (radx,rady) i.e. box size is (2*radx+1)x(2*rady+1)
	
	AlignedBuffer<float>* buffer = new AlignedBuffer<float> (W*H);
	float* temp = buffer->data;
	
	if (radx==0) {
		for (int row=0; row<H; row++) {
			int rr = MIN(H-1,MAX(0,row+dy));
			for (int col=0; col<H; col++) {
				int cc = MIN(W-1,MAX(0,col+dx));
				temp[row*H+col] = dy>0 ? (src[row*W+col])*(src[rr*W+cc]) : 0;
			}
		}
	} else {
		//horizontal blur
//OpenMP here		
		for (int row = 0; row < H; row++) {
			int len = radx + 1;
			int rr = MIN(H-1,MAX(0,row+dy));
			int cc = MIN(W-1,MAX(0,0+dx));
			temp[row*W+0] = ((float)src[row*W+0])*(src[rr*W+cc])/len;
			for (int j=1; j<=radx; j++) {
				int cc = MIN(W-1,MAX(0,j+dx));
				temp[row*W+0] += ((float)src[row*W+j])*(src[rr*W+cc])/len;
			}
			for (int col=1; col<=radx; col++) {
				int cc = MIN(W-1,MAX(0,col+dx+radx));
				temp[row*W+col] = (temp[row*W+col-1]*len + (src[row*W+col+radx])*(src[rr*W+cc]))/(len+1);
				len ++;
			}
			for (int col = radx+1; col < W-radx; col++) {
				int cc = MIN(W-1,MAX(0,col+dx+radx));
				int cc1 = MIN(W-1,MAX(0,col+dx-radx-1));
				temp[row*W+col] = temp[row*W+col-1] + ((float)((src[row*W+col+radx])*(src[rr*W+cc]) - 
															   (src[row*W+col-radx-1])*(src[rr*W+cc1])))/len;
			}
			for (int col=W-radx; col<W; col++) {
				int cc1 = MIN(W-1,MAX(0,col+dx-radx-1));
				temp[row*W+col] = (temp[row*W+col-1]*len - (src[row*W+col-radx-1])*(src[rr*W+cc1]))/(len-1);
				len --;
			}
		}
	}
	
	if (rady==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				dst[row*W+col] = temp[row*W+col];
			}
	} else {
		//vertical blur
//OpenMP here		
		for (int col = 0; col < W; col++) {
			int len = rady + 1;
			dst[0*W+col] = temp[0*W+col]/len;
			for (int i=1; i<=rady; i++) {
				dst[0*W+col] += temp[i*W+col]/len;
			}
			for (int row=1; row<=rady; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len + temp[(row+rady)*W+col])/(len+1);
				len ++;
			}
			for (int row = rady+1; row < H-rady; row++) {
				dst[row*W+col] = dst[(row-1)*W+col] + (temp[(row+rady)*W+col] - temp[(row-rady-1)*W+col])/len;
			}
			for (int row=H-rady; row<H; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len - temp[(row-rady-1)*W+col])/(len-1);
				len --;
			}
		}
	}
	
	delete buffer;
	
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<class T, class A> void boxabsblur (T* src, A* dst, int radx, int rady, int W, int H) {
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//box blur image; box range = (radx,rady) i.e. box size is (2*radx+1)x(2*rady+1)
	
	AlignedBuffer<float>* buffer = new AlignedBuffer<float> (W*H);
	float* temp = buffer->data;
	
	if (radx==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				temp[row*H+col] = fabs(src[row*W+col]);
			}
	} else {
		//horizontal blur
//OpenMP here		
		for (int row = 0; row < H; row++) {
			int len = radx + 1;
			temp[row*W+0] = fabs((float)src[row*W+0])/len;
			for (int j=1; j<=radx; j++) {
				temp[row*W+0] += fabs((float)src[row*W+j])/len;
			}
			for (int col=1; col<=radx; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len + fabs(src[row*W+col+radx]))/(len+1);
				len ++;
			}
			for (int col = radx+1; col < W-radx; col++) {
				temp[row*W+col] = temp[row*W+col-1] + ((float)(fabs(src[row*W+col+radx]) - fabs(src[row*W+col-radx-1])))/len;
			}
			for (int col=W-radx; col<W; col++) {
				temp[row*W+col] = (temp[row*W+col-1]*len - fabs(src[row*W+col-radx-1]))/(len-1);
				len --;
			}
		}
	}
	
	if (rady==0) {
		for (int row=0; row<H; row++) 
			for (int col=0; col<H; col++) {
				dst[row*W+col] = temp[row*W+col];
			}
	} else {
		//vertical blur
//OpenMP here		
		for (int col = 0; col < W; col++) {
			int len = rady + 1;
			dst[0*W+col] = temp[0*W+col]/len;
			for (int i=1; i<=rady; i++) {
				dst[0*W+col] += temp[i*W+col]/len;
			}
			for (int row=1; row<=rady; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len + temp[(row+rady)*W+col])/(len+1);
				len ++;
			}
			for (int row = rady+1; row < H-rady; row++) {
				dst[row*W+col] = dst[(row-1)*W+col] + (temp[(row+rady)*W+col] - temp[(row-rady-1)*W+col])/len;
			}
			for (int row=H-rady; row<H; row++) {
				dst[row*W+col] = (dst[(row-1)*W+col]*len - temp[(row-rady-1)*W+col])/(len-1);
				len --;
			}
		}
	}
	
	delete buffer;
	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#endif /* _BOXBLUR_H_ */