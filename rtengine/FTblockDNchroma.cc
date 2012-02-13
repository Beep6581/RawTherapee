////////////////////////////////////////////////////////////////
//
//			CFA chroma denoise by FT filtering
//
//	copyright (c) 2008-2012  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: December 24, 2011
//
//	FTblockDNchroma.cc is free software: you can redistribute it and/or modify
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
	
	
#include <math.h>
#include <fftw3.h>

//#include "bilateral2.h"
//#include "gauss.h"

#include "rtengine.h"
#include "improcfun.h"
#include "LUT.h"
#include "array2D.h"
#include "iccmatrices.h"
#include "boxblur.h"

#ifdef _OPENMP
#include <omp.h>
#endif


#include "EdgePreserveLab.h"
#include "cplx_wavelet_dec.h"


#define SQR(x) ((x)*(x))
//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
//#define LIM(x,min,max) MAX(min,MIN(x,max))
//#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
//#define CLIP(x) LIM(x,0,65535)

#define TS 256		// Tile size
#define offset (TS/4)	// shift between tiles
#define fTS ((TS/2+1))	// second dimension of Fourier tiles
//#define eps 0.01f/(TS*TS) //tolerance

namespace rtengine {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
 Structure of the algorithm: 
 
 1. Compute a high pass filter of the image via bilateral filter with user input range
 2. Decompose the image into TSxTS size tiles, shifting by 'offset' each step (so roughly 
	each pixel is in (TS/offset)^2 tiles); Fourier transform the tiles after applying a mask
	to prevent long range tails in the FT data due to boundary discontinuities.
 3. Compute the average size of Fourier coefficients.
 4. Damp the FT data of the tile by a Wiener filter factor
			(image_variance)/(image_variance + noise_control)
	where noise_control is the user specified noise reduction amount.  
	Noise_control is altered according to neighbor average.
 6. Inverse FT the denoised tile data and combine the tiles into a denoised output
	image.
 
*/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGB_InputTransf(Imagefloat * src, LabImage * dst, const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe) {
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// gamma transform input channel data
		float gam = dnparams.gamma;
		float gamthresh = 0.03;
		float gamslope = exp(log((double)gamthresh)/gam)/gamthresh;
		
		LUTf gamcurve(65536,0);
		
		for (int i=0; i<65536; i++) {
			gamcurve[i] = (gamma((double)i/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)) * 32768.0f;
		}
		
		//srand((unsigned)time(0));//test with random data
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<src->height; i++) 
			for (int j=0; j<src->width; j++) {
				
				float X = src->r[i][j];//xyz_prophoto[0][0]*src->r[i][j] + xyz_prophoto[0][1]*src->g[i][j] + xyz_prophoto[0][2]*src->b[i][j];
				float Y = src->g[i][j];//xyz_prophoto[1][0]*src->r[i][j] + xyz_prophoto[1][1]*src->g[i][j] + xyz_prophoto[1][2]*src->b[i][j];
				float Z = src->b[i][j];//xyz_prophoto[2][0]*src->r[i][j] + xyz_prophoto[2][1]*src->g[i][j] + xyz_prophoto[2][2]*src->b[i][j];
				
				X = X<65535.0f ? gamcurve[X] : (gamma((double)X/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
				Y = Y<65535.0f ? gamcurve[Y] : (gamma((double)Y/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
				Z = Z<65535.0f ? gamcurve[Z] : (gamma((double)Z/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
				
				dst->L[i][j] = Y;
				dst->a[i][j] = 0.2f*(X-Y);
				dst->b[i][j] = 0.2f*(Y-Z);
				
				//Y = 0.05+0.1*((float)rand()/(float)RAND_MAX);//test with random data
				//dst->L[i][j] = gamcurve[65535.0f*Y];
				
			}

		
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGB_OutputTransf(LabImage * src, Imagefloat * dst, const procparams::DirPyrDenoiseParams & dnparams) {
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// gamma transform output channel data
		float gam = dnparams.gamma;
		float gamthresh = 0.03;
		float gamslope = exp(log((double)gamthresh)/gam)/gamthresh;
		float igam = 1/gam;
		float igamthresh = gamthresh*gamslope;
		float igamslope = 1/gamslope;
		
		LUTf igamcurve(65536,0);
		
		for (int i=0; i<65536; i++) {
			igamcurve[i] = (gamma((float)i/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
		}
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<src->H; i++) {
			float X,Y,Z;
			for (int j=0; j<src->W; j++) {
				//input normalized to (0,1)
				//Y = igamcurveL[ src->L[i][j] ];
				Y = src->L[i][j];
				X = (5.0f*(src->a[i][j])) + Y;
				Z = Y - (5.0f*(src->b[i][j]));
				
				X = X<32768.0f ? igamcurve[X] : (gamma((float)X/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
				Y = Y<32768.0f ? igamcurve[Y] : (gamma((float)Y/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
				Z = Z<32768.0f ? igamcurve[Z] : (gamma((float)Z/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
				
				//Y = 65535.0f*(0.05+0.1*((float)rand()/(float)RAND_MAX));//test with random data
				
				dst->r[i][j] = X;//prophoto_xyz[0][0]*X + prophoto_xyz[0][1]*Y + prophoto_xyz[0][2]*Z;
				dst->g[i][j] = Y;//prophoto_xyz[1][0]*X + prophoto_xyz[1][1]*Y + prophoto_xyz[1][2]*Z;
				dst->b[i][j] = Z;//prophoto_xyz[2][0]*X + prophoto_xyz[2][1]*Y + prophoto_xyz[2][2]*Z;
				
			}
		}
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGB_denoise(Imagefloat * src, Imagefloat * dst, /*int Roffset,*/ const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe)
	{	
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		/*if (plistener) {
		 plistener->setProgressStr ("Block FT Luma Denoise...");
		 plistener->setProgress (0.0);
		 }*/
		
		volatile double progress = 0.0;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		const short int height=src->height, width=src->width;
		const short int hfh=(height+1)/2, hfw=(width+1)/2;
		//const short int hfh=(height+1), hfw=(width+1);
		
		const int blkrad=1;
		float noisevar_L = SQR(dnparams.luma * TS * 10.0f);
		float noisevar_ab = SQR(dnparams.chroma * TS * 100.0f);
		
		//int dxr=Roffset&1, dyr=(Roffset&2)/2, dxb=(1-dxr), dyb=(1-dyr);
		//int rdx, rdy, bdx, bdy;
		
		// calculation for tiling
		const int numblox_W = ceil(((float)(width))/(offset))+2*blkrad;
		const int numblox_H = ceil(((float)(height))/(offset))+2*blkrad;
		//const int nrtiles = numblox_W*numblox_H;
		// end of tiling calc
		
		//const float eps = 1.0f;
		
		
		
		array2D<float> tilemask_in(TS,TS);
		array2D<float> tilemask_out(TS,TS);
		
		array2D<float> totwt(width,height,ARRAY2D_CLEAR_DATA);
		
		for (int i=0; i<TS; i++) {
			float i1 = abs((i>TS/2 ? i-TS : i));
			float vmask = (i1<8 ? SQR(sin(M_PI*i1/16.0f)) : 1.0f);
			float vmask2 = (i1<32 ? SQR(sin(M_PI*i1/64.0f)) : 1.0f);
			for (int j=0; j<TS; j++) {
				float j1 = abs((j>TS/2 ? j-TS : j));
				//tilemask_in[i][j] = vmask * (j1<4 ? SQR(sin(M_PI*(float)j1/8.0f)) : 1.0f);
				tilemask_in[i][j] = (exp(-(SQR(i-TS/2-0.5f)+SQR(j-TS/2-0.5f))/(2*SQR(TS/4.0f))) * \
									 vmask * (j1<8 ? SQR(sin(M_PI*(float)j1/16.0f)) : 1.0f));
				
				//tilemask_out[i][j] = (SQR(MAX(0.0f, sin(M_PI*(float)(i-8.0f)/(TS-17) ))) * \
				SQR(MAX(0.0f, sin(M_PI*(float)(j-8.0f)/(TS-17) ))));
				tilemask_out[i][j] = vmask2 * (j1<32 ? SQR(sin(M_PI*(float)j1/64.0f)) : 1.0f);
				
				//tilemask_out[i][j] = exp(-(SQR(i-TS/2-0.5f)+SQR(j-TS/2-0.5f))/(SQR(TS/4.0f)));
				
			}
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		LabImage * labin = new LabImage(width,height);
		LabImage * labblur = new LabImage(width,height);
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// transform RGB input to ersatz Lab
		
		RGB_InputTransf(src, labin, dnparams, defringe);
		
		int datalen = labin->W * labin->H;
		
		wavelet_decomposition Ldecomp(labin->data, labin->W, labin->H, 5/*maxlevels*/, 0/*subsampling*/ );
		wavelet_decomposition adecomp(labin->data+datalen, labin->W, labin->H, 5, 1 );//last args are maxlevels, subsampling
		wavelet_decomposition bdecomp(labin->data+2*datalen, labin->W, labin->H, 5, 1 );//last args are maxlevels, subsampling

		float noisevarL	 = SQR(dnparams.Lamt/25.0f);//TODO: clean up naming confusion about params
		float noisevarab = SQR(dnparams.chroma/25.0f);
		
		//WaveletDenoise(Ldecomp, SQR((float)dnparams.Lamt/25.0f));
		WaveletDenoiseAll(Ldecomp, adecomp, bdecomp, noisevarL, noisevarab);

		Ldecomp.reconstruct(labblur->data);
		adecomp.reconstruct(labblur->data+datalen);
		bdecomp.reconstruct(labblur->data+2*datalen);

		
		//impulse_nr (dst, 50.0f/20.0f);
		//PF_correct_RT(dst, dst, defringe.radius, defringe.threshold);
		
		
		//dirpyr_ab(labin, labblur, dnparams);//use dirpyr here if using it to blur ab channels only
		//dirpyrLab_denoise(labin, labblur, dnparams);//use dirpyr here if using it to blur ab channels only

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// initialize FT and output data structures
		LabImage * labdn = new LabImage(width,height);
		for (int i=0; i<3*labdn->W*labdn->H; i++) {
			labdn->data[i] = 0.0f;
		}
		float ** Lblox = new float *[8] ;
		float ** RLblox = new float *[8] ;
		float ** BLblox = new float *[8] ;
		
		fftwf_complex ** fLblox = new fftwf_complex *[8] ;	//for FT
		//float ** fLblox = new float *[8] ;				//for DCT
		
		fftwf_complex ** fRLblox = new fftwf_complex *[8] ;
		fftwf_complex ** fBLblox = new fftwf_complex *[8] ;
		
		for( int i = 0 ; i < 8 ; i++ ) {
			Lblox[i] = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));
			//RLblox[i] = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));
			//BLblox[i] = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));
			
			fLblox[i] = (fftwf_complex *) fftwf_malloc (numblox_W*TS*fTS * sizeof (fftwf_complex));	//for FT
			//fLblox[i] = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));				//for DCT
			
			//fRLblox[i] = (fftwf_complex *) fftwf_malloc (numblox_W*TS*fTS * sizeof (fftwf_complex));
			//fBLblox[i] = (fftwf_complex *) fftwf_malloc (numblox_W*TS*fTS * sizeof (fftwf_complex));
		}
		
		//make a plan for FFTW
		fftwf_plan plan_forward_blox, plan_backward_blox;
		
		int nfwd[2]={TS,TS};
		
		//for FT:
		plan_forward_blox  = fftwf_plan_many_dft_r2c(2, nfwd, numblox_W, Lblox[0], NULL, 1, TS*TS, \
													 fLblox[0], NULL, 1, TS*fTS, FFTW_ESTIMATE );
		plan_backward_blox = fftwf_plan_many_dft_c2r(2, nfwd, numblox_W, fLblox[0], NULL, 1, TS*fTS, \
													 Lblox[0], NULL, 1, TS*TS, FFTW_ESTIMATE );
		
		//for DCT:
		//const fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
		//const fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};
		
		//plan_forward_blox  = fftwf_plan_many_r2r(2, nfwd, numblox_W, Lblox[0], NULL, 1, TS*TS, \
		fLblox[0], NULL, 1, TS*TS, fwdkind, FFTW_ESTIMATE );
		//plan_backward_blox = fftwf_plan_many_r2r(2, nfwd, numblox_W, fLblox[0], NULL, 1, TS*TS, \
		Lblox[0], NULL, 1, TS*TS, bwdkind, FFTW_ESTIMATE );
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Main algorithm: Tile loop
#pragma omp parallel for schedule(dynamic)
		//int vblock=0, hblock=0;
		for (int vblk=0; vblk<numblox_H; vblk++) {
			//printf("vblock=%d",vblk);
			int vblkmod = vblk%8;
			
			int top = (vblk-blkrad)*offset;
			
			for (int hblk=0; hblk<numblox_W; hblk++) {
				int left = (hblk-blkrad)*offset;
				int bottom = MIN( top+TS,height);
				int right  = MIN(left+TS, width);
				int imin = MIN(TS-1,MAX(0,-top));
				int jmin = MIN(TS-1,MAX(0,-left));
				int imax = MAX(0,bottom - top);
				int jmax = MAX(0,right - left);
				
				int indx = (hblk)*TS;//index of block in malloc
				
				// load Lab high pass data
				
				for (int i=imin; i<imax; i++)
					for (int j=jmin; j<jmax; j++) {
						
						Lblox[vblkmod][(indx + i)*TS+j]  = tilemask_in[i][j]*(labin->L[top+i][left+j]-labblur->L[top+i][left+j]);// luma data
						//RLblox[vblkmod][(indx + i)*TS+j] = tilemask_in[i][j]*(labin->a[top+i][left+j]-labblur->a[top+i][left+j]);// high pass chroma data
						//BLblox[vblkmod][(indx + i)*TS+j] = tilemask_in[i][j]*(labin->b[top+i][left+j]-labblur->b[top+i][left+j]);// high pass chroma data
						
						totwt[top+i][left+j] += tilemask_in[i][j]*tilemask_out[i][j];
					}
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				//pad the image to size TS on both sides
				
				if (imin>0) {
					for (int i=0; i<imin; i++) 
						for (int j=jmin; j<jmax; j++) {
							Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[2*imin-i-1][left+j]-labblur->L[2*imin-i-1][left+j]);
							//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[2*imin-i-1][left+j]-labblur->a[2*imin-i-1][left+j]);
							//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[2*imin-i-1][left+j]-labblur->b[2*imin-i-1][left+j]);
						}
					if (jmin>0) {
						for (int i=0; i<imin; i++) 
							for (int j=0; j<jmin; j++) {
								Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[2*imin-i-1][2*jmin-j-1]-labblur->L[2*imin-i-1][2*jmin-j-1]);
								//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[2*imin-i-1][2*jmin-j-1]-labblur->a[2*imin-i-1][2*jmin-j-1]);
								//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[2*imin-i-1][2*jmin-j-1]-labblur->b[2*imin-i-1][2*jmin-j-1]);
							}
					}
				}
				
				if (imax<TS) {
					for (int i=imax; i<TS; i++) 
						for (int j=jmin; j<jmax; j++) {
							Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[height+imax-i-1][left+j]-labblur->L[height+imax-i-1][left+j]);
							//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[height+imax-i-1][left+j]-labblur->a[height+imax-i-1][left+j]);
							//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[height+imax-i-1][left+j]-labblur->b[height+imax-i-1][left+j]);
						}
					if (jmax<TS) {
						for (int i=imax; i<TS; i++) 
							for (int j=jmax; j<TS; j++) {
								Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[height+imax-i-1][width+jmax-j-1]-labblur->L[height+imax-i-1][width+jmax-j-1]);
								//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[height+imax-i-1][width+jmax-j-1]-labblur->a[height+imax-i-1][width+jmax-j-1]);
								//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[height+imax-i-1][width+jmax-j-1]-labblur->b[height+imax-i-1][width+jmax-j-1]);
							}
					}
				}
				
				if (jmin>0) {
					for (int j=0; j<jmin; j++) 
						for (int i=imin; i<imax; i++) {
							Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[top+i][2*jmin-j-1]-labblur->L[top+i][2*jmin-j-1]);
							//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[top+i][2*jmin-j-1]-labblur->a[top+i][2*jmin-j-1]);
							//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[top+i][2*jmin-j-1]-labblur->b[top+i][2*jmin-j-1]);
						}
					if (imax<TS) {
						for (int j=0; j<jmin; j++) 
							for (int i=imax; i<TS; i++) {
								Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[height+imax-i-1][2*jmin-j-1]-labblur->L[height+imax-i-1][2*jmin-j-1]);
								//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[height+imax-i-1][2*jmin-j-1]-labblur->a[height+imax-i-1][2*jmin-j-1]);
								//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[height+imax-i-1][2*jmin-j-1]-labblur->b[height+imax-i-1][2*jmin-j-1]);
							}
					}
				}
				
				
				if (jmax<TS) {
					for (int j=jmax; j<TS; j++) 
						for (int i=imin; i<imax; i++) {
							Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[top+i][width+jmax-j-1]-labblur->L[top+i][width+jmax-j-1]);
							//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[top+i][width+jmax-j-1]-labblur->a[top+i][width+jmax-j-1]);
							//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[top+i][width+jmax-j-1]-labblur->b[top+i][width+jmax-j-1]);
						}
					if (imin>0) {
						for (int i=0; i<imin; i++)
							for (int j=jmax; j<TS; j++) {
								Lblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->L[2*imin-i-1][width+jmax-j-1]-labblur->L[2*imin-i-1][width+jmax-j-1]);
								//RLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->a[2*imin-i-1][width+jmax-j-1]-labblur->a[2*imin-i-1][width+jmax-j-1]);
								//BLblox[vblkmod][(indx + i)*TS+j]=tilemask_in[i][j]*(labin->b[2*imin-i-1][width+jmax-j-1]-labblur->b[2*imin-i-1][width+jmax-j-1]);
							}
					}
				}
				//end of tile padding
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				//Lblox[vblkmod][(indx + TS/2)*TS+TS/2]=32768.0f;//testing: locate block centers
				
			}//end of filling block row
			
			//fftwf_print_plan (plan_forward_blox);
			fftwf_execute_dft_r2c(plan_forward_blox,Lblox[vblkmod],fLblox[vblkmod]);	// FT an entire row of tiles
			//fftwf_execute_r2r(plan_forward_blox,Lblox[vblkmod],fLblox[vblkmod]);		// DCT an entire row of tiles
			
			//fftwf_execute_dft_r2c(plan_forward_blox,RLblox[vblkmod],fRLblox[vblkmod]);// FT an entire row of tiles
			//fftwf_execute_dft_r2c(plan_forward_blox,BLblox[vblkmod],fBLblox[vblkmod]);// FT an entire row of tiles
			
			if (vblk<blkrad) continue;
			
			int vblproc = (vblk-blkrad);
			int vblprocmod = vblproc%8;
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// now process the vblproc row of tiles for noise reduction
			for (int hblk=0; hblk<numblox_W; hblk++) {
				
				int hblproc = hblk;
				
				RGBtile_denoise (fLblox, fRLblox, fBLblox, vblproc, hblproc, blkrad, numblox_H, numblox_W, 
								 noisevar_L, noisevar_ab );
				
				if (vblk==(numblox_H-1)) {//denoise last blkrad rows
					for (vblproc=(vblk-blkrad+1); vblproc<numblox_H; vblproc++) {
						vblprocmod = vblproc%8;
						RGBtile_denoise (fLblox, fRLblox, fBLblox, vblproc, hblproc, blkrad, numblox_H, numblox_W, 
										 noisevar_L, noisevar_ab );
					}
				}
				
			}//end of horizontal tile loop
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			//now perform inverse FT of an entire row of tiles
			fftwf_execute_dft_c2r(plan_backward_blox,fLblox[vblprocmod],Lblox[vblprocmod]);	//for FT
			//fftwf_execute_r2r(plan_backward_blox,fLblox[vblprocmod],Lblox[vblprocmod]);	//for DCT
			
			//fftwf_execute_dft_c2r(plan_backward_blox,fRLblox[vblprocmod],RLblox[vblprocmod]);
			//fftwf_execute_dft_c2r(plan_backward_blox,fBLblox[vblprocmod],BLblox[vblprocmod]);
			
			int topproc = (vblproc-blkrad)*offset;
			
			//add row of tiles to output image
			RGBoutput_tile_row (Lblox[vblprocmod], RLblox[vblprocmod], BLblox[vblprocmod], labdn, \
								tilemask_out, height, width, topproc, blkrad );
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			if (vblk==(numblox_H-1)) {//inverse FT last blkrad rows
				for (int vblproc=(vblk-blkrad+1); vblproc<numblox_H; vblproc++) {
					topproc=(vblproc-blkrad)*offset;
					vblprocmod=vblproc%8;
					fftwf_execute_dft_c2r(plan_backward_blox,fLblox[vblprocmod],Lblox[vblprocmod]);	//for FT
					//fftwf_execute_r2r(plan_backward_blox,fLblox[vblprocmod],Lblox[vblprocmod]);	//for DCT
					
					//fftwf_execute_dft_c2r(plan_backward_blox,fRLblox[vblprocmod],RLblox[vblprocmod]);
					//fftwf_execute_dft_c2r(plan_backward_blox,fBLblox[vblprocmod],BLblox[vblprocmod]);
					
					RGBoutput_tile_row (Lblox[vblprocmod], RLblox[vblprocmod], BLblox[vblprocmod], labdn, \
										tilemask_out, height, width, topproc, blkrad );
					
				}
			}
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
		}//end of vertical tile loop
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		// clean up
		//#pragma omp single nowait
		fftwf_destroy_plan( plan_forward_blox );
		//#pragma omp single nowait
		fftwf_destroy_plan( plan_backward_blox );
		
		for( int i = 0 ; i < 8 ; i++ ) {
			fftwf_free ( Lblox[i]);
			//fftwf_free (RLblox[i]);
			//fftwf_free (BLblox[i]);
			fftwf_free ( fLblox[i]);
			//fftwf_free (fRLblox[i]);
			//fftwf_free (fBLblox[i]);
		}
		delete[]  Lblox;	
		//delete[] RLblox;
		//delete[] BLblox;
		delete[]  fLblox;
		//delete[] fRLblox;
		//delete[] fBLblox;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int i=0; i<height; i++) {
			for (int j=0; j<width; j++) {
				//may want to include masking threshold for large hipass data to preserve edges/detail
				float hpdn = labdn->L[i][j]/totwt[i][j];//note that labdn initially stores the denoised hipass data
				
				labdn->L[i][j] = labblur->L[i][j] + hpdn;
				//labdn->L[i][j] = -(hporig)+0.25;
				
				//hpdn = labdn->a[i][j]/totwt[i][j];
				//labdn->a[i][j] = labblur->a[i][j] + hpdn;
				
				//hpdn = labdn->b[i][j]/totwt[i][j];
				//labdn->b[i][j] = labblur->b[i][j] + hpdn;
				
				labdn->a[i][j] = labblur->a[i][j];
				labdn->b[i][j] = labblur->b[i][j];
			}
		}
		
		//dirpyr_ab(labdn, labdn, dnparams);//use dirpyr here if using it to blur ab channels only
		
		//dirpyrdenoise(labdn);//denoise ab channels using ImProcFns denoise (stripped to ab channels only)
		
		//Wavelet denoise of ab channels
		//int numpix = labdn->W*labdn->H;
		//cplx_wavelet_decomposition adecomp(labdn->data+numpix, labdn->W, labdn->H, 5 );//last arg is num levels
		//WaveletDenoise(adecomp, SQR((float)dnparams.chroma*100.0f));
		//adecomp.reconstruct(labdn->data+numpix);
		//cplx_wavelet_decomposition bdecomp(labdn->data+2*numpix, labdn->W, labdn->H, 5 );//last arg is num levels
		//WaveletDenoise(bdecomp, SQR((float)dnparams.chroma*100.0f));
		//bdecomp.reconstruct(labdn->data+2*numpix);
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// transform denoised "Lab" to output RGB
		
		RGB_OutputTransf(labdn, dst, dnparams);
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		delete labin;
		delete labdn;
		delete labblur;
		
	}//end of main RB_denoise
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	void ImProcFunctions::RGBtile_denoise (fftwf_complex ** fLblox, fftwf_complex ** fRLblox, fftwf_complex ** fBLblox, \
										   int vblproc, int hblproc, int blkrad, int numblox_H, int numblox_W, \
										   float noisevar_L, float noisevar_ab )	//for FT
	//void ImProcFunctions::RGBtile_denoise (float ** fLblox, fftwf_complex ** fRLblox, fftwf_complex ** fBLblox, \
	int vblproc, int hblproc, int blkrad, int numblox_H, int numblox_W, \
	float noisevar_L, float noisevar_ab )	//for DCT
	{
		int vblprocmod=vblproc%8;
		
		const float eps = 0.01f/(TS*TS); //tolerance
		const float cutoffsq = 8.0f;//frequency cutoff
		
		
		float RLblockvar=eps, BLblockvar=eps, Lblockvar=eps;
		
		for (int i=TS/4; i<3*TS/4; i++) 
			for (int j=TS/4; j<fTS; j++) {	//for FT
				//for (int j=TS/4; j<3*TS/4; j++) {	//for DCT
				
				Lblockvar += (SQR( fLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0])+SQR( fLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1]));	//for FT
				//Lblockvar += SQR( fLblox[vblprocmod][(hblproc*TS+i)*TS+j]);		//for DCT
				
				//RLblockvar += (SQR(fRLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0])+SQR(fRLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1]));
				//BLblockvar += (SQR(fBLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0])+SQR(fBLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1]));
				
			}
		
		Lblockvar /= (TS/2)*(fTS-TS/4);	//for FT
		//Lblockvar /= (TS/2)*(TS/2);	//for DCT
		
		//RLblockvar /= (TS/2)*(fTS-TS/4);
		//BLblockvar /= (TS/2)*(fTS-TS/4);
		Lblockvar = (3*Lblockvar);
		//RLblockvar = (3*RLblockvar);
		//BLblockvar = (3*BLblockvar);
		
		//printf("vblock=%d  hblock=%d  blockstddev=%f \n",vblproc,hblproc,sqrt(blockvar));
		
		//float wsqave=0;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for (int i=0; i<TS; i++) 
			for (int j=0; j<fTS; j++) {		//for FT
				//for (int j=0; j<TS; j++) {	//for DCT
				
				
				float Lcoeffre = fLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0];	//for FT
				float Lcoeffim = fLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1];
				
				//float Lcoeff = fLblox[vblprocmod][(hblproc*TS+i)*TS+j];		//for DCT
				
				//float RLcoeffre = fRLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0];
				//float RLcoeffim = fRLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1];
				
				//float BLcoeffre = fBLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0];
				//float BLcoeffim = fBLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1];
				
				/*double nbrave=0, nbrsqave=0, coeffsq;
				int vblnbrmod, hblnbrmod;
				for (int k=0; k<9; k++) {
					vblnbrmod = (vblproc+(k/3)+7)%8;
					hblnbrmod = MAX(0,hblproc+(k%3)-1);
					if (hblnbrmod==numblox_W) hblnbrmod=numblox_W-1;
					coeffsq = SQR(fLblox[vblnbrmod][(hblnbrmod*TS+i)*fTS+j][0])+SQR(fLblox[vblnbrmod][(hblnbrmod*TS+i)*fTS+j][1]);	//for FT
					//coeffsq = SQR(fLblox[vblnbrmod][(hblnbrmod*TS+i)*TS+j]);	//for DCT
					
					nbrave += sqrt(coeffsq);
					nbrsqave += coeffsq;
				}
				float nbrvar = (nbrsqave/9.0f)-SQR(nbrave/9.0f);*/
				
				float  Lwsq = eps+SQR( Lcoeffre)+SQR( Lcoeffim);	//for FT
				//float  Lwsq = eps+SQR( Lcoeff);					//for DCT
				
				//float RLwsq = eps+SQR(RLcoeffre)+SQR(RLcoeffim);
				//float BLwsq = eps+SQR(BLcoeffre)+SQR(BLcoeffim);
				
				//wsqave += Lwsq;
				//float Lfactor = (4*Lblockvar)/(eps+(Lwsq+nbrvar)+2*Lblockvar);
				//float Lfactor = expf(-Lwsq/(9*Lblockvar));
				float freqfactor = 1.0f-MAX((expf(-(SQR(i)+SQR(j))/cutoffsq)),(expf(-(SQR(TS-i)+SQR(j))/cutoffsq)));
				float  Lfactor = 1;//freqfactor;//*(2* Lblockvar)/(eps+ Lwsq+ Lblockvar);
				//float RLfactor = 1;//(2*RLblockvar)/(eps+RLwsq+RLblockvar);
				//float BLfactor = 1;//(2*BLblockvar)/(eps+BLwsq+BLblockvar);

				float  Lshrinkfactor = SQR(Lwsq/(Lwsq + noisevar_L * Lfactor*exp(-Lwsq/(3*noisevar_L))));
				//float RLshrinkfactor = RLwsq/(RLwsq+noisevar_ab*RLfactor*exp(-Lwsq/(3*noisevar_L)));
				//float BLshrinkfactor = BLwsq/(BLwsq+noisevar_ab*BLfactor*exp(-Lwsq/(3*noisevar_L)));
				
				//float shrinkfactor = (wsq<noisevar ? 0 : 1);//hard threshold
				
				fLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0] *= Lshrinkfactor;	//for FT
				fLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1] *= Lshrinkfactor;
				//fLblox[vblprocmod][(hblproc*TS+i)*TS+j] *= Lshrinkfactor;		//for DCT
				
				//fRLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0] *= RLshrinkfactor;
				//fRLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1] *= RLshrinkfactor;
				
				//fBLblox[vblprocmod][(hblproc*TS+i)*fTS+j][0] *= BLshrinkfactor;
				//fBLblox[vblprocmod][(hblproc*TS+i)*fTS+j][1] *= BLshrinkfactor;
				
			}//end of block denoise		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//printf("vblk=%d  hlk=%d  wsqave=%f   ||   ",vblproc,hblproc,wsqave);
		
	}//end of function tile_denoise
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGBoutput_tile_row (float *bloxrow_L, float *bloxrow_a, float *bloxrow_b, LabImage * labdn,
											  float ** tilemask_out, int height, int width, int top, int blkrad )
	{
		const int numblox_W = ceil(((float)(width))/(offset))+2*blkrad;
		
		//add row of tiles to output image
		for (int hblk=0; hblk < numblox_W; hblk++) {
			int left = (hblk-blkrad)*offset;
			int bottom = MIN( top+TS,height);
			int right  = MIN(left+TS, width);
			int imin = MAX(0,-top);
			int jmin = MAX(0,-left);
			int imax = bottom - top;
			int jmax = right - left;
			
			int indx = hblk*TS;
			
			for (int i=imin; i<imax; i++)
				for (int j=jmin; j<jmax; j++) {
					
					labdn->L[top+i][left+j] += tilemask_out[i][j]*bloxrow_L[(indx + i)*TS+j]/(TS*TS); 
					//labdn->a[top+i][left+j] += tilemask_out[i][j]*bloxrow_a[(indx + i)*TS+j]/(TS*TS); 
					//labdn->b[top+i][left+j] += tilemask_out[i][j]*bloxrow_b[(indx + i)*TS+j]/(TS*TS); 
					
				}
		}
		
	}
	
#undef TS
#undef fTS
#undef offset
	//#undef eps
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//experimental dirpyr low-pass filter

void ImProcFunctions::dirpyr_ab(LabImage * data_fine, LabImage * data_coarse, const procparams::DirPyrDenoiseParams & dnparams)
{
	int W = data_fine->W;
	int H = data_fine->H;
	float thresh_L = 10.0f*dnparams.luma;
	float threshsq_L = SQR(thresh_L);
	float thresh_ab = 10.0f*dnparams.chroma;
	float threshsq_ab = SQR(thresh_ab);
	LUTf rangefn_L(0x10000);
	LUTf rangefn_ab(0x10000);

	LabImage * dirpyrlo[2];
	
	//set up range functions
	
	for (int i=0; i<0x10000; i++) {
		rangefn_L[i] =  exp(-((float)i) / (1.0+thresh_L)) ;// * (1.0+thresh_L)/(((float)i) + thresh_L+1.0); 
		rangefn_ab[i] =  exp(-SQR((float)i) / (1.0+threshsq_ab)) ;// * (1.0+thresh_ab)/(((float)i) + thresh_ab+1.0); 
	}
	dirpyrlo[0] = new LabImage (W, H);
	dirpyrlo[1] = new LabImage (W, H);
	
	//int scale[4]={1,3,5,9/*1*/};
	int scale[5]={1,2,4,7,13/*1*/};

	int level=0;
	int indx=0;
	dirpyr_ablevel(data_fine, dirpyrlo[indx], W, H, rangefn_L,rangefn_ab, 0, scale[level] );
	level += 1;
	indx = 1-indx;
	while (level<3) {
		dirpyr_ablevel(dirpyrlo[1-indx], dirpyrlo[indx], W, H, rangefn_L,rangefn_ab, level, scale[level] );
		level += 1;
		indx = 1-indx;
	}
	
	dirpyr_ablevel(dirpyrlo[1-indx], data_coarse, W, H, rangefn_L,rangefn_ab, level, scale[level] );
	
	//delete dirpyrlo[0];//TODO: this seems to disable the NR ???
	//delete dirpyrlo[1];
}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

void ImProcFunctions::dirpyr_ablevel(LabImage * data_fine, LabImage * data_coarse, int width, int height, LUTf & rangefn_L, LUTf & rangefn_ab, int level, int scale)
{
	//scale is spacing of directional averaging weights
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// calculate weights, compute directionally weighted average
	
	//int domker[5][5] = {{1,1,1,1,1},{1,2,2,2,1},{1,2,2,2,1},{1,2,2,2,1},{1,1,1,1,1}};
	//int domker[5][5] = {{1, 2, 4, 2, 1}, {2, 4, 8, 4, 2}, {4, 8, 16, 8, 4}, {2, 4, 8, 4, 2}, {1, 2, 4, 2, 1}};
	float domker[5][5] = {{0.129923f, 0.279288f, 0.360448f, 0.279288f, 0.129923f}, \
		{0.279288f, 0.600373f, 0.774837f, 0.600373f, 0.279288f}, \
		{0.360448f, 0.774837f, 1.0f,      0.774837f, 0.360448f}, \
		{0.279288f, 0.600373f, 0.774837f, 0.600373f, 0.279288f}, \
		{0.129923f, 0.279288f, 0.360448f, 0.279288f, 0.129923f}};//Gaussian with sigma=1.4
	
	
	int scalewin = 2*scale;
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i = 0; i < height; i++) {
		for(int j = 0; j < width; j++)
		{
			float valL=0, vala=0, valb=0;
			float norm_L=0, norm_ab=0;
			
			for(int inbr=MAX(0,(i-scalewin)); inbr<=MIN(height-1,(i+scalewin)); inbr+=scale) {
				for (int jnbr=MAX(0,(j-scalewin)); jnbr<=MIN(width-1,(j+scalewin)); jnbr+=scale) {
					//it seems that weighting the blur by L (gamma=3) works better
					//than using the variable gamma source
					//float desat = 1-rangefn_ab[data_fine->L[i][j]+abs(data_fine->a[i][j])+abs(data_fine->b[i][j])];

					float nbrdiff_L = fabs(data_fine->L[inbr][jnbr]-data_fine->L[i][j])/level;
					float nbrdiff_ab = (fabs(data_fine->a[inbr][jnbr]-data_fine->a[i][j]) + \
										fabs(data_fine->b[inbr][jnbr]-data_fine->b[i][j]));
					float dirwt_L = ( domker[(inbr-i)/scale+2][(jnbr-j)/scale+2] * rangefn_L[nbrdiff_L] );
					float dirwt_ab = ( /*domker[(inbr-i)/scale+2][(jnbr-j)/scale+2] */ rangefn_ab[nbrdiff_ab] );
					//valL += dirwt_L *data_fine->L[inbr][jnbr];
					vala += dirwt_L*dirwt_ab*data_fine->a[inbr][jnbr];
					valb += dirwt_L*dirwt_ab*data_fine->b[inbr][jnbr];
					//norm_L += dirwt_L;
					norm_ab += dirwt_L*dirwt_ab;

				}
			}
			
			//data_coarse->L[i][j] = valL/norm_L; // low pass filter
			data_coarse->L[i][j] = data_fine->L[i][j];
			data_coarse->a[i][j] = vala/norm_ab; // low pass filter
			data_coarse->b[i][j] = valb/norm_ab; // low pass filter
			
			/*if (level!=3) {
				data_coarse->L[i][j] = valL/norm_L; // low pass filter
			} else {
				float valL=0, vala=0, valb=0;
				float norm=0;
				for(int inbr=MAX(0,(i-2)); inbr<=MIN(height-1,(i+2)); inbr++) {
					for (int jnbr=MAX(0,(j-2)); jnbr<=MIN(width-1,(j+2)); jnbr++) {
						//it seems that weighting the blur by Lab luminance (~gamma=3) 
						//works better than using the variable gamma source
						float nbrdiff = (fabs(data_fine->L[inbr][jnbr]-data_fine->L[i][j]) + \
										 fabs(data_fine->a[inbr][jnbr]-data_fine->a[i][j]) + \
										 fabs(data_fine->b[inbr][jnbr]-data_fine->b[i][j]))/(level);
						float dirwt = ( domker[(inbr-i)/scale+2][(jnbr-j)+2] * rangefn_L[nbrdiff] );
						valL += dirwt*data_fine->L[inbr][jnbr];
						vala += dirwt*data_fine->a[inbr][jnbr];
						valb += dirwt*data_fine->b[inbr][jnbr];
						norm += dirwt;
						
					}
				}
				data_coarse->L[i][j] = data_fine->L[i][j];//valL/norm;
				data_coarse->a[i][j] = vala/norm; 
				data_coarse->b[i][j] = valb/norm; 
			}*/

		}
	}
	
}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
/*
void ImProcFunctions::FixImpulse_ab(LabImage * src, LabImage * dst, double radius, int thresh) { 
 
	
	float threshsqr = SQR(thresh);
	int halfwin = ceil(2*radius)+1;
		
	// local variables
	int width=src->W, height=src->H;
	//temporary array to store chromaticity
	float *fringe = float * calloc ((height)*(width), sizeof *fringe);
	
	LabImage * tmp1; 
	tmp1 = new LabImage(width, height);
	
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(src->W,src->H));
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
			if (33*fringe[i*width+j]>thresh*chromave) {
				float atot=0;
				float btot=0;
				float norm=0;
				float wt;
				for (int i1=MAX(0,i-halfwin+1); i1<MIN(height,i+halfwin); i1++) 
					for (int j1=MAX(0,j-halfwin+1); j1<MIN(width,j+halfwin); j1++) {
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
*/
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::WaveletDenoise(cplx_wavelet_decomposition &DualTreeCoeffs, float noisevar ) 
	{
		int maxlvl = DualTreeCoeffs.maxlevel();
		int rad_stage1[8] = {3,2,1,1,1,1,1,1};
		int rad_stage2[8] = {2,1,1,1,1,1,1,1};
		
		for (int lvl=0; lvl<maxlvl-1; lvl++) {
			int Wlvl = DualTreeCoeffs.level_W(lvl,0);
			int Hlvl = DualTreeCoeffs.level_H(lvl,0);
			
			array2D<float> wiener1(Wlvl,Hlvl);

			for (int m=0; m<2; m++) {
				float ** ReCoeffs = DualTreeCoeffs.level_coeffs(lvl,0+m);
				float ** ImCoeffs = DualTreeCoeffs.level_coeffs(lvl,2+m);
				float ** ReParents = DualTreeCoeffs.level_coeffs(lvl+1,0+m);
				float ** ImParents = DualTreeCoeffs.level_coeffs(lvl+1,2+m);
				int ParentPadding = DualTreeCoeffs.level_pad(lvl+1,0+m);
				for (int dir=1; dir<4; dir++) {
					//FirstStageWiener (ReCoeffs[dir],ImCoeffs[dir],wiener1,Wlvl,Hlvl,rad_stage1[lvl], noisevar);
					//SecondStageWiener(ReCoeffs[dir],ImCoeffs[dir],wiener1,Wlvl,Hlvl,rad_stage2[lvl], noisevar);

					BiShrink(ReCoeffs[dir], ImCoeffs[dir], ReParents[dir], ImParents[dir], Wlvl, Hlvl, lvl, ParentPadding, noisevar);
				}
			}
		}
		
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::WaveletDenoise(wavelet_decomposition &WaveletCoeffs, float noisevar ) 
	{
		int maxlvl = WaveletCoeffs.maxlevel();
		int rad_stage1[8] = {3,2,1,1,1,1,1,1};
		int rad_stage2[8] = {2,1,1,1,1,1,1,1};
		
		for (int lvl=0; lvl<maxlvl/*-1*/; lvl++) {
			int Wlvl = WaveletCoeffs.level_W(lvl);
			int Hlvl = WaveletCoeffs.level_H(lvl);
			
			//array2D<float> wiener1(Wlvl,Hlvl);
			int ParentPadding;
			float ** WavParents;
			float ** WavCoeffs = WaveletCoeffs.level_coeffs(lvl);
			if (lvl<maxlvl-1) {
				WavParents = WaveletCoeffs.level_coeffs(lvl+1);
				ParentPadding = WaveletCoeffs.level_pad(lvl+1);
			} else {
				WavParents = WaveletCoeffs.level_coeffs(lvl);
				ParentPadding = 0;
			}
			
			float threshsq = noisevar;//*SQR(UniversalThresh(WaveletCoeffs.level_coeffs(lvl)[3], Wlvl*Hlvl));
				
			Shrink(WavCoeffs, Wlvl, Hlvl, lvl, threshsq/*noisevar*/);
		}
		
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::BiShrink(float * ReCoeffs, float * ImCoeffs, float * ReParents, float * ImParents, int W, int H, int level, int pad, float noisevar) 
	{	
		//bivariate shrinkage of Sendur & Selesnick
		float * sigma = new float[W*H]; 
		int rad = 3*(1<<level);
		boxvar(ReCoeffs,sigma,rad,rad,W,H);//box blur detail coeffs to estimate local variance
		const float root3 = sqrt(3);
		const int Wpar = (W+2*pad);
		//const int Hpar = (H+1+2*pad)/2;
		const float eps = 0.01f;
		
		for (int i=0; i<H; i++) {
			for (int j=0; j<W; j++) {
				
				float thresh = root3 * noisevar/sqrt(MAX(sigma[i*W+j]-noisevar, eps));
				int parentloc = ((i)-pad)*Wpar+(j)-pad;
				int coeffloc  = i*W+j;
				float mag = sqrt(SQR(ReCoeffs[coeffloc]) + SQR(ImCoeffs[coeffloc]) + SQR(ReParents[parentloc]) + SQR(ImParents[parentloc]));
				float shrinkfactor = MAX(0,mag-thresh);
				shrinkfactor /= (shrinkfactor+thresh+eps);
				//float shrinkfactor = mag/(mag+noisevar+eps);
				//float shrinkre = SQR(ReCoeffs[coeffloc])/(noisevar+ SQR(ReCoeffs[coeffloc]) +eps);
				//float shrinkim = SQR(ImCoeffs[coeffloc])/(noisevar+ SQR(ImCoeffs[coeffloc]) +eps);
				ReCoeffs[coeffloc] *= shrinkfactor;
				ImCoeffs[coeffloc] *= shrinkfactor;
				
			}
		}
		
		delete[] sigma;
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::Shrink(float ** WavCoeffs, int W, int H, int level, float noisevar) 
	{	
		//simple wavelet shrinkage
		float * sigma = new float[W*H]; 
		const float eps = 0.01f;
		printf("level=%d  ",level);
		for (int dir=1; dir<4; dir++) {
			float mad = SQR(UniversalThresh(WavCoeffs[dir], W*H));//*6*/(level+1);
			printf("  dir=%d  mad=%f	",dir,sqrt(mad));
			for (int i=0; i<H; i++) {
				for (int j=0; j<W; j++) {
					
					int coeffloc  = i*W+j;
					float mag = SQR(WavCoeffs[dir][coeffloc]);
					float shrinkfactor = mag/(mag+noisevar*mad*exp(-mag/(3*noisevar*mad))+eps);
					//float shrinkfactor = mag/(mag+noisevar*SQR(sigma[coeffloc])+eps);
					
					//WavCoeffs[dir][coeffloc] *= shrinkfactor;
					sigma[coeffloc] = shrinkfactor;
				}
			}
			
			boxblur(sigma, sigma, 1, 1, W, H);//increase smoothness by locally averaging shrinkage
			for (int i=0; i<W*H; i++) {
				float mag = SQR(WavCoeffs[dir][i]);
				float sf = mag/(mag+noisevar*mad+eps);
				
				//use smoothed shrinkage unless local shrinkage is much less
				WavCoeffs[dir][i] *= (SQR(sigma[i])+SQR(sf))/(sigma[i]+sf+eps);
				
				//the following is for testing
				//float wdn = WavCoeffs[dir][i]*sf;
				//float sf1 = mag/(mag+4*noisevar*mad+eps);
				//float wdn1 = WavCoeffs[dir][i]*sf1;

				//WavCoeffs[dir][i] = wdn-wdn1;
			}
		}
		printf("\n");
		delete[] sigma;
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	float ImProcFunctions::UniversalThresh(float * Detail_Coeffs, int datalen) {
		
		int * histo = new int[65536];
		//memset(histo, 0, 65536*sizeof(histo));
		for (int i=0; i<65536; i++) histo[i]=0;
		
		//calculate histogram of absolute values of HH wavelet coeffs
		for (int i=0; i<datalen; i++) {
			histo[MAX(0,MIN(65535,abs((int)Detail_Coeffs[i])))]++;
		}
		
		//find median of histogram
		int median=0, count=0;
		while (count<datalen/2) {
			count += histo[median];
			median++;
		}
		
		int count_ = count - histo[median-1];
		
		delete[] histo;
		
		// interpolate 
		return (( (median-1) + (datalen/2-count_)/((float)(count-count_)) )/0.6745);
		
	}

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::WaveletDenoiseAll(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a, 
											wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_ab ) 
	{
		int maxlvl = WaveletCoeffs_L.maxlevel();
		
		for (int lvl=0; lvl<maxlvl; lvl++) {
			
			int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
			int Hlvl_L = WaveletCoeffs_L.level_H(lvl);
			
			int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
			int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);
			
			float skip_L = WaveletCoeffs_L.level_stride(lvl);
			float skip_ab = WaveletCoeffs_a.level_stride(lvl);

			float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);
						
			ShrinkAll(WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, lvl, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,
					  skip_L, skip_ab, noisevar_L, noisevar_ab);
		}
		
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	void ImProcFunctions::ShrinkAll(float ** WavCoeffs_L, float ** WavCoeffs_a, float ** WavCoeffs_b, int level, 
								 int W_L, int H_L, int W_ab, int H_ab, int skip_L, int skip_ab, float noisevar_L, float noisevar_ab) 
	{	
		//simple wavelet shrinkage
		const float eps = 0.01f;
		float * sigma = new float[W_L*H_L]; 

		printf("\n level=%d  \n",level);
		
		for (int dir=1; dir<4; dir++) {
			float mad_L = SQR(UniversalThresh(WavCoeffs_L[dir], W_L*H_L));//*6*/(level+1);
			float mad_a = SQR(UniversalThresh(WavCoeffs_a[dir], W_ab*H_ab));//*6*/(level+1);
			float mad_b = SQR(UniversalThresh(WavCoeffs_b[dir], W_ab*H_ab));//*6*/(level+1);
			
			float thresh_L = sqrt(mad_L*noisevar_L);
			float thresh_a = sqrt(mad_a*noisevar_ab);
			float thresh_b = sqrt(mad_b*noisevar_ab);

			printf("  dir=%d  mad_L=%f	mad_a=%f	mad_b=%f	\n",dir,sqrt(mad_L),sqrt(mad_a),sqrt(mad_b));
			
			for (int i=0; i<H_ab; i++) {
				for (int j=0; j<W_ab; j++) {
					
					int coeffloc_ab = i*W_ab+j;
					int coeffloc_L	= ((i*skip_L)/skip_ab)*W_L + ((j*skip_L)/skip_ab);
					
					float mag_L = fabs(WavCoeffs_L[dir][coeffloc_L ])+eps;
					float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
					float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;
					
					float edgefactor = exp(-mag_L/(sqrt(noisevar_L*mad_L))) * exp(-mag_a/(3*noisevar_ab*mad_a)) * exp(-mag_b/(3*noisevar_ab*mad_b));

					//WavCoeffs_a[dir][coeffloc_ab] *= mag_a/(mag_a + noisevar_ab*mad_a*edgefactor + eps);
					//WavCoeffs_b[dir][coeffloc_ab] *= mag_b/(mag_b + noisevar_ab*mad_b*edgefactor + eps);
					
					float coeff_a = fabs(WavCoeffs_a[dir][coeffloc_ab]);
					float coeff_b = fabs(WavCoeffs_b[dir][coeffloc_ab]);

					// 'firm' threshold of chroma coefficients
					WavCoeffs_a[dir][coeffloc_ab] *= (coeff_a>2*thresh_a ? 1 : (coeff_a<thresh_a ? 0 : (coeff_a/thresh_a - 1)));
					WavCoeffs_b[dir][coeffloc_ab] *= (coeff_b>2*thresh_b ? 1 : (coeff_b<thresh_b ? 0 : (coeff_b/thresh_b - 1)));

					//WavCoeffs_b[dir][coeffloc_ab] *= (fabs(WavCoeffs_b[dir][coeffloc_ab])<thresh_a*noise_ab ? 0 : 1);


				}
			}
			
			for (int i=0; i<W_L*H_L; i++) {
				
				//float coeff_L = fabs(WavCoeffs_L[dir][i]);
				// 'firm' threshold of luma coefficients
				//float shrinkfactor = (coeff_L>2*thresh_L ? 1 : (coeff_L<thresh_L ? 0 : (coeff_L/thresh_L - 1)));
				
				float mag = SQR(WavCoeffs_L[dir][i]);
				float shrinkfactor = mag/(mag+noisevar_L*mad_L*exp(-mag/(3*noisevar_L*mad_L))+eps);
				
				//float shrinkfactor = mag/(mag+noisevar*SQR(sigma[coeffloc])+eps);
				
				//WavCoeffs_L[dir][i] *= shrinkfactor;
				sigma[i] = shrinkfactor;
			}
			
			boxblur(sigma, sigma, level+2, level+2, W_L, H_L);//increase smoothness by locally averaging shrinkage
			for (int i=0; i<W_L*H_L; i++) {
				
				//float coeff_L = fabs(WavCoeffs_L[dir][i]);
				// 'firm' threshold of chroma coefficients
				//float sf = (coeff_L>2*thresh_L ? 1 : (coeff_L<thresh_L ? 0 : (coeff_L/thresh_L - 1)));
				
				float mag = SQR(WavCoeffs_L[dir][i]);
				float sf = mag/(mag+noisevar_L*mad_L+eps);
				
				//use smoothed shrinkage unless local shrinkage is much less
				WavCoeffs_L[dir][i] *= (SQR(sigma[i])+SQR(sf))/(sigma[i]+sf+eps);
				
				//the following is for testing
				//float wdn = WavCoeffs[dir][i]*sf;
				//float sf1 = mag/(mag+4*noisevar*mad+eps);
				//float wdn1 = WavCoeffs[dir][i]*sf1;
				
				//WavCoeffs[dir][i] = wdn-wdn1;
			}//now luminance coeffs are denoised
			
			

		}
		delete[] sigma;
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
};
