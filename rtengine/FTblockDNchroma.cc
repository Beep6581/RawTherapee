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
#include "gauss.h"

#include "rtengine.h"
#include "improcfun.h"
#include "LUT.h"
#include "array2D.h"
#include "iccmatrices.h"
#include "boxblur.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cplx_wavelet_dec.h"


#define SQR(x) ((x)*(x))
//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
//#define LIM(x,min,max) MAX(min,MIN(x,max))
//#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
//#define CLIP(x) LIM(x,0,65535)

#define TS 64		// Tile size
#define offset 25	// shift between tiles
#define fTS ((TS/2+1))	// second dimension of Fourier tiles
#define blkrad 0	// radius of block averaging
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
		
		//int max;
		//float median = MadMax(src->data, max, src->width*src->height);
		//gain = sqrt(MAX(1.0f,(0.15f*65535.0f/median))*(65535.0f/max));//'gain' is public float allocated in improcfun.h;
		const float gain = pow (2.0, dnparams.expcomp);
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<src->height; i++) 
			for (int j=0; j<src->width; j++) {
				
				float X = gain*src->r[i][j];//xyz_prophoto[0][0]*src->r[i][j] + xyz_prophoto[0][1]*src->g[i][j] + xyz_prophoto[0][2]*src->b[i][j];
				float Y = gain*src->g[i][j];//xyz_prophoto[1][0]*src->r[i][j] + xyz_prophoto[1][1]*src->g[i][j] + xyz_prophoto[1][2]*src->b[i][j];
				float Z = gain*src->b[i][j];//xyz_prophoto[2][0]*src->r[i][j] + xyz_prophoto[2][1]*src->g[i][j] + xyz_prophoto[2][2]*src->b[i][j];
				
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
		
		const float gain = pow (2.0, dnparams.expcomp);
		
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
				
				dst->r[i][j] = X/gain;//prophoto_xyz[0][0]*X + prophoto_xyz[0][1]*Y + prophoto_xyz[0][2]*Z;
				dst->g[i][j] = Y/gain;//prophoto_xyz[1][0]*X + prophoto_xyz[1][1]*Y + prophoto_xyz[1][2]*Z;
				dst->b[i][j] = Z/gain;//prophoto_xyz[2][0]*X + prophoto_xyz[2][1]*Y + prophoto_xyz[2][2]*Z;
				
			}
		}
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGB_denoise(Imagefloat * src, Imagefloat * dst, /*int Roffset,*/ const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe)
	{	
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		/*if (plistener) {
		 plistener->setProgressStr ("Denoise...");
		 plistener->setProgress (0.0);
		 }*/
		
		volatile double progress = 0.0;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		const short int height=src->height, width=src->width;
		const short int hfh=(height+1)/2, hfw=(width+1)/2;
		
		if (dnparams.Lamt==0) {//nothing to do; copy src to dst
			for (int i=0; i<height; i++) {
				for (int j=0; j<width; j++) {
					dst->r[i][j] = src->r[i][j];
					dst->r[i][j] = src->r[i][j];
					dst->r[i][j] = src->r[i][j];
				}
			}
			return;
		}
		
		//const int blkrad=2;
		float noisevar_L = SQR((100-dnparams.luma) * TS * 100.0f);
		float noisevar_ab = SQR(dnparams.chroma * TS * 150.0f);
		
		// calculation for tiling
		const int numblox_W = ceil(((float)(width))/(offset))+2*blkrad;
		const int numblox_H = ceil(((float)(height))/(offset))+2*blkrad;
		//const int nrtiles = numblox_W*numblox_H;
		// end of tiling calc
				
		
		
		array2D<float> tilemask_in(TS,TS);
		array2D<float> tilemask_out(TS,TS);
		
		array2D<float> totwt(width,height,ARRAY2D_CLEAR_DATA);
				
		const int border = MAX(2,TS/16);
		
		for (int i=0; i<TS; i++) {
			float i1 = abs((i>TS/2 ? i-TS+1 : i));
			float vmask = (i1<border ? SQR(sin((M_PI*i1)/(2*border))) : 1.0f);
			float vmask2 = (i1<2*border ? SQR(sin((M_PI*i1)/(2*border))) : 1.0f);
			for (int j=0; j<TS; j++) {
				float j1 = abs((j>TS/2 ? j-TS+1 : j));
				tilemask_in[i][j] = (vmask * (j1<border ? SQR(sin((M_PI*j1)/(2*border))) : 1.0f));
				
				tilemask_out[i][j] = vmask2 * (j1<2*border ? SQR(sin((M_PI*j1)/(2*border))) : 1.0f);
				
			}
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		LabImage * labin = new LabImage(width,height);
		LabImage * labdn = new LabImage(width,height);
		array2D<float> Ldn(width,height);
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// transform RGB input to ersatz Lab
		
		RGB_InputTransf(src, labin, dnparams, defringe);
		
		memcpy (labdn->data, labin->data, 3*width*height*sizeof(float));

		impulse_nr (labdn, 50.0f/20.0f);
		
		int datalen = labin->W * labin->H;
		
		wavelet_decomposition Ldecomp(labin->data, labin->W, labin->H, 5/*maxlevels*/, 0/*subsampling*/ );
		wavelet_decomposition adecomp(labin->data+datalen, labin->W, labin->H, 5, 1 );//last args are maxlevels, subsampling
		wavelet_decomposition bdecomp(labin->data+2*datalen, labin->W, labin->H, 5, 1 );//last args are maxlevels, subsampling

		float noisevarL	 = SQR(dnparams.Lamt/25.0f);//TODO: clean up naming confusion about params
		float noisevarab = SQR(dnparams.chroma/10.0f);
		
		WaveletDenoiseAll_BiShrink(Ldecomp, adecomp, bdecomp, noisevarL, noisevarab);

		Ldecomp.reconstruct(labdn->data);
		adecomp.reconstruct(labdn->data+datalen);
		bdecomp.reconstruct(labdn->data+2*datalen);
		
		impulse_nr (labdn, 50.0f/20.0f);
		//PF_correct_RT(dst, dst, defringe.radius, defringe.threshold);

		
		float * Ldnptr = Ldn;
		memcpy (Ldnptr, labdn->data, width*height*sizeof(float));
		for (int i=0; i<labdn->W*labdn->H; i++) {
			labdn->data[i] = 0.0f;
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// allocate FT data structures
		
		float ** Lblox = new float *[8] ;
		
		float ** fLblox = new float *[8] ;				//for DCT
		
		for( int i = 0 ; i < 8 ; i++ ) {
			Lblox[i] = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));
			fLblox[i] = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));				//for DCT
			
		}
		
		//make a plan for FFTW
		fftwf_plan plan_forward_blox, plan_backward_blox;
		
		int nfwd[2]={TS,TS};
		
		
		//for DCT:
		const fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
		const fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};
		
		plan_forward_blox  = fftwf_plan_many_r2r(2, nfwd, numblox_W, Lblox[0], NULL, 1, TS*TS, fLblox[0], NULL, 1, TS*TS, fwdkind, FFTW_ESTIMATE );
		plan_backward_blox = fftwf_plan_many_r2r(2, nfwd, numblox_W, fLblox[0], NULL, 1, TS*TS, Lblox[0], NULL, 1, TS*TS, bwdkind, FFTW_ESTIMATE );
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Main algorithm: Tile loop
#pragma omp parallel for schedule(dynamic)
		//int vblock=0, hblock=0;
		for (int vblk=0; vblk<numblox_H; vblk++) {
			//printf("vblock=%d",vblk);
			int vblkmod = vblk%8;
			
			int top = (vblk-blkrad)*offset;

			float * buffer = new float [width + TS + 2*blkrad*offset];
			float * datarow = buffer+blkrad*offset;
			for (int i=0, row=top; i<TS; i++, row++) {
				
				int rr = row;
				if (row<0) {
					rr = MIN(-row,height-1);
				} else if (row>=height) {
					rr = MAX(0,2*height-2-row);
				} 
				
				for (int j=0; j<labin->W; j++) {
					datarow[j] = (labin->L[rr][j]-Ldn[rr][j]);
				}
				
				for (int j=-blkrad*offset; j<0; j++) {
					datarow[j] = datarow[MIN(-j,width-1)];
				}
				for (int j=width; j<width+(TS-(width%TS)-1)+blkrad*offset; j++) {
					datarow[j] = datarow[MAX(0,2*width-2-j)];
				}//now we have a padded data row
				
				//now fill this row of the tiles with Lab high pass data
				for (int hblk=0; hblk<numblox_W; hblk++) {
					int left = (hblk-blkrad)*offset;
					int indx = (hblk)*TS;//index of block in malloc
					
					for (int j=0; j<TS; j++) {
						
						Lblox[vblkmod][(indx + i)*TS+j]  = tilemask_in[i][j]*datarow[left+j];// luma data
						if (top+i>=0 && top+i<height && left+j>=0 && left+j<width) {
							totwt[top+i][left+j] += tilemask_in[i][j]*tilemask_out[i][j];
						}
					}
				}
			}//end of filling block row
			delete[] buffer;

			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			
			//fftwf_print_plan (plan_forward_blox);
			fftwf_execute_r2r(plan_forward_blox,Lblox[vblkmod],fLblox[vblkmod]);		// DCT an entire row of tiles

			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			if (vblk<blkrad) continue;
			
			int vblproc = (vblk-blkrad);
			int vblprocmod = vblproc%8;
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// now process the vblproc row of tiles for noise reduction
			for (int hblk=0; hblk<numblox_W; hblk++) {
				
				int hblproc = hblk;
				
				RGBtile_denoise (fLblox, vblproc, hblproc, numblox_H, numblox_W, noisevar_L, noisevar_ab );
				
				if (vblk==(numblox_H-1)) {//denoise last blkrad rows
					for (vblproc=(vblk-blkrad+1); vblproc<numblox_H; vblproc++) {
						vblprocmod = vblproc%8;
						RGBtile_denoise (fLblox, vblproc, hblproc, numblox_H, numblox_W, noisevar_L, noisevar_ab );
					}
				}
				
			}//end of horizontal tile loop
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			//now perform inverse FT of an entire row of tiles
			fftwf_execute_r2r(plan_backward_blox,fLblox[vblprocmod],Lblox[vblprocmod]);	//for DCT
			
			int topproc = (vblproc-blkrad)*offset;
			
			//add row of tiles to output image
			RGBoutput_tile_row (Lblox[vblprocmod], labdn, tilemask_out, height, width, topproc );
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			if (vblk==(numblox_H-1)) {//inverse FT last blkrad rows
				for (int vblproc=(vblk-blkrad+1); vblproc<numblox_H; vblproc++) {
					topproc=(vblproc-blkrad)*offset;
					vblprocmod=vblproc%8;
					fftwf_execute_r2r(plan_backward_blox,fLblox[vblprocmod],Lblox[vblprocmod]);	//for DCT
					
					RGBoutput_tile_row (Lblox[vblprocmod], labdn, tilemask_out, height, width, topproc );
					
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
			fftwf_free ( fLblox[i]);
		}
		delete[]  Lblox;
		delete[]  fLblox;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int i=0; i<height; i++) {
			for (int j=0; j<width; j++) {
				//may want to include masking threshold for large hipass data to preserve edges/detail
				float hpdn = labdn->L[i][j]/totwt[i][j];//note that labdn initially stores the denoised hipass data
				
				labdn->L[i][j] = Ldn[i][j] + hpdn;
				
			}
		}
				
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// transform denoised "Lab" to output RGB
		
		RGB_OutputTransf(labdn, dst, dnparams);
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		delete labin;
		delete labdn;
		
	}//end of main RGB_denoise
	

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGBtile_denoise (float ** fLblox, int vblproc, int hblproc, int numblox_H, int numblox_W, \
										   float noisevar_L, float noisevar_ab )	//for DCT
	{
		int vblprocmod=vblproc%8;
		
		float * nbrwt  = new float[TS*TS];	//for DCT
		
		int blkstart = hblproc*TS*TS;
		
		boxabsblur(fLblox[vblprocmod]+blkstart, nbrwt, 3, 3, TS, TS);//blur neighbor weights for more robust estimation	//for DCT
		
		for (int n=0; n<TS*TS; n++) {		//for DCT
			fLblox[vblprocmod][blkstart+n] *= (1-expf(-SQR(nbrwt[n])/noisevar_L));
		}//output neighbor averaged result
		
		delete[] nbrwt;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//printf("vblk=%d  hlk=%d  wsqave=%f   ||   ",vblproc,hblproc,wsqave);
		
	}//end of function tile_denoise
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGBoutput_tile_row (float *bloxrow_L, LabImage * labdn, float ** tilemask_out, int height, int width, int top )
	{
		const int numblox_W = ceil(((float)(width))/(offset))+2*blkrad;
		//const float FTnorm = 1.0f/(TS*TS); //for FT
		const float DCTnorm = 1.0f/(4*TS*TS); //for DCT

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
					
					labdn->L[top+i][left+j] += tilemask_out[i][j]*bloxrow_L[(indx + i)*TS+j]*DCTnorm; //for DCT
					
				}
		}
		
	}
	
#undef TS
#undef fTS
#undef offset
	//#undef eps
	
	
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
			
			float threshsq = noisevar;//*SQR(MadMax(WaveletCoeffs.level_coeffs(lvl)[3], max, Wlvl*Hlvl));
				
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
		int max;
		
		printf("level=%d  ",level);
		for (int dir=1; dir<4; dir++) {
			float mad = SQR(MadMax(WavCoeffs[dir], max, W*H));//*6*/(level+1);
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
				
			}
		}
		printf("\n");
		delete[] sigma;
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	float ImProcFunctions::MadMax(float * DataList, int & max, int datalen) {
		
		//computes Median Absolute Deviation and Maximum of DataList
		//DataList values should mostly have abs val < 65535
		
		int * histo = new int[65536];
		//memset(histo, 0, 65536*sizeof(histo));
		for (int i=0; i<65536; i++) histo[i]=0;
		
		//calculate histogram of absolute values of HH wavelet coeffs
		for (int i=0; i<datalen; i++) {
			histo[MAX(0,MIN(65535,abs((int)DataList[i])))]++;
		}
		
		//find median of histogram
		int median=0, count=0;
		while (count<datalen/2) {
			count += histo[median];
			median++;
		}
		
		//find max of histogram
		max=65535;
		while (histo[max]==0) {
			max--;
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
		float * sfave = new float[W_L*H_L];
		int max;

		printf("\n level=%d  \n",level);
		
		for (int dir=1; dir<4; dir++) {
			float madL = SQR(MadMax(WavCoeffs_L[dir], max, W_L*H_L));	
			float mada = SQR(MadMax(WavCoeffs_a[dir], max, W_ab*H_ab));
			float madb = SQR(MadMax(WavCoeffs_b[dir], max, W_ab*H_ab));
			
			//float thresh_L = sqrt(mad_L*noisevar_L);
			//float thresh_a = sqrt(mad_a*noisevar_ab);
			//float thresh_b = sqrt(mad_b*noisevar_ab);

			printf("  dir=%d  mad_L=%f	mad_a=%f	mad_b=%f	\n",dir,sqrt(madL),sqrt(mada),sqrt(madb));
			
			//float mad_L = noisevar_L *6/((level+2)*pow(2.0f,level));	
			//float mad_a = noisevar_ab;
			//float mad_b = noisevar_ab;
			float mad_L = madL*noisevar_L *6/(level+2);
			float mad_a = mada*noisevar_ab;
			float mad_b = madb*noisevar_ab;
			
			for (int i=0; i<H_ab; i++) {
				for (int j=0; j<W_ab; j++) {
					
					int coeffloc_ab = i*W_ab+j;
					int coeffloc_L	= ((i*skip_L)/skip_ab)*W_L + ((j*skip_L)/skip_ab);
					
					float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L ])+eps;
					float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
					float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;
					
					//float edgefactor = exp(-mag_L/(3*mad_L))) * exp(-mag_a/(3*mad_a)) * exp(-mag_b/(3*mad_b));
					//float edgefactor = 1-exp(-mag_L/(9*mad_L));
					
					//WavCoeffs_a[dir][coeffloc_ab] *= mag_a/(mag_a + noisevar_ab*mad_a*edgefactor + eps);
					//WavCoeffs_b[dir][coeffloc_ab] *= mag_b/(mag_b + noisevar_ab*mad_b*edgefactor + eps);
					
					//float coeff_a = fabs(WavCoeffs_a[dir][coeffloc_ab]);
					//float coeff_b = fabs(WavCoeffs_b[dir][coeffloc_ab]);

					// 'firm' threshold of chroma coefficients
					WavCoeffs_a[dir][coeffloc_ab] *= SQR(1-exp(-(mag_a/mad_a)-(mag_L/(9*madL))));//(coeff_a>2*thresh_a ? 1 : (coeff_a<thresh_a ? 0 : (coeff_a/thresh_a - 1)));
					WavCoeffs_b[dir][coeffloc_ab] *= SQR(1-exp(-(mag_b/mad_b)-(mag_L/(9*madL))));//(coeff_b>2*thresh_b ? 1 : (coeff_b<thresh_b ? 0 : (coeff_b/thresh_b - 1)));

					//WavCoeffs_b[dir][coeffloc_ab] *= (fabs(WavCoeffs_b[dir][coeffloc_ab])<thresh_a*noise_ab ? 0 : 1);


				}
			}
			
			for (int i=0; i<W_L*H_L; i++) {
				
				//float coeff_L = fabs(WavCoeffs_L[dir][i]);
				// 'firm' threshold of luma coefficients
				//float shrinkfactor = (coeff_L>2*thresh_L ? 1 : (coeff_L<thresh_L ? 0 : (coeff_L/thresh_L - 1)));
				
				float mag = SQR(WavCoeffs_L[dir][i]);
				float shrinkfactor = mag/(mag+mad_L*exp(-mag/(9*mad_L))+eps);
				//float shrinkfactor = SQR(1-exp(-(mag/(mad_L))));

				//float shrinkfactor = mag/(mag+noisevar*SQR(sfave[coeffloc])+eps);
				
				//WavCoeffs_L[dir][i] *= shrinkfactor;
				sfave[i] = shrinkfactor;
			}
			
			boxblur(sfave, sfave, level+2, level+2, W_L, H_L);//increase smoothness by locally averaging shrinkage
			for (int i=0; i<W_L*H_L; i++) {
				
				//float coeff_L = fabs(WavCoeffs_L[dir][i]);
				// 'firm' threshold of chroma coefficients
				//float sf = (coeff_L>2*thresh_L ? 1 : (coeff_L<thresh_L ? 0 : (coeff_L/thresh_L - 1)));
				
				float mag = SQR(WavCoeffs_L[dir][i]);
				float sf = mag/(mag+mad_L*exp(-mag/(9*mad_L))+eps);
				//float sf = SQR(1-exp(-(mag/(mad_L))));
				//float sf = mag/(mag+noisevar_L*mad_L);
				
				//use smoothed shrinkage unless local shrinkage is much less
				WavCoeffs_L[dir][i] *= (SQR(sfave[i])+SQR(sf))/(sfave[i]+sf+eps);
				
			}//now luminance coeffs are denoised
			
			

		}
		delete[] sfave;
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	void ImProcFunctions::WaveletDenoiseAll_BiShrink(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a, 
													 wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_ab ) 
	{
		int maxlvl = WaveletCoeffs_L.maxlevel();
		const float eps = 0.01f;
		int max;
		float parfrac = 0.05;
		
		float madL[8][3], mada[8][3], madb[8][3];
		
		for (int lvl=0; lvl<maxlvl; lvl++) {
			// compute median absolute deviation (MAD) of detail coefficients as robust noise estimator
			
			int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
			int Hlvl_L = WaveletCoeffs_L.level_H(lvl);
			
			int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
			int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);
			
			float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);
			
			for (int dir=1; dir<4; dir++) {
				madL[lvl][dir-1] = SQR(MadMax(WavCoeffs_L[dir], max, Wlvl_L*Hlvl_L));
				mada[lvl][dir-1] = SQR(MadMax(WavCoeffs_a[dir], max, Wlvl_ab*Hlvl_ab));
				madb[lvl][dir-1] = SQR(MadMax(WavCoeffs_b[dir], max, Wlvl_ab*Hlvl_ab));
			}
		}
		
		for (int lvl=maxlvl-1; lvl>=0; lvl--) {//for levels less than max, use level diff to make edge mask
			//for (int lvl=0; lvl<maxlvl; lvl++) {
			
			int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
			int Hlvl_L = WaveletCoeffs_L.level_H(lvl);
			
			int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
			int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);
			
			float skip_L = WaveletCoeffs_L.level_stride(lvl);
			float skip_ab = WaveletCoeffs_a.level_stride(lvl);
			
			float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);
			
			if (lvl==maxlvl-1) {
				ShrinkAll(WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, lvl, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,
						  skip_L, skip_ab, noisevar_L, noisevar_ab);//TODO: this implies redundant evaluation of MAD
			} else {
				
				float ** WavPars_L = WaveletCoeffs_L.level_coeffs(lvl+1);
				//float ** WavPars_a = WaveletCoeffs_a.level_coeffs(lvl+1);
				//float ** WavPars_b = WaveletCoeffs_b.level_coeffs(lvl+1);
				
				//simple wavelet shrinkage
				float * sfave = new float[Wlvl_L*Hlvl_L]; 
				//float * edge = new float[Wlvl_L*Hlvl_L];
				array2D<float> edge(Wlvl_L,Hlvl_L);
				AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(Wlvl_L,Hlvl_L));

				printf("\n level=%d  \n",lvl);
				
				for (int dir=1; dir<4; dir++) {
					float mad_L = madL[lvl][dir-1];
					float mad_a = noisevar_ab*mada[lvl][dir-1];
					float mad_b = noisevar_ab*madb[lvl][dir-1];
					//float mad_Lpar = madL[lvl+1][dir-1];
					//float mad_apar = mada[lvl+1][dir-1];
					//float mad_bpar = mada[lvl+1][dir-1];
					
					//float skip_ab_ratio = WaveletCoeffs_a.level_stride(lvl+1)/skip_ab;
					float skip_L_ratio =  WaveletCoeffs_L.level_stride(lvl+1)/skip_L;
					
					printf("  dir=%d  mad_L=%f		mad_a=%f		mad_b=%f	\n",dir,sqrt(mad_L),sqrt(mad_a),sqrt(mad_b));
					
					for (int i=0; i<Hlvl_ab; i++) {
						for (int j=0; j<Wlvl_ab; j++) {
							
							int coeffloc_ab = i*Wlvl_ab+j;
							//int coeffloc_abpar = (MAX(0,i-skip_ab)*Wlvl_ab+MAX(0,j-skip_ab))/skip_ab_ratio;
							
							int coeffloc_L	= ((i*skip_L)/skip_ab)*Wlvl_L + ((j*skip_L)/skip_ab);

							float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L ])+eps;
							float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
							float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;
							
							//float edgefactor = 1-exp(-mag_L/(9*mad_L));// * exp(-mag_a/(4*mad_a)) * exp(-mag_b/(4*mad_b));
							
							//float coeff_a = sqrt(SQR(WavCoeffs_a[dir][coeffloc_ab])/mad_a+SQR(parfrac*WavPars_a[dir][coeffloc_abpar])/mad_apar);
							//float coeff_b = sqrt(SQR(WavCoeffs_b[dir][coeffloc_ab])/mad_b+SQR(parfrac*WavPars_b[dir][coeffloc_abpar])/mad_bpar);
							
							// 'firm' threshold of chroma coefficients
							//WavCoeffs_a[dir][coeffloc_ab] *= edgefactor*(coeff_a>2 ? 1 : (coeff_a<1 ? 0 : (coeff_a - 1)));
							//WavCoeffs_b[dir][coeffloc_ab] *= edgefactor*(coeff_b>2 ? 1 : (coeff_b<1 ? 0 : (coeff_b - 1)));
							
							WavCoeffs_a[dir][coeffloc_ab] *= SQR(1-exp(-(mag_a/mad_a)-(mag_L/(9*mad_L))));
							WavCoeffs_b[dir][coeffloc_ab] *= SQR(1-exp(-(mag_b/mad_b)-(mag_L/(9*mad_L))));

						}
					}//now chrominance coefficients are denoised
					
					mad_L *= noisevar_L*5/(lvl+1);
										
					for (int i=0; i<Hlvl_L; i++) 
						for (int j=0; j<Wlvl_L; j++) {
							
							int coeffloc_L = i*Wlvl_L+j;
							int coeffloc_Lpar = (MAX(0,i-skip_L)*Wlvl_L+MAX(0,j-skip_L))/skip_L_ratio;
							
							float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
							//float mag_Lpar = SQR(parfrac*WavPars_L[dir][coeffloc_Lpar]);
							//float sf_L = SQR(1-expf(-(mag_L/mad_L)-(mag_Lpar/mad_L)));
							float sf_L = mag_L/(mag_L+mad_L*exp(-mag_L/(9*mad_L))+eps);

							sfave[coeffloc_L] = sf_L;
							
							edge[i][j] = (WavCoeffs_L[dir][coeffloc_L] - WavPars_L[dir][coeffloc_Lpar]);
						}
					
					//blur edge measure
					gaussHorizontal<float> (edge, edge, buffer, Wlvl_L, Hlvl_L, 1<<(lvl+1), false /*multiThread*/);
					gaussVertical<float>   (edge, edge, buffer, Wlvl_L, Hlvl_L, 1<<(lvl+1), false);
										
					boxblur(sfave, sfave, lvl+2, lvl+2, Wlvl_L, Hlvl_L);//increase smoothness by locally averaging shrinkage
					
					for (int i=0; i<Hlvl_L; i++) 
						for (int j=0; j<Wlvl_L; j++) {
							
							int coeffloc_L = i*Wlvl_L+j;
							
							float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
							//float sf_L = SQR(1-expf(-(mag_L/mad_L)-(mag_Lpar/mad_L)));
							
							float edgefactor = 1;//expf(-SQR(edge[i][j])/mad_L);
							
							float sf_L = mag_L/(mag_L + edgefactor*mad_L*exp(-mag_L/(9*mad_L))+eps);

							//use smoothed shrinkage unless local shrinkage is much less
							WavCoeffs_L[dir][coeffloc_L] *= (SQR(edgefactor*sfave[coeffloc_L])+SQR(sf_L))/(edgefactor*sfave[coeffloc_L]+sf_L+eps);
							
						}//now luminance coeffs are denoised
					
					
					
				}
				delete[] sfave;
				//delete[] edge;
				delete buffer;

			}
		}
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
};
