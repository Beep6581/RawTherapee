////////////////////////////////////////////////////////////////
//
//			CFA denoise by wavelet transform, FT filtering
//
//	copyright (c) 2008-2012  Emil Martinec <ejmartin@uchicago.edu>
//
//
//  code dated: March 9, 2012
//
//	FTblockDN.cc is free software: you can redistribute it and/or modify
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
#include "rt_math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cplx_wavelet_dec.h"

//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
//#define LIM(x,min,max) MAX(min,MIN(x,max))
//#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
//#define CLIP(x) LIM(x,0,65535)

#define TS 64		// Tile size
#define offset 25	// shift between tiles
#define fTS ((TS/2+1))	// second dimension of Fourier tiles
#define blkrad 1	// radius of block averaging

#define epsilon 0.001f/(TS*TS) //tolerance

namespace rtengine {
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*
	 Structure of the algorithm: 
	 
	 1. Compute an initial denoise of the image via undecimated wavelet transform 
	 and universal thresholding modulated by user input.
	 2. Decompose the residual image into TSxTS size tiles, shifting by 'offset' each step 
	 (so roughly each pixel is in (TS/offset)^2 tiles); Discrete Cosine transform the tiles.
	 3. Filter the DCT data to pick out patterns missed by the wavelet denoise
	 4. Inverse DCT the denoised tile data and combine the tiles into a denoised output image.
	 
	 */
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	
	
	
	void ImProcFunctions::RGB_denoise(Imagefloat * src, Imagefloat * dst, const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe)
	{	
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		/*if (plistener) {
		 plistener->setProgressStr ("Denoise...");
		 plistener->setProgress (0.0);
		 }*/
		
		volatile double progress = 0.0;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		const short int imheight=src->height, imwidth=src->width;
		
		if (dnparams.luma==0 && dnparams.chroma==0) {
            //nothing to do; copy src to dst
            memcpy(dst->data,src->data,dst->width*dst->height*3*sizeof(float));
			return;
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// gamma transform for input data
		float gam = dnparams.gamma;
		float gamthresh = 0.001;
		float gamslope = exp(log((double)gamthresh)/gam)/gamthresh;
		
		LUTf gamcurve(65536,0);
		
		for (int i=0; i<65536; i++) {
			gamcurve[i] = (Color::gamma((double)i/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)) * 32768.0f;
		}
		
		// inverse gamma transform for output data
		float igam = 1/gam;
		float igamthresh = gamthresh*gamslope;
		float igamslope = 1/gamslope;
		
		LUTf igamcurve(65536,0);
		
		for (int i=0; i<65536; i++) {
			igamcurve[i] = (Color::gamma((float)i/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		//srand((unsigned)time(0));//test with random data
		
		const float gain = pow (2.0, dnparams.expcomp);
		
		float noisevar_Ldetail = SQR((SQR(100-dnparams.Ldetail) + 50*(100-dnparams.Ldetail)) * TS * 0.5f);
		
		
		array2D<float> tilemask_in(TS,TS);
		array2D<float> tilemask_out(TS,TS);
		
		const int border = MAX(2,TS/16);

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<TS; i++) {
			float i1 = abs((i>TS/2 ? i-TS+1 : i));
			float vmask = (i1<border ? SQR(sin((M_PI*i1)/(2*border))) : 1.0f);
			float vmask2 = (i1<2*border ? SQR(sin((M_PI*i1)/(2*border))) : 1.0f);
			for (int j=0; j<TS; j++) {
				float j1 = abs((j>TS/2 ? j-TS+1 : j));
				tilemask_in[i][j] = (vmask * (j1<border ? SQR(sin((M_PI*j1)/(2*border))) : 1.0f)) + epsilon;
				
				tilemask_out[i][j] = (vmask2 * (j1<2*border ? SQR(sin((M_PI*j1)/(2*border))) : 1.0f)) + epsilon;
				
			}
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// begin tile processing of image
		
		//output buffer
		Imagefloat * dsttmp = new Imagefloat(imwidth,imheight);
		for (int n=0; n<3*imwidth*imheight; n++) dsttmp->data[n] = 0;
		
		const int tilesize = 1024;
		const int overlap = 128;

		int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;

		if (imwidth<tilesize) {
			numtiles_W = 1;
			tileWskip = imwidth;
			tilewidth = imwidth;
		} else {
			numtiles_W = ceil(((float)(imwidth))/(tilesize-overlap));
			tilewidth  = ceil(((float)(imwidth))/(numtiles_W))+overlap;
			tilewidth += (tilewidth&1);
			tileWskip = tilewidth-overlap;
		}
		if (imheight<tilesize) {
			numtiles_H = 1;
			tileHskip = imheight;
			tileheight = imheight;
		} else {
			numtiles_H = ceil(((float)(imheight))/(tilesize-overlap));
			tileheight = ceil(((float)(imheight))/(numtiles_H))+overlap;
			tileheight += (tileheight&1);
			tileHskip = tileheight-overlap;
		}
		//now we have tile dimensions, overlaps
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//adding omp here slows it down
		for (int tiletop=0; tiletop<imheight; tiletop+=tileHskip) {
			for (int tileleft=0; tileleft<imwidth; tileleft+=tileWskip) {
				
				int tileright = MIN(imwidth,tileleft+tilewidth);
				int tilebottom = MIN(imheight,tiletop+tileheight);
				int width  = tileright-tileleft;
				int height = tilebottom-tiletop;
				
				//input L channel
				array2D<float> Lin(width,height,ARRAY2D_CLEAR_DATA);
				//wavelet denoised image
				LabImage * labdn = new LabImage(width,height);
				//residual between input and denoised L channel
				array2D<float> Ldetail(width,height,ARRAY2D_CLEAR_DATA);
				//pixel weight
				array2D<float> totwt(width,height,ARRAY2D_CLEAR_DATA);//weight for combining DCT blocks

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//TODO: implement using AlignedBufferMP
				//fill tile from image; convert RGB to "luma/chroma"
				for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
					int i1 = i - tiletop;
					for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
						int j1 = j - tileleft;
						
						float X = gain*src->r[i][j];//xyz_prophoto[0][0]*src->r[i][j] + xyz_prophoto[0][1]*src->g[i][j] + xyz_prophoto[0][2]*src->b[i][j];
						float Y = gain*src->g[i][j];//xyz_prophoto[1][0]*src->r[i][j] + xyz_prophoto[1][1]*src->g[i][j] + xyz_prophoto[1][2]*src->b[i][j];
						float Z = gain*src->b[i][j];//xyz_prophoto[2][0]*src->r[i][j] + xyz_prophoto[2][1]*src->g[i][j] + xyz_prophoto[2][2]*src->b[i][j];
						
						X = X<65535.0f ? gamcurve[X] : (Color::gamma((double)X/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						Y = Y<65535.0f ? gamcurve[Y] : (Color::gamma((double)Y/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						Z = Z<65535.0f ? gamcurve[Z] : (Color::gamma((double)Z/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
						
						labdn->L[i1][j1] = Y;
						labdn->a[i1][j1] = (X-Y);
						labdn->b[i1][j1] = (Y-Z);
						
						Ldetail[i1][j1] = 0;
						Lin[i1][j1] = Y;
						totwt[i1][j1] = 0;
					}
				}
				
				//initial impulse denoise
				if (dnparams.luma>0.01) {
					impulse_nr (labdn, MIN(50.0f,dnparams.luma)/20.0f);
				}
				
				int datalen = labdn->W * labdn->H;
				
				//now perform basic wavelet denoise
				//last two arguments of wavelet decomposition are max number of wavelet decomposition levels; 
				//and whether to subsample the image after wavelet filtering.  Subsampling is coded as 
				//binary 1 or 0 for each level, eg subsampling = 0 means no subsampling, 1 means subsample 
				//the first level only, 7 means subsample the first three levels, etc.
				wavelet_decomposition Ldecomp(labdn->data, labdn->W, labdn->H, 5/*maxlevels*/, 0/*subsampling*/ );
				wavelet_decomposition adecomp(labdn->data+datalen, labdn->W, labdn->H, 5, 1 );
				wavelet_decomposition bdecomp(labdn->data+2*datalen, labdn->W, labdn->H, 5, 1 );
				
				float noisevarL	 = SQR((dnparams.luma/125.0f)*(1+ dnparams.luma/25.0f));
				float noisevarab = SQR(dnparams.chroma/10.0f);
				
				//WaveletDenoiseAll_BiShrink(Ldecomp, adecomp, bdecomp, noisevarL, noisevarab);
				WaveletDenoiseAll(Ldecomp, adecomp, bdecomp, noisevarL, noisevarab);

				Ldecomp.reconstruct(labdn->data);
				adecomp.reconstruct(labdn->data+datalen);
				bdecomp.reconstruct(labdn->data+2*datalen);
				
				//TODO: at this point wavelet coefficients storage can be freed
				
				//second impulse denoise
				if (dnparams.luma>0.01) {
					impulse_nr (labdn, MIN(50.0f,dnparams.luma)/20.0f);
				}				
				//PF_correct_RT(dst, dst, defringe.radius, defringe.threshold);
				
				//wavelet denoised L channel
				array2D<float> Lwavdn(width,height);
				float * Lwavdnptr = Lwavdn;
				memcpy (Lwavdnptr, labdn->data, width*height*sizeof(float));
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// now do detail recovery using block DCT to detect  
				// patterns missed by wavelet denoise
				// blocks are not the same thing as tiles!
				
				
				// allocate DCT data structures
				
				// calculation for detail recovery blocks
				const int numblox_W = ceil(((float)(width))/(offset))+2*blkrad;
				const int numblox_H = ceil(((float)(height))/(offset))+2*blkrad;
				//const int nrtiles = numblox_W*numblox_H;
				// end of tiling calc
				
				//DCT block data storage
				float * Lblox  = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));
				float * fLblox = (float *) fftwf_malloc (numblox_W*TS*TS * sizeof (float));
				
				
				//make a plan for FFTW
				fftwf_plan plan_forward_blox, plan_backward_blox;
				
				int nfwd[2]={TS,TS};
				
				//for DCT:
				const fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
				const fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};
				
				plan_forward_blox  = fftwf_plan_many_r2r(2, nfwd, numblox_W, Lblox, NULL, 1, TS*TS, fLblox, NULL, 1, TS*TS, fwdkind, FFTW_ESTIMATE );
				plan_backward_blox = fftwf_plan_many_r2r(2, nfwd, numblox_W, fLblox, NULL, 1, TS*TS, Lblox, NULL, 1, TS*TS, bwdkind, FFTW_ESTIMATE );
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// Main detail recovery algorithm: Block loop
//OpenMP here			
//adding omp here leads to artifacts	
				for (int vblk=0; vblk<numblox_H; vblk++) {
					//printf("vblock=%d",vblk);
					int vblkmod = vblk%8;
					
					int top = (vblk-blkrad)*offset;
					
					float * buffer = new float [width + TS + 2*blkrad*offset];
					float * datarow = buffer+blkrad*offset;

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//TODO: implement using AlignedBufferMP
					for (int i=0/*, row=top*/; i<TS; i++/*, row++*/) {
						int row = top + i;
						int rr = row;
						if (row<0) {
							rr = MIN(-row,height-1);
						} else if (row>=height) {
							rr = MAX(0,2*height-2-row);
						} 
						
						for (int j=0; j<labdn->W; j++) {
							datarow[j] = (Lin[rr][j]-Lwavdn[rr][j]);
						}
						
						for (int j=-blkrad*offset; j<0; j++) {
							datarow[j] = datarow[MIN(-j,width-1)];
						}
						for (int j=width; j<width+TS+blkrad*offset; j++) {
							datarow[j] = datarow[MAX(0,2*width-2-j)];
						}//now we have a padded data row
						
						//now fill this row of the blocks with Lab high pass data
//OMP here does not add speed, better handled on the outside loop
						for (int hblk=0; hblk<numblox_W; hblk++) {
							int left = (hblk-blkrad)*offset;
							int indx = (hblk)*TS;//index of block in malloc
							
							for (int j=0; j<TS; j++) {
								Lblox[(indx + i)*TS+j]  = tilemask_in[i][j]*datarow[left+j];// luma data
								if (top+i>=0 && top+i<height && left+j>=0 && left+j<width) {
									totwt[top+i][left+j] += tilemask_in[i][j]*tilemask_out[i][j];
								}
							}
						}
					}//end of filling block row
					delete[] buffer;
					
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					
					//fftwf_print_plan (plan_forward_blox);
					fftwf_execute_r2r(plan_forward_blox,Lblox,fLblox);		// DCT an entire row of tiles
					
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					// now process the vblk row of blocks for noise reduction
					for (int hblk=0; hblk<numblox_W; hblk++) {
						
						RGBtile_denoise (fLblox, vblk, hblk, numblox_H, numblox_W, noisevar_Ldetail );
						
					}//end of horizontal block loop
					
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					
					//now perform inverse FT of an entire row of blocks
					fftwf_execute_r2r(plan_backward_blox,fLblox,Lblox);	//for DCT
					
					int topproc = (vblk-blkrad)*offset;
					
					//add row of blocks to output image tile
					RGBoutput_tile_row (Lblox, Ldetail, tilemask_out, height, width, topproc );
					
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					
				}//end of vertical block loop
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				// clean up
				//#pragma omp single nowait
				fftwf_destroy_plan( plan_forward_blox );
				//#pragma omp single nowait
				fftwf_destroy_plan( plan_backward_blox );
								
				fftwf_free ( Lblox);
				fftwf_free ( fLblox);
				
				fftwf_cleanup();
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				for (int i=0; i<height; i++) {
					for (int j=0; j<width; j++) {
						//may want to include masking threshold for large hipass data to preserve edges/detail
						float hpdn = Ldetail[i][j]/totwt[i][j];//note that labdn initially stores the denoised hipass data
						
						labdn->L[i][j] = Lwavdn[i][j] + hpdn;
						
					}
				}
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// transform denoised "Lab" to output RGB
				
				//calculate mask for feathering output tile overlaps
				float * Vmask = new float [height];
				float * Hmask = new float [width];

				for (int i=0; i<height; i++) {
					Vmask[i] = 1;
				}
				for (int j=0; j<width; j++) {
					Hmask[j] = 1;
				}
				for (int i=0; i<overlap; i++) {
					float mask = SQR(sin((M_PI*i)/(2*overlap)));
					if (tiletop>0) Vmask[i] = mask;
					if (tilebottom<imheight) Vmask[height-1-i] = mask;
					if (tileleft>0) Hmask[i] = mask;
					if (tileright<imwidth) Hmask[width-1-i] = mask;
				}

#ifdef _OPENMP
#pragma omp parallel for
#endif
				//convert back to RGB and write to destination array
				for (int i=tiletop/* i1=0*/; i<tilebottom; i++/*, i1++*/){
					int i1 = i-tiletop;
					float X,Y,Z;
					for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
						int j1=j-tileleft;
						
						Y = labdn->L[i1][j1];
						X = (labdn->a[i1][j1]) + Y;
						Z = Y - (labdn->b[i1][j1]);
						
						X = X<32768.0f ? igamcurve[X] : (Color::gamma((float)X/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
						Y = Y<32768.0f ? igamcurve[Y] : (Color::gamma((float)Y/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
						Z = Z<32768.0f ? igamcurve[Z] : (Color::gamma((float)Z/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
						
						//Y = 65535.0f*(0.05+0.1*((float)rand()/(float)RAND_MAX));//test with random data
						
						float factor = Vmask[i1]*Hmask[j1]/gain;
						
						dsttmp->r[i][j] += factor*X;//prophoto_xyz[0][0]*X + prophoto_xyz[0][1]*Y + prophoto_xyz[0][2]*Z;
						dsttmp->g[i][j] += factor*Y;//prophoto_xyz[1][0]*X + prophoto_xyz[1][1]*Y + prophoto_xyz[1][2]*Z;
						dsttmp->b[i][j] += factor*Z;//prophoto_xyz[2][0]*X + prophoto_xyz[2][1]*Y + prophoto_xyz[2][2]*Z;
						
					}
				}
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				delete labdn;
				delete[] Vmask;
				delete[] Hmask;
				
			}//end of tile row
		}//end of tile loop
		
	//copy denoised image to output
	memcpy (dst->data, dsttmp->data, 3*dst->width*dst->height*sizeof(float));
	
	delete dsttmp;
	
}//end of main RGB_denoise


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	
	void ImProcFunctions::RGBtile_denoise (float * fLblox, int vblproc, int hblproc, int numblox_H, int numblox_W, float noisevar_Ldetail )	//for DCT
	{			
		float * nbrwt  = new float[TS*TS];	//for DCT
		
		int blkstart = hblproc*TS*TS;
		
		boxabsblur(fLblox+blkstart, nbrwt, 3, 3, TS, TS);//blur neighbor weights for more robust estimation	//for DCT

		for (int n=0; n<TS*TS; n++) {		//for DCT
			fLblox[blkstart+n] *= (1-expf(-SQR(nbrwt[n])/noisevar_Ldetail));
		}//output neighbor averaged result
		
		delete[] nbrwt;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//printf("vblk=%d  hlk=%d  wsqave=%f   ||   ",vblproc,hblproc,wsqave);
		
	}//end of function tile_denoise
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::RGBoutput_tile_row (float *bloxrow_L, float ** Ldetail, float ** tilemask_out, int height, int width, int top )
	{
		const int numblox_W = ceil(((float)(width))/(offset));
		const float DCTnorm = 1.0f/(4*TS*TS); //for DCT

#ifdef _OPENMP
#pragma omp parallel for
#endif			
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
					
					Ldetail[top+i][left+j] += tilemask_out[i][j]*bloxrow_L[(indx + i)*TS+j]*DCTnorm; //for DCT
					
				}
		}
		
	}
	
#undef TS
#undef fTS
#undef offset
#undef epsilon
	
	
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
	
	
	void ImProcFunctions::WaveletDenoiseAll_BiShrink(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a, 
													 wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_ab ) 
	{
		int maxlvl = WaveletCoeffs_L.maxlevel();
		const float eps = 0.01f;
		int max;
		float parfrac = 0.05;
		
		float madL[8][3], mada[8][3], madb[8][3];

//OpenMP here		
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
					
					if (noisevar_ab>0.01) {
						
						printf("  dir=%d  mad_L=%f		mad_a=%f		mad_b=%f	\n",dir,sqrt(mad_L),sqrt(mad_a),sqrt(mad_b));
						
//OpenMP here						
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
								
								//float satfactor_a = mad_a/(mad_a+0.5*SQR(WavCoeffs_a[0][coeffloc_ab]));
								//float satfactor_b = mad_b/(mad_b+0.5*SQR(WavCoeffs_b[0][coeffloc_ab]));
								
								WavCoeffs_a[dir][coeffloc_ab] *= SQR(1-exp(-(mag_a/mad_a)-(mag_L/(9*mad_L)))/*satfactor_a*/);
								WavCoeffs_b[dir][coeffloc_ab] *= SQR(1-exp(-(mag_b/mad_b)-(mag_L/(9*mad_L)))/*satfactor_b*/);
								
							}
						}//now chrominance coefficients are denoised
					}
					
					if (noisevar_L>0.01) {
						mad_L *= noisevar_L*5/(lvl+1);
//OpenMP here						
						for (int i=0; i<Hlvl_L; i++) 
							for (int j=0; j<Wlvl_L; j++) {
								
								int coeffloc_L = i*Wlvl_L+j;
								int coeffloc_Lpar = (MAX(0,i-skip_L)*Wlvl_L+MAX(0,j-skip_L))/skip_L_ratio;
								
								float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
								//float mag_Lpar = SQR(parfrac*WavPars_L[dir][coeffloc_Lpar]);
								//float sf_L = SQR(1-expf(-(mag_L/mad_L)-(mag_Lpar/mad_L)));
								float sf_L = mag_L/(mag_L+mad_L*exp(-mag_L/(9*mad_L))+eps);
								
								sfave[coeffloc_L] = sf_L;
								
								//edge[i][j] = (WavCoeffs_L[dir][coeffloc_L] - WavPars_L[dir][coeffloc_Lpar]);
							}
						
						//blur edge measure
						//gaussHorizontal<float> (edge, edge, buffer, Wlvl_L, Hlvl_L, 1<<(lvl+1), false /*multiThread*/);
						//gaussVertical<float>   (edge, edge, buffer, Wlvl_L, Hlvl_L, 1<<(lvl+1), false);
						
						boxblur(sfave, sfave, lvl+2, lvl+2, Wlvl_L, Hlvl_L);//increase smoothness by locally averaging shrinkage
//OpenMP here						
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
					
				}
				delete[] sfave;
				delete buffer;
				
			}
		}
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void ImProcFunctions::WaveletDenoiseAll(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a, 
											wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_ab ) 
	{
		int maxlvl = WaveletCoeffs_L.maxlevel();

#ifdef _OPENMP
#pragma omp parallel for
#endif
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
		float * sfavea = new float[W_L*H_L];
		float * sfaveb = new float[W_L*H_L];

		int max;
		
		printf("\n level=%d  \n",level);
	
		for (int dir=1; dir<4; dir++) {
			float madL = SQR(MadMax(WavCoeffs_L[dir], max, W_L*H_L));	
			float mada = SQR(MadMax(WavCoeffs_a[dir], max, W_ab*H_ab));
			float madb = SQR(MadMax(WavCoeffs_b[dir], max, W_ab*H_ab));
			
			
			printf("  dir=%d  mad_L=%f	mad_a=%f	mad_b=%f	\n",dir,sqrt(madL),sqrt(mada),sqrt(madb));

			float mad_L = madL*noisevar_L*5/(level+1);
			float mad_a = mada*noisevar_ab;
			float mad_b = madb*noisevar_ab;
			
			if (noisevar_ab>0.01) {
//OpenMP here

#ifdef _OPENMP
#pragma omp parallel for
#endif				
				for (int i=0; i<H_ab; i++) {
					for (int j=0; j<W_ab; j++) {
						
						int coeffloc_ab = i*W_ab+j;
						int coeffloc_L	= ((i*skip_L)/skip_ab)*W_L + ((j*skip_L)/skip_ab);
						
						float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L ])+eps;
						float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
						float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;
						
						sfavea[coeffloc_ab] = (1-exp(-(mag_a/mad_a)-(mag_L/(9*madL))));
						sfaveb[coeffloc_ab] = (1-exp(-(mag_b/mad_b)-(mag_L/(9*madL))));
						
						// 'firm' threshold of chroma coefficients
						//WavCoeffs_a[dir][coeffloc_ab] *= (1-exp(-(mag_a/mad_a)-(mag_L/(9*madL))));//(coeff_a>2*thresh_a ? 1 : (coeff_a<thresh_a ? 0 : (coeff_a/thresh_a - 1)));
						//WavCoeffs_b[dir][coeffloc_ab] *= (1-exp(-(mag_b/mad_b)-(mag_L/(9*madL))));//(coeff_b>2*thresh_b ? 1 : (coeff_b<thresh_b ? 0 : (coeff_b/thresh_b - 1)));
												
						
					}
				}//now chrominance coefficients are denoised
				
				boxblur(sfavea, sfavea, level+2, level+2, W_ab, H_ab);//increase smoothness by locally averaging shrinkage
				boxblur(sfaveb, sfaveb, level+2, level+2, W_ab, H_ab);//increase smoothness by locally averaging shrinkage

#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (int i=0; i<H_ab; i++) 
					for (int j=0; j<W_ab; j++) {
						
						int coeffloc_ab = i*W_ab+j;
						int coeffloc_L	= ((i*skip_L)/skip_ab)*W_L + ((j*skip_L)/skip_ab);
						
						float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L ])+eps;
						float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
						float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;
						
						float sfa = (1-exp(-(mag_a/mad_a)-(mag_L/(9*madL))));
						float sfb = (1-exp(-(mag_b/mad_b)-(mag_L/(9*madL))));
						
						//use smoothed shrinkage unless local shrinkage is much less
						WavCoeffs_a[dir][coeffloc_ab] *= (SQR(sfavea[coeffloc_ab])+SQR(sfa))/(sfavea[coeffloc_ab]+sfa+eps);
						WavCoeffs_b[dir][coeffloc_ab] *= (SQR(sfaveb[coeffloc_ab])+SQR(sfb))/(sfaveb[coeffloc_ab]+sfb+eps);
						
					}//now chrominance coefficients are denoised
			}
			
			if (noisevar_L>0.01) {
//OpenMP here				
#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (int i=0; i<W_L*H_L; i++) {

					float mag = SQR(WavCoeffs_L[dir][i]);
					float shrinkfactor = mag/(mag+mad_L*exp(-mag/(9*mad_L))+eps);
										
					//WavCoeffs_L[dir][i] *= shrinkfactor;
					sfave[i] = shrinkfactor;
				}
//OpenMP here				
				boxblur(sfave, sfave, level+2, level+2, W_L, H_L);//increase smoothness by locally averaging shrinkage
#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (int i=0; i<W_L*H_L; i++) {
					
					
					float mag = SQR(WavCoeffs_L[dir][i]);
					float sf = mag/(mag+mad_L*exp(-mag/(9*mad_L))+eps);
					
					//use smoothed shrinkage unless local shrinkage is much less
					WavCoeffs_L[dir][i] *= (SQR(sfave[i])+SQR(sf))/(sfave[i]+sf+eps);
					
				}//now luminance coefficients are denoised
			}
			
			
		}
		delete[] sfave;
		delete[] sfavea;
		delete[] sfaveb;

	}

	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
}

