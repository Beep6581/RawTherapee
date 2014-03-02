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
#include "../rtgui/threadutils.h"

//#include "bilateral2.h"
#include "gauss.h"

#include "rtengine.h"
#include "improcfun.h"
#include "LUT.h"
#include "array2D.h"
#include "iccmatrices.h"
#include "boxblur.h"
#include "rt_math.h"
#include "mytime.h"
#include "sleef.c"
#ifdef __SSE2__
	#include "sleefsseavx.c"
#endif

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

		extern const Settings* settings;




	void ImProcFunctions::RGB_denoise(Imagefloat * src, Imagefloat * dst, bool isRAW, const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe, const double expcomp)
	{
//#ifdef _DEBUG
//	MyTime t1e,t2e;
//	t1e.set();
//#endif
	
		static MyMutex FftwMutex;
		MyMutex::MyLock lock(FftwMutex);

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		/*if (plistener) {
		 plistener->setProgressStr ("Denoise...");
		 plistener->setProgress (0.0);
		 }*/

//		volatile double progress = 0.0;

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		const short int imheight=src->height, imwidth=src->width;

		if (dnparams.luma==0 && dnparams.chroma==0) {
            //nothing to do; copy src to dst
            memcpy(dst->data,src->data,dst->width*dst->height*3*sizeof(float));
			return;
		}
		perf=false;
		if(dnparams.dmethod=="RGB") perf=true;//RGB mode
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// gamma transform for input data
		float gam = dnparams.gamma;
		float gamthresh = 0.001f;
		if(!isRAW) {//reduce gamma under 1 for Lab mode ==> TIF and JPG
		 if(gam <1.9f) gam=1.f - (1.9f-gam)/3.f;//minimum gamma 0.7
		 else if (gam >= 1.9f && gam <= 3.f) gam=(1.4f/1.1f)*gam - 1.41818f;
		 }
		float gamslope = exp(log((double)gamthresh)/gam)/gamthresh;

		LUTf gamcurve(65536,0);
		if(perf) {
		for (int i=0; i<65536; i++) {
			gamcurve[i] = (Color::gamma((double)i/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)) * 32768.0f;
		}
		}
		else  {
		for (int i=0; i<65536; i++) {
			gamcurve[i] = (Color::gamman((double)i/65535.0,gam)) * 32768.0f;
		}
		}

		// inverse gamma transform for output data
		float igam = 1.f/gam;
		float igamthresh = gamthresh*gamslope;
		float igamslope = 1.f/gamslope;

		LUTf igamcurve(65536,0);
		if(perf) {
		for (int i=0; i<65536; i++) {
			igamcurve[i] = (Color::gamma((float)i/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
		}
		}
		else {
		for (int i=0; i<65536; i++) {
			igamcurve[i] = (Color::gamman((float)i/32768.0f,igam) * 65535.0f);
		}

		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		//srand((unsigned)time(0));//test with random data

		const float gain = pow (2.0f, float(expcomp));
		float incr=1.f;
		float noisevar_Ldetail = SQR((float)(SQR(100.-dnparams.Ldetail) + 50.*(100.-dnparams.Ldetail)) * TS * 0.5f * incr);
		bool enhance_denoise = dnparams.enhance;
		int gamlab = settings->denoiselabgamma;//gamma lab essentialy for Luminance detail
		if(gamlab > 2) gamlab=2;
		if(settings->verbose) printf("Denoise Lab=%i\n",gamlab);

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

        // According to FFTW-Doc 'it is safe to execute the same plan in parallel by multiple threads', so we now create 4 plans
        // outside the parallel region and use them inside the parallel region.

        // calculate max size of numblox_W.
        int max_numblox_W = ceil(((float)(MIN(imwidth,tilewidth)))/(offset))+2*blkrad;
        // calculate min size of numblox_W.
		int min_numblox_W = ceil(((float)((MIN(imwidth,((numtiles_W - 1) * tileWskip) + tilewidth) ) - ((numtiles_W - 1) * tileWskip)))/(offset))+2*blkrad;

        // these are needed only for creation of the plans and will be freed before entering the parallel loop
		float * Lbloxtmp;
		float * fLbloxtmp;
		Lbloxtmp  = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof (float));
		fLbloxtmp = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof (float));

		int nfwd[2]={TS,TS};

        //for DCT:
		const fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
		const fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};

        fftwf_plan plan_forward_blox[2];
        fftwf_plan plan_backward_blox[2];

        // Creating the plans with FFTW_MEASURE instead of FFTW_ESTIMATE speeds up the execute a bit
		plan_forward_blox[0]  = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, Lbloxtmp, NULL, 1, TS*TS, fLbloxtmp, NULL, 1, TS*TS, fwdkind, FFTW_MEASURE );
		plan_backward_blox[0] = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, fLbloxtmp, NULL, 1, TS*TS, Lbloxtmp, NULL, 1, TS*TS, bwdkind, FFTW_MEASURE );
		plan_forward_blox[1]  = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, Lbloxtmp, NULL, 1, TS*TS, fLbloxtmp, NULL, 1, TS*TS, fwdkind, FFTW_MEASURE );
		plan_backward_blox[1] = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, fLbloxtmp, NULL, 1, TS*TS, Lbloxtmp, NULL, 1, TS*TS, bwdkind, FFTW_MEASURE );
		fftwf_free ( Lbloxtmp );
		fftwf_free ( fLbloxtmp );

#ifdef _OPENMP
        // Calculate number of tiles. If less than omp_get_max_threads(), then limit num_threads to number of tiles
        int numtiles = numtiles_W * numtiles_H;
        int numthreads = MIN(numtiles,omp_get_max_threads());
        if(options.rgbDenoiseThreadLimit > 0) numthreads = MIN(numthreads,options.rgbDenoiseThreadLimit);
        // Issue 1887, overide setting of 1, if more than one thread is available. This way the inner omp-directives should become inactive
        if(numthreads == 1 && omp_get_max_threads() > 1)
			numthreads = 2;
#pragma omp parallel num_threads(numthreads)
#endif
        {
		//DCT block data storage
		float * Lblox;
		float * fLblox;
    TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);

    double wp[3][3] = {
        {wprof[0][0],wprof[0][1],wprof[0][2]},
		{wprof[1][0],wprof[1][1],wprof[1][2]},
		{wprof[2][0],wprof[2][1],wprof[2][2]}
    };

	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	//inverse matrix user select
	double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};



#ifdef _OPENMP
#pragma omp critical
#endif
        {
		Lblox  = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof(float));
		fLblox = (float*) fftwf_malloc(max_numblox_W*TS*TS*sizeof(float));
        }
#ifdef _OPENMP
#pragma omp for schedule(dynamic) collapse(2)
#endif
		for (int tiletop=0; tiletop<imheight; tiletop+=tileHskip) {
			for (int tileleft=0; tileleft<imwidth; tileleft+=tileWskip) {

				int tileright = MIN(imwidth,tileleft+tilewidth);
				int tilebottom = MIN(imheight,tiletop+tileheight);
				int width  = tileright-tileleft;
				int height = tilebottom-tiletop;
				//input L channel
				array2D<float> Lin(width,height);
				//wavelet denoised image
				LabImage * labdn = new LabImage(width,height);

				//residual between input and denoised L channel
				array2D<float> Ldetail(width,height,ARRAY2D_CLEAR_DATA);
				//pixel weight
				array2D<float> totwt(width,height,ARRAY2D_CLEAR_DATA);//weight for combining DCT blocks

				//
				//#ifdef _OPENMP
				//#pragma omp parallel for
				//#endif
				//TODO: implement using AlignedBufferMP
				//fill tile from image; convert RGB to "luma/chroma"
				if (isRAW) {//image is raw; use channel differences for chroma channels
					if(!perf){//lab mode
							//modification Jacques feb 2013					
					for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
						int i1 = i - tiletop;
						for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
							int j1 = j - tileleft;
							float R_ = gain*src->r(i,j);
							float G_ = gain*src->g(i,j);
							float B_ = gain*src->b(i,j);
							//modify arbitrary data for Lab..I have test : nothing, gamma 2.6 11 - gamma 4 5 - gamma 5.5 10
							//we can put other as gamma g=2.6 slope=11, etc.
							// but noting to do with real gamma !!!: it's only for data Lab # data RGB
							//finally I opted fot gamma55 and with options we can change
							if (gamlab == 0) {// options 12/2013
							R_ = Color::igammatab_26_11[R_];
							G_ = Color::igammatab_26_11[G_];
							B_ = Color::igammatab_26_11[B_];
							}
							else if (gamlab == 1) {
							//other new gamma 4 5 
							R_ = Color::igammatab_4[R_];
							G_ = Color::igammatab_4[G_];
							B_ = Color::igammatab_4[B_];
							}
							else if (gamlab == 2) {							
							//new gamma 5.5 10 better for detail luminance..it is a compromise...which depends on the image (distribution BL, ML, HL ...)
							R_ = Color::igammatab_55[R_];
							G_ = Color::igammatab_55[G_];
							B_ = Color::igammatab_55[B_];
							}
							//apply gamma noise	standard (slider)
							R_ = R_<65535.0f ? gamcurve[R_] : (Color::gamman((double)R_/65535.0, gam)*32768.0f);
							G_ = G_<65535.0f ? gamcurve[G_] : (Color::gamman((double)G_/65535.0, gam)*32768.0f);
							B_ = B_<65535.0f ? gamcurve[B_] : (Color::gamman((double)B_/65535.0, gam)*32768.0f);
							//true conversion xyz=>Lab
							float L,a,b;
							float X,Y,Z;
							Color::rgbxyz(R_,G_,B_,X,Y,Z,wp);

							//convert to Lab
							Color::XYZ2Lab(X, Y, Z, L, a, b);
							labdn->L[i1][j1] = L;
							labdn->a[i1][j1] = a;
							labdn->b[i1][j1] = b;
							Lin[i1][j1] = L;
//							totwt[i1][j1] = 0;
						}
					}
					}
					else {//RGB mode
					for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
						int i1 = i - tiletop;
						for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
							int j1 = j - tileleft;

							float X = gain*src->r(i,j);
							float Y = gain*src->g(i,j);
							float Z = gain*src->b(i,j);

							X = X<65535.0f ? gamcurve[X] : (Color::gamma((double)X/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
							Y = Y<65535.0f ? gamcurve[Y] : (Color::gamma((double)Y/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);
							Z = Z<65535.0f ? gamcurve[Z] : (Color::gamma((double)Z/65535.0, gam, gamthresh, gamslope, 1.0, 0.0)*32768.0f);

							labdn->L[i1][j1] = Y;
							labdn->a[i1][j1] = (X-Y);
							labdn->b[i1][j1] = (Y-Z);

//							Ldetail[i1][j1] = 0;
							Lin[i1][j1] = Y;
//							totwt[i1][j1] = 0;
						}
					}

					}
				} else {//image is not raw; use Lab parametrization
					for (int i=tiletop/*, i1=0*/; i<tilebottom; i++/*, i1++*/) {
						int i1 = i - tiletop;
						for (int j=tileleft/*, j1=0*/; j<tileright; j++/*, j1++*/) {
							int j1 = j - tileleft;
							float L,a,b;
							//use gamma sRGB, not good if TIF (JPG) Output profil not with gamma sRGB  (eg : gamma =1.0, or 1.8...)
							//very difficult to solve !
							// solution ==> save TIF with gamma sRGB and re open
							float rtmp = Color::igammatab_srgb[ src->r(i,j) ];
							float gtmp = Color::igammatab_srgb[ src->g(i,j) ];
							float btmp = Color::igammatab_srgb[ src->b(i,j) ];
							//modification Jacques feb 2013
							// gamma slider different from raw
							rtmp = rtmp<65535.0f ? gamcurve[rtmp] : (Color::gamman((double)rtmp/65535.0, gam)*32768.0f);
							gtmp = gtmp<65535.0f ? gamcurve[gtmp] : (Color::gamman((double)gtmp/65535.0, gam)*32768.0f);
							btmp = btmp<65535.0f ? gamcurve[btmp] : (Color::gamman((double)btmp/65535.0, gam)*32768.0f);


							float X,Y,Z;
							Color::rgbxyz(rtmp,gtmp,btmp,X,Y,Z,wp);

							//convert Lab
							Color::XYZ2Lab(X, Y, Z, L, a, b);
							labdn->L[i1][j1] = L;
							labdn->a[i1][j1] = a;
							labdn->b[i1][j1] = b;

//							Ldetail[i1][j1] = 0;
							Lin[i1][j1] = L;

//							totwt[i1][j1] = 0;
						}
					}
				}


				//initial impulse denoise
				if (dnparams.luma>0.01) {
					impulse_nr (labdn, float(MIN(50.0,dnparams.luma))/20.0f);
				}

				int datalen = labdn->W * labdn->H;

				//now perform basic wavelet denoise
				//last two arguments of wavelet decomposition are max number of wavelet decomposition levels;
				//and whether to subsample the image after wavelet filtering.  Subsampling is coded as
				//binary 1 or 0 for each level, eg subsampling = 0 means no subsampling, 1 means subsample
				//the first level only, 7 means subsample the first three levels, etc.
				float noisevarL	 = (float) (SQR((dnparams.luma/125.0)*(1.+ dnparams.luma/25.0)));
				
				float interm_med= (float) dnparams.chroma/10.0;
				float intermred, intermblue;
				if(dnparams.redchro > 0.) intermred=0.0014f* (float)SQR(dnparams.redchro); else intermred= (float) dnparams.redchro/7.0;//increase slower than linear for more sensit
				float intermred2=(float) dnparams.redchro/7.0;
				if(dnparams.bluechro > 0.) intermblue=0.0014f*(float) SQR(dnparams.bluechro); else intermblue= (float)dnparams.bluechro/7.0;//increase slower than linear		
				float intermblue2=(float) dnparams.bluechro/7.0;
				//adjust noise ab in function of sliders red and blue
				float realred = interm_med + intermred; if (realred < 0.f) realred=0.01f;
				float realred2 = interm_med + intermred2; if (realred2 < 0.f) realred2=0.01f;
				float noisevarab_r = SQR(realred);
				float realblue = interm_med + intermblue; if (realblue < 0.f) realblue=0.01f;
				float realblue2 = interm_med + intermblue2; if (realblue2 < 0.f) realblue2=0.01f;
				float noisevarab_b = SQR(realblue);
				
				
                { // enclosing this code in a block frees about 120 MB before allocating 20 MB after this block (measured with D700 NEF)
                wavelet_decomposition* Ldecomp;
                wavelet_decomposition* adecomp;
                wavelet_decomposition* bdecomp;

				int levwav=5;
				float maxreal = max(realred2, realblue2);
				//increase the level of wavelet if user increase much or very much sliders
				if( maxreal < 8.f) levwav=5;
				else if( maxreal < 10.f)levwav=6;
				else if( maxreal < 15.f)levwav=7;
				else levwav=8;//maximum ==> I have increase Maxlevel in cplx_wavelet_dec.h from 8 to 9
				
				
				//	if (settings->verbose) printf("levwavelet=%i  noisevarA=%f noisevarB=%f \n",levwav, noisevarab_r, noisevarab_b );
                Ldecomp = new wavelet_decomposition (labdn->data, labdn->W, labdn->H, levwav/*maxlevels*/, 0/*subsampling*/ );
                adecomp = new wavelet_decomposition (labdn->data+datalen, labdn->W, labdn->H,levwav, 1 );
                bdecomp = new wavelet_decomposition (labdn->data+2*datalen, labdn->W, labdn->H, levwav, 1 );

				if(enhance_denoise)	WaveletDenoiseAll_BiShrink(*Ldecomp, *adecomp, *bdecomp, noisevarL, noisevarab_r, noisevarab_b,labdn);//enhance mode
				else; WaveletDenoiseAll(*Ldecomp, *adecomp, *bdecomp, noisevarL, noisevarab_r, noisevarab_b,labdn);//

				Ldecomp->reconstruct(labdn->data);
				delete Ldecomp;
				adecomp->reconstruct(labdn->data+datalen);
				delete adecomp;
				bdecomp->reconstruct(labdn->data+2*datalen);
				delete bdecomp;
                }

				//TODO: at this point wavelet coefficients storage can be freed
				//Issue 1680: Done now

				//second impulse denoise
				if (dnparams.luma>0.01) {
					impulse_nr (labdn, MIN(50.0f,(float)dnparams.luma)/20.0f);
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



				// calculation for detail recovery blocks
				const int numblox_W = ceil(((float)(width))/(offset))+2*blkrad;
				const int numblox_H = ceil(((float)(height))/(offset))+2*blkrad;

				//const int nrtiles = numblox_W*numblox_H;
				// end of tiling calc
                {

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// Main detail recovery algorithm: Block loop
				//OpenMP here
				//adding omp here leads to artifacts
                AlignedBuffer<float> pBuf(width + TS + 2*blkrad*offset);
				for (int vblk=0; vblk<numblox_H; vblk++) {
					//printf("vblock=%d",vblk);

					int top = (vblk-blkrad)*offset;
					float * datarow = (float*)pBuf.data +blkrad*offset;

					//#ifdef _OPENMP
					//#pragma omp parallel for
					//#endif
					//TODO: implement using AlignedBufferMP
//					#pragma omp parallel for
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

					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					//fftwf_print_plan (plan_forward_blox);
					if(numblox_W == max_numblox_W)
                        fftwf_execute_r2r(plan_forward_blox[0],Lblox,fLblox);		// DCT an entire row of tiles
                    else
                        fftwf_execute_r2r(plan_forward_blox[1],Lblox,fLblox);		// DCT an entire row of tiles
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					// now process the vblk row of blocks for noise reduction
					for (int hblk=0; hblk<numblox_W; hblk++) {

						RGBtile_denoise (fLblox, vblk, hblk, numblox_H, numblox_W, noisevar_Ldetail );

					}//end of horizontal block loop

					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

					//now perform inverse FT of an entire row of blocks
					if(numblox_W == max_numblox_W)
                        fftwf_execute_r2r(plan_backward_blox[0],fLblox,Lblox);	//for DCT
                    else
                        fftwf_execute_r2r(plan_backward_blox[1],fLblox,Lblox);	//for DCT

					int topproc = (vblk-blkrad)*offset;

					//add row of blocks to output image tile
					RGBoutput_tile_row (Lblox, Ldetail, tilemask_out, height, width, topproc );

					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				}//end of vertical block loop
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                }
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
				float * Vmask = new float [height+1];
				float * Hmask = new float [width+1];

				for (int i=0; i<height; i++) {
					Vmask[i] = 1;
				}
				for (int j=0; j<width; j++) {
					Hmask[j] = 1;
				}
				for (int i=0; i<overlap; i++) {
					float mask = SQR(sin((M_PI*i)/(2*overlap)));
					if (tiletop>0) Vmask[i] = mask;
					if (tilebottom<imheight) Vmask[height-i] = mask;
					if (tileleft>0) Hmask[i] = mask;
					if (tileright<imwidth) Hmask[width-i] = mask;
				}

				//convert back to RGB and write to destination array
				if (isRAW) {
				if(!perf) {//Lab mode
#ifdef _OPENMP
//#pragma omp parallel for
#endif
					for (int i=tiletop; i<tilebottom; i++){
						int i1 = i-tiletop;
						float X,Y,Z,L,a,b;

						for (int j=tileleft; j<tileright; j++) {
							int j1=j-tileleft;
							//modification Jacques feb 2013
							//true conversion Lab==>xyz
							L = labdn->L[i1][j1];
							a = labdn->a[i1][j1];
							b = labdn->b[i1][j1];
							//convert XYZ
							Color::Lab2XYZ(L, a, b, X, Y, Z);
							//apply inverse gamma noise
							float r_,g_,b_;
							Color::xyz2rgb(X,Y,Z,r_,g_,b_,wip);
							//inverse gamma standard (slider)
							r_ = r_<32768.0f ? igamcurve[r_] : (Color::gamman((float)r_/32768.0f, igam) * 65535.0f);
							g_ = g_<32768.0f ? igamcurve[g_] : (Color::gamman((float)g_/32768.0f, igam) * 65535.0f);
							b_ = b_<32768.0f ? igamcurve[b_] : (Color::gamman((float)b_/32768.0f, igam) * 65535.0f);
							//readapt arbitrary gamma (inverse from beginning)
							if (gamlab == 0) {
							r_ = Color::gammatab_26_11[r_];
							g_ = Color::gammatab_26_11[g_];
							b_ = Color::gammatab_26_11[b_];
							}
							else if (gamlab == 1) {
							r_ = Color::gammatab_4[r_];
							g_ = Color::gammatab_4[g_];
							b_ = Color::gammatab_4[b_];
							}
							else if (gamlab == 2) {
							r_ = Color::gammatab_55[r_];
							g_ = Color::gammatab_55[g_];
							b_ = Color::gammatab_55[b_];
							}
							float factor = Vmask[i1]*Hmask[j1]/gain;

							dsttmp->r(i,j) += factor*r_;
							dsttmp->g(i,j) += factor*g_;
							dsttmp->b(i,j) += factor*b_;

						}
					}
					}
					else {//RGB mode
					for (int i=tiletop; i<tilebottom; i++){
						int i1 = i-tiletop;
						float X,Y,Z;
						for (int j=tileleft; j<tileright; j++) {
							int j1=j-tileleft;

							Y = labdn->L[i1][j1];
							X = (labdn->a[i1][j1]) + Y;
							Z = Y - (labdn->b[i1][j1]);

							X = X<32768.0f ? igamcurve[X] : (Color::gamma((float)X/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
							Y = Y<32768.0f ? igamcurve[Y] : (Color::gamma((float)Y/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);
							Z = Z<32768.0f ? igamcurve[Z] : (Color::gamma((float)Z/32768.0f, igam, igamthresh, igamslope, 1.0, 0.0) * 65535.0f);

							float factor = Vmask[i1]*Hmask[j1]/gain;

							dsttmp->r(i,j) += factor*X;
							dsttmp->g(i,j) += factor*Y;
							dsttmp->b(i,j) += factor*Z;

						}
					}

					}
				} else {
#ifdef _OPENMP
//#pragma omp parallel for
#endif
					for (int i=tiletop; i<tilebottom; i++){
						int i1 = i-tiletop;
						float X,Y,Z,L,a,b;
						for (int j=tileleft; j<tileright; j++) {
							int j1=j-tileleft;
							//modification Jacques feb 2013
							L = labdn->L[i1][j1];
							a = labdn->a[i1][j1];
							b = labdn->b[i1][j1];
							Color::Lab2XYZ(L, a, b, X, Y, Z);

							float factor = Vmask[i1]*Hmask[j1];
							float r_,g_,b_;
							Color::xyz2rgb(X,Y,Z,r_,g_,b_,wip);
							//gamma slider is different from Raw
							r_ = r_<32768.0f ? igamcurve[r_] : (Color::gamman((float)r_/32768.0f, igam) * 65535.0f);
							g_ = g_<32768.0f ? igamcurve[g_] : (Color::gamman((float)g_/32768.0f, igam) * 65535.0f);
							b_ = b_<32768.0f ? igamcurve[b_] : (Color::gamman((float)b_/32768.0f, igam) * 65535.0f);

							dsttmp->r(i,j) += factor*r_;
							dsttmp->g(i,j) += factor*g_;
							dsttmp->b(i,j) += factor*b_;

						}
					}
				}

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				delete labdn;
			//	delete noiseh;

				delete[] Vmask;
				delete[] Hmask;



			}//end of tile row
		}//end of tile loop
#ifdef _OPENMP
#pragma omp critical
#endif
{
		fftwf_free ( Lblox);
		fftwf_free ( fLblox);
}

        }
		//copy denoised image to output
		memcpy (dst->data, dsttmp->data, 3*dst->width*dst->height*sizeof(float));

		if (!isRAW) {//restore original image gamma
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (int i=0; i<3*dst->width*dst->height; i++) {
				dst->data[i] = Color::gammatab_srgb[ dst->data[i] ];
			}
		}

		delete dsttmp;

    // destroy the plans
	fftwf_destroy_plan( plan_forward_blox[0] );
	fftwf_destroy_plan( plan_backward_blox[0] );
	fftwf_destroy_plan( plan_forward_blox[1] );
	fftwf_destroy_plan( plan_backward_blox[1] );
	fftwf_cleanup();
//#ifdef _DEBUG
//	if (settings->verbose) {
//		t2e.set();
//		printf("Denoise performed in %d usec:\n", t2e.etime(t1e));
//	}
//#endif

	}//end of main RGB_denoise


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#if defined( __SSE2__ ) && defined( WIN32 )
__attribute__((force_align_arg_pointer)) void ImProcFunctions::RGBtile_denoise (float * fLblox, int vblproc, int hblproc, int numblox_H, int numblox_W, float noisevar_Ldetail )	//for DCT
#else
	void ImProcFunctions::RGBtile_denoise (float * fLblox, int vblproc, int hblproc, int numblox_H, int numblox_W, float noisevar_Ldetail )	//for DCT
#endif
	{
		float * nbrwt  = new float[TS*TS];	//for DCT
		int blkstart = hblproc*TS*TS;

		boxabsblur(fLblox+blkstart, nbrwt, 3, 3, TS, TS);//blur neighbor weights for more robust estimation	//for DCT
		
#ifdef __SSE2__
		__m128	tempv;
		__m128	noisevar_Ldetailv = _mm_set1_ps( noisevar_Ldetail );
		__m128	onev = _mm_set1_ps( 1.0f );
		for (int n=0; n<TS*TS; n+=4) {		//for DCT
			tempv  = onev - xexpf( -SQRV( LVFU(nbrwt[n]))/noisevar_Ldetailv);
			_mm_storeu_ps( &fLblox[blkstart+n], LVFU(fLblox[blkstart+n]) * tempv );
		}//output neighbor averaged result

#else
#pragma omp parallel for
		for (int n=0; n<TS*TS; n++) {		//for DCT
			fLblox[blkstart+n] *= (1-xexpf(-SQR(nbrwt[n])/noisevar_Ldetail));
		}//output neighbor averaged result
#endif
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

		int imin = MAX(0,-top);
		int bottom = MIN( top+TS,height);
		int imax = bottom - top;

		//add row of tiles to output image
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int hblk=0; hblk < numblox_W; hblk++) {
			int left = (hblk-blkrad)*offset;
			int right  = MIN(left+TS, width);
			int jmin = MAX(0,-left);
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
													 wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_abr, float noisevar_abb, LabImage * noi)
	{
		int maxlvl = WaveletCoeffs_L.maxlevel();
		const float eps = 0.01f;
		int max;
//		float parfrac = 0.05;

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
//			float skip_h;
			float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
			float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
			float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);

			if (lvl==maxlvl-1) {
			//	ShrinkAll(WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, lvl, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,10,10,
			//			  skip_L, skip_ab, skip_h, noisevar_L, noisevar_ab, NULL, NULL);//TODO: this implies redundant evaluation of MAD
				ShrinkAll(WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, lvl, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,
						  skip_L, skip_ab, noisevar_L, noisevar_abr, noisevar_abb, noi);//TODO: this implies redundant evaluation of MAD

			} else {

//				float ** WavPars_L = WaveletCoeffs_L.level_coeffs(lvl+1);
				//float ** WavPars_a = WaveletCoeffs_a.level_coeffs(lvl+1);
				//float ** WavPars_b = WaveletCoeffs_b.level_coeffs(lvl+1);

				//simple wavelet shrinkage
				float * sfave = new float[Wlvl_L*Hlvl_L];
				array2D<float> edge(Wlvl_L,Hlvl_L);

				//printf("\n level=%d  \n",lvl);

				for (int dir=1; dir<4; dir++) {
					float mad_L = madL[lvl][dir-1];
					float mad_a = noisevar_abr*mada[lvl][dir-1];
					float mad_b = noisevar_abb*madb[lvl][dir-1];
					//float mad_Lpar = madL[lvl+1][dir-1];
					//float mad_apar = mada[lvl+1][dir-1];
					//float mad_bpar = mada[lvl+1][dir-1];

					//float skip_ab_ratio = WaveletCoeffs_a.level_stride(lvl+1)/skip_ab;
					float skip_L_ratio =  WaveletCoeffs_L.level_stride(lvl+1)/skip_L;

					if (noisevar_abr>0.01f  || noisevar_abb>0.01f) {

						//printf("  dir=%d  mad_L=%f		mad_a=%f		mad_b=%f	\n",dir,sqrt(mad_L),sqrt(mad_a),sqrt(mad_b));

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

								WavCoeffs_a[dir][coeffloc_ab] *= SQR(1.f-xexpf(-(mag_a/mad_a)-(mag_L/(9.f*mad_L)))/*satfactor_a*/);
								WavCoeffs_b[dir][coeffloc_ab] *= SQR(1.f-xexpf(-(mag_b/mad_b)-(mag_L/(9.f*mad_L)))/*satfactor_b*/);

							}
						}//now chrominance coefficients are denoised
					}

					if (noisevar_L>0.01f) {
						mad_L *= noisevar_L*5/(lvl+1);
//OpenMP here
						for (int i=0; i<Hlvl_L; i++)
							for (int j=0; j<Wlvl_L; j++) {

								int coeffloc_L = i*Wlvl_L+j;
//								int coeffloc_Lpar = (MAX(0,i-skip_L)*Wlvl_L+MAX(0,j-skip_L))/skip_L_ratio;

								float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
								//float mag_Lpar = SQR(parfrac*WavPars_L[dir][coeffloc_Lpar]);
								//float sf_L = SQR(1-expf(-(mag_L/mad_L)-(mag_Lpar/mad_L)));
								float sf_L = mag_L/(mag_L+mad_L*xexpf(-mag_L/(9.f*mad_L))+eps);

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

								float sf_L = mag_L/(mag_L + edgefactor*mad_L*xexpf(-mag_L/(9.f*mad_L))+eps);

								//use smoothed shrinkage unless local shrinkage is much less
								WavCoeffs_L[dir][coeffloc_L] *= (SQR(edgefactor*sfave[coeffloc_L])+SQR(sf_L))/(edgefactor*sfave[coeffloc_L]+sf_L+eps);

							}//now luminance coeffs are denoised

					}

				}
				delete[] sfave;

			}
		}
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	void ImProcFunctions::WaveletDenoiseAll(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a,
											wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_abr, float noisevar_abb, LabImage * noi)//mod JD

	{
		int maxlvl = WaveletCoeffs_L.maxlevel();
//       printf("maxlevel = %d\n",maxlvl);
//omp_set_nested(true);
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

       //     printf("Hab : %d\n", Hlvl_ab);
        //    printf("Wab : %d\n", Wlvl_ab);
			ShrinkAll(WavCoeffs_L, WavCoeffs_a, WavCoeffs_b, lvl, Wlvl_L, Hlvl_L, Wlvl_ab, Hlvl_ab,
					  skip_L, skip_ab, noisevar_L, noisevar_abr, noisevar_abb, noi);

		}
//omp_set_nested(false);
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if defined( __SSE2__ ) && defined( WIN32 )
__attribute__((force_align_arg_pointer))	void ImProcFunctions::ShrinkAll(float ** WavCoeffs_L, float ** WavCoeffs_a, float ** WavCoeffs_b, int level,
									int W_L, int H_L, int W_ab, int H_ab,int skip_L, int skip_ab, float noisevar_L, float noisevar_abr,  float noisevar_abb, LabImage * noi)
#else
	void ImProcFunctions::ShrinkAll(float ** WavCoeffs_L, float ** WavCoeffs_a, float ** WavCoeffs_b, int level,
									int W_L, int H_L, int W_ab, int H_ab,int skip_L, int skip_ab, float noisevar_L, float noisevar_abr, float noisevar_abb, LabImage * noi)
#endif

									{
		//simple wavelet shrinkage
		const float eps = 0.01f;
		float * sfave = new float[W_L*H_L];
		float * sfavea = new float[W_L*H_L];
		float * sfaveb = new float[W_L*H_L];

		int max;

		//printf("\n level=%d  \n",level);

		for (int dir=1; dir<4; dir++) {
			float madL = SQR(MadMax(WavCoeffs_L[dir], max, W_L*H_L));
			float mada = SQR(MadMax(WavCoeffs_a[dir], max, W_ab*H_ab));
			float madb = SQR(MadMax(WavCoeffs_b[dir], max, W_ab*H_ab));


		//	printf("  dir=%d  mad_L=%f	mad_a=%f	mad_b=%f  skip_ab=%i	\n",dir,sqrt(madL),sqrt(mada),sqrt(madb), skip_ab);

			float mad_L = madL*noisevar_L*5/(level+1);
			float mad_a = mada*noisevar_abr;  // noisevar_abr between 0..2.25=default 100=middle value  ==> 582=max
			float mad_b = madb*noisevar_abb;

			if (noisevar_abr>0.01f  ||  noisevar_abb>0.01f) {
//OpenMP here

#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (int i=0; i<H_ab; i++) {
					for (int j=0; j<W_ab; j++) {
						float m_a,m_b;
						m_a=mad_a;
						m_b=mad_b;
						int coeffloc_ab = i*W_ab+j;
						int coeffloc_L	= ((i*skip_L)/skip_ab)*W_L + ((j*skip_L)/skip_ab);
						float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L ])+eps;
						float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
						float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;
						sfavea[coeffloc_ab] = (1.f-xexpf(-(mag_a/mad_a)-(mag_L/(9.f*madL))));
						sfaveb[coeffloc_ab] = (1.f-xexpf(-(mag_b/mad_b)-(mag_L/(9.f*madL))));
						mad_a=m_a;
						mad_b=m_b;
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
						float m_a,m_b;
						m_a=mad_a;
						m_b=mad_b;

						int coeffloc_ab = i*W_ab+j;
						int coeffloc_L	= ((i*skip_L)/skip_ab)*W_L + ((j*skip_L)/skip_ab);
						float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L ])+eps;
						float mag_a = SQR(WavCoeffs_a[dir][coeffloc_ab])+eps;
						float mag_b = SQR(WavCoeffs_b[dir][coeffloc_ab])+eps;

						float sfa = (1.f-xexpf(-(mag_a/mad_a)-(mag_L/(9.f*madL))));
						float sfb = (1.f-xexpf(-(mag_b/mad_b)-(mag_L/(9.f*madL))));

						//use smoothed shrinkage unless local shrinkage is much less
						WavCoeffs_a[dir][coeffloc_ab] *= (SQR(sfavea[coeffloc_ab])+SQR(sfa))/(sfavea[coeffloc_ab]+sfa+eps);
						WavCoeffs_b[dir][coeffloc_ab] *= (SQR(sfaveb[coeffloc_ab])+SQR(sfb))/(sfaveb[coeffloc_ab]+sfb+eps);
						mad_a=m_a;
						mad_b=m_b;

					}//now chrominance coefficients are denoised
			}

			if (noisevar_L>0.01f) {
#ifdef __SSE2__
				__m128	magv;
				__m128  mad_Lv = _mm_set1_ps( mad_L );
				__m128	ninev = _mm_set1_ps( 9.0f );
				__m128 	epsv = _mm_set1_ps( eps );
				for (int i=0; i<W_L*H_L-3; i+=4) {
					magv = SQRV(LVFU(WavCoeffs_L[dir][i]));
					_mm_storeu_ps( &sfave[i], magv / (magv + mad_Lv*xexpf(-magv/(ninev * mad_Lv)) + epsv));
				}
				for (int i=(W_L*H_L)-((W_L*H_L)%4); i<W_L*H_L; i++) {
					float mag = SQR(WavCoeffs_L[dir][i]);
					sfave[i] = mag/(mag+mad_L*xexpf(-mag/(9*mad_L))+eps);
				}
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (int i=0; i<W_L*H_L; i++) {

					float mag = SQR(WavCoeffs_L[dir][i]);
					float shrinkfactor = mag/(mag+mad_L*xexpf(-mag/(9*mad_L))+eps);

					//WavCoeffs_L[dir][i] *= shrinkfactor;
					sfave[i] = shrinkfactor;
				}
#endif
				boxblur(sfave, sfave, level+2, level+2, W_L, H_L);//increase smoothness by locally averaging shrinkage
#ifdef __SSE2__
				__m128	sfv;
				for (int i=0; i<W_L*H_L-3; i+=4) {
					magv = SQRV( LVFU(WavCoeffs_L[dir][i]));
					sfv = magv/(magv + mad_Lv * xexpf( -magv / (ninev * mad_Lv)) + epsv );
					//use smoothed shrinkage unless local shrinkage is much less
					_mm_storeu_ps( &WavCoeffs_L[dir][i], _mm_loadu_ps( &WavCoeffs_L[dir][i]) * (SQRV( LVFU(sfave[i] )) + SQRV(sfv)) / (LVFU(sfave[i])+sfv+epsv));
				}
				for (int i=(W_L*H_L)-((W_L*H_L)%4); i<W_L*H_L; i++) {
					float mag = SQR(WavCoeffs_L[dir][i]);
					float sf = mag/(mag+mad_L*xexpf(-mag/(9*mad_L))+eps);

					//use smoothed shrinkage unless local shrinkage is much less
					WavCoeffs_L[dir][i] *= (SQR(sfave[i])+SQR(sf))/(sfave[i]+sf+eps);
				}//now luminance coefficients are denoised

#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (int i=0; i<W_L*H_L; i++) {


					float mag = SQR(WavCoeffs_L[dir][i]);
					float sf = mag/(mag+mad_L*xexpf(-mag/(9.f*mad_L))+eps);

					//use smoothed shrinkage unless local shrinkage is much less
					WavCoeffs_L[dir][i] *= (SQR(sfave[i])+SQR(sf))/(sfave[i]+sf+eps);

				}//now luminance coefficients are denoised
#endif
			}


		}
		delete[] sfave;
		delete[] sfavea;
		delete[] sfaveb;

	}



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


}

