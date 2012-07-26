/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#ifndef _IMPROCFUN_H_
#define _IMPROCFUN_H_

#include "imagefloat.h"
#include "image16.h"
#include "image8.h"
#include "procparams.h"
#include "shmap.h"
#include "coord2d.h"
#include "color.h"
#include "labimage.h"
#include "LUT.h"
#include "lcp.h"
#include <fftw3.h>
#include "cplx_wavelet_dec.h"

namespace rtengine {

using namespace procparams;

class ImProcFunctions {

		static LUTf gamma2curve;

		cmsHTRANSFORM monitorTransform;

		const ProcParams* params;
		double scale;
		bool multiThread;
		float g;

        void calcVignettingParams(int oW, int oH, const VignettingParams& vignetting, double &w2, double &h2, double& maxRadius, double &v, double &b, double &mul);

		void transformPreview      (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, const LCPMapper *pLCPMap);
		void transformVignetteOnly (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int oW, int oH);
		void transformHighQuality  (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, const LCPMapper *pLCPMap, bool fullImage);

		void sharpenHaloCtrl    (LabImage* lab, float** blurmap, float** base, int W, int H);
		void firstAnalysisThread(Imagefloat* original, Glib::ustring wprofile, unsigned int* histogram, int row_from, int row_to);
		void dcdamping          (float** aI, float** aO, float damping, int W, int H);

		bool needsCA            ();
		bool needsDistortion    ();
		bool needsRotation      ();
		bool needsPerspective   ();
		bool needsVignetting    ();
        bool needsLCP           ();
 //   static cmsUInt8Number* Mempro = NULL;		


	public:
		static LUTf cachef;
		
		double lumimul[3];

		static void initCache ();
		static void cleanupCache ();
		
		ImProcFunctions       (const ProcParams* iparams, bool imultiThread=true)
			: monitorTransform(NULL), params(iparams), scale(1), multiThread(imultiThread) {}
		~ImProcFunctions      ();

		void setScale         (double iscale);

		bool needsTransform   ();

		void firstAnalysis    (Imagefloat* working, const ProcParams* params, LUTu & vhist16, double gamma);
		void rgbProc          (Imagefloat* working, LabImage* lab, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
		                       SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve);
		void rgbProc          (Imagefloat* working, LabImage* lab, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
		                       SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve,
		                       double expcomp, int hlcompr, int hlcomprthresh);
		void luminanceCurve   (LabImage* lold, LabImage* lnew, LUTf &curve);
		void chromiLuminanceCurve (LabImage* lold, LabImage* lnew, LUTf &acurve, LUTf &bcurve, LUTf & satcurve/*,LUTf & satbgcurve*/, LUTf &curve, bool utili, bool autili, bool butili, bool ccutili);
		void vibrance 		  (LabImage* lab);//Jacques' vibrance
		void colorCurve       (LabImage* lold, LabImage* lnew);
		void sharpening       (LabImage* lab, float** buffer);
		void transform        (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH,
                               double focalLen, double focalLen35mm, float focusDist, int rawRotationDeg, bool fullImage);
		void lab2monitorRgb   (LabImage* lab, Image8* image);
		void resize           (Image16* src, Image16* dst, float dScale);
		void deconvsharpening (LabImage* lab, float** buffer);
		void MLsharpen (LabImage* lab);// Manuel's clarity / sharpening
		void MLmicrocontrast(LabImage* lab ); //Manuel's microcontrast

		void impulsedenoise   (LabImage* lab);//Emil's impulse denoise
		void impulse_nr (LabImage* lab, double thresh);
		void dirpyrdenoise    (LabImage* src);//Emil's pyramid denoise
		void dirpyrequalizer  (LabImage* lab);//Emil's equalizer
		
		void EPDToneMap(LabImage *lab, unsigned int Iterates = 0, int skip = 1);
		
		// pyramid denoise
		procparams::DirPyrDenoiseParams dnparams;
		void dirpyrLab_denoise(LabImage * src, LabImage * dst, const procparams::DirPyrDenoiseParams & dnparams );//Emil's directional pyramid denoise
		void dirpyr           (LabImage* data_fine, LabImage* data_coarse, int level, LUTf &rangefn_L, LUTf &rangefn_ab,
							   int pitch, int scale, const int luma, int chroma );
		void idirpyr          (LabImage* data_coarse, LabImage* data_fine, int level, LUTf &rangefn_L, LUTf & nrwt_l, LUTf & nrwt_ab,
							   int pitch, int scale, const int luma, const int chroma/*, LUTf & Lcurve, LUTf & abcurve*/ );
		
		// FT denoise
		//void RGB_InputTransf(Imagefloat * src, LabImage * dst, const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe);
		//void RGB_OutputTransf(LabImage * src, Imagefloat * dst, const procparams::DirPyrDenoiseParams & dnparams);
		//void output_tile_row (float *Lbloxrow, float ** Lhipassdn, float ** tilemask, int height, int width, int top, int blkrad );
		void RGB_denoise(Imagefloat * src, Imagefloat * dst, const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe);
		void RGBtile_denoise (float * fLblox, int vblproc, int hblproc, int numblox_H, int numblox_W, float noisevar_L );	//for DCT
		void RGBoutput_tile_row (float *Lbloxrow, float ** Ldetail, float ** tilemask_out, int height, int width, int top );
		//void WaveletDenoise(cplx_wavelet_decomposition &DualTreeCoeffs, float noisevar );
		//void WaveletDenoise(wavelet_decomposition &WaveletCoeffs, float noisevar );
		void WaveletDenoiseAll(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a, 
							   wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_ab );
		void WaveletDenoiseAll_BiShrink(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_a, 
										wavelet_decomposition &WaveletCoeffs_b, float noisevar_L, float noisevar_ab );
		//void BiShrink(float * ReCoeffs, float * ImCoeffs, float * ReParents, float * ImParents, 
		//			  int W, int H, int level, int padding, float noisevar);
		//void Shrink(float ** WavCoeffs, int W, int H, int level, float noisevar);
		void ShrinkAll(float ** WavCoeffs_L, float ** WavCoeffs_a, float ** WavCoeffs_b, int level, 
					   int W_L, int H_L, int W_ab, int H_ab, int skip_L, int skip_ab, float noisevar_L, float noisevar_ab);
		float MadMax(float * HH_Coeffs, int &max, int datalen);
		
		// pyramid equalizer
		void dirpyr_equalizer    (float ** src, float ** dst, int srcwidth, int srcheight, const double * mult );//Emil's directional pyramid equalizer
		void dirpyr_channel      (float ** data_fine, float ** data_coarse, int width, int height, LUTf & rangefn, int level, int scale, const double * mult  );
		void idirpyr_eq_channel  (float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, const double * mult );

		void defringe         (LabImage* lab);
		void PF_correct_RT    (LabImage * src, LabImage * dst, double radius, int thresh);

		Image8*     lab2rgb   (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile);
		Image16*    lab2rgb16b (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile, Glib::ustring profi, Glib::ustring gam, bool freegamma, double gampos, double slpos, double &ga0, double &ga1, double &ga2, double &ga3, double &ga4, double &ga5, double &ga6);// for gamma output		
		Image16*    lab2rgb16 (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile);//without gamma ==>default

		bool transCoord       (int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef = -1, const LCPMapper *pLCPMap=NULL);
		bool transCoord       (int W, int H, const std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef = -1, const LCPMapper *pLCPMap=NULL);
		void getAutoExp       (LUTu & histogram, int histcompr, double defgain, double clip, double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh);
		static double getAutoDistor  (const Glib::ustring& fname, int thumb_size);	
		double getTransformAutoFill (int oW, int oH, const LCPMapper *pLCPMap=NULL);
};
}
#endif
