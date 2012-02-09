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
#include <fftw3.h>
#include "cplx_wavelet_dec.h"

namespace rtengine {
	
	using namespace procparams;
	
	class ImProcFunctions : public Color {
		
		cmsHTRANSFORM monitorTransform;
		
		const ProcParams* params;
		double scale;
		bool multiThread;
		float g;
		
		void simpltransform     (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH);
		void vignetting         (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int oW, int oH);
		void transformNonSep    (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH);
		void transformSep       (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH);
		void sharpenHaloCtrl    (LabImage* lab, float** blurmap, float** base, int W, int H);
		void firstAnalysisThread(Imagefloat* original, Glib::ustring wprofile, unsigned int* histogram, int row_from, int row_to);
		void dcdamping          (float** aI, float** aO, float damping, int W, int H);
		
		bool needsCA            ();
		bool needsDistortion    ();
		bool needsRotation      ();
		bool needsPerspective   ();
		bool needsVignetting    ();
		//   static cmsUInt8Number* Mempro = NULL;		
		
		
	public:
		// 195 LUTf for Munsell Lch correction
		static LUTf _4P10,_4P20,_4P30,_4P40,_4P50,_4P60;
		static LUTf _1P10,_1P20,_1P30,_1P40,_1P50,_1P60;
		static LUTf _5B40,_5B50,_5B60, _5B70,_5B80;
		static LUTf _7B40,_7B50,_7B60, _7B70,_7B80;
		static LUTf _9B40,_9B50,_9B60, _9B70,_9B80;
		static LUTf _10B40,_10B50,_10B60, _10B70,_10B80;
		static LUTf _05PB40,_05PB50,_05PB60, _05PB70,_05PB80;
		static LUTf _10PB10,_10PB20,_10PB30,_10PB40,_10PB50,_10PB60;
		static LUTf _9PB10,_9PB20,_9PB30,_9PB40,_9PB50,_9PB60,_9PB70,_9PB80;
		static LUTf _75PB10,_75PB20,_75PB30,_75PB40,_75PB50,_75PB60,_75PB70,_75PB80;
		static LUTf _6PB10,_6PB20,_6PB30,_6PB40,_6PB50,_6PB60,_6PB70,_6PB80;
		static LUTf _45PB10,_45PB20,_45PB30,_45PB40,_45PB50,_45PB60,_45PB70,_45PB80;
		static LUTf _3PB10,_3PB20,_3PB30,_3PB40,_3PB50,_3PB60,_3PB70,_3PB80;
		static LUTf _15PB10,_15PB20,_15PB30,_15PB40,_15PB50,_15PB60, _15PB70,_15PB80;
		static LUTf _10YR20, _10YR30, _10YR40,_10YR50,_10YR60,_10YR70,_10YR80,_10YR90;
		static LUTf _85YR20, _85YR30, _85YR40,_85YR50,_85YR60,_85YR70,_85YR80,_85YR90;
		static LUTf  _7YR30, _7YR40,_7YR50,_7YR60,_7YR70,_7YR80;
		static LUTf  _55YR30, _55YR40,_55YR50,_55YR60,_55YR70,_55YR80,_55YR90;
		static LUTf  _4YR30, _4YR40,_4YR50,_4YR60,_4YR70,_4YR80;
		static LUTf  _25YR30, _25YR40,_25YR50,_25YR60,_25YR70;
		static LUTf  _10R30, _10R40,_10R50,_10R60,_10R70;
		static LUTf  _9R30, _9R40,_9R50,_9R60,_9R70;
		static LUTf  _7R30, _7R40,_7R50,_7R60,_7R70;
		static LUTf  _5R10, _5R20,_5R30;
		static LUTf  _25R10, _25R20,_25R30;
		static LUTf  _10RP10, _10RP20,_10RP30;
		static LUTf  _7G30, _7G40,_7G50,_7G60,_7G70,_7G80;
		static LUTf  _5G30, _5G40,_5G50,_5G60,_5G70,_5G80;
		static LUTf  _25G30, _25G40,_25G50,_25G60,_25G70,_25G80;
		static LUTf  _1G30, _1G40,_1G50,_1G60,_1G70,_1G80;
		static LUTf  _10GY30, _10GY40,_10GY50,_10GY60,_10GY70,_10GY80;
		static LUTf  _75GY30, _75GY40,_75GY50,_75GY60,_75GY70,_75GY80;
		static LUTf  _5GY30, _5GY40,_5GY50,_5GY60,_5GY70,_5GY80;
		
		double lumimul[3];
		
		static void initMunsell ();
		
		ImProcFunctions       (const ProcParams* iparams, bool imultiThread=true)
		: monitorTransform(NULL), params(iparams), scale(1), multiThread(imultiThread) {}
		~ImProcFunctions      ();
		
		void setScale         (double iscale);
		
		bool needsTransform   ();
		
		void firstAnalysis    (Imagefloat* working, const ProcParams* params, LUTu & vhist16, double gamma);
		void rgbProc          (Imagefloat* working, LabImage* lab, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve, \
							   SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve);
		void luminanceCurve   (LabImage* lold, LabImage* lnew, LUTf &curve);
		void chrominanceCurve (LabImage* lold, LabImage* lnew, LUTf &acurve, LUTf &bcurve, LUTf & satcurve);
		void vibrance         (LabImage* lab);//Jacques' vibrance
		void skinsat          (float lum, float hue, float chrom, float &satreduc);//jacques Skin color
		void MunsellLch       (float lum, float hue, float chrom, float memChprov, float &correction, int zone);//jacques:  Munsell correction
		void colorCurve       (LabImage* lold, LabImage* lnew);
		void sharpening       (LabImage* lab, float** buffer);
		void transform        (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH);
		void lab2rgb          (LabImage* lab, Image8* image);
		void resize           (Image16* src, Image16* dst, float dScale);
		void resize			  (LabImage* src, LabImage* dst, float dScale);
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
		void dirpyr           (LabImage* data_fine, LabImage* data_coarse, int level, LUTf &rangefn_L, LUTf &rangefn_ab, \
							   int pitch, int scale, const int luma, int chroma );
		void idirpyr          (LabImage* data_coarse, LabImage* data_fine, int level, LUTf &rangefn_L, LUTf & nrwt_l, LUTf & nrwt_ab, \
							   int pitch, int scale, const int luma, const int chroma/*, LUTf & Lcurve, LUTf & abcurve*/ );
		
		// FT denoise
		//void L_denoise(Imagefloat * src, LabImage * dst, const procparams::DirPyrDenoiseParams & dnparams);//Emil's FT denoise
		//void tile_denoise (fftwf_complex ** fLblox, int vblproc, int hblproc, \
		int blkrad, int numblox_H, int numblox_W, float noisevar );
		void RGB_InputTransf(Imagefloat * src, LabImage * dst, const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe);
		void RGB_OutputTransf(LabImage * src, Imagefloat * dst, const procparams::DirPyrDenoiseParams & dnparams);
		void output_tile_row (float *Lbloxrow, float ** Lhipassdn, float ** tilemask, int height, int width, int top, int blkrad );
		void RGB_denoise(Imagefloat * src, Imagefloat * dst, const procparams::DirPyrDenoiseParams & dnparams, const procparams::DefringeParams & defringe);
		void RGBtile_denoise (fftwf_complex ** fLblox, fftwf_complex ** fRLblox, fftwf_complex ** fBLblox, \
							  int vblproc, int hblproc, int blkrad, int numblox_H, int numblox_W, float noisevar_L, float noisevar_ab );	//for FT
		//void RGBtile_denoise (float ** fLblox, fftwf_complex ** fRLblox, fftwf_complex ** fBLblox, \
							  int vblproc, int hblproc, int blkrad, int numblox_H, int numblox_W, float noisevar_L, float noisevar_ab );	//for DCT
		void RGBoutput_tile_row (float *Lbloxrow, float *RLbloxrow, float *BLbloxrow, LabImage * labdn, \
								 float ** tilemask_out, int height, int width, int top, int blkrad );
		void dirpyr_ab(LabImage * data_fine, LabImage * data_coarse, const procparams::DirPyrDenoiseParams & dnparams);
		void dirpyr_ablevel(LabImage * data_fine, LabImage * data_coarse, int width, int height, \
							LUTf &rangefn_L, LUTf &rangefn_ab, int level, int scale);
		void ImStats(float* src, float* dst, int H, int W, int box );
		void WaveletDenoise(cplx_wavelet_decomposition &DualTreeCoeffs, float noisevar );
		void WaveletDenoise(wavelet_decomposition &WaveletCoeffs, float noisevar );
		void BiShrink(float * ReCoeffs, float * ImCoeffs, float * ReParents, float * ImParents, \
					  int W, int H, int level, int padding, float noisevar);
		void BiShrink(float ** WavCoeffs, float ** WavParents, int W, int H, int level, int padding, float noisevar);
		void FirstStageWiener(float* ReCoeffs, float* ImCoeffs, float* wiener1, int W, int H, int rad, float noisevar);
		void SecondStageWiener(float* ReCoeffs, float* ImCoeffs, float* wiener1, int W, int H, int rad, float noisevar);
		void QCoeffs (float* srcre, float* srcim, float* wiener1, float* dst, int rad, int W, int H);
		float UniversalThresh(float * HH_Coeffs, int datalen);

		
		// pyramid equalizer
		void dirpyr_equalizer    (float ** src, float ** dst, int srcwidth, int srcheight, const double * mult );//Emil's directional pyramid equalizer
		void dirpyr_channel      (float ** data_fine, float ** data_coarse, int width, int height, LUTf & rangefn, int level, int scale, const double * mult  );
		void idirpyr_eq_channel  (float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, const double * mult );
		
		void defringe         (LabImage* lab);
		void PF_correct_RT    (LabImage * src, LabImage * dst, double radius, int thresh);
		
		Image8*     lab2rgb   (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile);
		Image16*    lab2rgb16b (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile, Glib::ustring profi, Glib::ustring gam, bool freegamma, double gampos, double slpos, double &ga0, double &ga1, double &ga2, double &ga3, double &ga4, double &ga5, double &ga6);// for gamma output		
		Image16*    lab2rgb16 (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile);//without gamma ==>default
		
		bool transCoord       (int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef = -1);
		bool transCoord       (int W, int H, std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef = -1);
		void getAutoExp       (LUTu & histogram, int histcompr, double defgain, double clip, double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh);
		static double getAutoDistor  (const Glib::ustring& fname, int thumb_size);	
		double getTransformAutoFill (int oW, int oH);
	};
}
#endif
