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

#ifndef _COLOR_H_
#define _COLOR_H_

#include "rt_math.h"
#include "LUT.h"
#include "labimage.h"
#include "iccstore.h"
#include "iccmatrices.h"

namespace rtengine {

#ifdef _DEBUG

class MunsellDebugInfo {
public:
	float maxdhuelum[4];
	float maxdhue[4];
	unsigned int depass;
	unsigned int depassLum;

	MunsellDebugInfo();
	void reinitValues();
};

#endif

class Color {

private:
	// Jacques' 195 LUTf for Munsell Lch correction
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

	// Separated from init() to keep the code clear
	static void initMunsell ();
    static double hue2rgb(double p, double q, double t);

public:
	const static double sRGBGamma;        // standard average gamma
	const static double sRGBGammaCurve;   // 2.4 in the curve
	const static double eps_max, kappa, epskap;
	const static float D50x, D50z;
	const static double u0, v0;

	static cmsToneCurve* linearGammaTRC;

	static LUTf cachef;
	static LUTf gamma2curve;

	// look-up tables for the standard srgb gamma and its inverse (filled by init())
	static LUTf igammatab_srgb;
	static LUTf gammatab_srgb;
//	static LUTf igammatab_709;
//	static LUTf gammatab_709;
	static LUTf igammatab_26_11;
	static LUTf gammatab_26_11;
	static LUTf igammatab_24_17;
	static LUTf gammatab_24_17a;
	
	// look-up tables for the simple exponential gamma
	static LUTf gammatab;


	static void init ();
	static void cleanup ();

	static float rgbLuminance(float r, float g, float b) { return r*float(xyz_sRGBd65[1][0]) + g*float(xyz_sRGBd65[1][1]) + b*float(xyz_sRGBd65[1][2]); }
	static double rgbLuminance(double r, double g, double b) { return r*xyz_sRGBd65[1][0] + g*xyz_sRGBd65[1][1] + b*xyz_sRGBd65[1][2]; }
	static void rgb2hsl (float r, float g, float b, float &h, float &s, float &l);
	static void hsl2rgb (float h, float s, float l, float &r, float &g, float &b);
	static void rgb2hsv (float r, float g, float b, float &h, float &s, float &v);
	static void hsv2rgb (float h, float s, float v, float &r, float &g, float &b);
	static void hsv2rgb (float h, float s, float v, int &r, int &g, int &b);
	static void hsv2rgb01 (float h, float s, float v, float &r, float &g, float &b);
	static void xyz2srgb (float x, float y, float z, float &r, float &g, float &b);
	static void xyz2Prophoto (float x, float y, float z, float &r, float &g, float &b);
	static void Prophotoxyz (float r, float g, float b, float &x, float &y, float &z);
	static void xyz2rgb (float x, float y, float z, float &r, float &g, float &b, double rgb_xyz[3][3]);
	static void xyz2rgb (float x, float y, float z, float &r, float &g, float &b, float rgb_xyz[3][3]);
	static void rgbxyz (float r, float g, float b, float &x, float &y, float &z, double xyz_rgb[3][3]);

	static void Lab2XYZ(float L, float a, float b, float &x, float &y, float &z);
	static void XYZ2Lab(float X, float Y, float Z, float &L, float &a, float &b);
	static void Lab2Yuv(float L, float a, float b, float &Y, float &u, float &v);
	static void Yuv2Lab(float Y, float u, float v, float &L, float &a, float &b, double wp[3][3]);
	static double f2xyz(double f);
	static void calcGamma (double pwr, double ts, int mode, int imax, double &gamma0, double &gamma1, double &gamma2, double &gamma3, double &gamma4,double &gamma5);
	static void trcGammaBW (float &r, float &g, float &b, float gammabwr, float gammabwg, float gammabwb);
	static void computeBWMixerConstants (const Glib::ustring &setting, const Glib::ustring &filter, const Glib::ustring &algo, float &mixerRed, float &mixerGreen,
										float &mixerBlue, float mixerOrange, float mixerYellow, float mixerCyan, float mixerPurple, float mixerMagenta,
										bool autoc, bool complement, float &kcorec, double &rrm, double &ggm, double &bbm);


	// standard srgb gamma and its inverse
	static inline double gamma2     (double x) {
											return x <= 0.003041 ? x*12.92 : 1.055011*exp(log(x)/sRGBGammaCurve)-0.055011;
									}
	static inline double igamma2    (double x) {
										return x <= 0.039293 ? x/12.92 : exp(log((x+0.055011)/1.055011)*sRGBGammaCurve);
									}
/*	static inline double gamma709     (double x) {
											return x <= 0.0176 ? x*4.5 : 1.0954*exp(log(x)/2.2)-0.0954;
									}
	static inline double igamma709    (double x) {
										return x <= 0.0795 ? x/4.5 : exp(log((x+0.0954)/1.0954)*2.2);
									}	
*/	
	static inline double gamma24_17     (double x) {
											return x <= 0.001867 ? x*17.0 : 1.044445*exp(log(x)/2.4)-0.044445;
									}
	static inline double igamma24_17    (double x) {
										return x <= 0.031746 ? x/17.0 : exp(log((x+0.044445)/1.044445)*2.4);
									}	

	static inline double gamma26_11     (double x) {
											return x <= 0.004921 ? x*11.0 : 1.086603*exp(log(x)/2.6)-0.086603;
									}
	static inline double igamma26_11    (double x) {
										return x <= 0.054127 ? x/11.0 : exp(log((x+0.086603)/1.086603)*2.6);
									}	

	// gamma function with adjustable parameters
	static inline double gamma      (double x, double gamma, double start, double slope, double mul, double add){
										return (x <= start ? x*slope : exp(log(x)/gamma)*mul-add);
									}
	static inline double igamma     (double x, double gamma, double start, double slope, double mul, double add){
										return (x <= start*slope ? x/slope : exp(log((x+add)/mul)*gamma) );
									}
	static inline double gamman      (double x, double gamma){//gamma standard without slope...
										return (x =exp(log(x)/gamma));
									}
	static inline double igamman     (double x, double gamma){//inverse gamma standard without slope...
										return (x = exp(log(x)*gamma) );
									}

	// gamma functions on [0,65535] based on look-up tables
	static inline float  gamma_srgb       (char x) { return gammatab_srgb[x]; }
	static inline float  gamma            (char x) { return gammatab[x]; }
	static inline float  igamma_srgb      (char x) { return igammatab_srgb[x]; }
	static inline float  gamma_srgb       (int x) { return gammatab_srgb[x]; }
	static inline float  gamma            (int x) { return gammatab[x]; }
	static inline float  igamma_srgb      (int x) { return igammatab_srgb[x]; }
	static inline float  gamma_srgb       (float x) { return gammatab_srgb[x]; }
	static inline float  gamma            (float x) { return gammatab[x]; }
	static inline float  igamma_srgb      (float x) { return igammatab_srgb[x]; }
	//static inline float  gamma_srgb       (double x) { return gammatab_srgb[x]; }
	//static inline float  gamma            (double x) { return gammatab[x]; }
	//static inline float  igamma_srgb      (double x) { return igammatab_srgb[x]; }

    //Jacques's Munsell correction
#ifdef _DEBUG
	static void AllMunsellLch (bool lumaMuns, float Lprov1, float Loldd, float HH, float Chprov1, float CC, float &correctionHueChroma, float &correctlum, MunsellDebugInfo* munsDbgInfo);
	static void gamutLchonly  (float HH, float &Lprov1, float &Chprov1, float &R, float &G, float &B, double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef, bool &neg, bool &more_rgb);
#else
	static void AllMunsellLch (bool lumaMuns, float Lprov1, float Loldd, float HH, float Chprov1, float CC, float &correctionHueChroma, float &correctlum);
	static void gamutLchonly  (float HH, float &Lprov1, float &Chprov1, float &R, float &G, float &B, double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef);
#endif
	static void LabGamutMunsell (LabImage *lab, float *Lold, float *Cold, bool corMunsell, bool lumaMuns, bool isHLEnabled, bool gamut, const Glib::ustring &working, bool multiThread );

	static void SkinSat         (float lum, float hue, float chrom, float &satreduc, int chromx);//jacques Skin color
	static void MunsellLch      (float lum, float hue, float chrom, float memChprov, float &correction, int zone, float &lbe, bool &correctL);//jacques:  Munsell correction
	// end Munsell
	static void scalered ( float rstprotection, float param, float limit, float HH, float deltaHH, float &scale, float &scaleext);
	static void transitred (float HH, float Chprov1, float dred, float factorskin, float protect_red, float factorskinext, float deltaHH, float factorsat, float &factor);
	static void skinred ( double J, double h, double sres, double Sp, float dred, float protect_red, int sk, float rstprotection, float ko, double &s);
	static void skinredfloat ( float J, float h, float sres, float Sp, float dred, float protect_red, int sk, float rstprotection, float ko, float &s);

	//void gamutmap(LabImage* );
	static void gamutmap(float &X, float &Y, float &Z, const double p[3][3]);

	static inline double huelab_to_huehsv2 (float HH){
					//hr=translate Hue Lab value  (-Pi +Pi) in approximative hr (hsv values) (0 1) [red 1/6 yellow 1/6 green 1/6 cyan 1/6 blue 1/6 magenta 1/6 ]
				// with multi linear correspondances (I expect there is no error !!)
				double hr;
				//allways put h between 0 and 1

				if      (HH>=0.f && HH < 0.6f) 	 hr=0.11666*(double) HH + 0.93;  //hr 0.93 1. full red
				else if (HH>=0.6f && HH < 1.4f)	 hr=0.1125*double(HH) - 0.0675;   //hr 0.0  0.09      red yellow orange          
				else if (HH>=1.4f && HH < 2.f)	 hr=0.2666*double(HH) - 0.2833;   //hr 0.09  0.25    orange yellow             
				else if (HH>=2.f && HH < 3.14159f) hr=0.1489*double(HH) - 0.04785;  //hr 0.25 0.42  yellow green green
				else if (HH>=-3.14159f && HH < -2.8f) hr=0.23419*double(HH) +1.1557;  //hr 0.42 0.5  green
				else if (HH>=-2.8f && HH < -2.3f) hr=0.16*double(HH) + 0.948;    //hr 0.5 0.58      cyan         
				else if (HH>=-2.3f && HH < -0.9f) hr=0.12143*double(HH)+ 0.85928;    //hr 0.58 0.75   blue blue-sky            
				else if (HH>=-0.9f && HH < -0.1f) hr=0.2125*double(HH) + 0.94125;        //hr 0.75  0.92    purple magenta          
				else if (HH>=-0.1f && HH < 0.f)   hr=0.1*double(HH) + 0.93;        //hr 0.92  0.93    red          
				// in case of !
				if     (hr<0.0) hr += 1.0;
				else if(hr>1.0) hr -= 1.0;
				return (hr);
	}
	
	
	static inline float f2xyz(float f) {
		const float epsilonExpInv3 = 6.0/29.0;
		const float kappaInv = 27.0/24389.0;  // inverse of kappa

		return (f > epsilonExpInv3) ? f*f*f : (116 * f - 16) * kappaInv;
	}

};

}

#endif
