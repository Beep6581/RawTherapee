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

#include <math.h>
#include "LUT.h"

namespace rtengine {

class Color {

public:
	const static double sRGBGamma;        // standard average gamma
    const static double sRGBGammaCurve;   // 2.4 in the curve
    const static double eps_max, kappa;
    const static float D50x, D50z;
    const static double u0, v0;

    static LUTf cachef;
	static LUTf gamma2curve;

    // look-up tables for the standard srgb gamma and its inverse (filled by init())
    static LUTf igammatab_srgb;
    static LUTf gammatab_srgb;
    // look-up tables for the simple exponential gamma
    static LUTf gammatab;


	static void init ();
	static void cleanup ();

	static void rgb2hsv (float r, float g, float b, float &h, float &s, float &v);
	static void hsv2rgb (float h, float s, float v, float &r, float &g, float &b);
	static void hsv2rgb (float h, float s, float v, int &r, int &g, int &b);
    static void hsv2rgb01 (float h, float s, float v, float &r, float &g, float &b);
	static void xyz2srgb (float x, float y, float z, float &r, float &g, float &b);
	static void xyz2rgb (float x, float y, float z, float &r, float &g, float &b, double rgb_xyz[3][3]);
    static void xyz2rgb (float x, float y, float z, float &r, float &g, float &b, float rgb_xyz[3][3]);
	static void Lab2XYZ(float L, float a, float b, float &x, float &y, float &z);
	static void XYZ2Lab(float X, float Y, float Z, float &L, float &a, float &b);
	static void Lab2Yuv(float L, float a, float b, float &Y, float &u, float &v);
	static void Yuv2Lab(float Y, float u, float v, float &L, float &a, float &b, double wp[3][3]);
    static double f2xyz(double f);
	static void calcGamma (double pwr, double ts, int mode, int imax, double &gamma0, double &gamma1, double &gamma2, double &gamma3, double &gamma4,double &gamma5);

	// standard srgb gamma and its inverse
	static inline double gamma2     (double x) {
											return x <= 0.00304 ? x*12.92 : 1.055*exp(log(x)/sRGBGammaCurve)-0.055;
									}
	static inline double igamma2    (double x) {
										return x <= 0.03928 ? x/12.92 : exp(log((x+0.055)/1.055)*sRGBGammaCurve);
									}
	// gamma function with adjustable parameters
	static inline double gamma      (double x, double gamma, double start, double slope, double mul, double add){
										return (x <= start ? x*slope : exp(log(x)/gamma)*mul-add);
									}
	static inline double igamma     (double x, double gamma, double start, double slope, double mul, double add){
										return (x <= start*slope ? x/slope : exp(log((x+add)/mul)*gamma) );
									}

	// gamma functions on [0,65535] based on look-up tables
	static inline float  gamma_srgb       (int x) { return gammatab_srgb[x]; }
	static inline float  gamma            (int x) { return gammatab[x]; }
	static inline float  igamma_srgb      (int x) { return igammatab_srgb[x]; }
	static inline float  gamma_srgb       (float x) { return gammatab_srgb[x]; }
	static inline float  gamma            (float x) { return gammatab[x]; }
	static inline float  igamma_srgb      (float x) { return igammatab_srgb[x]; }
	//static inline float  gamma_srgb       (double x) { return gammatab_srgb[x]; }
	//static inline float  gamma            (double x) { return gammatab[x]; }
	//static inline float  igamma_srgb      (double x) { return igammatab_srgb[x]; }

	//void gamutmap(LabImage* );
	static void gamutmap(float &X, float &Y, float &Z, const double p[3][3]);

	static inline float f2xyz(float f) {
		const float epsilonExpInv3 = 6.0/29.0;
		const float kappaInv = 27.0/24389.0;  // inverse of kappa

		return (f > epsilonExpInv3) ? f*f*f : (116 * f - 16) * kappaInv;
	}

};

}

#endif
