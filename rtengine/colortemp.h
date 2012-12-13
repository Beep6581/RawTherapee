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
#ifndef _COLORTEMP_
#define _COLORTEMP_

#include <gtkmm.h>
#include <cmath>
#include "colortemp.h"
#include "imagefloat.h"
#include "labimage.h"
#define MAXR(a,b) ((a) > (b) ? (a) : (b))

namespace rtengine {

#define MINTEMP 2000
#define MAXTEMP 25000
#define MINGREEN 0.02
#define MAXGREEN 5.0
#define INITIALBLACKBODY 4000

class ColorTemp {

	private:
		double temp;
		double green;
		Glib::ustring method;

		static void clip (double &temp, double &green);

	public:

		ColorTemp () : temp(-1), green(-1), method("Custom") {}
		ColorTemp (double t, double g, Glib::ustring m);
		ColorTemp (double mulr, double mulg, double mulb);

		inline double getTemp ()    { return temp;  }
		inline double getGreen ()   { return green; }

		void   getMultipliers (double &mulr, double &mulg, double &mulb) { temp2mul (temp, green, mulr, mulg, mulb); }

		void mul2temp (double rmul, double gmul, double bmul, double& temp, double& green);
		void temp2mul (double temp, double green, double& rmul, double& gmul, double& bmul);
		static void temp2mulxyz (double tem, double gree, Glib::ustring method, double &Xxyz, double &Zxyz);

		int XYZtoCorColorTemp(double x0,double y0 ,double z0, double &temp);
		static void cieCAT02(double Xw, double Yw, double Zw,double &CAM02BB00,double &CAM02BB01,double &CAM02BB02, double &CAM02BB10,double &CAM02BB11,double &CAM02BB12,double &CAM02BB20,double &CAM02BB21,double &CAM02BB22, double adap );
		//static	void CAT02 (Imagefloat* baseImg, const ProcParams* params);
		//static void ciecam_02 (LabImage* lab, const ProcParams* params);

		static double d_factor( double f, double la ) {
			return f * (1.0 - ((1.0 / 3.6) * exp((-la - 42.0) / 92.0)));
		}

		static double calculate_fl_from_la_ciecam02( double la ) {
			double la5 = la * 5.0;
			double k = 1.0 / (la5 + 1.0);

			/* Calculate k^4. */
			k = k * k;
			k = k * k;

			return (0.2 * k * la5) + (0.1 * (1.0 - k) * (1.0 - k) * pow(la5, 1.0 / 3.0));
		}

		static double achromatic_response_to_white( double x, double y, double z, double d, double fl, double nbb, int gamu ) {
			double r, g, b;
			double rc, gc, bc;
			double rp, gp, bp;
			double rpa, gpa, bpa;
			gamu=1;
			xyz_to_cat02( r, g, b, x, y, z, gamu );

			rc = r * (((y * d) / r) + (1.0 - d));
			gc = g * (((y * d) / g) + (1.0 - d));
			bc = b * (((y * d) / b) + (1.0 - d));

			cat02_to_hpe( rp, gp, bp, rc, gc, bc, gamu );
			if(gamu==1){//gamut correction M.H.Brill S.Susstrunk
				rp=MAXR(rp,0.0);
				gp=MAXR(gp,0.0);
				bp=MAXR(bp,0.0);
			}

			rpa = nonlinear_adaptation( rp, fl );
			gpa = nonlinear_adaptation( gp, fl );
			bpa = nonlinear_adaptation( bp, fl );

			return ((2.0 * rpa) + gpa + ((1.0 / 20.0) * bpa) - 0.305) * nbb;
		}

		static void xyz_to_cat02 ( double &r,  double &g,  double &b,  double x, double y, double z, int gamu );
		static void cat02_to_hpe ( double &rh, double &gh, double &bh, double r, double g, double b, int gamu );
		static void cat02_to_xyz ( double &x,  double &y,  double &z,  double r, double g, double b, int gamu );
		static void hpe_to_xyz   ( double &x,  double &y,  double &z,  double r, double g, double b );

		static double nonlinear_adaptation( double c, double fl ) {
		double p;
			if(c<0.0){ p = pow( (-1.0*fl * c) / 100.0, 0.42 );return ((-1.0*400.0 * p) / (27.13 + p)) + 0.1;}
			else {p = pow( (fl * c) / 100.0, 0.42 ); return ((400.0 * p) / (27.13 + p)) + 0.1;}
		}

		static double inverse_nonlinear_adaptation( double c, double fl ) {
		    int c1;
		    if(c-0.1 < 0.0) c1=-1; else c1=1;
			return c1*(100.0 / fl) * pow( (27.13 * fabs( c - 0.1 )) / (400.0 - fabs( c - 0.1 )), 1.0 / 0.42 );
		}
		static void curvecolor(double satind, double satval, double &sres, double parsat); 				
		static void curveJ (double br, double contr, int db, LUTf & outCurve , LUTu & histogram ) ;

		bool operator== (const ColorTemp& other) { return fabs(temp-other.temp)<1e-10 && fabs(green-other.green)<1e-10; }
		bool operator!= (const ColorTemp& other) { return !(*this==other); }

		static double blackbody_spect (double wavelength, double m1, double m2, double temp);
		static double daylight_spect  (double wavelength, double m1, double m2, double temp);
		static const double Cloudy6200_spect[97];
		static const double Daylight5300_spect[97];
		static const double Shade7600_spect[97];
		static const double A2856_spect[97];
		static const double FluoF1_spect[97];
		static const double FluoF2_spect[97];
		static const double FluoF3_spect[97];
		static const double FluoF4_spect[97];
		static const double FluoF5_spect[97];
		static const double FluoF6_spect[97];
		static const double FluoF7_spect[97];
		static const double FluoF8_spect[97];
		static const double FluoF9_spect[97];
		static const double FluoF10_spect[97];
		static const double FluoF11_spect[97];
		static const double FluoF12_spect[97];
		static const double HMI_spect[97];
		static const double GTI_spect[97];
		static const double JudgeIII_spect[97];
		static const double Solux3500_spect[97];
		static const double Solux4100_spect[97];
		static const double Solux4700_spect[97];
		static const double NG_Solux4700_spect[97];
		static const double NG_LEDLSI2040_spect[97];
		static const double NG_CRSSP12WWMR16_spect[97];
		static const double Flash5500_spect[97];
		static const double Flash6000_spect[97];
		static const double Flash6500_spect[97];

		//spectral data 8 color Colorchecker24
		static double get_spectral_color (double wavelength, const double* array) {
			int wlm = (int) ((wavelength -350.)/5.);
			return (array[wlm]);
		}

		static const double ColorchechredC3_spect[97];
		static const double ColorchechOraA2_spect[97];
		static const double ColorchechYelD3_spect[97];
		static const double ColorchechGreE2_spect[97];
		static const double ColorchechGreB3_spect[97];
		static const double ColorchechCyaF3_spect[97];
		static const double ColorchechPurD2_spect[97];
		static const double ColorchechMagE3_spect[97];
		static const double ColorchechSkiA138_13_14_spect[97];
		static const double ColorchechGraC4_67_spect[97];
		static const double ColorchechSkiB166_18_18_spect[97];
		static const double ColorchechBluC150_m5_m22_spect[97];
		static const double ColorchechDCBluN881_m7_m14_spect[97];//ColorChecker DC  N8
		static const double ColorchechSGSkiF763_14_26_spect[97];//ColorChecker SG  F7
		static const double ColorchechSGSkiK285_11_17_spect[97];//ColorChecker SG  K2
		static const double ColorchechWhiA496_spect[97];
		static const double ColorchechGreD1_spect[97];
		static const double ColorchechSGBlaN3_6_spect[97];//ColorChecker SG  N3
		static const double JDC468_GraK14_44_spect[97];//468  K14
		static const double JDC468_BluH10_spect[97];//468  H10
		static const double ColabSkin35_15_17_spect[97];//Skin L 35
		static const double ColabSkin57_22_18_spect[97];//Skin L 57
		static const double ColabSkin40_17_17_spect[97];//Skin L 40
		static const double ColabSkin91_4_14_spect[97];//Skin L 91
		static const double ColabSkin87_8_8_spect[97];//Skin L 87
		static const double ColabSkin89_8_21_spect[97];//Skin L 89
		static const double ColabSkin75_8_4_spect[97];//Skin L 75
		static const double ColabSkin75_10_33_spect[97];//Skin L 75
		static const double ColabSkin65_33_11_spect[97];//Skin L 65
		static const double ColabSkin65_7_24_spect[97];//Skin L 65
		static const double ColabSkin57_19_6_spect[97];//Skin L 57
		static const double ColabSkin57_4_19_spect[97];//Skin L 57
		static const double ColabSkin57_10_28_spect[97];//Skin L 57
		static const double ColabSkin40_7_19_spect[97];//Skin L 57
		static const double ColabSkin40_17_6_spect[97];//Skin L 40
		static const double ColabSkin40_4_11_spect[97];//Skin L 40
		static const double ColabSkin33_6_15_spect[97];//Skin L 33
		static const double ColabSkin33_15_5_spect[97];//Skin L 33
		static const double ColabSkin33_10_15_spect[97];//Skin L 33
		static const double ColabSkin24_5_6_spect[97];//Skin L 24
		static const double ColabSkin26_18_18_spect[97];//Skin L 26
		static const double ColabSkin24_7_5_spect[97];//Skin L 24
		static const double ColabSkin20_4_2_spect[97];//Skin L 20
		static const double ColabSkin98_m2_10_spect[97];//Skin L 98
		static const double ColabSkin90_m1_20_spect[97];//Skin L 90
		static const double ColabSkin95_0_4_spect[97];//Skin L 95
		static const double ColabSkin81_2_14_spect[97];//Skin L 81
		static const double ColabSkin87_3_10_spect[97];//Skin L 87
		static const double ColabSkin77_12_21_spect[97];//Skin L 77
		static const double ColabSkin70_7_32_spect[97];//Skin L 77
		static const double ColabSky60_0_m31_spect[97];//Sky L=60
		static const double ColabSky42_0_m24_spect[97];//Sky L=42

		static void spectrum_to_xyz_daylight  (double _m1, double _m2, double _temp, double &x, double &y, double &z);
		static void spectrum_to_xyz_blackbody (double _m1, double _m2, double _temp, double &x, double &y, double &z);
		static void spectrum_to_xyz_preset    (const double* spec_intens, double &x, double &y, double &z);

		static void spectrum_to_color_xyz_daylight  (const double* spec_color, double _m1, double _m2, double _temp, double &xx, double &yy, double &zz);
		static void spectrum_to_color_xyz_blackbody (const double* spec_color, double _m1, double _m2, double _temp, double &xx, double &yy, double &zz);
		static void spectrum_to_color_xyz_preset    (const double* spec_color, const double* spec_intens, double &xx, double &yy, double &zz);



/**
 * Inverse transform from CIECAM02 JCh to XYZ.
 */
		static void jch2xyz_ciecam02( double &x, double &y, double &z,
		                              double J, double C, double h,
		                              double xw, double yw, double zw,
		                              double yb, double la,
		                              double f, double c, double nc, bool doneinit2, int gamu );
/**
 * Forward transform from XYZ to CIECAM02 JCh.
 */
		static void xyz2jchqms_ciecam02( double &J, double &C, double &h,
		                                 double &Q, double &M, double &s,double &aw, double &fl, double &wh,
		                                 double x, double y, double z,
		                                 double xw, double yw, double zw,
		                                 double yb, double la,
		                                 double f, double c, double nc,  double pilotd, bool doneinit1, int gamu  );

};
inline void ColorTemp::xyz_to_cat02( double &r, double &g, double &b, double x, double y, double z, int gamu )
{
gamu=1;
if(gamu==0){
	r = ( 0.7328 * x) + (0.4296 * y) - (0.1624 * z);
	g = (-0.7036 * x) + (1.6975 * y) + (0.0061 * z);
	b = ( 0.0030 * x) + (0.0136 * y) + (0.9834 * z);
	}
else if (gamu==1) {//gamut correction M.H.Brill S.Susstrunk
	//r = ( 0.7328 * x) + (0.4296 * y) - (0.1624 * z);
	//g = (-0.7036 * x) + (1.6975 * y) + (0.0061 * z);
	//b = ( 0.0000 * x) + (0.0000 * y) + (1.0000 * z);
	r = ( 1.007245 * x) + (0.011136* y) - (0.018381 * z);//Changjun Li
	g = (-0.318061 * x) + (1.314589 * y) + (0.003471 * z);
	b = ( 0.0000 * x) + (0.0000 * y) + (1.0000 * z);
	
}
}

inline void ColorTemp::cat02_to_xyz( double &x, double &y, double &z, double r, double g, double b, int gamu )
{
gamu=1;
if(gamu==0) {
	x = ( 1.096124 * r) - (0.278869 * g) + (0.182745 * b);
	y = ( 0.454369 * r) + (0.473533 * g) + (0.072098 * b);
	z = (-0.009628 * r) - (0.005698 * g) + (1.015326 * b);
	}
else if(gamu==1) {//gamut correction M.H.Brill S.Susstrunk
	//x = ( 1.0978566 * r) - (0.277843 * g) + (0.179987 * b);
	//y = ( 0.455053 * r) + (0.473938 * g) + (0.0710096* b);
	//z = ( 0.000000 * r) - (0.000000 * g) + (1.000000 * b);
	x = ( 0.99015849 * r) - (0.00838772* g) + (0.018229217 * b);//Changjun Li
	y = ( 0.239565979 * r) + (0.758664642 * g) + (0.001770137* b);
	z = ( 0.000000 * r) - (0.000000 * g) + (1.000000 * b);
	
}	
}

inline void ColorTemp::hpe_to_xyz( double &x, double &y, double &z, double r, double g, double b )
{
	x = (1.910197 * r) - (1.112124 * g) + (0.201908 * b);
	y = (0.370950 * r) + (0.629054 * g) - (0.000008 * b);
	z = b;
	
}




inline void ColorTemp::cat02_to_hpe( double &rh, double &gh, double &bh, double r, double g, double b, int gamu )
{ gamu=1; 
 if(gamu==0){
	rh = ( 0.7409792 * r) + (0.2180250 * g) + (0.0410058 * b);
	gh = ( 0.2853532 * r) + (0.6242014 * g) + (0.0904454 * b);
	bh = (-0.0096280 * r) - (0.0056980 * g) + (1.0153260 * b);
	}
	else if (gamu==1) {//Changjun Li
	rh = ( 0.550930835 * r) + (0.519435987* g) - ( 0.070356303* b);
	gh = ( 0.055954056 * r) + (0.89973132 * g) + (0.044315524 * b);
	bh = (0.0 * r) - (0.0* g) + (1.0 * b);
	}
}

}
#endif
