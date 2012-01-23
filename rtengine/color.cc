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

#include "color.h"
#include "iccmatrices.h"

namespace rtengine {

#undef MAXVAL
#undef MAX
#undef MIN

#define MAXVAL  0xffff
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define SQR(x) ((x)*(x))

#define eps_max 580.40756 //(MAXVAL* 216.0f/24389.0);
#define kappa	903.29630 //24389.0/27.0;

LUTf Color::cachef ;
LUTf Color::gamma2curve = 0;

LUTf Color::gammatab;
LUTf Color::igammatab_srgb;
LUTf Color::gammatab_srgb;

// Wikipedia sRGB: Unlike most other RGB color spaces, the sRGB gamma cannot be expressed as a single numerical value.
// The overall gamma is approximately 2.2, consisting of a linear (gamma 1.0) section near black, and a non-linear section elsewhere involving a 2.4 exponent
// and a gamma (slope of log output versus log input) changing from 1.0 through about 2.3.
const double Color::sRGBGamma = 2.2;
const double Color::sRGBGammaCurve = 2.4;

void Color::init () {

	int maxindex = 65536;
	cachef(maxindex,0/*LUT_CLIP_BELOW*/);

	gamma2curve(maxindex,0);

	for (int i=0; i<maxindex; i++) {
		if (i>eps_max) {
			cachef[i] = 327.68*( exp(1.0/3.0 * log((double)i / MAXVAL) ));
		}
		else {
			cachef[i] = 327.68*((kappa*i/MAXVAL+16.0)/116.0);
		}
	}

	for (int i=0; i<maxindex; i++) {
		gamma2curve[i] = (gamma2(i/65535.0) * 65535.0);
	}

	/*******************************************/

	gammatab(65536,0);
	igammatab_srgb(65536,0);
	gammatab_srgb(65536,0);

	for (int i=0; i<65536; i++)
	gammatab_srgb[i] = (65535.0 * gamma2 (i/65535.0));
	for (int i=0; i<65536; i++)
	igammatab_srgb[i] = (65535.0 * igamma2 (i/65535.0));
	for (int i=0; i<65536; i++)
	gammatab[i] = (65535.0 * pow (i/65535.0, 0.454545));

	/*FILE* f = fopen ("c.txt", "wt");
	for (int i=0; i<256; i++)
		fprintf (f, "%g %g\n", i/255.0, clower (i/255.0, 2.0, 1.0));
	fclose (f);*/

}

void Color::cleanup () {

}

void Color::rgb2hsv (float r, float g, float b, float &h, float &s, float &v) {

	double var_R = r / 65535.0;
	double var_G = g / 65535.0;
	double var_B = b / 65535.0;

	double var_Min = MIN(MIN(var_R,var_G),var_B);
	double var_Max = MAX(MAX(var_R,var_G),var_B);
	double del_Max = var_Max - var_Min;
	v = var_Max;
	if (fabs(del_Max)<0.00001) {
		h = 0;
		s = 0;
	}
	else {
		s = del_Max/var_Max;

		if      ( var_R == var_Max ) h = (var_G - var_B)/del_Max;
		else if ( var_G == var_Max ) h = 2.0 + (var_B - var_R)/del_Max;
		else if ( var_B == var_Max ) h = 4.0 + (var_R - var_G)/del_Max;
		h /= 6.0;

		if ( h < 0 )  h += 1;
		if ( h > 1 )  h -= 1;
	}
}

void Color::rgb2hsv (int r, int g, int b, float &h, float &s, float &v) {

	double var_R = r / 65535.0;
	double var_G = g / 65535.0;
	double var_B = b / 65535.0;

	double var_Min = MIN(MIN(var_R,var_G),var_B);
	double var_Max = MAX(MAX(var_R,var_G),var_B);
	double del_Max = var_Max - var_Min;
	v = var_Max;
	if (fabs(del_Max)<0.00001) {
		h = 0;
		s = 0;
	}
	else {
		s = del_Max/var_Max;

		if      ( var_R == var_Max ) h = (var_G - var_B)/del_Max;
		else if ( var_G == var_Max ) h = 2.0 + (var_B - var_R)/del_Max;
		else if ( var_B == var_Max ) h = 4.0 + (var_R - var_G)/del_Max;
		h /= 6.0;

		if ( h < 0 )  h += 1;
		if ( h > 1 )  h -= 1;
	}
}

void Color::hsv2rgb (float h, float s, float v, float &r, float &g, float &b) {

	float h1 = h*6; // sector 0 to 5
	int i = floor( h1 );
	float f = h1 - i; // fractional part of h

	float p = v * ( 1 - s );
	float q = v * ( 1 - s * f );
	float t = v * ( 1 - s * ( 1 - f ) );

	float r1,g1,b1;

	if      (i==0) {r1 = v;  g1 = t;  b1 = p;}
	else if (i==1) {r1 = q;  g1 = v;  b1 = p;}
	else if (i==2) {r1 = p;  g1 = v;  b1 = t;}
	else if (i==3) {r1 = p;  g1 = q;  b1 = v;}
	else if (i==4) {r1 = t;  g1 = p;  b1 = v;}
	else if (i==5) {r1 = v;  g1 = p;  b1 = q;}

	r = ((r1)*65535.0);
	g = ((g1)*65535.0);
	b = ((b1)*65535.0);
}


// The same function but set float values instead if int
// Function copied for speed concerns
/* Not exactly the same as above ; this one return a result in the [0.0 ; 1.0] range
void Color::hsv2rgb (float h, float s, float v, float &r, float &g, float &b) {

	float h1 = h*6; // sector 0 to 5
	int i = floor( h1 );
	float f = h1 - i; // fractional part of h

	float p = v * ( 1 - s );
	float q = v * ( 1 - s * f );
	float t = v * ( 1 - s * ( 1 - f ) );

	if      (i==0) {r = v;  g = t;  b = p;}
	else if (i==1) {r = q;  g = v;  b = p;}
	else if (i==2) {r = p;  g = v;  b = t;}
	else if (i==3) {r = p;  g = q;  b = v;}
	else if (i==4) {r = t;  g = p;  b = v;}
	else if (i==5) {r = v;  g = p;  b = q;}
}
*/

void Color::hsv2rgb (float h, float s, float v, int &r, int &g, int &b) {

	float h1 = h*6; // sector 0 to 5
	int i = floor( h1 );
	float f = h1 - i; // fractional part of h

	float p = v * ( 1 - s );
	float q = v * ( 1 - s * f );
	float t = v * ( 1 - s * ( 1 - f ) );

	float r1,g1,b1;

	if (i==0) {r1 = v;  g1 = t;  b1 = p;}
	else if (i==1) {r1 = q;  g1 = v;  b1 = p;}
	else if (i==2) {r1 = p;  g1 = v;  b1 = t;}
	else if (i==3) {r1 = p;  g1 = q;  b1 = v;}
	else if (i==4) {r1 = t;  g1 = p;  b1 = v;}
	else if (i==5) {r1 = v;  g1 = p;  b1 = q;}

	r = (int)( r1 * 65535);
	g = (int)( g1 * 65535);
	b = (int)( b1 * 65535);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void Color::xyz2srgb (float x, float y, float z, float &r, float &g, float &b) {

	//Transform to output color.  Standard sRGB is D65, but internal representation is D50
	//Note that it is only at this point that we should have need of clipping color data

	/*float x65 = d65_d50[0][0]*x + d65_d50[0][1]*y + d65_d50[0][2]*z ;
	 float y65 = d65_d50[1][0]*x + d65_d50[1][1]*y + d65_d50[1][2]*z ;
	 float z65 = d65_d50[2][0]*x + d65_d50[2][1]*y + d65_d50[2][2]*z ;

	 r = sRGB_xyz[0][0]*x65 + sRGB_xyz[0][1]*y65 + sRGB_xyz[0][2]*z65;
	 g = sRGB_xyz[1][0]*x65 + sRGB_xyz[1][1]*y65 + sRGB_xyz[1][2]*z65;
	 b = sRGB_xyz[2][0]*x65 + sRGB_xyz[2][1]*y65 + sRGB_xyz[2][2]*z65;*/

	/*r = sRGBd65_xyz[0][0]*x + sRGBd65_xyz[0][1]*y + sRGBd65_xyz[0][2]*z ;
	 g = sRGBd65_xyz[1][0]*x + sRGBd65_xyz[1][1]*y + sRGBd65_xyz[1][2]*z ;
	 b = sRGBd65_xyz[2][0]*x + sRGBd65_xyz[2][1]*y + sRGBd65_xyz[2][2]*z ;*/

	r = ((sRGB_xyz[0][0]*x + sRGB_xyz[0][1]*y + sRGB_xyz[0][2]*z)) ;
	g = ((sRGB_xyz[1][0]*x + sRGB_xyz[1][1]*y + sRGB_xyz[1][2]*z)) ;
	b = ((sRGB_xyz[2][0]*x + sRGB_xyz[2][1]*y + sRGB_xyz[2][2]*z)) ;

}


void Color::xyz2rgb (float x, float y, float z, float &r, float &g, float &b, float rgb_xyz[3][3]) {

	//Transform to output color.  Standard sRGB is D65, but internal representation is D50
	//Note that it is only at this point that we should have need of clipping color data

	/*float x65 = d65_d50[0][0]*x + d65_d50[0][1]*y + d65_d50[0][2]*z ;
	 float y65 = d65_d50[1][0]*x + d65_d50[1][1]*y + d65_d50[1][2]*z ;
	 float z65 = d65_d50[2][0]*x + d65_d50[2][1]*y + d65_d50[2][2]*z ;

	 r = sRGB_xyz[0][0]*x65 + sRGB_xyz[0][1]*y65 + sRGB_xyz[0][2]*z65;
	 g = sRGB_xyz[1][0]*x65 + sRGB_xyz[1][1]*y65 + sRGB_xyz[1][2]*z65;
	 b = sRGB_xyz[2][0]*x65 + sRGB_xyz[2][1]*y65 + sRGB_xyz[2][2]*z65;*/

	/*r = sRGBd65_xyz[0][0]*x + sRGBd65_xyz[0][1]*y + sRGBd65_xyz[0][2]*z ;
	 g = sRGBd65_xyz[1][0]*x + sRGBd65_xyz[1][1]*y + sRGBd65_xyz[1][2]*z ;
	 b = sRGBd65_xyz[2][0]*x + sRGBd65_xyz[2][1]*y + sRGBd65_xyz[2][2]*z ;*/



	r = ((rgb_xyz[0][0]*x + rgb_xyz[0][1]*y + rgb_xyz[0][2]*z)) ;
	g = ((rgb_xyz[1][0]*x + rgb_xyz[1][1]*y + rgb_xyz[1][2]*z)) ;
	b = ((rgb_xyz[2][0]*x + rgb_xyz[2][1]*y + rgb_xyz[2][2]*z)) ;

}

void Color::calcGamma (double pwr, double ts, int mode, int imax, double &gamma0, double &gamma1, double &gamma2, double &gamma3, double &gamma4, double &gamma5) {
//from Dcraw (D.Coffin)
	int i;
	double g[6], bnd[2]={0,0}, r;

	g[0] = pwr;
	g[1] = ts;
	g[2] = g[3] = g[4] = 0;
	bnd[g[1] >= 1] = 1;
	if (g[1] && (g[1]-1)*(g[0]-1) <= 0) {
		for (i=0; i < 48; i++) {
			g[2] = (bnd[0] + bnd[1])/2;
			if (g[0])
				bnd[(pow(g[2]/g[1],-g[0]) - 1)/g[0] - 1/g[2] > -1] = g[2];
			else
				bnd[g[2]/exp(1-1/g[2]) < g[1]] = g[2];
		}
		g[3] = g[2] / g[1];
		if (g[0]) g[4] = g[2] * (1/g[0] - 1);
	}
	if (g[0])
		g[5] = 1 / (g[1]*SQR(g[3])/2 - g[4]*(1 - g[3]) + (1 - pow(g[3],1+g[0]))*(1 + g[4])/(1 + g[0])) - 1;
	else
		g[5] = 1 / (g[1]*SQR(g[3])/2 + 1 - g[2] - g[3] - g[2]*g[3]*(log(g[3]) - 1)) - 1;
	if (!mode--) {
		gamma0=g[0];gamma1=g[1];gamma2=g[2];gamma3=g[3];gamma4=g[4];gamma5=g[5];
		return;
	}
}

void Color::Lab2XYZ(float L, float a, float b, float &x, float &y, float &z) {
	float fy = (0.00862069 * L) + 0.137932; // (L+16)/116
	float fx = (0.002 * a) + fy;
	float fz = fy - (0.005 * b);

	x = 65535.0*f2xyz(fx)*D50x;
	y = 65535.0*f2xyz(fy);
	z = 65535.0*f2xyz(fz)*D50z;
}

void Color::XYZ2Lab(float X, float Y, float Z, float &L, float &a, float &b) {

	float X1 = X/D50x;
	float Z1 = Z/D50z;

	float fx = (X1<65535.0 ? cachef[X1] : (327.68*exp(log(X1/MAXVAL)/3.0 )));
	float fy = (Y<65535.0 ? cachef[Y] : (327.68*exp(log(Y/MAXVAL)/3.0 )));
	float fz = (Z1<65535.0 ? cachef[Z1] : (327.68*exp(log(Z1/MAXVAL)/3.0 )));

	L = (116.0 * fy - 5242.88); //5242.88=16.0*327.68;
	a = (500.0 * (fx - fy) );
	b = (200.0 * (fy - fz) );

}

void Color::Lab2Yuv(float L, float a, float b, float &Y, float &u, float &v) {
	float fy = (0.00862069 * L/327.68) + 0.137932; // (L+16)/116
	float fx = (0.002 * a/327.68) + fy;
	float fz = fy - (0.005 * b/327.68);

	float X = 65535.0*f2xyz(fx)*D50x;
	Y = 65535.0*f2xyz(fy);
	float Z = 65535.0*f2xyz(fz)*D50z;

	u = 4.0*X/(X+15*Y+3*Z)-u0;
	v = 9.0*Y/(X+15*Y+3*Z)-v0;

}

void Color::Yuv2Lab(float Yin, float u, float v, float &L, float &a, float &b, double wp[3][3]) {

	float u1 = u + u0;
	float v1 = v + v0;

	float Y = Yin;
	float X = (9*u1*Y)/(4*v1*D50x);
	float Z = (12 - 3*u1 - 20*v1)*Y/(4*v1*D50z);

	gamutmap(X,Y,Z,wp);

	float fx = (X<65535.0 ? cachef[X] : (327.68*exp(log(X/MAXVAL)/3.0 )));
	float fy = (Y<65535.0 ? cachef[Y] : (327.68*exp(log(Y/MAXVAL)/3.0 )));
	float fz = (Z<65535.0 ? cachef[Z] : (327.68*exp(log(Z/MAXVAL)/3.0 )));

	L = (116.0 * fy - 5242.88); //5242.88=16.0*327.68;
	a = (500.0 * (fx - fy) );
	b = (200.0 * (fy - fz) );

}

}
